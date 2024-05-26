mod pearson;
mod rm_lowentropy;
use rm_lowentropy::{cacu_by_entropy_gene, cacu_index_by_sort, core_threads_rm_low};
// mod spearman;
mod readfile;
use anyhow::{anyhow, Error, Ok, Result};
// use ndarray::{Array, Axis};
use pearson::core_threads_cacu_pearson;
use readfile::{read_expression_matrix, read_group_file, read_tpm_matrix};
use std::vec;
use std::{collections::HashMap, path::Path};
mod writerfile;
use writerfile::{
    generate_hist_file_name, write_2d_vec_to_file, writer_entropy, writer_entropy_from_hash,
    writer_filter_matrix, writer_final_hash,
};
// 计算运行时间
use std::time::Instant;
mod caculate_core;
// mod plot;
use caculate_core::{
    calculate_percentiles_for_pcc, generate_pcc_hashmap, remove_low_entropy_values,
};
// use plot::plot_histogram;
mod args;
use ansi_term::Colour;
use args::get_arg;
use rayon::prelude::*;

use crate::readfile::read_gz_file;

#[tokio::main]
async fn main() -> Result<(), Error> {
    let start_time = Instant::now();
    let matches = get_arg();

    let input = matches
        .value_of("input")
        .ok_or_else(|| anyhow!("No matrix input file specified"))?;
    let samples = matches
        .value_of("samples")
        .ok_or_else(|| anyhow!("No samples input file specified"))?;
    let output = matches.value_of("output").unwrap_or("cst_gene_net_out");
    let pearson_threshold = matches
        .value_of("pearson_threshold")
        .unwrap_or("0.8")
        .parse::<f64>()?;
    // 判断output是否存在不存在则创建
    let _ = std::fs::create_dir_all(output);
    let tpm_path = Path::new(input);
    let path = Path::new(samples);
    let _ = pearson_correlation_group(path, tpm_path, output, false, pearson_threshold).await?;
    let elapsed = start_time.elapsed();
    println!("{}: {:#?}", Colour::Green.paint("耗时"), elapsed);
    println!(
        "{}{}",
        Colour::Green.paint("Done!"),
        Colour::Red.paint("Goodbye!")
    );
    Ok(())
}

//根据分组计算各自分组的皮尔逊相关系数

async fn pearson_correlation_group(
    samples_path: &Path,
    matrix_path: &Path,
    output: &str,
    use_nd: bool,
    pearson_threshold: f64,
) -> Result<(), Error> {
    // let OUT_DIR = output;
    // 判断output/03_gene_matrix.csv是否存在如果存在
    // 则读取
    let (group_map, samples_tag) = read_group_file(samples_path)?;
    let (expression_vec, samples, gene_vec) = read_tpm_matrix(matrix_path, output)?;

    // println!("{:?}", group_map);
    // println!("{:?}", samples);
    let mut group_uniq_map = HashMap::new();
    // group_map分组
    for group in group_map.keys() {
        // println!("{}", group);
        let value = group_map
            .get(group)
            .ok_or_else(|| anyhow!("group not found"))?;
        // 计算group在samples中的index
        // println!("{:?}", value);
        let group_index = samples
            .iter()
            .position(|x| x == group)
            .ok_or_else(|| anyhow!("group value not found"))?;
        // println!("{:?}", group_index);
        if !group_uniq_map.contains_key(value) {
            group_uniq_map.insert(value.to_string(), vec![group_index]);
        } else {
            let mut tmp_value_vec = group_uniq_map.get_mut(value).unwrap().to_vec();
            tmp_value_vec.push(group_index);
            group_uniq_map.insert(value.to_string(), tmp_value_vec);
        }
    }
    let mut new_expression_vec = Vec::new();
    let mut new_gene_vec = Vec::new();

    if Path::new(&format!("{}/03_gene_matrix.csv", output)).exists() {
        new_expression_vec = expression_vec;
        new_gene_vec = gene_vec;
    } else {
        let entr_result = core_threads_rm_low(&group_uniq_map, &expression_vec, &gene_vec)?;
        let _ = writer_entropy(&entr_result, &format!("{}/01_entropy_value.tsv", output)).await?;
        println!(
            "{}:{}",
            Colour::Green.paint("Writer entropy value file"),
            Colour::Red.paint("01_entropy_value.tsv")
        );
        let entr_res_hash = cacu_index_by_sort(&entr_result)?;
        // 写出每种分组下的entropy值
        let _ = writer_entropy_from_hash(
            &entr_res_hash,
            &format!("{}/02_entropy_sorted_min_2_max.csv", output),
        )
        .await?;
        println!(
            "{}:{}",
            Colour::Green.paint("Writer entropy sort index file"),
            Colour::Red.paint("02_entropy_sorted_min_2_max.csv")
        );
        let entr_gene_hashmap = cacu_by_entropy_gene(&gene_vec, &entr_res_hash)?;
        // println!("{:#?}", entr_gene_hashmap);
        (new_expression_vec, new_gene_vec) =
            remove_low_entropy_values(&entr_gene_hashmap, &expression_vec, &gene_vec)?;
        // 再次写出表达矩阵
        let _ = writer_filter_matrix(
            &new_expression_vec,
            &new_gene_vec,
            &format!("{}/03_gene_matrix.csv", output),
            &samples,
        )
        .await?;
        println!(
            "{}:{}",
            Colour::Green.paint("Writer rm low entropy gene matrix file"),
            Colour::Red.paint("03_gene_matrix.csv")
        );
    }

    let res_vec: Vec<String> = core_threads_cacu_pearson(
        &group_uniq_map,
        new_expression_vec,
        &new_gene_vec,
        &samples,
        &samples_tag,
        use_nd,
        output,
        pearson_threshold,
    )?;
    // println!("{:#?}", res_vec);
    // 读取res_vec中的文件
    // res_vec.par_iter().for_each(|res| {
    //     let mut pearson_vec: Vec<f64> = Vec::new();
    //     if let std::result::Result::Ok(res_file) = read_gz_file(Path::new(&res)) {
    //         for line in res_file {
    //             if let std::result::Result::Ok(line) = line {
    //                 // 跳过空行
    //                 if line.trim().is_empty() {
    //                     continue;
    //                 } else {
    //                     pearson_vec.push(
    //                         line.trim().split("\t").collect::<Vec<&str>>()[2]
    //                             .parse::<f64>()
    //                             .expect("pearson value parse error"),
    //                     )
    //                 }
    //             }
    //         }
    //     }
    //     let pcc_min_max = calculate_percentiles_for_pcc(&pearson_vec).unwrap();
    //     let min5 = pcc_min_max[0];
    //     let max5 = pcc_min_max[1];
    //     println!(
    //         "{} group:{} min:{} max:{}",
    //         Colour::Green.paint("Caculate min pcc and max pcc"),
    //         res,
    //         min5,
    //         max5
    //     );
    //     // 最终结果文件
    //     let pearson_res = res.replace("rawPCC.txt", "resultPCC.txt");
    //     let mut res_vec_pcc = Vec::new();
    //     // println!("{}", res);
    //     if let std::result::Result::Ok(res_file) = read_gz_file(Path::new(&res)) {
    //         for line in res_file {
    //             if let std::result::Result::Ok(line) = line {
    //                 if line.trim().is_empty() {
    //                     continue;
    //                 }
    //                 let lines = line.trim().split("\t").collect::<Vec<&str>>(); // 提取文件
    //                 let gene_name_1 = lines[0];
    //                 let gene_name_2 = lines[1];
    //                 let pearson_value = lines[2].parse::<f64>().expect("");
    //                 let p_value = lines[3].parse::<f64>().expect("pearson value parse error");
    //                 let mut note = "";
    //                 if pearson_value <= min5 {
    //                     note = "negative";
    //                 } else if pearson_value >= max5 {
    //                     note = "positive";
    //                 } else {
    //                     note = "neutral";
    //                     continue;
    //                 }
    //                 res_vec_pcc.push(vec![
    //                     gene_name_1.to_string(),
    //                     gene_name_2.to_string(),
    //                     pearson_value.to_string(),
    //                     p_value.to_string(),
    //                     note.to_string(),
    //                 ]);
    //             }
    //         }
    //     }
    //     let _ = write_2d_vec_to_file(res_vec_pcc, &pearson_res, &"\t")
    //         .expect("write pearson result is error!");
    // });

    println!(
        "{} {}",
        Colour::Green.paint("Pearson caculate down!"),
        Colour::Red.paint("lalalala~~~~")
    );
    println!(
        "{} {}",
        Colour::Green.paint("Caculate PCC down!"),
        Colour::Red.paint("Start writer result")
    );
    // let _ = writer_final_hash(res_hash).await?;
    println!("{}", Colour::Green.paint("Writer result down!"));
    Ok(())
}

#[tokio::test]
async fn test_correlation() {
    let start = Instant::now();
    let tpm_path = Path::new("test_data/test.csv.gz");
    let path = Path::new("test_data/01_tissue_sample.csv.gz");
    let _ = pearson_correlation_group(path, tpm_path, "test_data", false, 0.8)
        .await
        .unwrap();
    // println!("{:?}", res);
    // write_2d_vec_to_file(res, "test_data/pearson_test.csv").unwrap();
    let elapsed = start.elapsed();
    println!("耗时: {:?}", elapsed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_vec() {
        let mut my_vec: Vec<Vec<i32>> = vec![vec![1, 2, 3], vec![4, 5, 6]];

        // 需要添加到现有数组后方的新元素
        let new_elements: Vec<i32> = (7..=10).collect();
        // 使用 par_iter_mut() 并行地对最后一个子数组进行操作，逐个添加新元素
        // 将新元素包装成一个子数组，并用 par_extend 方法并行地添加到现有数组后方
        my_vec.extend(vec![new_elements]);

        println!("{:?}", my_vec);
    }
}
