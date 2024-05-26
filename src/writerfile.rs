use crate::rm_lowentropy::MixedNumber;
use anyhow::{anyhow, Error, Ok, Result};
use std::collections::HashMap;
// use std::fs::File;
// use std::io::{self, Write};
use rayon::prelude::*;
use std::io::Write;
use tokio::fs::File;
use tokio::io::{self, AsyncWriteExt};

// 写出第一个文件每个基因的entropy值
pub async fn writer_entropy(
    entropy: &HashMap<String, Vec<Vec<MixedNumber>>>,
    file_name: &str,
) -> Result<(), Error> {
    let mut writer = tokio::fs::File::create(file_name).await?;
    for key in entropy.keys() {
        let value = entropy
            .get(key)
            .ok_or_else(|| anyhow!("Could not find key"))?;
        for v in value {
            let mut gene_name = "";
            let mut entropy_value = &0.0;
            for i in v {
                match i {
                    MixedNumber::String(s) => gene_name = &s,
                    MixedNumber::Float(f) => entropy_value = f,
                    _ => continue,
                }
            }
            // writeln!(writer, "{}\t{}\t{:#?}", key, gene_name, entropy_value)?;
            writer
                .write_all(format!("{}\t{}\t{:#?}\n", key, gene_name, entropy_value).as_bytes())
                .await?;
        }
    }
    // writer.flush()?;

    Ok(())
}

// 写出第2个文件根据entropy排序的基因列表
pub async fn writer_entropy_from_hash(
    entropy_from_hash: &HashMap<String, Vec<usize>>,
    file_name: &str,
) -> Result<(), Error> {
    let mut writer = File::create(file_name).await?;
    for (key, value) in entropy_from_hash {
        writer.write_all(key.as_bytes()).await?;
        // writer.write_all(b",").await?;
        // writer.write_all(value.as_slice())?;
        // usize写出
        for i in value {
            writer.write_all(b",").await?;
            writer.write_all(i.to_string().as_bytes()).await?;
        }
        writer.write_all(b"\n").await?;
    }
    Ok(())
}
// 写出过滤后的表达矩阵
pub async fn writer_filter_matrix(
    matrix: &Vec<Vec<f64>>,
    gene_vec: &Vec<String>,
    file_name: &str,
    samples: &Vec<String>,
) -> Result<(), Error> {
    let mut writer = File::create(file_name).await?;
    writer.write_all(b"gene_name").await?;
    for sample in samples {
        writer.write_all(b",").await?;
        writer.write_all(sample.as_bytes()).await?;
    }
    writer.write_all(b"\n").await?;
    for (i, gene) in gene_vec.iter().enumerate() {
        writer.write_all(gene.as_bytes()).await?;
        // writer.write_all(b",").await?;
        for j in matrix[i].iter() {
            writer.write_all(b",").await?;
            writer.write_all(j.to_string().as_bytes()).await?;
        }
        // writer.write_all(b",").await?;
        writer.write_all(b"\n").await?;
    }
    // writer.write_all(b"\n").await?;

    Ok(())
}

pub async fn writer_final_hash(
    final_res_hash: HashMap<String, Vec<Vec<String>>>,
) -> Result<(), Error> {
    // for path_file in final_res_hash.keys() {
    //     let mut writer = File::create(path_file).await?;
    //     let res_tmp_vec = final_res_hash
    //         .get(path_file)
    //         .ok_or_else(|| anyhow!("No data found for the given file"))?;
    //     writer
    //         .write_all(b"gene_a\tgene_b\tpearson_value\tnote\n")
    //         .await?;
    //     for res_tmp_vec in res_tmp_vec {
    //         writer
    //             .write_all(
    //                 format!(
    //                     "{}\t{}\t{}\t{}\n",
    //                     res_tmp_vec[0], res_tmp_vec[1], res_tmp_vec[2], res_tmp_vec[3]
    //                 )
    //                 .as_bytes(),
    //             )
    //             .await?;
    //     }
    //     // 关闭文件
    //     // writer.flush()?;
    //     println!("{} 写入完成", path_file);
    // }
    // 线程加异步
    let handles: Vec<_> = final_res_hash
        .into_iter()
        .map(|(path_file, res_tmp_vec)| {
            tokio::spawn(async move {
                let mut writer = File::create(path_file).await?;
                writer
                    .write_all(b"gene_a\tgene_b\tpearson_value\tnote\n")
                    .await?;

                for res_tmp in res_tmp_vec {
                    writer
                        .write_all(
                            format!(
                                "{}\t{}\t{}\t{}\n",
                                res_tmp[0], res_tmp[1], res_tmp[2], res_tmp[3]
                            )
                            .as_bytes(),
                        )
                        .await?;
                }

                Ok(())
            })
        })
        .collect();

    for handle in handles {
        handle.await??;
    }

    Ok(())
}

pub fn write_2d_vec_to_file(
    data: Vec<Vec<String>>,
    file_path: &str,
    sep: &str,
) -> Result<(), Error> {
    let mut file = std::fs::File::create(file_path)?;

    for row in data {
        let row_str = row.join(sep); // 定义列元素的分隔符
                                     // writeln!(file, "{}", row_str)?;
        file.write_all(row_str.as_bytes())?;
        file.write_all(b"\n")?;
    }
    println!("Data has been written to {}", file_path);
    Ok(())
}
// 根据两个Vec和输出路径生成hist文件名
pub fn generate_hist_file_name(
    samples_vec: &Vec<String>,
    labels_vec: &Vec<String>,
    output_path: &str,
    file_type: &str,
) -> Result<HashMap<String, String>, Error> {
    let mut file_names = HashMap::new();
    // 和samples_tag zip生成文件名
    for lable in labels_vec {
        let tags = lable.split(",").collect::<Vec<&str>>();
        let file_name = samples_vec
            .iter()
            .zip(tags.iter())
            .map(|(sample, tag)| format!("{}_{}", sample, tag))
            .collect::<Vec<String>>()
            .join("_");
        let path_file = format!("{}/{}.{}", output_path, file_name, file_type);
        file_names.insert(lable.to_string(), path_file);
    }
    // println!("{:#?}", file_names);
    Ok(file_names)
}
#[tokio::test]
async fn test_write_2d_vec_to_file() {
    let data = vec![
        vec!["1".to_string(), "2".to_string()],
        vec!["3".to_string(), "4".to_string()],
        vec!["5".to_string(), "6".to_string()],
    ];
    let file_path = "test_data/output.csv";
    write_2d_vec_to_file(data, file_path, &",").unwrap();
}
#[tokio::test]
async fn test_writer_entropy_from_hash() {
    let mut data = HashMap::new();
    let file_path = "test_data/01_entropy_output.csv";
    data.insert("gene1".to_string(), vec![5, 5, 5, 5]);
    data.insert("gene2".to_string(), vec![5, 2, 2, 1]);
    data.insert("gene3".to_string(), vec![5, 5, 5, 9]);
    writer_entropy_from_hash(&data, file_path).await.unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_hist_file_name() {
        let samples_vec = vec![
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
            "d".to_string(),
        ];
        let labels_vec = vec!["1,0,0,0".to_string(), "0,1,0,0".to_string()];
        let output_path = "test_data";
        let result =
            generate_hist_file_name(&samples_vec, &labels_vec, output_path, "svg").unwrap();
        // assert_eq!(result, expected_file_names);
        println!("{:?}", result);
    }
}
