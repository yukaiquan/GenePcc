use anyhow::{anyhow, Error, Result};
// use bgzip::{BGZFError, BGZFReader};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
// use std::fmt::Result;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Lines, Write};
use std::path::Path;

// //使用Box包装使下面函数可以读取压缩文件和非压缩文件 这个读取有问题 部分压缩文件只能读取一半
pub fn read_gz_file(path: &Path) -> Result<Lines<BufReader<Box<dyn std::io::Read>>>, Error> {
    let file = File::open(path)?;
    let reader: Box<dyn std::io::Read> = if path.extension().unwrap() == "gz" {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::new(reader).lines())
}
// //使用Box包装使下面函数可以读取压缩文件和非压缩文件
// pub fn read_bgz_file(path: &Path) -> Result<Lines<BufReader<Box<dyn std::io::Read>>>, Error> {
//     let file = File::open(path)?;
//     let reader: Box<dyn std::io::Read> = if path.extension().unwrap() == "gz" {
//         Box::new(BGZFReader::new(file)?)
//     } else {
//         Box::new(file)
//     };
//     Ok(BufReader::new(reader).lines())
//     // Ok(())
// }

// 读取分组文件
// 从计算机层面我们并不需要知道详细分组是什么，我们只需要一个分辨的符号即可
// 后续根据分组进行聚类
pub fn read_group_file(path: &Path) -> Result<(HashMap<String, String>, Vec<String>), Error> {
    let mut group_map = HashMap::new();
    let mut sampls_tag = Vec::new();
    let mut lines_num = 0;
    if let Ok(lines) = read_gz_file(path) {
        for line in lines {
            if let Ok(line) = line {
                // groups.push(group);
                // 如果是第一行
                if lines_num == 0 {
                    lines_num += 1;
                    sampls_tag = line
                        .split(",")
                        .map(|s| s.to_string())
                        // 跳过第一个
                        .skip(1)
                        .collect::<Vec<String>>();
                    continue;
                } else {
                    // let samples_vec: Vec<&str> = line.split(",").collect();
                    // 获得所有权
                    let samples_vec: Vec<String> = line.split(",").map(|s| s.to_string()).collect();
                    group_map.insert(
                        samples_vec[0].to_string(),
                        samples_vec[1..].join(",").to_string(),
                    );
                    lines_num += 1;
                }
            }
        }
    }
    println!(
        "read group file success, samples number is : {}",
        lines_num - 1
    );
    Ok((group_map, sampls_tag))
}

//读取去除批次效应的表达矩阵
pub fn read_expression_matrix(
    path: &Path,
) -> Result<(Vec<Vec<f64>>, Vec<String>, Vec<String>), Error> {
    let mut expression_vec = Vec::new();
    let mut samples = Vec::new();
    let mut gene_vec = Vec::new();
    let mut lines_num = 0;
    // let buf_reader = BufReader::new(reader);
    if let Ok(lines) = read_gz_file(path) {
        for line in lines {
            if let Ok(line) = line {
                // println!("{}", line);
                // 这里可以对每一行进行处理，而不是将其全部存储在内存中
                // 如果是第一行
                if lines_num == 0 {
                    let samples_vec: Vec<&str> = line.split(",").collect();
                    for sample in samples_vec.iter().skip(1) {
                        samples.push(sample.to_string());
                    }
                } else {
                    let samples_vec: Vec<String> =
                        line.split(",").map(|s| s.parse().unwrap()).collect();
                    gene_vec.push(samples_vec[0].to_string());
                    // expression_map.insert(samples_vec[0].to_string(), samples_vec[1..].to_vec());
                    // println!("{:#?}", &samples_vec[1..]);
                    expression_vec.push(
                        samples_vec[1..]
                            .to_vec()
                            .into_iter()
                            .map(|x| x.parse::<f64>().unwrap())
                            .collect(),
                    );
                }
                lines_num += 1;
            } else {
                // return Err("Error reading line".into());
                return Err(anyhow!("Error reading line"));
            }
        }
    }

    println!(
        "read expression matrix success, lines_num: {}",
        lines_num - 1
    );
    Ok((expression_vec, samples, gene_vec))
}
pub fn read_tpm_matrix(
    path: &Path,
    output: &str,
) -> Result<(Vec<Vec<f64>>, Vec<String>, Vec<String>), Error> {
    let mut expression_vec = Vec::new();
    let mut samples = Vec::new();
    let mut gene_vec = Vec::new();

    if Path::new(&format!("{}/03_gene_matrix.csv", output)).exists() {
        let result = read_expression_matrix(Path::new(&format!("{}/03_gene_matrix.csv", output)))?;
        expression_vec = result.0;
        samples = result.1;
        gene_vec = result.2;
    } else {
        let result = read_expression_matrix(path)?;
        expression_vec = result.0;
        samples = result.1;
        gene_vec = result.2;
    }

    Ok((expression_vec, samples, gene_vec))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_read_group_file() {
        let path = Path::new("test_data/01_tissue_sample.csv.gz");
        let (group_map, samples_tag) = read_group_file(path).unwrap();
        assert_eq!(group_map.len(), 258);
        // println!("{:?}", group_map);
        assert_eq!(samples_tag.len(), 9);
        println!("{:?}", samples_tag);
    }
    #[test]
    fn test_read_expression_file() {
        let path = Path::new("test_data/no_batch_tpm_matrix.csv.gz");
        let (expression_map, samples, gene_vec) = read_expression_matrix(path).unwrap();
        assert_eq!(expression_map.len(), 120352);
        assert_eq!(samples.len(), 258);
        // println!("{:?}", expression_map);
        // println!("{:?}", samples.len());
    }
    #[test]
    fn test_read_expression2_file() {
        let path = Path::new("test_data/test.csv.gz");
        let (expression_map, samples, gene_vec) = read_expression_matrix(path).unwrap();
        assert_eq!(expression_map.len(), 999);
        assert_eq!(samples.len(), 258);
        // println!("{:?}", expression_map);
        // println!("{:?}", samples.len());
    }
}
