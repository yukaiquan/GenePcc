// Butte, A. J., Tamayo, P., Slonim, D., Golub, T. R., & Kohane, I. S. (2000). Discovering functional relationships between RNA expression and chemotherapeutic susceptibility using relevance networks. Proceedings of the National Academy of Sciences of the United States of America, 97(22), 12182–12186. https://doi.org/10.1073/pnas.220392197
// Butte, A. J., & Kohane, I. S. (2000). Mutual information relevance networks: functional genomic clustering using pairwise entropy measurements. Pacific Symposium on Biocomputing. Pacific Symposium on Biocomputing, 418–429. https://doi.org/10.1142/9789814447331_0040

// where log2 is base 2 logarithm, and p(x) is the probability a value was within decile x of that feature. For example, a gene with the expression amounts 20, 22, 60, 80, and 90 would have deciles 7 units wide, with two values in the first decile, one in the sixth decile, and one in the ninth and 10th decile, making H = 1.92
use anyhow::{anyhow, Error, Ok};
use dashmap::DashMap;
use ndarray::{s, Array, Array2, Axis};
use rayon::prelude::*;
use std::collections::{btree_map::Keys, HashMap};
use std::sync::{Arc, Mutex};

fn calculate_entropy(gene_expression: &[f64], n: usize) -> f64 {
    // 计算分位数的宽度
    let min_value = *gene_expression
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let max_value = *gene_expression
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let decile_width = (max_value - min_value) / n as f64; // 计算分位数的宽度

    // 计算每个值在分位数中的位置
    let deciles_index: Vec<usize> = gene_expression
        .iter()
        .map(|&value| ((value - min_value) / decile_width).floor() as usize)
        .collect();

    // 计算每个表达量值的概率 p(x)
    let mut deciles_counts: HashMap<usize, usize> = HashMap::new();
    for &index in &deciles_index {
        *deciles_counts.entry(index).or_insert(0) += 1;
    }

    let total_values = gene_expression.len() as f64;
    let probabilities: HashMap<usize, f64> = deciles_counts
        .iter()
        .map(|(&index, &count)| (index, count as f64 / total_values))
        .collect();

    // 计算H值
    let h_value: f64 = probabilities
        .values()
        .filter_map(|&p_x| {
            if p_x > 0.0 {
                Some(p_x * p_x.log2())
            } else {
                None
            }
        })
        .map(|x| -x)
        .sum();

    h_value
}

use std::cmp::Ordering;
#[derive(Debug, PartialEq, Clone)]
pub enum MixedNumber {
    String(String),
    Float(f64),
}

impl PartialOrd for MixedNumber {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (Self::Float(a), Self::Float(b)) => a.partial_cmp(b),
            (Self::String(a), Self::String(b)) => a.partial_cmp(b),
            _ => None, // Return None for mismatched types
        }
    }
}

pub fn core_threads_rm_low(
    group_uniq_map: &HashMap<String, Vec<usize>>,
    expression_vec: &Vec<Vec<f64>>,
    gene_vec: &Vec<String>,
    // samples: &Vec<String>,
) -> Result<HashMap<String, Vec<Vec<MixedNumber>>>, Error> {
    let mut rm_low_result: DashMap<String, Vec<Vec<MixedNumber>>> = DashMap::new();
    // let mut res_h_value = Vec::new();
    // for key in group_uniq_map.keys() {
    //     // for (key,value) in group_uniq_map.iter(){
    //     let group_index = group_uniq_map
    //         .get(key)
    //         .ok_or_else(|| anyhow!("group not found"))?;
    //     for (index, exp) in expression_vec.iter().enumerate() {
    //         let target_exp_vec = exp
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(index, &value)| {
    //                 if group_index.contains(&index) {
    //                     Some(value)
    //                 } else {
    //                     None
    //                 }
    //             })
    //             .collect::<Vec<f64>>();
    //         // println!("{:#?}", target_exp_vec);
    //         let h_value = calculate_entropy(&target_exp_vec, 10);
    //         if let Some(tmp_vec) = rm_low_result.get_mut(key) {
    //             tmp_vec.push(vec![
    //                 MixedNumber::String(gene_vec[index].to_string()),
    //                 MixedNumber::Float(h_value),
    //             ]);
    //         } else {
    //             let entry = rm_low_result.entry(key.clone()).or_insert(Vec::new());
    //             entry.push(vec![
    //                 MixedNumber::String(gene_vec[index].to_string()),
    //                 MixedNumber::Float(h_value),
    //             ]);
    //         }
    //     }
    // }
    // 不需要分组

    // for (index, exp) in expression_vec.iter().enumerate() {
    //     let target_exp_vec = exp;
    //     // println!("{:#?}", target_exp_vec);
    //     let h_value = calculate_entropy(&target_exp_vec, 10);
    //     if let Some(tmp_vec) = rm_low_result.get_mut("all") {
    //         tmp_vec.push(vec![
    //             MixedNumber::String(gene_vec[index].to_string()),
    //             MixedNumber::Float(h_value),
    //         ]);
    //     } else {
    //         let entry = rm_low_result.entry("all".to_string()).or_insert(Vec::new());
    //         entry.push(vec![
    //             MixedNumber::String(gene_vec[index].to_string()),
    //             MixedNumber::Float(h_value),
    //         ]);
    //     }
    // }
    // 改为并行
    expression_vec
        .par_iter()
        .enumerate()
        .for_each(|(index, target_exp_vec)| {
            let h_value = calculate_entropy(target_exp_vec, 10);
            if let Some(mut tmp_vec) = rm_low_result.get_mut("all") {
                tmp_vec.push(vec![
                    MixedNumber::String(gene_vec[index].to_string()),
                    MixedNumber::Float(h_value),
                ]);
            } else {
                let mut entry = rm_low_result.entry("all".to_string()).or_insert(Vec::new());
                entry.push(vec![
                    MixedNumber::String(gene_vec[index].to_string()),
                    MixedNumber::Float(h_value),
                ]);
            }
        });

    // rm_low_result根据key不同对value排序并记录值
    // println!("{:#?}", rm_low_result);
    let mut sort_rm_low_result = HashMap::new();
    for pair in rm_low_result.iter() {
        let key = pair.key();
        let value_vec = rm_low_result
            .get(key)
            .ok_or_else(|| anyhow!("key is error!"))?;
        // Extract and sort the Float variants
        let mut mut_value = value_vec.clone();
        mut_value.sort_by(|a, b| {
            let a_floats: Vec<f64> = a
                .iter()
                .filter_map(|num| {
                    if let MixedNumber::Float(val) = num {
                        Some(*val)
                    } else {
                        None
                    }
                })
                .collect();

            let b_floats: Vec<f64> = b
                .iter()
                .filter_map(|num| {
                    if let MixedNumber::Float(val) = num {
                        Some(*val)
                    } else {
                        None
                    }
                })
                .collect();

            a_floats.partial_cmp(&b_floats).unwrap()
        });
        // 更新rm_low_result
        sort_rm_low_result.insert(key.to_string(), mut_value);
    }
    // println!("{:#?}", sort_rm_low_result);
    Ok(sort_rm_low_result)
}

// 根据排序的顺序得到该基因在此分组中的位置
// pub fn cacu_index_by_sort(
//     hash: &HashMap<String, Vec<Vec<MixedNumber>>>,
// ) -> Result<HashMap<String, Vec<usize>>, Error> {
//     let mut res_hash = HashMap::new();
//     for key in hash.keys() {
//         let values = hash.get(key).ok_or_else(|| anyhow!("key is error!"))?;
//         for (index, value) in values.iter().enumerate() {
//             for num in value {
//                 let mut gene_name = "";
//                 let mut entropy_value = 0.0;
//                 match num {
//                     MixedNumber::Float(val) => entropy_value = *val,
//                     MixedNumber::String(val) => gene_name = val,
//                     // Add handling for other MixedNumber variants if needed
//                 }
//                 if res_hash.contains_key(gene_name) {
//                     let mut res_tmp = Vec::new();
//                     let tmp: &mut Vec<usize> = res_hash
//                         .get_mut(gene_name)
//                         .ok_or_else(|| anyhow!("gene name is error!"))?;
//                     res_tmp.extend(tmp.clone());
//                     res_tmp.push(index);
//                     res_hash.insert(gene_name.to_string(), res_tmp);
//                 } else {
//                     res_hash.insert(gene_name.to_string(), vec![index]);
//                 }
//             }
//         }
//     }
//     Ok(res_hash)
// }
pub fn cacu_index_by_sort(
    hash: &HashMap<String, Vec<Vec<MixedNumber>>>,
) -> Result<HashMap<String, Vec<usize>>, Error> {
    // let res_hash: Arc<Mutex<HashMap<String, Vec<usize>>>> = Arc::new(Mutex::new(HashMap::new()));
    let res_hash: DashMap<String, Vec<usize>> = DashMap::new();

    hash.par_iter().for_each(|(key, values)| {
        values.iter().enumerate().for_each(|(index, value)| {
            value.iter().for_each(|num| {
                let mut gene_name = "";
                let mut entropy_value = 0.0;

                match num {
                    MixedNumber::Float(val) => entropy_value = *val,
                    MixedNumber::String(val) => gene_name = val,
                    // Add handling for other MixedNumber variants if needed
                }

                // let mut res_hash_locked = res_hash.lock().unwrap();
                // if res_hash_locked.contains_key(gene_name) {
                //     let mut res_tmp = Vec::new();
                //     if let Some(tmp) = res_hash_locked.get_mut(gene_name) {
                //         res_tmp.extend(tmp.clone());
                //         res_tmp.push(index);
                //         res_hash_locked.insert(gene_name.to_string(), res_tmp);
                //     }
                // } else {
                //     res_hash_locked.insert(gene_name.to_string(), vec![index]);
                // }
                if let Some(mut tmp) = res_hash.get_mut(gene_name) {
                    tmp.push(index);
                } else {
                    res_hash.insert(gene_name.to_string(), vec![index]);
                }
            });
        });
    });

    // // Extract the result from the Arc<Mutex<HashMap>>
    // let result = Arc::try_unwrap(res_hash)
    //     .unwrap()
    //     .into_inner()
    //     .map_err(|e| anyhow!("Failed to unwrap Arc: {}", e))?;
    let result: HashMap<_, _> = res_hash
        .iter()
        .map(|pair| (pair.key().clone(), pair.value().clone()))
        .collect();
    Ok(result)
}
// 根据基因vec和res_hash确定保留的基因
pub fn cacu_by_entropy_gene(
    gene: &Vec<String>,
    res_hash: &HashMap<String, Vec<usize>>,
) -> Result<HashMap<String, bool>, Error> {
    if gene.len() <= 10 {
        return Ok(res_hash
            .iter()
            // 跳过空值
            .filter(|(k, v)| !k.is_empty())
            .map(|(k, v)| (k.to_string(), true))
            .collect());
    }
    let mut res_gene_hash = HashMap::new();
    // 计算基因数量的5%用于去除
    let rm_gene_num = (gene.len() as f64 * 0.05).round() as usize;
    for (key, value) in res_hash.iter() {
        if key.is_empty() {
            continue;
        }
        if value.iter().all(|&x| x < rm_gene_num) {
            println!("Skipping gene: {}", key);
            continue; // Skip to the next iteration
        } else {
            res_gene_hash.insert(key.to_string(), true);
        }
    }
    let vaild_num = res_gene_hash.len();
    println!(
        "it have {:#?} genes entropy is very lower",
        gene.len() - vaild_num
    );

    Ok(res_gene_hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rm_low() {
        // Gene expression data
        let gene_expression = vec![20.0, 22.0, 60.0, 80.0, 90.0];

        // Calculate entropy for gene expression data
        let entropy = calculate_entropy(&gene_expression, 10);
        println!("Entropy for gene expression data: {:}", entropy);
    }
    #[test]
    fn test_core_threads_rm_low() {
        // Create test data for inputs
        let mut group_uniq_map: HashMap<String, Vec<usize>> = HashMap::new();
        group_uniq_map.insert("key1".to_string(), vec![0, 1, 2]);
        group_uniq_map.insert("key2".to_string(), vec![3, 4, 5]);

        let expression_vec: Vec<Vec<f64>> = vec![
            vec![2.0, 2.0, 3.0, 2.0, 3.0, 2.0],
            vec![6.0, 5.0, 6.0, 8.0, 8.0, 16.0],
            vec![5.0, 6.0, 6.0, 12.0, 10.0, 2.0],
        ];

        let gene_vec: Vec<String> = vec![
            "gene1".to_string(),
            "gene2".to_string(),
            "gene3".to_string(),
        ];
        // Call the function and assert the expected result
        let result = core_threads_rm_low(&group_uniq_map, &expression_vec, &gene_vec).unwrap();
        println!("{:#?}", result);
        let res_hash = cacu_index_by_sort(&result).unwrap();
        let res = cacu_by_entropy_gene(&gene_vec, &res_hash);
        println!("{:#?}", res);
        // assert!(result.is_ok()); // Check if the function returns Ok(())
    }
}
