use anyhow::{anyhow, Error, Ok, Result};
use std::collections::HashMap;

///根据计算得到的pearson相关系数生成每一个分组的皮尔逊相关系数哈希表\
pub fn generate_pcc_hashmap(
    res_vec: &Vec<Vec<String>>,
) -> Result<HashMap<String, Vec<f64>>, Error> {
    let mut tag_hash: HashMap<String, Vec<f64>> = HashMap::new();
    for string_value in res_vec {
        // 如果为NaN则跳过
        // println!("{:?}", string_value);
        let tag = &string_value[0];
        let value = &string_value[3];
        let parsed_number = value.parse::<f64>()?;
        // group
        if let Some(entry) = tag_hash.get_mut(tag) {
            entry.push(parsed_number);
            // entry.push(0.0);
        } else {
            tag_hash.insert(tag.to_string(), vec![parsed_number]);
        }
    }
    Ok(tag_hash)
}

///计算数组的左右各5%
pub fn calculate_percentiles_for_pcc(pcc_values: &[f64]) -> Result<Vec<f64>, Error> {
    if pcc_values.is_empty() {
        return Ok(vec![0.0, 0.0]);
    }

    let mut sorted_pcc_values = pcc_values.to_vec();
    sorted_pcc_values.sort_by(|a, b| a.partial_cmp(b).expect("NaN in PCC values"));

    let total_count = sorted_pcc_values.len();
    if total_count < 2 {
        return Ok(vec![0.0, 0.0]);
    }

    let mut per5_num = (total_count as f64 * 0.05).round() as usize;
    if per5_num == 0 {
        per5_num = 1;
    }
    // println!("per5_num: {}", per5_num);
    let pre_per5 = sorted_pcc_values[per5_num - 1];
    let last_per5 = sorted_pcc_values[total_count - per5_num];
    // Ok((pre_per5, last_per5))
    Ok(vec![pre_per5, last_per5])
}
// 根据hashmap去除low entropy的值
// rm_hash中存储的是需要包括的值
pub fn remove_low_entropy_values(
    rm_hash: &HashMap<String, bool>,
    expressions: &Vec<Vec<f64>>,
    gene_vec: &Vec<String>,
) -> Result<(Vec<Vec<f64>>, Vec<String>), Error> {
    let mut new_expressions = Vec::new();
    let mut new_gene_vec = Vec::new();
    for i in 0..gene_vec.len() {
        let gene_name = &gene_vec[i];
        if rm_hash.contains_key(gene_name) {
            new_expressions.push(expressions[i].clone());
            new_gene_vec.push(gene_name.clone());
        } else {
            continue;
        }
    }
    Ok((new_expressions, new_gene_vec))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_cacu_per() {
        // 假设有一组数据
        let data_f64 = vec![1.2, 3.4, 5.6, 7.8, 9.1, 10.11, 11.12];

        // 计算最小值和最大值的5%的临界值（f64类型）
        let m5_percent_value = calculate_percentiles_for_pcc(&data_f64).unwrap();

        println!(
            "f64 - Value corresponding to the 5th percentile: {:?}",
            m5_percent_value
        );
        println!(
            "f64 - Value corresponding to the 95th percentile: {:?}",
            m5_percent_value
        );
    }
}
