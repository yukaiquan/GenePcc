fn calculate_rank(data: &Vec<f64>) -> Vec<usize> {
    // 创建一个辅助向量，对数据进行排序并生成对应的秩次
    let mut indexed_data: Vec<(usize, &f64)> = data.iter().enumerate().collect();
    indexed_data.sort_by(|a, b| a.1.partial_cmp(b.1).unwrap());

    // 使用辅助向量记录排序后的索引顺序
    let ranks: Vec<usize> = indexed_data.iter().map(|&(i, _)| i).collect();

    // 创建一个保存秩次的向量
    let mut ranks_result = vec![0; data.len()];
    for (idx, &original_index) in ranks.iter().enumerate() {
        ranks_result[original_index] = idx + 1;
    }
    ranks_result
}

fn spearman_rank_correlation(x: &Vec<f64>, y: &Vec<f64>) -> f64 {
    // 计算秩次
    let ranks_x = calculate_rank(x);
    let ranks_y = calculate_rank(y);

    let n = x.len();

    // 计算皮尔逊相关系数
    let mean_rank_x: f64 = ranks_x.iter().sum::<usize>() as f64 / n as f64;
    let mean_rank_y: f64 = ranks_y.iter().sum::<usize>() as f64 / n as f64;

    let numerator: f64 = ranks_x
        .iter()
        .zip(ranks_y.iter())
        .map(|(rank_x, rank_y)| (*rank_x as f64 - mean_rank_x) * (*rank_y as f64 - mean_rank_y))
        .sum();

    let denominator_x: f64 = ranks_x
        .iter()
        .map(|&rank_x| (rank_x as f64 - mean_rank_x) * (rank_x as f64 - mean_rank_x))
        .sum();

    let denominator_y: f64 = ranks_y
        .iter()
        .map(|&rank_y| (rank_y as f64 - mean_rank_y) * (rank_y as f64 - mean_rank_y))
        .sum();

    let denominator = (denominator_x * denominator_y).sqrt();

    let spearman_corr = if denominator != 0.0 {
        numerator / denominator
    } else {
        0.0
    };

    spearman_corr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spearman_correlation() {
        let gene_x = vec![20.0, 22.0, 60.0, 80.0, 90.0];
        let gene_y = vec![18.0, 25.0, 58.0, 34.0, 22.0];

        // 计算斯皮尔曼等级相关系数
        let spearman_corr = spearman_rank_correlation(&gene_x, &gene_y);
        println!("Spearman's rank correlation coefficient: {}", spearman_corr);
        // assert_eq!(pearson_corr, 0.9982302625204867)
    }
}
