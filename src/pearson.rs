use crate::caculate_core::{
    calculate_percentiles_for_pcc, generate_pcc_hashmap, remove_low_entropy_values,
};
use crate::writerfile::{
    generate_hist_file_name, writer_entropy, writer_entropy_from_hash, writer_filter_matrix,
    writer_final_hash,
};
use ansi_term::Colour;
use anyhow::{anyhow, Error, Ok, Result};
use core::num;
use ndarray::{array, stack, Array2, ArrayView1};
use ndarray::{s, Array, Axis};
use ndarray_stats::CorrelationExt;
use num_cpus;
use rayon::prelude::*;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::{Continuous, StudentsT};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{mpsc, Mutex};
use std::sync::{Arc, RwLock};
use std::thread;

fn calculate_correlation<'a>(
    a_gene: ArrayView1<'a, f64>,
    b_gene: ArrayView1<'a, f64>,
    use_nd: bool,
) -> f64 {
    if use_nd {
        let gene = stack(Axis(0), &[a_gene, b_gene]).unwrap();
        if let Some(tmp) = gene.pearson_correlation().ok() {
            return if tmp[(0, 1)].is_nan() {
                0.0
            } else {
                tmp[(0, 1)]
            };
        } else {
            return 0.0;
        }
    } else {
        let a_gene_vec: Vec<f64> = a_gene.to_vec();
        let b_gene_vec: Vec<f64> = b_gene.to_vec();
        return pearson_correlation_vec(&a_gene_vec, &b_gene_vec);
    }
}

pub fn core_threads_cacu_pearson(
    group_uniq_map: &HashMap<String, Vec<usize>>,
    expression_vec: Vec<Vec<f64>>,
    gene_vec: &Vec<String>,
    samples: &Vec<String>,
    samples_tag: &Vec<String>,
    use_nd: bool,
    output: &str,
    pearson_threshold: f64,
) -> Result<Vec<String>, Error> {
    let flat_vec: Vec<f64> = expression_vec.into_iter().flatten().collect(); // 扁平化二维数组
                                                                             // let expression_vec_nd_array = Array::from_shape_vec((gene_vec.len(), samples.len()), flat_vec)?;
    let expression_vec_nd_array = Array::from_shape_vec((gene_vec.len(), samples.len()), flat_vec)?;

    let mut res: Vec<String> = Vec::new();
    let cpu_nums = num_cpus::get();
    // 1/8用于写出，剩下的用于并行计算
    let writer_nums = 1;
    let rayon_nums = cpu_nums - writer_nums;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(rayon_nums)
        .build()?;

    for key in group_uniq_map.keys() {
        if key.to_uppercase() == "CK" {
            continue;
        }
        if let Some(group_index) = group_uniq_map.get(key) {
            let file_name: String = key
                .split(',')
                .zip(samples_tag.iter())
                .map(|(s, tag)| format!("{}_{}", tag, s))
                .collect::<Vec<String>>()
                .join("_");
            // 判断文件是否存在
            let file_path = format!("{}/{}_rawPCC.txt", output, file_name);
            let file_gz = format!("{}/{}_rawPCC.txt.gz", output, file_name);
            if Path::new(&file_path).exists() | Path::new(&file_gz).exists() {
                println!("{} already exists!", file_path);
                res.push(format!("{}/{}_rawPCC.txt", output, file_name));
                continue;
            }
            res.push(format!("{}/{}_rawPCC.txt", output, file_name));
            let expression_vec_nd_array_ref =
                Arc::new(expression_vec_nd_array.select(Axis(1), group_index));
            // println!("{:#?}", expression_vec_nd_array_ref.shape());
            let num_rows = expression_vec_nd_array_ref.shape()[0];
            let num_col = expression_vec_nd_array_ref.shape()[1];
            // 多生产者dan消费者模型
            let (sender, receiver): (mpsc::Sender<String>, mpsc::Receiver<String>) =
                mpsc::channel();
            let sender_arc_c = Arc::new(sender);
            let gene_vec_c = Arc::new(gene_vec.clone());
            // let receiver = Arc::new(Mutex::new(receiver));
            // // 创建文件写入锁
            // let file_mutex = Arc::new(Mutex::new(
            //     File::create(format!("{}/{}_rawPCC.txt", output, file_name))
            //         .expect("Failed to create file"),
            // ));
            // // 创建计数器用于跟踪活动的消费者线程数量
            // let active_consumers = Arc::new(Mutex::new(writer_nums)); // n个消费者线程
            // let mut handles = Vec::new();
            // for _ in 0..writer_nums {
            //     let receiver_clone = Arc::clone(&receiver);
            //     let active_consumers_clone = Arc::clone(&active_consumers);
            //     let file_mutex_clone = Arc::clone(&file_mutex);
            //     let handle = thread::spawn(move || {
            //         // 消费者接收数据并写入文件
            //         loop {
            //             match receiver_clone.lock().unwrap().recv() {
            //                 std::result::Result::Ok(data) => {
            //                     write_data(data, &file_mutex_clone);
            //                 }
            //                 Err(_) => {
            //                     // 通道已关闭，减少活动消费者线程计数
            //                     let mut count = active_consumers_clone.lock().unwrap();
            //                     *count -= 1;

            //                     // 如果所有消费者线程都完成，退出循环
            //                     if *count == 0 {
            //                         break;
            //                     }
            //                 }
            //             }
            //         }
            //     });

            //     // 将消费者线程的 JoinHandle 存储到 handles 中
            //     handles.push(handle);
            // }

            let consumer_thread =
                spawn_consumer(receiver, format!("{}/{}_rawPCC.txt", output, file_name));
            // let sender_arc_c_c = Arc::clone(&sender_arc_c);
            pool.install(|| {
                (0..num_rows)
                    .into_par_iter()
                    .for_each_with(sender_arc_c, |sender_arc_c, i| {
                        // let expression_vec_nd_array_ref_ref = Arc::clone(&expression_vec_nd_array_ref);
                        let sender_arc_c_c = Arc::clone(&sender_arc_c);
                        (i + 1..num_rows).into_par_iter().for_each_with(
                            sender_arc_c_c,
                            |sender_arc_c_c, j| {
                                let sender_arc_c_c_c = Arc::clone(&sender_arc_c_c);
                                let gene1_arr = expression_vec_nd_array_ref.row(i);
                                let correlation_value = calculate_correlation(
                                    gene1_arr,
                                    expression_vec_nd_array_ref.row(j),
                                    use_nd,
                                );
                                if correlation_value.abs() <= pearson_threshold {
                                    return;
                                }
                                // println!("{:#?}", gene1_arr.shape());
                                // 计算 t 统计量
                                // let len_gene1 = gene1_arr.shape()[0];
                                let t_stat = t_statistic(correlation_value, num_col);
                                let p_value = p_value(t_stat, num_col);
                                if p_value > 0.01 {
                                    return;
                                }
                                let gene_name1 = &gene_vec_c[i];
                                let gene_name2 = &gene_vec_c[j];
                                let note;
                                if correlation_value <= 0.0 {
                                    note = "negative";
                                } else if correlation_value >= 0.0 {
                                    note = "positive";
                                } else {
                                    note = "neutral";
                                }
                                // vec![
                                //     key.clone(),
                                //     gene_name1.to_string(),
                                //     gene_name2.to_string(),
                                //     correlation_value.to_string(),
                                // ]
                                // println!("{}\t{}\t{}\n", gene_name1, gene_name2, correlation_value);
                                sender_arc_c_c_c
                                    .send(format!(
                                        "{}\t{}\t{}\t{}\t{}\n",
                                        gene_name1, gene_name2, correlation_value, p_value, note
                                    ))
                                    .expect("Error sending data");
                            },
                        )
                    });
            });

            //
            // spawn_producer.join().expect("Error joining thread");
            // let _ = receiver.iter().for_each(|s| {
            //     let mut file = File::create(&file_name).expect("Failed to create file");
            //     // 从通道接收数据并写入文件
            //     // if let Err(err) = file.write_all(data.as_bytes()) {
            //     //     eprintln!("Error writing to file: {}", err);
            //     // }
            //     let datas: Vec<&str> = s.trim().split("\t").collect();
            //     // println!("{:#?}", datas);
            //     if datas[2].parse::<f64>().unwrap() == 0.0 {
            //         return;
            //     } else {
            //         file.write_all(s.as_bytes()).expect("Error writing to file");
            //     }
            // });
            println!("{}/{}_rawPCC.txt", output, file_name);

            consumer_thread.join().expect("Error joining thread");

            // // 等待所有消费者线程结束
            // for handle in handles {
            //     handle.join().expect("Error joining thread");
            // }
        }
    }
    Ok(res)
}
fn write_data(data: String, file_mutex: &Arc<Mutex<File>>) {
    // 处理数据并写入文件
    let datas: Vec<&str> = data.trim().split("\t").collect();
    let correlation_value = datas[2].parse::<f64>().unwrap();

    if correlation_value == 0.0 {
        return;
    }

    let p_value = datas[3].parse::<f64>().unwrap();
    let mut note = "";

    if correlation_value <= 0.0 {
        note = "negative";
    } else if correlation_value >= 0.0 {
        note = "positive";
    } else {
        note = "neutral";
        return;
    }
    let mut file = file_mutex.lock().unwrap();
    file.write_all(
        format!(
            "{}\t{}\t{}\t{}\t{}\n",
            datas[0], datas[1], correlation_value, p_value, note
        )
        .as_bytes(),
    )
    .expect("Error writing to file");
}
fn spawn_consumer(receiver: mpsc::Receiver<String>, file_name: String) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        // 创建文件并将数据写入文件

        let mut file = BufWriter::new(File::create(file_name).expect("Failed to create file"));

        // 从通道接收数据并写入文件
        while let std::result::Result::Ok(data) = receiver.recv() {
            // if let Err(err) = file.write_all(data.as_bytes()) {
            //     eprintln!("Error writing to file: {}", err);
            // }
            file.write_all(format!("{}", data).as_bytes())
                .expect("Error writing to file");
            // let datas: Vec<&str> = data.trim().split("\t").collect();
            // let correlation_value = datas[2].parse::<f64>().unwrap();
            // // println!("{:#?}", datas);
            // if correlation_value == 0.0 {
            //     continue;
            // } else {
            //     let p_value = datas[3].parse::<f64>().unwrap();
            //     let mut note = "";
            //     if correlation_value <= 0.0 {
            //         note = "negative";
            //     } else if correlation_value >= 0.0 {
            //         note = "positive";
            //     } else {
            //         note = "neutral";
            //         continue;
            //     }

            //     // let p_value = 0.0;
            //     // file.write_all(
            //     //     format!(
            //     //         "{}\t{}\t{}\t{}\t{}\n",
            //     //         datas[0], datas[1], correlation_value, p_value, note
            //     //     )
            //     //     .as_bytes(),
            //     // )
            //     // .expect("Error writing to file");
            //     // file.write_all(data.as_bytes())
            //     //     .expect("Error writing to file");
            // }
        }
        file.flush().unwrap(); // 确保所有缓冲的数据都被写入文件
    })
}

// fn pearson_correlation_matrix(data1: &Array2<f64>, data2: &Array2<f64>) -> Result<f64, Error> {
//     // let num_vars = data.shape()[0];
//     // let mut correlation_matrix = Array2::zeros((num_vars, num_vars));

//     // for i in 0..num_vars {
//     //     for j in i..num_vars {
//     //         // let corr = data
//     //         //     .slice(s![i, ..])
//     //         //     .correlation(&data.slice(s![j, ..]))
//     //         //     .unwrap();
//     //         let tmp_data = data.select(Axis(0), &[i, j]);
//     //         let corr = tmp_data.pearson_correlation().unwrap();
//     //         correlation_matrix[[i, j]] = corr;
//     //         correlation_matrix[[j, i]] = corr;
//     //     }
//     // }
//     // 合并两个data
//     let corr = 0.0;
//     // let tmp_data = ndarray::stack(Axis(0), &[data1.view(), data2.view()])?;
//     // 合并为Array2
//     let tmp_data = ndarray::stack(Axis(0), &[data1.view(), data2.view()]).unwrap();
//     corr = tmp_data.pearson_correlation()?;
//     Ok(corr)
// }

// 计算平均值
fn mean(data: &[f64]) -> f64 {
    data.iter().sum::<f64>() / data.len() as f64
}

// 计算协方差
fn covariance(data1: &[f64], data2: &[f64]) -> f64 {
    let mean1 = mean(data1);
    let mean2 = mean(data2);

    let mut cov = 0.0;
    for i in 0..data1.len() {
        cov += (data1[i] - mean1) * (data2[i] - mean2);
    }
    if data1.len() as f64 == 0.0 {
        return 0.0;
    }
    cov / (data1.len() as f64)
}

// 计算标准差
fn standard_deviation(data: &[f64]) -> f64 {
    let m = mean(data);
    let variance = data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / data.len() as f64;
    variance.sqrt()
}

// 计算皮尔逊相关系数
pub fn pearson_correlation_vec(data1: &[f64], data2: &[f64]) -> f64 {
    let cov = covariance(data1, data2);
    let sd1 = standard_deviation(data1);
    let sd2 = standard_deviation(data2);
    if sd1 == 0.0 || sd2 == 0.0 {
        return 0.0;
    }
    cov / (sd1 * sd2)
}

// 计算p_value
pub fn t_statistic(correlation_coefficient: f64, sample_size: usize) -> f64 {
    let coe_powi = correlation_coefficient.powi(2);
    if coe_powi == 1.0 {
        return 0.0;
    } else {
        correlation_coefficient * ((sample_size as f64 - 2.0) / (1.0 - coe_powi)).sqrt()
    }
}

pub fn p_value(t_stat: f64, sample_size: usize) -> f64 {
    let df = sample_size - 2;
    let t_abs = t_stat.abs();

    let t_distribution =
        StudentsT::new(0.0, 1.0, df as f64).expect("Failed to create t distribution");

    // 使用 clamp 方法确保 t 统计量在合理范围内
    if t_abs.is_nan() {
        return 1.0;
    }

    // println!("{}", t_abs);
    // 计算 p-value
    2.0 * (1.0 - t_distribution.cdf(t_abs))
    // 对密度函数报错进行处理
    // if let Err(p_value) = t_distribution.cdf(t_abs_clamped) {
    //     return 2.0 * (1.0 - p_value);
    // } else {
    //     return 0.0;
    // }
}

#[cfg(test)]
mod tests {
    use super::*;
    // 计算运行时
    use ndarray::Array2;
    use ndarray_stats::CorrelationExt;
    // use statrs::statistics::{bivariate::pearson, Testing};
    use ndarray::Axis;
    use ndarray::Zip;
    use std::time::Instant;
    #[test]
    fn test_pearson_correlation() {
        let gene_x = vec![20.0, 22.0, 60.0, 80.0, 90.0];
        // let gene_x = Array2::from_shape_vec((5, 1), vec![20.0, 22.0, 60.0, 80.0, 90.0]).unwrap();
        let gene_y = vec![18.0, 25.0, 58.0, 78.0, 88.0];
        // let gene_y = Array2::from_shape_vec((5, 1), vec![18.0, 25.0, 58.0, 78.0, 88.0]).unwrap();
        // let gene_x = gene_x.row(0);
        // let gene_y = gene_y.row(0);
        println!("Gene X: {:?}", gene_x);
        println!("Gene Y: {:?}", gene_y);

        let pearson_corr = pearson_correlation_vec(&gene_x, &gene_y);
        println!("Gene X 和 Gene Y 之间的皮尔逊相关系数为：{}", pearson_corr);
        assert_eq!(pearson_corr, 0.9982302625204867)
    }
    #[test]
    fn test_nd_correlation() {
        let data = Array2::from_shape_vec(
            (2, 5),
            vec![20.0, 22.0, 60.0, 80.0, 90.0, 18.0, 25.0, 58.0, 78.0, 88.0],
        )
        .unwrap();

        // let correlation_matrix = pearson_correlation_matrix(&data);
        // println!("Pearson Correlation Matrix:");
        // println!("{:?}", correlation_matrix);
        // 创建两个示例数组
        // let arr1 = Array::from(vec![1, 2, 3]).into_shape((3, 1)).unwrap();
        // let arr2 = Array::from(vec![4, 5, 6]).into_shape((3, 1)).unwrap();

        // // 沿着指定轴合并数组
        // let merged = ndarray::stack(Axis(0), &[arr1.view(), arr2.view()]).unwrap();

        // println!("Merged array:");
        // println!("{:?}", merged);
        // let a = arr2(&[[1., 3., 5.], [2., 4., 6.]]);
        // let corr = a.pearson_correlation().unwrap();
        // let a = arr2(&[[1., 3., 5.], [2., 4., 6.]]);
        let covariance = data.pearson_correlation().unwrap();
        println!("Covariance: {:?}", covariance);
        // 创建两个示例的 Array2
        let arr1 = Array2::from_shape_vec((2, 3), vec![1, 2, 3, 4, 5, 6]).unwrap();
        let arr2 = Array2::from_shape_vec((2, 3), vec![7, 8, 9, 10, 11, 12]).unwrap();

        // 沿着指定轴（例如轴0）堆叠两个数组
        let merged = ndarray::stack(Axis(0), &[arr1.view(), arr2.view()]).unwrap();

        println!("Shape of merged array: {:?}", merged.shape());
    }
    #[test]
    fn test_core_threads_cacu_pearson() {
        let start = Instant::now();
        let mut group_uniq_map = HashMap::new();
        let expression_vec = vec![
            vec![20.0, 22.0, 60.0, 80.0, 90.0],
            vec![18.0, 25.0, 58.0, 78.0, 88.0],
        ];
        let gene_vec = vec!["gene1".to_string(), "gene2".to_string()];
        let samples = vec![
            "sample1".to_string(),
            "sample2".to_string(),
            "sample3".to_string(),
            "sample4".to_string(),
            "sample5".to_string(),
        ];
        group_uniq_map.insert("group1".to_string(), vec![0, 1, 2, 3, 4]);
        // group_uniq_map.insert("group2".to_string(), vec![2, 3, 4]);

        // let res =
        //     core_threads_cacu_pearson(&group_uniq_map, expression_vec, &gene_vec, &samples, false);
        // // let res =
        // //     core_threads_cacu_pearson1(&group_uniq_map, &expression_vec, &gene_vec, &samples, true)
        // //         .unwrap();
        // println!("{:?}", res);
        let duration = start.elapsed();
        println!("耗时: {:?}", duration);
    }
    #[test]
    fn test_nd_array_speed() {
        // 创建一个大小为 (rows, cols) 的二维数组
        let rows = 1000;
        let cols = 1000;
        let mut data = vec![0.0; rows * cols];
        for i in 0..rows {
            for j in 0..cols {
                data[i * cols + j] = (i * cols + j) as f64;
            }
        }
        let array = Array2::from_shape_vec((rows, cols), data).unwrap();

        // 测试使用 ndarray 进行列选择的性能
        let start_ndarray = std::time::Instant::now();
        let selected_column_ndarray = array.select(Axis(1), &[0, 1, 2, 3, 4]);
        let elapsed_ndarray = start_ndarray.elapsed();

        // 打印结果
        println!(
            "Time taken by ndarray column selection: {:?}",
            elapsed_ndarray
        );

        // 测试手写 Rust 代码进行列选择的性能
        let start_custom = std::time::Instant::now();
        let selected_columns_custom: Vec<_> = (0..rows).map(|i| array.slice(s![i, 0..4])).collect();
        let elapsed_custom = start_custom.elapsed();

        // 打印结果
        println!(
            "Time taken by custom Rust column selection: {:?}",
            elapsed_custom
        );
    }
    #[test]
    fn test_pearson_pvalue() {
        // 示例数据
        let x = vec![
            2.1721040,
            1.4957743,
            -0.1816372,
            0.6784108,
            -0.009159269,
            0.1361455,
        ];
        let y = vec![
            -0.2802266,
            -1.2365491,
            0.1526231,
            -0.4330048,
            -0.700769586,
            2.6310295,
        ];

        // 计算 Pearson 相关系数
        let correlation_coefficient = pearson_correlation_vec(&x, &y);

        // 计算 t 统计量
        let t_stat = t_statistic(correlation_coefficient, x.len());

        // 计算 p-value
        let p_value = p_value(t_stat, x.len());

        println!(
            "Pearson Correlation Coefficient: {}",
            correlation_coefficient
        );
        println!("T Statistic: {}", t_stat);
        println!("P-Value: {}", p_value);
    }
    #[test]
    fn test_get_cpu() {
        let max_cpus = num_cpus::get();

        println!("System's max core count: {}", max_cpus);
        let rayon_nums = (max_cpus as f64 * 62.0 / 64.0) as usize;
        println!("Use count: {}", rayon_nums);
    }
}
