use plotlib::page::Page;
use plotlib::repr::{Histogram, HistogramBins};
use plotlib::style::BoxStyle;
use plotlib::view::ContinuousView;

use rand::prelude::*;
use std::collections::HashMap;

pub fn plot_histogram(data: &[f64], file_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    // let data = [0.3, 0.5, 6.4, 5.3, 3.6, 3.6, 3.5, 7.5, 4.0];
    let h = Histogram::from_slice(&data, HistogramBins::Count(40))
        // .style(&BoxStyle::new().fill("burlywood"));
        // 换个好看的颜色
        .style(&BoxStyle::new().fill("steelblue"));
    // println!("{:?}", h);
    let v = ContinuousView::new()
        .add(h)
        .x_label("Value")
        .y_label("Frequency")
        .x_range(-1.0, 1.0)
        //增加刻度横坐标刻度
        .x_max_ticks(40);

    Page::single(&v).save(file_name).expect("saving svg");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_plot_histogram() {
        // 示例数据，随机生成一些数据作为示例
        let mut rng = rand::thread_rng();
        let random_points: Vec<f64> = (0..1000).map(|_| rng.gen_range(-1.0..1.0)).collect();

        let file_name = "test_data/histogram.svg"; // 图片文件名

        if let Err(e) = plot_histogram(&random_points, file_name) {
            eprintln!("Error: {}", e);
        }
    }
}
