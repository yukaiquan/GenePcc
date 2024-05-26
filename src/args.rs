use ansi_term::Colour;
use clap_v3::{App, Arg, ArgMatches};

pub fn get_arg() -> ArgMatches {
    println!(
        "{} {}",
        Colour::Green.paint("Welcome to use caculate PCC!"),
        Colour::Red.paint("Author: Yu kaiquan <1962568272@qq.com>")
    );
    App::new("cst_gene_net")
        .version("0.0.1")
        .author("Yu kaiquan <1962568272@qq.com>")
        .about("Caculate PCCs from TPM/FPKM/CPM matrix")
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("exp matrix file(csv)")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("samples")
                .short('s')
                .long("samples")
                .value_name("STRING")
                .help("samples file(csv)")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short('o')
                .long("output")
                .value_name("STRING")
                .help("output dir(default:cst_gene_net_out)")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::with_name("pearson_threshold")
                .short('p')
                .long("pearson threshold")
                .value_name("STRING")
                .help("pearson threshold(default:0.8)")
                .takes_value(true)
                .required(false),
        )
        .get_matches()
}
