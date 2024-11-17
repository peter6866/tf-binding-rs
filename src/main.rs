// use polars::prelude::*;
use tf_binding_rs::occupancy;

fn main() {
    // let df = occupancy::
    let seq = "GAGCCGGGTCATGAAAAAGGGGATCTTGTGTGTCTGTCCACGATAAGCACTATCACAAGGACTTTCTATAAACTCACAAGAAATTTCTGCCCACCCAGCACACAGTTTGTCCAGCTCATCCTGTAGGTGTCTCTATAATAGGACCTATCATAAAAAATTCCTCAAGACTGCAGCATTTCAGATAAGCCACCCTCACAAGA";
    let ewms = occupancy::read_pwm_to_ewm("tests/data/tdmMotifs.meme").unwrap();
    let landscape = occupancy::total_landscape(seq, &ewms, 8.0).unwrap();
    println!("{:?}", landscape);
}
