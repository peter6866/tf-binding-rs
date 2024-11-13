use polars::prelude::*;
use tf_binding_rs::occupancy;

fn main() {
    let df = occupancy::read_pwm_files("tests/data/tdmMotifs.meme").unwrap();
    println!("{:?}", df);
}
