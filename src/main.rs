// use polars::prelude::*;
use tf_binding_rs::occupancy;

fn main() {
    let df = occupancy::read_pwm_to_ewm("tests/data/tdmMotifs.meme").unwrap();

    let pwm1 = df.get("NRL_HUMAN.MA0842.1").unwrap();
    println!("{:?}", pwm1);
}
