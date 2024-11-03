use polars::prelude::*;
use tf_binding_rs::fasta;

fn main() {
    let df: DataFrame = fasta::read_fasta("tests/data/test1.fasta").unwrap();
    let restrictions = ["AGCTTTTTAATAGAGTCAGCAAAACTGAA", "TCACTGA"];
    let has_restriction_sites_df = fasta::has_restriction_sites(&df, &restrictions).unwrap();
    println!("{:?}", has_restriction_sites_df);
}
