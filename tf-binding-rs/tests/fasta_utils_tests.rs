use polars::prelude::*;
use tf_binding_rs::fasta;

#[test]
fn test_read_fasta() {
    let path = "tests/data/test1.fasta";
    let df = fasta::read_fasta(path).unwrap();
    assert_eq!(df.height(), 3);
    assert_eq!(df.width(), 2);

    // test file does not exist
    let result = fasta::read_fasta("tests/data/nonexistent.fasta");
    assert!(result.is_err());
}

#[test]
fn test_write_fasta() {
    let path = "tests/data/test1_out.fasta";
    let df: DataFrame = df!(
        "label" => ["chr1-4357766-4357930_CPPP_WT", "chr1-4357733-4357765_CPPP_WT", "chr1-4357712-4357732_CPPP_WT"],
        "sequence" => ["AGCTTTTTAATAGAGTCAGCAAAACTGAAGCCT", "TGCTTTTTTTTTGAGTCAGCAAAACTGAAGCCT", "CGCTTTTTAATAGAGTCAGCAAAACTGAAGCCT"],
    )
    .unwrap();

    fasta::write_fasta(&df, path).unwrap();

    let df_out = fasta::read_fasta(path).unwrap();
    assert_eq!(df_out.height(), 3);
    assert_eq!(df_out.width(), 2);

    // clean up
    std::fs::remove_file(path).unwrap();
}

#[test]
fn test_reverse_complement() {
    // Test basic reverse complement
    let sequence = "ATCG";
    assert_eq!(fasta::reverse_complement(sequence).unwrap(), "CGAT");

    // Test longer sequence
    let sequence = "AATTCCGG";
    assert_eq!(fasta::reverse_complement(sequence).unwrap(), "CCGGAATT");

    // Test palindromic sequence
    let sequence = "GCGC";
    assert_eq!(fasta::reverse_complement(sequence).unwrap(), "GCGC");

    // Test error case with invalid nucleotide
    let sequence = "ATCGX";
    assert!(fasta::reverse_complement(sequence).is_err());
}
