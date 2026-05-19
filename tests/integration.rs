/// Integration tests for the bedpull library crate.
/// These tests exercise the public API exactly as a downstream crate would use it.
use bedpull::{
    BamConfig, CigarOps, PafIndex, PafRecord, ReadCuts, ToCigarOps,
    utils::{calculate_qscore, get_read_cuts, read_bed, write_fasta_record, write_fastq_record},
};
use std::io::Write;

// --- re-exports reachable ---

#[test]
fn bam_config_accessible_from_crate_root() {
    let c = BamConfig::default();
    assert_eq!(c.min_mapq, 0);
}

#[test]
fn to_cigar_ops_trait_accessible_from_crate_root() {
    let ops: CigarOps = "4M2I3M".to_cigar_ops().expect("valid CIGAR");
    assert_eq!(ops.len(), 3);
}

#[test]
fn read_cuts_type_accessible() {
    let ops = "10M".to_cigar_ops().unwrap();
    let _cuts: ReadCuts = get_read_cuts(&ops, 1, 3, 7);
}

// --- PafIndex and PafRecord via lib root ---

fn temp_paf(contents: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", contents).unwrap();
    f
}

const PAF_LINE: &str = "q1\t100\t0\t50\t+\tchr1\t1000\t100\t150\t50\t50\t60\tcg:Z:50M";

#[test]
fn paf_index_build_and_query_via_lib() {
    let paf = temp_paf(&format!("{PAF_LINE}\n"));
    let idx = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
    let hits = idx.query("chr1", 100, 150);
    assert_eq!(hits.len(), 1);
}

#[test]
fn paf_record_from_line_via_lib() {
    let r = PafRecord::from_line(PAF_LINE).unwrap();
    assert_eq!(r.query_name, "q1");
    assert_eq!(r.cigar.as_deref(), Some("50M"));
}

// --- read_bed ---

fn temp_bed(contents: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", contents).unwrap();
    f
}

#[test]
fn read_bed_returns_correct_regions() {
    let f = temp_bed("chr1\t100\t200\nchr4\t39318077\t39318136\tRFC1\n");
    let regions = read_bed(f.path()).unwrap();
    assert_eq!(regions.len(), 2);
    assert_eq!(regions[1].1, "RFC1");
    assert_eq!(regions[1].2, "chr4");
}

// --- write_fasta_record / write_fastq_record ---

#[test]
fn write_fasta_produces_correct_output() {
    let mut buf = Vec::new();
    write_fasta_record(&mut buf, "read1|chr1:100-200|region", "ACGT").unwrap();
    assert_eq!(
        String::from_utf8(buf).unwrap(),
        ">read1|chr1:100-200|region\nACGT\n"
    );
}

#[test]
fn write_fastq_produces_correct_output() {
    let mut buf = Vec::new();
    write_fastq_record(&mut buf, "read1", "ACGT", "IIII").unwrap();
    assert_eq!(
        String::from_utf8(buf).unwrap(),
        "@read1\nACGT\n+\nIIII\n"
    );
}

// --- calculate_qscore ---

#[test]
fn qscore_uniform_phred40() {
    let q: String = std::iter::repeat('I').take(5).collect();
    assert!((calculate_qscore(&q) - 40.0).abs() < 0.01);
}

#[test]
fn qscore_below_threshold_detection() {
    // Phred 0 ('!') should give score of 0.0
    assert!((calculate_qscore("!") - 0.0).abs() < 0.01);
    // Phred 40 ('I') should pass a threshold of 20
    let q: String = std::iter::repeat('I').take(5).collect();
    assert!(calculate_qscore(&q) > 20.0);
}

// --- get_read_cuts end-to-end ---

#[test]
fn get_read_cuts_insertion_captured() {
    let ops = "3M5I4M".to_cigar_ops().unwrap();
    let cuts = get_read_cuts(&ops, 1, 3, 7);
    assert_eq!(cuts.read_end - cuts.read_start, 9); // 1M + 5I + 3M
}
