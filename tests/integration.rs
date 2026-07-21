/// Integration tests for the bedpull library crate.
/// These tests exercise the public API exactly as a downstream crate would use it.
use bedpull::{
    BamConfig, CigarOps, PafIndex, PafRecord, ReadCuts, ToCigarOps,
    utils::{calculate_qscore, get_read_cuts, read_bed, write_fasta_record, write_fastq_record},
};
use std::io::Write;
use std::path::Path;

// --- CRAM helpers ---

fn cram_query(
    cram_path: &Path,
    region_str: &str,
) -> (
    noodles::cram::io::IndexedReader<std::fs::File>,
    noodles::sam::Header,
    noodles::core::Region,
) {
    let region: noodles::core::Region = region_str.parse().unwrap();
    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .build_from_path(cram_path)
        .expect("failed to open test CRAM");
    let header = reader.read_header().expect("failed to read CRAM header");
    (reader, header, region)
}

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
    // 10M from align_start=1 → align_end (one-past) = 11.
    let _cuts: ReadCuts = get_read_cuts(&ops, 1, 11, 3, 7);
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
    let regions = read_bed(f.path(), false).unwrap();
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
    assert_eq!(String::from_utf8(buf).unwrap(), "@read1\nACGT\n+\nIIII\n");
}

// --- calculate_qscore ---

#[test]
fn qscore_uniform_phred40() {
    let q: String = "I".repeat(5);
    assert!((calculate_qscore(&q) - 40.0).abs() < 0.01);
}

#[test]
fn qscore_below_threshold_detection() {
    // Phred 0 ('!') should give score of 0.0
    assert!((calculate_qscore("!") - 0.0).abs() < 0.01);
    // Phred 40 ('I') should pass a threshold of 20
    let q: String = "I".repeat(5);
    assert!(calculate_qscore(&q) > 20.0);
}

// --- get_read_cuts end-to-end ---

#[test]
fn get_read_cuts_insertion_captured() {
    let ops = "3M5I4M".to_cigar_ops().unwrap();
    // 3M+4M consume 7 ref bases from align_start=1 → align_end (one-past) = 8.
    let cuts = get_read_cuts(&ops, 1, 8, 3, 7);
    assert_eq!(cuts.read_end - cuts.read_start, 9); // 1M + 5I + 3M
}

// --- BAM mode end-to-end ---

fn bam_query(
    bam_path: &Path,
    region_str: &str,
) -> (
    noodles::bam::io::indexed_reader::IndexedReader<
        impl noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
    >,
    noodles::sam::Header,
    noodles::core::Region,
) {
    let region: noodles::core::Region = region_str.parse().unwrap();
    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .expect("failed to open test BAM");
    let header = reader.read_header().expect("failed to read BAM header");
    (reader, header, region)
}

#[test]
fn bam_mode_rfc1_returns_reads() {
    let bam_path = Path::new("examples/rfc1_test.bam");
    let (mut reader, header, region) = bam_query(bam_path, "chr4:39318077-39318136");
    let query = reader.query(&header, &region).expect("BAM query failed");

    let (reads, _) = bedpull::get_bam_reads(&BamConfig::default(), query, &region, 0, 0)
        .expect("get_bam_reads failed");

    // 40 reads fully span chr4:39318077-39318136 in the test BAM
    assert_eq!(reads.len(), 40, "expected 40 spanning reads");
    for (name, seq, _qual, _rs, _re, _hap) in &reads {
        assert!(!name.is_empty());
        assert!(!seq.is_empty());
    }
}

#[test]
fn bam_mode_candidates_seen_counts_all_overlapping_reads_before_filtering() {
    // candidates_seen should count every read the query yields, including the 4 that
    // only partially overlap and get excluded by resolve_cuts in non-partial mode —
    // i.e. candidates_seen (44) > reads.len() (40), distinguishing "nothing overlapped"
    // from "some candidates were filtered out" for unmapped-region reporting.
    let bam_path = Path::new("examples/rfc1_test.bam");
    let (mut reader, header, region) = bam_query(bam_path, "chr4:39318077-39318136");
    let query = reader.query(&header, &region).expect("BAM query failed");

    let (reads, candidates_seen) =
        bedpull::get_bam_reads(&BamConfig::default(), query, &region, 0, 0)
            .expect("get_bam_reads failed");

    assert_eq!(reads.len(), 40);
    assert_eq!(candidates_seen, 44);
}

#[test]
fn bam_mode_partial_includes_more_reads() {
    let bam_path = Path::new("examples/rfc1_test.bam");
    let (mut reader, header, region) = bam_query(bam_path, "chr4:39318077-39318136");
    let query = reader.query(&header, &region).expect("BAM query failed");

    let config = BamConfig {
        partial: true,
        ..BamConfig::default()
    };
    let (reads, _) =
        bedpull::get_bam_reads(&config, query, &region, 0, 0).expect("get_bam_reads failed");

    // 44 reads overlap the region with --partial (40 spanning + 4 partial)
    assert_eq!(reads.len(), 44, "expected 44 reads with --partial");
    assert!(
        reads.len() > 40,
        "--partial should return more reads than spanning-only"
    );
}

#[test]
fn bam_mode_flanks_extend_extracted_sequence() {
    let bam_path = Path::new("examples/rfc1_test.bam");

    let (mut r0, h0, region0) = bam_query(bam_path, "chr4:39318077-39318136");
    let q0 = r0.query(&h0, &region0).unwrap();
    let (reads_no_flank, _) =
        bedpull::get_bam_reads(&BamConfig::default(), q0, &region0, 0, 0).unwrap();

    let (mut r1, h1, region1) = bam_query(bam_path, "chr4:39318077-39318136");
    let q1 = r1.query(&h1, &region1).unwrap();
    let (reads_with_flank, _) =
        bedpull::get_bam_reads(&BamConfig::default(), q1, &region1, 100, 100).unwrap();

    // Flanks can pull in reads whose effective window is clamped to alignment end,
    // so count can increase; sequences are always wider or equal per read.
    assert!(
        reads_with_flank.len() >= reads_no_flank.len(),
        "flanks should not reduce read count"
    );
    let total_no_flank: usize = reads_no_flank.iter().map(|(_, s, ..)| s.len()).sum();
    let total_with_flank: usize = reads_with_flank.iter().map(|(_, s, ..)| s.len()).sum();
    assert!(
        total_with_flank > total_no_flank,
        "flanked extraction ({total_with_flank} bp) should exceed unflanked ({total_no_flank} bp)"
    );
}

// --- CRAM mode end-to-end ---

#[test]
fn cram_mode_rfc1_returns_same_reads_as_bam() {
    let cram_path = Path::new("examples/rfc1_test.cram");
    let (mut reader, header, region) = cram_query(cram_path, "chr4:39318077-39318136");
    let query = reader.query(&header, &region).expect("CRAM query failed");

    let (reads, _) = bedpull::get_cram_reads(&BamConfig::default(), query, &region, 0, 0)
        .expect("get_cram_reads failed");

    // Should match the 40 spanning reads returned by BAM mode
    assert_eq!(reads.len(), 40, "expected 40 spanning reads from CRAM");
    for (name, seq, _qual, _rs, _re, _hap) in &reads {
        assert!(!name.is_empty());
        assert!(!seq.is_empty());
    }
}

#[test]
fn cram_mode_partial_includes_more_reads() {
    let cram_path = Path::new("examples/rfc1_test.cram");
    let (mut reader, header, region) = cram_query(cram_path, "chr4:39318077-39318136");
    let query = reader.query(&header, &region).expect("CRAM query failed");

    let config = BamConfig {
        partial: true,
        ..BamConfig::default()
    };
    let (reads, _) =
        bedpull::get_cram_reads(&config, query, &region, 0, 0).expect("get_cram_reads failed");

    assert_eq!(
        reads.len(),
        44,
        "expected 44 reads with --partial from CRAM"
    );
    assert!(
        reads.len() > 40,
        "--partial should return more reads than spanning-only"
    );
}

#[test]
fn cram_mode_flanks_extend_extracted_sequence() {
    let cram_path = Path::new("examples/rfc1_test.cram");

    let (mut r0, h0, region0) = cram_query(cram_path, "chr4:39318077-39318136");
    let q0 = r0.query(&h0, &region0).unwrap();
    let (reads_no_flank, _) =
        bedpull::get_cram_reads(&BamConfig::default(), q0, &region0, 0, 0).unwrap();

    let (mut r1, h1, region1) = cram_query(cram_path, "chr4:39318077-39318136");
    let q1 = r1.query(&h1, &region1).unwrap();
    let (reads_with_flank, _) =
        bedpull::get_cram_reads(&BamConfig::default(), q1, &region1, 100, 100).unwrap();

    assert!(
        reads_with_flank.len() >= reads_no_flank.len(),
        "flanks should not reduce read count"
    );
    let total_no_flank: usize = reads_no_flank.iter().map(|(_, s, ..)| s.len()).sum();
    let total_with_flank: usize = reads_with_flank.iter().map(|(_, s, ..)| s.len()).sum();
    assert!(
        total_with_flank > total_no_flank,
        "flanked CRAM extraction ({total_with_flank} bp) should exceed unflanked ({total_no_flank} bp)"
    );
}

#[test]
fn cram_mode_sequences_match_bam_mode() {
    let bam_path = Path::new("examples/rfc1_test.bam");
    let (mut br, bh, bregion) = bam_query(bam_path, "chr4:39318077-39318136");
    let bq = br.query(&bh, &bregion).unwrap();
    let (mut bam_reads, _) =
        bedpull::get_bam_reads(&BamConfig::default(), bq, &bregion, 0, 0).unwrap();

    let cram_path = Path::new("examples/rfc1_test.cram");
    let (mut cr, ch, cregion) = cram_query(cram_path, "chr4:39318077-39318136");
    let cq = cr.query(&ch, &cregion).unwrap();
    let (mut cram_reads, _) =
        bedpull::get_cram_reads(&BamConfig::default(), cq, &cregion, 0, 0).unwrap();

    assert_eq!(
        bam_reads.len(),
        cram_reads.len(),
        "read counts differ between BAM and CRAM"
    );

    bam_reads.sort_by(|a, b| a.0.cmp(&b.0));
    cram_reads.sort_by(|a, b| a.0.cmp(&b.0));

    for ((bname, bseq, ..), (cname, cseq, ..)) in bam_reads.iter().zip(cram_reads.iter()) {
        assert_eq!(bname, cname, "read names differ");
        assert_eq!(bseq, cseq, "sequences differ for read {bname}");
    }
}
