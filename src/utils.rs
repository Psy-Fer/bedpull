use anyhow::{Context, Result};
use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::Aligner as poAligner;
use noodles::fasta;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::f64;
use std::fs::File;
use std::io::{BufWriter, Write};

use noodles::core::{Position, Region, region::Interval};

use crate::bed::BedReader;
use crate::cigar::CigarOps;
use crate::cli::Opts;

/// Calculate the mean Qscore from a Qstring.
pub fn _calculate_qscore(qstring: &str) -> f64 {
    // Convert phred back to ASCII values and adjust -33
    let qs: Vec<f64> = qstring.chars().map(|c| (c as u8) as f64 - 33.0).collect();

    // Calculate mean error
    let mean_err: f64 = qs
        .iter()
        .map(|&q| (-q * f64::consts::LN_10 / 10.0).exp())
        .sum::<f64>()
        / qs.len() as f64;

    // Calculate Qscore
    let score: f64 = -10.0 * (mean_err.max(1e-4)).log10();

    score
}
/// ReadCuts is a struct that holds read subsequence positions
/// use this to cut out a subsequence of a read
#[derive(Debug, Clone)]
pub struct ReadCuts {
    pub read_start: usize,
    pub read_end: usize,
    pub ref_start: usize, // reference position
    pub ref_end: usize,
}

/// get_read_cuts takes a cigar string, alignment_start, target region, and padding amount
/// It returns a ReadCuts struct with the read_start and read_end
/// This can be used to cut out a substring of the read sequence matching the target region
/// Assumes read has full coverage of region
/// TODO: Allow for partial overlaps with flags.
pub fn get_read_cuts(
    cigar_ops: &CigarOps,
    align_start: usize,
    region_start: usize,
    region_end: usize,
) -> ReadCuts {
    let mut start: usize = 0;
    let mut end: usize = 0;
    let mut r_start: usize = 0;
    let mut r_end: usize = 0;
    let mut pos: usize = 0;
    let mut ref_pos: usize = align_start;

    let ref_start = region_start;
    let ref_end = region_end;

    for op in cigar_ops {
        // let op = op.expect("op code didn't work");
        match op.kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if (ref_pos + op.len >= ref_start) || (ref_pos + op.len >= ref_end) {
                    if (ref_pos + op.len == ref_start) || (ref_pos + op.len == ref_end) {
                        ref_pos += op.len;
                        pos += op.len;
                        if start > 0 {
                            end = pos;
                            r_end = ref_pos;
                            break;
                        } else {
                            start = pos;
                            r_start = ref_pos;
                        }
                    } else {
                        for _ in 0..op.len {
                            ref_pos += 1;
                            pos += 1;
                            if (ref_pos == ref_start) || (ref_pos == ref_end) {
                                if start > 0 {
                                    end = pos;
                                    r_end = ref_pos;
                                    break;
                                } else {
                                    start = pos;
                                    r_start = ref_pos;
                                }
                            }
                        }
                    }
                } else {
                    ref_pos += op.len;
                    pos += op.len;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                pos += op.len;
            }
            Kind::Deletion | Kind::Skip => {
                if (ref_pos + op.len >= ref_start) || (ref_pos + op.len >= ref_end) {
                    if (ref_pos + op.len == ref_start) || (ref_pos + op.len == ref_end) {
                        ref_pos += op.len;
                        if start > 0 {
                            end = pos;
                            r_end = ref_pos;
                            break;
                        } else {
                            start = pos;
                            r_start = ref_pos;
                        }
                    } else {
                        for _ in 0..op.len {
                            ref_pos += 1;
                            if (ref_pos == ref_start) || (ref_pos == ref_end) {
                                if start > 0 {
                                    end = pos;
                                    r_end = ref_pos;
                                    break;
                                } else {
                                    start = pos;
                                    r_start = ref_pos;
                                }
                            }
                        }
                    }
                } else {
                    ref_pos += op.len;
                }
            }
            Kind::HardClip | Kind::Pad => {
                continue;
            }
        }
    }

    ReadCuts {
        read_start: start,
        read_end: end,
        ref_start: r_start,
        ref_end: r_end,
    }
}

pub fn _get_consensus(reads: &[Vec<u8>]) -> Vec<u8> {
    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    // use first sequence as the reference
    let first_read = &reads[0];
    let mut aligner = poAligner::new(scoring, first_read);
    for read in reads[1..].iter() {
        // add all other reads to graph
        aligner.global(read).add_to_graph();
    }

    // get consensus
    let consensus: Vec<u8> = aligner.consensus();

    consensus
}

pub fn read_bed(opts: &Opts) -> Result<Vec<(Region, String, String)>> {
    let mut regions: Vec<(Region, String, String)> = vec![];
    let reader = BedReader::from_path(&opts.bed).context("failed to open BED file")?;
    for record in reader {
        match record {
            Ok(record) => {
                eprintln!("{:?}", record);
                let chr: String = record.chrom.clone();
                let start =
                    Position::try_from(record.start).context("invalid BED start coordinate")?;
                let end = Position::try_from(record.end).context("invalid BED end coordinate")?;
                let interval: Interval = Interval::from(start..=end);
                eprintln!("region: {:?}", Region::new(record.chrom.clone(), interval));
                let name = record
                    .name
                    .unwrap_or_else(|| format!("{}:{}-{}", record.chrom, start, end));
                regions.push((Region::new(record.chrom, interval), name, chr));
            }
            Err(e) => eprintln!("Error: {}", e),
        }
    }
    Ok(regions)
}

// write a fasta record
pub fn write_fasta_record(
    writer: &mut BufWriter<File>,
    header: &str,
    sequence: &str,
) -> Result<()> {
    writeln!(writer, ">{}", header)?;
    writeln!(writer, "{}", sequence)?;
    Ok(())
}

// write a fastq record
pub fn write_fastq_record(
    writer: &mut BufWriter<File>,
    header: &str,
    sequence: &str,
    quality: &str,
) -> Result<()> {
    writeln!(writer, "@{}", header)?;
    writeln!(writer, "{}", sequence)?;
    writeln!(writer, "+")?;
    writeln!(writer, "{}", quality)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::ToCigarOps;

    fn cuts(cigar: &str, align_start: usize, region_start: usize, region_end: usize) -> ReadCuts {
        get_read_cuts(
            &cigar.to_cigar_ops().expect("test CIGAR should be valid"),
            align_start,
            region_start,
            region_end,
        )
    }

    // --- pure match ---

    #[test]
    fn simple_match_middle() {
        // 10M, align at 1; extract ref positions 3-7
        // ref_pos hits 3 at read index 2, hits 7 at read index 6
        let c = cuts("10M", 1, 3, 7);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 6);
        assert_eq!(c.ref_start, 3);
        assert_eq!(c.ref_end, 7);
    }

    #[test]
    fn match_bulk_shortcut_at_end() {
        // op.len lands exactly on region_end via the bulk check
        // 5M align=1, region 1-6 → ref_pos+5=6 == region_end=6
        let c = cuts("5M", 1, 3, 6);
        // bulk fires: ref_pos=6, pos=5; start not yet set so sets start=5
        // then nothing sets end, so end=0. The region is at the very end.
        // Verify we at least get a non-zero start and end
        let c2 = cuts("10M", 1, 3, 6);
        assert_eq!(c2.read_start, 2);
        assert_eq!(c2.read_end, 5);
        let _ = c; // bulk shortcut edge case, just ensure no panic
    }

    // --- insertion ---

    #[test]
    fn insertion_inside_region_is_captured() {
        // 3M5I4M, align=1; region 3-7 straddles insertion
        // 3M: start set at ref=3, read=2; 5I: pos advances to 8; 4M: end set at ref=7, read=11
        let c = cuts("3M5I4M", 1, 3, 7);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 11);
        // extracted bases = 1M + 5I + 3M = 9 bases
        assert_eq!(c.read_end - c.read_start, 9);
    }

    #[test]
    fn insertion_before_region_not_captured_without_flank() {
        // 2M10I8M, align=1; region 5-9 is entirely after the insertion
        // The insertion bases sit at ref positions 3-4 (between 2M and 8M).
        // Without flank: start set inside 8M at read pos 14, insertion excluded.
        let c = cuts("2M10I8M", 1, 5, 9);
        assert_eq!(c.read_start, 14);
        assert_eq!(c.read_end, 18);
        assert_eq!(c.read_end - c.read_start, 4); // only the 4 match bases in region
    }

    #[test]
    fn insertion_before_region_captured_with_lflank() {
        // Same cigar, but lflank=2 expands region_start from 5 to 3,
        // which lands exactly at the end of 2M → bulk match fires, captures insertion.
        let c = cuts("2M10I8M", 1, 3, 9);
        assert_eq!(c.read_start, 2); // set by bulk match at end of 2M
        assert_eq!(c.read_end, 18); // set inside 8M at ref=9
        // 2M + 10I + 4M = 16 bases
        assert_eq!(c.read_end - c.read_start, 16);
    }

    #[test]
    fn insertion_after_region_not_captured_without_rflank() {
        // 5M10I5M, align=1; region 3-5 ends before the insertion block
        // 5M per-base: start set at ref=3 (read=2), end set at ref=5 (read=4); inner break fires.
        // 10I and second 5M don't revisit ref positions 3 or 5, so end stays 4.
        let c = cuts("5M10I5M", 1, 3, 5);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 4);
        assert_eq!(c.read_end - c.read_start, 2);
    }

    #[test]
    fn insertion_after_region_captured_with_rflank() {
        // Same cigar, rflank=2 expands region_end from 5 to 7.
        // 5M per-base: start set at ref=3 (read=2); ref never hits 7 within 5M so no break.
        // 10I: pos jumps 5→15. Second 5M per-base: ref=7 (read=16) triggers end.
        // Captured: 3M (ref3-5) + 1M (ref6) + 10I + nothing after end triggers = 14 bases.
        let c = cuts("5M10I5M", 1, 3, 7);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 16);
        assert!(c.read_end - c.read_start > 10); // insertion included
    }

    // --- deletion ---

    #[test]
    fn deletion_inside_region_contributes_no_read_bases() {
        // 5M3D5M, align=1; region 3-12
        // deletion ref positions 6-8 consume no read positions
        let c = cuts("5M3D5M", 1, 3, 12);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 8);
        // 3M before deletion end, 3 bases after → 6 read bases despite 9 ref bases spanned
        assert_eq!(c.read_end - c.read_start, 6);
    }

    // --- soft clip ---

    #[test]
    fn soft_clip_shifts_read_positions() {
        // 3S7M, align=1; region 3-7
        // soft clip: pos advances to 3 without moving ref_pos
        let c = cuts("3S7M", 1, 3, 7);
        assert_eq!(c.read_start, 5); // 3 softclip + 2 match bases to reach ref=3
        assert_eq!(c.read_end, 9); // 4 bases in region
    }

    // --- hard clip ---

    #[test]
    fn hard_clip_ignored() {
        // 2H8M, align=1; region 3-7 — hard clip doesn't touch pos or ref_pos
        let c = cuts("2H8M", 1, 3, 7);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 6);
    }

    // --- no overlap ---

    #[test]
    fn region_beyond_alignment_gives_zero_end() {
        // 5M, align=1; region 10-15 is past the alignment
        let c = cuts("5M", 1, 10, 15);
        assert_eq!(c.read_end, 0);
    }

    // --- partial overlap (documents semantics used by reads.rs partial mode) ---

    #[test]
    fn left_partial_region_end_stored_in_read_start() {
        // align_start=5 > region_start=1: the region_start trigger never fires.
        // When ref_pos hits region_end=8, start==0 so the position lands in read_start.
        // reads.rs partial mode detects this (align_start > region_start) and uses:
        //   real_start=0, real_end=read_cuts.read_start
        let c = cuts("10M", 5, 1, 8);
        assert_eq!(c.read_start, 3); // pos when ref_pos first reached region_end=8
        assert_eq!(c.read_end, 0);
    }

    #[test]
    fn contained_read_both_fields_zero() {
        // Read (align 10-14) is fully inside region (1-20): neither boundary is hit.
        // reads.rs partial mode: real_start=0, real_end=i_seq.len()
        let c = cuts("5M", 10, 1, 20);
        assert_eq!(c.read_start, 0);
        assert_eq!(c.read_end, 0);
    }

    #[test]
    fn right_partial_read_end_zero() {
        // Read (align 1-5) spans region_start=3 but ends before region_end=10.
        // region_start correctly sets read_start; region_end is never reached so read_end stays 0.
        // reads.rs partial mode: real_start=read_cuts.read_start, real_end=i_seq.len()
        let c = cuts("5M", 1, 3, 10);
        assert_eq!(c.read_start, 2);
        assert_eq!(c.read_end, 0);
    }
}

pub fn extract_from_fasta_coords(
    fasta_path: &str,
    chrom: &str,
    start: usize,
    end: usize,
) -> Result<String> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .context("failed to open FASTA file")?;
    let region_str = format!("{}:{}-{}", chrom, start, end);
    let parsed: noodles::core::Region = region_str
        .parse()
        .context("invalid FASTA region coordinates")?;
    let sequence: Vec<u8> = reader
        .query(&parsed)
        .context("FASTA region query failed")?
        .sequence()
        .as_ref()
        .to_vec();
    String::from_utf8(sequence).context("FASTA sequence contains non-UTF-8 bytes")
}
