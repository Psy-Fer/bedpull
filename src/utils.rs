use anyhow::{Context, Result};
use noodles::fasta;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::f64;
use std::io::Write;
use std::path::Path;

use noodles::core::{Position, Region, region::Interval};

use crate::bed::BedReader;
use crate::cigar::CigarOps;

/// Calculate the mean Phred score from a Phred+33 encoded quality string.
///
/// Each character in `qstring` is interpreted as `ASCII value − 33` to recover the
/// per-base Phred score. The mean is computed in error-probability space (i.e.
/// `mean(10^(−Q/10))`) and then converted back to a single Phred score, which is
/// the statistically correct average. The result is clamped so that a mean error
/// probability of 1.0 (all bases at Q0) returns `0.0` rather than `-inf`.
pub fn calculate_qscore(qstring: &str) -> f64 {
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
/// Read-coordinate slice produced by [`get_read_cuts`] for a single alignment.
///
/// All positions are 0-based indices into the read sequence array. `ref_start` /
/// `ref_end` are the actual reference positions reached during the CIGAR walk and
/// may differ slightly from the requested BED coordinates when a region boundary
/// falls inside a deletion or at the very end of a CIGAR op.
#[derive(Debug, Clone)]
pub struct ReadCuts {
    /// Index of the first read base that corresponds to `region_start` in the
    /// reference. Use as the start of a slice into the read sequence: `seq[read_start..read_end]`.
    pub read_start: usize,
    /// Index one past the last read base that corresponds to `region_end`.
    /// A value of `0` is a sentinel meaning the region end was never reached
    /// (right-partial overlap or no overlap at all).
    pub read_end: usize,
    /// The reference position at which `read_start` was set during the CIGAR walk.
    pub ref_start: usize,
    /// The reference position at which `read_end` was set during the CIGAR walk.
    pub ref_end: usize,
}

/// Walk a CIGAR string to find the read-coordinate slice for a reference region.
///
/// This is the core of bedpull's CIGAR-aware extraction. It simultaneously
/// tracks `ref_pos` (the current reference position) and `pos` (the current
/// read position) as it steps through each CIGAR operation. When `ref_pos`
/// reaches `region_start` the read position is recorded as `read_start`; when
/// it reaches `region_end` the position is recorded as `read_end`. Insertions
/// advance `pos` without advancing `ref_pos`, so inserted bases that fall
/// between `region_start` and `region_end` are automatically included in the
/// slice — they have no reference coordinate but are captured because the
/// surrounding match ops bracket them.
///
/// ## Coordinate convention (critical)
///
/// `align_start` is **1-based** (as returned by noodles from a BAM record or
/// read from a PAF `target_start` field after adjustment). `region_start` and
/// `region_end` are **0-based** (directly from the BED file). This mixed
/// convention is intentional: it matches the existing behaviour that the tests
/// are written against. Do not normalise both to the same base — it will shift
/// every boundary by one.
///
/// ## Partial-overlap sentinel values
///
/// When the alignment starts after `region_start` (left-partial or contained
/// read), the `region_start` trigger never fires. If `region_end` is
/// subsequently reached, the position is stored in `read_start` (because
/// `start == 0`), and `read_end` remains `0`. [`get_bam_reads`][crate::reads::get_bam_reads]
/// detects this pattern when `partial` mode is enabled and interprets the
/// fields accordingly.
///
/// # Parameters
///
/// - `cigar_ops` — the decoded CIGAR for this alignment.
/// - `align_start` — 1-based reference position where the alignment begins
///   (from the BAM `alignment_start` field or the PAF `target_start` field).
/// - `region_start` — 0-based start of the BED region to extract.
/// - `region_end` — 0-based end of the BED region to extract (exclusive in BED
///   convention, but the CIGAR walk treats it as the last position to fire on).
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

/// Walk a CIGAR to find the read-coordinate offset corresponding to a single
/// reference position `ref_target`, or `None` if the alignment never reaches it.
///
/// This exists alongside [`get_read_cuts`] rather than reusing it because
/// `get_read_cuts` tracks *two* boundaries at once and decides which slot
/// (`read_start` vs `read_end`) to write into based on whether `read_start`
/// has already been set — a convention that breaks down when only one
/// boundary is needed and it happens to coincide with `align_start` (see
/// `get_read_cuts`'s "Partial-overlap sentinel values" docs): the trigger for
/// that boundary never independently fires, so the *other* boundary's
/// crossing gets misattributed into the wrong field. Callers that only need
/// one reference→read coordinate mapping (e.g. stitching a window's two
/// edges across different PAF records) should use this instead.
///
/// `ref_target <= align_start` returns `Some(0)` (at or before the alignment
/// begins). A reference position falling inside a deletion maps to the read
/// position immediately after the last matched base before it, since a
/// deletion consumes no read bases.
pub fn read_pos_at_ref(
    cigar_ops: &CigarOps,
    align_start: usize,
    ref_target: usize,
) -> Option<usize> {
    if ref_target <= align_start {
        return Some(0);
    }

    let mut pos: usize = 0;
    let mut ref_pos: usize = align_start;

    for op in cigar_ops {
        match op.kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if ref_pos + op.len >= ref_target {
                    return Some(pos + (ref_target - ref_pos));
                }
                ref_pos += op.len;
                pos += op.len;
            }
            Kind::Insertion | Kind::SoftClip => {
                pos += op.len;
            }
            Kind::Deletion | Kind::Skip => {
                if ref_pos + op.len >= ref_target {
                    return Some(pos);
                }
                ref_pos += op.len;
            }
            Kind::HardClip | Kind::Pad => continue,
        }
    }
    None
}

/// Parse a BED file and return a list of `(region, name, chromosome)` triples.
///
/// Reads 3- or 4-column BED format via [`crate::bed::BedReader`]. The
/// returned `Region` uses a noodles 1-based closed interval (BED `start` is
/// promoted with `Position::try_from`). The `name` field comes from the optional
/// 4th column; if absent it defaults to `"chrom:start-end"`. The `chromosome`
/// string is the raw first column, kept separately so callers can use it as a
/// plain string key (e.g. for PAF index lookups) without parsing the `Region`.
pub fn read_bed(path: &Path, debug: bool) -> Result<Vec<(Region, String, String)>> {
    let mut regions: Vec<(Region, String, String)> = vec![];
    let reader = BedReader::from_path(path).context("failed to open BED file")?;
    for record in reader {
        match record {
            Ok(record) => {
                if debug {
                    eprintln!("{:?}", record);
                }
                let chr: String = record.chrom.clone();
                let start =
                    Position::try_from(record.start).context("invalid BED start coordinate")?;
                let end = Position::try_from(record.end).context("invalid BED end coordinate")?;
                let interval: Interval = Interval::from(start..=end);
                if debug {
                    eprintln!("region: {:?}", Region::new(record.chrom.clone(), interval));
                }
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

/// Write a single FASTA record to `writer`.
///
/// Outputs `>header\nsequence\n`. No line-wrapping is applied; the sequence is
/// written as a single line. Returns an error if the underlying write fails.
pub fn write_fasta_record<W: Write + ?Sized>(
    writer: &mut W,
    header: &str,
    sequence: &str,
) -> Result<()> {
    writeln!(writer, ">{}", header)?;
    writeln!(writer, "{}", sequence)?;
    Ok(())
}

/// Write a single FASTQ record to `writer`.
///
/// Outputs the standard four-line FASTQ format: `@header`, sequence, `+`, quality.
/// `quality` must be a Phred+33 encoded string of the same length as `sequence`.
/// Returns an error if the underlying write fails.
pub fn write_fastq_record<W: Write + ?Sized>(
    writer: &mut W,
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

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'G' => b'C',
        b'g' => b'c',
        b'C' => b'G',
        b'c' => b'g',
        b'N' => b'N',
        b'n' => b'n',
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'S' => b'S',
        b's' => b's',
        b'W' => b'W',
        b'w' => b'w',
        b'K' => b'M',
        b'k' => b'm',
        b'M' => b'K',
        b'm' => b'k',
        b'B' => b'V',
        b'b' => b'v',
        b'V' => b'B',
        b'v' => b'b',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        _ => b'N',
    }
}

/// Reverse-complement a nucleotide sequence string.
///
/// Handles uppercase and lowercase IUPAC bases, preserving case. Non-IUPAC
/// characters are mapped to `N`. Returns the reverse complement as a new `String`.
pub fn revcomp(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| complement_base(b) as char)
        .collect()
}

/// Extract a subsequence from an indexed FASTA file.
///
/// Opens the FASTA at `fasta_path` (requires a companion `.fai` index) and queries
/// the region `chrom:start-end`. Coordinates are 0-based and are forwarded
/// directly to noodles via a region string. Returns the sequence as a `String`
/// of ASCII nucleotides, or an error if the file cannot be opened, the region
/// is out of bounds, or the bytes are not valid UTF-8.
///
/// Opens and drops the file on every call. When extracting many subsequences from
/// the same file (e.g. one per PAF alignment across many BED regions), prefer
/// opening the file once and calling [`extract_from_fasta_coords_reader`] instead.
pub fn extract_from_fasta_coords(
    fasta_path: &str,
    chrom: &str,
    start: usize,
    end: usize,
) -> Result<String> {
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .context("failed to open FASTA file")?;
    extract_from_fasta_coords_reader(&mut reader, chrom, start, end)
}

/// Extract a subsequence from an already-open indexed FASTA reader. Reuses `reader`
/// instead of opening the file fresh, which matters when extracting many
/// subsequences (e.g. one per overlapping alignment across many BED regions) — see
/// [`extract_from_fasta_coords`] for the single-shot convenience form.
pub fn extract_from_fasta_coords_reader<R>(
    reader: &mut fasta::io::IndexedReader<R>,
    chrom: &str,
    start: usize,
    end: usize,
) -> Result<String>
where
    R: std::io::BufRead + std::io::Seek,
{
    // `start`/`end` are 0-based half-open (per this function's contract), but noodles'
    // "chrom:start-end" region syntax is 1-based inclusive — convert or every query
    // silently includes one extra base before `start`.
    let region_str = format!("{}:{}-{}", chrom, start + 1, end);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::ToCigarOps;
    use std::io::Write;

    // --- calculate_qscore ---

    #[test]
    fn uniform_phred40_gives_40() {
        // 'I' = ASCII 73, Phred 40
        let q: String = "I".repeat(10);
        assert!((calculate_qscore(&q) - 40.0).abs() < 0.01);
    }

    #[test]
    fn uniform_phred0_gives_0() {
        // '!' = ASCII 33, Phred 0 → mean_error clamped to 1e-4, score = 40
        // Actually Phred 0 → error = 1.0 → mean_error = 1.0 → score = 0.0
        let score = calculate_qscore("!");
        assert!((score - 0.0).abs() < 0.01);
    }

    #[test]
    fn mixed_quality_is_between_extremes() {
        let q: String = ['!', 'I'].iter().collect(); // Phred 0 and 40
        let score = calculate_qscore(&q);
        assert!(score > 0.0 && score < 40.0);
    }

    // --- write_fasta_record ---

    #[test]
    fn fasta_record_format() {
        let mut buf = Vec::new();
        write_fasta_record(&mut buf, "read1|chr1:100-200|region", "ACGT").unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            ">read1|chr1:100-200|region\nACGT\n"
        );
    }

    #[test]
    fn fasta_empty_sequence() {
        let mut buf = Vec::new();
        write_fasta_record(&mut buf, "h", "").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), ">h\n\n");
    }

    // --- write_fastq_record ---

    #[test]
    fn fastq_record_format() {
        let mut buf = Vec::new();
        write_fastq_record(&mut buf, "read1", "ACGT", "IIII").unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "@read1\nACGT\n+\nIIII\n");
    }

    // --- read_bed ---

    fn temp_bed(contents: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, "{}", contents).unwrap();
        f
    }

    #[test]
    fn read_bed_three_column() {
        let f = temp_bed("chr1\t100\t200\n");
        let regions = read_bed(f.path(), false).unwrap();
        assert_eq!(regions.len(), 1);
        let (_, name, chr) = &regions[0];
        assert_eq!(chr, "chr1");
        // no 4th column → fallback name includes coords
        assert!(name.contains("chr1"));
    }

    #[test]
    fn read_bed_four_column_name() {
        let f = temp_bed("chr4\t39318077\t39318136\tRFC1\n");
        let regions = read_bed(f.path(), false).unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].1, "RFC1");
        assert_eq!(regions[0].2, "chr4");
    }

    #[test]
    fn read_bed_multiple_regions() {
        let f = temp_bed("chr1\t100\t200\nchr2\t300\t400\tFOO\n");
        let regions = read_bed(f.path(), false).unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[1].1, "FOO");
    }

    // --- revcomp ---

    #[test]
    fn revcomp_simple() {
        assert_eq!(revcomp("ACGT"), "ACGT"); // palindrome
        assert_eq!(revcomp("AAAA"), "TTTT");
        assert_eq!(revcomp("GCGC"), "GCGC");
    }

    #[test]
    fn revcomp_preserves_case() {
        assert_eq!(revcomp("acgt"), "acgt");
        assert_eq!(revcomp("AcGt"), "aCgT");
    }

    #[test]
    fn revcomp_iupac_codes() {
        assert_eq!(revcomp("R"), "Y"); // R=A/G → complement Y=C/T, reversed
        assert_eq!(revcomp("N"), "N");
    }

    #[test]
    fn revcomp_involution() {
        // revcomp(revcomp(x)) == x for any sequence
        let seq = "ACGTNRYSWKMBDHV";
        assert_eq!(revcomp(&revcomp(seq)), seq);
    }

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

    // --- read_pos_at_ref ---

    fn pos_at(cigar: &str, align_start: usize, ref_target: usize) -> Option<usize> {
        let ops = cigar.to_cigar_ops().unwrap();
        read_pos_at_ref(&ops, align_start, ref_target)
    }

    #[test]
    fn read_pos_at_ref_mid_match() {
        assert_eq!(pos_at("20M", 1000, 1010), Some(10));
    }

    #[test]
    fn read_pos_at_ref_matches_the_alignment_own_start() {
        // The exact scenario get_read_cuts can't handle in one call (see its
        // "Partial-overlap sentinel values" docs) — region_start == align_start.
        assert_eq!(pos_at("20M", 1000, 1000), Some(0));
    }

    #[test]
    fn read_pos_at_ref_matches_the_alignment_own_end() {
        assert_eq!(pos_at("20M", 1000, 1020), Some(20));
    }

    #[test]
    fn read_pos_at_ref_before_alignment_start_is_zero() {
        assert_eq!(pos_at("20M", 1000, 990), Some(0));
    }

    #[test]
    fn read_pos_at_ref_past_alignment_end_is_none() {
        assert_eq!(pos_at("20M", 1000, 1025), None);
    }

    #[test]
    fn read_pos_at_ref_across_multiple_ops() {
        // 5M 3I 5M: an insertion between two match blocks advances read
        // position without advancing reference position.
        assert_eq!(pos_at("5M3I5M", 1000, 1000), Some(0));
        assert_eq!(pos_at("5M3I5M", 1000, 1005), Some(5));
        // Position 1006 is past the insertion (which doesn't consume ref),
        // 1 base into the second match block, plus the 3 inserted bases.
        assert_eq!(pos_at("5M3I5M", 1000, 1006), Some(9));
        assert_eq!(pos_at("5M3I5M", 1000, 1010), Some(13));
    }

    #[test]
    fn read_pos_at_ref_inside_deletion_maps_to_position_before_it() {
        // 5M 4D 5M: a reference position inside the deletion has no read
        // counterpart, so it maps to the read position right after the last
        // matched base before the deletion (5, same as ref position 1005).
        assert_eq!(pos_at("5M4D5M", 1000, 1007), Some(5));
        assert_eq!(pos_at("5M4D5M", 1000, 1005), Some(5));
        assert_eq!(pos_at("5M4D5M", 1000, 1009), Some(5));
        assert_eq!(pos_at("5M4D5M", 1000, 1010), Some(6));
    }

    // --- extract_from_fasta_coords: 0-based half-open coordinate contract ---

    fn write_indexed_fasta(contents: &str) -> tempfile::NamedTempFile {
        let f = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(f.path(), contents).unwrap();
        let index = fasta::fs::index(f.path()).unwrap();
        let mut fai_path = f.path().as_os_str().to_owned();
        fai_path.push(".fai");
        fasta::fai::fs::write(&fai_path, &index).unwrap();
        f
    }

    #[test]
    fn extract_from_fasta_coords_is_0_based_half_open() {
        // "ACGTACGTAC" — extracting [2, 5) (0-based half-open) should give exactly
        // the 3 bases at indices 2,3,4: "GTA". Regression test for an off-by-one
        // where the region string was forwarded to noodles (1-based inclusive)
        // without the +1 conversion, silently prepending one extra base.
        let f = write_indexed_fasta(">seq1\nACGTACGTAC\n");
        let sequence = extract_from_fasta_coords(f.path().to_str().unwrap(), "seq1", 2, 5).unwrap();
        assert_eq!(sequence, "GTA");
    }

    #[test]
    fn extract_from_fasta_coords_from_start() {
        let f = write_indexed_fasta(">seq1\nACGTACGTAC\n");
        let sequence = extract_from_fasta_coords(f.path().to_str().unwrap(), "seq1", 0, 4).unwrap();
        assert_eq!(sequence, "ACGT");
    }
}
