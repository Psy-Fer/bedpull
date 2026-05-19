use anyhow::{Context, Result};
use itertools::Itertools;
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;

pub use crate::cigar::ToCigarOps;
use crate::utils::{ReadCuts, calculate_qscore, get_read_cuts};

/// Filter and extraction settings for BAM reads.
///
/// This is the library equivalent of the CLI `Opts` fields that affect BAM extraction.
/// Construct with [`BamConfig::default`] for permissive defaults (no filters, spanning
/// reads only), then override individual fields as needed.
#[derive(Debug, Clone)]
pub struct BamConfig {
    /// Minimum mapping quality (MAPQ) a read must have to be included.
    /// Reads whose MAPQ is below this threshold are skipped. Default: `0` (no filter).
    pub min_mapq: u8,

    /// Whether to include secondary alignments (SAM flag `0x100`).
    /// Default: `false` (secondary alignments are skipped).
    pub include_secondary: bool,

    /// Whether to include supplementary alignments (SAM flag `0x800`).
    /// Default: `false` (supplementary alignments are skipped).
    pub include_supplementary: bool,

    /// Whether to include reads that only partially overlap the BED region.
    /// When `false` (the default), only reads whose alignment spans the entire
    /// region — from `region_start` to `region_end` — are returned. When `true`,
    /// reads that overlap any part of the region are included; the extracted
    /// subsequence is clipped to whatever portion the read covers.
    pub partial: bool,

    /// Minimum mean Phred quality of the extracted subsequence.
    /// After the CIGAR walk, the quality string of the extracted slice is scored
    /// with [`crate::utils::calculate_qscore`]; reads below this
    /// threshold are discarded. Default: `0.0` (no filter). Only meaningful when
    /// quality scores are present (BAM input with `--fastq`).
    pub min_region_quality: f64,
}

impl Default for BamConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: false,
            partial: false,
            min_region_quality: 0.0,
        }
    }
}

/// A single BAM read after CIGAR-aware extraction.
///
/// The tuple fields are `(name, sequence, quality_string, ref_start, ref_end, haplotype)`.
///
/// - `name` — read name from the BAM record.
/// - `sequence` — the extracted subsequence as raw bytes (ASCII nucleotides).
/// - `quality_string` — Phred+33 encoded quality string for the extracted slice.
/// - `ref_start` / `ref_end` — the reference coordinates actually covered by the
///   extracted slice (may differ from the requested BED coordinates when the alignment
///   does not perfectly span the region boundary).
/// - `haplotype` — value of the `HP` aux tag; `0` means the tag was absent (unphased).
pub type BamRead = (String, Vec<u8>, String, usize, usize, u8);

/// Extract subsequences from all BAM reads that overlap a given region.
///
/// Iterates `query` (a BAM region query result), applies the filters defined in
/// `config`, performs a CIGAR walk via [`crate::utils::get_read_cuts`]
/// to find the exact read-coordinate slice that corresponds to the reference region
/// (expanded by `lflank` / `rflank` bp on each side), and returns one [`BamRead`]
/// per passing read.
///
/// # Parameters
///
/// - `config` — filter and extraction settings; see [`BamConfig`].
/// - `query` — an open BAM region query produced by
///   `bam::io::IndexedReader::query`. Consumed by this function.
/// - `region` — the target region (0-based BED coordinates encoded as a noodles
///   `noodles::core::Region`; the interval bounds are read from
///   `region.interval()`).
/// - `lflank` — number of reference base pairs to extend the extraction window to
///   the left of `region_start`. Useful for capturing insertions that sit just
///   outside the annotated BED boundary.
/// - `rflank` — same as `lflank` but for the right side of `region_end`.
///
/// # Partial-overlap semantics
///
/// `get_read_cuts` signals a partial overlap through the values of `read_start` /
/// `read_end` (see that function's documentation). When `config.partial` is `true`,
/// this function interprets those sentinel values and slices the read from the
/// beginning or to the end of the sequence as appropriate.
pub fn get_bam_reads<R>(
    config: &BamConfig,
    query: bam::io::reader::Query<R>,
    region: &Region,
    lflank: usize,
    rflank: usize,
) -> Result<Vec<BamRead>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    let mut h0_subseq_vec: Vec<BamRead> = vec![]; // no hap assigned

    for result in query.records() {
        let record = result.context("failed to read BAM record")?;

        let map_quality = record.mapping_quality().map(u8::from).unwrap_or(255);
        if map_quality < config.min_mapq {
            continue;
        }

        let flags = record.flags();
        if flags.is_secondary() && !config.include_secondary {
            continue;
        }
        if flags.is_supplementary() && !config.include_supplementary {
            continue;
        }

        let align_start = usize::from(
            record
                .alignment_start()
                .ok_or_else(|| anyhow::anyhow!("BAM record has no alignment start"))?
                .context("invalid alignment start position")?,
        );
        let align_end = usize::from(
            record
                .alignment_end()
                .ok_or_else(|| anyhow::anyhow!("BAM record has no alignment end"))?
                .context("invalid alignment end position")?,
        );
        let name_bytes: &[u8] = record
            .name()
            .ok_or_else(|| anyhow::anyhow!("BAM record has no name"))?
            .as_ref();
        let name = String::from_utf8(name_bytes.to_vec())
            .context("BAM record name contains invalid UTF-8")?;
        let seq = record.sequence();
        let i_seq = seq.iter().collect_vec();
        let i_qual = record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&score| score + 33) // Adjust quality scores
            .collect::<Vec<_>>();

        // eprintln!("quality_scores adjusted: {:?}", i_qual);
        // now convert that to a String
        let quality_scores_str: String = String::from_utf8_lossy(&i_qual).into_owned();
        let cigar = record
            .cigar()
            .to_cigar_ops()
            .context("invalid CIGAR in BAM record")?;
        // Convert CIGAR operations to string by getting each kind, converting to a char, and going len|char and collecting
        // let cigar_string: String = cigar_to_string(&cigar);
        // get start and end position in read sequence coordinates using cigar string
        // take alignment start as 0 in read position. Work through cigar to get ref co-ord->read position for subsequence start
        // then continue through read until ref coord end-> read position for subsequence end
        // extract subsequence and chuck it into vector to be worked on
        // read_cuts = (start, end)

        // TODO: filter reads by mean read quality
        // let read_mean_qscore: f64 = calculate_qscore(&quality_scores_str);
        // if read_mean_qscore < config.min_read_quality {
        //     eprintln!("{} mean read quality of {} too low (min: {})", name, read_mean_qscore, config.min_read_quality);
        //     counter -= 1;
        //     continue;
        // }

        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());

        // Expand window by flanks; reads overlapping anywhere in the flanked window are included.
        let desired_start = region_start.saturating_sub(lflank);
        let desired_end = region_end + rflank;

        if (align_end < desired_start) || (align_start > desired_end) {
            continue;
        }

        // Clamp the desired window to this read's alignment span so get_read_cuts stays in bounds.
        let (eff_start, eff_end) = if lflank == 0 && rflank == 0 {
            (region_start, region_end)
        } else {
            (desired_start.max(align_start), desired_end.min(align_end))
        };

        let read_cuts: ReadCuts = get_read_cuts(&cigar, align_start, eff_start, eff_end);
        // get_read_cuts fires on ref_pos == region_start and ref_pos == region_end.
        // When align_start > region_start, the start trigger never fires; instead get_read_cuts
        // stores the region_end position in read_start (wrong field) and leaves read_end = 0.
        // Partial mode detects this and swaps the fields; see utils tests for exact semantics.
        let (read_start, read_end) = if config.partial && align_start > region_start {
            // Left-partial or contained: real start = 0 (beginning of read).
            // read_start holds region_end position if found; otherwise read is fully contained.
            let end = if read_cuts.read_start > 0 { read_cuts.read_start } else { i_seq.len() };
            (0, end)
        } else if read_cuts.read_end == 0 {
            // Right-partial (region_end not reached) or no overlap.
            if config.partial { (read_cuts.read_start, i_seq.len()) } else { continue }
        } else {
            (read_cuts.read_start, read_cuts.read_end)
        };
        let subseq = i_seq[read_start..read_end].to_vec();
        let subqual: String = quality_scores_str[read_start..read_end].to_string();

        if config.min_region_quality > 0.0 && calculate_qscore(&subqual) < config.min_region_quality {
            continue;
        }

        let hap: u8 = record
            .data()
            .get(b"HP")
            .and_then(|v| v.ok())
            .and_then(|v| v.as_int())
            .map(|i| i as u8)
            .unwrap_or(0);

        h0_subseq_vec.push((
            name,
            subseq.clone(),
            subqual,
            read_cuts.ref_start,
            read_cuts.ref_end,
            hap,
        ));
        // match hap {
        //     _ => eprintln!("multiple haplotypes detected. bedpull currently does not support phased data")
        // }

        // eprintln!("ref_id: {}", ref_id);
        // eprintln!("align_start: {}", align_start);
        // eprintln!("align_end: {}", align_end);
        // eprintln!("map_quality: {:?}", map_quality);
        // eprintln!("flags: {:?}", flags);
        // eprintln!("read span: {}", span);
        // eprintln!("subseq alignment span: {}", subseq_align_span);
        // eprintln!("subseq relative to reference: {}", subseq_align_span.saturating_sub(subseq.len() as isize));
        // // eprintln!("cigar: {:?}", cigar_string);
        // eprintln!("subseq len: {}", subseq.len());
        // eprintln!("subseq: {:?}", subseq_str);
        // // eprintln!("subqual: {:?}", subqual);
        // // eprintln!("all quality_scores_str: {:?}", quality_scores_str);
        // // eprintln!("data: {:?}", data);
        // eprintln!("HP tag: {:?}", hap);
        // eprintln!("subseq_vec: {:?}", subseq_vec);
    }
    // eprintln!("Number of reads in region: {}", counter);
    // eprintln!("Number of reads HAP1: {}", h1_subseq_vec.len());
    // eprintln!("Number of reads HAP2: {}", h2_subseq_vec.len());
    // eprintln!("Number of reads no HAP: {}", h0_subseq_vec.len());

    Ok(h0_subseq_vec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bam_config_default_has_no_filters() {
        let c = BamConfig::default();
        assert_eq!(c.min_mapq, 0);
        assert!(!c.include_secondary);
        assert!(!c.include_supplementary);
        assert!(!c.partial);
        assert_eq!(c.min_region_quality, 0.0);
    }

    #[test]
    fn bam_config_custom_values() {
        let c = BamConfig {
            min_mapq: 20,
            include_secondary: true,
            include_supplementary: false,
            partial: true,
            min_region_quality: 15.0,
        };
        assert_eq!(c.min_mapq, 20);
        assert!(c.include_secondary);
        assert!(!c.include_supplementary);
        assert!(c.partial);
        assert_eq!(c.min_region_quality, 15.0);
    }
}
