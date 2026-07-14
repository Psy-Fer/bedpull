use anyhow::{Context, Result};
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::Cigar as SamCigar;
use noodles::sam::alignment::record::data::field::Tag as SamTag;

pub use crate::cigar::ToCigarOps;
use crate::paf::{PafIndexEntry, read_paf_record_from_reader};
use crate::utils::{
    ReadCuts, calculate_qscore, extract_from_fasta_coords_reader, get_read_cuts, revcomp,
};

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

    /// Minimum fraction (0.0-1.0) of the requested (region ± flanks) window a
    /// `partial` read must cover to be included. Reads/alignments below this
    /// threshold are discarded, similar to liftOver's `-minMatch`. Default: `0.0`
    /// (no filter — any overlap is accepted). Only meaningful when `partial` is
    /// `true`; non-partial reads always cover the full window by construction, so
    /// this is a no-op otherwise.
    pub min_partial_coverage: f64,

    /// Minimum mean Phred quality of the extracted subsequence.
    /// After the CIGAR walk, the quality string of the extracted slice is scored
    /// with [`crate::utils::calculate_qscore`]; reads below this
    /// threshold are discarded. Default: `0.0` (no filter). Applies to BAM/CRAM
    /// input regardless of output format — the read's quality scores are always
    /// present, whether or not `--fastq` is used for output.
    pub min_region_quality: f64,
}

impl Default for BamConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: false,
            partial: false,
            min_partial_coverage: 0.0,
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

/// A single PAF alignment after CIGAR-aware extraction.
///
/// The tuple fields are `(sequence, query_name, query_start, query_end, strand, haplotype)`.
///
/// - `sequence` — the extracted subsequence, always in reference (target) orientation.
/// - `query_name` — name of the query contig the sequence came from.
/// - `query_start` / `query_end` — 0-based half-open coordinates on the query contig.
/// - `strand` — alignment strand (`'+'` or `'-'`).
/// - `haplotype` — value of the `hp:i:` tag; `0` means absent (unphased).
pub type PafRead = (String, String, usize, usize, char, u8);

/// Resolve the final `(read_start, read_end, ref_start, ref_end)` for a read from
/// its raw [`ReadCuts`], applying partial-overlap fallback rules.
///
/// `get_read_cuts` fires on `ref_pos == region_start` and `ref_pos == region_end`. When
/// the alignment starts after `region_start`, the `region_start` trigger never fires;
/// `read_cuts.ref_start`/`read_cuts.read_start` end up holding the `region_end` position
/// instead (see [`get_read_cuts`][crate::utils::get_read_cuts] docs), and `read_end`/`ref_end`
/// stay `0`. This function undoes that mislabelling so callers get the reference span
/// actually covered by the returned slice, not the raw (and sometimes swapped) fields.
/// Returns `None` when the read doesn't satisfy the region and `config.partial` is `false`
/// (the read should be skipped).
fn resolve_cuts(
    read_cuts: &ReadCuts,
    config: &BamConfig,
    align_start: usize,
    align_end: usize,
    region_start: usize,
    seq_len: usize,
) -> Option<(usize, usize, usize, usize)> {
    if config.partial && align_start > region_start {
        let (read_end, ref_end) = if read_cuts.read_start > 0 {
            (read_cuts.read_start, read_cuts.ref_start)
        } else {
            (seq_len, align_end)
        };
        Some((0, read_end, align_start, ref_end))
    } else if read_cuts.read_end == 0 {
        if config.partial {
            Some((
                read_cuts.read_start,
                seq_len,
                read_cuts.ref_start,
                align_end,
            ))
        } else {
            None
        }
    } else {
        Some((
            read_cuts.read_start,
            read_cuts.read_end,
            read_cuts.ref_start,
            read_cuts.ref_end,
        ))
    }
}

/// Returns `false` if a read's actual reference coverage (`ref_start..ref_end`) falls short of
/// `config.min_partial_coverage` as a fraction of the requested (region ± flanks) window
/// (`desired_start..desired_end`). Always `true` when `min_partial_coverage` is `0.0` (the
/// default) or the window is empty.
fn passes_min_partial_coverage(
    config: &BamConfig,
    desired_start: usize,
    desired_end: usize,
    ref_start: usize,
    ref_end: usize,
) -> bool {
    if config.min_partial_coverage <= 0.0 {
        return true;
    }
    let desired_len = desired_end.saturating_sub(desired_start);
    if desired_len == 0 {
        return true;
    }
    let covered_len = ref_end.saturating_sub(ref_start);
    (covered_len as f64) >= config.min_partial_coverage * (desired_len as f64)
}

/// Extract subsequences from all BAM reads that overlap a given region.
///
/// Iterates `query` (a BAM region query result), applies the filters defined in
/// `config`, performs a CIGAR walk via [`crate::utils::get_read_cuts`]
/// to find the exact read-coordinate slice that corresponds to the reference region
/// (expanded by `lflank` / `rflank` bp on each side), and returns one [`BamRead`]
/// per passing read.
///
/// Returns `(reads, candidates_seen)`, where `candidates_seen` is the number of BAM
/// records the query yielded before any filter was applied — a region with
/// `candidates_seen == 0` had no overlapping alignments at all, whereas
/// `candidates_seen > 0 && reads.is_empty()` means every candidate was excluded by
/// `config` (mapq, flags, spanning/partial coverage, or region quality). Useful for
/// distinguishing "nothing here" from "filtered out" when reporting unmapped regions.
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
) -> Result<(Vec<BamRead>, usize)>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    let mut results: Vec<BamRead> = Vec::new();
    let mut candidates_seen: usize = 0;

    for result in query.records() {
        let record = result.context("failed to read BAM record")?;
        candidates_seen += 1;

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
        let i_seq: Vec<u8> = seq.iter().collect();
        let i_qual = record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&score| score + 33)
            .collect::<Vec<_>>();
        let quality_scores_str: String = String::from_utf8_lossy(&i_qual).into_owned();
        let cigar = record
            .cigar()
            .to_cigar_ops()
            .context("invalid CIGAR in BAM record")?;

        let region_start = region
            .interval()
            .start()
            .map(usize::from)
            .ok_or_else(|| anyhow::anyhow!("BED region has unbounded start"))?;
        let region_end = region
            .interval()
            .end()
            .map(usize::from)
            .ok_or_else(|| anyhow::anyhow!("BED region has unbounded end"))?;

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
        let Some((read_start, read_end, ref_start, ref_end)) = resolve_cuts(
            &read_cuts,
            config,
            align_start,
            align_end,
            region_start,
            i_seq.len(),
        ) else {
            continue;
        };
        if !passes_min_partial_coverage(config, desired_start, desired_end, ref_start, ref_end) {
            continue;
        }
        let subseq = i_seq[read_start..read_end].to_vec();
        let subqual: String = quality_scores_str[read_start..read_end].to_string();

        if config.min_region_quality > 0.0 && calculate_qscore(&subqual) < config.min_region_quality
        {
            continue;
        }

        let hap: u8 = record
            .data()
            .get(b"HP")
            .and_then(|v| v.ok())
            .and_then(|v| v.as_int())
            .map(|i| i as u8)
            .unwrap_or(0);

        results.push((name, subseq, subqual, ref_start, ref_end, hap));
    }

    Ok((results, candidates_seen))
}

/// Extract subsequences from all CRAM records that overlap a given region.
///
/// Works identically to [`get_bam_reads`] but consumes any iterator that yields
/// `io::Result<sam::alignment::RecordBuf>` — the item type produced by
/// `cram::io::IndexedReader::query`. Returns `(reads, candidates_seen)`; see
/// [`get_bam_reads`] for what `candidates_seen` means.
pub fn get_cram_reads(
    config: &BamConfig,
    query: impl Iterator<Item = std::io::Result<noodles::sam::alignment::RecordBuf>>,
    region: &Region,
    lflank: usize,
    rflank: usize,
) -> Result<(Vec<BamRead>, usize)> {
    use noodles::sam::alignment::record_buf::data::field::Value as RecordBufValue;

    let mut results: Vec<BamRead> = Vec::new();
    let mut candidates_seen: usize = 0;

    let region_start = region
        .interval()
        .start()
        .map(usize::from)
        .ok_or_else(|| anyhow::anyhow!("BED region has unbounded start"))?;
    let region_end = region
        .interval()
        .end()
        .map(usize::from)
        .ok_or_else(|| anyhow::anyhow!("BED region has unbounded end"))?;

    for result in query {
        let record = result.context("failed to read CRAM record")?;
        candidates_seen += 1;

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
                .ok_or_else(|| anyhow::anyhow!("CRAM record has no alignment start"))?,
        );
        let align_end = usize::from(
            record
                .alignment_end()
                .ok_or_else(|| anyhow::anyhow!("CRAM record has no alignment end"))?,
        );

        let name_bytes: &[u8] = record
            .name()
            .ok_or_else(|| anyhow::anyhow!("CRAM record has no name"))?
            .as_ref();
        let name = String::from_utf8(name_bytes.to_vec())
            .context("CRAM record name contains invalid UTF-8")?;

        let i_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let i_qual = record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&score| score + 33)
            .collect::<Vec<_>>();
        let quality_scores_str: String = String::from_utf8_lossy(&i_qual).into_owned();

        let cigar = (record.cigar() as &dyn SamCigar)
            .to_cigar_ops()
            .context("invalid CIGAR in CRAM record")?;

        let desired_start = region_start.saturating_sub(lflank);
        let desired_end = region_end + rflank;

        if (align_end < desired_start) || (align_start > desired_end) {
            continue;
        }

        let (eff_start, eff_end) = if lflank == 0 && rflank == 0 {
            (region_start, region_end)
        } else {
            (desired_start.max(align_start), desired_end.min(align_end))
        };

        let read_cuts = get_read_cuts(&cigar, align_start, eff_start, eff_end);
        let Some((read_start, read_end, ref_start, ref_end)) = resolve_cuts(
            &read_cuts,
            config,
            align_start,
            align_end,
            region_start,
            i_seq.len(),
        ) else {
            continue;
        };
        if !passes_min_partial_coverage(config, desired_start, desired_end, ref_start, ref_end) {
            continue;
        }

        let subseq = i_seq[read_start..read_end].to_vec();
        let subqual: String = quality_scores_str[read_start..read_end].to_string();

        if config.min_region_quality > 0.0 && calculate_qscore(&subqual) < config.min_region_quality
        {
            continue;
        }

        let hp_tag = SamTag::new(b'H', b'P');
        let hap: u8 = record
            .data()
            .get(&hp_tag)
            .and_then(|v| {
                if let RecordBufValue::Int32(i) = v {
                    Some(*i as u8)
                } else {
                    v.as_int().map(|i| i as u8)
                }
            })
            .unwrap_or(0);

        results.push((name, subseq, subqual, ref_start, ref_end, hap));
    }
    Ok((results, candidates_seen))
}

/// Extract subsequences from PAF alignment records that overlap a given region.
///
/// For each [`PafIndexEntry`] in `entries`, reads the full PAF record via `paf_reader`,
/// performs a CIGAR walk to compute the query-coordinate slice that corresponds to
/// the reference region (expanded by `lflank`/`rflank`), extracts the subsequence
/// from `fasta_reader`, and reverse-complements it for minus-strand alignments so the
/// output is always in target (reference) orientation.
///
/// `paf_reader` and `fasta_reader` are caller-opened and reused across calls (e.g. once
/// per BED region, or once for an entire run across many regions) rather than being
/// reopened per record — reopening a file per alignment is the dominant cost once BED
/// files grow past a handful of regions.
///
/// Records without a `cg:Z:` CIGAR tag or with invalid cut coordinates are skipped
/// with a diagnostic message to stderr.
#[allow(clippy::too_many_arguments)]
pub fn get_paf_reads<R>(
    paf_reader: &mut std::io::BufReader<std::fs::File>,
    fasta_reader: &mut noodles::fasta::io::IndexedReader<R>,
    entries: &[&PafIndexEntry],
    region_start: usize,
    region_end: usize,
    lflank: usize,
    rflank: usize,
    debug: bool,
) -> Result<Vec<PafRead>>
where
    R: std::io::BufRead + std::io::Seek,
{
    let mut results: Vec<PafRead> = Vec::new();

    for entry in entries {
        let paf_record = read_paf_record_from_reader(paf_reader, entry.offset)
            .with_context(|| format!("failed to read PAF record at offset {}", entry.offset))?;

        if paf_record.target_start > region_start {
            eprintln!("Warning: Alignment starts after region start, may be incomplete");
        }
        if paf_record.target_end < region_end {
            eprintln!("Warning: Alignment ends before region end, may be incomplete");
        }

        let eff_start = region_start
            .saturating_sub(lflank)
            .max(paf_record.target_start);
        let eff_end = (region_end + rflank).min(paf_record.target_end);

        let cigar_str = match &paf_record.cigar {
            Some(s) => s.clone(),
            None => {
                eprintln!("Warning: PAF record has no CIGAR (cg:Z: tag), skipping");
                continue;
            }
        };

        let cigar_ops = cigar_str
            .as_str()
            .to_cigar_ops()
            .context("invalid CIGAR string in PAF record")?;
        let cuts = get_read_cuts(&cigar_ops, paf_record.target_start, eff_start, eff_end);

        if cuts.read_start == 0 && cuts.read_end == 0 {
            eprintln!("Warning: No valid overlap found, skipping");
            continue;
        }
        if cuts.read_start >= cuts.read_end {
            eprintln!(
                "Warning: Invalid coordinates (start {} >= end {}), skipping",
                cuts.read_start, cuts.read_end
            );
            continue;
        }

        // For '+' strand: cuts are offsets from paf_record.query_start.
        // For '-' strand: cuts are offsets into the reverse-complemented query,
        // so we flip them relative to paf_record.query_end.
        let (query_start, query_end) = if paf_record.strand == '+' {
            (
                paf_record.query_start + cuts.read_start,
                paf_record.query_start + cuts.read_end,
            )
        } else {
            (
                paf_record.query_end.saturating_sub(cuts.read_end),
                paf_record.query_end.saturating_sub(cuts.read_start),
            )
        };

        if debug {
            eprintln!(
                "Query coords: {}:{}-{} (strand {})",
                paf_record.query_name, query_start, query_end, paf_record.strand
            );
        }

        let sequence = extract_from_fasta_coords_reader(
            fasta_reader,
            &paf_record.query_name,
            query_start,
            query_end,
        )
        .with_context(|| format!("failed to extract sequence for {}", paf_record.query_name))?;

        let sequence = if paf_record.strand == '-' {
            revcomp(&sequence)
        } else {
            sequence
        };

        let hap = paf_record.haplotype.unwrap_or(0);
        results.push((
            sequence,
            paf_record.query_name,
            query_start,
            query_end,
            paf_record.strand,
            hap,
        ));
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn bam_config_default_has_no_filters() {
        let c = BamConfig::default();
        assert_eq!(c.min_mapq, 0);
        assert!(!c.include_secondary);
        assert!(!c.include_supplementary);
        assert!(!c.partial);
        assert_eq!(c.min_partial_coverage, 0.0);
        assert_eq!(c.min_region_quality, 0.0);
    }

    #[test]
    fn bam_config_custom_values() {
        let c = BamConfig {
            min_mapq: 20,
            include_secondary: true,
            include_supplementary: false,
            partial: true,
            min_partial_coverage: 0.5,
            min_region_quality: 15.0,
        };
        assert_eq!(c.min_mapq, 20);
        assert!(c.include_secondary);
        assert!(!c.include_supplementary);
        assert!(c.partial);
        assert_eq!(c.min_partial_coverage, 0.5);
        assert_eq!(c.min_region_quality, 15.0);
    }

    fn partial_config() -> BamConfig {
        BamConfig {
            partial: true,
            ..BamConfig::default()
        }
    }

    #[test]
    fn resolve_cuts_full_span_returns_read_cuts_unchanged() {
        let read_cuts = ReadCuts {
            read_start: 10,
            read_end: 50,
            ref_start: 100,
            ref_end: 140,
        };
        let resolved = resolve_cuts(&read_cuts, &BamConfig::default(), 90, 200, 100, 300);
        assert_eq!(resolved, Some((10, 50, 100, 140)));
    }

    #[test]
    fn resolve_cuts_non_partial_incomplete_overlap_is_skipped() {
        // align_start == region_start (not left-partial) but read_end == 0 (never reached
        // region_end) and partial is off: the read should be skipped entirely.
        let read_cuts = ReadCuts {
            read_start: 10,
            read_end: 0,
            ref_start: 100,
            ref_end: 0,
        };
        let resolved = resolve_cuts(&read_cuts, &BamConfig::default(), 100, 130, 100, 300);
        assert_eq!(resolved, None);
    }

    #[test]
    fn resolve_cuts_right_partial_takes_rest_of_read() {
        // Alignment starts at region_start but ends before region_end: right-partial.
        let read_cuts = ReadCuts {
            read_start: 10,
            read_end: 0,
            ref_start: 100,
            ref_end: 0,
        };
        let resolved = resolve_cuts(&read_cuts, &partial_config(), 100, 130, 100, 300);
        // ref_end falls back to align_end (130) since the region end was never reached.
        assert_eq!(resolved, Some((10, 300, 100, 130)));
    }

    #[test]
    fn resolve_cuts_left_partial_starts_from_zero() {
        // Alignment starts after region_start: left-partial. get_read_cuts stores the
        // region_end position in read_start/ref_start (see its docs).
        let read_cuts = ReadCuts {
            read_start: 40,
            read_end: 0,
            ref_start: 140,
            ref_end: 0,
        };
        let resolved = resolve_cuts(&read_cuts, &partial_config(), 110, 200, 100, 300);
        assert_eq!(resolved, Some((0, 40, 110, 140)));
    }

    #[test]
    fn resolve_cuts_left_partial_never_reaches_region_end() {
        // Left-partial and the read also ends before region_end is reached at all.
        let read_cuts = ReadCuts {
            read_start: 0,
            read_end: 0,
            ref_start: 0,
            ref_end: 0,
        };
        let resolved = resolve_cuts(&read_cuts, &partial_config(), 110, 150, 100, 300);
        assert_eq!(resolved, Some((0, 300, 110, 150)));
    }

    // --- passes_min_partial_coverage ---

    fn min_coverage_config(min_partial_coverage: f64) -> BamConfig {
        BamConfig {
            min_partial_coverage,
            ..BamConfig::default()
        }
    }

    #[test]
    fn min_partial_coverage_zero_accepts_any_overlap() {
        // Default (0.0): even 1 base out of a 100bp window passes.
        let config = min_coverage_config(0.0);
        assert!(passes_min_partial_coverage(&config, 100, 200, 150, 151));
    }

    #[test]
    fn min_partial_coverage_full_span_always_passes() {
        let config = min_coverage_config(1.0);
        assert!(passes_min_partial_coverage(&config, 100, 200, 100, 200));
    }

    #[test]
    fn min_partial_coverage_below_threshold_is_rejected() {
        // Window is 100bp (100..200); covered span is 40bp (100..140) = 40% coverage.
        let config = min_coverage_config(0.5);
        assert!(!passes_min_partial_coverage(&config, 100, 200, 100, 140));
    }

    #[test]
    fn min_partial_coverage_at_exactly_threshold_passes() {
        // Window is 100bp; covered span is exactly 50bp = 50% coverage.
        let config = min_coverage_config(0.5);
        assert!(passes_min_partial_coverage(&config, 100, 200, 100, 150));
    }

    #[test]
    fn min_partial_coverage_above_threshold_passes() {
        // Window is 100bp; covered span is 90bp = 90% coverage.
        let config = min_coverage_config(0.5);
        assert!(passes_min_partial_coverage(&config, 100, 200, 105, 195));
    }

    #[test]
    fn min_partial_coverage_empty_window_always_passes() {
        let config = min_coverage_config(0.9);
        assert!(passes_min_partial_coverage(&config, 100, 100, 100, 100));
    }

    // --- get_paf_reads: reused-reader plumbing ---
    //
    // Exercises get_paf_reads with caller-opened, reused reader handles (the fix for
    // reopening the PAF/FASTA file per record). Calls it twice on the same open
    // readers, as extract_from_paf now does once per BED region, to confirm seeking
    // back into an already-open PAF file and re-querying an already-open indexed
    // FASTA reader both work correctly across repeated calls.

    #[test]
    fn get_paf_reads_reuses_readers_across_repeated_calls() {
        // 20bp query, straight 20M alignment to chr1:100-120 — no indels, so the
        // extracted length should exactly match the requested reference span.
        let paf_line = "q1\t20\t0\t20\t+\tchr1\t1000\t100\t120\t20\t20\t60\tcg:Z:20M\n";
        let mut paf_file = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut paf_file, paf_line.as_bytes()).unwrap();

        let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut fasta_file, b">q1\nACGTACGTACGTACGTACGT\n").unwrap();
        let fasta_index = noodles::fasta::fs::index(fasta_file.path()).unwrap();
        let mut fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
            .set_index(fasta_index)
            .build_from_path(fasta_file.path())
            .unwrap();

        let mut paf_reader = BufReader::new(File::open(paf_file.path()).unwrap());
        let entry = PafIndexEntry {
            offset: 0,
            target_start: 100,
            target_end: 120,
        };
        let entries = [&entry];

        // First call: chr1:105-110 (0-based), 5bp span.
        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            105,
            110,
            0,
            0,
            false,
        )
        .unwrap();
        assert_eq!(reads.len(), 1);
        let (sequence, query_name, query_start, query_end, strand, hap) = &reads[0];
        assert_eq!(sequence.len(), 5);
        assert_eq!(query_name, "q1");
        assert_eq!((*query_start, *query_end), (5, 10));
        assert_eq!(*strand, '+');
        assert_eq!(*hap, 0);

        // Second call on the SAME open readers, different region — proves the PAF
        // reader seeks correctly and the FASTA reader re-queries correctly instead of
        // returning stale state from the first call.
        let reads2 = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            110,
            118,
            0,
            0,
            false,
        )
        .unwrap();
        assert_eq!(reads2.len(), 1);
        let (sequence2, _, query_start2, query_end2, _, _) = &reads2[0];
        assert_eq!(sequence2.len(), 8);
        assert_eq!((*query_start2, *query_end2), (10, 18));
    }
}
