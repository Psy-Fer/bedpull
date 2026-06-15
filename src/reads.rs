use anyhow::{Context, Result};
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::Cigar as SamCigar;
use noodles::sam::alignment::record::data::field::Tag as SamTag;

pub use crate::cigar::ToCigarOps;
use crate::paf::{PafIndexEntry, read_paf_record_at_offset};
use crate::utils::{ReadCuts, calculate_qscore, extract_from_fasta_coords, get_read_cuts, revcomp};

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
    let mut results: Vec<BamRead> = Vec::new();

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
        // get_read_cuts fires on ref_pos == region_start and ref_pos == region_end.
        // When align_start > region_start, the start trigger never fires; instead get_read_cuts
        // stores the region_end position in read_start (wrong field) and leaves read_end = 0.
        // Partial mode detects this and swaps the fields; see utils tests for exact semantics.
        let (read_start, read_end) = if config.partial && align_start > region_start {
            let end = if read_cuts.read_start > 0 {
                read_cuts.read_start
            } else {
                i_seq.len()
            };
            (0, end)
        } else if read_cuts.read_end == 0 {
            if config.partial {
                (read_cuts.read_start, i_seq.len())
            } else {
                continue;
            }
        } else {
            (read_cuts.read_start, read_cuts.read_end)
        };
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

        results.push((
            name,
            subseq,
            subqual,
            read_cuts.ref_start,
            read_cuts.ref_end,
            hap,
        ));
    }

    Ok(results)
}

/// Extract subsequences from all CRAM records that overlap a given region.
///
/// Works identically to [`get_bam_reads`] but consumes any iterator that yields
/// `io::Result<sam::alignment::RecordBuf>` — the item type produced by
/// `cram::io::IndexedReader::query`.
pub fn get_cram_reads(
    config: &BamConfig,
    query: impl Iterator<Item = std::io::Result<noodles::sam::alignment::RecordBuf>>,
    region: &Region,
    lflank: usize,
    rflank: usize,
) -> Result<Vec<BamRead>> {
    use noodles::sam::alignment::record_buf::data::field::Value as RecordBufValue;

    let mut results: Vec<BamRead> = Vec::new();

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
        let (read_start, read_end) = if config.partial && align_start > region_start {
            let end = if read_cuts.read_start > 0 {
                read_cuts.read_start
            } else {
                i_seq.len()
            };
            (0, end)
        } else if read_cuts.read_end == 0 {
            if config.partial {
                (read_cuts.read_start, i_seq.len())
            } else {
                continue;
            }
        } else {
            (read_cuts.read_start, read_cuts.read_end)
        };

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

        results.push((
            name,
            subseq,
            subqual,
            read_cuts.ref_start,
            read_cuts.ref_end,
            hap,
        ));
    }
    Ok(results)
}

/// Extract subsequences from PAF alignment records that overlap a given region.
///
/// For each [`PafIndexEntry`] in `entries`, reads the full PAF record from disk,
/// performs a CIGAR walk to compute the query-coordinate slice that corresponds to
/// the reference region (expanded by `lflank`/`rflank`), extracts the subsequence
/// from the query FASTA at `query_ref`, and reverse-complements it for minus-strand
/// alignments so the output is always in target (reference) orientation.
///
/// Records without a `cg:Z:` CIGAR tag or with invalid cut coordinates are skipped
/// with a diagnostic message to stderr.
pub fn get_paf_reads(
    paf_path: &str,
    query_ref: &str,
    entries: &[&PafIndexEntry],
    region_start: usize,
    region_end: usize,
    lflank: usize,
    rflank: usize,
) -> Result<Vec<PafRead>> {
    let mut results: Vec<PafRead> = Vec::new();

    for entry in entries {
        let paf_record = read_paf_record_at_offset(paf_path, entry.offset)
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

        eprintln!(
            "Query coords: {}:{}-{} (strand {})",
            paf_record.query_name, query_start, query_end, paf_record.strand
        );

        let sequence =
            extract_from_fasta_coords(query_ref, &paf_record.query_name, query_start, query_end)
                .with_context(|| {
                    format!("failed to extract sequence for {}", paf_record.query_name)
                })?;

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
