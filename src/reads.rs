use anyhow::{Context, Result};
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::Cigar as SamCigar;
use noodles::sam::alignment::record::data::field::Tag as SamTag;

pub use crate::cigar::ToCigarOps;
use crate::paf::{PafIndexEntry, PafRecord, read_paf_record_from_reader};
use crate::utils::{
    ReadCuts, calculate_qscore, extract_from_fasta_coords_reader, get_read_cuts, read_pos_at_ref,
    revcomp,
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
/// - `ref_start` / `ref_end` — the reference coordinates (0-based, matching BED)
///   actually covered by the extracted slice (may differ from the requested BED
///   coordinates when the alignment does not perfectly span the region boundary).
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

/// Cross-record stitching settings for PAF extraction — see [`get_paf_reads`].
///
/// A single PAF record only ever covers a contiguous block of the alignment;
/// a large structural variant (most often a big novel insertion with no
/// homologous target sequence) frequently causes the aligner to emit *two or
/// more separate, chained records* instead of one record with a large indel
/// operation. Without stitching, a BED window spanning that boundary is only
/// ever partially extracted from whichever single record happens to overlap
/// it. Stitching looks for a chain of records sharing a query contig and
/// strand, contiguous in target space, that together span the window, and
/// extracts a single sequence across the whole chain in one slice (see
/// [`get_paf_reads`] docs for how the query-side "gap" between chain members
/// is handled).
#[derive(Debug, Clone, Copy)]
pub struct StitchConfig {
    /// When `false` (the default), stitching is never attempted; a window not
    /// fully covered by any single PAF record is only ever partially
    /// extracted (or produces multiple separate partial hits, one per
    /// overlapping record).
    pub enabled: bool,
    /// Maximum allowed *target*-space gap (bp) between consecutive records in
    /// a candidate chain for them to still be considered part of the same
    /// split alignment. This bounds the target-side gap only — the
    /// reconstructed query-side splice (e.g. a large novel insertion with no
    /// target homolog) can be much larger and is not separately bounded.
    pub max_gap: usize,
}

impl Default for StitchConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            max_gap: 10_000,
        }
    }
}

/// Find a chain of 2+ records sharing a query contig and strand, contiguous
/// in target space (consecutive gaps no larger than `max_gap`), whose
/// combined target span `[first.target_start, last.target_end)` covers
/// `[desired_start, desired_end)`. Returns the original indices into
/// `records` belonging to the winning chain, sorted by `target_start`, or
/// `None` if no such chain exists. When more than one candidate chain
/// qualifies, the first one found wins, where "first" means the
/// `(query_name, strand)` group whose key first appears in `records` — a
/// deterministic order tracked explicitly, not `HashMap` iteration order
/// (which is randomized per-process and would otherwise make the result
/// vary run-to-run on identical input). Good enough for the common case of
/// one real split alignment per region; pathological inputs with multiple
/// overlapping candidate chains aren't specially disambiguated beyond that.
fn find_stitch_chain(
    records: &[PafRecord],
    desired_start: usize,
    desired_end: usize,
    max_gap: usize,
) -> Option<Vec<usize>> {
    // HashMap iteration order is randomized per-process, so groups must be
    // walked in a separately-tracked, deterministic order (first appearance
    // in `records`) rather than however the map happens to enumerate them —
    // otherwise which chain wins when more than one qualifies would vary
    // run-to-run on the exact same input.
    let mut groups: std::collections::HashMap<(&str, char), Vec<usize>> =
        std::collections::HashMap::new();
    let mut key_order: Vec<(&str, char)> = Vec::new();
    for (i, r) in records.iter().enumerate() {
        let key = (r.query_name.as_str(), r.strand);
        if !groups.contains_key(&key) {
            key_order.push(key);
        }
        groups.entry(key).or_default().push(i);
    }

    for key in &key_order {
        let idxs = groups
            .get_mut(key)
            .expect("key_order only holds keys already inserted into groups");
        idxs.sort_by_key(|&i| records[i].target_start);

        let mut run_start = 0usize;
        for i in 1..idxs.len() {
            let prev = &records[idxs[i - 1]];
            let cur = &records[idxs[i]];
            let gap = cur.target_start.saturating_sub(prev.target_end);
            if gap > max_gap {
                if let Some(span) =
                    winning_run(records, &idxs[run_start..i], desired_start, desired_end)
                {
                    return Some(span);
                }
                run_start = i;
            }
        }
        if let Some(span) = winning_run(records, &idxs[run_start..], desired_start, desired_end) {
            return Some(span);
        }
    }
    None
}

/// A contiguous run qualifies if it has 2+ members (a single record is
/// already handled by the normal per-record path) and its combined target
/// span covers `[desired_start, desired_end)`.
fn winning_run(
    records: &[PafRecord],
    run_idxs: &[usize],
    desired_start: usize,
    desired_end: usize,
) -> Option<Vec<usize>> {
    if run_idxs.len() < 2 {
        return None;
    }
    let first = &records[run_idxs[0]];
    let last = &records[run_idxs[run_idxs.len() - 1]];
    if first.target_start <= desired_start && last.target_end >= desired_end {
        Some(run_idxs.to_vec())
    } else {
        None
    }
}

/// Extract a single sequence spanning a whole stitched chain.
///
/// Only the *first* record's CIGAR (to locate `desired_start`) and the
/// *last* record's CIGAR (to locate `desired_end`) are walked — the chain
/// shares one query contig, so once both boundary query coordinates are
/// known, the result is a single contiguous slice of that contig. Any
/// members between the first and last, and any "gap" between chain members
/// with no aligned counterpart at all, fall automatically inside that slice
/// — which is the entire point: that gap is exactly the sequence a large
/// indel/SV split the alignment around, spliced back in from the raw query
/// rather than reconstructed from an alignment that was never made for it.
fn build_stitched_read<R>(
    records: &[PafRecord],
    chain: &[usize],
    desired_start: usize,
    desired_end: usize,
    fasta_reader: &mut noodles::fasta::io::IndexedReader<R>,
    debug: bool,
) -> Result<Option<PafRead>>
where
    R: std::io::BufRead + std::io::Seek,
{
    let first = &records[chain[0]];
    let last = &records[chain[chain.len() - 1]];

    let (Some(first_cigar_str), Some(last_cigar_str)) = (&first.cigar, &last.cigar) else {
        return Ok(None);
    };

    // A malformed CIGAR in one candidate chain is this chain's problem, not the
    // whole run's — decline to stitch rather than aborting every other region.
    let (Ok(first_ops), Ok(last_ops)) = (
        first_cigar_str.as_str().to_cigar_ops(),
        last_cigar_str.as_str().to_cigar_ops(),
    ) else {
        eprintln!("Warning: invalid CIGAR string in a stitch chain member, skipping stitch");
        return Ok(None);
    };

    // Deliberately not get_read_cuts: the two window edges live on *different*
    // chained records here (desired_start on the first record's CIGAR,
    // desired_end on the last's), so each is resolved independently with
    // read_pos_at_ref. get_read_cuts models a single alignment's two boundaries
    // and doesn't fit a boundary-per-record split.
    let (Some(read_pos_first), Some(read_pos_last)) = (
        read_pos_at_ref(&first_ops, first.target_start, desired_start),
        read_pos_at_ref(&last_ops, last.target_start, desired_end),
    ) else {
        eprintln!("Warning: stitched chain produced no valid overlap, skipping");
        return Ok(None);
    };

    let strand = first.strand;
    let (query_start, query_end) = if strand == '+' {
        (
            first.query_start + read_pos_first,
            last.query_start + read_pos_last,
        )
    } else {
        (
            last.query_end.saturating_sub(read_pos_last),
            first.query_end.saturating_sub(read_pos_first),
        )
    };

    // query_start == query_end is a valid (if unusual) zero-length result — the
    // whole stitched span collapsed to a single query position, meaning the
    // window has no corresponding query bases at all. Only a genuinely
    // inverted span (query_start > query_end) is an error.
    if query_start > query_end {
        eprintln!(
            "Warning: stitched chain produced invalid coordinates (start {} > end {}), skipping",
            query_start, query_end
        );
        return Ok(None);
    }

    if debug {
        eprintln!(
            "Stitched {} chained record(s): {}:{}-{} (strand {})",
            chain.len(),
            first.query_name,
            query_start,
            query_end,
            strand
        );
    }

    let sequence = match extract_from_fasta_coords_reader(
        fasta_reader,
        &first.query_name,
        query_start,
        query_end,
    ) {
        Ok(seq) => seq,
        Err(e) => {
            eprintln!(
                "Warning: failed to extract stitched sequence for {}: {:#}, skipping",
                first.query_name, e
            );
            return Ok(None);
        }
    };
    let sequence = if strand == '-' {
        revcomp(&sequence)
    } else {
        sequence
    };

    let hap = first.haplotype.unwrap_or(0);
    Ok(Some((
        sequence,
        first.query_name.clone(),
        query_start,
        query_end,
        strand,
        hap,
    )))
}

/// Resolve the final `(read_start, read_end, ref_start, ref_end)` for a read from
/// its raw [`ReadCuts`].
///
/// [`get_read_cuts`][crate::utils::get_read_cuts] now clamps the fire-on boundaries
/// to the alignment span itself, so `read_start`/`read_end` are already correct read
/// offsets and `ref_start`/`ref_end` already report the covered sub-span — no
/// sentinel un-swapping is needed. This function is left as the one place the
/// spanning/partial *policy* is applied:
///
/// - `read_end == 0` means the alignment never reached the window at all → skip.
/// - `spans` is `true` when the alignment covers the whole requested region
///   (`align_start <= region_start && align_end >= region_end`). In non-partial
///   mode a read that does not fully span the region is skipped; in partial mode
///   any real overlap is accepted and the covered sub-span is returned.
fn resolve_cuts(
    read_cuts: &ReadCuts,
    config: &BamConfig,
    spans: bool,
) -> Option<(usize, usize, usize, usize)> {
    if read_cuts.read_end == 0 {
        return None;
    }
    if !config.partial && !spans {
        return None;
    }
    Some((
        read_cuts.read_start,
        read_cuts.read_end,
        read_cuts.ref_start,
        read_cuts.ref_end,
    ))
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

        // noodles alignment_end() is the inclusive last aligned position; get_read_cuts
        // wants the exclusive one-past-end boundary, in the same frame as region_end.
        let align_end_excl = align_end + 1;

        // get_read_cuts clamps its fire-on boundaries to [align_start, align_end_excl]
        // internally, so the desired (flank-expanded) window is passed straight through.
        let read_cuts: ReadCuts = get_read_cuts(
            &cigar,
            align_start,
            align_end_excl,
            desired_start,
            desired_end,
        );
        // Spanning is judged against the requested region (not the flank window): a read
        // must cover region_start..region_end to count as full-length; flanks are captured
        // opportunistically when the alignment reaches into them.
        let spans = align_start <= region_start && align_end_excl >= region_end;
        let Some((read_start, read_end, ref_start, ref_end)) =
            resolve_cuts(&read_cuts, config, spans)
        else {
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

        // The CIGAR walk runs in the same 1-based frame as noodles' alignment_start
        // (and read_bed's +1-shifted region bounds), so ref_start/ref_end come back
        // 1-based. Callers (and the |missing_left/right header suffix) work in 0-based
        // BED coordinates, so normalise here — otherwise every fully-spanning read is
        // mislabelled `missing_left=1bp`. Lengths are unchanged, so the coverage check
        // above is unaffected.
        let ref_start = ref_start.saturating_sub(1);
        let ref_end = ref_end.saturating_sub(1);

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

        // noodles alignment_end() is inclusive; get_read_cuts wants one-past-end.
        let align_end_excl = align_end + 1;

        let read_cuts = get_read_cuts(
            &cigar,
            align_start,
            align_end_excl,
            desired_start,
            desired_end,
        );
        let spans = align_start <= region_start && align_end_excl >= region_end;
        let Some((read_start, read_end, ref_start, ref_end)) =
            resolve_cuts(&read_cuts, config, spans)
        else {
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

        // Normalise ref_start/ref_end from the 1-based walk frame to 0-based BED
        // coordinates, as in get_bam_reads (keeps the |missing_left/right suffix honest).
        let ref_start = ref_start.saturating_sub(1);
        let ref_end = ref_end.saturating_sub(1);

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
///
/// When `stitch.enabled` and no single record fully covers the requested
/// region (± flanks), a chain of same-contig, same-strand, target-contiguous
/// records that collectively do is searched for (see [`StitchConfig`]); if
/// found, one stitched [`PafRead`] replaces the fragmentary per-record output
/// those chain members would otherwise have produced. Any other overlapping
/// records not part of a winning chain are still extracted individually as
/// before.
#[allow(clippy::too_many_arguments)]
pub fn get_paf_reads<R>(
    paf_reader: &mut std::io::BufReader<std::fs::File>,
    fasta_reader: &mut noodles::fasta::io::IndexedReader<R>,
    entries: &[&PafIndexEntry],
    region_start: usize,
    region_end: usize,
    lflank: usize,
    rflank: usize,
    stitch: StitchConfig,
    debug: bool,
) -> Result<Vec<PafRead>>
where
    R: std::io::BufRead + std::io::Seek,
{
    let mut results: Vec<PafRead> = Vec::new();

    // Read every overlapping record once up front — reused for both the normal
    // per-record extraction below and (if enabled) cross-record stitching.
    let mut records: Vec<PafRecord> = Vec::with_capacity(entries.len());
    for entry in entries {
        let paf_record = read_paf_record_from_reader(paf_reader, entry.offset)
            .with_context(|| format!("failed to read PAF record at offset {}", entry.offset))?;
        records.push(paf_record);
    }

    let desired_start = region_start.saturating_sub(lflank);
    let desired_end = region_end + rflank;

    let mut stitched_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();
    if stitch.enabled {
        let already_covered = records
            .iter()
            .any(|r| r.target_start <= desired_start && r.target_end >= desired_end);
        if !already_covered
            && let Some(chain) =
                find_stitch_chain(&records, desired_start, desired_end, stitch.max_gap)
            && let Some(stitched) = build_stitched_read(
                &records,
                &chain,
                desired_start,
                desired_end,
                fasta_reader,
                debug,
            )?
        {
            results.push(stitched);
            stitched_indices.extend(chain.iter().copied());
        }
    }

    for (i, paf_record) in records.iter().enumerate() {
        if stitched_indices.contains(&i) {
            continue;
        }

        if paf_record.target_start > region_start {
            eprintln!("Warning: Alignment starts after region start, may be incomplete");
        }
        if paf_record.target_end < region_end {
            eprintln!("Warning: Alignment ends before region end, may be incomplete");
        }

        let eff_start = desired_start.max(paf_record.target_start);
        let eff_end = desired_end.min(paf_record.target_end);

        let cigar_str = match &paf_record.cigar {
            Some(s) => s.clone(),
            None => {
                eprintln!("Warning: PAF record has no CIGAR (cg:Z: tag), skipping");
                continue;
            }
        };

        // A malformed CIGAR in one record is this record's problem, not the
        // whole run's — skip it and keep processing the rest.
        let cigar_ops = match cigar_str.as_str().to_cigar_ops() {
            Ok(ops) => ops,
            Err(e) => {
                eprintln!(
                    "Warning: invalid CIGAR string in PAF record: {:#}, skipping",
                    e
                );
                continue;
            }
        };

        // read_pos_at_ref resolves each window edge to a query offset independently,
        // which is what the PAF path wants: it handles the strand flip and zero-length
        // (fully-inside-a-deletion) cases below by comparing the two offsets directly,
        // and eff_start routinely equals paf_record.target_start exactly (any record
        // whose window extends past its own left edge), which single-boundary mapping
        // handles cleanly.
        let (Some(read_start), Some(read_end)) = (
            read_pos_at_ref(&cigar_ops, paf_record.target_start, eff_start),
            read_pos_at_ref(&cigar_ops, paf_record.target_start, eff_end),
        ) else {
            eprintln!("Warning: No valid overlap found, skipping");
            continue;
        };

        // read_start == read_end is a valid (if unusual) zero-length result — the
        // whole requested span falls inside a deletion, so it has no
        // corresponding query bases at all. Only a genuinely inverted span
        // (read_start > read_end) is an error.
        if read_start > read_end {
            eprintln!(
                "Warning: Invalid coordinates (start {} > end {}), skipping",
                read_start, read_end
            );
            continue;
        }

        // For '+' strand: offsets are from paf_record.query_start.
        // For '-' strand: offsets are into the reverse-complemented query,
        // so they're flipped relative to paf_record.query_end.
        let (query_start, query_end) = if paf_record.strand == '+' {
            (
                paf_record.query_start + read_start,
                paf_record.query_start + read_end,
            )
        } else {
            (
                paf_record.query_end.saturating_sub(read_end),
                paf_record.query_end.saturating_sub(read_start),
            )
        };

        if debug {
            eprintln!(
                "Query coords: {}:{}-{} (strand {})",
                paf_record.query_name, query_start, query_end, paf_record.strand
            );
        }

        // A single bad/mismatched contig name (e.g. a PAF built against a
        // different FASTA than --query_ref) shouldn't abort every other
        // region's extraction — skip just this record instead.
        let sequence = match extract_from_fasta_coords_reader(
            fasta_reader,
            &paf_record.query_name,
            query_start,
            query_end,
        ) {
            Ok(seq) => seq,
            Err(e) => {
                eprintln!(
                    "Warning: failed to extract sequence for {}: {:#}, skipping",
                    paf_record.query_name, e
                );
                continue;
            }
        };

        let sequence = if paf_record.strand == '-' {
            revcomp(&sequence)
        } else {
            sequence
        };

        let hap = paf_record.haplotype.unwrap_or(0);
        results.push((
            sequence,
            paf_record.query_name.clone(),
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

    fn cuts_of(read_start: usize, read_end: usize, ref_start: usize, ref_end: usize) -> ReadCuts {
        ReadCuts {
            read_start,
            read_end,
            ref_start,
            ref_end,
            softclip_lead_start: read_start,
            softclip_trail_end: read_end,
        }
    }

    #[test]
    fn resolve_cuts_spanning_returns_read_cuts_unchanged() {
        // get_read_cuts already produced correct, clamped offsets — a spanning read
        // passes straight through in either mode.
        let read_cuts = cuts_of(10, 50, 100, 140);
        let resolved = resolve_cuts(&read_cuts, &BamConfig::default(), true);
        assert_eq!(resolved, Some((10, 50, 100, 140)));
    }

    #[test]
    fn resolve_cuts_no_overlap_is_skipped_even_in_partial_mode() {
        // read_end == 0 is the "alignment never reached the window" sentinel: nothing to
        // extract, so it's skipped regardless of partial mode.
        let read_cuts = cuts_of(10, 0, 100, 0);
        assert_eq!(resolve_cuts(&read_cuts, &BamConfig::default(), false), None);
        assert_eq!(resolve_cuts(&read_cuts, &partial_config(), false), None);
    }

    #[test]
    fn resolve_cuts_non_partial_non_spanning_is_skipped() {
        // A real overlap (read_end != 0) that doesn't fully span the region is dropped in
        // non-partial mode.
        let read_cuts = cuts_of(0, 40, 110, 150);
        assert_eq!(resolve_cuts(&read_cuts, &BamConfig::default(), false), None);
    }

    #[test]
    fn resolve_cuts_partial_non_spanning_returns_covered_subspan() {
        // Same overlap, partial mode: the covered sub-span (already clamped by
        // get_read_cuts) is returned as-is.
        let read_cuts = cuts_of(0, 40, 110, 150);
        assert_eq!(
            resolve_cuts(&read_cuts, &partial_config(), false),
            Some((0, 40, 110, 150))
        );
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
            StitchConfig::default(),
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
            StitchConfig::default(),
            false,
        )
        .unwrap();
        assert_eq!(reads2.len(), 1);
        let (sequence2, _, query_start2, query_end2, _, _) = &reads2[0];
        assert_eq!(sequence2.len(), 8);
        assert_eq!((*query_start2, *query_end2), (10, 18));
    }

    // --- cross-record stitching ---
    //
    // Simulates the real failure mode this feature targets: a large novel
    // insertion (no target homolog) causes the aligner to emit two separate,
    // target-adjacent records instead of one record with a big `I` op. A
    // window straddling the boundary is only ever partially covered by
    // either record alone.

    /// Query contig: [0..20) left arm, [20..50) unaligned 30bp "insertion",
    /// [50..70) right arm. Two records: A covers target [1000,1020) against
    /// query [0,20), B covers target [1020,1040) against query [50,70) —
    /// contiguous in target space, with the 30bp query gap between them
    /// standing in for the insertion neither record aligned. CIGARs are
    /// written as four 5M ops rather than one 20M so that a region boundary
    /// landing mid-record (as in the "disabled" test below) doesn't also
    /// trip the unrelated pre-existing `get_read_cuts` quirk where two
    /// boundaries resolving inside the *same* op get misattributed (see
    /// `read_pos_at_ref`'s docs) — orthogonal to what these tests exercise.
    fn stitch_test_paf_plus() -> (&'static str, &'static str) {
        (
            "q1\t70\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:5M5M5M5M\n\
             q1\t70\t50\t70\t+\tchr1\t2000\t1020\t1040\t20\t20\t60\tcg:Z:5M5M5M5M\n",
            ">q1\nCCCCCCCCCCGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAACCCCCCCCCC\n",
        )
    }

    fn setup_stitch_reader(
        paf_contents: &str,
        fasta_contents: &str,
    ) -> (
        tempfile::NamedTempFile,
        BufReader<File>,
        noodles::fasta::io::IndexedReader<noodles::fasta::io::BufReader<File>>,
    ) {
        let mut paf_file = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut paf_file, paf_contents.as_bytes()).unwrap();
        let paf_reader = BufReader::new(File::open(paf_file.path()).unwrap());

        let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut fasta_file, fasta_contents.as_bytes()).unwrap();
        let fasta_index = noodles::fasta::fs::index(fasta_file.path()).unwrap();
        let fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
            .set_index(fasta_index)
            .build_from_path(fasta_file.path())
            .unwrap();

        (paf_file, paf_reader, fasta_reader)
    }

    #[test]
    fn stitch_disabled_only_produces_partial_fragments() {
        let (paf_contents, fasta_contents) = stitch_test_paf_plus();
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_a = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_b = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1020,
            target_end: 1040,
        };
        let entries = [&entry_a, &entry_b];

        // Window [1010, 1030) straddles the boundary; neither record alone covers it.
        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1010,
            1030,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        // Two partial fragments — record A's [10,20) ("GGGGGGGGGG") and
        // record B's [50,60) ("AAAAAAAAAA") — not the correct 50bp stitched
        // result spanning the insertion between them. Record B's clamped
        // window start lands exactly on its own align_start
        // (region_start == align_start), which per-record extraction now
        // resolves correctly via read_pos_at_ref instead of misattributing
        // it as invalid (see read_pos_at_ref's docs) — but a correct partial
        // fragment is still not the same as the correct whole-window answer,
        // which is exactly what stitching is for.
        assert_eq!(reads.len(), 2);
        let mut by_seq: std::collections::HashMap<&str, (usize, usize)> =
            std::collections::HashMap::new();
        for (sequence, _, query_start, query_end, ..) in &reads {
            by_seq.insert(sequence.as_str(), (*query_start, *query_end));
        }
        assert_eq!(by_seq.get("GGGGGGGGGG"), Some(&(10, 20)));
        assert_eq!(by_seq.get("AAAAAAAAAA"), Some(&(50, 60)));
    }

    #[test]
    fn stitch_enabled_bridges_split_insertion_plus_strand() {
        let (paf_contents, fasta_contents) = stitch_test_paf_plus();
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_a = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_b = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1020,
            target_end: 1040,
        };
        let entries = [&entry_a, &entry_b];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1010,
            1030,
            0,
            0,
            StitchConfig {
                enabled: true,
                max_gap: 100,
            },
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        let (sequence, query_name, query_start, query_end, strand, _hap) = &reads[0];
        assert_eq!(query_name, "q1");
        assert_eq!(*strand, '+');
        assert_eq!((*query_start, *query_end), (10, 60));
        let expected = "GGGGGGGGGG".to_string() + &"T".repeat(30) + "AAAAAAAAAA";
        assert_eq!(*sequence, expected);
    }

    #[test]
    fn stitch_enabled_bridges_split_insertion_minus_strand() {
        // Same shape, but strand='-': target still ascends A -> B, forward-query
        // descends (B's arm is the *lower* forward coordinates, A's arm the
        // *higher* ones) — the reverse of the plus-strand case.
        let paf_contents = "q1\t70\t50\t70\t-\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:5M5M5M5M\n\
                             q1\t70\t0\t20\t-\tchr1\t2000\t1020\t1040\t20\t20\t60\tcg:Z:5M5M5M5M\n";
        let fasta_contents =
            ">q1\nCCCCCCCCCCGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAACCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_a = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_b = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1020,
            target_end: 1040,
        };
        let entries = [&entry_a, &entry_b];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1010,
            1030,
            0,
            0,
            StitchConfig {
                enabled: true,
                max_gap: 100,
            },
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        let (sequence, _query_name, _query_start, _query_end, strand, _hap) = &reads[0];
        assert_eq!(*strand, '-');
        // Forward-strand slice would be the same 50bp pattern as the plus-strand
        // case; the stitched result is that slice's reverse complement.
        let forward = "GGGGGGGGGG".to_string() + &"T".repeat(30) + "AAAAAAAAAA";
        let expected = revcomp(&forward);
        assert_eq!(*sequence, expected);
    }

    #[test]
    fn stitch_no_chain_when_window_exceeds_chain_coverage() {
        let (paf_contents, fasta_contents) = stitch_test_paf_plus();
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_a = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_b = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1020,
            target_end: 1040,
        };
        let entries = [&entry_a, &entry_b];

        // desired_end=1050 exceeds B's own target_end (1040), so the [A, B]
        // chain doesn't cover the full requested window — no stitch should
        // be emitted, confirming stitching cleanly declines rather than
        // fabricating a wrong answer when a chain is found but insufficient.
        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1010,
            1050,
            0,
            0,
            StitchConfig {
                enabled: true,
                max_gap: 100,
            },
            false,
        )
        .unwrap();

        // Falls back to the normal per-record loop — no stitch, so no 50bp
        // result spanning the insertion. Both records still contribute their
        // own ordinary partial fragments: A's 10bp [10,20) and B's full 20bp
        // [50,70) (desired_end=1050 clamps to B's own target_end=1040, which
        // maps to B's own full query span since it covers the whole record).
        assert!(reads.iter().all(|(sequence, ..)| sequence.len() < 30));
        assert_eq!(reads.len(), 2);
        let mut by_seq: std::collections::HashMap<&str, (usize, usize)> =
            std::collections::HashMap::new();
        for (sequence, _, query_start, query_end, ..) in &reads {
            by_seq.insert(sequence.as_str(), (*query_start, *query_end));
        }
        assert_eq!(by_seq.get("GGGGGGGGGG"), Some(&(10, 20)));
        assert_eq!(by_seq.get("AAAAAAAAAACCCCCCCCCC"), Some(&(50, 70)));
    }

    // --- get_read_cuts sentinel fix: region_start == align_start ---
    //
    // Before this fix, a window whose clamped left edge landed exactly on a
    // record's own align_start (very common: it's what eff_start clamps to
    // for any record whose window extends past the record's own left edge)
    // was silently discarded as "invalid coordinates", regardless of
    // whether stitching was in play at all. This is a single-record,
    // no-chain-involved case, isolating that this is fixed for ordinary
    // partial-overlap extraction, not just as a side effect of stitching.

    #[test]
    fn window_start_exactly_at_align_start_now_extracts_correctly() {
        // Single 20M record, target [1000,1020), query [0,20). Window
        // [1000,1010) starts exactly at align_start and ends 10bp in.
        let paf_contents = "q1\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:20M\n";
        let fasta_contents = ">q1\nAAAAAAAAAACCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entries = [&entry];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1000,
            1010,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        let (sequence, _, query_start, query_end, ..) = &reads[0];
        assert_eq!((*query_start, *query_end), (0, 10));
        assert_eq!(sequence, "AAAAAAAAAA");
    }

    #[test]
    fn window_end_exactly_at_align_end_now_extracts_correctly() {
        // Mirror case: window ends exactly at the record's own target_end.
        let paf_contents = "q1\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:20M\n";
        let fasta_contents = ">q1\nAAAAAAAAAACCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entries = [&entry];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1010,
            1020,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        let (sequence, _, query_start, query_end, ..) = &reads[0];
        assert_eq!((*query_start, *query_end), (10, 20));
        assert_eq!(sequence, "CCCCCCCCCC");
    }

    // --- window entirely inside a deletion: valid zero-length result ---

    #[test]
    fn window_entirely_inside_deletion_returns_empty_sequence_not_skipped() {
        // target [1000,1020): 5M (1000-1005), 10D (1005-1015), 5M (1015-1020);
        // query [0,10). A window fully inside the 10bp deletion has no
        // corresponding query bases — read_start == read_end is the correct
        // answer here, not "invalid coordinates".
        let paf_contents = "q1\t10\t0\t10\t+\tchr1\t2000\t1000\t1020\t10\t20\t60\tcg:Z:5M10D5M\n";
        let fasta_contents = ">q1\nAAAAACCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entries = [&entry];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1007,
            1012,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        let (sequence, _, query_start, query_end, ..) = &reads[0];
        assert_eq!(sequence, "");
        assert_eq!((*query_start, *query_end), (5, 5));
    }

    // --- graceful degradation: one bad record shouldn't error the whole call ---

    #[test]
    fn invalid_cigar_in_one_record_is_skipped_not_an_error() {
        let paf_contents = "q1\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:not_a_cigar\n\
                             q2\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:20M\n";
        let fasta_contents = ">q1\nAAAAAAAAAAAAAAAAAAAA\n>q2\nCCCCCCCCCCCCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_bad = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_good = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1000,
            target_end: 1020,
        };
        let entries = [&entry_bad, &entry_good];

        // Must not error — the bad record is skipped, the good one still succeeds.
        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1000,
            1020,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].1, "q2");
    }

    #[test]
    fn missing_fasta_contig_in_one_record_is_skipped_not_an_error() {
        // q_missing isn't in the FASTA at all (e.g. a PAF built against a
        // different assembly than --query_ref) — must not abort the whole run.
        let paf_contents = "q_missing\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:20M\n\
                             q2\t20\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:20M\n";
        let fasta_contents = ">q2\nCCCCCCCCCCCCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let entry_bad = PafIndexEntry {
            offset: 0,
            target_start: 1000,
            target_end: 1020,
        };
        let entry_good = PafIndexEntry {
            offset: paf_contents.lines().next().unwrap().len() as u64 + 1,
            target_start: 1000,
            target_end: 1020,
        };
        let entries = [&entry_bad, &entry_good];

        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entries,
            1000,
            1020,
            0,
            0,
            StitchConfig::default(),
            false,
        )
        .unwrap();

        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].1, "q2");
    }

    // --- deterministic chain selection when multiple candidate chains qualify ---

    #[test]
    fn stitch_picks_the_chain_whose_key_appears_first_in_records_order() {
        // Two independent, equally-valid 2-record chains (qA and qB) both
        // cover the same window. Whichever's key appears first in the PAF
        // (and therefore in `records`) must win, deterministically — not
        // whichever HashMap grouping happens to iterate first.
        let paf_contents = "qA\t70\t0\t20\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:5M5M5M5M\n\
                             qA\t70\t50\t70\t+\tchr1\t2000\t1020\t1040\t20\t20\t60\tcg:Z:5M5M5M5M\n\
                             qB\t70\t100\t120\t+\tchr1\t2000\t1000\t1020\t20\t20\t60\tcg:Z:5M5M5M5M\n\
                             qB\t70\t150\t170\t+\tchr1\t2000\t1020\t1040\t20\t20\t60\tcg:Z:5M5M5M5M\n";
        let fasta_contents = ">qA\nCCCCCCCCCCGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAACCCCCCCCCC\n\
                              >qB\nCCCCCCCCCCGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAACCCCCCCCCC\n";
        let (_paf_file, mut paf_reader, mut fasta_reader) =
            setup_stitch_reader(paf_contents, fasta_contents);

        let mut offset = 0u64;
        let mut offsets = Vec::new();
        for line in paf_contents.lines() {
            offsets.push(offset);
            offset += line.len() as u64 + 1;
        }
        let entries: Vec<PafIndexEntry> = offsets
            .iter()
            .map(|&o| PafIndexEntry {
                offset: o,
                target_start: 1000,
                target_end: 1040,
            })
            .collect();
        let entry_refs: Vec<&PafIndexEntry> = entries.iter().collect();

        // qA's records come first in the PAF (and therefore in `records`),
        // so qA's chain must win — deterministically, not by HashMap luck.
        let reads = get_paf_reads(
            &mut paf_reader,
            &mut fasta_reader,
            &entry_refs,
            1010,
            1030,
            0,
            0,
            StitchConfig {
                enabled: true,
                max_gap: 100,
            },
            false,
        )
        .unwrap();

        // qB's individual records still overlap the window (just not as the
        // winning chain), so they still contribute their own ordinary partial
        // fragments via the normal per-record loop — only the *stitched*
        // (50bp, spanning the full insertion) result is the one whose winner
        // this test cares about.
        let stitched: Vec<_> = reads.iter().filter(|(seq, ..)| seq.len() == 50).collect();
        assert_eq!(
            stitched.len(),
            1,
            "expected exactly one stitched (50bp) result"
        );
        assert_eq!(
            stitched[0].1, "qA",
            "qA's chain appears first in records and must win deterministically"
        );
    }
}
