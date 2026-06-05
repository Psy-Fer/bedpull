mod cli;

use anyhow::{Context, Result};
use bedpull::ToCigarOps;
use bedpull::paf::PafIndex;
use bedpull::paf::read_paf_record_at_offset;
use bedpull::reads::{BamConfig, get_bam_reads};
use bedpull::utils::{
    extract_from_fasta_coords, get_read_cuts, read_bed, revcomp, write_fasta_record,
    write_fastq_record,
};
use clap::Parser;
use noodles::bam;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

/// Replace empty `VN:` fields in @PG and @RG header lines.
/// Some samtools versions write `VN:\t` (empty value), which noodles rejects.
fn sanitize_bam_header_text(text: &str) -> String {
    text.lines()
        .map(|line| {
            if !line.starts_with("@PG") && !line.starts_with("@RG") {
                return line.to_string();
            }
            line.split('\t')
                .map(|f| if f == "VN:" { "VN:unknown" } else { f })
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Open the BAM file with a raw BGZF reader, extract and sanitize the SAM header
/// text, then parse it into a `sam::Header`.  Used as a fallback when noodles'
/// strict parser rejects the header (e.g. empty `VN:` in `@PG` records).
fn read_bam_header_lenient(path: &Path) -> Result<noodles::sam::Header> {
    use std::io::Read;

    let file = File::open(path).context("failed to open BAM file")?;
    let mut reader = noodles::bgzf::io::Reader::new(file);

    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .context("failed to read BAM magic bytes")?;
    if &magic != b"BAM\x01" {
        anyhow::bail!("not a valid BAM file (bad magic)");
    }

    let mut len_buf = [0u8; 4];
    reader
        .read_exact(&mut len_buf)
        .context("failed to read BAM header length")?;
    let l_text = u32::from_le_bytes(len_buf) as usize;

    let mut raw_text = vec![0u8; l_text];
    reader
        .read_exact(&mut raw_text)
        .context("failed to read BAM header text")?;

    let text = std::str::from_utf8(&raw_text)
        .context("BAM header text is not valid UTF-8")?
        .trim_end_matches('\0');

    let sanitized = sanitize_bam_header_text(text);

    sanitized
        .parse::<noodles::sam::Header>()
        .map_err(|e| anyhow::anyhow!("failed to parse sanitized BAM header: {e}"))
}

use cli::Opts;

fn effective_flanks(opts: &Opts) -> (usize, usize) {
    cli::resolve_flanks(opts.flanks, opts.lflank, opts.rflank)
}

fn bam_config(opts: &Opts) -> BamConfig {
    BamConfig {
        min_mapq: opts.min_mapq,
        include_secondary: opts.include_secondary,
        include_supplementary: opts.include_supplementary,
        partial: opts.partial,
        min_region_quality: opts.min_region_quality,
    }
}

fn hap_output_path(base: &Path, hap: u8) -> PathBuf {
    let mut name = base
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    name.push_str(&format!(".h{hap}"));
    if let Some(ext) = base.extension() {
        name.push('.');
        name.push_str(&ext.to_string_lossy());
    }
    base.parent().unwrap_or(Path::new(".")).join(name)
}

fn open_writer(path: &Path) -> Result<BufWriter<File>> {
    let f = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)
        .with_context(|| format!("failed to open output file: {}", path.display()))?;
    Ok(BufWriter::new(f))
}

fn main() -> Result<()> {
    let opts: Opts = Opts::parse();
    eprintln!("{:#?}", opts);
    crate::cli::check_option_values(&opts)?;
    crate::cli::check_inputs_exist(&opts)?;

    eprintln!("Reading bed file");
    let regions = read_bed(&opts.bed)?;

    let output_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&opts.output)?;

    let mut read_writer: BufWriter<File> = BufWriter::new(output_file);

    // if reference input
    // eprintln!("Reference mode");
    // eprintln!("Extracting sequences");
    // for region in bed
    // cut out sequence
    // write to fasta

    // if bam
    if opts.bam.to_str() != Some("None") {
        eprintln!("Bam mode");
        eprintln!("Extracting sequences");
        extract_from_bam(&opts, regions, &mut read_writer)?;
    }
    // if paf
    else if opts.paf.to_str() != Some("None") && opts.query_ref.to_str() != Some("None") {
        eprintln!("paf mode");
        eprintln!("Extracting sequences");
        extract_from_paf(&opts, regions, &mut read_writer)?;
    }

    eprintln!("Done");
    Ok(())
}

pub fn extract_from_bam(
    opts: &Opts,
    regions: Vec<(noodles::core::Region, String, String)>,
    read_writer: &mut BufWriter<File>,
) -> Result<()> {
    // Open per-haplotype writers when --hap_split is set.
    let mut hap_writers: Option<[BufWriter<File>; 3]> = None;
    if opts.hap_split {
        hap_writers = Some([
            open_writer(&hap_output_path(&opts.output, 0))?,
            open_writer(&hap_output_path(&opts.output, 1))?,
            open_writer(&hap_output_path(&opts.output, 2))?,
        ]);
    }

    for (region, region_name, chr) in regions.iter() {
        eprintln!("===============================");
        eprintln!("Analysing region: {}, {}", region, region_name);
        eprintln!("===============================");

        if region.name().contains(&b'#') {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }

        // open bam
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&opts.bam)
            .context("failed to open BAM file")?;
        let header = reader
            .read_header()
            .or_else(|e| {
                if e.kind() == std::io::ErrorKind::InvalidData {
                    eprintln!(
                        "Warning: BAM header has non-standard fields (e.g. empty VN: in @PG \
                         records — produced by some samtools versions). Retrying with lenient parser."
                    );
                    read_bam_header_lenient(&opts.bam)
                        .map_err(|ae| std::io::Error::new(std::io::ErrorKind::InvalidData, ae))
                } else {
                    Err(e)
                }
            })
            .context("failed to read BAM header")?;
        let query = reader
            .query(&header, region)
            .context("BAM region query failed")?;

        // find all reads that map to region
        // apply filters (full length, quality, etc)
        // cut out sequence (optionally qstring too and do quality calculation)
        let (lflank, rflank) = effective_flanks(opts);
        let overlapping_reads = get_bam_reads(&bam_config(opts), query, region, lflank, rflank)?;
        if overlapping_reads.is_empty() {
            eprintln!(
                "No reads found for region in bam file. Skipping region: {}",
                region_name
            );
            continue;
        }
        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());
        // write to fasta or fastq
        for (name, subseq, subqual, _ref_start, _ref_end, hap) in overlapping_reads {
            let head = format!(
                "{}|{}:{:?}-{:?}|{}",
                name, chr, region_start, region_end, region_name
            );
            let seq_str =
                std::str::from_utf8(&subseq).context("BAM sequence contains invalid UTF-8")?;
            let writer: &mut BufWriter<File> = match hap_writers.as_mut() {
                Some(writers) => match hap {
                    1 => &mut writers[1],
                    2 => &mut writers[2],
                    _ => {
                        if hap > 2 {
                            eprintln!("Warning: unexpected HP tag value {hap}, routing to h0");
                        }
                        &mut writers[0]
                    }
                },
                None => read_writer,
            };
            if opts.fastq {
                write_fastq_record(writer, &head, seq_str, &subqual)
                    .context("failed to write FASTQ record")?;
            } else {
                write_fasta_record(writer, &head, seq_str)
                    .context("failed to write FASTA record")?;
            }
        }
        // if consensus: generate consensus
        // write to consensus fasta (potential fastq using mean q score per base?)
    }
    Ok(())
}

pub fn extract_from_paf(
    opts: &Opts,
    regions: Vec<(noodles::core::Region, String, String)>,
    read_writer: &mut BufWriter<File>,
) -> Result<()> {
    // Build or load index
    let paf_path = opts
        .paf
        .to_str()
        .context("PAF path contains invalid UTF-8")?;
    let query_ref = opts
        .query_ref
        .to_str()
        .context("query_ref path contains invalid UTF-8")?;
    let index_path = format!("{}.idx", paf_path);
    let index = if opts.use_paf_index {
        if std::path::Path::new(&index_path).exists() {
            eprintln!("Loading PAF index from {}", index_path);
            PafIndex::load(&index_path)
                .with_context(|| format!("failed to load PAF index from {}", index_path))?
        } else {
            eprintln!("Building PAF index...");
            let index = PafIndex::build(paf_path).context("failed to build PAF index")?;
            index
                .save(&index_path)
                .with_context(|| format!("failed to save PAF index to {}", index_path))?;
            eprintln!("Index saved to {}", index_path);
            index
        }
    } else {
        PafIndex::build(paf_path).context("failed to build PAF index")?
    };

    // for each region, get paf regions and extract sequences
    for (region, region_name, chr) in regions.iter() {
        eprintln!("===============================");
        eprintln!("Analysing region: {}, {}", region, region_name);
        eprintln!("===============================");

        if region.name().contains(&b'#') {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }

        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());
        let (lflank, rflank) = effective_flanks(opts);

        // Query index for overlapping entries
        let overlapping_entries = index.query(chr, region_start, region_end);
        eprintln!("Found {} overlapping alignments", overlapping_entries.len());

        // TODO: move this stuff to the reads.rs file under get_paf_reads
        // Read actual PAF records from file
        for entry in overlapping_entries {
            let paf_record = read_paf_record_at_offset(paf_path, entry.offset)
                .with_context(|| format!("failed to read PAF record at offset {}", entry.offset))?;

            // Warn only when the core region (not the flanks) isn't fully covered.
            if paf_record.target_start > region_start {
                eprintln!("Warning: Alignment starts after region start, may be incomplete");
            }
            if paf_record.target_end < region_end {
                eprintln!("Warning: Alignment ends before region end, may be incomplete");
            }

            // Expand by flanks then clamp to what this alignment actually covers.
            let eff_start = region_start
                .saturating_sub(lflank)
                .max(paf_record.target_start);
            let eff_end = (region_end + rflank).min(paf_record.target_end);

            if let Some(cigar_str) = &paf_record.cigar {
                // Convert CIGAR and calculate query coordinates
                let cigar_ops = cigar_str
                    .as_str()
                    .to_cigar_ops()
                    .context("invalid CIGAR string in PAF record")?;
                let cuts = get_read_cuts(&cigar_ops, paf_record.target_start, eff_start, eff_end);

                // Validate cuts
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

                // Calculate actual query coordinates.
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

                // Extract from query FASTA
                let sequence = extract_from_fasta_coords(
                    query_ref,
                    &paf_record.query_name,
                    query_start,
                    query_end,
                )
                .with_context(|| {
                    format!("failed to extract sequence for {}", paf_record.query_name)
                })?;

                // Reverse-complement for minus-strand alignments so the output is
                // always in the reference (target) orientation.
                let sequence = if paf_record.strand == '-' {
                    revcomp(&sequence)
                } else {
                    sequence
                };

                // Write fasta output
                let header = format!(
                    "{}|ref_{}:{}-{}|query_{}:{}-{}",
                    paf_record.query_name,
                    region_name,
                    region_start,
                    region_end,
                    paf_record.query_name,
                    query_start,
                    query_end
                );
                write_fasta_record(read_writer, &header, &sequence)
                    .context("failed to write FASTA record")?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::cli::resolve_flanks;
    use bedpull::ToCigarOps;
    use bedpull::paf::{PafIndex, read_paf_record_at_offset};
    use bedpull::utils::get_read_cuts;

    const PAF_PATH: &str = "examples/hg002pat_to_hs1.rfc1_only.paf";

    // RFC1 BED region: chr4 39318077-39318136 (59 bp reference span)
    // This alignment has a 520 bp insertion inside the region, so bedpull
    // should extract 579 bp from the query — see README for context.
    const RFC1_TARGET_START: usize = 31058861; // PAF field 8  (align_start)
    const RFC1_REGION_START: usize = 39318077; // BED start
    const RFC1_REGION_END: usize = 39318136; // BED end
    const RFC1_EXPECTED_BP: usize = 579;

    #[test]
    fn paf_index_build_finds_rfc1_region() {
        let idx = PafIndex::build(PAF_PATH).expect("failed to build index");
        let hits = idx.query("chr4", RFC1_REGION_START, RFC1_REGION_END);
        assert_eq!(hits.len(), 1, "expected exactly one alignment over RFC1");
    }

    #[test]
    fn paf_record_at_offset_zero_parses_correctly() {
        let r = read_paf_record_at_offset(PAF_PATH, 0).expect("failed to read record");
        assert_eq!(r.query_name, "chr4_PATERNAL");
        assert_eq!(r.target_name, "chr4");
        assert_eq!(r.target_start, RFC1_TARGET_START);
        assert!(r.cigar.is_some(), "expected cg:Z: tag");
    }

    #[test]
    fn get_read_cuts_rfc1_captures_insertion() {
        let r = read_paf_record_at_offset(PAF_PATH, 0).unwrap();
        let cigar_str = r.cigar.as_deref().expect("no CIGAR");
        let ops = cigar_str.to_cigar_ops().expect("valid CIGAR from PAF file");
        let cuts = get_read_cuts(&ops, RFC1_TARGET_START, RFC1_REGION_START, RFC1_REGION_END);
        let extracted_len = cuts.read_end - cuts.read_start;
        assert_eq!(
            extracted_len, RFC1_EXPECTED_BP,
            "expected {RFC1_EXPECTED_BP} bp (59 bp ref span + 520 bp insertion); got {extracted_len}"
        );
    }

    #[test]
    fn resolve_flanks_zero_is_identity() {
        assert_eq!(resolve_flanks(0, 0, 0), (0, 0));
    }
}
