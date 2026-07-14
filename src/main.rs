mod cli;

use anyhow::{Context, Result};
use bedpull::paf::PafIndex;
use bedpull::reads::{BamConfig, get_bam_reads, get_cram_reads, get_paf_reads};
use bedpull::utils::{read_bed, write_fasta_record, write_fastq_record};
use clap::Parser;
use noodles::bam;
use noodles::cram;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader as FastaIndexedReader;
use std::collections::HashSet;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
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
    if opts.debug {
        eprintln!("{:#?}", opts);
    }
    crate::cli::check_option_values(&opts)?;
    crate::cli::check_inputs_exist(&opts)?;

    if opts.debug {
        eprintln!("Reading bed file");
    }
    let regions = read_bed(&opts.bed, opts.debug)?;

    let mut read_writer: Box<dyn std::io::Write> = if cli::is_stdout(&opts.output) {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        let output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&opts.output)
            .with_context(|| format!("failed to open output file: {}", opts.output.display()))?;
        Box::new(BufWriter::new(output_file))
    };

    if opts.bam.to_str() != Some("None") {
        eprintln!("BAM mode");
        eprintln!("Extracting sequences");
        extract_from_bam(&opts, regions, read_writer.as_mut())?;
    } else if opts.cram.to_str() != Some("None") {
        eprintln!("CRAM mode");
        eprintln!("Extracting sequences");
        extract_from_cram(&opts, regions, read_writer.as_mut())?;
    } else if opts.paf.to_str() != Some("None") && opts.query_ref.to_str() != Some("None") {
        eprintln!("PAF mode");
        eprintln!("Extracting sequences");
        extract_from_paf(&opts, regions, read_writer.as_mut())?;
    }

    eprintln!("Done");
    Ok(())
}

pub fn extract_from_bam(
    opts: &Opts,
    regions: Vec<(noodles::core::Region, String, String)>,
    read_writer: &mut dyn std::io::Write,
) -> Result<()> {
    let mut seen: HashSet<String> = HashSet::new();

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
        if opts.debug {
            eprintln!("===============================");
            eprintln!("Analysing region: {}, {}", region, region_name);
            eprintln!("===============================");
        }

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
        let region_start =
            region.interval().start().map(usize::from).ok_or_else(|| {
                anyhow::anyhow!("BED region '{}' has unbounded start", region_name)
            })?;
        let region_end = region
            .interval()
            .end()
            .map(usize::from)
            .ok_or_else(|| anyhow::anyhow!("BED region '{}' has unbounded end", region_name))?;
        // write to fasta or fastq
        for (name, subseq, subqual, _ref_start, _ref_end, hap) in overlapping_reads {
            if opts.dedup && !seen.insert(name.clone()) {
                continue;
            }
            let hap_suffix = if hap > 0 {
                format!("|h{}", hap)
            } else {
                String::new()
            };
            let head = format!(
                "{}|{}:{}-{}|{}{}",
                name, chr, region_start, region_end, region_name, hap_suffix
            );
            let seq_str =
                std::str::from_utf8(&subseq).context("BAM sequence contains invalid UTF-8")?;
            let writer: &mut dyn std::io::Write = match hap_writers.as_mut() {
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
    }
    Ok(())
}

pub fn extract_from_cram(
    opts: &Opts,
    regions: Vec<(noodles::core::Region, String, String)>,
    read_writer: &mut dyn std::io::Write,
) -> Result<()> {
    let mut seen: HashSet<String> = HashSet::new();

    // Open per-haplotype writers when --hap_split is set.
    let mut hap_writers: Option<[BufWriter<File>; 3]> = None;
    if opts.hap_split {
        hap_writers = Some([
            open_writer(&hap_output_path(&opts.output, 0))?,
            open_writer(&hap_output_path(&opts.output, 1))?,
            open_writer(&hap_output_path(&opts.output, 2))?,
        ]);
    }

    // Build a reference sequence repository from --reference if provided.
    let reference_repo = if opts.reference.to_str() != Some("None") {
        let indexed = fasta::io::indexed_reader::Builder::default()
            .build_from_path(&opts.reference)
            .with_context(|| {
                format!(
                    "failed to open reference FASTA: {}",
                    opts.reference.display()
                )
            })?;
        fasta::Repository::new(FastaIndexedReader::new(indexed))
    } else {
        fasta::Repository::default()
    };

    for (region, region_name, chr) in regions.iter() {
        if opts.debug {
            eprintln!("===============================");
            eprintln!("Analysing region: {}, {}", region, region_name);
            eprintln!("===============================");
        }

        if region.name().contains(&b'#') {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }

        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(reference_repo.clone())
            .build_from_path(&opts.cram)
            .context("failed to open CRAM file")?;
        let header = reader.read_header().context("failed to read CRAM header")?;
        let query = reader
            .query(&header, region)
            .context("CRAM region query failed")?;

        let (lflank, rflank) = effective_flanks(opts);
        let overlapping_reads = get_cram_reads(&bam_config(opts), query, region, lflank, rflank)?;
        if overlapping_reads.is_empty() {
            eprintln!(
                "No reads found for region in CRAM file. Skipping region: {}",
                region_name
            );
            continue;
        }
        let region_start =
            region.interval().start().map(usize::from).ok_or_else(|| {
                anyhow::anyhow!("BED region '{}' has unbounded start", region_name)
            })?;
        let region_end = region
            .interval()
            .end()
            .map(usize::from)
            .ok_or_else(|| anyhow::anyhow!("BED region '{}' has unbounded end", region_name))?;

        for (name, subseq, subqual, _ref_start, _ref_end, hap) in overlapping_reads {
            if opts.dedup && !seen.insert(name.clone()) {
                continue;
            }
            let hap_suffix = if hap > 0 {
                format!("|h{}", hap)
            } else {
                String::new()
            };
            let head = format!(
                "{}|{}:{}-{}|{}{}",
                name, chr, region_start, region_end, region_name, hap_suffix
            );
            let seq_str =
                std::str::from_utf8(&subseq).context("CRAM sequence contains invalid UTF-8")?;
            let writer: &mut dyn std::io::Write = match hap_writers.as_mut() {
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
    }
    Ok(())
}

pub fn extract_from_paf(
    opts: &Opts,
    regions: Vec<(noodles::core::Region, String, String)>,
    read_writer: &mut dyn std::io::Write,
) -> Result<()> {
    let mut seen: HashSet<String> = HashSet::new();

    // Open per-haplotype writers when --hap_split is set.
    let mut hap_writers: Option<[BufWriter<File>; 3]> = None;
    if opts.hap_split {
        hap_writers = Some([
            open_writer(&hap_output_path(&opts.output, 0))?,
            open_writer(&hap_output_path(&opts.output, 1))?,
            open_writer(&hap_output_path(&opts.output, 2))?,
        ]);
    }

    // Open --bed_out writer if requested.
    let mut bed_out_writer: Option<Box<dyn std::io::Write>> =
        if opts.bed_out.to_str() != Some("None") {
            if cli::is_stdout(&opts.bed_out) {
                Some(Box::new(BufWriter::new(std::io::stdout())))
            } else {
                let f = OpenOptions::new()
                    .write(true)
                    .create(true)
                    .truncate(true)
                    .open(&opts.bed_out)
                    .with_context(|| {
                        format!("failed to open bed_out file: {}", opts.bed_out.display())
                    })?;
                Some(Box::new(BufWriter::new(f)))
            }
        } else {
            None
        };

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
        if opts.debug {
            eprintln!("===============================");
            eprintln!("Analysing region: {}, {}", region, region_name);
            eprintln!("===============================");
        }

        if region.name().contains(&b'#') {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }

        let region_start =
            region.interval().start().map(usize::from).ok_or_else(|| {
                anyhow::anyhow!("BED region '{}' has unbounded start", region_name)
            })?;
        let region_end = region
            .interval()
            .end()
            .map(usize::from)
            .ok_or_else(|| anyhow::anyhow!("BED region '{}' has unbounded end", region_name))?;
        let (lflank, rflank) = effective_flanks(opts);

        // Query index for overlapping entries
        let overlapping_entries = index.query(chr, region_start, region_end);
        if opts.debug {
            eprintln!("Found {} overlapping alignments", overlapping_entries.len());
        }

        let reads = get_paf_reads(
            paf_path,
            query_ref,
            &overlapping_entries,
            region_start,
            region_end,
            lflank,
            rflank,
            opts.debug,
        )?;

        for (sequence, query_name, query_start, query_end, strand, hap) in reads {
            if opts.dedup && !seen.insert(query_name.clone()) {
                continue;
            }
            if let Some(bed_writer) = bed_out_writer.as_mut() {
                writeln!(
                    bed_writer,
                    "{}\t{}\t{}\t{}\t0\t{}",
                    query_name, query_start, query_end, region_name, strand
                )
                .context("failed to write BED record")?;
            }

            let hap_suffix = if hap > 0 {
                format!("|h{}", hap)
            } else {
                String::new()
            };
            let header = format!(
                "{}|{}:{}-{}|{}|{}:{}-{}|{}{}",
                query_name,
                chr,
                region_start,
                region_end,
                region_name,
                query_name,
                query_start,
                query_end,
                strand,
                hap_suffix
            );
            let writer: &mut dyn std::io::Write = match hap_writers.as_mut() {
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
            write_fasta_record(writer, &header, &sequence)
                .context("failed to write FASTA record")?;
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
