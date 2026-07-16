use anyhow::{Context, Result, bail};
use clap::{Parser, crate_version};
use noodles::{bam, cram, fasta};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
// #[command(author, version, about, long_about = None)]
#[clap(name = "bedpull",
       version = concat!("v", crate_version!()),
       arg_required_else_help = true,
       about = "bedpull - Pull the query sequence from bam or fasta references using a bed file\nhttps://github.com/Psy-Fer/bedpull",
       before_help = concat!(r#"     __             __            ____"#, "\n",
                             r#"    / /_  ___  ____/ /___  __  __/ / /"#, "\n",
                             r#"   / __ \/ _ \/ __  / __ \/ / / / / / "#, "\n",
                             r#"  / /_/ /  __/ /_/ / /_/ / /_/ / / /  "#, "\n",
                             r#" /_.___/\___/\__,_/ .___/\__,_/_/_/   "#, "\n",
                             r#"                 /_/                  "#, "\n",)
        )
        ]
pub struct Opts {
    /// Aligned BAM file
    #[clap(
        short = 'b',
        long = "bam",
        parse(from_os_str),
        default_value = "None",
        display_order = 1
    )]
    pub bam: PathBuf,

    /// Aligned CRAM file (requires --reference for reference-compressed CRAMs)
    #[clap(
        long = "cram",
        parse(from_os_str),
        default_value = "None",
        display_order = 1
    )]
    pub cram: PathBuf,

    /// Reference FASTA (required for CRAM decoding; used with --cram)
    #[clap(
        short = 'f',
        long = "reference",
        parse(from_os_str),
        default_value = "None",
        display_order = 2
    )]
    pub reference: PathBuf,

    /// Bed file of targets
    #[clap(
        short = 'r',
        long = "bed",
        parse(from_os_str),
        required = true,
        display_order = 3
    )]
    pub bed: PathBuf,

    /// paf file - ideally used for assembly to reference mapping
    #[clap(
        long = "paf",
        parse(from_os_str),
        default_value = "None",
        display_order = 3
    )]
    pub paf: PathBuf,

    /// query reference file (used with paf for extracting sequence)
    #[clap(
        long = "query_ref",
        parse(from_os_str),
        default_value = "None",
        display_order = 3
    )]
    pub query_ref: PathBuf,

    /// Output file path. Use '-' or omit to write to stdout.
    #[clap(short = 'o', long = "output", default_value = "-", display_order = 4)]
    pub output: PathBuf,

    /// Write lifted-over query coordinates as BED6 to this file (PAF mode only). Use '-' to write to stdout.
    #[clap(long = "bed_out", default_value = "None", display_order = 5)]
    pub bed_out: PathBuf,

    /// Write input BED regions that produced no output to this file, each preceded by a
    /// '#reason' comment (all modes). Use '-' to write to stdout.
    #[clap(long = "unmapped", default_value = "None", display_order = 5)]
    pub unmapped: PathBuf,

    /// Use paf index
    #[clap(long = "use_paf_index", default_value = "true", display_order = 8)]
    pub use_paf_index: bool,

    /// Write FASTQ instead of FASTA (BAM only)
    #[clap(long = "fastq", display_order = 8)]
    pub fastq: bool,

    /// Minimum mapping quality (MAPQ) to include a read (BAM/CRAM only; 0 = no filter)
    #[clap(long = "min_mapq", default_value = "0", display_order = 8)]
    pub min_mapq: u8,

    /// Include secondary alignments (BAM/CRAM only; default: skip)
    #[clap(long = "include_secondary", display_order = 8)]
    pub include_secondary: bool,

    /// Include supplementary alignments (BAM/CRAM only; default: skip)
    #[clap(long = "include_supplementary", display_order = 8)]
    pub include_supplementary: bool,

    /// Include reads that only partially overlap a BED region (BAM/CRAM only; default: spanning reads only)
    #[clap(long = "partial", display_order = 9)]
    pub partial: bool,

    /// Minimum fraction (0.0-1.0) of the requested region (± flanks) a --partial read must cover to be included (BAM/CRAM only, requires --partial; 0 = no filter, any overlap accepted)
    #[clap(long = "min_partial_coverage", default_value = "0", display_order = 9)]
    pub min_partial_coverage: f64,

    /// Minimum mean Phred quality of the extracted region to include a read (BAM/CRAM only, regardless of --fastq; 0 = no filter)
    #[clap(long = "min_region_quality", default_value = "0", display_order = 10)]
    pub min_region_quality: f64,

    /// Split output by haplotype tag (HP in BAM, hp:i: in PAF) into separate files; reads/alignments without the tag go to h0
    #[clap(long = "hap_split", display_order = 11)]
    pub hap_split: bool,

    /// Deduplicate output: if the same read/contig name is seen more than once (e.g. spanning multiple BED regions), emit it only the first time
    #[clap(long = "dedup", display_order = 12)]
    pub dedup: bool,

    /// Symmetric reference flank to add on both sides before CIGAR walk (bp)
    #[clap(long = "flanks", default_value = "0", display_order = 9)]
    pub flanks: usize,

    /// Left-side reference flank in bp (overrides --flanks for the left side)
    #[clap(long = "lflank", default_value = "0", display_order = 9)]
    pub lflank: usize,

    /// Right-side reference flank in bp (overrides --flanks for the right side)
    #[clap(long = "rflank", default_value = "0", display_order = 9)]
    pub rflank: usize,

    /// Attempt to stitch a region across multiple chained PAF records when no single
    /// record fully spans it (PAF mode only) — recovers cases where a large structural
    /// variant caused the aligner to split one alignment into two or more records
    /// instead of a single record with a large indel operation
    #[clap(long = "stitch_records", display_order = 9)]
    pub stitch_records: bool,

    /// Maximum target-space gap in bp between consecutive records for --stitch_records
    /// to still treat them as part of the same split alignment (PAF mode only, requires
    /// --stitch_records)
    #[clap(long = "max_stitch_gap", default_value = "10000", display_order = 9)]
    pub max_stitch_gap: usize,

    /// Print verbose per-region/per-read diagnostic output (parsed args, BED records, region banners, query coordinates)
    #[clap(long = "debug", display_order = 13)]
    pub debug: bool,
}

fn check_if_file_exists(filename: &PathBuf) -> Result<()> {
    if !Path::new(filename).exists() {
        bail!("file not found: {}", filename.display());
    }
    Ok(())
}

/// Returns the path to a BAM index if one exists, checking both `<bam>.bai`
/// and `<bam>.bam.bai` conventions. Returns `None` if neither is present.
fn bam_index_path(bam: &Path) -> Option<PathBuf> {
    // Most common: foo.bam.bai (samtools index default)
    let bam_bai = bam.with_extension("bam.bai");
    if bam_bai.exists() {
        return Some(bam_bai);
    }
    // Alternative: foo.bai
    let bai = bam.with_extension("bai");
    if bai.exists() {
        return Some(bai);
    }
    None
}

/// Returns the path to a FASTA index (`<fasta>.fai`), or `None` if absent.
fn fasta_index_path(fasta: &Path) -> Option<PathBuf> {
    let mut s = fasta.as_os_str().to_owned();
    s.push(".fai");
    let fai = PathBuf::from(s);
    if fai.exists() { Some(fai) } else { None }
}

/// Returns the CRAI path (`<cram>.crai`) if it exists, `None` otherwise.
fn cram_index_path(cram: &Path) -> Option<PathBuf> {
    let mut s = cram.as_os_str().to_owned();
    s.push(".crai");
    let crai = PathBuf::from(s);
    if crai.exists() { Some(crai) } else { None }
}

/// Build and write a `<bam>.bai` index (the `samtools index` default convention).
fn build_bam_index(path: &Path) -> Result<()> {
    eprintln!(
        "BAM index not found for '{}' — building it (this may take a moment for large files)...",
        path.display()
    );
    let index = bam::fs::index(path).with_context(|| {
        format!(
            "failed to auto-build BAM index for '{}'\n\
             You can build it manually with:\n  samtools index {}",
            path.display(),
            path.display()
        )
    })?;
    let dst = path.with_extension("bam.bai");
    bam::bai::fs::write(&dst, &index).with_context(|| {
        format!(
            "failed to write BAM index to '{}'\n\
             You can build it manually with:\n  samtools index {}",
            dst.display(),
            path.display()
        )
    })?;
    eprintln!("BAM index written to {}", dst.display());
    Ok(())
}

/// Build and write a `<cram>.crai` index (the `samtools index` default convention).
fn build_cram_index(path: &Path) -> Result<()> {
    eprintln!(
        "CRAM index not found for '{}' — building it (this may take a moment for large files)...",
        path.display()
    );
    let index = cram::fs::index(path).with_context(|| {
        format!(
            "failed to auto-build CRAM index for '{}'\n\
             You can build it manually with:\n  samtools index {}",
            path.display(),
            path.display()
        )
    })?;
    let mut dst = path.as_os_str().to_owned();
    dst.push(".crai");
    let dst = PathBuf::from(dst);
    cram::crai::fs::write(&dst, &index).with_context(|| {
        format!(
            "failed to write CRAM index to '{}'\n\
             You can build it manually with:\n  samtools index {}",
            dst.display(),
            path.display()
        )
    })?;
    eprintln!("CRAM index written to {}", dst.display());
    Ok(())
}

/// Build and write a `<fasta>.fai` index (the `samtools faidx` convention).
fn build_fasta_index(path: &Path) -> Result<()> {
    eprintln!(
        "FASTA index not found for '{}' — building it (this may take a moment for large files)...",
        path.display()
    );
    let index = fasta::fs::index(path).with_context(|| {
        format!(
            "failed to auto-build FASTA index for '{}'\n\
             You can build it manually with:\n  samtools faidx {}",
            path.display(),
            path.display()
        )
    })?;
    let mut dst = path.as_os_str().to_owned();
    dst.push(".fai");
    let dst = PathBuf::from(dst);
    fasta::fai::fs::write(&dst, &index).with_context(|| {
        format!(
            "failed to write FASTA index to '{}'\n\
             You can build it manually with:\n  samtools faidx {}",
            dst.display(),
            path.display()
        )
    })?;
    eprintln!("FASTA index written to {}", dst.display());
    Ok(())
}

pub fn check_inputs_exist(opts: &Opts) -> Result<()> {
    if opts.bam.to_str() != Some("None") {
        check_if_file_exists(&opts.bam)?;
        if bam_index_path(&opts.bam).is_none() {
            build_bam_index(&opts.bam)?;
        }
    }
    if opts.cram.to_str() != Some("None") {
        check_if_file_exists(&opts.cram)?;
        if cram_index_path(&opts.cram).is_none() {
            build_cram_index(&opts.cram)?;
        }
    }
    if opts.reference.to_str() != Some("None") {
        check_if_file_exists(&opts.reference)?;
        if fasta_index_path(&opts.reference).is_none() {
            build_fasta_index(&opts.reference)?;
        }
    }
    if opts.paf.to_str() != Some("None") {
        check_if_file_exists(&opts.paf)?;
    }
    if opts.query_ref.to_str() != Some("None") {
        check_if_file_exists(&opts.query_ref)?;
        if fasta_index_path(&opts.query_ref).is_none() {
            build_fasta_index(&opts.query_ref)?;
        }
    }
    if opts.bed.to_str() != Some("None") {
        check_if_file_exists(&opts.bed)?;
    }
    Ok(())
}

/// Resolve the three flank args into (lflank, rflank).
/// Per-side values take precedence over the symmetric --flanks shorthand.
pub fn resolve_flanks(flanks: usize, lflank: usize, rflank: usize) -> (usize, usize) {
    let left = if lflank > 0 { lflank } else { flanks };
    let right = if rflank > 0 { rflank } else { flanks };
    (left, right)
}

pub fn is_stdout(output: &std::path::Path) -> bool {
    output.to_str() == Some("-")
}

/// Verify `dir` is writable by actually creating and removing a probe file in it.
/// More reliable than checking permission bits, which don't account for ACLs,
/// read-only mounts, or ownership mismatches.
fn check_dir_writable(dir: &Path) -> Result<()> {
    let probe = dir.join(format!(".bedpull_write_test_{}", std::process::id()));
    std::fs::File::create(&probe)
        .with_context(|| format!("output directory is not writable: {}", dir.display()))?;
    let _ = std::fs::remove_file(&probe);
    Ok(())
}

pub fn check_option_values(opts: &Opts) -> Result<()> {
    let has_bam = opts.bam.to_str() != Some("None");
    let has_cram = opts.cram.to_str() != Some("None");
    let has_paf = opts.paf.to_str() != Some("None");
    let has_query_ref = opts.query_ref.to_str() != Some("None");

    // Conflicting modes
    if has_bam && has_cram {
        bail!("--bam and --cram are mutually exclusive; provide only one");
    }
    if (has_bam || has_cram) && has_paf {
        bail!("--paf cannot be combined with --bam or --cram");
    }

    // PAF requires both halves
    if has_paf && !has_query_ref {
        bail!("--paf requires --query_ref <fasta>");
    }
    if has_query_ref && !has_paf {
        bail!("--query_ref requires --paf <paf>");
    }

    // At least one mode must be specified
    if !has_bam && !has_cram && !has_paf {
        bail!("no input mode specified — provide --bam, --cram, or --paf + --query_ref");
    }

    if opts.fastq && !has_bam && !has_cram {
        bail!("--fastq requires --bam or --cram");
    }
    if opts.min_region_quality > 0.0 && !has_bam && !has_cram {
        bail!(
            "--min_region_quality requires --bam or --cram (quality scores are only available from alignment input)"
        );
    }
    if opts.min_mapq > 0 && !has_bam && !has_cram {
        bail!("--min_mapq requires --bam or --cram (PAF records have no mapping quality)");
    }
    if opts.include_secondary && !has_bam && !has_cram {
        bail!("--include_secondary requires --bam or --cram (not applicable to PAF mode)");
    }
    if opts.include_supplementary && !has_bam && !has_cram {
        bail!("--include_supplementary requires --bam or --cram (not applicable to PAF mode)");
    }
    if opts.partial && !has_bam && !has_cram {
        bail!("--partial requires --bam or --cram (not yet supported in PAF mode)");
    }
    if opts.stitch_records && !has_paf {
        bail!("--stitch_records requires --paf + --query_ref (not applicable to BAM/CRAM mode)");
    }
    if opts.min_partial_coverage > 0.0 && !opts.partial {
        bail!(
            "--min_partial_coverage requires --partial (no-op otherwise: non-partial reads always cover the full requested window)"
        );
    }
    if opts.min_partial_coverage < 0.0 || opts.min_partial_coverage > 1.0 {
        bail!(
            "--min_partial_coverage must be between 0.0 and 1.0 (got {})",
            opts.min_partial_coverage
        );
    }
    if opts.hap_split && is_stdout(&opts.output) {
        bail!("--hap_split requires --output <file> (cannot split haplotypes to stdout)");
    }

    let has_bed_out = opts.bed_out.to_str() != Some("None");
    if has_bed_out && !has_paf {
        bail!("--bed_out is only valid in PAF mode (--paf + --query_ref)");
    }
    if has_bed_out && is_stdout(&opts.bed_out) && is_stdout(&opts.output) {
        bail!(
            "--bed_out and --output cannot both be stdout ('-'); use a file path for one of them"
        );
    }
    if has_bed_out && !is_stdout(&opts.bed_out) {
        let parent = opts.bed_out.parent().unwrap_or(Path::new("."));
        let parent = if parent.as_os_str().is_empty() {
            Path::new(".")
        } else {
            parent
        };
        if !parent.exists() {
            bail!(
                "bed_out directory does not exist: {}\n\
                 Create it with:\n  mkdir -p {}",
                parent.display(),
                parent.display()
            );
        }
        check_dir_writable(parent)?;
    }

    let has_unmapped = opts.unmapped.to_str() != Some("None");
    if has_unmapped && is_stdout(&opts.unmapped) && is_stdout(&opts.output) {
        bail!(
            "--unmapped and --output cannot both be stdout ('-'); use a file path for one of them"
        );
    }
    if has_unmapped && has_bed_out && is_stdout(&opts.unmapped) && is_stdout(&opts.bed_out) {
        bail!(
            "--unmapped and --bed_out cannot both be stdout ('-'); use a file path for one of them"
        );
    }
    if has_unmapped && !is_stdout(&opts.unmapped) {
        let parent = opts.unmapped.parent().unwrap_or(Path::new("."));
        let parent = if parent.as_os_str().is_empty() {
            Path::new(".")
        } else {
            parent
        };
        if !parent.exists() {
            bail!(
                "unmapped directory does not exist: {}\n\
                 Create it with:\n  mkdir -p {}",
                parent.display(),
                parent.display()
            );
        }
        check_dir_writable(parent)?;
    }

    // Output parent directory must exist
    if !is_stdout(&opts.output) {
        let parent = opts.output.parent().unwrap_or(Path::new("."));
        let parent = if parent.as_os_str().is_empty() {
            Path::new(".")
        } else {
            parent
        };
        if !parent.exists() {
            bail!(
                "output directory does not exist: {}\n\
                 Create it with:\n  mkdir -p {}",
                parent.display(),
                parent.display()
            );
        }
        check_dir_writable(parent)?;
    }

    Ok(())
}

// pub fn get_opts() -> Opts{

//     let opts: Opts = Opts::parse();

//     check_option_values(&opts);
//     check_inputs_exist(&opts);

//     opts
// }

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;
    use std::path::PathBuf;

    fn make_opts(fastq: bool, bam: &str) -> Opts {
        Opts {
            bam: PathBuf::from(bam),
            cram: PathBuf::from("None"),
            reference: PathBuf::from("None"),
            bed: PathBuf::from("test.bed"),
            paf: PathBuf::from("None"),
            query_ref: PathBuf::from("None"),
            output: PathBuf::from("out.fasta"),
            bed_out: PathBuf::from("None"),
            unmapped: PathBuf::from("None"),
            use_paf_index: true,
            fastq,
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: false,
            partial: false,
            min_partial_coverage: 0.0,
            flanks: 0,
            lflank: 0,
            rflank: 0,
            min_region_quality: 0.0,
            hap_split: false,
            dedup: false,
            stitch_records: false,
            max_stitch_gap: 10_000,
            debug: false,
        }
    }

    #[test]
    fn fastq_without_bam_is_error() {
        assert!(check_option_values(&make_opts(true, "None")).is_err());
    }

    #[test]
    fn fastq_with_bam_is_ok() {
        assert!(check_option_values(&make_opts(true, "reads.bam")).is_ok());
    }

    #[test]
    fn no_mode_specified_is_error() {
        assert!(check_option_values(&make_opts(false, "None")).is_err());
    }

    #[test]
    fn bam_mode_is_ok() {
        assert!(check_option_values(&make_opts(false, "reads.bam")).is_ok());
    }

    fn make_opts_paf(paf: &str, query_ref: &str) -> Opts {
        let mut o = make_opts(false, "None");
        o.paf = PathBuf::from(paf);
        o.query_ref = PathBuf::from(query_ref);
        o
    }

    #[test]
    fn paf_without_query_ref_is_error() {
        let mut o = make_opts(false, "None");
        o.paf = PathBuf::from("aln.paf");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn query_ref_without_paf_is_error() {
        let mut o = make_opts(false, "None");
        o.query_ref = PathBuf::from("asm.fa");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn paf_with_query_ref_is_ok() {
        assert!(check_option_values(&make_opts_paf("aln.paf", "asm.fa")).is_ok());
    }

    #[test]
    fn bam_and_cram_together_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.cram = PathBuf::from("reads.cram");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn bam_and_paf_together_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.paf = PathBuf::from("aln.paf");
        o.query_ref = PathBuf::from("asm.fa");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn min_mapq_in_paf_mode_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.min_mapq = 20;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn min_mapq_in_bam_mode_is_ok() {
        let mut o = make_opts(false, "reads.bam");
        o.min_mapq = 20;
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn include_secondary_in_paf_mode_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.include_secondary = true;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn include_supplementary_in_paf_mode_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.include_supplementary = true;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn partial_in_paf_mode_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.partial = true;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn partial_in_bam_mode_is_ok() {
        let mut o = make_opts(false, "reads.bam");
        o.partial = true;
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn stitch_records_in_bam_mode_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.stitch_records = true;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn stitch_records_in_paf_mode_is_ok() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.stitch_records = true;
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn min_partial_coverage_without_partial_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.min_partial_coverage = 0.5;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn min_partial_coverage_with_partial_is_ok() {
        let mut o = make_opts(false, "reads.bam");
        o.partial = true;
        o.min_partial_coverage = 0.5;
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn min_partial_coverage_above_one_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.partial = true;
        o.min_partial_coverage = 1.5;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn min_partial_coverage_negative_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.partial = true;
        o.min_partial_coverage = -0.1;
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn min_partial_coverage_zero_with_partial_off_is_ok() {
        // 0.0 is the "disabled" sentinel, so it's fine even without --partial.
        let o = make_opts(false, "reads.bam");
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn output_dir_not_writable_is_error() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o500)).unwrap();
        let mut o = make_opts(false, "reads.bam");
        o.output = dir.path().join("out.fasta");
        let result = check_option_values(&o);
        // restore permissions so the tempdir can be cleaned up
        std::fs::set_permissions(dir.path(), std::fs::Permissions::from_mode(0o700)).unwrap();
        assert!(result.is_err());
    }

    #[test]
    fn output_missing_parent_dir_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.output = PathBuf::from("/nonexistent_dir_bedpull_test/out.fasta");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn output_existing_parent_dir_is_ok() {
        let dir = tempfile::tempdir().unwrap();
        let mut o = make_opts(false, "reads.bam");
        o.output = dir.path().join("out.fasta");
        assert!(check_option_values(&o).is_ok());
    }

    // --- index detection ---

    #[test]
    fn bam_index_found_as_bam_bai() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        let bai = dir.path().join("reads.bam.bai");
        std::fs::write(&bam, b"").unwrap();
        std::fs::write(&bai, b"").unwrap();
        assert!(bam_index_path(&bam).is_some());
    }

    #[test]
    fn bam_index_found_as_bai() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        let bai = dir.path().join("reads.bai");
        std::fs::write(&bam, b"").unwrap();
        std::fs::write(&bai, b"").unwrap();
        assert!(bam_index_path(&bam).is_some());
    }

    #[test]
    fn bam_index_missing_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        std::fs::write(&bam, b"").unwrap();
        assert!(bam_index_path(&bam).is_none());
    }

    #[test]
    fn fasta_index_found() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("assembly.fa");
        let fai = dir.path().join("assembly.fa.fai");
        std::fs::write(&fa, b"").unwrap();
        std::fs::write(&fai, b"").unwrap();
        assert!(fasta_index_path(&fa).is_some());
    }

    #[test]
    fn fasta_index_missing_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("assembly.fa");
        std::fs::write(&fa, b"").unwrap();
        assert!(fasta_index_path(&fa).is_none());
    }

    #[test]
    fn missing_bam_index_is_auto_built() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        std::fs::copy("examples/rfc1_test.bam", &bam).unwrap();
        let bed = dir.path().join("test.bed");
        std::fs::write(&bed, b"chr4\t39318077\t39318136\tRFC1\n").unwrap();
        let mut opts = make_opts(false, "None");
        opts.bam = bam.clone();
        opts.bed = bed;

        assert!(bam_index_path(&bam).is_none());
        check_inputs_exist(&opts).unwrap();
        assert!(bam_index_path(&bam).is_some());
    }

    #[test]
    fn invalid_bam_index_build_error_contains_samtools_command() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("reads.bam");
        std::fs::write(&bam, b"").unwrap();
        let mut opts = make_opts(false, "None");
        opts.bam = bam.clone();
        let err = check_inputs_exist(&opts).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("samtools index"),
            "error should suggest samtools index: {msg}"
        );
        assert!(
            msg.contains("reads.bam"),
            "error should mention the BAM path: {msg}"
        );
    }

    #[test]
    fn missing_query_ref_index_is_auto_built() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("assembly.fa");
        std::fs::write(&fa, b">chr1\nACGTACGTACGT\n").unwrap();
        let bed = dir.path().join("test.bed");
        std::fs::write(&bed, b"chr1\t0\t4\n").unwrap();
        let mut opts = make_opts(false, "None");
        opts.query_ref = fa.clone();
        opts.bed = bed;

        assert!(fasta_index_path(&fa).is_none());
        check_inputs_exist(&opts).unwrap();
        assert!(fasta_index_path(&fa).is_some());
    }

    #[test]
    fn invalid_fasta_index_build_error_contains_samtools_command() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("assembly.fa");
        std::fs::write(&fa, b"not a fasta file\n").unwrap();
        let mut opts = make_opts(false, "None");
        opts.query_ref = fa.clone();
        let err = check_inputs_exist(&opts).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("samtools faidx"),
            "error should suggest samtools faidx: {msg}"
        );
        assert!(
            msg.contains("assembly.fa"),
            "error should mention the FASTA path: {msg}"
        );
    }

    #[test]
    fn missing_cram_index_is_auto_built() {
        let dir = tempfile::tempdir().unwrap();
        let cram = dir.path().join("reads.cram");
        std::fs::copy("examples/rfc1_test.cram", &cram).unwrap();
        let bed = dir.path().join("test.bed");
        std::fs::write(&bed, b"chr4\t39318077\t39318136\tRFC1\n").unwrap();
        let mut opts = make_opts(false, "None");
        opts.cram = cram.clone();
        opts.bed = bed;

        assert!(cram_index_path(&cram).is_none());
        check_inputs_exist(&opts).unwrap();
        assert!(cram_index_path(&cram).is_some());
    }

    #[test]
    fn invalid_cram_index_build_error_contains_samtools_command() {
        let dir = tempfile::tempdir().unwrap();
        let cram = dir.path().join("reads.cram");
        std::fs::write(&cram, b"").unwrap();
        let mut opts = make_opts(false, "None");
        opts.cram = cram.clone();
        let err = check_inputs_exist(&opts).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("samtools index"),
            "error should suggest samtools index: {msg}"
        );
        assert!(
            msg.contains("reads.cram"),
            "error should mention the CRAM path: {msg}"
        );
    }

    // --- bed_out validation ---

    #[test]
    fn bed_out_without_paf_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.bed_out = PathBuf::from("out.bed");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn bed_out_with_paf_is_ok() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.bed_out = PathBuf::from("out.bed");
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn bed_out_and_output_both_stdout_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.output = PathBuf::from("-");
        o.bed_out = PathBuf::from("-");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn bed_out_stdout_with_output_file_is_ok() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.output = PathBuf::from("out.fasta");
        o.bed_out = PathBuf::from("-");
        assert!(check_option_values(&o).is_ok());
    }

    // --- --unmapped validation ---

    #[test]
    fn unmapped_valid_in_bam_mode() {
        let mut o = make_opts(false, "reads.bam");
        o.unmapped = PathBuf::from("unmapped.bed");
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn unmapped_valid_in_paf_mode() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.unmapped = PathBuf::from("unmapped.bed");
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn unmapped_and_output_both_stdout_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.output = PathBuf::from("-");
        o.unmapped = PathBuf::from("-");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn unmapped_and_bed_out_both_stdout_is_error() {
        let mut o = make_opts_paf("aln.paf", "asm.fa");
        o.output = PathBuf::from("out.fasta");
        o.bed_out = PathBuf::from("-");
        o.unmapped = PathBuf::from("-");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn unmapped_stdout_with_output_file_is_ok() {
        let mut o = make_opts(false, "reads.bam");
        o.output = PathBuf::from("out.fasta");
        o.unmapped = PathBuf::from("-");
        assert!(check_option_values(&o).is_ok());
    }

    #[test]
    fn unmapped_missing_parent_dir_is_error() {
        let mut o = make_opts(false, "reads.bam");
        o.unmapped = PathBuf::from("/nonexistent_dir_bedpull_test/unmapped.bed");
        assert!(check_option_values(&o).is_err());
    }

    #[test]
    fn all_zero_gives_zero() {
        assert_eq!(resolve_flanks(0, 0, 0), (0, 0));
    }

    #[test]
    fn symmetric_flanks_applied_to_both_sides() {
        assert_eq!(resolve_flanks(100, 0, 0), (100, 100));
    }

    #[test]
    fn lflank_overrides_left() {
        assert_eq!(resolve_flanks(100, 50, 0), (50, 100));
    }

    #[test]
    fn rflank_overrides_right() {
        assert_eq!(resolve_flanks(100, 0, 75), (100, 75));
    }

    #[test]
    fn per_side_overrides_both() {
        assert_eq!(resolve_flanks(100, 30, 40), (30, 40));
    }

    #[test]
    fn flanks_zero_per_side_nonzero() {
        assert_eq!(resolve_flanks(0, 10, 20), (10, 20));
    }
}
