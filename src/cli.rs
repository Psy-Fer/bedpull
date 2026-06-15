use anyhow::{Result, bail};
use clap::{Parser, crate_version};
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

    /// Use paf index
    #[clap(long = "use_paf_index", default_value = "true", display_order = 8)]
    pub use_paf_index: bool,

    /// Write FASTQ instead of FASTA (BAM only)
    #[clap(long = "fastq", display_order = 8)]
    pub fastq: bool,

    /// Minimum mapping quality (MAPQ) to include a read (BAM only; 0 = no filter)
    #[clap(long = "min_mapq", default_value = "0", display_order = 8)]
    pub min_mapq: u8,

    /// Include secondary alignments (BAM only; default: skip)
    #[clap(long = "include_secondary", display_order = 8)]
    pub include_secondary: bool,

    /// Include supplementary alignments (BAM only; default: skip)
    #[clap(long = "include_supplementary", display_order = 8)]
    pub include_supplementary: bool,

    /// Include reads that only partially overlap a BED region (default: spanning reads only)
    #[clap(long = "partial", display_order = 9)]
    pub partial: bool,

    /// Minimum mean Phred quality of the extracted region to include a read (BAM --fastq only; 0 = no filter)
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

pub fn check_inputs_exist(opts: &Opts) -> Result<()> {
    if opts.bam.to_str() != Some("None") {
        check_if_file_exists(&opts.bam)?;
        if bam_index_path(&opts.bam).is_none() {
            bail!(
                "BAM index not found for '{}'\n\
                 Create it with:\n  samtools index {}",
                opts.bam.display(),
                opts.bam.display()
            );
        }
    }
    if opts.cram.to_str() != Some("None") {
        check_if_file_exists(&opts.cram)?;
        if cram_index_path(&opts.cram).is_none() {
            bail!(
                "CRAM index not found for '{}'\n\
                 Create it with:\n  samtools index {}",
                opts.cram.display(),
                opts.cram.display()
            );
        }
    }
    if opts.reference.to_str() != Some("None") {
        check_if_file_exists(&opts.reference)?;
        if fasta_index_path(&opts.reference).is_none() {
            bail!(
                "FASTA index not found for '{}'\n\
                 Create it with:\n  samtools faidx {}",
                opts.reference.display(),
                opts.reference.display()
            );
        }
    }
    if opts.paf.to_str() != Some("None") {
        check_if_file_exists(&opts.paf)?;
    }
    if opts.query_ref.to_str() != Some("None") {
        check_if_file_exists(&opts.query_ref)?;
        if fasta_index_path(&opts.query_ref).is_none() {
            bail!(
                "FASTA index not found for '{}'\n\
                 Create it with:\n  samtools faidx {}",
                opts.query_ref.display(),
                opts.query_ref.display()
            );
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
            use_paf_index: true,
            fastq,
            min_mapq: 0,
            include_secondary: false,
            include_supplementary: false,
            partial: false,
            flanks: 0,
            lflank: 0,
            rflank: 0,
            min_region_quality: 0.0,
            hap_split: false,
            dedup: false,
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
    fn missing_bam_index_error_contains_samtools_command() {
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
    fn missing_query_ref_index_error_contains_samtools_command() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("assembly.fa");
        std::fs::write(&fa, b"").unwrap();
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
