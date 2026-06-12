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
    /// Aligned bam file
    #[clap(
        short = 'b',
        long = "bam",
        parse(from_os_str),
        default_value = "None",
        display_order = 1
    )]
    pub bam: PathBuf,

    /// Reference *.fa/fasta
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

    // /// Create a consensus sequence from extracted sequences for each region
    // #[clap(short = 'c', long = "consensus", display_order = 5)]
    // pub consensus: bool,

    // /// Consensus fasta file output
    // #[clap(long = "c_output", default_value = "consensus.bedpull.fasta", display_order = 5)]
    // pub c_output: PathBuf,

    // /// Minimum number of reads in a haplotype group to build a consensus
    // #[clap(long = "min_read_count", default_value = "3", display_order = 6)]
    // pub min_read_count: usize,

    // /// Split sequences based on haplotype tag HP
    // #[clap(short = 'h', long = "hap_split", display_order = 7)]
    // pub hap_split: bool,

    // /// Use 4th column of bed to add name to output
    // #[clap(short = 'n', long = "name", display_order = 8)]
    // pub name: bool,
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

    /// Split output by haplotype HP tag into separate files (BAM only; reads without HP tag go to h0)
    #[clap(long = "hap_split", display_order = 11)]
    pub hap_split: bool,

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
    if opts.reference.to_str() != Some("None") {
        check_if_file_exists(&opts.reference)?;
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
    if opts.fastq && opts.bam.to_str() == Some("None") {
        bail!("--fastq requires --bam");
    }
    if opts.min_region_quality > 0.0 && opts.bam.to_str() == Some("None") {
        bail!(
            "--min_region_quality requires --bam (quality scores are only available from BAM input)"
        );
    }
    if opts.hap_split && is_stdout(&opts.output) {
        bail!("--hap_split requires --output <file> (cannot split haplotypes to stdout)");
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
            reference: PathBuf::from("None"),
            bed: PathBuf::from("test.bed"),
            paf: PathBuf::from("None"),
            query_ref: PathBuf::from("None"),
            output: PathBuf::from("out.fasta"),
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
    fn no_fastq_no_bam_is_ok() {
        assert!(check_option_values(&make_opts(false, "None")).is_ok());
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
