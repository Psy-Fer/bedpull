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
    #[clap(short = 'b', long = "bam", parse(from_os_str), default_value = "None", display_order = 1)]
    pub bam: PathBuf,
    
    /// Reference *.fa/fasta
    #[clap(short = 'f', long = "reference", parse(from_os_str), default_value = "None", display_order = 2)]
    pub reference: PathBuf,

    /// Bed file of targets
    #[clap(short = 'r', long = "bed", parse(from_os_str), required=true, display_order = 3)]
    pub bed: PathBuf,

    /// paf file - ideally used for assembly to reference mapping
    #[clap(long = "paf", parse(from_os_str), display_order = 3)]
    pub paf: PathBuf,

    /// query reference file (used with paf for extracting sequence)
    #[clap(long = "query_ref", parse(from_os_str), display_order = 3)]
    pub query_ref: PathBuf,

    /// Write a fasta or optionally fastq (bam required) file with extracted query sequences
    #[clap(short = 'o', long = "output", required=true, display_order = 4)]
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
    #[clap(long = "use_paf_index", default_value="true", display_order = 8)]
    pub use_paf_index: bool
}

fn quit_with_error(text: &str) {
    eprintln!("\n\nError: {}", text);
    std::process::exit(1);
}

fn check_if_file_exists(filename: &PathBuf) {
    if !Path::new(filename).exists() {
        let error_msg = format!("The file: {:?} does not exist", filename);
        quit_with_error(&error_msg);
    }
}

pub fn check_inputs_exist(opts: &Opts) {
    if opts.bam.to_str() != Some("None") {
        check_if_file_exists(&opts.bam);
    }
    if opts.reference.to_str() != Some("None") {
        check_if_file_exists(&opts.reference);
    }
    if opts.bam.to_str() != Some("None") {
        let mut bai = opts.bam.clone();
        bai.set_extension("bam.bai");
        check_if_file_exists(&bai);
    }
    if opts.bed.to_str() != Some("None") {
        check_if_file_exists(&opts.bed);
    }
}


pub fn check_option_values(_opts: &Opts) {
    // if opts.min_read_count < 1 {
    //     quit_with_error("--min_read_count must be > 0")
    // }
}


// pub fn get_opts() -> Opts{

//     let opts: Opts = Opts::parse();

//     check_option_values(&opts);
//     check_inputs_exist(&opts);

//     opts
// }