mod cli;
mod utils;
mod bed;
mod reads;

use std::fs::File;
use std::io::BufWriter;
use std::fs::OpenOptions;
use anyhow::Result;

use anyhow::Ok;
use clap::Parser;
use noodles::bam;

use cli::Opts;
use reads::get_reads;
use utils::{read_bed, write_fasta_record};

fn main() -> Result<()>{
    // use bedpull::cli::Opts;
    let opts: Opts = Opts::parse();
    eprintln!("{:#?}", opts);
    crate::cli::check_option_values(&opts);
    crate::cli::check_inputs_exist(&opts);

    
    eprintln!("Reading bed file");
    let regions = read_bed(&opts);

    let output_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true) 
        .open(&opts.output)?;

    let mut read_writer = BufWriter::new(output_file);

    // if reference input
    eprintln!("Reference mode");
    eprintln!("Extracting sequences");
    // for region in bed
    // cut out sequence
    // write to fasta
    
    // if bam
    eprintln!("Bam mode");
    eprintln!("Extracting sequences");
    // for region in bed
    for (region, region_name, chr) in regions.iter() {
        eprintln!("===============================");
        eprintln!("Analysing region: {}, {}",region, region_name);
        eprintln!("===============================");

        if region.name().contains(&('#' as u8)) {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }
    
    // open bam
    let mut reader= bam::io::indexed_reader::Builder::default().build_from_path(&opts.bam).expect("Couldn't read bam");
    let header: noodles::sam::Header = reader.read_header().expect("Couldn't read header");
    let query: bam::io::reader::Query<File> = reader.query(&header, &region).expect("Couldn't find query");
    
    // find all reads that map to region
    // apply filters (full length, quality, etc)
    // cut out sequence (optionally qstring too and do quality calculation)
    let overlapping_reads: Vec<(String, Vec<u8>, String, usize, usize)> = get_reads(&opts, query, &region);
    if overlapping_reads.len() == 0 {
        eprintln!("No reads found for region in bam file. Skipping region: {}", region_name);
        continue;
    }
    // write to fasta or fastq
    for (name, subseq, _subqual, _ref_start, _ref_end) in overlapping_reads {
        let head = format!("{}|{}:{:?}-{:?}|{}", name, chr, region.start(), region.end(), region_name);
        write_fasta_record(&mut read_writer, &head, str::from_utf8(&subseq.to_vec()).expect("unexpected utf8 in sequence")).expect("Couldn't write fasta record")
    }
    // if consensus: generate consensus
    // write to consensus fasta (potential fastq using mean q score per base?)
    }

    eprintln!("Done");
    Ok(())
}