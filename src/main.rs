mod cli;
mod utils;
mod bed;
mod reads;
mod cigar;
mod paf;

use std::fs::File;
use std::io::BufWriter;
use std::fs::OpenOptions;
use anyhow::Result;

use anyhow::Ok;
use clap::Parser;
use noodles::bam;

use cli::Opts;
use reads::get_bam_reads;
use utils::{read_bed, write_fasta_record};
use paf::PafIndex;

use crate::paf::read_paf_record_at_offset;
use crate::reads::ToCigarOps;
use crate::utils::extract_from_fasta_coords;
use crate::utils::get_read_cuts;

fn main() -> Result<()>{
    let opts: Opts = Opts::parse();
    eprintln!("{:#?}", opts);
    crate::cli::check_option_values(&opts);
    crate::cli::check_inputs_exist(&opts);

    
    eprintln!("Reading bed file");
    let regions: Vec<(noodles::core::Region, String, String)> = read_bed(&opts);

    let output_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true) 
        .open(&opts.output)?;

    let mut read_writer: BufWriter<File> = BufWriter::new(output_file);

    // if reference input
    eprintln!("Reference mode");
    eprintln!("Extracting sequences");
    // for region in bed
    // cut out sequence
    // write to fasta
    
    // if bam
    if opts.bam.to_str() != Some("None"){
        eprintln!("Bam mode");
        eprintln!("Extracting sequences");
        extract_from_bam(&opts, regions, &mut read_writer);
    } 
    
    // if paf
    else if opts.paf.to_str() != Some("None") &&  opts.query_ref.to_str() == Some("None"){
        eprintln!("paf mode");
        eprintln!("Extracting sequences");
        extract_from_paf(&opts, regions, &mut read_writer);
               
    }  

    eprintln!("Done");
    Ok(())
}

pub fn extract_from_bam(opts: &Opts, regions: Vec<(noodles::core::Region, String, String)>, read_writer: &mut BufWriter<File>) {
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
        let overlapping_reads: Vec<(String, Vec<u8>, String, usize, usize)> = get_bam_reads(&opts, query, &region);
        if overlapping_reads.len() == 0 {
            eprintln!("No reads found for region in bam file. Skipping region: {}", region_name);
            continue;
        }
        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());
        // write to fasta or fastq
        for (name, subseq, _subqual, _ref_start, _ref_end) in overlapping_reads {
            let head = format!("{}|{}:{:?}-{:?}|{}", name, chr, region_start, region_end, region_name);
            write_fasta_record(read_writer, &head, &std::str::from_utf8(&subseq).expect("unexpected utf8 in sequence")).expect("Couldn't write fasta record");
        }
        // if consensus: generate consensus
        // write to consensus fasta (potential fastq using mean q score per base?)
        }

}


pub fn extract_from_paf(opts: &Opts, regions: Vec<(noodles::core::Region, String, String)>, read_writer: &mut BufWriter<File>) {
    // Build or load index
    let paf_path: &str = opts.paf.to_str().expect("couldn't get paf string");
    let query_ref = opts.query_ref.to_str().expect("Couldn't get query_ref string");
    let index_path = format!("{}.idx", paf_path);
    let index = if opts.use_paf_index {
        if std::path::Path::new(&index_path).exists() {
            eprintln!("Loading PAF index from {}", index_path);
            PafIndex::load(&index_path).expect("couldn't load index")
        // if can't load it, build it
        } else {
            eprintln!("Building PAF index...");
            let index = PafIndex::build(paf_path).expect("couldn't build paf index");
            index.save(&index_path).expect("Failed to save paf index");
            eprintln!("Index saved to {}", index_path);
            index
        }
    // else build it
    } else {
        PafIndex::build(paf_path).expect("Couldn't build paf index")
    };

    // for each region, get paf regions and extract sequences
    for (region, region_name, chr) in regions.iter() {
        eprintln!("===============================");
        eprintln!("Analysing region: {}, {}",region, region_name);
        eprintln!("===============================");

        if region.name().contains(&('#' as u8)) {
            eprintln!("Region {} has a #, skipping", region_name);
            continue;
        }
        
        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());
        // Query index for overlapping entries
        let overlapping_entries = index.query(chr, region_start, region_end);
        eprintln!("Found {} overlapping alignments", overlapping_entries.len());
        
        // TODO: move this stuff to the reads.rs file under get_paf_reads
        // Read actual PAF records from file
        for entry in overlapping_entries {
            let paf_record = read_paf_record_at_offset(paf_path, entry.offset).unwrap();
            
            if let Some(cigar_str) = &paf_record.cigar {
                // Convert CIGAR and calculate query coordinates
                let cigar_ops = cigar_str.as_str().to_cigar_ops();
                let cuts = get_read_cuts(
                    &cigar_ops,
                    paf_record.target_start,
                    region_start,
                    region_end,
                );
                
                // Calculate actual query coordinates
                let query_start = paf_record.query_start + cuts.read_start;
                let query_end = paf_record.query_start + cuts.read_end;
                
                eprintln!("Query coords: {}:{}-{}", paf_record.query_name, query_start, query_end);
                
                // Extract from query FASTA
                let sequence = extract_from_fasta_coords(query_ref, &paf_record.query_name, query_start, query_end).expect("couldn't extract fasta sequence");
                
                // Write fasta output
                let header = format!(">{}|ref_{}:{}-{}|query_{}:{}-{}", 
                                                paf_record.query_name,
                                                region_name, 
                                                region_start, 
                                                region_end,
                                                paf_record.query_name,
                                                query_start,
                                                query_end);
                write_fasta_record(read_writer, &header, &sequence).expect("Failed to write fasta");
            }
        }
    }
}