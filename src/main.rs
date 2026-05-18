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
use crate::utils::write_fastq_record;

fn effective_flanks(opts: &Opts) -> (usize, usize) {
    cli::resolve_flanks(opts.flanks, opts.lflank, opts.rflank)
}

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
    // eprintln!("Reference mode");
    // eprintln!("Extracting sequences");
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
    else if opts.paf.to_str() != Some("None") && opts.query_ref.to_str() != Some("None"){
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
        let (lflank, rflank) = effective_flanks(&opts);
        let overlapping_reads: Vec<(String, Vec<u8>, String, usize, usize)> = get_bam_reads(&opts, query, &region, lflank, rflank);
        if overlapping_reads.len() == 0 {
            eprintln!("No reads found for region in bam file. Skipping region: {}", region_name);
            continue;
        }
        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());
        // write to fasta or fastq
        for (name, subseq, subqual, _ref_start, _ref_end) in overlapping_reads {
            let head = format!("{}|{}:{:?}-{:?}|{}", name, chr, region_start, region_end, region_name);
            if opts.fastq {
                write_fastq_record(read_writer, &head, &std::str::from_utf8(&subseq).expect("unexpected utf8 in sequence"), &subqual).expect("Couldn't write fastq record");
            }
            else {
                write_fasta_record(read_writer, &head, &std::str::from_utf8(&subseq).expect("unexpected utf8 in sequence")).expect("Couldn't write fasta record");
            }
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
        let (lflank, rflank) = effective_flanks(opts);

        // Query index for overlapping entries
        let overlapping_entries = index.query(chr, region_start, region_end);
        eprintln!("Found {} overlapping alignments", overlapping_entries.len());

        // TODO: move this stuff to the reads.rs file under get_paf_reads
        // Read actual PAF records from file
        for entry in overlapping_entries {
            let paf_record = read_paf_record_at_offset(paf_path, entry.offset).unwrap();

            // Warn only when the core region (not the flanks) isn't fully covered.
            if paf_record.target_start > region_start {
                eprintln!("Warning: Alignment starts after region start, may be incomplete");
            }
            if paf_record.target_end < region_end {
                eprintln!("Warning: Alignment ends before region end, may be incomplete");
            }

            // Expand by flanks then clamp to what this alignment actually covers.
            let eff_start = region_start.saturating_sub(lflank).max(paf_record.target_start);
            let eff_end = (region_end + rflank).min(paf_record.target_end);

            if let Some(cigar_str) = &paf_record.cigar {
                // Convert CIGAR and calculate query coordinates
                let cigar_ops = cigar_str.as_str().to_cigar_ops();
                let cuts = get_read_cuts(&cigar_ops, paf_record.target_start, eff_start, eff_end);
                
                // Validate cuts
                if cuts.read_start == 0 && cuts.read_end == 0 {
                    eprintln!("Warning: No valid overlap found, skipping");
                    continue;
                }
                
                if cuts.read_start >= cuts.read_end {
                    eprintln!("Warning: Invalid coordinates (start {} >= end {}), skipping", 
                        cuts.read_start, cuts.read_end);
                    continue;
                }
                
                // Calculate actual query coordinates
                let query_start = paf_record.query_start + cuts.read_start;
                let query_end = paf_record.query_start + cuts.read_end;
                
                eprintln!("Query coords: {}:{}-{}", paf_record.query_name, query_start, query_end);
                
                // Extract from query FASTA
                let sequence = extract_from_fasta_coords(query_ref, &paf_record.query_name, query_start, query_end).expect("couldn't extract fasta sequence");
                
                // Write fasta output
                let header = format!("{}|ref_{}:{}-{}|query_{}:{}-{}", 
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

#[cfg(test)]
mod tests {
    use crate::paf::{PafIndex, read_paf_record_at_offset};
    use crate::utils::get_read_cuts;
    use crate::cigar::ToCigarOps;
    use crate::cli::resolve_flanks;

    const PAF_PATH: &str = "examples/hg002pat_to_hs1.rfc1_only.paf";

    // RFC1 BED region: chr4 39318077-39318136 (59 bp reference span)
    // This alignment has a 520 bp insertion inside the region, so bedpull
    // should extract 579 bp from the query — see README for context.
    const RFC1_TARGET_START: usize = 31058861; // PAF field 8  (align_start)
    const RFC1_REGION_START: usize = 39318077; // BED start
    const RFC1_REGION_END:   usize = 39318136; // BED end
    const RFC1_EXPECTED_BP:  usize = 579;

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
        let ops = cigar_str.to_cigar_ops();
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