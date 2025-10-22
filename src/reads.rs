use itertools::Itertools;
use std::fs::File;
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use crate::cli::Opts;

use crate::utils::{get_read_cuts, ReadCuts};
pub use crate::cigar::ToCigarOps;

// For bam reading
pub fn get_reads(_opts: &Opts, query: bam::io::reader::Query<File> , region: &Region) -> Vec<(String, Vec<u8>, String, usize, usize)>{

    let mut h0_subseq_vec: Vec<(String, Vec<u8>, String, usize, usize)> = vec![]; // no hap assigned
    
    // let mut counter = 0;
    for result in query {
        // counter += 1;
        let record = result.expect("Couldn't read result");
        let ref_id = record.reference_sequence_id().unwrap().expect("ref_id had an error");
        // this is a 1-based position
        let align_start = usize::from(record.alignment_start().unwrap().expect("Couldn't get align_start"));
        // TODO: can filter on this using arg
        let map_quality = u8::from(record.mapping_quality().unwrap());
        let flags = record.flags();
        let data = record.data();
        let span: usize = record.alignment_span().unwrap().expect("couldn't get alignment span");
        let align_end = usize::from(record.alignment_end().unwrap().expect("couldn't get alignment end"));
        // turn this to bytes, then a vec, then convert to string from utf8
        let name = String::from_utf8(record.name().unwrap().as_bytes().to_vec()).expect("unexpected utf in name");
        let seq = record.sequence();
        let i_seq = seq.iter().collect_vec();
        // let i_qual =  record.quality_scores().as_ref()
        //                                             .iter()
        //                                             .map(|&score| (score + 33) as u8) // Adjust quality scores
        //                                             .collect::<Vec<_>>();
        // eprintln!("record.quality_scores: {:?}", record.quality_scores());
        // let i_qual = record.quality_scores().as_ref()
        //                             .iter()
        //                             .map(|&score| (score + 33) as u8) // Adjust quality scores
        //                             .collect::<Vec<_>>();

        // eprintln!("quality_scores adjusted: {:?}", i_qual);
        // now convert that to a String
        // let quality_scores_str: String = String::from_utf8_lossy(&i_qual).into_owned();
        let cigar = record.cigar().to_cigar_ops();
        // Convert CIGAR operations to string by getting each kind, converting to a char, and going len|char and collecting
        // let cigar_string: String = cigar_to_string(&cigar);
        // get start and end position in read sequence coordinates using cigar string
        // take alignment start as 0 in read position. Work through cigar to get ref co-ord->read position for subsequence start
        // then continue through read until ref coord end-> read position for subsequence end
        // extract subsequence and chuck it into vector to be worked on
        // read_cuts = (start, end)

        // TODO: filter reads by map_quality
        // if map_quality < opts.min_map_score {
        //     eprintln!("{} read map score {} too low (min: {})", name, map_quality, opts.min_map_score);
        //     counter -= 1;
        //     continue;
        // }

        // TODO: filter reads by mean read quality
        // let read_mean_qscore: f64 = calculate_qscore(&quality_scores_str);
        // if read_mean_qscore < opts.min_read_quality {
        //     eprintln!("{} mean read quality of {} too low (min: {})", name, read_mean_qscore, opts.min_read_quality);
        //     counter -= 1;
        //     continue;
        // }

        // filter reads that don't map across the full region
        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());

        // Skip reads that don't overlap
        if (align_end < region_start) || (align_start > region_end) {
            eprintln!("{} Doesn't fit in region", name);
            eprintln!("A_start/start: {} / {}: {}", align_start, usize::from(region.interval().start().unwrap()), align_start as i64 - usize::from(region.interval().start().unwrap()) as i64);
            eprintln!("A_end/end: {} / {}: {}", align_end, usize::from(region.interval().end().unwrap()), align_end as i64 - usize::from(region.interval().end().unwrap()) as i64);
            continue;
        }

        eprintln!("A_start/start: {} / {}: {}", align_start, usize::from(region.interval().start().unwrap()), align_start as i64 - usize::from(region.interval().start().unwrap()) as i64);
        eprintln!("A_end/end: {} / {}: {}", align_end, usize::from(region.interval().end().unwrap()), align_end as i64 - usize::from(region.interval().end().unwrap()) as i64);
        let read_cuts: ReadCuts = get_read_cuts(&cigar, align_start, usize::from(region.interval().start().unwrap()), usize::from(region.interval().end().unwrap()));
        eprintln!("read_cuts: {:?}", read_cuts);
        if read_cuts.read_end == 0 {
            eprintln!("{} read_cut end is zero", name);
            // counter -= 1;
            continue;
        }
        eprintln!("----------------------------------------------------------");
        eprintln!("read_cuts: {:?}", read_cuts);
        let subseq = i_seq[read_cuts.read_start..read_cuts.read_end].to_vec();
        // let subqual: String = quality_scores_str[read_cuts.read_start..read_cuts.read_end].to_string();
        let subseq_str = String::from_utf8(subseq.to_vec()).expect("unexpected utf8 in sequence");
        // let subqual_str: String = String::from_utf8_lossy(&subqual).into_owned();
        let subseq_align_span: isize = (read_cuts.ref_end as isize).saturating_sub(read_cuts.ref_start as isize);
        eprintln!("name: {:?}", name);
        let mut hap: u8 = 0;
        match data.get(b"HP") {
            Some(Ok(value)) => hap = value.as_int().unwrap() as u8,
            Some(Err(e)) => eprintln!("Error occurred: {}", e),
            None => eprintln!("Tag not found"),
        }
        match hap {
            0 => h0_subseq_vec.push((name, subseq.clone(), ".".to_string(), read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
            // 2 => h2_subseq_vec.push((name, subseq.clone(), subqual.clone(), read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
            // 1 => h1_subseq_vec.push((name, subseq.clone(), subqual.clone(), read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
            _ => eprintln!("More than 3 hap groups detected. BladeRunner currently does not support more than diploid")
        }
        

        eprintln!("ref_id: {}", ref_id);
        eprintln!("align_start: {}", align_start);
        eprintln!("align_end: {}", align_end);
        eprintln!("map_quality: {:?}", map_quality);
        eprintln!("flags: {:?}", flags);
        eprintln!("read span: {}", span);
        eprintln!("subseq alignment span: {}", subseq_align_span);
        eprintln!("subseq relative to reference: {}", subseq_align_span.saturating_sub(subseq.len() as isize));
        // eprintln!("cigar: {:?}", cigar_string);
        eprintln!("subseq len: {}", subseq.len());
        eprintln!("subseq: {:?}", subseq_str);
        // eprintln!("subqual: {:?}", subqual);
        // eprintln!("all quality_scores_str: {:?}", quality_scores_str);
        // eprintln!("data: {:?}", data);
        eprintln!("HP tag: {:?}", hap);
        // eprintln!("subseq_vec: {:?}", subseq_vec);

    }
    // eprintln!("Number of reads in region: {}", counter);
    // eprintln!("Number of reads HAP1: {}", h1_subseq_vec.len());
    // eprintln!("Number of reads HAP2: {}", h2_subseq_vec.len());
    // eprintln!("Number of reads no HAP: {}", h0_subseq_vec.len());

    h0_subseq_vec
}