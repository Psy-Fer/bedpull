use std::f64;
use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use noodles::sam::alignment::record::cigar::op::Kind;
use bio::alignment::poa::Aligner as poAligner;
use bio::alignment::pairwise::Scoring;

use noodles::core::{Region, region::Interval, Position};

use crate::bed::BedReader;
use crate::cli::Opts;
use crate::cigar::CigarOps;


/// Calculate the mean Qscore from a Qstring.
pub fn _calculate_qscore(qstring: &String) -> f64 {
    // Convert phred back to ASCII values and adjust -33
    let qs: Vec<f64> = qstring
        .chars()
        .map(|c| (c as u8) as f64 - 33.0)
        .collect();

    // Calculate mean error
    let mean_err: f64 = qs
        .iter()
        .map(|&q| (-q * f64::consts::LN_10 / 10.0).exp())
        .sum::<f64>()
        / qs.len() as f64;

    // Calculate Qscore
    let score: f64 = -10.0 * (mean_err.max(1e-4)).log10();

    score
}
/// ReadCuts is a struct that holds read subsequence positions
/// use this to cut out a subsequence of a read
#[derive(Debug, Clone)]
pub struct ReadCuts {
    pub read_start: usize,
    pub read_end: usize,
    pub ref_start: usize, // reference position
    pub ref_end: usize,
}

/// get_read_cuts takes a cigar string, alignment_start, target region, and padding amount
/// It returns a ReadCuts struct with the read_start and read_end
/// This can be used to cut out a substring of the read sequence matching the target region
/// Assumes read has full coverage of region
/// TODO: Allow for partial overlaps with flags.
pub fn get_read_cuts(cigar_ops: &CigarOps, align_start: usize,region_start: usize, region_end: usize) -> ReadCuts {
    let mut start: usize = 0;
    let mut end: usize = 0;
    let mut r_start: usize = 0;
    let mut r_end: usize = 0;
    let mut pos: usize = 0;
    let mut ref_pos: usize = align_start;
    
    let ref_start = region_start;
    let ref_end = region_end;

    for op in cigar_ops {
        // let op = op.expect("op code didn't work");
        match op.kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if (ref_pos + op.len >= ref_start) || (ref_pos + op.len >= ref_end) {
                    if (ref_pos + op.len == ref_start) || (ref_pos + op.len == ref_end) {
                        ref_pos += op.len;
                        pos += op.len;
                        if start > 0 {
                            end = pos;
                            r_end = ref_pos;
                            break;
                        } else {
                            start = pos;
                            r_start = ref_pos;
                        }
                    } else {
                        for _ in 0..op.len {
                            ref_pos += 1;
                            pos += 1;
                            if (ref_pos == ref_start) || (ref_pos == ref_end) {
                                if start > 0 {
                                    end = pos;
                                    r_end = ref_pos;
                                    break;
                                } else {
                                    start = pos;
                                    r_start = ref_pos;
                                }
                            }
                        }
                    }
                } else {
                    ref_pos += op.len;
                    pos += op.len;
                }
            },
            Kind::Insertion | Kind::SoftClip => {
                pos += op.len;
            }
            Kind::Deletion | Kind::Skip => {
                if (ref_pos + op.len >= ref_start) || (ref_pos + op.len >= ref_end) {
                    if (ref_pos + op.len == ref_start) || (ref_pos + op.len == ref_end) {
                        ref_pos += op.len;
                        if start > 0 {
                            end = pos;
                            r_end = ref_pos;
                            break;
                        } else {
                            start = pos;
                            r_start = ref_pos;
                        }
                    } else {
                        for _ in 0..op.len {
                            ref_pos += 1;
                            if (ref_pos == ref_start) || (ref_pos == ref_end) {
                                if start > 0 {
                                    end = pos;
                                    r_end = ref_pos;
                                    break;
                                } else {
                                    start = pos;
                                    r_start = ref_pos;
                                }
                            }
                        }
                    }
                } else {
                    ref_pos += op.len;
                }
            },
            Kind::HardClip | Kind::Pad => {continue;},
        }
    }
    
    ReadCuts { read_start: start, read_end: end, ref_start: r_start, ref_end: r_end }
}


pub fn _get_consensus(reads: &Vec<Vec<u8>>) -> Vec<u8>{

    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    // use first sequence as the reference
    let first_read = &reads[0];
    let mut aligner = poAligner::new(scoring,  first_read);
    for read in reads[1..].iter(){
        // add all other reads to graph
        aligner.global(read).add_to_graph();
    }
    
    // get consensus
    let consensus: Vec<u8> = aligner.consensus();

    consensus
}


pub fn read_bed(opts: &Opts,) -> Vec<(Region, String, String)> {
    let mut regions: Vec<(Region, String, String)> = vec![];
    let reader = BedReader::from_path(&opts.bed).expect("failed to read bed file");
    for record in reader {
        match record {
            Ok(record) => {
                eprintln!("{:?}", record);
                let chr: String = record.chrom.clone();
                let start = Position::try_from(record.start).expect("Couldn't get start position");
                let end = Position::try_from(record.end).expect("Couldn't get end position");
                let interval = Interval::from(start..=end);
                eprintln!("region: {:?}", Region::new(record.chrom.clone(), interval));
                let name = record.name.unwrap_or(format!("{}:{}-{}", record.chrom, start, end));
                regions.push((Region::new(record.chrom, interval), name, chr));
            },
            Err(e) => eprintln!("Error: {}", e),
        }
    }
    regions
    // println!("{:?}", regions);
}


// write a fasta record
pub fn write_fasta_record(writer: &mut BufWriter<File>, header: &str, sequence: &str) -> Result<()> {
    writeln!(writer, ">{}", header)?;
    writeln!(writer, "{}", sequence)?;
    Ok(())
}

// write a fastq record
pub fn _write_fastq_record(writer: &mut BufWriter<File>, header: &str, sequence: &str, quality: &str) -> Result<()> {
    writeln!(writer, "@{}", header)?;
    writeln!(writer, "{}", sequence)?;
    writeln!(writer, "+")?;
    writeln!(writer, "{}", quality)?;
    Ok(())
}