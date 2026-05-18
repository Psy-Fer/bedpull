use anyhow::{Context, Result};
use itertools::Itertools;
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::Record;
use crate::cli::Opts;

use crate::utils::{get_read_cuts, ReadCuts};
pub use crate::cigar::ToCigarOps;

// For bam reading
pub fn get_bam_reads<R>(_opts: &Opts, query: bam::io::reader::Query<R>, region: &Region, lflank: usize, rflank: usize) -> Result<Vec<(String, Vec<u8>, String, usize, usize)>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{

    let mut h0_subseq_vec: Vec<(String, Vec<u8>, String, usize, usize)> = vec![]; // no hap assigned

    for result in query.records() {
        let record = result.context("failed to read BAM record")?;
        let align_start = usize::from(
            record.alignment_start()
                .ok_or_else(|| anyhow::anyhow!("BAM record has no alignment start"))?
                .context("invalid alignment start position")?
        );
        let align_end = usize::from(
            record.alignment_end()
                .ok_or_else(|| anyhow::anyhow!("BAM record has no alignment end"))?
                .context("invalid alignment end position")?
        );
        let name_bytes: &[u8] = record.name()
            .ok_or_else(|| anyhow::anyhow!("BAM record has no name"))?
            .as_ref();
        let name = String::from_utf8(name_bytes.to_vec()).context("BAM record name contains invalid UTF-8")?;
        let seq = record.sequence();
        let i_seq = seq.iter().collect_vec();
        let i_qual =  record.quality_scores().as_ref()
                                                    .iter()
                                                    .map(|&score| (score + 33) as u8) // Adjust quality scores
                                                    .collect::<Vec<_>>();

        // eprintln!("quality_scores adjusted: {:?}", i_qual);
        // now convert that to a String
        let quality_scores_str: String = String::from_utf8_lossy(&i_qual).into_owned();
        let cigar = record.cigar().to_cigar_ops().context("invalid CIGAR in BAM record")?;
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

        let region_start = usize::from(region.interval().start().unwrap());
        let region_end = usize::from(region.interval().end().unwrap());

        // Expand window by flanks; reads overlapping anywhere in the flanked window are included.
        let desired_start = region_start.saturating_sub(lflank);
        let desired_end = region_end + rflank;

        if (align_end < desired_start) || (align_start > desired_end) {
            continue;
        }

        // Clamp the desired window to this read's alignment span so get_read_cuts stays in bounds.
        let (eff_start, eff_end) = if lflank == 0 && rflank == 0 {
            (region_start, region_end)
        } else {
            (desired_start.max(align_start), desired_end.min(align_end))
        };

        let read_cuts: ReadCuts = get_read_cuts(&cigar, align_start, eff_start, eff_end);
        // eprintln!("read_cuts: {:?}", read_cuts);
        if read_cuts.read_end == 0 {
            // eprintln!("{} read_cut end is zero", name);
            // counter -= 1;
            continue;
        }
        eprintln!("----------------------------------------------------------");
        // eprintln!("read_cuts: {:?}", read_cuts);
        let subseq = i_seq[read_cuts.read_start..read_cuts.read_end].to_vec();
        let subqual: String = quality_scores_str[read_cuts.read_start..read_cuts.read_end].to_string();
        // let subseq_str = String::from_utf8(subseq.to_vec()).expect("unexpected utf8 in sequence");
        // let subqual_str: String = String::from_utf8_lossy(&subqual).into_owned();
        // let subseq_align_span: isize = (read_cuts.ref_end as isize).saturating_sub(read_cuts.ref_start as isize);
        // eprintln!("name: {:?}", name);
        // let mut hap: u8 = 0;
        // match data.get(b"HP") {
        //     Some(Ok(value)) => hap = value.as_int().unwrap() as u8,
        //     Some(Err(e)) => eprintln!("Error occurred: {}", e),
        //     None => eprintln!("Tag not found"),
        // }
        h0_subseq_vec.push((name, subseq.clone(), subqual, read_cuts.ref_start.clone(), read_cuts.ref_end.clone()));
        // match hap {
        //     0 => h0_subseq_vec.push((name, subseq.clone(), subqual, read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
        //     // 2 => h2_subseq_vec.push((name, subseq.clone(), subqual.clone(), read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
        //     // 1 => h1_subseq_vec.push((name, subseq.clone(), subqual.clone(), read_cuts.ref_start.clone(), read_cuts.ref_end.clone())),
        //     _ => eprintln!("multiple haplotypes detected. bedpull currently does not support phased data")
        // }
        

        // eprintln!("ref_id: {}", ref_id);
        // eprintln!("align_start: {}", align_start);
        // eprintln!("align_end: {}", align_end);
        // eprintln!("map_quality: {:?}", map_quality);
        // eprintln!("flags: {:?}", flags);
        // eprintln!("read span: {}", span);
        // eprintln!("subseq alignment span: {}", subseq_align_span);
        // eprintln!("subseq relative to reference: {}", subseq_align_span.saturating_sub(subseq.len() as isize));
        // // eprintln!("cigar: {:?}", cigar_string);
        // eprintln!("subseq len: {}", subseq.len());
        // eprintln!("subseq: {:?}", subseq_str);
        // // eprintln!("subqual: {:?}", subqual);
        // // eprintln!("all quality_scores_str: {:?}", quality_scores_str);
        // // eprintln!("data: {:?}", data);
        // eprintln!("HP tag: {:?}", hap);
        // eprintln!("subseq_vec: {:?}", subseq_vec);

    }
    // eprintln!("Number of reads in region: {}", counter);
    // eprintln!("Number of reads HAP1: {}", h1_subseq_vec.len());
    // eprintln!("Number of reads HAP2: {}", h2_subseq_vec.len());
    // eprintln!("Number of reads no HAP: {}", h0_subseq_vec.len());

    Ok(h0_subseq_vec)
}
