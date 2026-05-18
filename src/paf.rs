use anyhow::{anyhow, bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;
use std::io::{Seek, SeekFrom};


// TODO: remove this later
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub target_name: String,
    pub target_length: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub num_matches: usize,
    pub alignment_length: usize,
    pub mapping_quality: u8,
    pub cigar: Option<String>,  // cg:Z: tag
}

impl PafRecord {
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            bail!("invalid PAF line: too few fields (got {})", fields.len());
        }

        // Extract CIGAR from cg tag ("cg:Z:" is 5 bytes)
        let cigar = fields.iter()
            .skip(12)
            .find(|f| f.starts_with("cg:Z:"))
            .map(|f| f[5..].to_string());

        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().context("invalid PAF query length")?,
            query_start: fields[2].parse().context("invalid PAF query start")?,
            query_end: fields[3].parse().context("invalid PAF query end")?,
            strand: fields[4].chars().next().ok_or_else(|| anyhow!("PAF strand field is empty"))?,
            target_name: fields[5].to_string(),
            target_length: fields[6].parse().context("invalid PAF target length")?,
            target_start: fields[7].parse().context("invalid PAF target start")?,
            target_end: fields[8].parse().context("invalid PAF target end")?,
            num_matches: fields[9].parse().context("invalid PAF num matches")?,
            alignment_length: fields[10].parse().context("invalid PAF alignment length")?,
            mapping_quality: fields[11].parse().context("invalid PAF mapping quality")?,
            cigar,
        })
    }
    
    // pub fn overlaps_region(&self, chrom: &str, start: usize, end: usize) -> bool {
    //     self.target_name == chrom 
    //         && self.target_start < end 
    //         && self.target_end > start
    // }
}

// paf indexing
// Store file offset for each alignment
#[derive(Debug, Clone)]
pub struct PafIndexEntry {
    pub offset: u64,           // file byte offset
    pub target_start: usize,   // overlap checks can use these
    pub target_end: usize,
}

pub struct PafIndex {
    // chromosome: sorted list of index entries
    pub entries: HashMap<String, Vec<PafIndexEntry>>,
}

impl PafIndex {
    pub fn new() -> Self {
        PafIndex {
            entries: HashMap::new(),
        }
    }
    
    // Build index from PAF file
    pub fn build(paf_path: &str) -> Result<Self> {
        let file = File::open(paf_path)?;
        let mut reader = BufReader::new(file);
        let mut index = PafIndex::new();
        
        let mut offset: u64 = 0;
        let mut line = String::new();
        
        while reader.read_line(&mut line)? > 0 {
            let line_len = line.len() as u64;
            
            if !line.starts_with('#') && !line.is_empty() {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 9 {
                    let target_name = fields[5].to_string();
                    let target_start: usize = fields[7].parse()?;
                    let target_end: usize = fields[8].parse()?;
                    
                    let entry = PafIndexEntry {
                        offset,
                        target_start,
                        target_end,
                    };
                    
                    index.entries
                        .entry(target_name)
                        .or_insert_with(Vec::new)
                        .push(entry);
                }
            }
            
            offset += line_len;
            line.clear();
        }
        
        // Sort entries by start position for each chromosome
        for entries in index.entries.values_mut() {
            entries.sort_by_key(|e| e.target_start);
        }
        
        Ok(index)
    }
    
    // Save index to file
    pub fn save(&self, index_path: &str) -> Result<()> {
        let mut file = File::create(index_path)?;
        
        for (chrom, entries) in &self.entries {
            for entry in entries {
                writeln!(file, "{}\t{}\t{}\t{}", 
                    chrom, entry.offset, entry.target_start, entry.target_end)?;
            }
        }
        
        Ok(())
    }
    
    // Load index from file
    pub fn load(index_path: &str) -> Result<Self> {
        let file = File::open(index_path)?;
        let reader = BufReader::new(file);
        let mut index = PafIndex::new();
        
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let chrom = fields[0].to_string();
                let entry = PafIndexEntry {
                    offset: fields[1].parse()?,
                    target_start: fields[2].parse()?,
                    target_end: fields[3].parse()?,
                };
                
                index.entries
                    .entry(chrom)
                    .or_insert_with(Vec::new)
                    .push(entry);
            }
        }
        
        Ok(index)
    }
    
    // Query index for overlapping entries
    pub fn query(&self, chrom: &str, start: usize, end: usize) -> Vec<&PafIndexEntry> {
        if let Some(entries) = self.entries.get(chrom) {
            entries.iter()
                .filter(|e| e.target_start < end && e.target_end > start)
                .collect()
        } else {
            Vec::new()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const PAF_LINE_A: &str = "q1\t100\t0\t50\t+\tchr1\t1000\t100\t150\t50\t50\t60\tcg:Z:50M";
    const PAF_LINE_B: &str = "q2\t200\t0\t80\t+\tchr1\t1000\t200\t280\t80\t80\t60\tcg:Z:80M";
    const PAF_LINE_NO_CIGAR: &str = "q3\t100\t0\t50\t+\tchr2\t1000\t100\t150\t50\t50\t60";

    // --- PafRecord::from_line ---

    #[test]
    fn parse_record_with_cigar() {
        let r = PafRecord::from_line(PAF_LINE_A).unwrap();
        assert_eq!(r.query_name, "q1");
        assert_eq!(r.query_length, 100);
        assert_eq!(r.query_start, 0);
        assert_eq!(r.query_end, 50);
        assert_eq!(r.strand, '+');
        assert_eq!(r.target_name, "chr1");
        assert_eq!(r.target_start, 100);
        assert_eq!(r.target_end, 150);
        assert_eq!(r.mapping_quality, 60);
        assert_eq!(r.cigar.as_deref(), Some("50M"));
    }

    #[test]
    fn parse_record_without_cigar() {
        let r = PafRecord::from_line(PAF_LINE_NO_CIGAR).unwrap();
        assert!(r.cigar.is_none());
    }

    #[test]
    fn parse_record_too_few_fields_is_error() {
        assert!(PafRecord::from_line("q1\t100\t0").is_err());
    }

    // --- PafIndex::build + query ---

    fn temp_paf(contents: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, "{}", contents).unwrap();
        f
    }

    #[test]
    fn build_and_query_overlapping() {
        let paf = temp_paf(&format!("{}\n{}\n", PAF_LINE_A, PAF_LINE_B));
        let idx = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
        // A: chr1 100-150, B: chr1 200-280
        let hits = idx.query("chr1", 120, 130);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].target_start, 100);
    }

    #[test]
    fn query_no_overlap_returns_empty() {
        let paf = temp_paf(&format!("{}\n", PAF_LINE_A));
        let idx = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
        assert!(idx.query("chr1", 200, 300).is_empty());
    }

    #[test]
    fn query_unknown_chrom_returns_empty() {
        let paf = temp_paf(&format!("{}\n", PAF_LINE_A));
        let idx = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
        assert!(idx.query("chrX", 100, 150).is_empty());
    }

    #[test]
    fn query_boundary_exclusive() {
        // target_end=150; query start=150 should NOT overlap (filter: target_end > start)
        let paf = temp_paf(&format!("{}\n", PAF_LINE_A));
        let idx = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
        assert!(idx.query("chr1", 150, 200).is_empty());
    }

    // --- save / load round-trip ---

    #[test]
    fn save_and_load_round_trip() {
        let paf = temp_paf(&format!("{}\n{}\n", PAF_LINE_A, PAF_LINE_B));
        let idx_file = tempfile::NamedTempFile::new().unwrap();
        let idx_path = idx_file.path().to_str().unwrap();

        let original = PafIndex::build(paf.path().to_str().unwrap()).unwrap();
        original.save(idx_path).unwrap();

        let loaded = PafIndex::load(idx_path).unwrap();
        let orig_hits  = original.query("chr1", 120, 130);
        let load_hits  = loaded.query("chr1", 120, 130);
        assert_eq!(orig_hits.len(), load_hits.len());
        assert_eq!(orig_hits[0].target_start, load_hits[0].target_start);
        assert_eq!(orig_hits[0].target_end,   load_hits[0].target_end);
    }

    // --- read_paf_record_at_offset ---

    #[test]
    fn read_record_at_offset_zero() {
        let paf = temp_paf(&format!("{}\n{}\n", PAF_LINE_A, PAF_LINE_B));
        let r = read_paf_record_at_offset(paf.path().to_str().unwrap(), 0).unwrap();
        assert_eq!(r.query_name, "q1");
    }

    #[test]
    fn read_record_at_second_line_offset() {
        let first_line = format!("{}\n", PAF_LINE_A);
        let offset = first_line.len() as u64;
        let paf = temp_paf(&format!("{}{}\n", first_line, PAF_LINE_B));
        let r = read_paf_record_at_offset(paf.path().to_str().unwrap(), offset).unwrap();
        assert_eq!(r.query_name, "q2");
    }
}

// get the paf record using an index offset
pub fn read_paf_record_at_offset(paf_path: &str, offset: u64) -> Result<PafRecord> {
    let mut file = File::open(paf_path)?; // can I open this once and move the seek backwards?
    file.seek(SeekFrom::Start(offset))?;
    
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    reader.read_line(&mut line)?;
    
    PafRecord::from_line(&line)
}