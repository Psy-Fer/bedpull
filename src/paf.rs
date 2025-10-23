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
    pub fn from_line(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() < 12 {
            return Err("Invalid PAF line: too few fields".into());
        }
        
        // Extract CIGAR from cg tag
        let cigar = fields.iter()
            .skip(12)
            .find(|f| f.starts_with("cg:Z:"))
            .map(|f| f.strip_prefix("cg:Z:").unwrap().to_string());
        
        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse()?,
            query_start: fields[2].parse()?,
            query_end: fields[3].parse()?,
            strand: fields[4].chars().next().unwrap(),
            target_name: fields[5].to_string(),
            target_length: fields[6].parse()?,
            target_start: fields[7].parse()?,
            target_end: fields[8].parse()?,
            num_matches: fields[9].parse()?,
            alignment_length: fields[10].parse()?,
            mapping_quality: fields[11].parse()?,
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
    pub fn build(paf_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
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
    pub fn save(&self, index_path: &str) -> Result<(), Box<dyn std::error::Error>> {
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
    pub fn load(index_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
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

// get the paf record using an index offset
pub fn read_paf_record_at_offset(paf_path: &str, offset: u64) -> Result<PafRecord, Box<dyn std::error::Error>> {
    let mut file = File::open(paf_path)?; // can I open this once and move the seek backwards?
    file.seek(SeekFrom::Start(offset))?;
    
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    reader.read_line(&mut line)?;
    
    PafRecord::from_line(&line)
}