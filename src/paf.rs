use anyhow::{Context, Result, anyhow, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::io::{Seek, SeekFrom};

/// A parsed record from a PAF (Pairwise mApping Format) alignment file.
///
/// Fields correspond 1-to-1 with the 12 mandatory PAF columns. The `cigar` field
/// is populated from the optional `cg:Z:` tag; it is `None` for alignments that
/// were run without `--cs`/`--eqx` or that do not emit a full CIGAR. bedpull
/// requires a CIGAR to perform coordinate-accurate extraction; records without
/// one are skipped with a warning.
// TODO: remove this later
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PafRecord {
    /// Name of the query sequence (PAF column 1).
    pub query_name: String,
    /// Total length of the query sequence in bp (PAF column 2).
    pub query_length: usize,
    /// 0-based start of the aligned block on the query (PAF column 3).
    pub query_start: usize,
    /// 0-based, exclusive end of the aligned block on the query (PAF column 4).
    pub query_end: usize,
    /// Strand of the alignment relative to the target: `'+'` or `'-'` (PAF column 5).
    pub strand: char,
    /// Name of the target (reference) sequence (PAF column 6).
    pub target_name: String,
    /// Total length of the target sequence in bp (PAF column 7).
    pub target_length: usize,
    /// 0-based start of the aligned block on the target (PAF column 8).
    pub target_start: usize,
    /// 0-based, exclusive end of the aligned block on the target (PAF column 9).
    pub target_end: usize,
    /// Number of matching bases in the alignment (PAF column 10).
    pub num_matches: usize,
    /// Total number of bases (including gaps) in the alignment block (PAF column 11).
    pub alignment_length: usize,
    /// Mapping quality score, 255 if unavailable (PAF column 12).
    pub mapping_quality: u8,
    /// Full CIGAR string from the `cg:Z:` optional tag, if present.
    pub cigar: Option<String>,
}

impl PafRecord {
    /// Parse a single tab-delimited PAF line into a [`PafRecord`].
    ///
    /// Returns an error if the line has fewer than 12 fields or if any
    /// mandatory numeric field cannot be parsed.
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            bail!("invalid PAF line: too few fields (got {})", fields.len());
        }

        // Extract CIGAR from cg tag ("cg:Z:" is 5 bytes)
        let cigar = fields
            .iter()
            .skip(12)
            .find(|f| f.starts_with("cg:Z:"))
            .map(|f| f[5..].to_string());

        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().context("invalid PAF query length")?,
            query_start: fields[2].parse().context("invalid PAF query start")?,
            query_end: fields[3].parse().context("invalid PAF query end")?,
            strand: fields[4]
                .chars()
                .next()
                .ok_or_else(|| anyhow!("PAF strand field is empty"))?,
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

/// A single entry in a [`PafIndex`], recording where one alignment record lives on disk.
///
/// Storing `target_start` and `target_end` alongside the byte offset allows
/// [`PafIndex::query`] to filter by genomic interval without seeking to and
/// parsing each record.
#[derive(Debug, Clone)]
pub struct PafIndexEntry {
    /// Byte offset from the start of the PAF file to the beginning of this record's line.
    pub offset: u64,
    /// 0-based start of the alignment on the target (reference) sequence.
    pub target_start: usize,
    /// 0-based, exclusive end of the alignment on the target sequence.
    pub target_end: usize,
}

/// A byte-offset index over a PAF file, partitioned by target chromosome.
///
/// The index maps each chromosome name to a list of [`PafIndexEntry`] values sorted
/// by `target_start`. Building the index requires a single linear pass over the PAF
/// file; subsequent region queries are O(n) in the number of entries on that
/// chromosome but avoid re-reading or re-parsing the entire file because only the
/// matching records are seeked to with [`read_paf_record_at_offset`].
///
/// The on-disk format is a plain TSV (`<paf>.idx`):
/// ```text
/// <chrom>\t<byte_offset>\t<target_start>\t<target_end>
/// ```
pub struct PafIndex {
    /// Per-chromosome lists of index entries, sorted by `target_start`.
    pub entries: HashMap<String, Vec<PafIndexEntry>>,
}

impl Default for PafIndex {
    fn default() -> Self {
        Self::new()
    }
}

impl PafIndex {
    pub fn new() -> Self {
        PafIndex {
            entries: HashMap::new(),
        }
    }

    /// Build a [`PafIndex`] by scanning the PAF file at `paf_path`.
    ///
    /// Reads the file line by line, recording the byte offset of each alignment
    /// record along with its target coordinates. Entries are sorted by
    /// `target_start` per chromosome before returning. Returns an error if the
    /// file cannot be opened or if target coordinate fields cannot be parsed.
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

                    index.entries.entry(target_name).or_default().push(entry);
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

    /// Persist the index to a TSV file at `index_path`.
    ///
    /// Each line has the format `<chrom>\t<offset>\t<target_start>\t<target_end>`.
    /// The file can be reloaded with [`PafIndex::load`]. Returns an error if the
    /// file cannot be created or written.
    pub fn save(&self, index_path: &str) -> Result<()> {
        let mut file = File::create(index_path)?;

        for (chrom, entries) in &self.entries {
            for entry in entries {
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}",
                    chrom, entry.offset, entry.target_start, entry.target_end
                )?;
            }
        }

        Ok(())
    }

    /// Load a previously saved index from `index_path`.
    ///
    /// Reads the TSV written by [`PafIndex::save`] and reconstructs the in-memory
    /// index. Note that the loaded index is not re-sorted; the saved file is
    /// assumed to be ordered correctly. Returns an error if the file cannot be
    /// read or if any field cannot be parsed.
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

                index.entries.entry(chrom).or_default().push(entry);
            }
        }

        Ok(index)
    }

    /// Return all index entries whose alignment overlaps `[start, end)` on `chrom`.
    ///
    /// Uses the standard half-open overlap test: `entry.target_start < end && entry.target_end > start`.
    /// Returns an empty `Vec` if the chromosome is not present in the index.
    /// The returned references are valid for the lifetime of `self`.
    pub fn query(&self, chrom: &str, start: usize, end: usize) -> Vec<&PafIndexEntry> {
        if let Some(entries) = self.entries.get(chrom) {
            entries
                .iter()
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
        let orig_hits = original.query("chr1", 120, 130);
        let load_hits = loaded.query("chr1", 120, 130);
        assert_eq!(orig_hits.len(), load_hits.len());
        assert_eq!(orig_hits[0].target_start, load_hits[0].target_start);
        assert_eq!(orig_hits[0].target_end, load_hits[0].target_end);
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

/// Seek to `offset` bytes into `paf_path` and parse the line there as a [`PafRecord`].
///
/// This is the complement of the index: given a [`PafIndexEntry::offset`] returned
/// by [`PafIndex::query`], this function retrieves the full record without scanning
/// the entire file. Returns an error if the file cannot be opened, the seek fails,
/// or the line cannot be parsed.
pub fn read_paf_record_at_offset(paf_path: &str, offset: u64) -> Result<PafRecord> {
    let mut file = File::open(paf_path)?; // can I open this once and move the seek backwards?
    file.seek(SeekFrom::Start(offset))?;

    let mut reader = BufReader::new(file);
    let mut line = String::new();
    reader.read_line(&mut line)?;

    PafRecord::from_line(&line)
}
