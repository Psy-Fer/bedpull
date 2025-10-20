use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;



#[derive(Debug)]
pub struct BedRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub name: Option<String>,
}


impl BedRecord {
    /// Parses a tab-delimited string into a `BedRecord`.
    ///
    /// ## Arguments
    ///
    /// * `record` - A string slice that holds a single line from the BED file.
    ///
    /// ## Returns
    ///
    /// * `Result<Self, &'static str>` - Returns an instance of `BedRecord` if successful,
    ///   or an error message if the parsing fails.
    fn from_str(record: &str) -> Result<Self, &'static str> {
        let fields: Vec<&str> = record.split('\t').collect();
        if fields.len() < 3 || fields.len() > 4 {
            return Err("Incorrect number of fields | 3 or 4 fields allowed.");
        }
        Ok(Self {
            chrom: fields[0].to_string(),
            start: fields[1].parse().map_err(|_| "Invalid start value")?,
            end: fields[2].parse().map_err(|_| "Invalid end value")?,
            name: if fields.len() > 3 && fields[3] != "." {
                Some(fields[3].to_string())
            } else {
                None
            },
        })
    }
}


/// A struct representing a custom BED file reader.
///
/// This reader is designed to parse a BED file with custom fields. Each line in the file
/// should be tab-delimited and contain at least the following fields:
/// 
/// 1. `chrom` - Chromosome name (String)
/// 2. `start` - Start position (usize)
/// 3. `end` - End position (usize)
///
/// Additionally, it can contain the following optional fields:
///
/// 4. `name` - Name (Option<String>) the name of the repeat region
/// 
/// Empty fields can be set to None with a dot: `.`
///
/// The reader validates the types of these fields when reading the file.
/// 
/// ## custom bed file format
/// 
/// chr | start | end | name (optional)
/// 
/// chr1 | 895325 | 895345 | HMNR7_VWA1
/// 
pub struct BedReader<R: BufRead> {
    reader: R,
}

impl BedReader<BufReader<File>> {
    pub fn from_path(path: &PathBuf) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Ok(BedReader { reader })
    }
}

impl<R: BufRead> Iterator for BedReader<R> {
     /// Reads the next `BedRecord` from the file.
    ///
    /// ## Returns
    ///
    /// * `Some(Result<BedRecord, &'static str>)` - Returns the next record wrapped in `Some` if successful,
    ///   or `None` if the end of the file is reached. If there's an error reading the line or parsing the record,
    ///   an error message is returned wrapped in `Some`.
    type Item = Result<BedRecord, &'static str>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => None, // EOF reached
            Ok(_) => Some(BedRecord::from_str(line.trim_end())),
            Err(_) => Some(Err("Error reading line")),
        }
    }
}