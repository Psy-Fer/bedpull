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

impl<R: BufRead> BedReader<R> {
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn from_reader(reader: R) -> Self {
        BedReader { reader }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn parse(s: &str) -> Result<BedRecord, &'static str> {
        BedRecord::from_str(s)
    }

    // --- BedRecord::from_str ---

    #[test]
    fn three_column_no_name() {
        let r = parse("chr1\t100\t200").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert!(r.name.is_none());
    }

    #[test]
    fn four_column_with_name() {
        let r = parse("chr4\t39318077\t39318136\tRFC1").unwrap();
        assert_eq!(r.chrom, "chr4");
        assert_eq!(r.start, 39318077);
        assert_eq!(r.end, 39318136);
        assert_eq!(r.name.as_deref(), Some("RFC1"));
    }

    #[test]
    fn dot_name_becomes_none() {
        let r = parse("chr1\t0\t100\t.").unwrap();
        assert!(r.name.is_none());
    }

    #[test]
    fn too_few_fields_is_error() {
        assert!(parse("chr1\t100").is_err());
    }

    #[test]
    fn too_many_fields_is_error() {
        assert!(parse("chr1\t100\t200\tname\textra").is_err());
    }

    #[test]
    fn invalid_start_is_error() {
        assert!(parse("chr1\tabc\t200").is_err());
    }

    #[test]
    fn invalid_end_is_error() {
        assert!(parse("chr1\t100\tabc").is_err());
    }

    // --- BedReader iterator ---

    fn reader_from(s: &str) -> BedReader<BufReader<&[u8]>> {
        BedReader::from_reader(BufReader::new(s.as_bytes()))
    }

    #[test]
    fn reader_yields_multiple_records() {
        let data = "chr1\t100\t200\nchr2\t300\t400\tregion2\n";
        let records: Vec<_> = reader_from(data).collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].as_ref().unwrap().chrom, "chr1");
        assert_eq!(records[1].as_ref().unwrap().chrom, "chr2");
        assert_eq!(records[1].as_ref().unwrap().name.as_deref(), Some("region2"));
    }

    #[test]
    fn reader_empty_input_yields_nothing() {
        let records: Vec<_> = reader_from("").collect();
        assert!(records.is_empty());
    }

    #[test]
    fn reader_propagates_parse_error() {
        let data = "chr1\tbad\t200\n";
        let records: Vec<_> = reader_from(data).collect();
        assert_eq!(records.len(), 1);
        assert!(records[0].is_err());
    }
}