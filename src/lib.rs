//! # bedpull
//!
//! CIGAR-aware sequence extraction from BAM and PAF alignment files using BED coordinates.
//!
//! Standard coordinate-lifting tools (liftOver, `samtools faidx`) work purely in reference
//! space: an inserted sequence has no reference coordinate, so it is silently dropped when
//! you lift a BED interval. `bedpull` avoids this by walking the CIGAR string of each
//! alignment, tracking both the reference position and the read position simultaneously.
//! Insertions that fall inside — or immediately adjacent to — a BED region are captured
//! automatically.
//!
//! ## Extraction modes
//!
//! | Mode | Required flags | Description |
//! |------|----------------|-------------|
//! | BAM  | `--bam`        | Query an indexed BAM for reads overlapping each BED region. |
//! | PAF  | `--paf --query_ref` | Build/load a byte-offset index of a PAF file, find overlapping alignments, extract from the query FASTA. |
//!
//! ## Library usage
//!
//! ```no_run
//! use bedpull::{BamConfig, PafIndex, get_read_cuts};
//! ```
//!
//! See the individual module docs for full API details.

pub mod bed;
pub mod cigar;
pub mod paf;
pub mod reads;
pub mod utils;

pub use cigar::{CigarOp, CigarOps, ToCigarOps};
pub use paf::{PafIndex, PafIndexEntry, PafRecord, read_paf_record_at_offset};
pub use reads::{BamConfig, BamRead, get_bam_reads};
pub use utils::{
    ReadCuts, calculate_qscore, extract_from_fasta_coords, get_read_cuts, read_bed, revcomp,
    write_fasta_record, write_fastq_record,
};
