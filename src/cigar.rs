use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::bam::record::Cigar as BamCigar;
use noodles::sam::alignment::record::cigar::Cigar as SamCigar;

// create struct to hold a generic cigar format.
// convert any cigar format into this generic one to be analysed by get_read_cuts
#[derive(Debug, Clone)]
pub struct CigarOp {
    pub kind: Kind,
    pub len: usize,
}

pub type CigarOps = Vec<CigarOp>;


pub trait ToCigarOps {
    fn to_cigar_ops(&self) -> CigarOps;
}

// handle cigar from noodles
impl<'a> ToCigarOps for BamCigar<'a> {
    fn to_cigar_ops(&self) -> CigarOps {
        self.iter()
            .map(|op| {
                let op = op.expect("Invalid CIGAR operation");
                CigarOp {
                    kind: op.kind(),
                    len: op.len(),
                }
            })
            .collect()
    }
}

impl ToCigarOps for dyn SamCigar {
    fn to_cigar_ops(&self) -> CigarOps {
        self.iter()
            .map(|op| {
                let op = op.expect("Invalid CIGAR operation");
                CigarOp {
                    kind: op.kind(),
                    len: op.len(),
                }
            })
            .collect()
    }
}
// handle cigar from paf
impl ToCigarOps for &str {
    fn to_cigar_ops(&self) -> CigarOps {
        let mut ops = Vec::new();
        let mut num = String::new();
        
        for ch in self.chars() {
            if ch.is_numeric() {
                num.push(ch);
            } else {
                let len = num.parse::<usize>().expect("Invalid CIGAR length");
                let kind = match ch {
                    'M' => Kind::Match,
                    '=' => Kind::SequenceMatch,
                    'X' => Kind::SequenceMismatch,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    'N' => Kind::Skip,
                    'S' => Kind::SoftClip,
                    'H' => Kind::HardClip,
                    'P' => Kind::Pad,
                    _ => panic!("Unknown CIGAR operation: {}", ch),
                };
                ops.push(CigarOp { kind, len });
                num.clear();
            }
        }
        ops
    }
}