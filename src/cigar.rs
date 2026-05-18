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

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> CigarOps {
        s.to_cigar_ops()
    }

    #[test]
    fn single_match() {
        let ops = parse("10M");
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].kind, Kind::Match);
        assert_eq!(ops[0].len, 10);
    }

    #[test]
    fn all_op_types() {
        let ops = parse("1M2I3D4N5S6H7P8=9X");
        let expected = [
            (Kind::Match, 1),
            (Kind::Insertion, 2),
            (Kind::Deletion, 3),
            (Kind::Skip, 4),
            (Kind::SoftClip, 5),
            (Kind::HardClip, 6),
            (Kind::Pad, 7),
            (Kind::SequenceMatch, 8),
            (Kind::SequenceMismatch, 9),
        ];
        assert_eq!(ops.len(), expected.len());
        for (op, (kind, len)) in ops.iter().zip(expected.iter()) {
            assert_eq!(op.kind, *kind);
            assert_eq!(op.len, *len);
        }
    }

    #[test]
    fn multi_op_mixed() {
        let ops = parse("5M2I3D");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].len, 5);
        assert_eq!(ops[1].len, 2);
        assert_eq!(ops[2].len, 3);
    }

    #[test]
    fn large_lengths() {
        let ops = parse("100000M999I");
        assert_eq!(ops[0].len, 100000);
        assert_eq!(ops[1].len, 999);
    }

    #[test]
    #[should_panic]
    fn unknown_op_panics() {
        parse("5Z");
    }
}