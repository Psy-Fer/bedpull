use anyhow::{Context, Result, bail};
use noodles::bam::record::Cigar as BamCigar;
use noodles::sam::alignment::record::cigar::Cigar as SamCigar;
use noodles::sam::alignment::record::cigar::op::Kind;

#[derive(Debug, Clone)]
pub struct CigarOp {
    pub kind: Kind,
    pub len: usize,
}

pub type CigarOps = Vec<CigarOp>;

pub trait ToCigarOps {
    fn to_cigar_ops(&self) -> Result<CigarOps>;
}

impl<'a> ToCigarOps for BamCigar<'a> {
    fn to_cigar_ops(&self) -> Result<CigarOps> {
        self.iter()
            .map(|op| {
                let op = op.context("invalid BAM CIGAR operation")?;
                Ok(CigarOp {
                    kind: op.kind(),
                    len: op.len(),
                })
            })
            .collect()
    }
}

impl ToCigarOps for dyn SamCigar {
    fn to_cigar_ops(&self) -> Result<CigarOps> {
        self.iter()
            .map(|op| {
                let op = op.context("invalid SAM CIGAR operation")?;
                Ok(CigarOp {
                    kind: op.kind(),
                    len: op.len(),
                })
            })
            .collect()
    }
}

impl ToCigarOps for &str {
    fn to_cigar_ops(&self) -> Result<CigarOps> {
        let mut ops = Vec::new();
        let mut num = String::new();

        for ch in self.chars() {
            if ch.is_numeric() {
                num.push(ch);
            } else {
                let len = num.parse::<usize>().context("invalid CIGAR length")?;
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
                    _ => bail!("unknown CIGAR operation: {}", ch),
                };
                ops.push(CigarOp { kind, len });
                num.clear();
            }
        }
        Ok(ops)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> CigarOps {
        s.to_cigar_ops().expect("test CIGAR should be valid")
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
    fn unknown_op_returns_error() {
        assert!("5Z".to_cigar_ops().is_err());
    }
}
