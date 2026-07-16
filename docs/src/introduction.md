# Introduction

`bedpull` extracts sequences from BAM, CRAM, or PAF alignment files using BED coordinates. It is designed for situations where standard coordinate-lifting tools (`samtools faidx`, `bedtools getfasta`, `liftOver`) give you the wrong answer because the region of interest contains structural variation relative to the reference.

## The problem with reference-only extraction

BED coordinates are defined in reference space. A tool that naively slices a reference FASTA at those coordinates returns only the bases that exist in the reference. A tool that lifts BED intervals to query coordinates using only the alignment start/end positions makes the same mistake at the alignment level: it maps a reference interval to a query interval by applying a simple offset, which is correct for pure matches but wrong whenever the CIGAR string contains insertions.

An **insertion** in a CIGAR string means the read (or assembly) carries bases that have no reference coordinate at all. Those bases sit between two match blocks in read space, but they are invisible in reference space. Any tool that converts reference coordinates to read coordinates using only the alignment's start position will produce a slice that skips over the inserted bases entirely.

`bedpull` avoids this by walking the CIGAR string one operation at a time. It tracks both the reference position and the read position simultaneously as it steps through each operation. When it reaches the start of a BED region it records the read position, continues walking (accumulating read positions through any insertions it encounters), and when it reaches the end of the BED region it records the read position again. The slice `seq[read_start..read_end]` includes every base (matched, mismatched, and inserted) that the alignment places between those two reference coordinates.

## A concrete example: RFC1 VNTR

The RFC1 locus on chromosome 4 contains a variable-number tandem repeat (VNTR) that is relevant to spinocerebellar ataxia. A BED region spanning the annotated RFC1 repeat is 59 bp in reference coordinates (`chr4:39318077-39318136` on hg38/hs1). A long-read assembly of a pathological haplotype aligns to this region with a 520 bp insertion inside those 59 reference bases; the VNTR expansion is carried entirely in the insertion block.

Without CIGAR-aware extraction you get 59 bp of flanking sequence and miss the expansion entirely. With `bedpull`, walking the CIGAR correctly produces a 579 bp extracted sequence (59 reference-matching bases plus 520 inserted bases). This is the actual unit of biological interest, and it cannot be obtained from any tool that works purely in reference coordinates.

## When to use bedpull

Use `bedpull` when:

- Your target regions may contain insertions relative to the reference (VNTRs, STRs, mobile element insertions, structural variants).
- You are working with long-read sequencing data aligned to a reference and need the full read sequence over a locus, not just the reference-matching bases.
- You have an assembly-to-reference PAF alignment and need to extract the corresponding assembly sequence for a set of reference BED regions.

If your regions are simple SNP/indel loci with no large insertions, standard tools will give equivalent results. `bedpull` adds value specifically when the query sequence diverges structurally from the reference within the BED interval.
