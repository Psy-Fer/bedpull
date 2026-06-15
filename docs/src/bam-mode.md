# BAM Mode

BAM mode queries an indexed BAM file for reads that overlap each BED region, walks each read's CIGAR string to calculate exact read-coordinate boundaries for the region, and extracts the corresponding subsequence (and optionally the quality string).

## Basic usage

```bash
bedpull \
    --bam reads.bam \
    --bed regions.bed \
    --output extracted.fasta
```

The BAM must have a companion `.bam.bai` index in the same directory. The output is FASTA by default; use `--fastq` to include quality scores.

Omit `--output` (or pass `-`) to write to stdout:

```bash
bedpull --bam reads.bam --bed regions.bed | seqkit stats
```

## FASTQ output

```bash
bedpull \
    --bam reads.bam \
    --bed regions.bed \
    --output extracted.fastq \
    --fastq
```

## Filtering reads

### Mapping quality

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --min_mapq 20
```

Reads with MAPQ below `20` are skipped before any CIGAR processing. The default is `0` (all reads pass).

### Secondary and supplementary alignments

By default only primary alignments are returned. Include secondary or supplementary alignments explicitly:

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --include_secondary \
    --include_supplementary
```

### Quality of the extracted region

When using FASTQ output you can filter on the mean Phred score of the extracted slice:

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fastq \
    --fastq \
    --min_region_quality 20
```

The mean Phred score is computed in error-probability space (not a simple arithmetic mean of the Phred values), which is the statistically correct approach.

## Partial overlaps

By default only reads whose alignment spans the entire BED region are returned (both the left and right boundaries are crossed). To include reads that overlap only part of the region:

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --partial
```

Partial reads are clipped to whatever portion of the region they cover. Reads that start after `region_start` are extracted from read position 0; reads that end before `region_end` are extracted to the end of the read.

## Flanks

To extend the extraction window beyond the BED coordinates — for example to capture an insertion that sits just outside the annotated boundary — use `--flanks`, `--lflank`, or `--rflank`:

```bash
# Symmetric 500 bp flank on both sides
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --flanks 500

# Asymmetric: 1000 bp on the left, 200 bp on the right
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --lflank 1000 \
    --rflank 200
```

Per-side values (`--lflank`, `--rflank`) take precedence over `--flanks` for their respective side. Flanks are applied in reference space before the CIGAR walk, so any insertions within the extended window are captured.

## Haplotype splitting

If reads carry the `HP` aux tag (set by `whatshap haplotag` or similar), you can split the output by haplotype:

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fasta \
    --hap_split
```

This creates three output files alongside the specified `--output` path:

| File | Contents |
|------|----------|
| `out.h0.fasta` | Reads with no HP tag (unphased) |
| `out.h1.fasta` | Reads with HP tag = 1 |
| `out.h2.fasta` | Reads with HP tag = 2 |

## Deduplication

A read whose alignment spans two BED regions will appear in the output once per region by default, which can be useful when you want region-specific context. To emit each read only once across all regions:

```bash
bedpull --bam reads.bam --bed regions.bed --output out.fasta --dedup
```

## Output header format

Each FASTA/FASTQ record has a header of the form:

```
>read_name|chr:region_start-region_end|region_name
```

When `--hap_split` is active and the read carries an HP tag, a haplotype suffix is appended:

```
>read_name|chr:region_start-region_end|region_name|h1
```

`region_name` comes from the 4th column of the BED file, or defaults to `chr:start-end` if absent.

## Example: RFC1 VNTR from a long-read BAM

```bash
bedpull \
    --bam HG002.GRCh38.bam \
    --bed rfc1_vntr.bed \
    --output rfc1_reads.fasta \
    --min_mapq 20 \
    --flanks 200
```

Where `rfc1_vntr.bed` contains:
```
chr4	39318077	39318136	RFC1_VNTR
```
