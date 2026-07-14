# CRAM Mode

CRAM mode works identically to BAM mode but reads from an indexed CRAM file. All filtering flags, output formats, haplotype splitting, flanks, and deduplication work the same way. The only difference is the input format and the optional `--reference` flag required for reference-compressed CRAMs.

## Requirements

- A CRAM file. If the `<file>.cram.crai` index is missing, bedpull builds it automatically the
  first time you run against that file — no `samtools` required.
- For reference-compressed CRAMs: the reference FASTA used during compression. If its `.fai`
  index is missing, it's auto-built the same way.
- For CRAMs with embedded sequences (created with `samtools view -C` from a FASTA source with no external reference, or with `embed_ref=2`): no reference is needed.

## Basic usage

```bash
bedpull \
    --cram reads.cram \
    --bed regions.bed \
    --output extracted.fasta
```

## Reference-compressed CRAMs

If the CRAM was compressed against an external reference (the default for most pipelines), provide the same reference FASTA so bedpull can decode the sequences:

```bash
bedpull \
    --cram reads.cram \
    --reference GRCh38.fa \
    --bed regions.bed \
    --output extracted.fasta
```

If the reference FASTA's `.fai` index is missing, bedpull builds it automatically before extraction.

## FASTQ output

```bash
bedpull \
    --cram reads.cram \
    --reference GRCh38.fa \
    --bed regions.bed \
    --output extracted.fastq \
    --fastq
```

## Filtering reads

CRAM mode supports the same read filters as BAM mode:

```bash
# Mapping quality filter
bedpull --cram reads.cram --bed regions.bed --output out.fasta \
    --min_mapq 20

# Include secondary and supplementary alignments
bedpull --cram reads.cram --bed regions.bed --output out.fasta \
    --include_secondary \
    --include_supplementary

# Mean quality of the extracted region (works regardless of --fastq)
bedpull --cram reads.cram --bed regions.bed --output out.fastq \
    --fastq \
    --min_region_quality 20

# Include partial overlaps
bedpull --cram reads.cram --bed regions.bed --output out.fasta \
    --partial

# Include partial overlaps, but require at least 90% coverage of the requested window
# (see BAM mode's "Minimum partial coverage" section for details)
bedpull --cram reads.cram --bed regions.bed --output out.fasta \
    --partial \
    --min_partial_coverage 0.9
```

## Flanks

```bash
bedpull --cram reads.cram --bed regions.bed --output out.fasta \
    --flanks 500
```

## Haplotype splitting

If reads carry the `HP` aux tag, split output by haplotype:

```bash
bedpull \
    --cram reads.cram \
    --bed regions.bed \
    --output out.fasta \
    --hap_split
```

This creates `out.h0.fasta`, `out.h1.fasta`, and `out.h2.fasta`. `--hap_split` requires a file `--output` (not stdout).

## Deduplication

```bash
bedpull --cram reads.cram --bed regions.bed --output out.fasta --dedup
```

## Output header format

Identical to BAM mode:

```
>read_name|chr:region_start-region_end|region_name
```

When `--hap_split` is active and the read carries an HP tag:

```
>read_name|chr:region_start-region_end|region_name|h1
```

With `--partial`, a read that doesn't fully span the requested window gets a
`|missing_left=Nbp`/`|missing_right=Nbp` suffix — see the [BAM mode](./bam-mode.md#output-header-format) page for details.

## Diagnostics

By default bedpull only prints the mode banner and a completion line. Pass `--debug` to see
verbose per-region diagnostics: the parsed CLI args, each BED record as it's read, and a banner
for every region as it's processed.

## stdout

Omit `--output` (or pass `-`) to write to stdout:

```bash
bedpull --cram reads.cram --reference ref.fa --bed regions.bed | seqkit stats
```

## Converting a BAM to an embedded-reference CRAM

If you do not have the original reference available at runtime, create a self-contained CRAM with embedded sequences:

```bash
samtools view -C --no-PG reads.bam -o reads.cram
```

The `samtools index` step is optional — bedpull builds the `.crai` itself on first use if it's missing.

The resulting CRAM stores the decoded sequences internally; no `--reference` flag is needed.
