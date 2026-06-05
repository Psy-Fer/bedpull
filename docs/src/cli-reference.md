# CLI Reference

## Synopsis

```
bedpull [OPTIONS] --bed <BED> --output <OUTPUT>
```

At least one of `--bam` or (`--paf` + `--query_ref`) must be provided.

## Flags

| Flag | Short | Default | Applies to | Description |
|------|-------|---------|------------|-------------|
| `--bam <FILE>` | `-b` | — | BAM | Indexed BAM file. Requires a companion `.bam.bai` in the same directory. |
| `--reference <FILE>` | `-f` | — | — | Reference FASTA (reserved; not yet fully implemented). |
| `--bed <FILE>` | `-r` | required | both | BED file of target regions. 3- or 4-column format; the 4th column is used as the region name in output headers. |
| `--paf <FILE>` | — | — | PAF | PAF alignment file. Must include the `cg:Z:` CIGAR tag. |
| `--query_ref <FILE>` | — | — | PAF | Query FASTA to extract sequence from when using `--paf`. Must be indexed with `samtools faidx`. |
| `--output <FILE>` | `-o` | required | both | Output FASTA or FASTQ file path. |
| `--fastq` | — | false | BAM | Write FASTQ instead of FASTA. Quality scores are Phred+33 encoded. |
| `--min_mapq <N>` | — | `0` | BAM | Minimum mapping quality. Reads with MAPQ below `N` are skipped. `0` disables the filter. |
| `--min_region_quality <F>` | — | `0` | BAM + `--fastq` | Minimum mean Phred quality of the extracted subsequence. Reads below this threshold are discarded. `0` disables the filter. |
| `--include_secondary` | — | false | BAM | Include secondary alignments (SAM flag `0x100`). Excluded by default. |
| `--include_supplementary` | — | false | BAM | Include supplementary alignments (SAM flag `0x800`). Excluded by default. |
| `--partial` | — | false | BAM | Include reads that only partially overlap the BED region. By default only reads spanning the entire region are returned. |
| `--flanks <N>` | — | `0` | both | Expand the extraction window by `N` bp on both sides of each BED region before the CIGAR walk. |
| `--lflank <N>` | — | `0` | both | Left-side flank in bp. Overrides `--flanks` for the left side when non-zero. |
| `--rflank <N>` | — | `0` | both | Right-side flank in bp. Overrides `--flanks` for the right side when non-zero. |
| `--hap_split` | — | false | BAM | Split output by HP haplotype tag into `<output>.h0.<ext>`, `<output>.h1.<ext>`, `<output>.h2.<ext>`. Reads without an HP tag go to `h0`. |
| `--use_paf_index` | — | `true` | PAF | Build and/or load a byte-offset index (`<paf>.idx`). Set to `false` to rebuild in memory every run without saving. |

## Flank precedence

When both `--flanks` and a per-side flag are set, the per-side value takes precedence:

```
effective_lflank = lflank  (if lflank > 0)  else  flanks
effective_rflank = rflank  (if rflank > 0)  else  flanks
```

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| non-zero | An error was encountered; details are written to stderr. |

## Diagnostics

bedpull writes progress and diagnostic messages to **stderr** so that stdout (or the output file) contains only sequence data. Messages include region names as they are processed, overlap counts, and warnings for partial alignments or unexpected HP tag values.
