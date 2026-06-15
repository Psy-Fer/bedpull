# CLI Reference

## Synopsis

```
bedpull [OPTIONS] --bed <BED>
```

At least one input mode must be specified: `--bam`, `--cram`, or `--paf` + `--query_ref`. `--bam` and `--cram` are mutually exclusive. `--paf` cannot be combined with `--bam` or `--cram`.

## Input flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--bam <FILE>` | `-b` | — | BAM | Indexed BAM file. Requires a companion `.bam.bai` or `.bai` index. |
| `--cram <FILE>` | — | — | CRAM | Indexed CRAM file. Requires a companion `.cram.crai` index. |
| `--reference <FILE>` | `-f` | — | CRAM | Reference FASTA used during CRAM compression. Required for reference-compressed CRAMs; not needed for CRAMs with embedded sequences. Must be indexed with `samtools faidx`. |
| `--paf <FILE>` | — | — | PAF | PAF alignment file. Must include the `cg:Z:` CIGAR tag (`minimap2 -c`). Must be used together with `--query_ref`. |
| `--query_ref <FILE>` | — | — | PAF | Query FASTA to extract sequence from. Must be indexed with `samtools faidx`. |
| `--bed <FILE>` | `-r` | required | all | BED file of target regions. 3- or 4-column format; the 4th column is used as the region name in output headers. |

## Output flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--output <FILE>` | `-o` | `-` (stdout) | all | Output FASTA or FASTQ file. Pass `-` to write to stdout. |
| `--bed_out <FILE>` | — | — | PAF | Write lifted-over query coordinates as BED6 alongside the FASTA output. Pass `-` for stdout (but not simultaneously with `--output -`). |
| `--fastq` | — | false | BAM, CRAM | Write FASTQ instead of FASTA. Quality scores are Phred+33 encoded. |

## Filtering flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--min_mapq <N>` | — | `0` | BAM, CRAM | Minimum mapping quality. Reads with MAPQ below `N` are skipped. `0` disables the filter. Reads with unavailable MAPQ (255) always pass. |
| `--min_region_quality <F>` | — | `0` | BAM, CRAM | Minimum mean Phred quality of the extracted subsequence. Reads below this threshold are discarded. `0` disables the filter. Computed in error-probability space. |
| `--include_secondary` | — | false | BAM, CRAM | Include secondary alignments (SAM flag `0x100`). Excluded by default. |
| `--include_supplementary` | — | false | BAM, CRAM | Include supplementary alignments (SAM flag `0x800`). Excluded by default. |
| `--partial` | — | false | BAM, CRAM | Include reads that only partially overlap the BED region. By default only reads spanning the entire region are returned. |
| `--dedup` | — | false | all | Deduplicate output: if the same read or contig name is seen more than once across BED regions, emit it only the first time. |

## Extraction flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--flanks <N>` | — | `0` | all | Expand the extraction window by `N` bp on both sides of each BED region before the CIGAR walk. |
| `--lflank <N>` | — | `0` | all | Left-side flank in bp. Overrides `--flanks` for the left side when non-zero. |
| `--rflank <N>` | — | `0` | all | Right-side flank in bp. Overrides `--flanks` for the right side when non-zero. |
| `--hap_split` | — | false | BAM, CRAM, PAF | Split output by haplotype tag into `<output>.h0.<ext>`, `<output>.h1.<ext>`, `<output>.h2.<ext>`. BAM/CRAM use the `HP` aux tag; PAF uses the `hp:i:` optional field. Records without the tag go to `h0`. Requires a file `--output` (not stdout). |
| `--use_paf_index` | — | `true` | PAF | Build and/or load a byte-offset index (`<paf>.idx`). Set to `false` to rebuild in memory every run without saving. |

## Flank precedence

When both `--flanks` and a per-side flag are set, the per-side value takes precedence:

```
effective_lflank = lflank  (if lflank > 0)  else  flanks
effective_rflank = rflank  (if rflank > 0)  else  flanks
```

## Output header formats

### BAM and CRAM

```
>read_name|chr:region_start-region_end|region_name
```

With `--hap_split` and a non-zero HP tag:

```
>read_name|chr:region_start-region_end|region_name|h1
```

### PAF

```
>contig_name|chr:ref_start-ref_end|bed_name|contig_name:query_start-query_end|strand
```

With `--hap_split` and a non-zero `hp:i:` tag:

```
>contig_name|chr:ref_start-ref_end|bed_name|contig_name:query_start-query_end|strand|h1
```

`region_name` / `bed_name` comes from the 4th column of the BED file, or `chr:start-end` if the column is absent.

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| non-zero | An error was encountered; details are written to stderr. |

## Diagnostics

bedpull writes progress and diagnostic messages to **stderr** so that stdout (or the output file) contains only sequence data. Messages include region names as they are processed, overlap counts, and warnings for partial alignments or unexpected HP tag values.
