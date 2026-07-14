# CLI Reference

## Synopsis

```
bedpull [OPTIONS] --bed <BED>
```

At least one input mode must be specified: `--bam`, `--cram`, or `--paf` + `--query_ref`. `--bam` and `--cram` are mutually exclusive. `--paf` cannot be combined with `--bam` or `--cram`.

## Input flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--bam <FILE>` | `-b` | — | BAM | Coordinate-sorted BAM file. A missing `.bam.bai`/`.bai` index is built automatically. |
| `--cram <FILE>` | — | — | CRAM | CRAM file. A missing `.cram.crai` index is built automatically. |
| `--reference <FILE>` | `-f` | — | CRAM | Reference FASTA used during CRAM compression. Required for reference-compressed CRAMs; not needed for CRAMs with embedded sequences. A missing `.fai` index is built automatically. |
| `--paf <FILE>` | — | — | PAF | PAF alignment file. Must include the `cg:Z:` CIGAR tag (`minimap2 -c`). Must be used together with `--query_ref`. |
| `--query_ref <FILE>` | — | — | PAF | Query FASTA to extract sequence from. A missing `.fai` index is built automatically. |
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
| `--min_mapq <N>` | — | `0` | BAM, CRAM | Minimum mapping quality. Reads with MAPQ below `N` are skipped. `0` disables the filter. Reads with unavailable MAPQ (255) always pass. Errors if combined with `--paf`. |
| `--min_region_quality <F>` | — | `0` | BAM, CRAM | Minimum mean Phred quality of the extracted subsequence. Reads below this threshold are discarded. `0` disables the filter. Computed in error-probability space. Applies regardless of `--fastq`. |
| `--include_secondary` | — | false | BAM, CRAM | Include secondary alignments (SAM flag `0x100`). Excluded by default. Errors if combined with `--paf`. |
| `--include_supplementary` | — | false | BAM, CRAM | Include supplementary alignments (SAM flag `0x800`). Excluded by default. Errors if combined with `--paf`. |
| `--partial` | — | false | BAM, CRAM | Include reads that only partially overlap the BED region. By default only reads spanning the entire region are returned. Errors if combined with `--paf` — PAF mode always returns whatever portion of the region an overlapping alignment covers. |
| `--dedup` | — | false | all | Deduplicate output: if the same read or contig name is seen more than once across BED regions, emit it only the first time. |

## Extraction flags

| Flag | Short | Default | Mode | Description |
|------|-------|---------|------|-------------|
| `--flanks <N>` | — | `0` | all | Expand the extraction window by `N` bp on both sides of each BED region before the CIGAR walk. |
| `--lflank <N>` | — | `0` | all | Left-side flank in bp. Overrides `--flanks` for the left side when non-zero. |
| `--rflank <N>` | — | `0` | all | Right-side flank in bp. Overrides `--flanks` for the right side when non-zero. |
| `--hap_split` | — | false | BAM, CRAM, PAF | Split output by haplotype tag into `<output>.h0.<ext>`, `<output>.h1.<ext>`, `<output>.h2.<ext>`. BAM/CRAM use the `HP` aux tag; PAF uses the `hp:i:` optional field. Records without the tag go to `h0`. Requires a file `--output` (not stdout). |
| `--use_paf_index` | — | `true` | PAF | Build and/or load a byte-offset index (`<paf>.idx`). Set to `false` to rebuild in memory every run without saving. |
| `--debug` | — | false | all | Print verbose per-region/per-read diagnostic output: parsed CLI args, raw BED records, region banners, PAF overlap counts and per-alignment query coordinates. See [Diagnostics](#diagnostics). |

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

With `--partial` and a read that doesn't fully span the requested (region ± flank) window:

```
>read_name|chr:region_start-region_end|region_name|missing_left=12bp|missing_right=8bp
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

bedpull writes progress and diagnostic messages to **stderr** so that stdout (or the output file) contains only sequence data.

By default this is minimal: a mode banner, index-build/load status (PAF mode), and a completion
line — plus warnings that indicate a real data issue (an unexpected HP tag value, a PAF
alignment that doesn't fully cover a region, a region skipped because it contains `#`, or a
region with no overlapping reads/alignments). These warnings are always shown, regardless of
`--debug`.

Pass `--debug` to additionally see: the fully parsed `Opts` struct, each BED record as it's
parsed, a banner for every region as it's processed, and (PAF mode) the overlapping-alignment
count per region plus per-alignment query coordinates.

## Input validation

Before extraction starts, bedpull validates the input files and argument combinations and fails
fast with an actionable message rather than partway through processing. This includes: file
existence, output/`--bed_out` directory existence *and* writability, mutually exclusive modes
(`--bam`+`--cram`, alignment modes+`--paf`), `--paf`/`--query_ref` must be used together,
`--fastq`/`--min_region_quality`/`--min_mapq`/`--include_secondary`/`--include_supplementary`/`--partial`
requiring `--bam`/`--cram`, and `--hap_split` requiring a file `--output` (not stdout).
