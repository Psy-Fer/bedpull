# PAF Mode

PAF mode extracts sequences from a query FASTA using coordinates derived from a PAF (Pairwise mApping Format) alignment file. The typical use case is assembly-to-reference mapping: you have aligned an assembly to a reference with `minimap2`, you have a set of reference BED regions, and you want to pull out the corresponding assembly sequence — including any insertions the assembly carries within those regions.

## Requirements

- A PAF file with the `cg:Z:` CIGAR tag. Produce this with `minimap2 -c` or `minimap2 --cs=long`.
- The query FASTA (the assembly). If its `.fai` index is missing, bedpull builds it
  automatically the first time you run against that file — no `samtools` required.

`--partial` (see [BAM mode](./bam-mode.md#partial-overlaps)) is BAM/CRAM-only and errors if combined with `--paf`; PAF mode always extracts whatever portion of the region an overlapping alignment covers (see [Alignment coverage warnings](#alignment-coverage-warnings) below).

## Basic usage

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output extracted.fasta
```

Output is always FASTA in PAF mode (quality scores are not available from a FASTA). Omit `--output` (or pass `-`) to write to stdout.

## The PAF index

Scanning a multi-gigabyte PAF file for every BED region would be slow. bedpull builds a byte-offset index the first time you run on a PAF file and saves it alongside the PAF as `<paf_file>.idx`. On subsequent runs the index is loaded instead of rebuilt.

### Index format

The index is a plain TSV with four columns:

```
<chrom>  <byte_offset>  <target_start>  <target_end>
```

Each row records where one PAF record lives on disk (`byte_offset` is the byte position of the start of that line) plus its target chromosome coordinates so that overlap queries can be answered without seeking to the record.

### Index lifecycle

| Condition | Action |
|-----------|--------|
| `<paf>.idx` does not exist | Index is built and saved automatically |
| `<paf>.idx` exists | Index is loaded from disk |
| `--use_paf_index false` | Index is built in memory every run; never saved |

If you modify the PAF file, delete the `.idx` file to force a rebuild.

### Disabling the index

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output extracted.fasta \
    --use_paf_index false
```

## Flanks

As with BAM mode, `--flanks`, `--lflank`, and `--rflank` expand the extraction window in reference coordinates before the CIGAR walk:

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output extracted.fasta \
    --flanks 500
```

The effective window is clamped to the actual extent of each overlapping alignment, so a flank that extends beyond the alignment boundary does not cause an error — it is silently truncated.

## BED liftover output

PAF mode can simultaneously write a BED6 file containing the lifted-over query coordinates for every extracted region. This is useful for the two-step liftover workflow: extract assembly coordinates from the PAF, then use those coordinates to query reads aligned to the assembly.

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed ref_regions.bed \
    --output extracted.fasta \
    --bed_out assembly_regions.bed
```

The output BED6 has the form:

```
query_name  query_start  query_end  bed_name  0  strand
```

`--bed_out` can be combined with `--output` independently. Pass `-` to write the BED to stdout (but not at the same time as `--output -`).

## Haplotype splitting

If PAF records carry the `hp:i:` optional tag (set by assemblers or phasing tools that emit phased contigs), you can split the output by haplotype:

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output out.fasta \
    --hap_split
```

This creates three output files:

| File | Contents |
|------|----------|
| `out.h0.fasta` | Alignments with no `hp:i:` tag |
| `out.h1.fasta` | Alignments with `hp:i:1` |
| `out.h2.fasta` | Alignments with `hp:i:2` |

## Deduplication

A query contig that aligns to two BED regions will appear in the output once per region by default. To emit each contig only once across all regions:

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output out.fasta \
    --dedup
```

## Output header format

Each FASTA record produced in PAF mode has a header of the form:

```
>contig_name|chr:ref_start-ref_end|bed_name|contig_name:query_start-query_end|strand
```

When the alignment carries an `hp:i:` tag, a haplotype suffix is appended:

```
>contig_name|chr:ref_start-ref_end|bed_name|contig_name:query_start-query_end|strand|h1
```

The query coordinates (`query_start`, `query_end`) are the 0-based half-open positions within the query contig that were extracted. These are the same values written to `--bed_out`.

## Alignment coverage warnings

If an overlapping alignment does not fully span the requested BED region, bedpull emits a warning but still extracts whatever portion of the region the alignment covers:

```
Warning: Alignment starts after region start, may be incomplete
Warning: Alignment ends before region end, may be incomplete
```

These warnings are always printed. A per-alignment `Query coords: contig:start-end (strand X)`
line is also available but only shown with `--debug` (see below), since it's noisy for PAF
files with many overlapping alignments.

## Cross-record stitching

A single PAF record only ever covers a contiguous block of one alignment. A large structural
variant — most often a big novel insertion with no homologous target sequence — frequently
causes the aligner to emit *two or more separate, chained records* instead of one record with a
large indel operation. Without help, a BED region spanning that boundary is only ever partially
extracted from whichever single record happens to overlap it (or produces two separate partial
hits, one per record).

`--stitch_records` looks for a chain of records that share a query contig and strand, are
contiguous in target space (within `--max_stitch_gap` bp), and together span the requested
region — then extracts one sequence across the whole chain in a single slice, rather than
walking each record's CIGAR separately. Any "gap" in query space between chain members (the
inserted sequence itself, which has no aligned counterpart in either record) is spliced in
directly from the raw query FASTA, since that gap *is* the structural variant the aligner split
around.

```bash
bedpull --paf assembly.paf --query_ref assembly.fa --bed regions.bed \
    --stitch_records --max_stitch_gap 10000 --output output.fasta
```

`--max_stitch_gap` (default `10000`) bounds only the *target*-side gap between consecutive
chain members — the reconstructed query-side splice (the insertion itself) can be much larger
and isn't separately bounded. Off by default; when a single record already fully covers a
region, stitching is skipped even if enabled, since there's nothing to bridge.

## Diagnostics

By default bedpull only prints index-build/load status, the mode banner, and a completion
line. Pass `--debug` to see verbose per-region diagnostics: the parsed CLI args, each BED record
as it's read, a banner for every region as it's processed, the overlapping-alignment count per
region, and per-alignment query coordinates. With `--stitch_records`, a successful stitch also
prints `Stitched N chained record(s): contig:start-end (strand X)`.

## Example: extracting RFC1 from a paternal assembly

```bash
minimap2 -a -c chm13v2.0.fa HG002_paternal.fasta > hg002pat_to_chm13.paf

bedpull \
    --paf hg002pat_to_chm13.paf \
    --query_ref HG002_paternal.fasta \
    --bed rfc1_vntr.bed \
    --output rfc1_paternal.fasta \
    --bed_out rfc1_paternal_coords.bed
```

The first run builds `hg002pat_to_chm13.paf.idx`. The second run loads it in milliseconds. `rfc1_paternal.fasta` contains the paternal assembly sequence over the RFC1 region, including inserted bases invisible to reference-only tools. `rfc1_paternal_coords.bed` contains the corresponding query coordinates, ready to use as input to a subsequent BAM-mode run against reads aligned to the paternal assembly.
