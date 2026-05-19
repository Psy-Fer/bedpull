# PAF Mode

PAF mode extracts sequences from a query FASTA using coordinates derived from a PAF (Pairwise mApping Format) alignment file. The typical use case is assembly-to-reference mapping: you have aligned an assembly to a reference with `minimap2`, you have a set of reference BED regions, and you want to pull out the corresponding assembly sequence — including any insertions the assembly carries within those regions.

## Requirements

- A PAF file with the `cg:Z:` CIGAR tag. Produce this with `minimap2 -c` or `minimap2 --cs=long`.
- The query FASTA (the assembly) indexed with `samtools faidx`.

## Basic usage

```bash
bedpull \
    --paf assembly_to_ref.paf \
    --query_ref assembly.fasta \
    --bed regions.bed \
    --output extracted.fasta
```

Output is always FASTA in PAF mode (quality scores are not available from a FASTA).

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

## Output header format

Each FASTA record produced in PAF mode has a header of the form:

```
>query_name|ref_region_name:ref_start-ref_end|query_query_name:query_start-query_end
```

The query coordinates are computed as `paf_record.query_start + cuts.read_start` and `paf_record.query_start + cuts.read_end`, giving the precise position within the assembly sequence.

## Alignment coverage warnings

If an overlapping alignment does not fully span the requested BED region, bedpull emits a warning but still extracts whatever portion of the region the alignment covers:

```
Warning: Alignment starts after region start, may be incomplete
Warning: Alignment ends before region end, may be incomplete
```

## Example: extracting RFC1 from a paternal assembly

```bash
minimap2 -a -c chm13v2.0.fa HG002_paternal.fasta > hg002pat_to_chm13.paf

bedpull \
    --paf hg002pat_to_chm13.paf \
    --query_ref HG002_paternal.fasta \
    --bed rfc1_vntr.bed \
    --output rfc1_paternal.fasta
```

The first run builds `hg002pat_to_chm13.paf.idx`. The second run loads it in milliseconds. The output FASTA contains the paternal assembly sequence corresponding to the RFC1 region, including any inserted bases that are not present in the CHM13 reference.
