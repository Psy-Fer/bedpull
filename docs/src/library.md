# Library Usage

bedpull is published as a Rust library crate as well as a binary. You can use it programmatically to integrate CIGAR-aware extraction into your own pipelines.

## Adding to Cargo.toml

```toml
[dependencies]
bedpull = "0.3"
```

## Import paths

All public items are re-exported from the crate root:

```rust
use bedpull::{
    // CIGAR types
    CigarOp, CigarOps, ToCigarOps,
    // PAF index and record access
    PafIndex, PafIndexEntry, PafRecord, read_paf_record_at_offset, read_paf_record_from_reader,
    // BAM/CRAM/PAF extraction
    BamConfig, BamRead, PafRead, StitchConfig, get_bam_reads, get_cram_reads, get_paf_reads,
    // Core utilities
    ReadCuts, get_read_cuts, calculate_qscore, revcomp,
    read_bed, write_fasta_record, write_fastq_record,
    extract_from_fasta_coords, extract_from_fasta_coords_reader,
};
```

`read_paf_record_at_offset` and `extract_from_fasta_coords` each open their file fresh on
every call, which is fine for a one-off lookup. `read_paf_record_from_reader` and
`extract_from_fasta_coords_reader` take an already-open reader instead, so you can open the
PAF/FASTA once and reuse it across many records; see [Example 3](#example-3-paf-based-extraction).

## Example 1: BAM extraction with BamConfig

```rust
use anyhow::Result;
use bedpull::{BamConfig, get_bam_reads};
use noodles::bam;
use noodles::core::Region;

fn extract_region(bam_path: &str, region_str: &str) -> Result<()> {
    let config = BamConfig {
        min_mapq: 20,
        include_secondary: false,
        include_supplementary: false,
        partial: false,
        min_region_quality: 0.0,
    };

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region: Region = region_str.parse()?;
    let query = reader.query(&header, &region)?;

    // No flanks; extract exactly the BED region.
    let reads = get_bam_reads(&config, query, &region, 0, 0)?;

    for (name, seq, _qual, ref_start, ref_end, hap) in reads {
        let seq_str = String::from_utf8(seq)?;
        println!(
            ">{} ref={}:{}-{} hp={}",
            name, region_str, ref_start, ref_end, hap
        );
        println!("{}", seq_str);
    }

    Ok(())
}
```

`get_bam_reads` returns a `Vec<BamRead>`, where each element is a tuple
`(name, sequence, quality_string, ref_start, ref_end, haplotype)`. The `sequence`
field contains only the bases that fall within the requested region after the CIGAR
walk: inserted bases are included, deleted bases are absent.

## Example 2: CIGAR coordinate math with get_read_cuts

This example shows how to use `get_read_cuts` directly to calculate the read slice
for a known CIGAR string. This is useful if you are integrating bedpull's coordinate
logic with data from another source.

```rust
use anyhow::Result;
use bedpull::{ToCigarOps, get_read_cuts};

fn show_cuts(
    cigar: &str,
    align_start: usize,
    align_end: usize,
    region_start: usize,
    region_end: usize,
) -> Result<()> {
    // Parse a raw CIGAR string from a PAF cg:Z: tag or any str source.
    let ops = cigar.to_cigar_ops()?;

    // align_start / align_end and region_start / region_end must all be in the same
    // coordinate frame (the walk is relative). align_end is exclusive — one past the
    // last reference base. region_start / region_end are the desired window (expand
    // them yourself if you want flanks); get_read_cuts clamps its fire-on boundaries
    // to the alignment span internally.
    let cuts = get_read_cuts(&ops, align_start, align_end, region_start, region_end);

    println!(
        "read slice: [{}..{}]  ({} bases)",
        cuts.read_start,
        cuts.read_end,
        cuts.read_end - cuts.read_start,
    );
    println!(
        "reference covered: [{}-{}]",
        cuts.ref_start, cuts.ref_end,
    );

    Ok(())
}

fn main() -> Result<()> {
    // A read aligned at ref position 1, with a 5-base insertion inside the region.
    // CIGAR: 3M5I4M consumes 7 reference bases, so align_end = 1 + 7 = 8.
    // Region: ref bases 3 to 7.
    //
    // Expected: the insertion is captured, so the slice is 9 bases
    // (1 match + 5 inserted + 3 match), not 4 reference bases.
    show_cuts("3M5I4M", 1, 8, 3, 7)?;
    Ok(())
}
```

Output:
```
read slice: [2..11]  (9 bases)
reference covered: [3-7]
```

## Example 3: PAF-based extraction

`get_paf_reads` is the high-level entry point: it drives the CIGAR walk, the FASTA extraction,
and the minus-strand reverse-complement for you. It takes an already-open PAF reader and an
already-open indexed FASTA reader so you can call it once per region (or even across many
regions) without reopening either file each time. Reopening per record is the single biggest
cost in PAF-mode extraction once you're processing more than a handful of alignments.

```rust
use anyhow::Result;
use bedpull::{PafIndex, StitchConfig, get_paf_reads, write_fasta_record};
use noodles::fasta;
use std::fs::File;
use std::io::{self, BufReader};

fn extract_from_paf(
    paf_path: &str,
    query_fasta: &str,
    chrom: &str,
    region_start: usize,
    region_end: usize,
) -> Result<()> {
    // Build or load the byte-offset index.
    let index_path = format!("{}.idx", paf_path);
    let index = if std::path::Path::new(&index_path).exists() {
        PafIndex::load(&index_path)?
    } else {
        let idx = PafIndex::build(paf_path)?;
        idx.save(&index_path)?;
        idx
    };

    // Open both files once; reuse them across every region you process.
    let mut paf_reader = BufReader::new(File::open(paf_path)?);
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(query_fasta)?;

    let entries = index.query(chrom, region_start, region_end);
    let stdout = io::stdout();
    let mut out = stdout.lock();

    let reads = get_paf_reads(
        &mut paf_reader,
        &mut fasta_reader,
        &entries,
        region_start,
        region_end,
        0, // lflank
        0, // rflank
        StitchConfig::default(), // no cross-record stitching
        false, // debug
    )?;

    for (sequence, query_name, query_start, query_end, strand, _hap) in reads {
        let header = format!(
            "{}|{}:{}-{}|{}:{}-{}|{}",
            query_name, chrom, region_start, region_end, query_name, query_start, query_end, strand
        );
        write_fasta_record(&mut out, &header, &sequence)?;
    }

    Ok(())
}
```

To process many regions, open `paf_reader` and `fasta_reader` once outside a loop over regions
and call `get_paf_reads` once per region, exactly what `bedpull`'s own CLI does internally.

## Notes on the coordinate convention

`get_read_cuts` is **frame-agnostic**: `align_start`, `align_end`, `region_start`, and
`region_end` must all be in the *same* coordinate frame, because the walk is relative
(`ref_pos` starts at `align_start` and advances). Only the *difference* between the boundaries
determines the returned read offsets; the returned `ref_start`/`ref_end` inherit whichever frame
you passed in. `align_end` is **exclusive** (one past the last reference base) — with noodles'
1-based inclusive `alignment_end()`, add 1.

In BAM/CRAM mode, `get_bam_reads`/`get_cram_reads` supply all four in noodles' 1-based frame and
then normalise the returned `ref_start`/`ref_end` back to 0-based (BED) before handing them to you.
See the [Concepts](concepts.md) chapter for the full derivation. The test suite in `src/utils.rs`
documents the exact expected behaviour for all edge cases.

`extract_from_fasta_coords`/`extract_from_fasta_coords_reader` take `start`/`end` as **0-based
half-open**, the same convention as everywhere else in the crate, and convert internally to
noodles' 1-based-inclusive region syntax. Pass `read_cuts`/`get_paf_reads` output straight
through without adjusting it yourself.
