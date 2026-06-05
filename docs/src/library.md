# Library Usage

bedpull is published as a Rust library crate as well as a binary. You can use it programmatically to integrate CIGAR-aware extraction into your own pipelines.

## Adding to Cargo.toml

```toml
[dependencies]
bedpull = "0.1"
```

## Import paths

All public items are re-exported from the crate root:

```rust
use bedpull::{
    // CIGAR types
    CigarOp, CigarOps, ToCigarOps,
    // PAF index and record access
    PafIndex, PafIndexEntry, PafRecord, read_paf_record_at_offset,
    // BAM extraction
    BamConfig, BamRead, get_bam_reads,
    // Core utilities
    ReadCuts, get_read_cuts, calculate_qscore,
    read_bed, write_fasta_record, write_fastq_record,
    extract_from_fasta_coords,
};
```

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
walk — inserted bases are included, deleted bases are absent.

## Example 2: CIGAR coordinate math with get_read_cuts

This example shows how to use `get_read_cuts` directly to calculate the read slice
for a known CIGAR string. This is useful if you are integrating bedpull's coordinate
logic with data from another source.

```rust
use anyhow::Result;
use bedpull::{ToCigarOps, get_read_cuts};

fn show_cuts(cigar: &str, align_start: usize, region_start: usize, region_end: usize) -> Result<()> {
    // Parse a raw CIGAR string from a PAF cg:Z: tag or any str source.
    let ops = cigar.to_cigar_ops()?;

    // align_start is 1-based (as from a BAM record or PAF target_start field).
    // region_start / region_end are 0-based (as from a BED file).
    let cuts = get_read_cuts(&ops, align_start, region_start, region_end);

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
    // CIGAR: 3M5I4M
    // Region: ref bases 3 to 7 (0-based BED coordinates)
    //
    // Expected: the insertion is captured, so the slice is 9 bases
    // (1 match + 5 inserted + 3 match), not 4 reference bases.
    show_cuts("3M5I4M", 1, 3, 7)?;
    Ok(())
}
```

Output:
```
read slice: [2..11]  (9 bases)
reference covered: [3-7]
```

## Example 3: PAF-based extraction

```rust
use anyhow::Result;
use bedpull::{PafIndex, read_paf_record_at_offset, ToCigarOps, get_read_cuts, extract_from_fasta_coords, write_fasta_record};
use std::io;

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

    let entries = index.query(chrom, region_start, region_end);
    let stdout = io::stdout();
    let mut out = stdout.lock();

    for entry in entries {
        let record = read_paf_record_at_offset(paf_path, entry.offset)?;
        if let Some(cigar_str) = &record.cigar {
            let ops = cigar_str.as_str().to_cigar_ops()?;
            let cuts = get_read_cuts(&ops, record.target_start, region_start, region_end);

            let query_start = record.query_start + cuts.read_start;
            let query_end = record.query_start + cuts.read_end;

            let sequence = extract_from_fasta_coords(
                query_fasta,
                &record.query_name,
                query_start,
                query_end,
            )?;

            let header = format!(
                "{}|{}:{}-{}",
                record.query_name, chrom, region_start, region_end
            );
            write_fasta_record(&mut out, &header, &sequence)?;
        }
    }

    Ok(())
}
```

## Notes on the coordinate convention

`get_read_cuts` uses a **mixed coordinate system** that matches the conventions of its two callers:

- `align_start` is **1-based** (noodles `Position` from BAM; PAF `target_start` is 0-based but is used with the CIGAR walk in a way that matches 1-based BAM behaviour — see the [Concepts](concepts.md) chapter for details).
- `region_start` and `region_end` are **0-based** (directly from BED).

Do not normalise both to the same base before calling `get_read_cuts`. The test suite in `src/utils.rs` documents the exact expected behaviour for all edge cases.
