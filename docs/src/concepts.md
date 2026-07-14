# Concepts

This chapter explains the technical foundations that underpin bedpull's extraction logic. It is aimed at readers who want to understand why the tool behaves the way it does, or who are integrating the library into their own code.

## CIGAR strings

A CIGAR string encodes how an aligned sequence maps to a reference. It is a series of (length, operation) pairs. The operations relevant to coordinate extraction are:

| Operation | Symbol | Advances ref? | Advances read? |
|-----------|--------|---------------|----------------|
| Match/mismatch | `M` | yes | yes |
| Sequence match | `=` | yes | yes |
| Sequence mismatch | `X` | yes | yes |
| Insertion | `I` | **no** | yes |
| Deletion | `D` | yes | **no** |
| Skip (intron) | `N` | yes | **no** |
| Soft clip | `S` | **no** | yes |
| Hard clip | `H` | **no** | **no** |
| Padding | `P` | **no** | **no** |

The key asymmetry: **insertions advance the read position without advancing the reference position**. Inserted bases exist in the read (or assembly) but have no corresponding reference coordinate. This is why reference-coordinate slicing misses them.

## Why insertions have no reference coordinate

Reference coordinates number the bases of the reference sequence. An insertion is a sequence of bases present in the read that does not correspond to any reference base. In the CIGAR string, an `I` operation of length `n` means: "the next `n` bases in the read have no reference counterpart."

If you convert a BED interval `[start, end)` to read coordinates by taking `align_start + (start - align_start)` and `align_start + (end - align_start)` — which is what a simple offset-based tool does — you are only counting reference positions. You skip over inserted bases because they are not counted in the reference arithmetic.

To include them, you must walk the CIGAR and count read positions and reference positions simultaneously. When you encounter an `I` operation, you increment only the read position counter. The read position counter therefore "runs ahead" of the reference position counter, and the difference between them at any point is the total number of inserted bases seen so far. When you eventually reach `region_end` in reference coordinates, the read position cursor is already past all insertions that occurred between `region_start` and `region_end`, and they are included in the slice.

This is exactly what `get_read_cuts` in `src/utils.rs` does.

## The mixed 0-based / 1-based coordinate convention

`get_read_cuts` takes three coordinate arguments with different bases:

- `align_start` — **1-based**
- `region_start` — **0-based** (from BED)
- `region_end` — **0-based** (from BED)

This is intentional and non-obvious. Here is why.

### BAM coordinates are 1-based in noodles

The noodles library represents BAM alignment positions as `Position`, which is a 1-based type. When bedpull reads `record.alignment_start()` and converts it with `usize::from(...)`, the result is a 1-based integer. This value is passed directly to `get_read_cuts` as `align_start`.

### BED coordinates are 0-based

The BED format is 0-based, half-open: a BED record with `start=100, end=200` refers to reference bases 100–199 (0-based). bedpull reads BED coordinates directly from the file and passes them to `get_read_cuts` as `region_start` and `region_end`.

### The offset cancels correctly

Inside `get_read_cuts`, the walk begins with `ref_pos = align_start` (1-based). The trigger condition is `ref_pos == region_start` or `ref_pos == region_end` (0-based). These conditions fire when the 1-based running position equals the 0-based BED coordinate. Because BED `start` is 0-based and the alignment start is 1-based, the trigger fires exactly one position earlier in 0-based space than it would if both were 1-based — which is what you want: BED `start=100` means "the first base of the region is at 0-based position 100", and the 1-based ref_pos of `100` is that same base.

In short: the mixed convention is not a bug. It is the natural consequence of BED being 0-based and noodles BAM positions being 1-based, and the arithmetic works out correctly because the two offsets cancel. **Do not normalise `align_start` to 0-based before calling `get_read_cuts`.**

### PAF coordinates

PAF `target_start` is 0-based. However, bedpull passes it to `get_read_cuts` as `align_start` (the 1-based slot). In practice this means the PAF-mode walk starts at one position "before" the true alignment start relative to what the BAM walk does. The PAF integration (`get_paf_reads` in `src/reads.rs`) was developed and validated against the RFC1 test case, and the extracted lengths are correct (`579 bp = 59 ref bp + 520 inserted bp`). The test in `src/main.rs::tests::get_read_cuts_rfc1_captures_insertion` pins this behaviour at the CIGAR-cut level. Any change to the coordinate handling must keep that test passing.

Note that the CIGAR-cut math and the final extracted FASTA sequence are two separate steps:
`get_paf_reads` converts the 0-based half-open cut coordinates into query-contig coordinates and
extracts via `extract_from_fasta_coords_reader`, which itself converts them into noodles'
1-based-inclusive region syntax internally. Both conversions matter — a bug in either one
silently shifts the extracted sequence by one base without changing the *reported* length.

## Soft clips

Soft-clipped bases (`S`) are present in the read sequence but are not aligned to the reference. They advance the read position without advancing the reference position — similar to insertions. `get_read_cuts` handles them the same way as `I` operations, so soft-clipped bases at the start of a read shift the read-position cursor before the first aligned base.

## Partial overlaps

When a read's alignment starts after `region_start`, the `region_start` trigger in `get_read_cuts` never fires (because `ref_pos` is initialised to `align_start`, which is already past `region_start`). If `region_end` is subsequently reached, the position is stored in `read_start` (because `start == 0` at that point, indicating the start sentinel has not been set). `read_end` remains `0`.

`resolve_cuts` in `src/reads.rs` — a private helper shared by `get_bam_reads` and
`get_cram_reads` — detects this sentinel pattern when `BamConfig::partial` is `true` and remaps
the fields into the actual `(read_start, read_end, ref_start, ref_end)` covered by the returned
slice:

- If `align_start > region_start` (left-partial or contained): the read is extracted from position `0` to `read_cuts.read_start` (where the `region_end` was reached).
- If `read_cuts.read_end == 0` and partial mode is on (right-partial): the read is extracted from `read_cuts.read_start` to the end of the read.
- If `read_cuts.read_end == 0` and partial mode is off: the read is skipped.

The corrected `ref_start`/`ref_end` this produces are what drive the `missing_left=Nbp` /
`missing_right=Nbp` header annotations for `--partial` reads that don't fully cover the
requested window.

The unit tests for `get_read_cuts`'s underlying cut semantics are in `src/utils.rs` under the
`// --- partial overlap ---` section; the tests for `resolve_cuts`'s field-remapping on top of
that are in `src/reads.rs` (`resolve_cuts_*`).
