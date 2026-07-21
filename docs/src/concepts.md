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

If you convert a BED interval `[start, end)` to read coordinates by taking `align_start + (start - align_start)` and `align_start + (end - align_start)`, which is what a simple offset-based tool does, you are only counting reference positions. You skip over inserted bases because they are not counted in the reference arithmetic.

To include them, you must walk the CIGAR and count read positions and reference positions simultaneously. When you encounter an `I` operation, you increment only the read position counter. The read position counter therefore "runs ahead" of the reference position counter, and the difference between them at any point is the total number of inserted bases seen so far. When you eventually reach `region_end` in reference coordinates, the read position cursor is already past all insertions that occurred between `region_start` and `region_end`, and they are included in the slice.

This is exactly what `get_read_cuts` in `src/utils.rs` does.

## Coordinate frames

`get_read_cuts` takes `align_start`, `align_end`, `region_start`, and `region_end`. It is
**frame-agnostic**: all four must be in the *same* coordinate frame. The walk is relative —
`ref_pos` is initialised to `align_start` and advances one step per reference-consuming base — so
only the *difference* between the boundaries affects the returned read offsets. The returned
`ref_start`/`ref_end` come back in whatever frame you passed in. `align_end` is **exclusive**:
one position past the last reference base.

### How BAM/CRAM mode keeps the frame consistent

noodles represents BAM/CRAM alignment positions as 1-based `Position` values. `get_bam_reads` /
`get_cram_reads` pass `alignment_start()` as `align_start` and `alignment_end() + 1` as the
exclusive `align_end`. The region bounds come from the noodles `Region` built by `read_bed`, which
stores each BED coordinate shifted by `+1` — a workaround for `Position` being `NonZeroUsize`
(which cannot represent a BED `start` of `0`). That same `+1` shift puts the region into the same
1-based frame as the alignment positions, so the walk is internally consistent. Before returning,
both functions normalise `ref_start`/`ref_end` back to 0-based, so callers — and the
`|missing_left`/`|missing_right` header suffix — see plain BED coordinates. (Skipping that
normalisation is what caused the historical `missing_left=1bp`-on-every-read bug.)

Do not put `align_start` and the region bounds in different frames — mixing them shifts every
boundary by one.

### PAF coordinates

PAF-mode extraction does **not** go through `get_read_cuts`. Because a single BED window's two
edges can fall on two *different* chained PAF records, `get_paf_reads` resolves each edge
independently with `read_pos_at_ref`, a single-boundary reference→read mapper. PAF `target_start`
is 0-based and is used consistently with 0-based region bounds throughout that path. It is
validated against the RFC1 test case (`579 bp = 59 ref bp + 520 inserted bp`):
`src/main.rs::tests::get_read_cuts_rfc1_captures_insertion` pins the CIGAR-cut math and the
`get_paf_reads` tests in `src/reads.rs` pin the end-to-end extraction.

Note that the CIGAR-cut math and the final extracted FASTA sequence are two separate steps:
`get_paf_reads` converts the 0-based half-open cut coordinates into query-contig coordinates and
extracts via `extract_from_fasta_coords_reader`, which itself converts them into noodles'
1-based-inclusive region syntax internally. Both conversions matter: a bug in either one
silently shifts the extracted sequence by one base without changing the *reported* length.

## Soft clips

Soft-clipped bases (`S`) are present in the read sequence but are not aligned to the reference. Like insertions, they advance the read position without advancing the reference position. `get_read_cuts` handles them the same way as `I` operations, so soft-clipped bases at the start of a read shift the read-position cursor before the first aligned base.

`ReadCuts` also reports `softclip_lead_start` and `softclip_trail_end` — the read offsets of any leading/trailing soft-clip runs. A caller can use these to optionally extend an extraction into clipped "expansion" sequence that sits just outside the aligned window (bladerunner uses this for flanking-expansion detection). When there is no such extension they equal `read_start`/`read_end` respectively.

## Partial overlaps

`get_read_cuts` clamps its fire-on boundaries to the alignment's own span:
`ref_start = max(region_start, align_start)` and `ref_end = min(region_end, align_end)`. A read
that only partially overlaps the window therefore still produces a real slice — `ref_start` /
`ref_end` report the covered sub-span — instead of a sentinel. `read_end == 0` now means exactly
one thing: the alignment never reached the window at all.

The start boundary is tracked with an explicit `found_start` flag (not a `start > 0` sentinel) and
is recorded at op entry when `ref_pos == ref_start`. This is what makes `read_start == 0` a valid
result: a read whose alignment begins exactly at the window start correctly yields `read_start = 0`
instead of being dropped (the boundary-coincidence bug fixed in 0.3.0).

`resolve_cuts` in `src/reads.rs` (a private helper shared by `get_bam_reads` and `get_cram_reads`)
then applies the spanning policy:

- `read_end == 0` (no overlap): the read is skipped.
- Non-partial mode and the alignment does not fully cover the region (i.e.
  `align_start <= region_start && align_end >= region_end` is false): the read is skipped.
- Otherwise: the clamped `(read_start, read_end, ref_start, ref_end)` is returned unchanged.

The resulting `ref_start`/`ref_end` (normalised to 0-based) drive the `missing_left=Nbp` /
`missing_right=Nbp` header annotations for `--partial` reads that don't fully cover the requested
window.

The unit tests for `get_read_cuts`'s cut semantics are in `src/utils.rs`; the tests for
`resolve_cuts`'s spanning policy are in `src/reads.rs` (`resolve_cuts_*`).
