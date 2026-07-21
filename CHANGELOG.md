# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project uses [Semantic Versioning](https://semver.org/).

---

## [0.3.0] — unreleased

Coordinate-math correctness pass, driven by a differential test against
[bladerunner](https://github.com/Psy-Fer/bladerunner) and aimed at letting bladerunner link
`bedpull` as a library.

### Fixed

- **Boundary-coincidence bug in `get_read_cuts`** — a read whose alignment began (or ended)
  exactly at a region boundary was mis-sliced or silently dropped. The `start > 0` sentinel
  conflated "start not found yet" with "start found at read offset 0". Replaced with an explicit
  `found_start` flag plus an op-entry guard, and the op-level `== ref_end` shortcut no longer
  skips a mid-op `ref_start` crossing. This affected BAM/CRAM extraction (the PAF path already
  routed around it via `read_pos_at_ref`).

### Changed

- **`get_read_cuts` signature** is now
  `get_read_cuts(cigar_ops, align_start, align_end, region_start, region_end)`. It takes the
  alignment end and clamps its fire-on boundaries to the alignment span internally
  (`ref_start = max(region_start, align_start)`, `ref_end = min(region_end, align_end)`), so a
  partially-overlapping read yields a real slice instead of a `read_end == 0` sentinel.
  `region_start`/`region_end` are the desired (flank-expanded) window; the caller no longer
  pre-clamps. **`align_end` is exclusive (one past the last reference base)** — callers using
  noodles' inclusive `alignment_end()` must add 1.
- **`ReadCuts` gains `softclip_lead_start` and `softclip_trail_end`** (read offsets of the
  leading/trailing soft-clip runs; equal to `read_start`/`read_end` when no extension is
  available). Layout is now field-identical to bladerunner's `ReadCuts`.
- **noodles bumped `0.110` → `0.111`** to match bladerunner across the crate boundary.

### Added

- Regression battery for the boundary-coincidence cases and the soft-clip extension fields.

## [0.2.0] — 2026-07-16

### Added

- **Library crate** — core logic (`bed`, `cigar`, `paf`, `reads`, `utils`) is now a public
  library crate. Downstream tools can depend on `bedpull` directly without going through the CLI.
  `BamConfig` replaces `&Opts` in `get_bam_reads`; `read_bed` takes `&Path`.

- **`--flanks` / `--lflank` / `--rflank`** — expand the reference extraction window by N bases
  on each side before the CIGAR walk. Captures insertions that sit just outside BED coordinates
  (common in repeat regions). The window is clamped to the alignment boundaries.

- **`--partial`** — include reads that only partially overlap a BED region (default: spanning
  reads only). The extracted subsequence is clipped to whatever portion the read covers; the
  header records the actual covered coordinates.

- **`--min_mapq`** — skip reads with MAPQ below the threshold (default: 0, no filter). Reads
  with unavailable MAPQ (255) always pass.

- **`--fastq`** — write output as FASTQ instead of FASTA, preserving per-base quality scores
  (BAM input only).

- **`--min_region_quality`** — filter reads by the mean Phred quality of the extracted
  subsequence (BAM `--fastq` only; default: 0, no filter).

- **`--include_secondary` / `--include_supplementary`** — opt in to secondary (flag `0x100`)
  or supplementary (flag `0x800`) alignments; both are skipped by default.

- **`--hap_split`** — split output by haplotype HP aux tag into separate files
  `<output>.h0.<ext>`, `<output>.h1.<ext>`, `<output>.h2.<ext>`. Reads without an HP tag go
  to h0. HP values > 2 warn and route to h0 (BAM only).

- **Rustdoc** — `///` documentation on all public API items; zero warnings from `cargo doc`.

- **mdBook site** under `docs/` — introduction, installation, BAM mode, PAF mode, CLI
  reference, library usage, and concepts pages.

- **Test suite** — 54 unit tests across all modules plus 12 integration tests exercising the
  public library API, including a BAM-mode end-to-end test against
  `examples/rfc1_test.bam` (hs1 coords, chr4 RFC1 VNTR region).

- **CRAM support (`--cram`)** — `extract_from_cram` mirrors BAM extraction (mapq, secondary/
  supplementary, partial, flanks, hap_split, region quality all work identically); optional
  `--reference` builds a lazy `fasta::Repository` for reference-compressed CRAMs.

- **Automatic index creation** — missing `.bai`/`.crai`/`.fai` indexes are now built
  automatically (pure-Rust via `noodles`, no `samtools` dependency) instead of erroring.
  Falls back to an actionable `samtools index`/`faidx` error message if auto-build itself fails
  (e.g. a corrupt input file).

- **`--dedup`** — deduplicate output by read/contig name when the same one is emitted more than
  once (e.g. spanning multiple BED regions).

- **`--bed_out <file>`** (PAF mode only) — writes lifted-over query coordinates as BED6
  alongside the FASTA/FASTQ output.

- **`--unmapped <file>`** — write input BED regions that produced no output to a file, each
  preceded by a `#reason` comment (similar to liftOver's `-unmapped`), across all three modes.
  Distinguishes "no overlapping reads/alignments" from "candidates existed but were all filtered
  out" from "all matches already emitted for another region (--dedup)" from "region explicitly
  skipped (chromosome name contains '#')". `get_bam_reads`/`get_cram_reads` now return
  `(reads, candidates_seen)` instead of just `reads` to support this.

- **`--min_partial_coverage`** — minimum fraction (0.0-1.0) of the requested (region ± flanks)
  window a `--partial` read must cover to be included, similar to liftOver's `-minMatch`.
  Default `0.0` (no filter, any overlap accepted); requires `--partial`.

- **`--debug`** — gates verbose per-region/per-read diagnostic output (parsed args, raw BED
  record dumps, region banners, PAF query coordinates) behind an opt-in flag. Default output is
  now just the mode banner and a completion line.

- **`--stitch_records` / `--max_stitch_gap <N>`** (PAF mode) — when no single PAF record fully
  spans a region, look for a chain of records sharing a query contig and strand, contiguous in
  target space, whose combined span covers the region, and extract one sequence across the
  whole chain. Recovers large structural variants (typically big novel insertions) that cause
  the aligner to split one alignment into two or more records instead of a single record with a
  large indel operation. Off by default, PAF-mode only.

- **Partial-overlap header info** — `--partial` reads that don't fully span the requested
  (region ± `--flanks`) window now get `|missing_left=Nbp` / `|missing_right=Nbp` appended to
  their header, showing how many reference bases were not covered.

- **PAF mode zero-output notice** — PAF mode now reports when a region produces no output,
  matching BAM/CRAM's existing "No reads found... skipping" behavior.

- **Arg-combination and permission validation** — `--min_mapq`, `--include_secondary`,
  `--include_supplementary`, `--partial`, and `--stitch_records` now error clearly if combined
  with an incompatible mode (BAM/CRAM-only flags with `--paf`, and vice versa); output,
  `--bed_out`, and `--unmapped` directories are checked for actual writability, not just
  existence.

### Fixed

- **PAF minus-strand extraction** — query coordinates are now correctly flipped relative to
  `query_end` for `-` strand alignments, and the extracted sequence is reverse-complemented.
  Previously, minus-strand records produced the wrong subsequence.

- **Cryptic BAM header panic** ([#4](https://github.com/Psy-Fer/bedpull/issues/4)) — samtools
  ≥ 1.23.1 writes `@PG` records with an empty `VN:` field which noodles 0.110 rejects with
  `InvalidProgram(InvalidOther(Other("VN"), Missing))`. bedpull now catches `InvalidData`
  header errors, re-reads the raw BAM header bytes via bgzf, sanitises empty `VN:` fields to
  `VN:unknown`, and continues normally. A warning is printed so the user is aware.

- **Error handling** — replaced `.unwrap()` / `.expect()` / `process::exit` calls throughout
  with `anyhow`-based `?` propagation and `.context()` messages. Errors now surface as
  human-readable messages rather than Rust backtraces.

- **PAF off-by-one in sequence extraction** — `extract_from_fasta_coords` documented 0-based
  half-open coordinates but forwarded them to noodles' 1-based-inclusive region syntax without
  converting, so every PAF-mode extraction silently included one extra base before the intended
  start. Went undetected because no test previously exercised the FASTA-extraction step directly.

- **PAF CIGAR parsing crash on trailing newline** — `PafRecord::from_line` didn't trim the
  line ending before splitting fields, so a line whose last tab-field was `cg:Z:...` (common
  `minimap2` output without `--cs=long`) had a `\n` appended to the CIGAR string, crashing
  `to_cigar_ops()`. Only masked in prior tests because the bundled example PAF has `cg:Z:` as a
  middle field, not the last one.

- **BED `start == 0` crashed the entire run, all modes** — `read_bed` converted BED's 0-based
  coordinates straight into noodles `Position` (`NonZeroUsize`), and a region starting at the
  very first base of a chromosome has no valid `NonZeroUsize` representation, so it aborted the
  whole program, taking every other (valid) region in the same BED file down with it. Fixed by
  storing both bounds shifted by `+1` in `read_bed` and undoing the shift in one shared
  `region_bounds` helper used by all three extraction modes.

- **PAF mode `get_read_cuts` sentinel bug** — `get_paf_reads` used to call `get_read_cuts` per
  record, which has a documented sentinel convention where a boundary crossing gets
  misattributed whenever it coincides exactly with a record's own `align_start`. BAM/CRAM's
  `--partial` mode knows how to reinterpret this; PAF mode didn't, so any record whose clamped
  window edge landed exactly on its own start or end was silently discarded as invalid. Fixed by
  replacing the per-record `get_read_cuts` call with a new ambiguity-free `read_pos_at_ref`
  helper (`src/utils.rs`).

- **Flank-blind index query** — `index.query()` was called with the raw BED region instead of
  the flank-expanded window, so a `--stitch_records` chain member (or plain `--flanks` target)
  that only overlapped within the flank zone never reached extraction at all.

- **Hard crash on any single bad PAF record** — a malformed CIGAR or a contig name missing from
  the query FASTA aborted the entire run via `?` propagation, losing every other region's
  results. Downgraded to the same skip-with-warning pattern already used for "no valid overlap."

- **`read_start == read_end` treated as invalid** — a BED window falling entirely inside a
  deletion (no query bases correspond to it at all) was skipped as if invalid, even though a
  0bp extraction is the correct answer there. Narrowed the check to `read_start > read_end` only.

- **Stitch-chain selection nondeterminism** — `find_stitch_chain` grouped candidate chains in a
  `HashMap`, whose iteration order is randomized per-process, so on inputs where more than one
  chain could satisfy a region, which one won could vary run-to-run on identical input. Fixed by
  walking groups in a deterministic, first-seen-in-file order instead.

- **Reference-compressed CRAM without `--reference` decoded silently wrong** — whether a given
  CRAM needs an external reference isn't reliably detectable up front, so rather than guess,
  `--cram` without `--reference` now prints a warning that sequences will decode incorrectly
  (rather than error) if the CRAM turns out to be reference-compressed.

### Changed

- **noodles upgraded from v0.66 → v0.110** — required three API fixes: `fasta::indexed_reader`
  → `fasta::io::indexed_reader`, `Query<File>` → `Query<R>` with bgzf trait bounds,
  `for result in query` → `for result in query.records()`.

- **clap upgraded from v3 → v4** — removed the six clap-3-only `parse(from_os_str)` attributes
  in `src/cli.rs`; clap 4's derive picks a `PathBuf` value parser automatically. No other
  behavioural changes.

- **Consensus building removed.** The unused, `bio`-crate-based `_get_consensus` stub has been
  dropped along with the `bio` dependency. Use the external
  [poa-consensus](https://github.com/Psy-Fer/poa-consensus) tool instead — pair it with
  `--hap_split` to build a consensus per haplotype from bedpull's extracted reads (see README).

- **Code reorganisation** — `get_paf_reads`/`PafRead` moved to `reads.rs` alongside
  `get_bam_reads`/`get_cram_reads`; dead code removed; `extract_from_paf` simplified.

- **Reduced per-item I/O overhead** — BAM/CRAM readers are now opened once per run instead of
  once per BED region; the PAF file and query FASTA are opened once per run instead of once per
  alignment record. This was the dominant cost once a BED file has more than a handful of
  regions, independent of any future threading work.

---

## [0.1.1] — 2025

### Changed

- Updated crate categories for crates.io.

---

## [0.1.0] — 2025

Initial release.

- BAM mode: extract sequences from indexed BAM files using BED coordinates.
- PAF mode: build/load a byte-offset index of a PAF file and extract from a query FASTA.
- CIGAR-aware extraction via `get_read_cuts` — correctly handles insertions and deletions that
  coordinate-lifting tools miss.
