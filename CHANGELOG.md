# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project uses [Semantic Versioning](https://semver.org/).

---

## [0.2.0] — 2026-06-05

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

### Changed

- **noodles upgraded from v0.66 → v0.110** — required three API fixes: `fasta::indexed_reader`
  → `fasta::io::indexed_reader`, `Query<File>` → `Query<R>` with bgzf trait bounds,
  `for result in query` → `for result in query.records()`.

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
