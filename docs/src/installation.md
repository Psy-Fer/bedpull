# Installation

## From crates.io

```bash
cargo install bedpull
```

## Build from source

```bash
git clone https://github.com/Psy-Fer/bedpull
cd bedpull
cargo build --release
./target/release/bedpull --help
```

## Static binary for HPC (musl)

Cluster environments often run older glibc versions that are incompatible with binaries built on a modern workstation. Building with the `x86_64-unknown-linux-musl` target produces a fully static binary with no shared library dependencies.

```bash
rustup target add x86_64-unknown-linux-musl
RUSTFLAGS='-C link-arg=-s' cargo build --release --target x86_64-unknown-linux-musl
```

The binary is at `target/x86_64-unknown-linux-musl/release/bedpull`. Copy it to your HPC login node or scratch space; it requires no runtime dependencies.

## Requirements

- Rust 1.75 or later (install via [rustup](https://rustup.rs)).
- An indexed BAM (`.bam` + `.bam.bai`) for BAM mode.
- An indexed FASTA (`.fa` + `.fa.fai`) for PAF mode query extraction.
- A PAF file for PAF mode; the PAF must include the `cg:Z:` CIGAR tag (produced by `minimap2 -c` or `--cs=long`).
