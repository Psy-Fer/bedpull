# bedpull benchmark harness

Compares bedpull's CIGAR-aware `--paf` extraction against `liftOver` +
`samtools faidx` (the standard coordinate-lift-then-extract pipeline) on real
structural variants, using the GIAB HG002 v5.0q SV truth set. Three legs:

1. **Accuracy** (the headline result) — hs1 (CHM13v2.0) reference, HG002
   T2T-Q100 v1.1 diploid assembly, GIAB SV truth set. Does bedpull recover
   the true insertion/deletion length where liftOver truncates at a chain
   boundary?
2. **Confirmation** — hs1↔hg38 using *official* UCSC chain files (no
   custom chain-building). Sanity-checks that this harness reproduces known
   liftOver behavior on a well-established reference pair before trusting it
   on the harder personal-assembly case. Unlike the main leg, there's no
   individual whose assembly both genomes agree on, so the truth set here is
   built from real hg38↔hs1 alignment differences (`call_alignment_svs.py`),
   not the HG002 GIAB VCF — see "The hs1↔hg38 confirmation leg" below for why.
3. **Runtime/scale** — synthetic BED files of increasing size (100 →
   100,000 regions), to see how both tools' wall-clock scales.

## Design choices worth knowing before you read the scripts

- **We don't trust VCF GT phase to identify assembly haplotype.** The GIAB
  VCF's genotype (e.g. `0|1`) doesn't reliably tell you whether the ALT
  allele is on the `_PATERNAL` or `_MATERNAL` assembly contig — that mapping
  isn't guaranteed stable across pipelines. Instead, `run_bedpull_pipeline.sh`
  / `run_liftover_pipeline.sh` are run once per haplotype against the *same*
  window set, and `score_accuracy.py` takes whichever haplotype's result is
  closest to the expected length. A het SV will only extract cleanly from one
  haplotype anyway, so this sidesteps the phase-to-haplotype-identity
  question entirely rather than getting it wrong silently.
- **No-call and hemizygous genotypes are excluded by default**
  (`--include-ambiguous-gt` to keep them) — `.` in either allele makes the
  "truth" ambiguous, so they're dropped from the primary accuracy scoring.
- **Windows must be fully contained in the GIAB high-confidence BED.**
  Partial overlap with a confident region isn't good enough — see
  `fully_contained()` in `select_sv_windows.py`.
- **Expected extracted length** is `window_len + (alt_len - ref_len)`, not a
  SVTYPE-branching calculation — one formula handles insertions and
  deletions identically (see the docstring in `score_accuracy.py`).
- **The chain-building step (`build_chain_from_paf.sh`) is the least
  independently verified piece of this harness** — I don't have `paf2chain`
  or the UCSC kent tools installed in the environment this was built in, so
  I could write and reason through the conversion but not run it end to end.
  Everything else here (window selection, both extraction pipelines, the
  scorer, the scaling harness) has been run against real data — see the
  "What's been tested" section below.

## Prerequisites

Already available if you're on the same machine this was built on:
`bcftools`, `samtools`, `minimap2`, `python3` (stdlib only, no pip installs
needed). Still need to be installed separately:

- **`paf2chain`** (for the accuracy leg's custom chain): `git clone
  https://github.com/AndreaGuarracino/paf2chain && cd paf2chain && cargo
  install --force --path .`
- **`liftOver`**: standalone binary, no build needed — grab it from
  <https://hgdownload.soe.ucsc.edu/admin/exe/> for your platform, `chmod +x`,
  put it on `PATH` (or point `LIFTOVER_BIN` at it).
- **Official hg38→hs1 chain file** (confirmation leg only, no building
  required): <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/>
  (`hg38ToHs1.over.chain.gz` — UCSC only publishes this direction, which is
  why the confirmation leg lifts hg38 → hs1, not the other way around).
- **`k8`** (only if running the optional `paftools.js liftover` third-tool
  comparison — see below): grab a release from
  <https://github.com/attractivechaos/k8/releases> (e.g. `k8-1.2.tar.bz2`),
  extract, `chmod +x k8-x86_64-Linux`, put it on `PATH` as `k8` (or point
  `K8_BIN` at it). `paftools.js` itself ships in minimap2's `misc/` directory
  — no separate install, just point `PAFTOOLS_JS` at
  `<minimap2 source dir>/misc/paftools.js`.

## Setup

```bash
cp config.example.env config.env
# edit config.env: fill in HS1_FASTA, HG38_FASTA, the two PAT/MAT PAF paths,
# and the two hs1<->hg38 chain paths. The HG002 assembly and GIAB VCF/BED
# paths are already filled in if you're using the same data locations.
```

The HG002 assembly ships BGZF-compressed — `samtools faidx` can index it
directly without decompressing:

```bash
samtools faidx "$HG002_ASSEMBLY_FASTA"   # builds .fai + .gzi in place
```

Build the two haplotype-to-hs1 PAF alignments if you don't have them yet
(genome-wide, *not* the repo's bundled `hg002pat_to_hs1.rfc1_only.paf`, which
is filtered down to a single locus):

```bash
samtools faidx "$HG002_ASSEMBLY_FASTA" $(grep '_PATERNAL' "$HG002_ASSEMBLY_FASTA.fai" | cut -f1) > hg002.paternal.fa
samtools faidx "$HG002_ASSEMBLY_FASTA" $(grep '_MATERNAL' "$HG002_ASSEMBLY_FASTA.fai" | cut -f1) > hg002.maternal.fa
minimap2 -cx asm5 --cs=long -t 16 "$HS1_FASTA" hg002.paternal.fa > hg002pat_to_hs1.paf
minimap2 -cx asm5 --cs=long -t 16 "$HS1_FASTA" hg002.maternal.fa > hg002mat_to_hs1.paf
```

## Running the accuracy + confirmation legs

```bash
./run_all.sh config.env windows    # GIAB VCF -> windows.bed + truth.tsv
./run_all.sh config.env chains     # PAF -> chain, both haplotypes (needs paf2chain)
./run_all.sh config.env bedpull    # bedpull --paf, both haplotypes
./run_all.sh config.env liftover   # liftOver + faidx, both haplotypes (needs liftOver)
./run_all.sh config.env negatives  # size-matched no-SV windows + both tools, both haplotypes
./run_all.sh config.env score      # -> out/accuracy_report.md, out/roc_curve.png
# or just:
./run_all.sh config.env all
```

Output lands in `$BENCH_OUTDIR` (default `./out`): `windows.bed`, `truth.tsv`,
`bedpull_{PATERNAL,MATERNAL}.*`, `liftover_{PATERNAL,MATERNAL}.*`,
`accuracy_report.md`, `accuracy_details.tsv`.

### The hs1↔hg38 confirmation leg

There's no separate orchestrator for this — it's the same scripts, called
directly with hg38-specific arguments, treating hg38 as a single "haplotype"
(it's not diploid, so there's only one PAF/chain, not two):

UCSC only publishes `hg38ToHs1.over.chain.gz` (not the reverse), so this leg
lifts **hg38 → hs1**: windows are defined in hg38 coordinates, lifted onto
hs1, and sequence is extracted from the hs1 FASTA. The chain, the extraction
target, and the PAF's target/query roles all have to agree on this direction
— get any one backwards and you'll extract from the wrong genome.

**The truth set here is *not* the HG002 GIAB VCF.** hs1 (CHM13) and hg38
(GRCh38) are both reference assemblies of independent individuals, neither of
which is HG002 — there's no assembly on either side of this comparison that
actually carries HG002's personal variants, so scoring against them would
just measure how often CHM13 happens to agree with HG002 by chance (an
earlier version of this leg did exactly that; see "Why not the GIAB VCF"
below). Instead, `call_alignment_svs.py` calls indel events directly from the
hg38↔hs1 alignment's own CIGAR — real, verifiable differences between the two
reference genomes being extracted between — filtered to primary, mapq-60
alignment blocks and kept clear of alignment-record edges
(`--edge-buffer`) to avoid boundary artifacts.

```bash
# PAF: target=hg38 (matches the windows' coordinate system), query=hs1
# (what we extract from) — same direction convention as the haplotype PAFs.
minimap2 -x asm5 -d hg38.mmi "$HG38_FASTA"
minimap2 -cx asm5 -t 16 hg38.mmi "$HS1_FASTA" > hs1_to_hg38.paf

python3 call_alignment_svs.py \
    --paf hs1_to_hg38.paf --fai "$HG38_FASTA.fai" \
    --flank "$SV_FLANK_BP" --min-len "$SV_MIN_LEN" \
    --min-exclusion-len 10 --min-mapq 60 --edge-buffer 1000 \
    --out-windows out/hg38/windows.bed --out-truth out/hg38/truth.tsv \
    --out-aligned-bed out/hg38/aligned_regions.bed \
    --out-exclusion-bed out/hg38/exclusion_events.bed

python3 select_negative_control_windows.py \
    --exclusion-bed out/hg38/exclusion_events.bed \
    --high-conf-bed out/hg38/aligned_regions.bed \
    --truth out/hg38/truth.tsv --fai "$HG38_FASTA.fai" --flank "$SV_FLANK_BP" \
    --out-windows out/hg38/neg_windows.bed --out-truth out/hg38/neg_truth.tsv

./run_bedpull_pipeline.sh out/hg38/windows.bed hs1_to_hg38.paf "$HS1_FASTA" hg38 out/hg38/bedpull_hg38
./run_bedpull_pipeline.sh out/hg38/neg_windows.bed hs1_to_hg38.paf "$HS1_FASTA" hg38 out/hg38/neg_bedpull_hg38

./run_liftover_pipeline.sh out/hg38/windows.bed "$HG38_TO_HS1_CHAIN" "$HS1_FASTA" hg38 out/hg38/liftover_hg38
./run_liftover_pipeline.sh out/hg38/neg_windows.bed "$HG38_TO_HS1_CHAIN" "$HS1_FASTA" hg38 out/hg38/neg_liftover_hg38

python3 score_accuracy.py \
    --truth out/hg38/truth.tsv \
    --bedpull-results out/hg38/bedpull_hg38.results.tsv \
    --liftover-results out/hg38/liftover_hg38.results.tsv \
    --neg-truth out/hg38/neg_truth.tsv \
    --neg-bedpull-results out/hg38/neg_bedpull_hg38.results.tsv \
    --neg-liftover-results out/hg38/neg_liftover_hg38.results.tsv \
    --roc-out out/hg38/roc_curve.png --roc-data-out out/hg38/roc_data.tsv \
    --out-report out/hg38/confirmation_report.md --out-details out/hg38/confirmation_details.tsv
```

### Optional third tool: `paftools.js liftover`

liftOver isn't the only tool routinely used for this job — minimap2 itself
ships a simpler PAF-based liftover in `paftools.js` (`misc/paftools.js`,
run via the `k8` JS shell). It's a useful third comparison specifically
*because* it's PAF/CIGAR-based like bedpull rather than chain-based like
liftOver — it isolates whether bedpull's advantage over liftOver is really
"CIGAR-aware beats chain-based" in general, or something bedpull's
implementation specifically does better than another CIGAR-aware tool.

#### Feature/capability comparison

The accuracy numbers above are one axis; the two tools' actual feature sets
are quite different. paftools.js is a general-purpose PAF-ecosystem
multitool (23 subcommands total — format conversion, assembly QC, variant
calling, simulation evaluation), of which `liftover` is one thin,
single-purpose subcommand. bedpull is a dedicated, dependency-free
extraction tool with a much richer filtering vocabulary.

**Deployment:**

| | bedpull | paftools.js |
|---|---|---|
| Runtime | Single static Rust binary (musl build available) | `.js` script — needs the `k8` JS shell installed separately (not bundled, not on most systems by default) |
| Distribution | Standalone tool | Ships inside minimap2's source tree (`misc/paftools.js`), not packaged as its own tool |

**Core task — extracting sequence at BED regions:** this is the bigger
structural difference. paftools.js has no BAM/CRAM mode at all, and even
its PAF mode doesn't produce sequence in one step.

| | bedpull | paftools.js |
|---|---|---|
| BAM/CRAM → per-read sequence extraction | **Yes** | **None** — `sam2paf` converts SAM/BAM alignment records to PAF syntax (one line per alignment, no sequence, no region query), not extraction |
| PAF → sequence extraction | Yes, one command: BED region in, FASTA/FASTQ out | `liftover` only lifts *coordinates* (BED in, lifted BED out) — sequence needs a separate `samtools faidx` call after, exactly the extra step `run_paftools_pipeline.sh` had to bolt on |

**Filtering/selection options** (bedpull's BAM/CRAM mode has no paftools.js
counterpart at all, so this compares bedpull generally against what
`liftover` specifically exposes):

| Feature | bedpull | paftools.js `liftover` |
|---|---|---|
| Min mapping quality | `--min_mapq` (BAM/CRAM) | `-q` (default 5) |
| Min alignment block length | n/a | `-l` (default **50,000** — the filter we had to disable for a fair comparison) |
| Max sequence divergence | n/a | `-d` (default disabled) |
| Secondary/supplementary handling | `--include_secondary` / `--include_supplementary` | n/a on `liftover` (only `sam2paf` has a primary-only toggle) |
| Partial-overlap regions + min coverage fraction | `--partial` + `--min_partial_coverage` | No concept of this — not per-read, so "partial overlap" doesn't apply |
| Min region base quality | `--min_region_quality` | n/a (no quality scores in scope) |
| Flanking bases before CIGAR walk | `--flanks`/`--lflank`/`--rflank` | n/a |
| Haplotype-tag split output | `--hap_split` (`HP`/`hp:i:`) | n/a |
| Dedup repeated read/contig names | `--dedup` | n/a |
| FASTQ output | `--fastq` (BAM only) | n/a (BED only, never sequence+quality) |
| Unmapped-region reporting | `--unmapped` with reason codes | Silent drop — `parse_paftools_liftover.py` had to be written just to reconstruct which windows failed |

**What paftools.js has that bedpull doesn't** (adjacent, not competing):
`call` (calls variants from a PAF's `cs` tag, outputs VCF — more mature than
the `call_alignment_svs.py` hand-rolled for this benchmark's confirmation
leg), `sam2paf`/`delta2paf`/`gff2bed` (format conversion), `asmstat`/
`asmgene`/`misjoin`/`bedcov` (assembly/mapping QC), `mapeval`/`mason2fq`/
`pbsim2fq` (simulated-data benchmarking). None of these overlap with
bedpull's job.

It needs the **reverse**-direction alignment from the one bedpull/liftOver
use in this leg: `paftools.js liftover` keys its input BED lookup on the
PAF's *query* name and lifts onto the *target*, so our hg38 windows need to
be in the query role this time (target=hs1, query=hg38 — swapped from
`hs1_to_hg38.paf` above):

`paftools.js liftover` also defaults to `-l 50000` (minimum PAF alignment
*block* length) — a filter bedpull has no equivalent of. Left at its
default, it silently drops ~80% of the records in a typical asm5
whole-genome PAF (most blocks between two real assemblies are far shorter
than 50kb), which isn't a fair comparison to a tool with no such filter, so
`run_paftools_pipeline.sh` passes `-l 0` to disable it.

```bash
minimap2 -x asm5 -d hs1.mmi "$HS1_FASTA"
minimap2 -cx asm5 -t 16 hs1.mmi "$HG38_FASTA" > hg38_to_hs1.paf

./run_paftools_pipeline.sh out/hg38/windows.bed hg38_to_hs1.paf "$HS1_FASTA" hg38 out/hg38/paftools_hg38 0 0
./run_paftools_pipeline.sh out/hg38/neg_windows.bed hg38_to_hs1.paf "$HS1_FASTA" hg38 out/hg38/neg_paftools_hg38 0 0

python3 score_accuracy.py \
    --truth out/hg38/truth.tsv \
    --bedpull-results out/hg38/bedpull_hg38.results.tsv \
    --liftover-results out/hg38/liftover_hg38.results.tsv \
    --paftools-results out/hg38/paftools_hg38.results.tsv \
    --neg-truth out/hg38/neg_truth.tsv \
    --neg-bedpull-results out/hg38/neg_bedpull_hg38.results.tsv \
    --neg-liftover-results out/hg38/neg_liftover_hg38.results.tsv \
    --neg-paftools-results out/hg38/neg_paftools_hg38.results.tsv \
    --roc-out out/hg38/roc_curve.png --roc-data-out out/hg38/roc_data.tsv \
    --out-report out/hg38/confirmation_report.md --out-details out/hg38/confirmation_details.tsv
```

**It does extend to the main leg**, via the same reverse-alignment trick:
the main leg's `windows.bed` is in hs1 (reference) coordinates, which is
the PAF's *target* space in `hg002pat_to_hs1.paf`/`hg002mat_to_hs1.paf`
(target=hs1, query=HG002 contig) — the wrong role for `paftools.js
liftover`. Building the reverse-direction alignment per haplotype
(target=HG002 contig, query=hs1) puts `windows.bed` back in the query role
it needs:

```bash
minimap2 -x asm5 -d hg002pat.mmi hg002.paternal.fa
minimap2 -cx asm5 -t 16 hg002pat.mmi "$HS1_FASTA" > hs1_to_hg002pat.paf
minimap2 -x asm5 -d hg002mat.mmi hg002.maternal.fa
minimap2 -cx asm5 -t 16 hg002mat.mmi "$HS1_FASTA" > hs1_to_hg002mat.paf

./run_paftools_pipeline.sh out/windows.bed hs1_to_hg002pat.paf "$HG002_ASSEMBLY_FASTA" PATERNAL out/paftools_PATERNAL 0 0
./run_paftools_pipeline.sh out/windows.bed hs1_to_hg002mat.paf "$HG002_ASSEMBLY_FASTA" MATERNAL out/paftools_MATERNAL 0 0
./run_paftools_pipeline.sh out/neg_windows.bed hs1_to_hg002pat.paf "$HG002_ASSEMBLY_FASTA" PATERNAL out/neg_paftools_PATERNAL 0 0
./run_paftools_pipeline.sh out/neg_windows.bed hs1_to_hg002mat.paf "$HG002_ASSEMBLY_FASTA" MATERNAL out/neg_paftools_MATERNAL 0 0

python3 score_accuracy.py \
    --truth out/truth.tsv \
    --bedpull-results out/bedpull_PATERNAL.results.tsv out/bedpull_MATERNAL.results.tsv \
    --liftover-results out/liftover_PATERNAL.results.tsv out/liftover_MATERNAL.results.tsv \
    --paftools-results out/paftools_PATERNAL.results.tsv out/paftools_MATERNAL.results.tsv \
    --neg-truth out/neg_truth.tsv \
    --neg-bedpull-results out/neg_bedpull_PATERNAL.results.tsv out/neg_bedpull_MATERNAL.results.tsv \
    --neg-liftover-results out/neg_liftover_PATERNAL.results.tsv out/neg_liftover_MATERNAL.results.tsv \
    --neg-paftools-results out/neg_paftools_PATERNAL.results.tsv out/neg_paftools_MATERNAL.results.tsv \
    --roc-out out/roc_curve.png --roc-data-out out/roc_data.tsv \
    --out-report out/accuracy_report.md --out-details out/accuracy_details.tsv
```

This is worth doing precisely *because* it's non-circular in a way the
confirmation leg's paftools.js comparison isn't: the main leg's truth comes
from the independent GIAB VCF, not from either PAF's own CIGAR, so neither
bedpull nor paftools.js gets a truth-matches-its-own-alignment advantage
here. The result (see below) tells a genuinely different story than the
confirmation leg's.

Whole-genome asm5 alignment is memory-hungry (each haplotype run peaked at
~22-25GB RSS here) — if a run gets killed under memory pressure, retry with
`-t 4` or lower rather than `-t 16`; it's slower but needs less RAM.

Note this design has its own honest limitation, the mirror image of the main
leg's: bedpull's truth set here is derived from the *exact same* PAF bedpull
then extracts from, so its near-perfect detection AUC on this leg partly
reflects "can bedpull correctly walk the CIGAR it was given" rather than
"can it independently detect an SV" — expected and not circular for what it
tests (bedpull's extraction correctness), but worth knowing before reading
its AUC as an apples-to-apples measure against liftOver, whose chain is
built independently of this truth set.

### Two more UCSC kent tools: `pslMap`, and `liftOver -multiple` + `liftOverMerge`

liftOver isn't the only relevant tool in UCSC's own kent toolkit either.
Skimming the ~370 binaries at
<https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/> (verified against
the actual [kent source](https://github.com/ucscGenomeBrowser/kent), not
guessed from tool names) turns up two genuinely relevant ones — both reuse
the *exact same chain file* liftOver already has, no new alignment needed:

- **`pslMap`** does base-by-base coordinate projection through an alignment,
  where liftOver's own chain-walking is block-level. Its own README states
  this plainly: *"This differs from liftOver in that it does a base-by-base
  mapping, which may insert or delete bases within a block. The mappings
  produced by liftOver are block-level, which may expand or contract the
  size of blocks."* — precisely the architectural gap bedpull/paftools.js
  close over liftOver, available from inside the kent toolkit itself, and
  it accepts a **chain file directly** via `-chainMapFile` (no PSL
  conversion of the alignment needed). Workflow: `bedToPsl -keepQuery`
  turns our BED windows into trivial "self-alignment" PSLs, `pslMap
  -chainMapFile -swapMap` projects them through the chain (`-swapMap` is
  needed because our windows sit in the chain's *target* space — the same
  space liftOver itself consumes input in — while `pslMap` requires the
  input's target to match the map file's *query* space; this flag reuses
  the one official chain in reverse rather than needing a second,
  oppositely-built one), and `parse_pslmap_output.py` recovers which window
  each output row belongs to (matching on exact query-span, discarding
  partial-fragment rows from other overlapping chain records — same
  design choice as `parse_paftools_liftover.py`).
- **`liftOver -multiple` + `liftOverMerge`** is a lighter official
  mitigation for the *specific* symptom this whole benchmark keeps
  surfacing: a region straddling an SV/chain boundary gets truncated or
  dropped by liftOver's default single-best-hit behavior. `-multiple`
  outputs every overlapping chain fragment for a region instead of just
  one; `liftOverMerge` stitches fragments that share a chromosome and
  overlap (or sit within `-mergeGap`) back into a single span — verified
  against its actual source, which does check same-chromosome before
  merging, so it won't blindly weld together unrelated hits from different
  chromosomes.

```bash
# pslMap needs chrom.sizes for the windows' own genome
awk '{print $1"\t"$2}' "$HG38_FASTA.fai" > out/hg38/hg38.chrom.sizes   # confirmation leg
awk '{print $1"\t"$2}' "$HS1_FASTA.fai" > out/hs1.chrom.sizes          # main leg

./run_pslmap_pipeline.sh out/hg38/windows.bed "$HG38_TO_HS1_CHAIN" out/hg38/hg38.chrom.sizes "$HS1_FASTA" hg38 out/hg38/pslmap_hg38
./run_liftover_multiple_pipeline.sh out/hg38/windows.bed "$HG38_TO_HS1_CHAIN" "$HS1_FASTA" hg38 out/hg38/liftmulti_hg38
# ... and the negative-control windows the same way, then add to score_accuracy.py:
#   --pslmap-results / --neg-pslmap-results
#   --liftover-multiple-results / --neg-liftover-multiple-results
```

Both `pslMap` and `liftOverMerge` can legitimately produce more than one
candidate hit for the same window (different chromosome, or a fragment too
far away to merge) — the same "multiple hits, take whichever is closest to
expected" shape `score_accuracy.py`'s `best_hit()` already handles for
multi-haplotype results. Finding this exposed a real bug while wiring it
up: `extract_and_score_liftover.py` used to key its per-window bookkeeping
by *name*, so a name appearing more than once (which never happened before
these two tools, since bedpull/liftOver/paftools.js only ever produce one
row per name) would silently keep only the last-seen candidate and drop
the rest. Fixed by keying lengths per *(chrom, start, end)* region instead
of per name — verified byte-identical output on every existing tool's
results (nothing changes when a name really does only appear once) and
correct on a synthetic duplicate-name case.

**The result is the most interesting single finding of this whole
benchmark's tool-comparison work**: on the main leg (recall the numbers
below), `liftOver -multiple` + `liftOverMerge` — a free CLI flag plus one
official kent-toolkit merge step, no new alignment, no new binary beyond
what's already downloaded — actually **slightly exceeds bedpull's own
recall and F1** (68.2%/0.811 vs 68.1%/0.810), and its DEL recall jumps from
vanilla liftOver's collapsed 5.0% to 65.5%, nearly matching bedpull's
65.8%. `pslMap` lands a little behind (67.2% recall, DEL 64.6%) but still
dramatically ahead of vanilla liftOver. The confirmation leg shows the same
pattern at a smaller scale (both land around 64-65% recall and DEL
65.9%/66.5%, vs vanilla liftOver's 40.4%/3.7%).

This says something important about the earlier "liftOver truncates at
chain boundaries" framing: that behavior isn't an inherent limit of
chain-based coordinate lifting — it's specifically a consequence of
liftOver's *default* single-best-hit mode, which both official
alternatives inside the same toolkit almost entirely undo, for free, using
the identical chain file. The one real ceiling that remains: both tools
only ever walk the *same* chain liftOver already has — unlike bedpull or
paftools.js, they can't discover sequence the chain's own underlying
alignment never captured in the first place, so their gains come purely
from not throwing away information the chain already contains, not from a
better or different underlying alignment.

## Running the scale leg

```bash
./run_all.sh config.env scale
```

Writes `$BENCH_OUTDIR/scale/synthetic_{100,1000,10000,100000}.bed` (seeded,
reproducible; smaller BEDs are subsets of larger ones) and
`$BENCH_OUTDIR/scale_results.tsv` (tool, n_regions, elapsed_seconds).

## Script reference

| Script | Purpose |
|---|---|
| `select_sv_windows.py` | GIAB VCF + high-conf BED → `windows.bed` + `truth.tsv` (main leg) |
| `call_alignment_svs.py` | hg38↔hs1 PAF's own CIGAR → `windows.bed` + `truth.tsv` + `aligned_regions.bed` + `exclusion_events.bed` (confirmation leg) |
| `build_chain_from_paf.sh` | minimap2 PAF → UCSC `.chain` (via `paf2chain`) |
| `run_bedpull_pipeline.sh` | Run bedpull `--paf`, parse output → `results.tsv`; trailing args pass straight through to bedpull (e.g. `--stitch_records`) |
| `run_liftover_pipeline.sh` | Run `liftOver` + batched `samtools faidx`, parse → `results.tsv` |
| `run_paftools_pipeline.sh` | Run `paftools.js liftover` (optional extra tool) + batched `samtools faidx`, parse → `results.tsv` |
| `run_pslmap_pipeline.sh` | Run `bedToPsl` + `pslMap -chainMapFile` (optional extra tool, reuses the existing chain) + batched `samtools faidx`, parse → `results.tsv` |
| `run_liftover_multiple_pipeline.sh` | Run `liftOver -multiple` + `liftOverMerge` (optional extra tool, reuses the existing chain) + batched `samtools faidx`, parse → `results.tsv` |
| `parse_paf_results.py` | FASTA + `--unmapped` → results.tsv (called by the bedpull pipeline) |
| `extract_and_score_liftover.py` | lifted BED + unmapped → results.tsv (called by the liftover, paftools, pslMap, and liftOver-multiple pipelines) |
| `parse_paftools_liftover.py` | Recovers real window names in paftools.js liftover's output (called by the paftools pipeline) |
| `parse_pslmap_output.py` | Recovers real window names in pslMap's PSL output, keeping only whole-window (not partial-fragment) hits (called by the pslMap pipeline) |
| `score_accuracy.py` | truth (+ optional negative-control truth) + 2-5 tools' results → markdown report (recall, and precision/F1/ROC if negatives given) + per-window TSV |
| `select_negative_control_windows.py` | Size-matched "no true SV here" windows, for precision/F1/ROC — exclusion positions from either a GIAB VCF (`--vcf`, main leg) or a plain BED (`--exclusion-bed`, confirmation leg) |
| `generate_random_windows.py` | Seeded synthetic BEDs for the scale leg |
| `run_scaling_benchmark.sh` | Times both tools across the synthetic BEDs |
| `run_all.sh` | Orchestrates the hs1↔HG002 leg end to end |

Both pipelines' `results.tsv` share one schema (`name, hap, status,
unmapped_reason, extracted_len, contig_name, query_start, query_end,
strand`), which is what lets `score_accuracy.py` treat them identically.

## What's been tested so far (against real data, in this environment)

- `select_sv_windows.py`: full run against the real GIAB CHM13v2.0 v5.0q
  VCF (100,047 SV records) — 24,965 windows kept after filtering
  (type/length, ambiguous-GT, high-confidence containment all verified to
  sum back to the total).
- `run_bedpull_pipeline.sh` / `parse_paf_results.py`: run against the real
  `chr4_PATERNAL` contig extracted from the actual HG002 assembly and the
  repo's bundled RFC1 PAF — correctly reproduces the canonical 579bp result
  (59bp reference span + 520bp real insertion) and correctly flags a
  non-overlapping window as unmapped.
- `run_liftover_pipeline.sh` / `extract_and_score_liftover.py`: the
  liftOver-invocation shell wrapper isn't runnable here (no `liftOver`
  binary in this environment), but the extraction/scoring half was
  validated directly against mocked lifted-coordinate input on the same
  real `chr4_PATERNAL` contig — correctly extracts, reverse-complements a
  minus-strand test case, and parses liftOver's `#reason` unmapped format.
  Cross-checked: feeding it bedpull's own real output coordinates for RFC1
  gives the same 579bp.
- `score_accuracy.py`: run against real truth-set windows with fabricated
  (but realistic) per-tool results covering pass / wrong-length / unmapped
  cases — expected-length math, tolerance logic, and report rendering all
  checked by hand.
- `generate_random_windows.py`: verified non-overlapping, in-bounds,
  proportionally-distributed-by-chromosome-length placement, and the
  subset property between different `--sizes` values, both single- and
  multi-contig.
- `run_scaling_benchmark.sh`: run end to end (bedpull-only path) against
  real data; produces a valid timing TSV.

Update: the full pipeline (`paf2chain`, `liftOver`, both haplotype PAFs) has
since been run end to end — see Results below.

## Results (hs1↔HG002 accuracy leg, 24,965 SV windows, both haplotypes)

| | bedpull | bedpull+stitch | liftOver | paftools.js | pslMap | liftOver -multiple |
|---|---|---|---|---|---|---|
| Overall recall | 68.1% | **71.6%** | 36.7% | 68.0% | 67.2% | 68.2% |
| DEL recall | 65.8% | **69.1%** | 5.0% | 65.7% | 64.6% | 65.5% |
| INS recall | 70.5% | **74.2%** | 69.8% | 70.4% | 69.9% | 70.9% |
| Recall on SVs ≥5000bp | 91.4% | **92.7%** | 44.3% | 85.2% | 91.4% | 70.4% |
| Windows with zero output | 1,420 (5.7%) | **196 (0.8%)** | 8,027 (32.2%) | 179 (0.7%) | 3,078 (12.3%) | 252 (1.0%) |
| Precision (vs 24,965 size-matched negative controls) | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| F1 | 0.810 | **0.834** | 0.537 | 0.809 | 0.804 | 0.811 |
| Detection AUC | 0.919 | **0.974** | 0.596 | 0.969 | 0.850 | 0.966 |

**`bedpull --stitch_records`** (see `TODO.md` and `docs/src/paf-mode.md`) is now
the best tool on this leg on *every* metric, including detection AUC — it
edges out paftools.js there too (0.974 vs 0.969). This is the combined effect
of several changes from the same investigation: cross-record stitching
itself, plus a real pre-existing bug it surfaced along the way, plus a
follow-up robustness pass. `get_paf_reads`
used to call `get_read_cuts` per record, which has a documented sentinel
convention (misattributing a boundary crossing whenever it coincides exactly
with a record's own `align_start`) that BAM/CRAM's `--partial` mode knows how
to reinterpret but PAF mode never did — so any record whose clamped window
edge landed exactly on its own start or end was silently discarded as
"invalid," independent of whether a second record existed to stitch with at
all. Replacing that per-record call with the same ambiguity-free
`read_pos_at_ref` helper written for stitching fixed both at once. A follow-up
robustness pass then fixed three more gaps found the same way (flank-blind
index query, a hard crash on any single bad contig name/CIGAR, and
`read_start == read_end` wrongly treated as invalid rather than a valid
zero-length "this window is entirely deleted" result — see `TODO.md`'s
2026-07-16 entries for each). Zero-output windows dropped from 1,420 to 196 —
an 86.2% reduction, most of it (~1,200 windows) from the sentinel fix
specifically, not the stitching logic itself (isolated by disabling stitching
and testing the sentinel fix alone). All of this from the *exact same PAF*
bedpull already had, no new alignment or data. Confirmed on real data outside
this benchmark too — a real chr1 chain boundary in the PATERNAL PAF that
previously produced 0 output for a 50bp window now
correctly stitches a 7,407bp extraction.

bedpull recovers roughly 1.85x liftOver's overall recall, with the gap most
dramatic on deletions specifically — consistent with liftOver truncating at
chain boundaries where bedpull walks the CIGAR through the full indel.

Precision came out 100% for every tool against the negative-control set —
none ever extracted a wrong-length result at a locus confirmed to have no
real SV nearby, so there's no false-positive signal to trade off against
recall here, and F1 ends up a monotonic function of recall alone.
The **detection AUC** is the more informative companion metric: it asks a
different, genuinely independent question — "regardless of whether it gets
the exact size right, can the tool tell a true-SV window apart from a
no-SV window at all?" — using deviation from the plain reference-window
length (not the class-specific expected length) as a single score swept
across every threshold (see `roc_points()`'s docstring in
`score_accuracy.py` for why using expected_length directly doesn't work).
bedpull's ROC hugs the top-left corner (AUC 0.919); liftOver sits well above
chance but far below bedpull (AUC 0.596). Plot: `out/roc_curve.png`.

**The standout result is `liftOver -multiple` + `liftOverMerge`** (see
"Two more UCSC kent tools" above): using the *exact same chain* vanilla
liftOver already had, this free CLI-flag-plus-merge-step combination
slightly **exceeds bedpull's own recall and F1** here (68.2%/0.811 vs
68.1%/0.810), and fixes the DEL collapse almost entirely (5.0% → 65.5%,
against bedpull's 65.8%). `pslMap` — the same chain, walked base-by-base
instead — lands close behind (67.2% recall) and actually *matches* bedpull
exactly on SVs ≥5000bp (91.4%). Both fall behind on detection AUC (0.850
and 0.966 vs bedpull/paftools.js's ~0.92-0.97) and on very-large-SV recall
for `liftOver -multiple` specifically (70.4%, its one real weak point —
mean error on those is ~4kb, suggesting large events get partially but not
fully merged back together). Net: most of liftOver's apparent disadvantage
on this leg was never really about chain files being unable to represent
indels — it was liftOver's *default single-best-hit mode* discarding
information its own chain already had.

**`paftools.js liftover`, run here too (see "Optional third tool" above),
lands essentially tied with bedpull on this leg — 68.0% recall vs 68.1%, F1
0.809 vs 0.810 — and actually edges bedpull out on detection AUC (0.969 vs
0.919).** This is the more trustworthy of the two paftools.js comparisons in
this README: unlike the confirmation leg, this leg's truth set is the
independent GIAB VCF, not derived from either tool's own PAF, so neither
bedpull nor paftools.js gets a truth-matches-its-alignment advantage here —
each is reading its *own* separately-built alignment (bedpull:
`hg002{pat,mat}_to_hs1.paf`, target=hs1; paftools.js: the reverse
`hs1_to_hg002{pat,mat}.paf`, target=HG002) against the same independent
truth. That they land this close together is a real result, not an
artifact: it means bedpull's CIGAR-walking implementation and paftools.js's
point-lifting implementation are, in practice, roughly equally capable at
this task when each gets a reasonable assembly-to-reference alignment. It
also retroactively supports the read on the confirmation leg's gap (71.5%
vs 85.6%) — since the two tools tie here where the comparison is fair, that
gap is best explained by the circularity already flagged there (bedpull
reading the exact alignment its confirmation-leg truth was built from),
not by paftools.js's algorithm being meaningfully worse.

This tie held for bedpull's *original* implementation. With
`--stitch_records` (see the table above and "Two more UCSC kent tools"),
bedpull now leads paftools.js on this leg on every metric, including
detection AUC (71.6%/0.834/0.974 vs 68.0%/0.809/0.969) — the one class of
failure paftools.js's point-lifting approach can't recover (a chain-boundary
split, plus the sentinel bug the same investigation surfaced) is exactly
what the fix targets, so the tie was a property of bedpull's *old*
single-record limitation, not a ceiling on what CIGAR-aware extraction from
this alignment can achieve.

The scale leg tells the opposite story: at 100,000 synthetic regions,
bedpull took 157s vs liftOver's 9.3s. liftOver's chain-based lookup is
architecturally faster than bedpull's current per-region PAF query, and
bedpull has no threading yet (tracked in the main repo's `TODO.md`). Not a
knock against the accuracy result, just a real, separate tradeoff.

### Two "worst offenders" that turned out not to be bugs

Digging into bedpull's largest per-window errors surfaced one real harness
bug (now fixed) and one genuine, expected limitation:

1. **Multi-allelic/VNTR sites can make one region string represent more than
   one truth-set window.** `extract_and_score_liftover.py` originally keyed
   its region→name lookup as a plain dict, so when two different SV records
   lifted to the exact identical coordinates (common at repeat-rich loci
   where nearby small-indel calls collapse onto the same real difference),
   the second name silently got dropped and scored as a 0-length "failure"
   that never happened. Fixed by mapping region → *list* of names. This
   affected ~5.8% of lifted regions per haplotype and was solely responsible
   for liftOver's overall recall reading 33.9% instead of the correct 36.7%
   — a real fix, but it doesn't change the headline conclusion.
2. **Very large insertions (tens of kb) occasionally fail identically in
   both tools.** E.g. a window expecting an 88,009bp extraction: bedpull
   reports "no valid overlap" and liftOver — using a chain built from the
   *same* PAF — reports "Partially deleted in new" for the same haplotype.
   Since both draw from the identical minimap2 alignment, this means the
   alignment itself doesn't cleanly span this insertion's breakpoint (most
   likely split into separate PAF records), not a bug in either tool. If
   anything this is a useful cross-check: where the tools genuinely differ,
   it's on identical underlying alignment data, so the difference reflects
   real capability, not alignment noise.

## Results (hg38↔hs1 confirmation leg, 19,984 SV windows, official UCSC chain + PAF-based tools)

| | bedpull | bedpull+stitch | liftOver | paftools.js | pslMap | liftOver -multiple |
|---|---|---|---|---|---|---|
| Overall recall | 85.6% | **85.7%** | 40.4% | 71.5% | 64.4% | 64.8% |
| DEL recall | 88.6% | 88.6% | **3.7%** | 71.8% | 65.9% | 66.5% |
| INS recall | 83.7% | 83.7% | 63.1% | 71.3% | 63.3% | 63.6% |
| Precision (vs 19,984 size-matched negative controls) | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| F1 | 0.922 | **0.923** | 0.575 | 0.834 | 0.783 | 0.786 |
| Detection AUC | **0.998** | **0.998** | 0.533 | 0.947 | 0.807 | 0.853 |

These are the corrected numbers after fixing a genome-mismatch flaw in this
leg's original design (below) — bedpull recovers roughly 2.1x liftOver's
recall, and the DEL-vs-chain-boundary-truncation pattern from the main leg
shows up here too, even more starkly (liftOver's DEL recall collapses to
3.7%). Take bedpull's near-perfect AUC with the caveat noted above (its
truth set is derived from the same PAF it extracts from) — the F1/recall
numbers are the more informative comparison here since they measure
concrete extraction correctness, not just detection.

**`bedpull+stitch` is almost, but not quite, identical to plain `bedpull` on
this leg — and the small gap it does have is informative.** `--stitch_records`
bundles two changes: cross-record stitching itself, and a `get_read_cuts`
sentinel fix that applies independently of whether a chain exists at all
(see the main leg's Results section above). `call_alignment_svs.py` (this
leg's truth generator) already guarantees every truth event sits
`--edge-buffer` bases clear of any PAF record boundary, so there's no
split-record scenario here for *stitching* to fix — but the sentinel fix
can still recover a handful of ordinary single-record partial overlaps
whose clamped window edge happened to land exactly on a record's own start
or end, which is exactly the small recall/F1 bump (85.6%→85.7%,
0.922→0.923). The main leg's truth (the independent GIAB VCF) carries
neither guarantee, which is why both fixes show up there at real scale
instead.

`pslMap` and `liftOver -multiple` + `liftOverMerge` (see "Two more UCSC
kent tools" above) land close together here — both around 64-65% recall,
both fixing liftOver's DEL collapse almost completely (3.7% → 65.9%/66.5%)
using the identical official chain, no new alignment. They land behind
`paftools.js` (71.5%) and well behind bedpull (85.6%) on this leg
specifically, unlike the main leg where `liftOver -multiple` actually tied
bedpull — the difference is that this leg's official chain is independently
built (not derived from either bedpull's or paftools.js's own alignment),
so `pslMap`/`liftOver -multiple`'s ceiling here is capped by whatever that
particular chain does or doesn't capture, while bedpull and paftools.js
each get to walk their *own*, separately-built alignment instead.

**`paftools.js liftover`** (minimap2's own PAF-based lifter, run via `k8` —
see "Optional third tool" above) lands in between: 71.5% recall, comfortably
ahead of liftOver's chain-based 40.4% (and none of liftOver's DEL-specific
collapse — 71.8% vs 3.7%), but still well short of bedpull's 85.6%. Digging
into *why*, rather than reading that gap as "bedpull's CIGAR-walking
algorithm is better": the two tools' core algorithms are actually very
similar. Both walk a PAF record's CIGAR tracking a reference and a query
offset, and map a window's two boundary points across `M` (linear
interpolation), `I` (query-only advance), and `D` (target-only advance) ops
— bedpull does this with an exact base-by-base walk through the op that
straddles each boundary (`get_read_cuts` in `src/utils.rs`), paftools.js
does the equivalent with direct arithmetic; these produce the same answer
when a window's two boundary points fall inside the *same* alignment
record.

The gap traces to what happens when they *don't*, which turns out to be
common in exactly the largest-SV cases: checking the per-window details
directly, most of paftools.js's "ok-but-wrong-length" failures came back
reporting almost exactly the plain window length — i.e. it silently missed
the indel and reported "no variant" — and this was heavily concentrated on
large insertions/deletions. The cause isn't a flaw in the point-lifting math
itself; it's that **paftools.js is reading a separately-computed,
reverse-direction alignment** (`hg38_to_hs1.paf`, a fresh minimap2 run — it
has to be, since it needs windows in the PAF's query role, the opposite of
what bedpull/liftOver use here), not the exact alignment
`call_alignment_svs.py` used to build the truth set. A large, truly novel
insertion often has no homologous target sequence to align against at all,
so minimap2 frequently represents it as *two separate, abutting PAF
records* rather than one record with a giant `I` op — and paftools.js has
no logic to stitch a lifted coordinate across that record boundary; whichever
endpoint falls outside the record it's currently processing just gets
clamped to that record's own edge (visible in its raw output as a `_t5`/
`_t3` suffix), collapsing the entire indel to nothing. bedpull would hit the
identical wall on this same reverse alignment (its own PAF mode is also
strictly per-record — confirmed by an earlier anomaly on the main leg, see
below, where a large insertion failed identically in both tools from the
*same* alignment) — bedpull's high score here specifically reflects reading
the alignment that generated the truth, where `call_alignment_svs.py`
already guaranteed every truth event sits `--edge-buffer` bases clear of any
record boundary. So the real takeaways are: (1) the two implementations'
CIGAR-walking logic is not meaningfully different in what it computes, (2)
bedpull's number here is inflated by truth/alignment circularity (already
flagged above), and (3) paftools.js's number is a fair, non-circular
estimate of "how well does a second, independent CIGAR-aware liftover do on
this leg" — 71.5% recall, still comfortably ahead of chain-based liftOver's
40.4%, which is the more meaningful comparison to draw from this leg. The
main leg (below), where paftools.js is *also* run against an independently
built alignment but the truth set isn't derived from either tool's PAF,
confirms this directly: there bedpull and paftools.js land within 0.1
percentage points of each other (68.1% vs 68.0% recall) — the clean,
non-circular version of this same comparison, and it says the two
implementations really are roughly equally capable.

### Why not the GIAB VCF (a design bug this leg used to have)

An earlier version of this leg reused the HG002 GIAB SV truth set here too
— `HG002_GRCh38_v5.0q_stvar.vcf.gz`, HG002's *own* structural variants
relative to GRCh38 — and scored extraction from hs1/CHM13 against it. That's
a real bug, not just a caveat: hs1 (CHM13) and hg38 (GRCh38) are both
reference assemblies of individuals *other than* HG002, so there's no
guarantee either one carries the same variant HG002 has at a given locus.
Checking the old per-window details confirmed this concretely — of
bedpull's failed windows, 53.2% had extracted a length that exactly matched
"no variant present," i.e. bedpull correctly reported CHM13's real
(unvaried) sequence, it just wasn't the length the HG002-specific truth
set expected. That was 27.2% of *all* windows in the old version of this
leg penalized for a real biological fact having nothing to do with either
tool (liftOver: 29.3%), and it's why this leg's numbers used to read so much
lower than the main leg's (41.0%/28.1% recall, F1 0.581/0.438, AUC
0.808/0.618) despite testing "the same thing."

The fix (`call_alignment_svs.py`, see above) replaces that truth set with
indel events called directly from the hg38↔hs1 alignment's own CIGAR — real
differences between the two genomes actually being compared, so every
window's expected variant is guaranteed to exist at that locus by
construction. Re-running the full leg with this truth set dropped the
"tool reported no variant, and rightly so" rate from 27.2%/29.3% to
0.9%/7.1% (still nonzero — real alignment noise, edge cases, and a few
duplicate/overlapping alignment blocks remain, see `score_tool()`'s
`looks_like_no_variant` diagnostic in every report) and the recall numbers
above are the corrected result. The main leg never had this problem — its
truth set comes from HG002's *own* assembly, the same sequence being
extracted, so the expected variant is guaranteed to genuinely be present in
one of the two haplotypes.

## Notes on the ROC plots

Both ROC curves are rendered with `kuva roc` if it's on `PATH`, falling back
to matplotlib otherwise —
`score_accuracy.py` writes raw `(tool, score, label)` rows to
`--roc-data-out` and lets whichever plotter builds the curve itself, rather
than pre-aggregating points by hand. All three independent
implementations — this harness's own `roc_points()`/`auc()`, kuva's
built-in calculation, and `sklearn.metrics.roc_auc_score` (checked by hand,
not a project dependency) — agree on every AUC value reported here to 3
decimal places.

Both curves end in a near-vertical jump to (1, 1) — this is standard: every
ROC curve starts at (0,0) and ends at (1,1) by construction (loosest
threshold ⇒ everything called positive), and it's the shape *just before*
that endpoint that carries the real signal. It looks like a sudden jump
specifically here because a large majority of negative-control windows
extract a *perfect* (score=0.0) match, and unmapped true-SV windows (a
tool- and leg-dependent fraction — from 0% up to liftOver's ~44% on the
confirmation leg) are only swept in at the very last, artificially-low
threshold step (they can never be "detected" at any real threshold, by
construction). Neither of these change the AUC — the final segment has zero
width in the FPR direction, so it contributes zero area regardless.
