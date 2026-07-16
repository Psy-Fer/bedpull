#!/usr/bin/env bash
# Run `paftools.js liftover` (minimap2's own PAF-based liftover, a k8 script)
# over a set of windows, and turn the output into a results.tsv comparable
# to the bedpull/liftOver pipelines' (see parse_paftools_liftover.py /
# extract_and_score_liftover.py for the shared schema).
#
# Unlike UCSC liftOver's chain files, paftools.js liftover reads a PAF
# directly, so it needs the *reverse*-direction alignment from the one
# bedpull/liftOver use: windows must be defined in the PAF's *query* space
# (paftools.js keys its input BED lookup on PAF column 1, the query name),
# and it lifts onto the PAF's target space. For the hg38<->hs1 confirmation
# leg (windows in hg38 coordinates, extracting from hs1) that means a PAF
# built with target=hs1, query=hg38 — the opposite of the PAF bedpull/
# liftOver use in that same leg.
#
# Note: paftools.js liftover defaults to -l 50000 (minimum PAF alignment
# *block* length) — a filter bedpull has no equivalent of. Left at its
# default, it silently drops ~80% of records in a typical asm5 whole-genome
# PAF (most alignment blocks between two real assemblies are much shorter
# than 50kb), which isn't a fair comparison to a tool with no such filter.
# We pass -l 0 below to disable it.
#
# Usage:
#   run_paftools_pipeline.sh <windows.bed> <paf> <target_fasta> <hap_tag> <out_prefix> [min_mapq] [min_aln_len]
#
# Writes:
#   <out_prefix>.raw_lifted.bed   paftools.js liftover's raw stdout
#   <out_prefix>.lifted.bed       recovered-name BED6 (see parse_paftools_liftover.py)
#   <out_prefix>.unmapped.bed     windows paftools.js never lifted
#   <out_prefix>.results.tsv      one row per window
#   <out_prefix>.timing.tsv       hap, n_windows, liftover/extract/total seconds

set -euo pipefail

if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
    echo "usage: $(basename "$0") <windows.bed> <paf> <target_fasta> <hap_tag> <out_prefix> [min_mapq] [min_aln_len]" >&2
    exit 1
fi

windows_bed="$1"
paf="$2"
target_fasta="$3"
hap_tag="$4"
out_prefix="$5"
min_mapq="${6:-0}"
min_aln_len="${7:-0}"

: "${K8_BIN:=k8}"
: "${PAFTOOLS_JS:=paftools.js}"
: "${SAMTOOLS_BIN:=samtools}"

for f in "$windows_bed" "$paf" "$target_fasta"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

if ! command -v "$K8_BIN" >/dev/null 2>&1 && [ ! -x "$K8_BIN" ]; then
    echo "error: k8 binary not found (checked: $K8_BIN). Grab a release from" \
         "https://github.com/attractivechaos/k8/releases, chmod +x the" \
         "k8-x86_64-Linux binary, and point K8_BIN at it." >&2
    exit 1
fi
if [ ! -f "$PAFTOOLS_JS" ]; then
    echo "error: paftools.js not found at: $PAFTOOLS_JS (ships in minimap2's misc/ dir)" >&2
    exit 1
fi

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(wc -l < "$windows_bed")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

start_ts=$(date +%s.%N)
"$K8_BIN" "$PAFTOOLS_JS" liftover -q "$min_mapq" -l "$min_aln_len" "$paf" "$windows_bed" >"${out_prefix}.raw_lifted.bed"
end_ts=$(date +%s.%N)
lift_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

python3 "$script_dir/parse_paftools_liftover.py" \
    --windows "$windows_bed" \
    --raw-lifted "${out_prefix}.raw_lifted.bed" \
    --out-lifted "${out_prefix}.lifted.bed" \
    --out-unmapped "${out_prefix}.unmapped.bed"

start_ts=$(date +%s.%N)
python3 "$script_dir/extract_and_score_liftover.py" \
    --lifted "${out_prefix}.lifted.bed" \
    --unmapped "${out_prefix}.unmapped.bed" \
    --target-fasta "$target_fasta" \
    --hap "$hap_tag" \
    --samtools "$SAMTOOLS_BIN" \
    --out "${out_prefix}.results.tsv"
end_ts=$(date +%s.%N)
extract_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

total_elapsed=$(awk -v a="$lift_elapsed" -v b="$extract_elapsed" 'BEGIN { printf "%.3f", a + b }')
printf "hap\tn_windows\tliftover_seconds\textract_seconds\telapsed_seconds\n%s\t%s\t%s\t%s\t%s\n" \
    "$hap_tag" "$n_windows" "$lift_elapsed" "$extract_elapsed" "$total_elapsed" >"${out_prefix}.timing.tsv"

echo "paftools.js liftover ($hap_tag): $n_windows windows in ${total_elapsed}s -> ${out_prefix}.results.tsv" >&2
