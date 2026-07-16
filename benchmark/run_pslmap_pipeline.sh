#!/usr/bin/env bash
# Run UCSC kent's pslMap over a set of windows through an existing chain
# file, doing base-by-base (not block-level) coordinate projection — see
# README.md's "Other UCSC kent tools" section for why this is architecturally
# closer to bedpull/paftools.js than liftOver's own block-level chain walk.
#
# Reuses the *same* chain file liftOver would use for these windows — no
# separate alignment needed. -swapMap is required because our windows are
# always in the chain's target ("t") space (the same space liftOver takes
# input in), while pslMap requires the input PSL's target to match the
# mapFile's *query* ("q") space; -swapMap tells pslMap to use the chain in
# reverse without needing a second, differently-built chain file.
#
# Usage:
#   run_pslmap_pipeline.sh <windows.bed> <chain> <chrom_sizes> <target_fasta> <hap_tag> <out_prefix>
#
# Writes:
#   <out_prefix>.windows.psl    bedToPsl -keepQuery output (trivial self-alignment)
#   <out_prefix>.mapped.psl     pslMap's output
#   <out_prefix>.lifted.bed     recovered-name BED6 (see parse_pslmap_output.py)
#   <out_prefix>.unmapped.bed   windows pslMap never lifted end-to-end
#   <out_prefix>.results.tsv    one row per window
#   <out_prefix>.timing.tsv     hap, n_windows, map/extract/total seconds

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "usage: $(basename "$0") <windows.bed> <chain> <chrom_sizes> <target_fasta> <hap_tag> <out_prefix>" >&2
    exit 1
fi

windows_bed="$1"
chain="$2"
chrom_sizes="$3"
target_fasta="$4"
hap_tag="$5"
out_prefix="$6"

: "${BEDTOPSL_BIN:=bedToPsl}"
: "${PSLMAP_BIN:=pslMap}"
: "${SAMTOOLS_BIN:=samtools}"

for f in "$windows_bed" "$chain" "$chrom_sizes" "$target_fasta"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(wc -l < "$windows_bed")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

start_ts=$(date +%s.%N)
"$BEDTOPSL_BIN" -keepQuery "$chrom_sizes" "$windows_bed" "${out_prefix}.windows.psl"
"$PSLMAP_BIN" -chainMapFile -swapMap "${out_prefix}.windows.psl" "$chain" "${out_prefix}.mapped.psl" 2>"${out_prefix}.pslmap.log"
end_ts=$(date +%s.%N)
map_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

python3 "$script_dir/parse_pslmap_output.py" \
    --windows "$windows_bed" \
    --mapped-psl "${out_prefix}.mapped.psl" \
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

total_elapsed=$(awk -v a="$map_elapsed" -v b="$extract_elapsed" 'BEGIN { printf "%.3f", a + b }')
printf "hap\tn_windows\tmap_seconds\textract_seconds\telapsed_seconds\n%s\t%s\t%s\t%s\t%s\n" \
    "$hap_tag" "$n_windows" "$map_elapsed" "$extract_elapsed" "$total_elapsed" >"${out_prefix}.timing.tsv"

echo "pslMap ($hap_tag): $n_windows windows in ${total_elapsed}s -> ${out_prefix}.results.tsv" >&2
