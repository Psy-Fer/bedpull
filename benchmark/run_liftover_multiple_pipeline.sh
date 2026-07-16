#!/usr/bin/env bash
# Run liftOver in -multiple mode (emit every overlapping chain fragment for
# a region, instead of just the single best one) and stitch same-name
# fragments back together with liftOverMerge — a lighter, official-UCSC-
# toolchain mitigation for the "region straddles a chain/SV boundary and
# gets truncated" symptom the main liftOver pipeline demonstrates.
#
# liftOverMerge only merges fragments that share a chromosome and overlap
# (or sit within -mergeGap); genuinely separate candidate loci for the same
# window are left as separate output rows, one per name — exactly the
# "multiple hits, pick the best" shape score_accuracy.py's best_hit() (and
# extract_and_score_liftover.py, since its recent refactor) already handle
# for multi-haplotype results.
#
# Usage:
#   run_liftover_multiple_pipeline.sh <windows.bed> <chain> <target_fasta> <hap_tag> <out_prefix> [min_match]
#
# Writes:
#   <out_prefix>.windows6.bed   BED6 version of the input
#   <out_prefix>.multi.bed      liftOver -multiple's raw fragment output
#   <out_prefix>.merged.bed     liftOverMerge's stitched-fragment output
#   <out_prefix>.unmapped.bed   liftOver's own unmapped report
#   <out_prefix>.results.tsv    one row per (window, surviving fragment)
#   <out_prefix>.timing.tsv     hap, n_windows, liftover/merge/extract/total seconds

set -euo pipefail

if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
    echo "usage: $(basename "$0") <windows.bed> <chain> <target_fasta> <hap_tag> <out_prefix> [min_match]" >&2
    exit 1
fi

windows_bed="$1"
chain="$2"
target_fasta="$3"
hap_tag="$4"
out_prefix="$5"
min_match="${6:-0.1}"

: "${LIFTOVER_BIN:=liftOver}"
: "${LIFTOVERMERGE_BIN:=liftOverMerge}"
: "${SAMTOOLS_BIN:=samtools}"

for f in "$windows_bed" "$chain" "$target_fasta"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(wc -l < "$windows_bed")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,"+"}' "$windows_bed" >"${out_prefix}.windows6.bed"

start_ts=$(date +%s.%N)
"$LIFTOVER_BIN" -multiple -minMatch="$min_match" "${out_prefix}.windows6.bed" "$chain" "${out_prefix}.multi.bed" "${out_prefix}.unmapped.bed"
end_ts=$(date +%s.%N)
lift_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

start_ts=$(date +%s.%N)
"$LIFTOVERMERGE_BIN" "${out_prefix}.multi.bed" "${out_prefix}.merged.bed"
end_ts=$(date +%s.%N)
merge_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

# liftOverMerge's output is BED5 (no strand) — fine here since scoring only
# ever compares extracted *length*, never sequence content/orientation, so
# defaulting every fragment to '+' doesn't affect any reported number.
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"+"}' "${out_prefix}.merged.bed" >"${out_prefix}.merged6.bed"

start_ts=$(date +%s.%N)
python3 "$script_dir/extract_and_score_liftover.py" \
    --lifted "${out_prefix}.merged6.bed" \
    --unmapped "${out_prefix}.unmapped.bed" \
    --target-fasta "$target_fasta" \
    --hap "$hap_tag" \
    --samtools "$SAMTOOLS_BIN" \
    --out "${out_prefix}.results.tsv"
end_ts=$(date +%s.%N)
extract_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

total_elapsed=$(awk -v a="$lift_elapsed" -v b="$merge_elapsed" -v c="$extract_elapsed" 'BEGIN { printf "%.3f", a + b + c }')
printf "hap\tn_windows\tliftover_seconds\tmerge_seconds\textract_seconds\telapsed_seconds\n%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$hap_tag" "$n_windows" "$lift_elapsed" "$merge_elapsed" "$extract_elapsed" "$total_elapsed" >"${out_prefix}.timing.tsv"

echo "liftOver -multiple + liftOverMerge ($hap_tag): $n_windows windows in ${total_elapsed}s -> ${out_prefix}.results.tsv" >&2
