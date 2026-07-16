#!/usr/bin/env bash
# Run liftOver + samtools faidx over a set of windows through a chain file,
# and turn the output into a results.tsv comparable to the bedpull pipeline's
# (see parse_paf_results.py / extract_and_score_liftover.py for the shared
# schema).
#
# Usage:
#   run_liftover_pipeline.sh <windows.bed> <chain> <target_fasta> <hap_tag> <out_prefix>
#
# Writes:
#   <out_prefix>.windows6.bed   BED6 version of the input (liftOver needs a
#                               strand column to report chain-flips)
#   <out_prefix>.lifted.bed     liftOver's successfully-lifted coordinates
#   <out_prefix>.unmapped.bed   liftOver's own '#reason' + BED unmapped report
#   <out_prefix>.results.tsv    one row per window
#   <out_prefix>.timing.tsv     hap, n_windows, liftover/extract/total seconds

set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "usage: $(basename "$0") <windows.bed> <chain> <target_fasta> <hap_tag> <out_prefix>" >&2
    exit 1
fi

windows_bed="$1"
chain="$2"
target_fasta="$3"
hap_tag="$4"
out_prefix="$5"

: "${LIFTOVER_BIN:=liftOver}"
: "${SAMTOOLS_BIN:=samtools}"

for f in "$windows_bed" "$chain" "$target_fasta"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

if ! command -v "$LIFTOVER_BIN" >/dev/null 2>&1; then
    cat >&2 <<EOF
error: liftOver binary not found (checked: $LIFTOVER_BIN).
Download the standalone binary for your platform from
https://hgdownload.soe.ucsc.edu/admin/exe/ (no build/install needed — just
chmod +x and put it on PATH, or set LIFTOVER_BIN).
EOF
    exit 1
fi

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(wc -l < "$windows_bed")

# liftOver needs at least a BED6 to report strand flips; our windows are
# unstranded BED4, so treat every input as '+' and let liftOver tell us
# whether the chain flips it — that's what determines whether the extracted
# sequence needs reverse-complementing.
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,"+"}' "$windows_bed" >"${out_prefix}.windows6.bed"

start_ts=$(date +%s.%N)
"$LIFTOVER_BIN" "${out_prefix}.windows6.bed" "$chain" "${out_prefix}.lifted.bed" "${out_prefix}.unmapped.bed"
end_ts=$(date +%s.%N)
lift_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

echo "liftOver ($hap_tag): $n_windows windows in ${total_elapsed}s -> ${out_prefix}.results.tsv" >&2
