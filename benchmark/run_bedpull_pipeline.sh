#!/usr/bin/env bash
# Run bedpull --paf mode over a set of windows for one HG002 haplotype, and
# turn the output into a results.tsv comparable to the liftOver pipeline's.
#
# Usage:
#   run_bedpull_pipeline.sh <windows.bed> <paf> <query_ref.fa> <hap_tag> <out_prefix> [extra bedpull args...]
#
# Extra trailing args are passed straight through to bedpull — e.g.
# --stitch_records --max_stitch_gap 10000 to enable cross-record stitching.
#
# Writes:
#   <out_prefix>.fasta            raw bedpull FASTA output
#   <out_prefix>.unmapped.bed     bedpull's --unmapped report
#   <out_prefix>.bed_out.bed      lifted query coordinates (BED6)
#   <out_prefix>.results.tsv      one row per (window, hit) — see parse_paf_results.py
#   <out_prefix>.timing.tsv       hap, n_windows, elapsed_seconds

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "usage: $(basename "$0") <windows.bed> <paf> <query_ref.fa> <hap_tag> <out_prefix> [extra bedpull args...]" >&2
    exit 1
fi

windows_bed="$1"
paf="$2"
query_ref="$3"
hap_tag="$4"
out_prefix="$5"
shift 5
extra_args=("$@")

: "${BEDPULL_BIN:=bedpull}"

for f in "$windows_bed" "$paf" "$query_ref"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(wc -l < "$windows_bed")

start_ts=$(date +%s.%N)
"$BEDPULL_BIN" \
    --paf "$paf" \
    --query_ref "$query_ref" \
    --bed "$windows_bed" \
    --output "${out_prefix}.fasta" \
    --unmapped "${out_prefix}.unmapped.bed" \
    --bed_out "${out_prefix}.bed_out.bed" \
    "${extra_args[@]}" \
    >"${out_prefix}.stderr.log" 2>&1
end_ts=$(date +%s.%N)
elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

printf "hap\tn_windows\telapsed_seconds\n%s\t%s\t%s\n" \
    "$hap_tag" "$n_windows" "$elapsed" >"${out_prefix}.timing.tsv"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "$script_dir/parse_paf_results.py" \
    --fasta "${out_prefix}.fasta" \
    --unmapped "${out_prefix}.unmapped.bed" \
    --hap "$hap_tag" \
    --out "${out_prefix}.results.tsv"

echo "bedpull ($hap_tag): $n_windows windows in ${elapsed}s -> ${out_prefix}.results.tsv" >&2
