#!/usr/bin/env bash
# Time bedpull (and, if a chain is given, liftOver) against synthetic BED
# files of increasing size, appending results to a single TSV.
#
# Usage:
#   run_scaling_benchmark.sh <bed_dir> <paf> <query_ref> <out_tsv> [<chain> <target_fasta>]
#
# <bed_dir> must contain files named synthetic_<N>.bed (as written by
# generate_random_windows.py). Runs each size once; re-run and concatenate
# if you want repeats for variance.

set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "usage: $(basename "$0") <bed_dir> <paf> <query_ref> <out_tsv> [<chain> <target_fasta>]" >&2
    exit 1
fi

bed_dir="$1"
paf="$2"
query_ref="$3"
out_tsv="$4"
chain="${5:-}"
target_fasta="${6:-}"

: "${BEDPULL_BIN:=bedpull}"
: "${LIFTOVER_BIN:=liftOver}"
: "${SAMTOOLS_BIN:=samtools}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$(dirname "$out_tsv")"
tmp_work="$(mktemp -d)"
trap 'rm -rf "$tmp_work"' EXIT

echo -e "tool\tn_regions\telapsed_seconds" >"$out_tsv"

shopt -s nullglob
bed_files=("$bed_dir"/synthetic_*.bed)
if [ "${#bed_files[@]}" -eq 0 ]; then
    echo "error: no synthetic_*.bed files found in $bed_dir" >&2
    exit 1
fi

# Sort by region count (the numeric part of the filename) so results.tsv
# reads smallest-to-largest regardless of glob order.
IFS=$'\n' bed_files=($(printf '%s\n' "${bed_files[@]}" | sed -E 's/.*synthetic_([0-9]+)\.bed/\1 &/' | sort -n | cut -d' ' -f2-))
unset IFS

for bed in "${bed_files[@]}"; do
    n=$(wc -l <"$bed")

    start_ts=$(date +%s.%N)
    "$BEDPULL_BIN" --paf "$paf" --query_ref "$query_ref" --bed "$bed" \
        --output "$tmp_work/bedpull_out.fasta" >"$tmp_work/bedpull.log" 2>&1
    end_ts=$(date +%s.%N)
    elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')
    echo -e "bedpull\t$n\t$elapsed" >>"$out_tsv"
    echo "bedpull: $n regions in ${elapsed}s" >&2

    if [ -n "$chain" ] && [ -n "$target_fasta" ]; then
        if ! command -v "$LIFTOVER_BIN" >/dev/null 2>&1; then
            echo "skipping liftOver timing: $LIFTOVER_BIN not found" >&2
            continue
        fi
        start_ts=$(date +%s.%N)
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,"+"}' "$bed" >"$tmp_work/windows6.bed"
        "$LIFTOVER_BIN" "$tmp_work/windows6.bed" "$chain" "$tmp_work/lifted.bed" "$tmp_work/unmapped.bed" 2>>"$tmp_work/liftover.log"
        python3 "$script_dir/extract_and_score_liftover.py" \
            --lifted "$tmp_work/lifted.bed" \
            --unmapped "$tmp_work/unmapped.bed" \
            --target-fasta "$target_fasta" \
            --hap scaling \
            --samtools "$SAMTOOLS_BIN" \
            --out "$tmp_work/liftover_results.tsv" >>"$tmp_work/liftover.log" 2>&1
        end_ts=$(date +%s.%N)
        elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')
        echo -e "liftover\t$n\t$elapsed" >>"$out_tsv"
        echo "liftOver: $n regions in ${elapsed}s" >&2
    fi
done

echo "wrote $out_tsv" >&2
