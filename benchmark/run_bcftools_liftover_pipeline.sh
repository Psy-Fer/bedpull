#!/usr/bin/env bash
# Run bcftools +liftover over a set of windows through a chain file and turn
# the output into a results.tsv comparable to every other tool in this
# harness.
#
# bcftools +liftover lifts VCF records, not BED intervals, so each window is
# expressed as two single-base anchor records (its first and last base) by
# make_anchor_vcf.py, lifted, then paired back into an interval by
# bcftools_lifted_to_bed.py. See make_anchor_vcf.py's module docstring for
# why the windows are anchored rather than the SV records themselves being
# lifted directly — lifting the SV record would reconstruct the answer from
# the allele the tool was handed and measure nothing.
#
# The final extract-and-tabulate step reuses extract_and_score_liftover.py
# unchanged, so this tool is measured with identical extraction code to
# liftOver/pslMap/liftOverMerge.
#
# Usage:
#   run_bcftools_liftover_pipeline.sh <windows.bed> <chain> <src_fasta> <target_fasta> <hap_tag> <out_prefix>
#
# Environment:
#   BCFTOOLS_BIN                   bcftools >= 1.20 with the liftover plugin
#   BCFTOOLS_PLUGINS               directory containing liftover.so
#   SAMTOOLS_BIN                   samtools
#   BCFTOOLS_LIFTOVER_ALT_STYLE    ref-only (default) | snv
#   BCFTOOLS_LIFTOVER_EXTRA_OPTS   extra plugin options, appended verbatim after '--'
#                                  (e.g. "--max-indel-inc 5000"; note that has no
#                                  effect on single-base anchors and only matters
#                                  if you adapt this to lift real indel records)
#
# Writes:
#   <out_prefix>.anchors.vcf     two anchor records per window, sorted
#   <out_prefix>.lifted.vcf      +liftover output, sorted
#   <out_prefix>.reject.vcf      +liftover --reject output (reason in FILTER)
#   <out_prefix>.lifted.bed      BED6 of paired-up lifted intervals
#   <out_prefix>.unmapped.bed    liftOver-style '#reason' + BED report
#   <out_prefix>.results.tsv     one row per window
#   <out_prefix>.timing.tsv      hap, n_windows, anchor/lift/extract/total seconds

set -euo pipefail

if [ "$#" -ne 6 ]; then
    echo "usage: $(basename "$0") <windows.bed> <chain> <src_fasta> <target_fasta> <hap_tag> <out_prefix>" >&2
    exit 1
fi

windows_bed="$1"
chain="$2"
src_fasta="$3"
target_fasta="$4"
hap_tag="$5"
out_prefix="$6"

: "${BCFTOOLS_BIN:=bcftools}"
: "${SAMTOOLS_BIN:=samtools}"
: "${BCFTOOLS_LIFTOVER_ALT_STYLE:=ref-only}"
: "${BCFTOOLS_LIFTOVER_EXTRA_OPTS:=}"

for f in "$windows_bed" "$chain" "$src_fasta" "$target_fasta"; do
    if [ ! -f "$f" ]; then
        echo "error: input file not found: $f" >&2
        exit 1
    fi
done

if [ ! -f "${src_fasta}.fai" ]; then
    echo "error: source FASTA index not found: ${src_fasta}.fai (run: $SAMTOOLS_BIN faidx $src_fasta)" >&2
    exit 1
fi

if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1 && [ ! -x "$BCFTOOLS_BIN" ]; then
    echo "error: bcftools binary not found (checked: $BCFTOOLS_BIN)" >&2
    exit 1
fi

if ! "$BCFTOOLS_BIN" plugin -l 2>/dev/null | grep -qx "liftover"; then
    cat >&2 <<'MSG'
error: the bcftools 'liftover' plugin is not available.

It does NOT ship with bcftools — it lives in Giulio Genovese's separate
score plugin collection and must be compiled against the bcftools source:

    git clone https://github.com/freeseek/score
    cp score/liftover.c score/score.h /path/to/bcftools-<version>/plugins/
    cd /path/to/bcftools-<version> && make plugins/liftover.so
    export BCFTOOLS_PLUGINS=/path/to/bcftools-<version>/plugins

Then re-run. Reference: Genovese et al., Bioinformatics 40(2), 2024.
MSG
    exit 1
fi

mkdir -p "$(dirname "$out_prefix")"
n_windows=$(grep -cve '^[[:space:]]*$' "$windows_bed")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- 1. windows -> anchor VCF -------------------------------------------------
start_ts=$(date +%s.%N)
python3 "$script_dir/make_anchor_vcf.py" \
    --windows "$windows_bed" \
    --src-fasta "$src_fasta" \
    --out-vcf "${out_prefix}.anchors.unsorted.vcf" \
    --samtools "$SAMTOOLS_BIN" \
    --alt-style "$BCFTOOLS_LIFTOVER_ALT_STYLE"
"$BCFTOOLS_BIN" sort -Ov -o "${out_prefix}.anchors.vcf" "${out_prefix}.anchors.unsorted.vcf" 2>/dev/null
rm -f "${out_prefix}.anchors.unsorted.vcf"
end_ts=$(date +%s.%N)
anchor_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

# --- 2. lift ------------------------------------------------------------------
# shellcheck disable=SC2086
start_ts=$(date +%s.%N)
"$BCFTOOLS_BIN" +liftover -Ou "${out_prefix}.anchors.vcf" -- \
    -s "$src_fasta" \
    -f "$target_fasta" \
    -c "$chain" \
    --reject "${out_prefix}.reject.vcf" \
    -O v \
    --write-src \
    --write-reject \
    $BCFTOOLS_LIFTOVER_EXTRA_OPTS \
  | "$BCFTOOLS_BIN" sort -Ov -o "${out_prefix}.lifted.vcf" 2>/dev/null
end_ts=$(date +%s.%N)
lift_elapsed=$(awk -v a="$start_ts" -v b="$end_ts" 'BEGIN { printf "%.3f", b - a }')

# --- 3. pair anchors back into intervals --------------------------------------
python3 "$script_dir/bcftools_lifted_to_bed.py" \
    --lifted-vcf "${out_prefix}.lifted.vcf" \
    --reject-vcf "${out_prefix}.reject.vcf" \
    --windows "$windows_bed" \
    --out-lifted-bed "${out_prefix}.lifted.bed" \
    --out-unmapped-bed "${out_prefix}.unmapped.bed"

# --- 4. extract + tabulate (shared with the liftOver pipeline) ----------------
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

total_elapsed=$(awk -v a="$anchor_elapsed" -v b="$lift_elapsed" -v c="$extract_elapsed" 'BEGIN { printf "%.3f", a + b + c }')
printf "hap\tn_windows\tanchor_seconds\tliftover_seconds\textract_seconds\telapsed_seconds\n%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$hap_tag" "$n_windows" "$anchor_elapsed" "$lift_elapsed" "$extract_elapsed" "$total_elapsed" \
    >"${out_prefix}.timing.tsv"

echo "bcftools +liftover ($hap_tag): $n_windows windows in ${total_elapsed}s -> ${out_prefix}.results.tsv" >&2
