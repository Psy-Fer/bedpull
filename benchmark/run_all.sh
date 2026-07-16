#!/usr/bin/env bash
# Orchestrate the hs1<->HG002 benchmark leg (the main event — see README.md
# for the hs1<->hg38 confirmation leg, which reuses these same scripts with
# different config values rather than a separate orchestrator).
#
# Usage:
#   run_all.sh <config.env> <stage>
#   stages: windows | chains | bedpull | liftover | negatives | score | scale | all

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $(basename "$0") <config.env> <stage>" >&2
    echo "  stages: windows | chains | bedpull | liftover | negatives | score | scale | all" >&2
    exit 1
fi

config="$1"
stage="$2"

if [ ! -f "$config" ]; then
    echo "error: config file not found: $config (copy config.example.env and fill it in)" >&2
    exit 1
fi
# shellcheck disable=SC1090
source "$config"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$BENCH_OUTDIR"

run_windows() {
    python3 "$script_dir/select_sv_windows.py" \
        --vcf "$GIAB_HS1_STVAR_VCF" \
        --high-conf-bed "$GIAB_HS1_STVAR_BED" \
        --fai "${HS1_FASTA}.fai" \
        --flank "$SV_FLANK_BP" \
        --min-len "$SV_MIN_LEN" \
        --max-windows "$SV_MAX_WINDOWS" \
        --out-windows "$BENCH_OUTDIR/windows.bed" \
        --out-truth "$BENCH_OUTDIR/truth.tsv" \
        --bcftools "$BCFTOOLS_BIN"
}

run_chains() {
    bash "$script_dir/build_chain_from_paf.sh" "$HG002_PATERNAL_TO_HS1_PAF" "$HG002_PATERNAL_TO_HS1_CHAIN"
    bash "$script_dir/build_chain_from_paf.sh" "$HG002_MATERNAL_TO_HS1_PAF" "$HG002_MATERNAL_TO_HS1_CHAIN"
}

run_bedpull() {
    for hap in PATERNAL MATERNAL; do
        paf_var="HG002_${hap}_TO_HS1_PAF"
        BEDPULL_BIN="$BEDPULL_BIN" bash "$script_dir/run_bedpull_pipeline.sh" \
            "$BENCH_OUTDIR/windows.bed" "${!paf_var}" "$HG002_ASSEMBLY_FASTA" "$hap" \
            "$BENCH_OUTDIR/bedpull_${hap}"
    done
}

run_liftover() {
    for hap in PATERNAL MATERNAL; do
        chain_var="HG002_${hap}_TO_HS1_CHAIN"
        LIFTOVER_BIN="$LIFTOVER_BIN" SAMTOOLS_BIN="$SAMTOOLS_BIN" bash "$script_dir/run_liftover_pipeline.sh" \
            "$BENCH_OUTDIR/windows.bed" "${!chain_var}" "$HG002_ASSEMBLY_FASTA" "$hap" \
            "$BENCH_OUTDIR/liftover_${hap}"
    done
}

run_negatives() {
    python3 "$script_dir/select_negative_control_windows.py" \
        --vcf "$GIAB_HS1_STVAR_VCF" \
        --high-conf-bed "$GIAB_HS1_STVAR_BED" \
        --truth "$BENCH_OUTDIR/truth.tsv" \
        --fai "${HS1_FASTA}.fai" \
        --flank "$SV_FLANK_BP" \
        --out-windows "$BENCH_OUTDIR/neg_windows.bed" \
        --out-truth "$BENCH_OUTDIR/neg_truth.tsv" \
        --bcftools "$BCFTOOLS_BIN"
    for hap in PATERNAL MATERNAL; do
        paf_var="HG002_${hap}_TO_HS1_PAF"
        BEDPULL_BIN="$BEDPULL_BIN" bash "$script_dir/run_bedpull_pipeline.sh" \
            "$BENCH_OUTDIR/neg_windows.bed" "${!paf_var}" "$HG002_ASSEMBLY_FASTA" "$hap" \
            "$BENCH_OUTDIR/neg_bedpull_${hap}"
        chain_var="HG002_${hap}_TO_HS1_CHAIN"
        LIFTOVER_BIN="$LIFTOVER_BIN" SAMTOOLS_BIN="$SAMTOOLS_BIN" bash "$script_dir/run_liftover_pipeline.sh" \
            "$BENCH_OUTDIR/neg_windows.bed" "${!chain_var}" "$HG002_ASSEMBLY_FASTA" "$hap" \
            "$BENCH_OUTDIR/neg_liftover_${hap}"
    done
}

run_score() {
    python3 "$script_dir/score_accuracy.py" \
        --truth "$BENCH_OUTDIR/truth.tsv" \
        --bedpull-results "$BENCH_OUTDIR/bedpull_PATERNAL.results.tsv" "$BENCH_OUTDIR/bedpull_MATERNAL.results.tsv" \
        --liftover-results "$BENCH_OUTDIR/liftover_PATERNAL.results.tsv" "$BENCH_OUTDIR/liftover_MATERNAL.results.tsv" \
        --neg-truth "$BENCH_OUTDIR/neg_truth.tsv" \
        --neg-bedpull-results "$BENCH_OUTDIR/neg_bedpull_PATERNAL.results.tsv" "$BENCH_OUTDIR/neg_bedpull_MATERNAL.results.tsv" \
        --neg-liftover-results "$BENCH_OUTDIR/neg_liftover_PATERNAL.results.tsv" "$BENCH_OUTDIR/neg_liftover_MATERNAL.results.tsv" \
        --roc-out "$BENCH_OUTDIR/roc_curve.png" \
        --roc-data-out "$BENCH_OUTDIR/roc_data.tsv" \
        --out-report "$BENCH_OUTDIR/accuracy_report.md" \
        --out-details "$BENCH_OUTDIR/accuracy_details.tsv"
    echo "report: $BENCH_OUTDIR/accuracy_report.md" >&2
}

run_scale() {
    mkdir -p "$BENCH_OUTDIR/scale"
    python3 "$script_dir/generate_random_windows.py" \
        --fai "${HS1_FASTA}.fai" \
        --sizes 100,1000,10000,100000 \
        --window-size 200 \
        --out-prefix "$BENCH_OUTDIR/scale/synthetic"
    BEDPULL_BIN="$BEDPULL_BIN" LIFTOVER_BIN="$LIFTOVER_BIN" SAMTOOLS_BIN="$SAMTOOLS_BIN" \
        bash "$script_dir/run_scaling_benchmark.sh" \
        "$BENCH_OUTDIR/scale" "$HG002_PATERNAL_TO_HS1_PAF" "$HG002_ASSEMBLY_FASTA" \
        "$BENCH_OUTDIR/scale_results.tsv" \
        "$HG002_PATERNAL_TO_HS1_CHAIN" "$HG002_ASSEMBLY_FASTA"
}

case "$stage" in
windows) run_windows ;;
chains) run_chains ;;
bedpull) run_bedpull ;;
liftover) run_liftover ;;
negatives) run_negatives ;;
score) run_score ;;
scale) run_scale ;;
all)
    run_windows
    run_chains
    run_bedpull
    run_liftover
    run_negatives
    run_score
    run_scale
    ;;
*)
    echo "error: unknown stage '$stage'" >&2
    echo "  stages: windows | chains | bedpull | liftover | negatives | score | scale | all" >&2
    exit 1
    ;;
esac
