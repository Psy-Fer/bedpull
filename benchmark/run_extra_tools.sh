#!/usr/bin/env bash
# Orchestrate the four "extra" tool comparisons (paftools.js, pslMap,
# liftOver -multiple, bedpull+stitch) for either leg, and produce the full
# 6-tool score_accuracy.py report — the exact commands used to generate the
# 6-tool tables in README.md, which previously existed only as ad hoc
# one-off commands run by hand during the investigation, not as a single
# reproducible script (a real gap flagged in an independent pre-release
# review — see TODO.md's 2026-07-16 entries).
#
# Usage:
#   run_extra_tools.sh <config.env> <leg> <stage>
#   leg:   main | confirmation
#   stage: pslmap | liftover_multiple | bedpull_stitch | paftools | score | all
#
# Prerequisites beyond run_all.sh's (see README.md's "Two more UCSC kent
# tools" / "Optional third tool" sections):
#   - BEDTOPSL_BIN, PSLMAP_BIN, LIFTOVERMERGE_BIN (UCSC kent tools)
#   - K8_BIN, PAFTOOLS_JS (k8 + minimap2's misc/paftools.js)
#   - main leg only: HS1_TO_HG002_PATERNAL_PAF / HS1_TO_HG002_MATERNAL_PAF
#     (reverse-direction PAFs, target=HG002 haplotype, query=hs1 — built
#     below if missing, but expensive: ~15-25 min each of whole-genome
#     asm5 alignment). confirmation leg's equivalent (HS1_TO_HG38_PAF's
#     reverse) must already exist as out/hg38/hg38_to_hs1.paf (see the
#     "Optional third tool" section for the exact build command) since it's
#     shared with that leg's own paftools.js writeup.
#
# The `main` leg's pslMap/liftover_multiple/bedpull_stitch stages are cheap
# (reuse the existing paf2chain-built chain and the existing PAT/MAT PAFs —
# seconds, not minutes). Only `paftools` needs a new alignment, and only for
# the main leg (the confirmation leg's reverse PAF is a one-time prerequisite
# already covered by that leg's own paftools.js section).

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "usage: $(basename "$0") <config.env> <leg: main|confirmation> <stage>" >&2
    echo "  stages: pslmap | liftover_multiple | bedpull_stitch | paftools | score | all" >&2
    exit 1
fi

config="$1"
leg="$2"
stage="$3"

if [ ! -f "$config" ]; then
    echo "error: config file not found: $config" >&2
    exit 1
fi
# shellcheck disable=SC1090
source "$config"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$leg" != "main" ] && [ "$leg" != "confirmation" ]; then
    echo "error: leg must be 'main' or 'confirmation', got: $leg" >&2
    exit 1
fi

: "${MINIMAP2_BIN:=minimap2}"
: "${SAMTOOLS_BIN:=samtools}"

if [ "$leg" = "main" ]; then
    outdir="$BENCH_OUTDIR"
    windows="$outdir/windows.bed"
    neg_windows="$outdir/neg_windows.bed"
    hap_tags=(PATERNAL MATERNAL)
    chain_var_prefix="HG002"
    target_fasta="$HG002_ASSEMBLY_FASTA"
    chrom_sizes="$outdir/hs1.chrom.sizes"
    windows_genome_fasta="$HS1_FASTA"
else
    outdir="$BENCH_OUTDIR/hg38"
    windows="$outdir/windows.bed"
    neg_windows="$outdir/neg_windows.bed"
    hap_tags=(hg38)
    target_fasta="$HS1_FASTA"
    chrom_sizes="$outdir/hg38.chrom.sizes"
    windows_genome_fasta="$HG38_FASTA"
fi
mkdir -p "$outdir"

run_pslmap_stage() {
    awk '{print $1"\t"$2}' "${windows_genome_fasta}.fai" > "$chrom_sizes"
    for hap in "${hap_tags[@]}"; do
        if [ "$leg" = "main" ]; then
            chain_var="${chain_var_prefix}_${hap}_TO_HS1_CHAIN"
            chain="${!chain_var}"
        else
            chain="$HG38_TO_HS1_CHAIN"
        fi
        bash "$script_dir/run_pslmap_pipeline.sh" "$windows" "$chain" "$chrom_sizes" "$target_fasta" "$hap" "$outdir/pslmap_${hap}"
        bash "$script_dir/run_pslmap_pipeline.sh" "$neg_windows" "$chain" "$chrom_sizes" "$target_fasta" "$hap" "$outdir/neg_pslmap_${hap}"
    done
}

run_liftover_multiple_stage() {
    for hap in "${hap_tags[@]}"; do
        if [ "$leg" = "main" ]; then
            chain_var="${chain_var_prefix}_${hap}_TO_HS1_CHAIN"
            chain="${!chain_var}"
        else
            chain="$HG38_TO_HS1_CHAIN"
        fi
        bash "$script_dir/run_liftover_multiple_pipeline.sh" "$windows" "$chain" "$target_fasta" "$hap" "$outdir/liftmulti_${hap}"
        bash "$script_dir/run_liftover_multiple_pipeline.sh" "$neg_windows" "$chain" "$target_fasta" "$hap" "$outdir/neg_liftmulti_${hap}"
    done
}

run_bedpull_stitch_stage() {
    for hap in "${hap_tags[@]}"; do
        if [ "$leg" = "main" ]; then
            paf_var="HG002_${hap}_TO_HS1_PAF"
            paf="${!paf_var}"
        else
            paf="$HS1_TO_HG38_PAF"
        fi
        bash "$script_dir/run_bedpull_pipeline.sh" "$windows" "$paf" "$target_fasta" "$hap" "$outdir/bedpullstitch_${hap}" --stitch_records
        bash "$script_dir/run_bedpull_pipeline.sh" "$neg_windows" "$paf" "$target_fasta" "$hap" "$outdir/neg_bedpullstitch_${hap}" --stitch_records
    done
}

run_paftools_stage() {
    for hap in "${hap_tags[@]}"; do
        if [ "$leg" = "main" ]; then
            reverse_paf="$outdir/hs1_to_hg002${hap,,}.paf"
            if [ ! -f "$reverse_paf" ]; then
                mmi="$outdir/hg002${hap,,}.mmi"
                fasta_var="HG002_ASSEMBLY_FASTA"
                haplotype_fasta="$outdir/hg002.${hap,,}.fa"
                if [ ! -f "$haplotype_fasta" ]; then
                    echo "extracting ${hap} contigs from ${!fasta_var}..." >&2
                    "$SAMTOOLS_BIN" faidx "${!fasta_var}" $(grep "_${hap}" "${!fasta_var}.fai" | cut -f1) > "$haplotype_fasta"
                fi
                echo "building reverse-direction alignment for paftools.js (target=${hap}, query=hs1) — this is a whole-genome asm5 alignment, expect ~15-25 minutes..." >&2
                "$MINIMAP2_BIN" -x asm5 -d "$mmi" "$haplotype_fasta"
                "$MINIMAP2_BIN" -cx asm5 -t 16 "$mmi" "$HS1_FASTA" > "$reverse_paf"
                rm -f "$mmi"
            fi
        else
            reverse_paf="$outdir/hg38_to_hs1.paf"
            if [ ! -f "$reverse_paf" ]; then
                echo "error: $reverse_paf not found — build it first (see README.md's" \
                     "'Optional third tool' section for the exact minimap2 command)" >&2
                exit 1
            fi
        fi
        bash "$script_dir/run_paftools_pipeline.sh" "$windows" "$reverse_paf" "$target_fasta" "$hap" "$outdir/paftools_${hap}" 0 0
        bash "$script_dir/run_paftools_pipeline.sh" "$neg_windows" "$reverse_paf" "$target_fasta" "$hap" "$outdir/neg_paftools_${hap}" 0 0
    done
}

run_score_stage() {
    bedpull_results=() liftover_results=() paftools_results=() pslmap_results=() liftmulti_results=() bedpullstitch_results=()
    neg_bedpull_results=() neg_liftover_results=() neg_paftools_results=() neg_pslmap_results=() neg_liftmulti_results=() neg_bedpullstitch_results=()
    for hap in "${hap_tags[@]}"; do
        bedpull_results+=("$outdir/bedpull_${hap}.results.tsv")
        liftover_results+=("$outdir/liftover_${hap}.results.tsv")
        paftools_results+=("$outdir/paftools_${hap}.results.tsv")
        pslmap_results+=("$outdir/pslmap_${hap}.results.tsv")
        liftmulti_results+=("$outdir/liftmulti_${hap}.results.tsv")
        bedpullstitch_results+=("$outdir/bedpullstitch_${hap}.results.tsv")
        neg_bedpull_results+=("$outdir/neg_bedpull_${hap}.results.tsv")
        neg_liftover_results+=("$outdir/neg_liftover_${hap}.results.tsv")
        neg_paftools_results+=("$outdir/neg_paftools_${hap}.results.tsv")
        neg_pslmap_results+=("$outdir/neg_pslmap_${hap}.results.tsv")
        neg_liftmulti_results+=("$outdir/neg_liftmulti_${hap}.results.tsv")
        neg_bedpullstitch_results+=("$outdir/neg_bedpullstitch_${hap}.results.tsv")
    done

    if [ "$leg" = "main" ]; then
        truth="$outdir/truth.tsv"; neg_truth="$outdir/neg_truth.tsv"
        report="$outdir/accuracy_report.md"; details="$outdir/accuracy_details.tsv"
        roc="$outdir/roc_curve.png"; roc_data="$outdir/roc_data.tsv"
    else
        truth="$outdir/truth.tsv"; neg_truth="$outdir/neg_truth.tsv"
        report="$outdir/confirmation_report.md"; details="$outdir/confirmation_details.tsv"
        roc="$outdir/roc_curve.png"; roc_data="$outdir/roc_data.tsv"
    fi

    python3 "$script_dir/score_accuracy.py" \
        --truth "$truth" \
        --bedpull-results "${bedpull_results[@]}" \
        --liftover-results "${liftover_results[@]}" \
        --paftools-results "${paftools_results[@]}" \
        --pslmap-results "${pslmap_results[@]}" \
        --liftover-multiple-results "${liftmulti_results[@]}" \
        --bedpull-stitch-results "${bedpullstitch_results[@]}" \
        --neg-truth "$neg_truth" \
        --neg-bedpull-results "${neg_bedpull_results[@]}" \
        --neg-liftover-results "${neg_liftover_results[@]}" \
        --neg-paftools-results "${neg_paftools_results[@]}" \
        --neg-pslmap-results "${neg_pslmap_results[@]}" \
        --neg-liftover-multiple-results "${neg_liftmulti_results[@]}" \
        --neg-bedpull-stitch-results "${neg_bedpullstitch_results[@]}" \
        --roc-out "$roc" --roc-data-out "$roc_data" \
        --out-report "$report" --out-details "$details"
    echo "full 6-tool report: $report" >&2
}

case "$stage" in
pslmap) run_pslmap_stage ;;
liftover_multiple) run_liftover_multiple_stage ;;
bedpull_stitch) run_bedpull_stitch_stage ;;
paftools) run_paftools_stage ;;
score) run_score_stage ;;
all)
    run_pslmap_stage
    run_liftover_multiple_stage
    run_bedpull_stitch_stage
    run_paftools_stage
    run_score_stage
    ;;
*)
    echo "error: unknown stage '$stage'" >&2
    echo "  stages: pslmap | liftover_multiple | bedpull_stitch | paftools | score | all" >&2
    exit 1
    ;;
esac
