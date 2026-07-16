#!/usr/bin/env python3
"""Select 'no true SV here' negative-control windows, matched in size to a
positive truth.tsv (from select_sv_windows.py), for computing precision/F1/
ROC alongside the recall-only accuracy report.

Without negative controls, "precision" is meaningless here: every window
select_sv_windows.py produces is a window we already know has a real SV, so
a tool can never be observed producing a spurious result — there's nothing
for it to be spurious *about*. This script generates windows we're
confident do NOT contain a real SV (far from every SVTYPE-classified VCF
record, inside the high-confidence regions), sized to match the real SVs'
REF-span-implied window sizes one-for-one, so extraction difficulty is
comparable between the positive and negative sets.

Expected behavior for a negative control is "extract exactly window_len" —
same expected-length formula as the positive set, just with alt_len==ref_len
(no length delta).
"""

import argparse
import bisect
import csv
import random
import subprocess
import sys
from collections import defaultdict


def load_fai_lengths(fai_path):
    chrom_lengths = []
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom_lengths.append((fields[0], int(fields[1])))
    return chrom_lengths


def load_high_conf_regions(bed_path):
    regions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            regions[fields[0]].append((int(fields[1]), int(fields[2])))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def fully_contained(regions_for_chrom, start, end):
    for r_start, r_end in regions_for_chrom:
        if r_start <= start and end <= r_end:
            return True
        if r_start > start:
            break
    return False


def merge_intervals(raw):
    """{chrom: [(start, end), ...]} -> same shape, sorted and merged."""
    merged = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        out = []
        for start, end in intervals:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = out
    return merged


def load_exclusion_intervals_from_vcf(vcf_path, bcftools_bin, pad):
    """All SVTYPE-classified variant positions, padded by `pad` and merged
    per chromosome, as a sorted list of (start, end) tuples for bisecting.
    """
    raw = defaultdict(list)
    cmd = [bcftools_bin, "query", "-f", "%CHROM\t%POS\t%REF\n", "-i", 'INFO/SVTYPE!="."', vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3:
            continue
        chrom, pos, ref = fields[0], int(fields[1]), fields[2]
        start = pos - 1 - pad
        end = pos - 1 + len(ref) + pad
        raw[chrom].append((max(0, start), end))
    if proc.wait() != 0:
        sys.exit("bcftools query failed while loading exclusion intervals")
    return merge_intervals(raw)


def load_exclusion_intervals_from_bed(bed_path, pad):
    """A plain 3-column BED of event spans (e.g. call_alignment_svs.py's
    exclusion_events.bed), padded by `pad` and merged per chromosome.
    """
    raw = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            raw[chrom].append((max(0, start - pad), end + pad))
    return merge_intervals(raw)


def too_close_to_any_sv(exclusions_for_chrom, start, end):
    """True if [start, end) overlaps any padded exclusion interval."""
    if not exclusions_for_chrom:
        return False
    starts = [iv[0] for iv in exclusions_for_chrom]
    i = bisect.bisect_right(starts, end)
    # Check the interval immediately before i (could still overlap) and i itself.
    for j in (i - 1, i):
        if 0 <= j < len(exclusions_for_chrom):
            iv_start, iv_end = exclusions_for_chrom[j]
            if iv_start < end and start < iv_end:
                return True
    return False


def load_truth_svlens(truth_path):
    svlens = []
    with open(truth_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            svlens.append(int(row["svlen"]))
    return svlens


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", help="GIAB stvar VCF (for exclusion positions) — mutually exclusive with --exclusion-bed")
    p.add_argument("--exclusion-bed", help="Plain 3-col BED of event spans to exclude near — mutually exclusive with --vcf")
    p.add_argument("--high-conf-bed", required=True, help="Benchmark high-confidence (or trusted-alignment) regions BED")
    p.add_argument("--truth", required=True, help="Positive truth.tsv from select_sv_windows.py — sizes negatives to match")
    p.add_argument("--fai", required=True, help="Reference .fai (for chromosome lengths and bounds)")
    p.add_argument("--min-distance", type=int, default=10000, help="Minimum distance from any known SV (default 10000bp)")
    p.add_argument("--flank", type=int, default=500, help="Must match the --flank used for the positive windows")
    p.add_argument("--max-attempts-per-window", type=int, default=200)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out-windows", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--bcftools", default="bcftools")
    args = p.parse_args()

    if bool(args.vcf) == bool(args.exclusion_bed):
        sys.exit("error: exactly one of --vcf or --exclusion-bed must be given")

    rng = random.Random(args.seed)
    chrom_lengths = load_fai_lengths(args.fai)
    high_conf = load_high_conf_regions(args.high_conf_bed)
    if args.vcf:
        exclusions = load_exclusion_intervals_from_vcf(args.vcf, args.bcftools, args.min_distance)
    else:
        exclusions = load_exclusion_intervals_from_bed(args.exclusion_bed, args.min_distance)
    svlens = load_truth_svlens(args.truth)

    kept = []
    n_failed = 0
    for i, svlen in enumerate(svlens):
        placed = False
        for _ in range(args.max_attempts_per_window):
            chrom, chrom_len = rng.choices(chrom_lengths, weights=[l for _, l in chrom_lengths])[0]
            if svlen + 2 * args.min_distance >= chrom_len:
                continue
            bed_start = rng.randint(0, chrom_len - svlen)
            bed_end = bed_start + svlen

            if too_close_to_any_sv(exclusions.get(chrom, []), bed_start, bed_end):
                continue
            if not fully_contained(high_conf.get(chrom, []), bed_start, bed_end):
                continue

            window_start = max(0, bed_start - args.flank)
            window_end = min(bed_end + args.flank, chrom_len)
            name = f"NEG_{chrom}_{bed_start}_{svlen}"
            kept.append(
                {
                    "name": name,
                    "chrom": chrom,
                    "svtype": "NONE",
                    "svlen": 0,
                    "matched_svlen": svlen,
                    "ref_len": svlen,
                    "alt_len": svlen,
                    "bed_start": bed_start,
                    "bed_end": bed_end,
                    "window_start": window_start,
                    "window_end": window_end,
                }
            )
            placed = True
            break
        if not placed:
            n_failed += 1

    kept.sort(key=lambda r: (r["chrom"], r["window_start"]))

    with open(args.out_windows, "w") as f:
        for r in kept:
            f.write(f"{r['chrom']}\t{r['window_start']}\t{r['window_end']}\t{r['name']}\n")

    truth_cols = [
        "name",
        "chrom",
        "svtype",
        "svlen",
        "matched_svlen",
        "ref_len",
        "alt_len",
        "bed_start",
        "bed_end",
        "window_start",
        "window_end",
    ]
    with open(args.out_truth, "w") as f:
        f.write("\t".join(truth_cols) + "\n")
        for r in kept:
            f.write("\t".join(str(r[c]) for c in truth_cols) + "\n")

    print(
        f"placed {len(kept)}/{len(svlens)} negative-control windows "
        f"({n_failed} could not find a valid spot within "
        f"{args.max_attempts_per_window} attempts)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
