#!/usr/bin/env python3
"""Select SV windows from a GIAB structural-variant VCF for the bedpull benchmark.

Reads a sequence-resolved SV VCF (SVTYPE/SVLEN populated only on the actual SV
records — GIAB's v5.0q "stvar" files also contain small variants with these
fields unset), restricts to records fully contained in the high-confidence BED,
and emits:

  - windows.bed: chrom, window_start, window_end, name  (0-based half-open,
    the SV's REF span expanded by --flank on each side)
  - truth.tsv: one row per window with the fields needed to score extraction
    accuracy later (svtype, svlen, genotype, window bounds, etc.)

By design this does NOT try to map VCF genotype phase (0|1 vs 1|0) to a
specific assembly haplotype (_PATERNAL vs _MATERNAL) — that mapping isn't
guaranteed stable across pipelines. Downstream, run_bedpull_pipeline.sh and
run_liftover_pipeline.sh try both haplotypes for every window and the scorer
takes the best match. Genotype is still recorded so ambiguous/no-call sites
can be filtered out (the default) or included for inspection.
"""

import argparse
import os
import random
import subprocess
import sys
from collections import defaultdict

GENOTYPE_SEP_CHARS = "|/"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="GIAB stvar VCF (bgzipped, tabix-indexed)")
    p.add_argument("--high-conf-bed", required=True, help="Benchmark high-confidence regions BED")
    p.add_argument("--fai", help="Reference .fai, used to clip windows to chromosome bounds")
    p.add_argument("--flank", type=int, default=500, help="Bases of context on each side of the SV")
    p.add_argument("--min-len", type=int, default=50, help="Minimum |SVLEN| to include")
    p.add_argument("--max-len", type=int, default=0, help="Maximum |SVLEN| to include (0 = no cap)")
    p.add_argument(
        "--svtype",
        default="INS,DEL",
        help="Comma-separated SVTYPE values to include (default: INS,DEL)",
    )
    p.add_argument(
        "--include-ambiguous-gt",
        action="store_true",
        help="Also emit windows for no-call/hemizygous genotypes (excluded by default)",
    )
    p.add_argument("--max-windows", type=int, default=0, help="Cap on emitted windows (0 = no cap)")
    p.add_argument("--seed", type=int, default=1, help="RNG seed used when --max-windows subsamples")
    p.add_argument("--out-windows", required=True, help="Output windows.bed path")
    p.add_argument("--out-truth", required=True, help="Output truth.tsv path")
    p.add_argument("--bcftools", default="bcftools", help="bcftools binary (default: on PATH)")
    return p.parse_args()


def load_high_conf_regions(bed_path):
    """Load the high-confidence BED into {chrom: [(start, end), ...]}, sorted."""
    regions = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            regions[chrom].append((start, end))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def load_fai_lengths(fai_path):
    lengths = {}
    if not fai_path:
        return lengths
    if not os.path.isfile(fai_path):
        print(f"warning: --fai path not found, skipping chromosome-bound clipping: {fai_path}", file=sys.stderr)
        return lengths
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def fully_contained(regions_for_chrom, start, end):
    """True if [start, end) sits entirely inside one high-confidence interval.

    regions_for_chrom is sorted by start, so a linear scan is fine for a
    one-off containment check; callers only invoke this once per SV.
    """
    for r_start, r_end in regions_for_chrom:
        if r_start <= start and end <= r_end:
            return True
        if r_start > start:
            break
    return False


def classify_genotype(gt):
    """Return (allele1, allele2, genotype_class) for a bcftools GT string.

    genotype_class is one of: het, hom, no_call, hemizygous.
    """
    alleles = None
    for sep in GENOTYPE_SEP_CHARS:
        if sep in gt:
            alleles = gt.split(sep)
            break
    if alleles is None:
        # No separator: a single value, e.g. chrX/chrY hemizygous regions.
        return gt, None, "hemizygous"

    a1, a2 = alleles[0], alleles[1]
    if a1 == "." or a2 == ".":
        return a1, a2, "no_call"
    if a1 == a2:
        return a1, a2, "hom"
    return a1, a2, "het"


def query_vcf_records(vcf_path, bcftools_bin):
    """Yield (chrom, pos, ref, alt, svtype, svlen, gt) for SV-classified records."""
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\n"
    cmd = [bcftools_bin, "query", "-f", fmt, "-i", 'INFO/SVTYPE!="."', vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 7:
            continue
        chrom, pos, ref, alt, svtype, svlen, gt = fields
        yield chrom, int(pos), ref, alt, svtype, int(svlen), gt
    ret = proc.wait()
    if ret != 0:
        sys.exit(f"bcftools query failed with exit code {ret}")


def main():
    args = parse_args()
    wanted_svtypes = set(args.svtype.split(","))
    high_conf = load_high_conf_regions(args.high_conf_bed)
    chrom_lengths = load_fai_lengths(args.fai)

    kept = []
    seen_names = set()
    n_seen = 0
    n_wrong_type_or_len = 0
    n_ambiguous_gt = 0
    n_not_contained = 0

    for chrom, pos, ref, alt, svtype, svlen, gt in query_vcf_records(args.vcf, args.bcftools):
        n_seen += 1
        abs_len = abs(svlen)
        if svtype not in wanted_svtypes or abs_len < args.min_len:
            n_wrong_type_or_len += 1
            continue
        if args.max_len > 0 and abs_len > args.max_len:
            n_wrong_type_or_len += 1
            continue

        a1, a2, genotype_class = classify_genotype(gt)
        if genotype_class in ("no_call", "hemizygous") and not args.include_ambiguous_gt:
            n_ambiguous_gt += 1
            continue

        # 0-based half-open span of the REF allele (covers the full deleted
        # region for a DEL; for an INS this is just the anchor base(s)).
        bed_start = pos - 1
        bed_end = bed_start + len(ref)

        if not fully_contained(high_conf.get(chrom, []), bed_start, bed_end):
            n_not_contained += 1
            continue

        window_start = max(0, bed_start - args.flank)
        window_end = bed_end + args.flank
        chrom_len = chrom_lengths.get(chrom)
        if chrom_len is not None:
            window_end = min(window_end, chrom_len)

        name = f"{svtype}_{chrom}_{pos}_{abs_len}"
        if name in seen_names:
            name = f"{name}_{n_seen}"
        seen_names.add(name)

        kept.append(
            {
                "name": name,
                "chrom": chrom,
                "svtype": svtype,
                "svlen": abs_len,
                "ref_len": len(ref),
                "alt_len": len(alt),
                "genotype": gt,
                "genotype_class": genotype_class,
                "bed_start": bed_start,
                "bed_end": bed_end,
                "window_start": window_start,
                "window_end": window_end,
            }
        )

    if args.max_windows > 0 and len(kept) > args.max_windows:
        rng = random.Random(args.seed)
        kept = rng.sample(kept, args.max_windows)

    kept.sort(key=lambda r: (r["chrom"], r["window_start"]))

    with open(args.out_windows, "w") as bed_out:
        for r in kept:
            bed_out.write(f"{r['chrom']}\t{r['window_start']}\t{r['window_end']}\t{r['name']}\n")

    truth_cols = [
        "name",
        "chrom",
        "svtype",
        "svlen",
        "ref_len",
        "alt_len",
        "genotype",
        "genotype_class",
        "bed_start",
        "bed_end",
        "window_start",
        "window_end",
    ]
    with open(args.out_truth, "w") as truth_out:
        truth_out.write("\t".join(truth_cols) + "\n")
        for r in kept:
            truth_out.write("\t".join(str(r[c]) for c in truth_cols) + "\n")

    print(
        f"scanned {n_seen} SV-classified records; kept {len(kept)}; "
        f"dropped {n_wrong_type_or_len} (type/length), {n_ambiguous_gt} (ambiguous GT), "
        f"{n_not_contained} (outside high-confidence regions)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
