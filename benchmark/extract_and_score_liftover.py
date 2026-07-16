#!/usr/bin/env python3
"""Extract sequence at liftOver's output coordinates and write a results.tsv
in the same shape as parse_paf_results.py, for direct comparison.

Batches the samtools faidx calls (one per strand, via -r region-file) rather
than one process per region — the same reopen-per-item cost bedpull itself
had to fix applies just as much here once you're extracting thousands of
windows.

liftOver's own unmapped file uses the same convention we borrowed for
bedpull's --unmapped: a '#reason' comment line followed by the original
input BED record.

A name can legitimately appear more than once in --lifted — e.g.
liftOverMerge can emit several disjoint merged fragments for one input
region (different chromosome, or too far apart to merge), and pslMap can
likewise produce multiple candidate mappings per window. Every occurrence
must be scored (score_accuracy.py's best_hit() already picks whichever is
closest to the expected length across all rows sharing a name — the same
mechanism that already handles multi-haplotype results), so lengths are
looked up per (chrom, start, end) region rather than collapsed into a
single per-name value.
"""

import argparse
import subprocess
import sys
from collections import defaultdict


def parse_bed6(path):
    """Yield (chrom, start, end, name, score, strand) from a BED6 file."""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            strand = fields[5] if len(fields) > 5 else "+"
            yield chrom, start, end, name, strand


def parse_unmapped(path):
    """Yield (reason, name) from a liftOver-style '#reason' + BED unmapped file."""
    reason = None
    try:
        f = open(path)
    except FileNotFoundError:
        return
    with f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                reason = line[1:]
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            yield reason, fields[3]


def batch_faidx(samtools_bin, target_fasta, regions, reverse_complement, workdir_prefix):
    """Run one samtools faidx -r call for all given regions (sharing a strand).

    Returns {region_string: sequence_length}. Keyed by region rather than
    name — a name can map to more than one region (see module docstring),
    but a given region string always has one unambiguous length.
    """
    if not regions:
        return {}

    region_file = f"{workdir_prefix}.regions.txt"
    with open(region_file, "w") as f:
        for region in regions:
            f.write(region + "\n")

    cmd = [samtools_bin, "faidx", "-r", region_file, "--mark-strand", "no"]
    if reverse_complement:
        cmd.append("-i")
    cmd.append(target_fasta)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"samtools faidx failed: {proc.stderr}")

    lengths = {}

    def flush(header, seq_len):
        if header is None:
            return
        lengths[header] = seq_len

    header = None
    seq_len = 0
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            flush(header, seq_len)
            header = line[1:].strip()
            seq_len = 0
        else:
            seq_len += len(line.strip())
    flush(header, seq_len)
    return lengths


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lifted", required=True, help="liftOver's lifted-coordinates BED6 output")
    p.add_argument("--unmapped", required=True, help="liftOver's unmapped output (may not exist)")
    p.add_argument("--target-fasta", required=True, help="Reference/assembly FASTA to extract from")
    p.add_argument("--hap", required=True, help="haplotype tag to stamp on every row")
    p.add_argument("--samtools", default="samtools", help="samtools binary (default: on PATH)")
    p.add_argument("--out", required=True, help="output results.tsv path")
    args = p.parse_args()

    plus_regions = set()
    minus_regions = set()
    lifted_rows = []
    for chrom, start, end, name, strand in parse_bed6(args.lifted):
        region = f"{chrom}:{start + 1}-{end}"
        lifted_rows.append((name, chrom, start, end, strand, region))
        if strand == "-":
            minus_regions.add(region)
        else:
            plus_regions.add(region)

    lengths = {}
    lengths.update(
        batch_faidx(args.samtools, args.target_fasta, plus_regions, False, f"{args.out}.plus")
    )
    lengths.update(
        batch_faidx(args.samtools, args.target_fasta, minus_regions, True, f"{args.out}.minus")
    )

    cols = [
        "name",
        "hap",
        "status",
        "unmapped_reason",
        "extracted_len",
        "contig_name",
        "query_start",
        "query_end",
        "strand",
    ]
    rows = []
    for name, chrom, start, end, strand, region in lifted_rows:
        rows.append(
            {
                "name": name,
                "hap": args.hap,
                "status": "ok",
                "unmapped_reason": "",
                "extracted_len": lengths.get(region, 0),
                "contig_name": chrom,
                "query_start": start,
                "query_end": end,
                "strand": strand,
            }
        )

    n_unmapped = 0
    for reason, name in parse_unmapped(args.unmapped):
        n_unmapped += 1
        rows.append(
            {
                "name": name,
                "hap": args.hap,
                "status": "unmapped",
                "unmapped_reason": reason or "",
                "extracted_len": 0,
                "contig_name": "",
                "query_start": "",
                "query_end": "",
                "strand": "",
            }
        )

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_lifted_names = len({r["name"] for r in rows if r["status"] == "ok"})
    print(f"{args.hap}: {n_lifted_names} lifted ({len(lifted_rows)} region(s)), {n_unmapped} unmapped", file=sys.stderr)


if __name__ == "__main__":
    main()
