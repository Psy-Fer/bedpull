#!/usr/bin/env python3
"""Turn a bedpull --paf mode FASTA + --unmapped output into a results.tsv.

PAF-mode headers look like:
  contig_name|chr:region_start-region_end|region_name|contig_name:query_start-query_end|strand
  contig_name|chr:region_start-region_end|region_name|contig_name:query_start-query_end|strand|h1

(--partial, and its missing_left/missing_right header suffix, is BAM/CRAM-only
and never appears in PAF-mode output, so this parser doesn't need to handle it.)

A window (region_name) can produce zero, one, or several FASTA records — zero
if --unmapped logged it, several if more than one alignment overlapped the
window. This script keeps every hit; the scorer picks the best one per name
across haplotype runs.
"""

import argparse
import sys


def parse_fasta_lengths(fasta_path):
    """Yield (header, sequence_length) for each record in a FASTA file."""
    header = None
    seq_len = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_len
                header = line[1:]
                seq_len = 0
            else:
                seq_len += len(line)
    if header is not None:
        yield header, seq_len


def parse_unmapped_bed(unmapped_path):
    """Yield (reason, chrom, start, end, name) from a bedpull --unmapped file."""
    reason = None
    try:
        f = open(unmapped_path)
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
            yield reason, fields[0], int(fields[1]), int(fields[2]), fields[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True, help="bedpull FASTA output")
    p.add_argument("--unmapped", required=True, help="bedpull --unmapped output (may not exist)")
    p.add_argument("--hap", required=True, help="haplotype tag to stamp on every row, e.g. PATERNAL")
    p.add_argument("--out", required=True, help="output results.tsv path")
    args = p.parse_args()

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
    for header, seq_len in parse_fasta_lengths(args.fasta):
        fields = header.split("|")
        if len(fields) < 5:
            print(f"warning: unrecognised header shape, skipping: {header}", file=sys.stderr)
            continue
        contig_name = fields[0]
        region_name = fields[2]
        query_coords = fields[3]
        strand = fields[4][:1]
        query_start, query_end = "", ""
        if ":" in query_coords:
            _, coord_part = query_coords.split(":", 1)
            if "-" in coord_part:
                query_start, query_end = coord_part.split("-", 1)
        rows.append(
            {
                "name": region_name,
                "hap": args.hap,
                "status": "ok",
                "unmapped_reason": "",
                "extracted_len": seq_len,
                "contig_name": contig_name,
                "query_start": query_start,
                "query_end": query_end,
                "strand": strand,
            }
        )

    for reason, _chrom, _start, _end, name in parse_unmapped_bed(args.unmapped):
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

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_unmapped = sum(1 for r in rows if r["status"] == "unmapped")
    print(f"{args.hap}: {n_ok} hit(s), {n_unmapped} unmapped window(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
