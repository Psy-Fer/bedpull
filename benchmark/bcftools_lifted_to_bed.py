#!/usr/bin/env python3
"""Convert bcftools +liftover's output back into liftOver-shaped BED files.

Pairs the two anchor records emitted per window by make_anchor_vcf.py back
into a single lifted interval, and writes exactly the two files
liftOver itself produces:

  * a BED6 of successfully lifted intervals, and
  * a '#reason' + BED unmapped report.

Emitting liftOver's own output shape (rather than a results.tsv directly) is
deliberate: it lets this pipeline hand off to the existing
extract_and_score_liftover.py for the extract-and-tabulate half, so bcftools
+liftover is scored through *identical* extraction code to liftOver, pslMap
and liftOverMerge — same samtools batching, same minus-strand handling, same
per-region length keying. Any difference in the final numbers is then a
difference in lifting, not in how we measured it.

A window is reported as lifted only when *both* of its anchors lifted and
landed on the same target contig. The interval is [min(pos)-1, max(pos)) in
BED coordinates; anchors arriving in reversed order mean the chain flipped
the window, which is recorded as strand '-' so the extractor
reverse-complements it (matching what liftOver reports in BED6 column 6).

Reject reasons come from the FILTER column of +liftover's --reject output
(the plugin sets MissingContig, UnmappedAnchors, UnmappedAnchor5,
UnmappedAnchor3, MismatchAnchors, ApartAnchors or MissingFasta there).
Two reasons are ours rather than the plugin's, and are prefixed to keep that
distinction visible in the report:

  * local:DifferentContigs — both anchors lifted, but to different contigs.
  * local:NoAnchorRecord   — an anchor is in neither the lifted nor the
                             reject output (e.g. dropped by make_anchor_vcf.py
                             for having no valid reference base).
"""

import argparse
import sys
from collections import defaultdict

ID_SEP = "#"


def read_windows(path):
    """Yield (name, original_line) from the BED4 windows file, in file order."""
    with open(path) as f:
        for line in f:
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            fields = stripped.split("\t")
            if len(fields) < 4:
                continue
            yield fields[3], stripped


def split_anchor_id(anchor_id):
    """Split 'NAME#L' into ('NAME', 'L'); return (None, None) if malformed."""
    if ID_SEP not in anchor_id:
        return None, None
    name, _, tag = anchor_id.rpartition(ID_SEP)
    if tag not in ("L", "R"):
        return None, None
    return name, tag


def parse_vcf(path):
    """Yield (chrom, pos, id, filter) for each record in a (plain) VCF."""
    try:
        f = open(path)
    except FileNotFoundError:
        return
    with f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            yield fields[0], int(fields[1]), fields[2], fields[6]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lifted-vcf", required=True, help="bcftools +liftover output VCF (sorted)")
    p.add_argument("--reject-vcf", required=True, help="bcftools +liftover --reject VCF (may not exist)")
    p.add_argument("--windows", required=True, help="the original windows.bed handed to make_anchor_vcf.py")
    p.add_argument("--out-lifted-bed", required=True, help="output BED6 of lifted intervals")
    p.add_argument("--out-unmapped-bed", required=True, help="output liftOver-style unmapped report")
    args = p.parse_args()

    anchors = defaultdict(dict)
    n_dupe = 0
    for chrom, pos, anchor_id, _filt in parse_vcf(args.lifted_vcf):
        name, tag = split_anchor_id(anchor_id)
        if name is None:
            continue
        if tag in anchors[name]:
            # +liftover emits one record per input record, so a repeat means
            # the same anchor lifted more than once. Keep the first and count.
            n_dupe += 1
            continue
        anchors[name][tag] = (chrom, pos)

    rejects = defaultdict(dict)
    for _chrom, _pos, anchor_id, filt in parse_vcf(args.reject_vcf):
        name, tag = split_anchor_id(anchor_id)
        if name is None:
            continue
        rejects[name][tag] = filt or "."

    n_lifted = 0
    n_unmapped = 0
    reason_counts = defaultdict(int)

    with open(args.out_lifted_bed, "w") as lifted_out, open(args.out_unmapped_bed, "w") as unmapped_out:
        for name, original_line in read_windows(args.windows):
            got = anchors.get(name, {})
            left, right = got.get("L"), got.get("R")

            if left and right and left[0] == right[0]:
                chrom = left[0]
                strand = "+" if left[1] <= right[1] else "-"
                start = min(left[1], right[1]) - 1
                end = max(left[1], right[1])
                lifted_out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
                n_lifted += 1
                continue

            if left and right:
                reason = "local:DifferentContigs"
            else:
                missing = [tag for tag in ("L", "R") if tag not in got]
                reasons = []
                for tag in missing:
                    reasons.append(rejects.get(name, {}).get(tag) or "local:NoAnchorRecord")
                # Distinct reasons, order-stable, so a window failing both
                # anchors for the same reason reports it once.
                reason = ",".join(dict.fromkeys(reasons))

            reason_counts[reason] += 1
            n_unmapped += 1
            unmapped_out.write(f"#{reason}\n{original_line}\n")

    if n_dupe:
        print(f"warning: {n_dupe} duplicate anchor record(s) ignored", file=sys.stderr)
    summary = ", ".join(f"{r}={c}" for r, c in sorted(reason_counts.items(), key=lambda kv: -kv[1]))
    print(
        f"bcftools +liftover: {n_lifted} windows lifted, {n_unmapped} unmapped"
        + (f" ({summary})" if summary else ""),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
