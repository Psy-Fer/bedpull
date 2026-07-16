#!/usr/bin/env python3
"""Call structural-variant-like indel events directly from a whole-genome PAF
alignment's own CIGAR, for use as the confirmation leg's truth set.

Why this exists: the confirmation leg extracts sequence from the hs1
(CHM13v2.0) *reference* itself, not from any individual's assembly. The GIAB
HG002 SV truth set (used by the main leg) describes differences between
HG002 and a reference — those variants belong to HG002 the person, and there
is no guarantee CHM13 (an independent, different individual's genome) carries
the same variant at the same locus. Scoring "did we extract the HG002-sized
sequence from hs1" is therefore not a real accuracy measurement; it just
measures how often CHM13 happens to agree with HG002.

The fix: build truth from genuine hg38<->hs1 differences, i.e. indel events
present in the alignment between the two reference genomes actually being
extracted between. This makes the confirmation leg internally consistent —
both bedpull (PAF mode) and liftOver (official chain) are scored against
"does the extracted hs1 sequence have the length the alignment itself says it
should", not against a third genome's variant calls.

Reads a PAF with target=hg38, query=hs1 (matches the windows' hg38 coordinate
system and the direction the official hg38->hs1 chain lifts in). Walks each
qualifying record's cg:Z: CIGAR, merges consecutive I/D runs into single
events, and keeps ones far enough from record boundaries and big enough to
matter.

Emits three files:
  - windows.bed / truth.tsv: same schema as select_sv_windows.py's output,
    so score_accuracy.py needs no changes to consume it.
  - aligned_regions.bed: merged, edge-trimmed spans of every qualifying PAF
    record — the confirmation leg's analogue of the GIAB high-confidence BED,
    used both to keep windows away from alignment-block edges and as the
    containment region for select_negative_control_windows.py.
  - exclusion_events.bed: every indel event >= --min-exclusion-len (a lower
    bar than --min-len), for select_negative_control_windows.py to avoid when
    placing "no SV here" negative controls.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paf", required=True, help="PAF with target=hg38, query=hs1, cg:Z: CIGAR tag required")
    p.add_argument("--fai", help="Target (hg38) .fai, used to clip windows to chromosome bounds")
    p.add_argument("--flank", type=int, default=500, help="Bases of context on each side of the event")
    p.add_argument("--min-len", type=int, default=50, help="Minimum |svlen| to keep as positive truth")
    p.add_argument("--min-exclusion-len", type=int, default=10, help="Minimum event size to record for negative-control exclusion")
    p.add_argument("--min-mapq", type=int, default=60, help="Minimum PAF mapq to trust a record")
    p.add_argument("--edge-buffer", type=int, default=1000, help="Bases from a record's target start/end to stay clear of (alignment-boundary artifacts)")
    p.add_argument("--max-windows", type=int, default=0, help="Cap on emitted windows (0 = no cap)")
    p.add_argument("--seed", type=int, default=1, help="RNG seed used when --max-windows subsamples")
    p.add_argument("--out-windows", required=True)
    p.add_argument("--out-truth", required=True)
    p.add_argument("--out-aligned-bed", required=True)
    p.add_argument("--out-exclusion-bed", required=True)
    return p.parse_args()


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


def parse_paf_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def iter_qualifying_records(paf_path, min_mapq):
    n_seen = 0
    n_low_mapq = 0
    n_not_primary = 0
    n_no_cigar = 0
    with open(paf_path) as f:
        for line in f:
            n_seen += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            mapq = int(fields[11])
            tags = parse_paf_tags(fields)
            tp = tags.get("tp")
            if tp is not None and tp != "P":
                n_not_primary += 1
                continue
            if mapq < min_mapq:
                n_low_mapq += 1
                continue
            cigar = tags.get("cg")
            if not cigar:
                n_no_cigar += 1
                continue
            yield {
                "target_name": fields[5],
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "cigar": cigar,
            }
    print(
        f"PAF records: {n_seen} total; dropped {n_not_primary} (not primary), "
        f"{n_low_mapq} (mapq < {min_mapq}), {n_no_cigar} (no cg:Z: tag)",
        file=sys.stderr,
    )


def call_events(record):
    """Yield (bed_start, bed_end, ref_len, alt_len) for each merged I/D run
    in the record's CIGAR, in target coordinate space.
    """
    t_pos = record["target_start"]
    event_start = None
    cur_ref = 0
    cur_alt = 0

    def finalize():
        if event_start is None:
            return None
        return (event_start, event_start + cur_ref, cur_ref, cur_alt)

    for length_str, op in CIGAR_OP_RE.findall(record["cigar"]):
        length = int(length_str)
        if op in "M=X":
            ev = finalize()
            if ev is not None:
                yield ev
            event_start = None
            cur_ref = 0
            cur_alt = 0
            t_pos += length
        elif op == "D":
            if event_start is None:
                event_start = t_pos
            cur_ref += length
            t_pos += length
        elif op == "I":
            if event_start is None:
                event_start = t_pos
            cur_alt += length
        elif op == "N":
            ev = finalize()
            if ev is not None:
                yield ev
            event_start = None
            cur_ref = 0
            cur_alt = 0
            t_pos += length
        # S/H/P don't consume target and shouldn't appear mid-alignment in a
        # PAF cg:Z: string; ignore rather than let them merge into an event.

    ev = finalize()
    if ev is not None:
        yield ev


def main():
    args = parse_args()
    chrom_lengths = load_fai_lengths(args.fai)

    aligned_intervals = defaultdict(list)
    exclusion_events = []
    truth_candidates = []
    n_events_seen = 0
    n_edge_dropped = 0

    for record in iter_qualifying_records(args.paf, args.min_mapq):
        chrom = record["target_name"]
        lo = record["target_start"] + args.edge_buffer
        hi = record["target_end"] - args.edge_buffer
        if hi > lo:
            aligned_intervals[chrom].append((lo, hi))

        for bed_start, bed_end, ref_len, alt_len in call_events(record):
            n_events_seen += 1
            if bed_start < lo or bed_end > hi:
                n_edge_dropped += 1
                continue

            svlen = alt_len - ref_len
            abs_len = abs(svlen)
            if abs_len >= args.min_exclusion_len:
                exclusion_events.append((chrom, bed_start, bed_end))
            if abs_len >= args.min_len:
                truth_candidates.append(
                    {
                        "chrom": chrom,
                        "svtype": "INS" if svlen > 0 else "DEL",
                        "svlen": abs_len,
                        "ref_len": ref_len,
                        "alt_len": alt_len,
                        "bed_start": bed_start,
                        "bed_end": bed_end,
                    }
                )

    # Merge aligned (trusted) intervals per chromosome.
    merged_aligned = {}
    for chrom, intervals in aligned_intervals.items():
        intervals.sort()
        out = []
        for start, end in intervals:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged_aligned[chrom] = out

    # Drop truth candidates whose locus was seen more than once (e.g.
    # overlapping alignment blocks disagreeing on the same target span) —
    # keep the first (largest-mapq-sorted input order is not guaranteed, so
    # this is a conservative de-dup on exact position collisions only).
    seen_positions = set()
    kept = []
    n_duplicate = 0
    for c in sorted(truth_candidates, key=lambda r: (r["chrom"], r["bed_start"])):
        key = (c["chrom"], c["bed_start"], c["bed_end"])
        if key in seen_positions:
            n_duplicate += 1
            continue
        seen_positions.add(key)
        kept.append(c)

    for r in kept:
        chrom_len = chrom_lengths.get(r["chrom"])
        window_start = max(0, r["bed_start"] - args.flank)
        window_end = r["bed_end"] + args.flank
        if chrom_len is not None:
            window_end = min(window_end, chrom_len)
        r["window_start"] = window_start
        r["window_end"] = window_end
        r["name"] = f"{r['svtype']}_{r['chrom']}_{r['bed_start']}_{r['svlen']}"
        r["genotype"] = "1/1"
        r["genotype_class"] = "hom"

    if args.max_windows > 0 and len(kept) > args.max_windows:
        import random

        rng = random.Random(args.seed)
        kept = rng.sample(kept, args.max_windows)

    kept.sort(key=lambda r: (r["chrom"], r["window_start"]))

    with open(args.out_windows, "w") as f:
        for r in kept:
            f.write(f"{r['chrom']}\t{r['window_start']}\t{r['window_end']}\t{r['name']}\n")

    truth_cols = [
        "name", "chrom", "svtype", "svlen", "ref_len", "alt_len",
        "genotype", "genotype_class", "bed_start", "bed_end",
        "window_start", "window_end",
    ]
    with open(args.out_truth, "w") as f:
        f.write("\t".join(truth_cols) + "\n")
        for r in kept:
            f.write("\t".join(str(r[c]) for c in truth_cols) + "\n")

    with open(args.out_aligned_bed, "w") as f:
        for chrom in sorted(merged_aligned):
            for start, end in merged_aligned[chrom]:
                f.write(f"{chrom}\t{start}\t{end}\n")

    exclusion_events.sort()
    with open(args.out_exclusion_bed, "w") as f:
        for chrom, start, end in exclusion_events:
            f.write(f"{chrom}\t{start}\t{end}\n")

    print(
        f"events seen: {n_events_seen}; dropped {n_edge_dropped} (near alignment-block "
        f"edge), {n_duplicate} (duplicate locus); kept {len(kept)} truth windows "
        f"(>= {args.min_len}bp) and {len(exclusion_events)} exclusion events "
        f"(>= {args.min_exclusion_len}bp)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
