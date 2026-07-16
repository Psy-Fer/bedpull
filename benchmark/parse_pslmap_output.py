#!/usr/bin/env python3
"""Recover which of our windows pslMap successfully lifted end-to-end, and
turn its PSL output into the same lifted.bed/unmapped.bed shape
extract_and_score_liftover.py already consumes.

pslMap (given a whole-window "self-alignment" PSL from `bedToPsl -keepQuery`
projected through a chain) emits one output row per *overlapping chain
record*, not one row per input window — a window whose true extent isn't
covered by a single chain record produces several short fragment rows
instead of one whole-window row (or none, if it's badly enough split). We
only trust rows whose query span exactly matches one of our windows'
(chrom, start, end) — that's pslMap having found one chain record covering
the entire window in one going (matching how bedpull/liftOver/paftools.js
all try to produce a single whole-window answer). Partial-fragment rows are
silently discarded rather than guessed at; a window with no exact match is
reported unmapped, exactly like any other tool that couldn't do it in one
piece.

PSL columns (0-based): 0 matches, 1 misMatches, 2 repMatches, 3 nCount,
4 qNumInsert, 5 qBaseInsert, 6 tNumInsert, 7 tBaseInsert, 8 strand,
9 qName, 10 qSize, 11 qStart, 12 qEnd, 13 tName, 14 tSize, 15 tStart,
16 tEnd, 17 blockCount, 18 blockSizes, 19 qStarts, 20 tStarts.
"""

import argparse
import sys
from collections import defaultdict


def load_windows(path):
    by_key = defaultdict(list)
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            by_key[(chrom, start, end)].append(name)
            rows.append((chrom, start, end, name))
    return by_key, rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--windows", required=True, help="Original input windows.bed (chrom, start, end, name)")
    p.add_argument("--mapped-psl", required=True, help="pslMap's output PSL")
    p.add_argument("--out-lifted", required=True, help="BED6 output with real names recovered")
    p.add_argument("--out-unmapped", required=True, help="liftOver-style '#reason' + BED unmapped output")
    args = p.parse_args()

    by_key, rows = load_windows(args.windows)
    seen_keys = set()
    written_names = set()

    with open(args.out_lifted, "w") as out:
        with open(args.mapped_psl) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 17:
                    continue
                strand, qName, qStart, qEnd = fields[8], fields[9], int(fields[11]), int(fields[12])
                tName, tStart, tEnd = fields[13], int(fields[15]), int(fields[16])
                key = (qName, qStart, qEnd)
                names = by_key.get(key)
                if not names:
                    continue  # partial fragment, not a whole-window match — see module docstring
                seen_keys.add(key)
                for name in names:
                    if name in written_names:
                        continue  # already emitted from an earlier overlapping chain record
                    written_names.add(name)
                    out.write(f"{tName}\t{tStart}\t{tEnd}\t{name}\t0\t{strand}\n")

    with open(args.out_unmapped, "w") as out:
        for chrom, start, end, name in rows:
            if (chrom, start, end) in seen_keys:
                continue
            out.write("#no single chain record spans this window end-to-end\n")
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    n_unmapped = len(rows) - len(written_names)
    print(f"recovered names for {len(written_names)}/{len(rows)} windows; {n_unmapped} unmapped", file=sys.stderr)


if __name__ == "__main__":
    main()
