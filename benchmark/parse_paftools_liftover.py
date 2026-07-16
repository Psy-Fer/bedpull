#!/usr/bin/env python3
"""Recover our original window names in paftools.js liftover's output, and
report the windows it silently dropped, in the same lifted.bed/unmapped.bed
shape extract_and_score_liftover.py already consumes for UCSC liftOver.

paftools.js liftover (misc/paftools.js in the minimap2 repo) doesn't read or
propagate a BED's name column — it only uses chrom/start/end from the input
query.bed, and fabricates its own output name from the query coordinates:
`{query_chrom}_{query_start}_{query_end}`, optionally suffixed with `_t5`
and/or `_t3` if a boundary landed outside the alignment record. This script
reverses that back to our real window names by matching on (chrom, start,
end), which paftools.js reproduces exactly from its input. Chromosome names
containing underscores (e.g. chrUn_KI270751v1) are handled by matching the
*last* two underscore-separated numeric fields as start/end via regex,
rather than a naive split("_").

Any window from --windows that never appears in --raw-lifted (wrong mapq,
too-short alignment, no overlapping record, or divergence over the paftools
-d threshold) is written to --out-unmapped instead.
"""

import argparse
import re
import sys
from collections import defaultdict

NAME_RE = re.compile(r"^(.*)_(\d+)_(\d+)$")


def load_windows(path):
    """Return {(chrom, start, end): [name, ...]} and the ordered list of
    (chrom, start, end, name) rows (for reporting unmapped windows in
    input order).

    Distinct truth-set windows can legitimately share identical coordinates
    (multi-allelic/VNTR sites are common enough not to be an edge case — see
    the near-identical bug fixed in extract_and_score_liftover.py earlier
    this session), so by_key must map to a *list* of names, not a single
    name — otherwise every name but the last sharing a key silently
    disappears from both the lifted and unmapped output.
    """
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


def parse_autoname(raw):
    """('chrom', start, end, hit_t5, hit_t3) from paftools.js's fabricated
    name, or None if it doesn't match the expected shape.
    """
    s = raw
    t3 = s.endswith("_t3")
    if t3:
        s = s[: -len("_t3")]
    t5 = s.endswith("_t5")
    if t5:
        s = s[: -len("_t5")]
    m = NAME_RE.match(s)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3)), t5, t3


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--windows", required=True, help="Original input windows.bed (chrom, start, end, name)")
    p.add_argument("--raw-lifted", required=True, help="paftools.js liftover's raw stdout output")
    p.add_argument("--out-lifted", required=True, help="BED6 output with real names recovered")
    p.add_argument("--out-unmapped", required=True, help="liftOver-style '#reason' + BED unmapped output")
    args = p.parse_args()

    by_key, rows = load_windows(args.windows)
    seen_keys = set()
    written_names = set()
    n_unparsed = 0

    with open(args.out_lifted, "w") as out:
        with open(args.raw_lifted) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                chrom, start, end, autoname = fields[0], fields[1], fields[2], fields[3]
                score = fields[4] if len(fields) > 4 else "0"
                strand = fields[5] if len(fields) > 5 else "+"
                parsed = parse_autoname(autoname)
                if parsed is None:
                    n_unparsed += 1
                    continue
                key = parsed[:3]
                names = by_key.get(key)
                if not names:
                    n_unparsed += 1
                    continue
                seen_keys.add(key)
                for name in names:
                    if name in written_names:
                        continue  # already emitted from an earlier overlapping PAF record
                    written_names.add(name)
                    out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

    with open(args.out_unmapped, "w") as out:
        for chrom, start, end, name in rows:
            if (chrom, start, end) in seen_keys:
                continue
            out.write("#no qualifying alignment overlap (mapq/length/divergence filter, or no primary record covers this locus)\n")
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    n_unmapped = len(rows) - len(written_names)
    print(
        f"recovered names for {len(written_names)}/{len(rows)} windows; {n_unmapped} unmapped"
        + (f"; {n_unparsed} output line(s) couldn't be matched back to an input window" if n_unparsed else ""),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
