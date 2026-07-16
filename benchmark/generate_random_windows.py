#!/usr/bin/env python3
"""Generate synthetic BED files of increasing size for the runtime/scale benchmark.

Windows are placed by dividing the (optionally chromosome-filtered) genome
into N evenly-sized slots and dropping one fixed-size window at a random
offset inside each slot — guaranteed non-overlapping and spread across the
whole genome, without the O(N^2) rejection-sampling an independent-random-
placement approach would need at N=100000.

Smaller sizes are an evenly-spaced subset of the largest requested size, so
e.g. the 100-region BED is a subset of the 1000-region BED — keeps the scale
series internally consistent rather than being N unrelated random draws.
"""

import argparse
import random


def load_fai_lengths(fai_path, exclude_chroms):
    lengths = []
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom, length = fields[0], int(fields[1])
            if chrom in exclude_chroms:
                continue
            lengths.append((chrom, length))
    return lengths


def place_windows(chrom_lengths, n, window_size, seed):
    total = sum(length for _, length in chrom_lengths)
    if n * window_size > total:
        raise SystemExit(
            f"error: {n} windows of {window_size}bp ({n * window_size}bp) "
            f"won't fit in a {total}bp genome"
        )

    rng = random.Random(seed)
    slot_size = total / n
    windows = []
    # Cumulative per-chromosome start offsets in the flattened coordinate space.
    cum_offsets = []
    running = 0
    for chrom, length in chrom_lengths:
        cum_offsets.append((running, running + length, chrom, length))
        running += length

    for i in range(n):
        slot_start = i * slot_size
        slot_end = (i + 1) * slot_size
        usable = max(0.0, (slot_end - slot_start) - window_size)
        flat_pos = int(slot_start + rng.uniform(0, usable))

        for offset_start, offset_end, chrom, length in cum_offsets:
            if offset_start <= flat_pos < offset_end:
                local_start = flat_pos - offset_start
                local_start = min(local_start, length - window_size)
                windows.append((chrom, local_start, local_start + window_size, f"synthetic_{i}"))
                break

    return windows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fai", required=True, help="Reference .fai (for chromosome lengths)")
    p.add_argument(
        "--sizes",
        default="100,1000,10000,100000",
        help="Comma-separated window counts to generate (default: 100,1000,10000,100000)",
    )
    p.add_argument("--window-size", type=int, default=200, help="Window width in bp (default 200)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument(
        "--exclude-chroms",
        default="",
        help="Comma-separated chromosome names to exclude (e.g. chrM,chrY_PATERNAL)",
    )
    p.add_argument("--out-prefix", required=True, help="Output prefix; writes <prefix>_<n>.bed")
    args = p.parse_args()

    exclude = set(c for c in args.exclude_chroms.split(",") if c)
    chrom_lengths = load_fai_lengths(args.fai, exclude)
    sizes = sorted(int(s) for s in args.sizes.split(","))
    max_n = sizes[-1]

    all_windows = place_windows(chrom_lengths, max_n, args.window_size, args.seed)
    all_windows.sort(key=lambda w: (w[0], w[1]))

    for n in sizes:
        # Evenly-spaced subset of the max-N windows, so smaller BEDs nest
        # inside larger ones rather than being independent random draws.
        stride = max_n / n
        subset = [all_windows[int(i * stride)] for i in range(n)]
        subset.sort(key=lambda w: (w[0], w[1]))
        out_path = f"{args.out_prefix}_{n}.bed"
        with open(out_path, "w") as f:
            for chrom, start, end, name in subset:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
        print(f"wrote {len(subset)} windows to {out_path}")


if __name__ == "__main__":
    main()
