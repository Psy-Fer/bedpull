#!/usr/bin/env python3
"""Turn a windows.bed into an anchor VCF that bcftools +liftover can lift.

## Why anchors, and not the SV records themselves

bcftools +liftover is a VCF->VCF tool: it moves a variant record into a new
coordinate system and carries REF/ALT across as *payload*. That makes the
obvious integration (lift the GIAB SV record, read its new POS, rebuild the
window from the carried ALT length) worthless as a measurement — the
extracted length would be reconstructed from the allele the tool was handed,
so it would equal truth.tsv's expected_length by construction every time the
lift succeeded. It would look like a 100% score and measure nothing.

Instead we ask bcftools the same question every other tool in this harness is
asked: *where do these two reference coordinates land in the target
assembly?* Each window contributes two single-base anchor records, one at the
window's first base and one at its last. The lifted interval between them is
an emergent property of the chain alignment, not of any allele we supplied —
if the assembly carries a 563bp insertion between the anchors, the lifted
interval comes out 563bp longer, and if the tool truncates at a chain
boundary it doesn't. That is exactly the quantity liftOver, pslMap and
bedpull are already scored on, so the comparison is like-for-like.

This is also not an abuse of the tool: +liftover's own internal mechanism is
anchor-based (see its UnmappedAnchor5/UnmappedAnchor3/MismatchAnchors/
ApartAnchors reject reasons), so lifting explicit anchors exercises the same
machinery its indel path uses.

## Why ALT defaults to '.'

An anchor exists to carry a coordinate, not an allele. With a concrete ALT
base, every anchor that happens to land on a real hs1-vs-HG002 SNV engages
the plugin's allele bookkeeping — SWAP when the assembly carries our
arbitrarily chosen ALT, or a reference-mismatch outcome when it carries a
third base — and could in principle be rejected for a reason that has
nothing to do with whether the coordinate lifted. ALT='.' (no alternate
allele) sidesteps that path entirely and isolates pure coordinate lifting.

Measured on the full 24,965-window PATERNAL set, the two styles are near
identical: both lift 21,662 windows and reject 3,303, differing on 5
windows' coordinates (0.02%). `--alt-style snv` does trip 22 swaps and 175
reference-added anchors, so the allele path really is engaged — it just
doesn't change the result here. ref-only stays the default as the cleaner
control, not because the choice is load-bearing.

Writes an unsorted VCF; the caller is expected to pipe it through
`bcftools sort` (run_bcftools_liftover_pipeline.sh does).
"""

import argparse
import subprocess
import sys

# Anchor IDs are "<window name>#L" / "<window name>#R". '#' is legal mid-field
# in a VCF ID and is not used by this harness's window names (which are
# built from SVTYPE/chrom/pos/len joined with '_'), so it round-trips through
# bcftools as an unambiguous separator.
ID_SEP = "#"

# Transversion, so the ALT can never coincide with the REF base.
TRANSVERSION = {"A": "C", "C": "A", "G": "T", "T": "G"}


def read_windows(path):
    """Yield (chrom, start, end, name) from a BED4 windows file."""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                sys.exit(f"windows BED needs 4 columns (name required), got: {line!r}")
            yield fields[0], int(fields[1]), int(fields[2]), fields[3]


def anchor_positions(start, end):
    """1-based VCF positions of a window's first and last base.

    BED is 0-based half-open, so the first base is `start` (1-based start+1)
    and the last base is `end - 1` (1-based `end`). Anchoring on the last
    base rather than one past it keeps both anchors inside the window, which
    matters because a position one past the window end can fall outside the
    alignment block the window sits in.
    """
    return start + 1, end


def fetch_bases(samtools_bin, fasta, regions, tmp_prefix):
    """Batch one `samtools faidx -r` call for all anchor regions.

    Returns {region_string: uppercased base}. Soft-masked (lowercase)
    reference bases are uppercased because the VCF REF must match the
    FASTA case-insensitively but bcftools compares literally.
    """
    region_file = f"{tmp_prefix}.anchor_regions.txt"
    with open(region_file, "w") as f:
        for region in regions:
            f.write(region + "\n")

    proc = subprocess.run(
        [samtools_bin, "faidx", "-r", region_file, fasta],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        sys.exit(f"samtools faidx failed on the source FASTA: {proc.stderr}")

    bases = {}
    header = None
    seq = []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if header is not None:
                bases[header] = "".join(seq).upper()
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if header is not None:
        bases[header] = "".join(seq).upper()
    return bases


def read_fai(fai_path):
    """Return [(contig, length), ...] in FASTA order, for the VCF header."""
    contigs = []
    with open(fai_path) as f:
        for line in f:
            fields = line.split("\t")
            contigs.append((fields[0], fields[1]))
    return contigs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--windows", required=True, help="windows.bed (BED4, name in column 4)")
    p.add_argument("--src-fasta", required=True, help="source/reference FASTA the windows are in (needs .fai)")
    p.add_argument("--out-vcf", required=True, help="output anchor VCF path")
    p.add_argument("--samtools", default="samtools", help="samtools binary (default: on PATH)")
    p.add_argument(
        "--alt-style",
        choices=["ref-only", "snv"],
        default="ref-only",
        help="ALT for each anchor: '.' (default, isolates coordinate lifting) or a transversion SNV",
    )
    args = p.parse_args()

    windows = list(read_windows(args.windows))
    if not windows:
        sys.exit(f"no windows read from {args.windows}")

    wanted = []
    seen = set()
    for chrom, start, end, _name in windows:
        for pos in anchor_positions(start, end):
            region = f"{chrom}:{pos}-{pos}"
            if region not in seen:
                seen.add(region)
                wanted.append(region)

    bases = fetch_bases(args.samtools, args.src_fasta, wanted, args.out_vcf)

    contigs = read_fai(f"{args.src_fasta}.fai")
    contig_order = {c: i for i, (c, _n) in enumerate(contigs)}

    records = []
    n_skipped = 0
    for chrom, start, end, name in windows:
        left, right = anchor_positions(start, end)
        for tag, pos in (("L", left), ("R", right)):
            ref = bases.get(f"{chrom}:{pos}-{pos}", "")
            if len(ref) != 1 or ref not in "ACGTN":
                # An anchor we can't give a valid REF for can't be lifted; drop
                # it and let bcftools_lifted_to_bed.py report the window as
                # unmapped (one missing anchor is enough to fail the window).
                n_skipped += 1
                continue
            if args.alt_style == "snv":
                alt = TRANSVERSION.get(ref, "A")
            else:
                alt = "."
            records.append((contig_order.get(chrom, len(contigs)), pos, chrom, f"{name}{ID_SEP}{tag}", ref, alt))

    records.sort(key=lambda r: (r[0], r[1]))

    with open(args.out_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig, length in contigs:
            out.write(f"##contig=<ID={contig},length={length}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _order, pos, chrom, anchor_id, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t{anchor_id}\t{ref}\t{alt}\t.\t.\t.\n")

    msg = f"{len(windows)} windows -> {len(records)} anchor records ({args.alt_style})"
    if n_skipped:
        msg += f", {n_skipped} anchor(s) skipped (no valid reference base)"
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
