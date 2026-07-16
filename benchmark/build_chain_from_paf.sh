#!/usr/bin/env bash
# Convert a minimap2 PAF alignment (with cg:Z: CIGAR tags — the same input
# bedpull's --paf mode consumes) into a UCSC .chain file, for the liftOver
# comparison leg of the benchmark.
#
# Uses paf2chain (https://github.com/AndreaGuarracino/paf2chain), a small
# purpose-built Rust tool with exactly one relevant flag. This is deliberately
# the *only* supported path here rather than a hand-rolled UCSC
# axtChain/chainNet/chainSort fallback: that classic pipeline expects
# axt/PSL input (not PAF) and getting the conversion flags right needs
# hands-on verification we can't do without the tool installed. If you need
# that path, build it yourself and drop the resulting .chain where this
# script would have written one — everything downstream just wants a
# .chain(.gz) file, not this specific tool.
#
# Usage:
#   build_chain_from_paf.sh <in.paf> <out.chain>
#
# Install paf2chain if missing:
#   git clone https://github.com/AndreaGuarracino/paf2chain
#   cd paf2chain && cargo install --force --path .

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $(basename "$0") <in.paf> <out.chain>" >&2
    exit 1
fi

paf_path="$1"
chain_path="$2"

if [ ! -f "$paf_path" ]; then
    echo "error: PAF file not found: $paf_path" >&2
    exit 1
fi

if ! command -v paf2chain >/dev/null 2>&1; then
    cat >&2 <<'EOF'
error: paf2chain not found on PATH.

Install it with:
  git clone https://github.com/AndreaGuarracino/paf2chain
  cd paf2chain && cargo install --force --path .

Then re-run this script.
EOF
    exit 1
fi

mkdir -p "$(dirname "$chain_path")"
paf2chain -i "$paf_path" > "$chain_path"

n_chains=$(grep -c '^chain ' "$chain_path" || true)
echo "wrote $n_chains chain(s) to $chain_path" >&2
