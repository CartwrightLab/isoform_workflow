#!/usr/bin/env bash

set -euo pipefail

# Script arguments and defaults
input_fasta="${1:-/dev/stdin}"
output_delta="${2:-/dev/stdout}"

# Find nucmer binary. Prefer `numcer` and then `mummer-nucmer`.
# Use `nucmer` if neither is found for readable error message.
NUCMER_BIN=$(command -v nucmer || command -v mummer-nucmer || echo "nucmer")
FILTER_BIN=$(command -v delta-filter || command -v mummer-delta-filter || echo "delta-filter")

# Create temporary file and delete it on exit.
ref_fasta=$(mktemp --suffix=.fasta)
trap 'rm -f "$ref_fasta"' EXIT

# Extract the first sequence as the reference.
awk '/^>/{n++} n>1{exit} {print}' "$input_fasta" > "$ref_fasta"

# Align the whole file (all sequences, including the reference itself) against it
# Filter the resulting delta file for global best hit for each query
"$NUCMER_BIN" -f --delta=/dev/stdout "$ref_fasta" "$input_fasta" \
  | "$FILTER_BIN" -q -g /dev/stdin > "$output_delta"


#  lastz ref.fasta AREG.fasta --gfextend --chain --gapped --strand=plus --inner=1000 --format=maf