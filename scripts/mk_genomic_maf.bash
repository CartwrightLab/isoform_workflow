#!/usr/bin/env bash

set -euo pipefail

# Script arguments and defaults
input_fasta="${1:-/dev/stdin}"
output_maf="${2:-/dev/stdout}"

# Create temporary file and delete it on exit.
ref_fasta=$(mktemp --suffix=.fasta)
trap 'rm -f "$ref_fasta"' EXIT

# Extract the first sequence as the reference.
awk '/^>/{n++} n>1{exit} {print}' "$input_fasta" > "$ref_fasta"

# Align the whole file (all sequences, including the reference itself) against it
lastz "$ref_fasta" "$input_fasta" \
  --gfextend \
  --chain \
  --gapped \
  --strand=plus \
  --inner=1000 \
  --format=maf \
  > "$output_maf"
