#!/usr/bin/env bash

# Download genomic seqeunces for a set of homologous genes. This batches genes
# using a POST request, but is limited to 50 genes at a time.

set -euo pipefail

input_csv="${1:-/dev/stdin}"
output_fasta="${2:-/dev/stdout}"

# Construct POST request body
body=$(awk -F',' 'NR > 1 && NF && $2 { print "\"" $2 "\"" }' "$input_csv" \
	| jq -s '{ids: .}')

curl -f -s "https://rest.ensembl.org/sequence/id?type=genomic;format=fasta" \
	-H 'Content-type:application/json' \
	-H 'Accept:text/x-fasta' \
	-X POST -d "$body" -o "$output_fasta"
