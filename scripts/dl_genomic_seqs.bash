#!/usr/bin/env bash

set -euo pipefail

input_csv="${1:-/dev/stdin}"
output_fasta="${2:-/dev/stdout}"

# Construct POST request body
body=$(awk -F',' 'NR > 1 && NF && $2 { print "\"" $2 "\"" }' "$input_csv" \
	| jq -s '{ids: .}')

curl -s "https://rest.ensembl.org/sequence/id?type=genomic;format=fasta" \
	-H 'Content-type:application/json' \
	-H 'Accept:text/x-fasta' \
	-X POST -d "$body" -o "$output_fasta"
