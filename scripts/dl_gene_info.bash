#!/usr/bin/env bash

# Download information for a set of human genes from ensembl using its REST API.
# Use a POST request to download information in bulk.

set -euo pipefail

input_gene="$1"
output_json="${2:-/dev/stdout}"

curl -f -s "https://rest.ensembl.org/lookup/symbol/homo_sapiens/${input_gene}?expand=1" \
	-H 'Content-type:application/json' \
	-H 'Accept:application/json' \
	-o "$output_json"
