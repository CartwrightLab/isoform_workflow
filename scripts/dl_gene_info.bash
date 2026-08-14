#!/usr/bin/env bash

set -euo pipefail

input_list="${1:-/dev/stdin}"
output_json="${2:-/dev/stdout}"

# REST API only supports 1000 requests at a time. TODO: refactor
test $(wc -l < ${input_list}) -le 1000 || exit 1

# Construct POST request body
body=$(jq -R -s -c 'split("\n") | map(select(length > 0)) | {symbols: ., expand: 1}' "${input_list}")

curl -s 'https://rest.ensembl.org/lookup/symbol/homo_sapiens' \
	-H 'Content-type:application/json' \
	-H 'Accept:application/json' \
	-X POST -d "$body" -o "$output_json"
