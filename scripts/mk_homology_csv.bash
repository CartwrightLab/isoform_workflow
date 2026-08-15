#!/usr/bin/env bash

set -euo pipefail

input_json="${1:-/dev/stdin}"
output_csv="${2:-/dev/stdout}"

# This jq filter transforms the downloaded homology information into a CSV that
# contains a species name and a ENSEMBL gene ID 
jq -r '.data[0] as $d 
  | ($d.homologies[0].source | [.species, $d.id]) as $human
  | [$human] + [$d.homologies[] | select(.type == "ortholog_one2one") | [.target.species, .target.id]] 
  | (["species","gene_id"], .[]) 
  | join(",")' "${input_json}" > "${output_csv}"

