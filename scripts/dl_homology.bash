#!/usr/bin/env bash

# Download information about primate orthologus for a human gene from ENSEMBL.

set -euo pipefail

input_gene="$1"
output_json="${2:-/dev/stdout}"

url="https://rest.ensembl.org/homology/symbol/homo_sapiens/${input_gene}"
url="${url}?target_taxon=9443;type=orthologues;sequence=none;format=full"

curl -f -s -H 'Content-type:application/json' "${url}" -o "$output_json"
