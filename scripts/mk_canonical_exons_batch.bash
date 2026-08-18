#!/usr/bin/env bash

# Identify the canonical transcript for a human gene and extract data about
# its translation and exons.

set -euo pipefail

input_gene="${1}"
input_json="${2:-/dev/stdin}"
output_json="${3:-/dev/stdout}"

jq -c ".${input_gene} as \$g
  | (\$g.Transcript[] | select(.is_canonical == 1)) as \$t
  | \$t.Exon as \$e
  | \$g
  | {id, version, seq_region_name, start, end, strand,
    translation: (\$t.Translation | {id, version, start, end}),
    exons: {
      id: (\$e | map(.id)),
      version: (\$e | map(.version)),
      start: (\$e | map(.start)),
      end: (\$e | map(.end))
    }
  }" \
  "${input_json}" > "${output_json}"
