#!/usr/bin/env bash

# Identify the canonical transcript for a human gene and extract data about
# its translation and exons.

set -euo pipefail

input_json="${1:-/dev/stdin}"
output_json="${2:-/dev/stdout}"

jq -c "(.Transcript[] | select(.is_canonical == 1)) as \$t
  | \$t.Exon as \$e
  | {
      id: .id,
      version: .version,
      seq_region_name: .seq_region_name,
      start: .start,
      end: .end,
      strand: .strand,
      translation: (
        \$t.Translation
        | {
            id: .id,
            version: .version,
            start: .start,
            end: .end
          }
      ),
      exons: {
        id: (\$e | map(.id)),
        version: (\$e | map(.version)),
        start: (\$e | map(.start)),
        end: (\$e | map(.end))
      }
    }" \
  "${input_json}" > "${output_json}"
