# Isoform Workflow

## Requirements

- GNU Make
- Bash Shell
- jq
- curl
- R (packages: jsonlite)
- CESAR 2.0 (https://github.com/hillerlab/CESAR2.0)

## Installation

This workflow can be downloaded from GitHub using git clone.

```shell
git clone https://github.com/CartwrightLab/isoform_workflow.git
```

## Usage

This workflow is built around GNU Make.

```shell
cd isoform_workflow
make data/fasta_cds/AREG.fasta
```

This command downloads information from Ensembl for a human gene. It will also
download the genomic sequences of the human gene and all the one-to-one
orthologs identified in other primates. It then constructs an alignment of
the human canonical CDS with the other species.

## Directory structure

- `fasta_cds`: Multiple sequence alignments of CDS sequences.
- `fasta_raw`: Raw genomic sequences downloaded from Ensembl.
- `gene_info`: Information about canonical isoforms.
- `gene_info_raw`: Raw gene information downloaded from Ensembl.
- `homology`: Information about one-to-one orthologous genes.
- `homology_raw`: Raw homology information downloaded from Ensembl.

## Potential Issues and Improvements

The CESAR alignments are currently used as is and low quality "hits" are not
filtered out.

CESAR can take a while to run, but provides better results on our problem
than lastz, nucmer, and minimap2. It's possible that exonerate might be able
to replace CESAR, or that we could rethink our approach to using lastz.

The results in fasta_cds are quick, reference-based alignments. A full
tree-based, codon-aware multiple sequence aligner would produce more optimal
results. But would likely require masking of "frameshifts" relative to the
reference first. Without writing a new aligner, we can probably do 
fasta_cds -> coati -> masking -> prank. However, coati might be redundant
because CESAR is already codon-aware.
