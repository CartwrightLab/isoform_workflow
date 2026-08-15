CURL_BIN=curl

default: all

all:

.PHONY: default all

.SECONDARY:

data/gene_info.json: config/gene_list.txt
	bash scripts/dl_gene_info.bash $< $@

data/gene_info/%.json: data/gene_info.json
	bash scripts/mk_canonical_exons.bash $* $< $@

# Download homology information for a human gene symbol from all of primates
data/homology_raw/%.json:
	$(CURL_BIN) -H 'Content-type:application/json' -o $@ \
	'https://rest.ensembl.org/homology/symbol/homo_sapiens/$*?target_taxon=9443;type=orthologues;sequence=none;format=full'

# Extract ENSEMBL IDs from homology information
data/homology/%.csv: data/homology_raw/%.json
	@mkdir -p $(@D)
	bash scripts/mk_homology_csv.bash $< $@

# Download unaligned data for each species
data/fasta_raw/%.fasta: data/homology/%.csv
	@mkdir -p $(@D)
	bash scripts/dl_genomic_seqs.bash $< $@

# Align gene regions using lastz
data/genomic_maf/%.maf: data/fasta_raw/%.fasta
	@mkdir -p $(@D)
	bash scripts/mk_genomic_maf.bash $< $@

# Produce aligned blocks from the delta alignments
data/genomic_blocks/%.json: data/fasta_raw/%.fasta data/genomic_maf/%.maf data/homology/%.csv
	@mkdir -p $(@D)
	Rscript --vanilla scripts/mk_genomic_blocks.R $^ $@

data/cds_aln/%.fasta: data/gene_info/%.json data/genomic_aln/%.json
	@mkdir -p $(@D)
	Rscript --vanilla scripts/mk_cds_aln.R $^ $@
