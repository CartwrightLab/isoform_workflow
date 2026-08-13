CURL_BIN=curl

default: all

all:

.PHONY: default all

.SECONDARY:

data/gene_list.json: config/gene_list.txt
	test $$(wc -l < $<) -le 1000 # REST API only supports 1000 requests at a time. TODO: refactor
	jq -R -s -c 'split("\n") | map(select(length > 0)) | {symbols: ., expand: 1}' $< \
	 | $(CURL_BIN) -s -H 'Content-type:application/json' -H 'Accept:application/json' \
	  -X POST -d @- 'https://rest.ensembl.org/lookup/symbol/homo_sapiens' -o $@

# Download homology information for a human gene symbol from all of primates
data/homology_raw/%.json:
	$(CURL_BIN) -H 'Content-type:application/json' -o $@ \
	'https://rest.ensembl.org/homology/symbol/homo_sapiens/$*?target_taxon=9443;type=orthologues;sequence=none;format=full'

# Extract ENSEMBL IDs from homology information
data/homology/%.csv: data/homology_raw/%.json
	bash scripts/mk_homology_csv.bash $< $@

# Download unaligned data for each species
data/fasta_raw/%.fasta: data/homology/%.csv
	bash scripts/dl_genomic_seqs.bash $< $@

data/genomic_aln_delta/%.delta: data/fasta_raw/%.fasta
	bash scripts/mk_genomic_delta.bash $< $@


data/genomic_aln/%.json : data/fasta_raw/%.fasta data/genomic_aln_delta/%.delta data/homology/%.csv
	Rscript --vanilla scripts/mk_genomic_blocks.R $^ $@
