default: all

all:

.PHONY: default all

.SECONDARY:

## Comment out downloading gene info as a batch.
# data/gene_info.json: config/gene_list.txt
# 	@mkdir -p $(@D)
# 	bash scripts/dl_gene_info_batch.bash $< $@

# data/gene_info/%.json: data/gene_info.json
# 	@mkdir -p $(@D)
# 	bash scripts/mk_canonical_exons_batch.bash $* $< $@

# Download information about a human gene
data/gene_info_raw/%.json:
	@mkdir -p $(@D)
	bash scripts/dl_gene_info.bash $* $@	

# Download homology information for a human gene symbol from all of primates
data/homology_raw/%.json:
	@mkdir -p $(@D)
	bash scripts/dl_homology.bash $* $@

# Download unaligned data for each species
data/fasta_raw/%.fasta: data/homology/%.csv
	@mkdir -p $(@D)
	bash scripts/dl_genomic_seqs.bash $< $@

# Identify the human canonical sequence, its translation, and its exons.
data/gene_info/%.json: data/gene_info_raw/%.json
	@mkdir -p $(@D)
	bash scripts/mk_canonical_exons.bash $< $@

# Extract ENSEMBL IDs from homology information
data/homology/%.csv: data/homology_raw/%.json
	@mkdir -p $(@D)
	bash scripts/mk_homology_csv.bash $< $@

# Align CDS sequences using CESAR
data/fasta_cds/%.fasta: data/fasta_raw/%.fasta data/gene_info/%.json data/homology/%.csv
	@mkdir -p $(@D)
	Rscript --vanilla scripts/mk_cds_aln.R $^ $@
