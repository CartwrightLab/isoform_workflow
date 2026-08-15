default: all

all:

.PHONY: default all

.SECONDARY:

data/gene_info.json: config/gene_list.txt
	@mkdir -p $(@D)
	bash scripts/dl_gene_info.bash $< $@

data/gene_info/%.json: data/gene_info.json
	@mkdir -p $(@D)
	bash scripts/mk_canonical_exons.bash $* $< $@

# Download homology information for a human gene symbol from all of primates
data/homology_raw/%.json:
	@mkdir -p $(@D)
	bash scripts/dl_homology.bash $* $@

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

data/cds_aln/%.fasta: data/gene_info/%.json data/genomic_blocks/%.json
	@mkdir -p $(@D)
	Rscript --vanilla scripts/mk_cds_aln.R $^ $@
