#!/usr/bin/env Rscript --vanilla

# Read the alignment data from a MAF file and output it in a JSON format. Also
# add the complete, unaligned genomic sequences to the output, as well as
# additional metadata.

library(seqinr)
library(jsonlite)

parse_maf <- function(path) {
  body <- trimws(readLines(path, warn = FALSE))
  body <- body[nzchar(body)]
  body <- body[!startsWith(body, "#")]
  fields <- strsplit(body, "\\s+")

  op <- vapply(fields, `[`, "", 1L)

  is_alnheader <- op == "a"
  is_seq <- op == "s"
  stopifnot(all(is_alnheader | is_seq))

  aln_fields <- fields[is_alnheader]
  score <- as.numeric(sub("^score=", "", vapply(aln_fields, `[`, "", 2L)))

  aln_id <- cumsum(is_alnheader)
  stopifnot(all(aln_id[is_seq] > 0L))

  seq_fields <- fields[is_seq]
  id <- vapply(seq_fields, `[`, "", 2L)
  start <- as.numeric(vapply(seq_fields, `[`, "", 3L)) + 1L
  len <- as.numeric(vapply(seq_fields, `[`, "", 4L))
  end <- start + len - 1L
  aligned <- vapply(seq_fields, `[`, "", 7L)
  score <- score[aln_id[is_seq]]

  data <- data.frame(
    aln_id = aln_id[is_seq],
    score = score,
    seq_id = id,
    start = start,
    end = end,
    aligned = aligned
  )
  data
}

make_blocks_main <- function(input_fasta, input_maf, input_csv) {
  seqs <- read.fasta(input_fasta,
    seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE,
    set.attributes = FALSE
  )
  csv <- read.csv(input_csv, header = TRUE)
  species_tab <- setNames(csv[["species"]], csv[["gene_id"]])

  maf_blocks <- parse_maf(input_maf)

  aln_id <- maf_blocks$aln_id
  block_score <- maf_blocks$score

  aln_blocks <- maf_blocks[setdiff(names(maf_blocks), c("score", "aln_id"))]
  aln_blocks <- split(aln_blocks, aln_id)
  names(aln_blocks) <- NULL

  dat <- list(
    seq_ids = names(seqs),
    species = setNames(species_tab[sub("[.]\\d+$", "", names(seqs))], NULL),
    lengths = nchar(seqs),
    sequences = setNames(seqs, NULL)
  )

  list(
    data = dat,
    blocks = aln_blocks,
    block_scores = block_score[!duplicated(aln_id)]
  )
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) >= 4)

  output <- make_blocks_main(args[1], args[2], args[3])
  write_json(output, args[4], auto_unbox = TRUE, pretty = 2)
}
