#!/usr/bin/env Rscript --vanilla

library(seqinr)
library(jsonlite)

build_global_alignment <- function(ref_seq, aln) {
  len <- nchar(ref_seq)

  S1 <- sapply(aln$blocks, \(x) x$reference$start)
  E1 <- sapply(aln$blocks, \(x) x$reference$end)
  A1 <- sapply(aln$blocks, \(x) x$reference$aligned)
  S2 <- sapply(aln$blocks, \(x) x$query$start)
  E2 <- sapply(aln$blocks, \(x) x$query$end)
  A2 <- sapply(aln$blocks, \(x) x$query$aligned)

  # Spacers between blocks
  Sz <- c(1L, E1 + 1L)
  Ez <- c(S1 - 1L, len)
  Az <- substring(ref_seq, Sz, Ez)
  Nz <- strrep("n", nchar(Az))

  # build reference alignment
  ref <- paste0(c("", A1), Az, collapse = "")
  qry <- paste0(c("", A2), Nz, collapse = "")
  stopifnot(nchar(ref) == nchar(qry))

  c(ref = ref, qry = qry)
}

pos_to_col <- function(x, seq) {
  p <- gregexpr("-", seq)[[1]]
  if (p[1] == -1) {
    return(x)
  }
  v <- p - seq_along(p)
  x + findInterval(x - 1, v)
}

subset_alignment <- function(aln, start, stop) {
  start <- pos_to_col(start, aln[1])
  stop <- pos_to_col(stop, aln[1])
  sapply(aln, \(x) paste0(substring(x, start, stop), collapse = ""))
}

normalize_alignment <- function(aln) {
  aln <- strsplit(aln, "", fixed = TRUE)
  ref0 <- aln[[1]]
  qry0 <- aln[[2]]
  stopifnot(length(ref0) == length(qry0))
  n <- length(ref0)
  o <- ref0 != "-" & qry0 != "-"
  ref <- ref0[ref0 != "-"]
  qry <- qry0[qry0 != "-"]

  ref_m <- ref0[o]
  qry_m <- qry0[o]
  len_r <- length(ref)
  len_q <- length(qry)

  index_r <- integer(len_r)
  index_q <- integer(len_q)
  r <- q <- k <- 1L
  for (m in seq_along(ref_m)) {
    while (ref[r] != ref_m[m]) {
      index_r[r] <- k
      r <- r + 1L
      k <- k + 1L
    }
    while (qry[q] != qry_m[m]) {
      index_q[q] <- k
      q <- q + 1L
      k <- k + 1L
    }
    # Column k corresponds to match m.
    index_r[r] <- k
    index_q[q] <- k
    r <- r + 1L
    q <- q + 1L
    k <- k + 1L
  }
  # Handle any trailing gaps.
  if (r <= len_r) {
    nn <- len_r - r + 1L
    index_r[r:len_r] <- seq.int(k, length.out = nn)
    k <- k + nn
  }
  if (q <= len_q) {
    nn <- len_q - q + 1L
    index_q[q:len_q] <- seq.int(k, length.out = nn)
    k <- k + nn
  }
  stopifnot(k == n + 1L)

  # Map sequences to alignment columns
  ref_o <- rep.int("-", n)
  ref_o[index_r] <- ref
  ref_o <- paste0(ref_o, collapse = "")
  qry_o <- rep.int("-", n)
  qry_o[index_q] <- qry
  qry_o <- paste0(qry_o, collapse = "")

  setNames(c(ref_o, qry_o), names(aln))
}

shatter_alignment <- function(aln) {
  aln <- strsplit(aln, "")
  # Assign each column in qry to a position in ref. Insertions are assigned
  # to the closest "match" upstream (to the left) of the insertion.
  p <- cumsum(c("N", aln[[1]]) != "-")
  s <- split(c("N", aln[[2]]), p)
  # Join the columns in alt together to make "alleles".
  sapply(s, paste0, collapse = "")
}

make_cds_main <- function(input_info, input_blocks) {
  info <- read_json(input_info, simplifyVector = TRUE)
  blocks <- read_json(input_blocks, simplifyVector = TRUE, simplifyDataFrame = FALSE)

  ref_id <- blocks$data$ids[1]
  ref_seq <- blocks$data$sequences[1]

  seq_start <- info$start
  seq_length <- info$end - info$start + 1L
  cds_start <- info$translation$start - seq_start + 1L
  cds_end <- info$translation$end - seq_start + 1L
  exon_starts <- pmax(info$exons$start - seq_start + 1L, cds_start)
  exon_ends <- pmin(info$exons$end - seq_start + 1L, cds_end)

  # identify exons that overlap with the CDS
  o <- exon_starts <= exon_ends

  if (info$strand == -1) {
    exon_lefts <- seq_length - exon_ends[o] + 1L
    exon_rights <- seq_length - exon_starts[o] + 1L
  } else {
    exon_lefts <- exon_starts[o]
    exon_rights <- exon_ends[o]
  }
  # Sanity checks
  stopifnot(length(exon_lefts) == length(exon_rights) && length(exon_lefts) > 0L)
  stopifnot(!is.unsorted(exon_lefts, strictly = TRUE))
  stopifnot(!is.unsorted(exon_rights, strictly = TRUE))
  stopifnot(all(exon_lefts <= exon_rights))

  cds_len <- sum(exon_rights - exon_lefts + 1L)

  query_ids <- sapply(blocks$alignments, \(x) x$query$id)
  species_ids <- blocks$data$species[match(query_ids, blocks$data$ids)]
  seq_names <- paste(species_ids, query_ids)

  cds_matrix <- matrix("", nrow = cds_len + 1L, ncol = length(seq_names))

  for (i in seq_along(blocks$alignments)) {
    aln <- blocks$alignments[[i]]
    stopifnot(aln$reference$id == ref_id)
    a <- build_global_alignment(ref_seq, aln)
    cds <- subset_alignment(a, exon_lefts, exon_rights)
    # cds <- normalize_alignment(cds)
    cds_matrix[, i] <- shatter_alignment(cds)
  }

  # Remove immortal link
  cds_matrix[1, ] <- sub("^N", "", cds_matrix[1, ])

  # Pad alleles
  nz <- nchar(cds_matrix)
  ng <- apply(nz, 1, max) - nz
  pg <- strrep("-", ng)
  cds_matrix[] <- paste0(cds_matrix, pg)
  cds <- apply(cds_matrix, 2, paste0, collapse = "")

  # Move human to the top
  o <- which(species_ids == "homo_sapiens")
  seq_names <- c(seq_names[o], seq_names[-o])
  cds <- c(cds[o], cds[-o])

  setNames(cds, seq_names)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) >= 3)

  seqs <- make_cds_main(args[1], args[2])

  write.fasta(as.list(seqs), names(seqs), args[3], as.string = TRUE)
}
