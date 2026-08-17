#!/usr/bin/env Rscript --vanilla

# Construct a draft multiple sequence alignment of a CDS region based on
# alignment blocks produced by lastz.
#
# ## Steps
#
# 1. Blocks are grouped by query sequence ids.
# 2. Blocks are sorted by score and overlaps are discarded.
#   - The higher scoring block is kept.
#   - Overlap is only checked on the reference strand.
# 3. A global alignment is constructed from the alignment blocks.
#   - Missing sections in the aligned reference are replaced by the reference
#   - Missing sections in the aligned query are replaced by "n"
# 4. The pairwise alignment is trimmed to just contain the coding exons.
# 5. A "reference-based" MSA is constructed assuming star-tree.

library(seqinr)
library(jsonlite)

filter_overlaps <- function(start, end) {
  n <- length(start)
  stopifnot(length(end) == n)

  ov <- outer(start, end, `<=`) & outer(end, start, `>=`)

  keep <- logical(n)
  for (i in seq_len(n)) {
    keep[i] <- !any(ov[i, keep])
  }

  keep
}

build_global_alignment <- function(ref_seq, blocks, scores) {
  len <- nchar(ref_seq)

  # create wide table
  S1 <- sapply(blocks, \(x) x$start[1])
  E1 <- sapply(blocks, \(x) x$end[1])
  A1 <- sapply(blocks, \(x) x$aligned[1])
  S2 <- sapply(blocks, \(x) x$start[2])
  E2 <- sapply(blocks, \(x) x$end[2])
  A2 <- sapply(blocks, \(x) x$aligned[2])

  tab <- data.frame(
    S1 = S1, E1 = E1, A1 = A1,
    S2 = S2, E2 = E2, A2 = A2,
    score = scores
  )

  # reorder blocks by scores
  o <- order(tab$score, decreasing = TRUE)
  tab <- tab[o, ]

  # identify the primary blocks, subset, and reorder
  p <- filter_overlaps(tab$S1, tab$E1) # & filter_overlaps(tab$S2, tab$E2)
  tab <- tab[p, ]
  o <- order(tab$S1)
  tab <- tab[o, ]

  # Spacers between blocks
  Sz <- c(1L, tab$E1 + 1L)
  Ez <- c(tab$S1 - 1L, len)
  Az <- substring(ref_seq, Sz, Ez)
  Nz <- strrep("n", nchar(Az))

  # build reference alignment
  ref <- paste0(c("", tab$A1), Az, collapse = "")
  qry <- paste0(c("", tab$A2), Nz, collapse = "")
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
  blocks <- read_json(input_blocks, simplifyVector = TRUE)

  ref_id <- blocks$data$seq_ids[1]
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
  cds_len <- sum(exon_rights - exon_lefts + 1L)

  # Sanity checks
  stopifnot(length(exon_lefts) == length(exon_rights) && length(exon_lefts) > 0L)
  stopifnot(!is.unsorted(exon_lefts, strictly = TRUE))
  stopifnot(!is.unsorted(exon_rights, strictly = TRUE))
  stopifnot(all(exon_lefts <= exon_rights))

  # Group blocks by query_id
  query_ids <- sapply(blocks$blocks, \(x) x$seq_id[2])
  grouped_blocks <- split(blocks$blocks, query_ids)
  grouped_scores <- split(blocks$block_scores, query_ids)
  species_ids <- blocks$data$species[match(names(grouped_blocks), blocks$data$seq_ids)]
  seq_names <- paste(species_ids, names(grouped_blocks))

  cds_matrix <- matrix("", nrow = cds_len + 1L, ncol = length(seq_names))

  for (i in seq_along(grouped_blocks)) {
    block_group <- grouped_blocks[[i]]
    score_group <- grouped_scores[[i]]

    # sanity check
    group_ref_id <- sapply(block_group, \(x) x$seq_id[1])
    stopifnot(all(group_ref_id == ref_id))

    a <- build_global_alignment(ref_seq, block_group, score_group)
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
