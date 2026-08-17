#!/usr/bin/env Rscript --vanilla

library(jsonlite)

# read a fasta file into a character vector
read_fasta <- function(path, text) {
  lines <- trimws(readLines(path))
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "#")]
  is_header <- grepl("^>", lines)
  seq_id <- cumsum(is_header)

  # sanity check
  stopifnot(all(seq_id[!is_header] > 0L))

  # extract names
  n <- sub("^>(\\S+).*", "\\1", lines[is_header])

  # build sequences
  g <- split(lines[!is_header], seq_id[!is_header])
  g <- vapply(g, paste0, "", collapse = "", USE.NAMES = FALSE)
  g <- gsub("\\s", "", g)

  setNames(g, n)
}

run_cesar <- function(exons, seqs) {
  exons <- toupper(exons)

  g <- rep(seq_along(exons), nchar(exons))
  cds <- unlist(strsplit(exons, "", fixed = TRUE))

  # Identify and mark split codons
  y <- seq(1, length(g), 3)
  b <- (g[y] == g[y + 1] & g[y] == g[y + 2])
  o <- which(!b)
  p <- 3 * (o - 1)
  cds[c(p + 1, p + 2, p + 3)] <- tolower(cds[c(p + 1, p + 2, p + 3)])
  exons <- sapply(split(cds, g), paste0, collapse = "")

  # Write input file for CESAR
  ref_text <- sprintf(">%s%d\n%s", "exon", seq_along(exons), exons)
  qry_text <- sprintf(">%s\n%s", names(seqs), seqs)
  temp_fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_fasta))
  writeLines(c(ref_text, "####", qry_text), temp_fasta)

  # Call CESAR
  output <- system2("bin/cesar", args = temp_fasta, stdout = TRUE, stderr = FALSE)
  exit_status <- attr(output, "status")
  stopifnot(is.null(exit_status))

  # Parse output
  output <- gsub(" ", "?", output)
  output_fasta <- read_fasta(textConnection(output))
  output_fasta
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

align_cds_main <- function(input_fasta, input_info, input_csv) {
  fasta <- read_fasta(input_fasta)
  info <- read_json(input_info, simplifyVector = TRUE)
  csv <- read.csv(input_csv, header = TRUE)
  species_tab <- setNames(csv[["species"]], csv[["gene_id"]])

  ref_id <- names(fasta)[1]
  ref_seq <- fasta[[1]]

  # Identify locations of coding exons in ref_seq
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

  # Extract Exons
  exons <- substring(ref_seq, exon_lefts, exon_rights)

  # Call CESAR
  blocks <- run_cesar(exons, fasta)

  # Process CESAR output
  aligned_ref <- toupper(blocks[seq.int(1, length(blocks), 2)])
  aligned_qry <- toupper(blocks[seq.int(2, length(blocks), 2)])
  spacer <- strrep("n", nchar(exons))

  segment_pos <- gregexpr("[^?]+", aligned_ref)
  segment_ref <- regmatches(aligned_ref, segment_pos)
  segment_qry <- regmatches(aligned_qry, segment_pos)

  # add spacers for mixing exons
  for (i in seq_along(segment_qry)) {
    ref <- gsub("-", "", segment_ref[[i]], fixed = TRUE)
    o <- match(ref, exons)

    qry <- spacer
    qry[o] <- segment_qry[[i]]
    segment_qry[[i]] <- qry

    ref <- exons
    ref[o] <- segment_ref[[i]]
    segment_ref[[i]] <- ref
  }

  # construct pairwise alignments
  refs <- sapply(segment_ref, paste0, collapse = "")
  qrys <- sapply(segment_qry, paste0, collapse = "")

  species_ids <- species_tab[sub("[.]\\d+$", "", names(qrys))]
  seq_names <- paste(species_ids, names(qrys))

  # construct multiple sequence alignments
  cds_matrix <- matrix("", nrow = cds_len + 1L, ncol = length(seq_names))
  for (i in seq_along(qrys)) {
    cds_matrix[, i] <- shatter_alignment(c(refs[[i]], qrys[[i]]))
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
  stopifnot(length(args) >= 4)

  # Get aligned sequences
  text <- align_cds_main(args[1], args[2], args[3])

  # Wrap sequences
  text <- gsub("(.{80})", "\\1\n", text)

  # Build Lines
  lines <- sprintf(">%s\n%s", names(text), text)

  # Save output
  writeLines(lines, args[4])
}
