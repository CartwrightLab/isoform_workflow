#!/usr/bin/env Rscript --vanilla

library(seqinr)

parse_delta <- function(path) {
  body <- trimws(readLines(path, warn = FALSE)[-(1:2)])
  body <- body[nzchar(body)]
  fields <- strsplit(body, "\\s+")
  nfields <- lengths(fields)

  is_seqheader <- startsWith(body, ">")
  is_coord <- !is_seqheader & nfields == 7L
  is_deltaval <- !is_seqheader & nfields == 1L
  stopifnot(is_seqheader | is_coord | is_deltaval)

  seq_id <- cumsum(is_seqheader)
  aln_id <- cumsum(is_coord)
  stopifnot(all(seq_id[!is_seqheader] > 0L), all(aln_id[is_deltaval] > 0L))

  # One entry per sequence pair, indexed by seq_id
  seq_fields <- fields[is_seqheader]
  ref_names <- sub("^>", "", vapply(seq_fields, `[`, "", 1L))
  qry_names <- vapply(seq_fields, `[`, "", 2L)
  ref_lens <- as.numeric(vapply(seq_fields, `[`, "", 3L))
  qry_lens <- as.numeric(vapply(seq_fields, `[`, "", 4L))

  # One entry per alignment block, indexed by aln_id
  # Coordinates are 1-based
  coord_fields <- fields[is_coord]
  S1 <- as.numeric(vapply(coord_fields, `[`, "", 1L))
  E1 <- as.numeric(vapply(coord_fields, `[`, "", 2L))
  S2 <- as.numeric(vapply(coord_fields, `[`, "", 3L))
  E2 <- as.numeric(vapply(coord_fields, `[`, "", 4L))
  stopifnot(all(S1 < E1), all(S2 < E2))

  # The sequence pair that each alignment block belongs to
  coord_seq_id <- seq_id[is_coord]

  # Delta values grouped by alignment id, terminating zero dropped
  delta_vals <- as.integer(unlist(fields[is_deltaval]))
  stopifnot(!anyNA(delta_vals))
  delta_aln_id <- aln_id[is_deltaval]
  keep <- delta_vals != 0L
  deltas_per_block <- split(
    delta_vals[keep],
    factor(delta_aln_id[keep], levels = seq_len(sum(is_coord)))
  )

  # Return a list of blocks
  data <- list(
    ref = ref_names[coord_seq_id],
    qry = qry_names[coord_seq_id],
    reflen = ref_lens[coord_seq_id],
    qrylen = qry_lens[coord_seq_id],
    S1 = S1,
    E1 = E1,
    S2 = S2,
    E2 = E2,
    deltas = deltas_per_block
  )
  do.call(mapply, c(FUN = list("list"), data, SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

insert_gaps <- function(chars, at) {
  n <- length(chars)
  if (length(at) == 0L) {
    return(chars)
  }
  gaps_le <- findInterval(seq_len(n), at) # num of `at` values <= k
  final_pos <- seq_len(n) + gaps_le # each char's final position
  result <- rep.int("-", n + length(at))
  result[final_pos] <- chars
  result
}

# Slide all gaps to the right if they don't change the matches in the alignment
normalize_alignment <- function(ref_seq, qry_seq) {}

build_alignment <- function(ref_seq, qry_seq, block) {
  ref_chars <- strsplit(substr(ref_seq, block$S1, block$E1), "", fixed = TRUE)[[1]]
  qry_chars <- strsplit(substr(qry_seq, block$S2, block$E2), "", fixed = TRUE)[[1]]

  d <- block$deltas
  ref_adv <- ifelse(d > 0L, d, -d - 1L) # ref chars consumed by each delta
  qry_adv <- ifelse(d > 0L, d - 1L, -d) # query chars consumed by each delta
  ref_next <- 1L + cumsum(ref_adv) # ref column just after each delta
  qry_next <- 1L + cumsum(qry_adv) # query column just after each delta

  gaps_ref <- ref_next[d < 0L] # gaps go into the ref row
  gaps_qry <- qry_next[d > 0L] # gaps go into the query row

  ref_aln <- insert_gaps(ref_chars, gaps_ref)
  qry_aln <- insert_gaps(qry_chars, gaps_qry)
  stopifnot(length(ref_aln) == length(qry_aln))
}


mk_blocks_main <- function(input_fasta, input_delta, input_csv, output_json) {
  seqs <- read.fasta(input_fasta,
    seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE,
    set.attributes = FALSE
  )

  blocks <- parse_delta(input_delta)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) >= 4)

  mk_blocks_main(args[1], args[2], args[3], args[4])
}
