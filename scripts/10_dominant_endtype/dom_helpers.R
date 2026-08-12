#!/usr/bin/env Rscript
# Stage 17 helpers: fine-grained 5'-end classification and the dominance
# cascade. Pure functions only, so that test/dominance/*.R can source this file
# without touching BAMs or the filesystem.

# ── 5'-end classification ─────────────────────────────────────────────────────

# Length of the 5' soft clip, strand-aware.
#   + strand: leading soft clip  ("2S48M" -> 2)
#   - strand: trailing soft clip ("48M2S" -> 2)
# Carried over unchanged from scripts/05_tss_end_types/01_build_ctss.R, where it
# is already validated.
clip5_len <- function(cigar, is_plus) {
  len <- integer(length(cigar))
  len[ is_plus] <- as.integer(ifelse(
    grepl("^\\d+S", cigar[is_plus]),
    sub("^(\\d+)S.*", "\\1", cigar[is_plus]), "0"))
  len[!is_plus] <- as.integer(ifelse(
    grepl("\\d+S$", cigar[!is_plus]),
    sub(".*?(\\d+)S$", "\\1", cigar[!is_plus]), "0"))
  len
}

# Reverse complement of a character vector, without pulling in Biostrings.
.revcomp <- function(x) {
  vapply(x, function(s) {
    if (nchar(s) == 0) return("")
    ch <- rev(strsplit(s, "")[[1]])
    paste(c(A = "T", T = "A", G = "C", C = "G", N = "N")[ch], collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

# The 5' soft clip in the read's ORIGINAL 5'->3' orientation.
#   + strand: the first clen stored bases, as stored;
#   - strand: the last clen stored bases, reverse-complemented, because the BAM
#     stores the reverse complement of the sequenced read.
clip5_seq <- function(seq_str, clen, is_plus) {
  out <- character(length(seq_str))
  has <- clen > 0
  p <- has &  is_plus
  m <- has & !is_plus
  if (any(p)) out[p] <- substr(seq_str[p], 1, clen[p])
  if (any(m)) {
    n <- nchar(seq_str[m])
    out[m] <- .revcomp(substr(seq_str[m], n - clen[m] + 1, n))
  }
  out
}

# Fine type label. "M" for no clip; "<len>S<seq>" in lower case otherwise;
# clips longer than max_clip collapse to "longclip"; clips containing an
# ambiguous base collapse to "Nclip".
fine_type <- function(clen, clipseq, max_clip) {
  type <- paste0(clen, "S", tolower(clipseq))
  type[clen == 0]                <- "M"
  type[clen > max_clip]          <- "longclip"
  type[grepl("N", clipseq, fixed = TRUE)] <- "Nclip"
  type
}

# ── clustering ────────────────────────────────────────────────────────────────

# distclu: merge positions within max_dist bp on the same chromosome and strand.
# This reproduces scripts/05_tss_end_types/02_cluster_annotate.R lines 125-133
# exactly, including the cluster_id naming, so that the two cluster sets are
# comparable. `total_ctss` must already have the min_total floor applied,
# because Stage 8 applies it BEFORE clustering.
derive_clusters <- function(total_ctss, max_dist) {
  dt <- data.table::copy(total_ctss)
  data.table::setorder(dt, chr, strand, pos)
  dt[, cluster_local := {
    if (.N == 1L) 1L else cumsum(c(Inf, diff(pos)) > max_dist) + 1L
  }, by = .(chr, strand)]
  dt[, cluster_id := paste(chr, strand, cluster_local, sep = ":")]
  dt[, cluster_local := NULL]
  dt[]
}

# ── dominance cascade ─────────────────────────────────────────────────────────

# `counts`: long data.table with cluster_id, sample_id, type, n.
# `bg`: named numeric vector of genome-wide background shares, summing to <= 1.
# Returns one row per cluster. Each pass_* column records that step in
# isolation; `dominated` is the conjunction of all four, so a cluster failing an
# early step can never be rescued by a later one.
call_dominance <- function(counts, bg, min_depth, threshold, min_samples,
                           min_reads_per_sample, fdr) {
  stopifnot(all(c("cluster_id", "sample_id", "type", "n") %in% names(counts)))

  pooled <- counts[, .(n = sum(n)), by = .(cluster_id, type)]
  tot    <- pooled[, .(n_total = sum(n)), by = cluster_id]
  top    <- pooled[pooled[, .I[which.max(n)], by = cluster_id]$V1]
  data.table::setnames(top, c("type", "n"), c("dom_type", "n_dom"))
  res <- merge(top, tot, by = "cluster_id")
  res[, frac := n_dom / n_total]

  res[, pass_depth := n_total >= min_depth]
  res[, pass_frac  := frac >= threshold]

  # One-sided binomial test: is the dominant type's share above its background?
  # Types absent from `bg` (longclip, Nclip, rare clips) get the smallest
  # observed background, which is the conservative choice for a rare type.
  bg_floor <- min(bg)
  res[, bg_p := ifelse(dom_type %in% names(bg), bg[dom_type], bg_floor)]
  res[, p_binom := stats::pbinom(n_dom - 1L, n_total, bg_p, lower.tail = FALSE)]
  res[, q_binom := stats::p.adjust(p_binom, method = "BH")]
  res[, pass_binom := q_binom < fdr]

  # Reproducibility: the same type is the per-sample argmax in at least
  # min_samples samples that have at least min_reads_per_sample reads.
  per <- counts[, .(n = sum(n)), by = .(cluster_id, sample_id, type)]
  dep <- per[, .(depth = sum(n)), by = .(cluster_id, sample_id)]
  per <- merge(per, dep, by = c("cluster_id", "sample_id"))
  per <- per[depth >= min_reads_per_sample]
  arg <- per[per[, .I[which.max(n)], by = .(cluster_id, sample_id)]$V1]
  arg <- merge(arg, res[, .(cluster_id, dom_type)], by = "cluster_id")
  sup <- arg[type == dom_type, .(n_support = .N), by = cluster_id]
  res <- merge(res, sup, by = "cluster_id", all.x = TRUE)
  res[is.na(n_support), n_support := 0L]
  res[, pass_reprod := n_support >= min_samples]

  res[, dominated := pass_depth & pass_frac & pass_binom & pass_reprod]
  res[, bg_p := NULL]
  data.table::setcolorder(res, c("cluster_id", "dom_type", "n_dom", "n_total",
                                 "frac", "p_binom", "q_binom", "n_support",
                                 "pass_depth", "pass_frac", "pass_binom",
                                 "pass_reprod", "dominated"))
  res[]
}

# ── supplementary aggregate-unusual test ─────────────────────────────────────

# `call_dominance()` above takes the argmax over the FINE type vocabulary, so a
# start site whose reads are overwhelmingly unusual but split across sibling
# clips (for example 3Sggt, 3Sgga and 4Sggta at the same position) can never
# reach the fraction threshold on any single label and is therefore invisible to
# that cascade. This supplementary test pools every type other than the two
# usual ones and applies the same statistical filters to the pooled share, so
# that such sites are counted rather than silently suppressed. It does not
# replace the per-type cascade and it does not feed the dominance classes; it is
# reported alongside them.
#
# `counts`: long data.table with cluster_id, sample_id, type, n.
# `bg`: named numeric vector of genome-wide background shares; the pooled
# background is the sum of the shares of every type outside `usual`.
# Reproducibility is defined here as the pooled unusual share reaching the same
# threshold in at least `min_samples` samples, which is the natural analogue of
# the per-sample argmax used by the per-type cascade.
call_aggregate_unusual <- function(counts, bg, min_depth, threshold,
                                   min_samples, min_reads_per_sample, fdr,
                                   usual = c("1Sg", "M")) {
  stopifnot(all(c("cluster_id", "sample_id", "type", "n") %in% names(counts)))

  res <- counts[, .(n_unusual = sum(n[!type %in% usual]), n_total = sum(n)),
                by = cluster_id]
  res[, frac_unusual := n_unusual / n_total]

  res[, pass_depth := n_total >= min_depth]
  res[, pass_frac  := frac_unusual >= threshold]

  bg_unusual <- sum(bg[!names(bg) %in% usual])
  res[, p_binom := stats::pbinom(n_unusual - 1L, n_total, bg_unusual,
                                 lower.tail = FALSE)]
  res[, q_binom := stats::p.adjust(p_binom, method = "BH")]
  res[, pass_binom := q_binom < fdr]

  per <- counts[, .(n_unusual = sum(n[!type %in% usual]), depth = sum(n)),
                by = .(cluster_id, sample_id)]
  per <- per[depth >= min_reads_per_sample]
  sup <- per[n_unusual / depth >= threshold, .(n_support = .N), by = cluster_id]
  res <- merge(res, sup, by = "cluster_id", all.x = TRUE)
  res[is.na(n_support), n_support := 0L]
  res[, pass_reprod := n_support >= min_samples]

  res[, aggregate_unusual := pass_depth & pass_frac & pass_binom & pass_reprod]
  data.table::setcolorder(res, c("cluster_id", "n_unusual", "n_total",
                                 "frac_unusual", "p_binom", "q_binom",
                                 "n_support", "pass_depth", "pass_frac",
                                 "pass_binom", "pass_reprod",
                                 "aggregate_unusual"))
  res[]
}
