# Shared definition of the cluster-defined promoter, sourced by
# 12_promoter_clusters.R (which writes the table the localization stages read)
# and 13_promoter_window_sweep.R (which compares window settings and draws the
# definition), so that the two can never disagree on what a promoter is.
#
# The definition, in gene orientation:
#   1. a window of w bases ends strictly upstream of the annotated gene start,
#      so the start base itself is outside it, and only TSS clusters on the
#      gene's own strand are considered;
#   2. the windows are tried in increasing order (e.g. 15, 20, 25, 30), and the
#      first one that a cluster overlaps decides: every cluster overlapping that
#      window is a promoter cluster of the gene, and a cluster reachable only by
#      a longer window is not. A gene that no window reaches has no promoter;
#   3. the promoter of a gene-cluster pair runs from the cluster's upstream edge
#      to the base just before the gene start, so it covers the gap between the
#      cluster and the start; a cluster that runs past the start keeps its full
#      extent.
# A single fixed window is the special case of a one-element windows vector.
#
# Requires data.table and GenomicRanges, and derive_clusters() from
# scripts/10_dominant_endtype/dom_helpers.R.

# Re-derives the Stage 8 cluster extents, which Stage 8 does not store, and
# asserts that they reproduce Stage 8's cluster set exactly.
load_cluster_extents <- function(ctss_f, mat_f, min_total, max_dist) {
  ctss   <- data.table::fread(ctss_f, select = c("chr", "pos", "strand", "n_total"))
  pooled <- ctss[, .(n_total = sum(n_total)), by = .(chr, pos, strand)]
  members <- derive_clusters(pooled[n_total >= min_total], max_dist)

  stage8 <- unique(data.table::fread(mat_f, select = c("cluster_id", "chr", "tss_pos", "strand")))
  stopifnot(!anyDuplicated(stage8$cluster_id))
  if (!setequal(unique(members$cluster_id), stage8$cluster_id)) {
    stop("Re-derived cluster set differs from ", mat_f, ": ",
         data.table::uniqueN(members$cluster_id), " here vs ", nrow(stage8), " there. ",
         "Check that a8_min_total and a8_max_dist equal Stage 8's settings.")
  }
  extent <- members[, .(start = min(pos), end = max(pos), reads = sum(n_total)),
                    by = cluster_id]
  cl <- merge(stage8, extent, by = "cluster_id")
  stopifnot(nrow(cl) == nrow(stage8), all(cl$tss_pos >= cl$start & cl$tss_pos <= cl$end))
  list(clusters = cl, members = members)
}

# Reads a BED6 gene file and adds the 1-based gene start. BED is 0-based
# half-open, so the start of a + gene is start + 1 and that of a - gene is end.
load_gene_starts <- function(genes_f) {
  genes <- data.table::fread(genes_f, header = FALSE, select = 1:6,
                             col.names = c("chr", "start", "end", "gene", "score", "strand"))
  genes[, gene_start := ifelse(strand == "+", start + 1L, end)]
  genes[]
}

# Every same-strand gene-cluster pair within max(windows) bases upstream, with
# `gap`, the distance in bases from the gene start to the nearest base of the
# cluster that lies upstream of it (1 = the base just before the start). A
# cluster overlaps the window of length w exactly when gap <= w.
candidate_pairs <- function(cl, genes, max_window) {
  window_gr <- GRanges(genes$chr,
                       IRanges(pmax(1L, ifelse(genes$strand == "+",
                                               genes$gene_start - max_window,
                                               genes$gene_start + 1L)),
                               ifelse(genes$strand == "+",
                                      genes$gene_start - 1L,
                                      genes$gene_start + max_window)),
                       strand = genes$strand)
  cl_gr <- GRanges(cl$chr, IRanges(cl$start, cl$end), strand = cl$strand)
  # findOverlaps respects strand by default, so antisense pairs never appear.
  hits <- findOverlaps(cl_gr, window_gr)
  pairs <- cbind(cl[queryHits(hits), .(cluster_id, chr, cluster_start = start,
                                       cluster_end = end, strand, tss_pos, reads)],
                 genes[subjectHits(hits), .(gene, gene_start)])
  pairs[, gap := ifelse(strand == "+",
                        gene_start - pmin(cluster_end, gene_start - 1L),
                        pmax(cluster_start, gene_start + 1L) - gene_start)]
  stopifnot(all(pairs$gap >= 1L & pairs$gap <= max_window))
  pairs[]
}

# Applies the growing-window rule to the candidate pairs and returns the
# promoter pairs, with the deciding `window` and the promoter span in `start`
# and `end` (1-based, inclusive), which are the columns load_promoter_clusters()
# in scripts/utils/localization.R reads.
define_promoters <- function(pairs, windows) {
  windows <- sort(unique(as.integer(windows)))
  if (length(windows) == 0L || anyNA(windows) || any(windows < 1L))
    stop("promoter windows must be positive integers")
  p <- data.table::copy(pairs[gap <= max(windows)])
  # The deciding window of a gene is the shortest one its nearest cluster
  # reaches; findInterval with left.open picks the first window >= min(gap).
  p[, window := windows[findInterval(min(gap), windows, left.open = TRUE) + 1L],
    by = gene]
  p <- p[gap <= window]
  p[, `:=`(start = ifelse(strand == "+", cluster_start, pmin(cluster_start, gene_start + 1L)),
           end   = ifelse(strand == "+", pmax(cluster_end, gene_start - 1L), cluster_end))]
  p[, dist_to_gene_start := ifelse(strand == "+", gene_start - tss_pos,
                                   tss_pos - gene_start)]
  data.table::setorder(p, chr, start, strand, gene)
  p[]
}
