#!/usr/bin/env Rscript
# Stage 8 step 12: define promoters from the data, as the distclu TSS clusters
# that lie just upstream of a gene.
#
# The definition lives in promoter_helpers.R: windows of the lengths listed in
# promoter_cluster_windows end strictly upstream of each gene's annotated start
# and are tried shortest first; the clusters on the gene's strand that overlap
# the first window reaching any cluster are its promoter clusters, and the
# promoter runs from such a cluster's upstream edge to the gene start. The
# localization stages (10 and 19) label a position "Promoter" from this table
# when promoter_source is "clusters" in the params file -- see
# scripts/utils/localization.R. 13_promoter_window_sweep.R compares window
# settings and draws the definition.
#
# Stage 8 stores only the dominant TSS of a cluster, so the extents are
# re-derived with derive_clusters() and asserted equal to Stage 8's cluster
# set, exactly as 11_export_cluster_bed.R does. Every Stage 8 cluster is
# eligible, not only those Stage 9 kept for being present in two samples.
#
# promoter_clusters.tsv has one row per gene-cluster pair: cluster_start and
# cluster_end are the cluster's extent, start and end are the promoter span,
# window is the window that decided the gene, and gap is the distance from the
# gene start to the cluster. dist_to_gene_start is measured from the dominant
# TSS and is positive upstream of the gene start; it is zero or negative for
# clusters whose dominant TSS sits at or past the start.
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/12_promoter_clusters.R \
#       [--ctss     results/tss/typed_ctss.tsv.gz] \
#       [--matrix   results/tss/tss_matrix.tsv]    \
#       [--genes    data/sacCer3.genes.bed]        \
#       [--outdir   results/tss]                   \
#       [--trackdir results/tracks]                \
#       [--params   config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag)
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))
source(file.path(script_dir, "promoter_helpers.R"))

params   <- yaml.load_file(get_arg("--params", "config/params.yaml"))
ctss_f   <- get_arg("--ctss",     "results/tss/typed_ctss.tsv.gz")
mat_f    <- get_arg("--matrix",   "results/tss/tss_matrix.tsv")
genes_f  <- get_arg("--genes",    params$annotation_bed)
outdir   <- get_arg("--outdir",   "results/tss")
trackdir <- get_arg("--trackdir", "results/tracks")
force    <- has_flag("--force")

# The Stage 8 floor and distance, read from where 11_export_cluster_bed.R reads
# them.
min_total <- as.integer(params$a8_min_total)
max_dist  <- as.integer(params$a8_max_dist)
windows   <- sort(as.integer(unlist(params$promoter_cluster_windows)))
if (length(windows) == 0L || anyNA(windows) || any(windows < 1L))
  stop("promoter_cluster_windows must be a list of positive integers in the params file")

out_tsv <- file.path(outdir,   "promoter_clusters.tsv")
out_bed <- file.path(trackdir, "promoter_clusters.bed")
if (!force && file.exists(out_tsv) && file.exists(out_bed)) {
  message("Outputs exist, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}
dir.create(outdir,   recursive = TRUE, showWarnings = FALSE)
dir.create(trackdir, recursive = TRUE, showWarnings = FALSE)

cl <- load_cluster_extents(ctss_f, mat_f, min_total, max_dist)$clusters
message("Cluster set matches Stage 8: ", nrow(cl), " clusters")
genes <- load_gene_starts(genes_f)

res <- define_promoters(candidate_pairs(cl, genes, max(windows)), windows)
fwrite(res[, .(cluster_id, chr, cluster_start, cluster_end, strand, tss_pos, reads,
               gene, gene_start, window, gap, start, end, dist_to_gene_start)],
       out_tsv, sep = "\t")

n_cl <- uniqueN(res$cluster_id)
message("Windows tried: ", paste(windows, collapse = ", "), " bp")
message("Promoter clusters: ", n_cl, " of ", nrow(cl), " (",
        round(100 * sum(cl$reads[cl$cluster_id %in% res$cluster_id]) / sum(cl$reads), 1),
        " % of clustered reads)")
message("Genes with a promoter: ", uniqueN(res$gene), " of ", nrow(genes),
        "; decided by window: ",
        paste(res[, .(w = window[1]), by = gene][, .N, by = w][order(w)][
          , sprintf("%d bp %d", w, N)], collapse = ", "))
message("Clusters assigned to more than one gene: ",
        sum(res[, .N, by = cluster_id]$N > 1L))
message("Written: ", out_tsv)

# ── track ─────────────────────────────────────────────────────────────────────
# One box per promoter span, named by its gene or genes; the score is the pooled
# depth of the cluster capped at 1000 as BED requires.
bed <- res[, .(gene = paste(gene, collapse = ",")),
           by = .(chr, start, end, strand, reads)]
setorder(bed, chr, start, strand)
con <- file(out_bed, "w")
writeLines(sprintf('track name="Promoter clusters" description="TSS cluster to gene start, windows %s bp" visibility=pack',
                   paste(windows, collapse = "/")), con)
close(con)
fwrite(bed[, .(chr, start - 1L, end, gene, pmin(1000L, reads), strand)], out_bed,
       sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
message("Written: ", out_bed, " (", nrow(bed), " promoters)")

sessionInfo()
