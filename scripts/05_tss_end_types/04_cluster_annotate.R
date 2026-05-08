#!/usr/bin/env Rscript
# Stage 9 step 1: hierarchical clustering of TSS by end-type vector +
# strand-aware BED annotation (overlap → sense/antisense genic; else intergenic).
#
# Clustering: average pct_1Sg/2Sg/M/other across all samples per TSS →
# Euclidean distance → Ward.D2 linkage → cut at --k_clusters (default 6).
#
# Annotation: for each dominant TSS position, intersect with the gene annotation
# BED (strand-aware). TSS overlapping a gene on the same strand → "sense_genic";
# opposite strand → "antisense_genic"; no overlap → "intergenic" (nearest gene
# distance reported).
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/04_cluster_annotate.R \
#       [--matrix     results/tss/tss_matrix.tsv]    \
#       [--samples    config/samples.tsv]             \
#       [--outdir     results/tss]                    \
#       [--params     config/params.yaml]             \
#       [--annot_bed  /path/to/annotation.bed]        \
#       [--k_clusters 6]                              \
#       [--min_samples 2]                             \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(readr)
  library(yaml)
})

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

matrix_f    <- get_arg("--matrix",     "results/tss/tss_matrix.tsv")
samples_f   <- get_arg("--samples",    "config/samples.tsv")
outdir      <- get_arg("--outdir",     "results/tss")
params_f    <- get_arg("--params",     "config/params.yaml")
annot_bed   <- get_arg("--annot_bed",  NULL)
k_clusters  <- as.integer(get_arg("--k_clusters",  "6"))
min_samples <- as.integer(get_arg("--min_samples", "2"))
force       <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_tsv <- file.path(outdir, "tss_clustered.tsv")
out_rds <- file.path(outdir, "tss_hclust.rds")

if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

params <- yaml.load_file(params_f)

# Resolve annotation BED: CLI arg overrides params.yaml
if (is.null(annot_bed) || annot_bed == "") {
  annot_bed <- params$annotation_bed_yeast
}

# ── load data ─────────────────────────────────────────────────────────────────
message("Stage 9 step 1: cluster + annotate")
mat <- fread(matrix_f)
message("  clusters: ", uniqueN(mat$cluster_id),
        "  samples: ", uniqueN(mat$sample_id))

# ── build per-TSS average end-type vector ─────────────────────────────────────
# Filter: keep TSS present (n_total > 0) in at least min_samples samples
present_counts <- mat[n_total > 0, .(n_present = .N), by = cluster_id]
keep_ids <- present_counts[n_present >= min_samples, cluster_id]
mat_f <- mat[cluster_id %in% keep_ids]

message("  clusters with >= ", min_samples, " samples: ", length(keep_ids))

# Average pct values across all samples (ignoring NA / zero-coverage)
avg <- mat_f[n_total > 0, .(
  mean_pct_1Sg   = mean(pct_1Sg,   na.rm = TRUE),
  mean_pct_2Sg   = mean(pct_2Sg,   na.rm = TRUE),
  mean_pct_M     = mean(pct_M,     na.rm = TRUE),
  mean_pct_other = mean(pct_other, na.rm = TRUE),
  mean_total     = mean(n_total,   na.rm = TRUE)
), by = .(cluster_id, chr, tss_pos, strand, initiator)]

# ── hierarchical clustering ───────────────────────────────────────────────────
# Euclidean distance on 4-component end-type vector; Ward.D2 linkage.
feat_cols <- c("mean_pct_1Sg", "mean_pct_2Sg", "mean_pct_M", "mean_pct_other")
feat_mat  <- as.matrix(avg[, ..feat_cols])
rownames(feat_mat) <- avg$cluster_id

# Remove rows with any NA (can arise if a cluster has 0 coverage in all samples)
complete_rows <- complete.cases(feat_mat)
if (any(!complete_rows)) {
  message("  dropping ", sum(!complete_rows), " clusters with incomplete pct vectors")
  feat_mat <- feat_mat[complete_rows, , drop = FALSE]
  avg      <- avg[cluster_id %in% rownames(feat_mat)]
}

message("  clustering ", nrow(feat_mat), " TSS (Euclidean + Ward.D2)")
d   <- dist(feat_mat, method = "euclidean")
hc  <- hclust(d, method = "ward.D2")

# Save hclust object for dendrogram plotting in next script
saveRDS(hc, out_rds)
message("  hclust object saved: ", out_rds)

# Cut tree
hclust_group <- cutree(hc, k = k_clusters)
avg$hclust_group <- paste0("C", hclust_group[avg$cluster_id])

# ── BED annotation ────────────────────────────────────────────────────────────
if (!is.null(annot_bed) && annot_bed != "" && file.exists(annot_bed)) {
  message("  annotating against: ", annot_bed)

  genes <- import(annot_bed, format = "BED")

  tss_gr <- GRanges(
    seqnames = avg$chr,
    ranges   = IRanges(start = avg$tss_pos, width = 1),
    strand   = avg$strand
  )
  names(tss_gr) <- avg$cluster_id

  # Same-strand overlap
  hits_sense <- findOverlaps(tss_gr, genes, ignore.strand = FALSE)
  # Any-strand overlap (to catch antisense)
  hits_any   <- findOverlaps(tss_gr, genes, ignore.strand = TRUE)
  hits_anti  <- hits_any[!(queryHits(hits_any) %in% queryHits(hits_sense))]

  annotation <- rep("intergenic", length(tss_gr))
  annotation[queryHits(hits_anti)]  <- "antisense_genic"
  annotation[queryHits(hits_sense)] <- "sense_genic"  # sense wins over antisense

  # Nearest gene distance for intergenic TSS
  intergenic_idx <- which(annotation == "intergenic")
  dist_nearest   <- rep(NA_integer_, length(tss_gr))
  nearest_gene   <- rep(NA_character_, length(tss_gr))

  if (length(intergenic_idx) > 0) {
    nn <- distanceToNearest(tss_gr[intergenic_idx], genes, ignore.strand = TRUE)
    dist_nearest[intergenic_idx[queryHits(nn)]] <- mcols(nn)$distance
    if (!is.null(genes$name)) {
      nearest_gene[intergenic_idx[queryHits(nn)]] <-
        genes$name[subjectHits(nn)]
    }
  }

  avg$annotation    <- annotation
  avg$dist_nearest  <- dist_nearest
  avg$nearest_gene  <- nearest_gene

} else {
  if (!is.null(annot_bed) && annot_bed != "") {
    warning("Annotation BED not found: ", annot_bed, " — skipping annotation")
  } else {
    message("  no annotation BED provided (set annotation_bed_yeast in params.yaml)")
  }
  avg$annotation   <- NA_character_
  avg$dist_nearest <- NA_integer_
  avg$nearest_gene <- NA_character_
}

# ── write output ──────────────────────────────────────────────────────────────
setDT(avg)
fwrite(avg, out_tsv, sep = "\t")
message("Written: ", out_tsv, "  (", nrow(avg), " TSS, ", k_clusters, " hclust groups)")
message("  group sizes:")
print(avg[, .N, by = hclust_group][order(hclust_group)])

sessionInfo()
