#!/usr/bin/env Rscript
# Stage 8 step 11: export the distclu TSS clusters as genome-browser tracks.
#
# Stage 8 keeps one row per cluster and sample (results/tss/tss_matrix.tsv)
# carrying only the dominant TSS, so the cluster's own extent is not stored
# anywhere. It is recovered here by re-running the same distclu rule with
# derive_clusters() from scripts/10_dominant_endtype/dom_helpers.R, which
# reproduces Stage 8's cluster ids exactly; the resulting cluster set is
# asserted equal to Stage 8's before anything is written, so a mismatch in
# min_total or max_dist fails loudly rather than drawing shifted boxes.
#
# The two tracks cover different cluster sets. Stage 9 groups only the clusters
# present in at least two samples (results/tss/tss_clustered.tsv), so the
# hclust-group track is that subset, while Stage 17 calls dominance for every
# Stage 8 cluster, so the dominant-type track draws all of them.
#
# Two BED9 tracks come out, identical except for their colour and name:
#
#   tss_clusters_by_hclust_group.bed   coloured by the Stage 9 hclust_group,
#                                      i.e. by the cluster's end-type profile;
#   tss_clusters_by_dominant_type.bed  coloured by the Stage 17 dominance call
#                                      at the primary threshold, with every
#                                      cluster that is not dominated in grey.
#
# The score is the cluster's pooled depth divided by the number of samples,
# capped at 1000 as BED requires.
#
# Each box spans the cluster's first to last CTSS, and the thick part marks the
# dominant TSS. A small legend TSV is written beside each track, because a BED file
# cannot carry one.
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/11_export_cluster_bed.R \
#       [--ctss      results/tss/typed_ctss.tsv.gz]   \
#       [--matrix    results/tss/tss_matrix.tsv]      \
#       [--clustered results/tss/tss_clustered.tsv]   \
#       [--dominance results/dominance/dominance_calls.tsv] \
#       [--outdir    results/tracks]                  \
#       [--params    config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(yaml)
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
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))

ctss_f  <- get_arg("--ctss",      "results/tss/typed_ctss.tsv.gz")
mat_f   <- get_arg("--matrix",    "results/tss/tss_matrix.tsv")
clus_f  <- get_arg("--clustered", "results/tss/tss_clustered.tsv")
dom_f   <- get_arg("--dominance", "results/dominance/dominance_calls.tsv")
outdir  <- get_arg("--outdir",    "results/tracks")
params  <- yaml.load_file(get_arg("--params", "config/params.yaml"))
force   <- has_flag("--force")

# The Stage 8 floor and distance. The A-8 stage keeps its own copies of both,
# documented as having to equal Stage 8's, so they are read from there.
min_total <- as.integer(params$a8_min_total)
max_dist  <- as.integer(params$a8_max_dist)
dom_thr   <- as.numeric(params$dom_thresholds[1])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_grp <- file.path(outdir, "tss_clusters_by_hclust_group.bed")
out_dom <- file.path(outdir, "tss_clusters_by_dominant_type.bed")

if (!force && file.exists(out_grp) && file.exists(out_dom)) {
  message("Outputs exist, skipping (use --force to rerun): ", out_grp)
  quit(save = "no", status = 0)
}

# ── cluster extents ───────────────────────────────────────────────────────────
ctss   <- fread(ctss_f, select = c("chr", "pos", "strand", "sample_id", "n_total"))
n_samples <- uniqueN(ctss$sample_id)
pooled <- ctss[, .(n_total = sum(n_total)), by = .(chr, pos, strand)]
members <- derive_clusters(pooled[n_total >= min_total], max_dist)

stage8 <- unique(fread(mat_f, select = c("cluster_id", "chr", "tss_pos", "strand",
                                         "initiator")))
stopifnot(!anyDuplicated(stage8$cluster_id))
if (!setequal(unique(members$cluster_id), stage8$cluster_id)) {
  stop("Re-derived cluster set differs from ", mat_f, ": ",
       uniqueN(members$cluster_id), " here vs ", nrow(stage8), " there. ",
       "Check that a8_min_total and a8_max_dist equal Stage 8's settings.")
}
message("Cluster set matches Stage 8: ", nrow(stage8), " clusters")

extent <- members[, .(start = min(pos), end = max(pos),
                      mean_depth = sum(n_total) / n_samples), by = cluster_id]
bed <- merge(stage8, extent, by = "cluster_id")
stopifnot(nrow(bed) == nrow(stage8),
          all(bed$tss_pos >= bed$start & bed$tss_pos <= bed$end))

clus <- fread(clus_f)
stopifnot(all(clus$cluster_id %in% bed$cluster_id))
bed <- merge(bed, clus[, .(cluster_id, hclust_group)], by = "cluster_id", all.x = TRUE)
message("Clusters with an hclust group: ", sum(!is.na(bed$hclust_group)),
        " of ", nrow(bed))

hex_to_rgb <- function(hex) apply(grDevices::col2rgb(hex), 2, paste, collapse = ",")

# ── colour 1: hclust group ────────────────────────────────────────────────────
# Each group is labelled by its two most abundant end types and their mean
# shares, taken from the table itself, so a re-run of Stage 9 that reorders or
# redefines the groups relabels them rather than leaving stale names.
pct_cols <- grep("^mean_pct_", names(clus), value = TRUE)
profile  <- clus[, lapply(.SD, mean), by = hclust_group, .SDcols = pct_cols]
setorder(profile, hclust_group)
profile[, label := apply(.SD, 1, function(r) {
  top <- order(-r)[1:2]
  paste(sprintf("%s %.0f%%", sub("^mean_pct_", "", pct_cols[top]), r[top]),
        collapse = ", ")
}), .SDcols = pct_cols]

# Okabe-Ito, which stays distinguishable under the common colour-vision
# deficiencies; recycled only if Stage 9 ever cuts more than eight groups.
GROUP_PAL <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00",
               "#CC79A7", "#D55E00", "#F0E442", "#000000")
profile[, colour := GROUP_PAL[(seq_len(.N) - 1L) %% length(GROUP_PAL) + 1L]]
bed <- merge(bed, profile[, .(hclust_group, colour_grp = colour)],
             by = "hclust_group", all.x = TRUE)

# ── colour 2: dominant end type ───────────────────────────────────────────────
calls <- fread(dom_f)[threshold == dom_thr, .(cluster_id, dom_type, frac, dominated)]
bed <- merge(bed, calls, by = "cluster_id", all.x = TRUE)
bed[is.na(dominated), dominated := FALSE]

NOT_DOMINATED  <- "not dominated"
OFF_VOCABULARY <- "other clip"
bed[, dom_label := fifelse(!dominated, NOT_DOMINATED,
                   fifelse(dom_type %in% END_TYPES_REAL, dom_type, OFF_VOCABULARY))]

dom_pal <- c(END_TYPE_COLORS_BRIGHT[END_TYPES_REAL], "#495057", "#c8c8c8")
names(dom_pal) <- c(END_TYPES_REAL, OFF_VOCABULARY, NOT_DOMINATED)
bed[, colour_dom := dom_pal[dom_label]]

# ── write ─────────────────────────────────────────────────────────────────────
# BED is 0-based half-open: a 1-based CTSS at p covers [p - 1, p).
setorder(bed, chr, start, strand)

write_track <- function(path, name, description, rows, item_name, colour) {
  b <- bed[rows]
  out <- data.table(chrom = b$chr, chromStart = b$start - 1L, chromEnd = b$end,
                    name = item_name[rows],
                    score = pmin(1000L, as.integer(round(b$mean_depth))),
                    strand = b$strand,
                    thickStart = b$tss_pos - 1L, thickEnd = b$tss_pos,
                    itemRgb = hex_to_rgb(colour[rows]))
  con <- file(path, "w")
  writeLines(sprintf('track name="%s" description="%s" itemRgb="On" visibility=full',
                     name, description), con)
  close(con)
  fwrite(out, path, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
  message("Written: ", path, " (", nrow(out), " clusters)")
}

write_track(out_grp, "TSS clusters by hclust group",
            "distclu TSS clusters coloured by end-type profile group (Stage 9)",
            which(!is.na(bed$hclust_group)),
            paste(bed$hclust_group, bed$initiator, sep = "|"), bed$colour_grp)

write_track(out_dom, "TSS clusters by dominant type",
            paste0("distclu TSS clusters coloured by dominant 5' end type at ",
                   dom_thr, " (Stage 17); grey = not dominated"),
            seq_len(nrow(bed)),
            # UCSC splits a custom-track line on any whitespace, so an item name
            # must not contain a space.
            ifelse(bed$dominated,
                   sprintf("%s_%.0f%%|%s", bed$dom_type, 100 * bed$frac, bed$initiator),
                   paste(gsub(" ", "_", NOT_DOMINATED), bed$initiator, sep = "|")),
            bed$colour_dom)

legend_grp <- merge(profile[, .(hclust_group, label, colour)],
                    bed[!is.na(hclust_group), .(n_clusters = .N), by = hclust_group],
                    by = "hclust_group")
fwrite(legend_grp, file.path(outdir, "tss_clusters_by_hclust_group_legend.tsv"), sep = "\t")

legend_dom <- bed[, .(n_clusters = .N), by = .(dom_label, colour = colour_dom)][order(-n_clusters)]
fwrite(legend_dom, file.path(outdir, "tss_clusters_by_dominant_type_legend.tsv"), sep = "\t")

print(legend_grp)
print(legend_dom)

sessionInfo()
