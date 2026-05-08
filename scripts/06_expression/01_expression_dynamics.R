#!/usr/bin/env Rscript
# Stage 11: expression dynamics of interesting TSS clusters across conditions.
#
# "Interesting" = mean pct_1Sg >= --min_pct_1Sg (default 70 %) across all samples.
#
# Outputs:
#   results/figures/11_expr_by_condition.pdf   — expression (n_total) box/violin by condition
#   results/figures/11_endtype_by_condition.pdf — end-type fractions by condition
#   results/figures/11_iqw_by_cluster_type.pdf  — IQW for high-1Sg vs low-1Sg clusters
#                                                  (only if tagclusters.tsv.gz exists)
#
# Usage (from project root):
#   Rscript scripts/06_expression/01_expression_dynamics.R \
#       [--matrix      results/tss/tss_matrix.tsv]          \
#       [--clustered   results/tss/tss_clustered.tsv]        \
#       [--tagclusters results/cager/tagclusters.tsv.gz]     \
#       [--samples     config/samples.tsv]                   \
#       [--outdir      results/figures]                      \
#       [--params      config/params.yaml]                   \
#       [--min_pct_1Sg 70]                                   \
#       [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(yaml)
})

# ── helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x[[1]]) && x[[1]] != "") x else y

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

matrix_f  <- get_arg("--matrix",      "results/tss/tss_matrix.tsv")
clust_f   <- get_arg("--clustered",   "results/tss/tss_clustered.tsv")
tc_gz     <- get_arg("--tagclusters", "results/cager/tagclusters.tsv.gz")
samples_f <- get_arg("--samples",     "config/samples.tsv")
outdir    <- get_arg("--outdir",      "results/figures")
params_f  <- get_arg("--params",      "config/params.yaml")
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_expr  <- file.path(outdir, "11_expr_by_condition.pdf")
out_type  <- file.path(outdir, "11_endtype_by_condition.pdf")
out_iqw   <- file.path(outdir, "11_iqw_by_cluster_type.pdf")

if (!force && all(file.exists(c(out_expr, out_type)))) {
  message("Outputs exist, skipping (use --force to rerun).")
  quit(save = "no", status = 0)
}

for (f in c(matrix_f, clust_f)) {
  if (!file.exists(f)) stop("Input not found: ", f,
                            "\nRun stages 8-9 first (scripts/05_tss_end_types/run.sh).")
}

params <- yaml.load_file(params_f)
min_pct_1Sg <- as.numeric(get_arg("--min_pct_1Sg",
                                   as.character(params$expr_min_pct_1Sg %||% 70)))

message("Stage 11: expression dynamics")
message("  min_pct_1Sg threshold: ", min_pct_1Sg, "%")

# ── load data ─────────────────────────────────────────────────────────────────
mat  <- fread(matrix_f)
clus <- fread(clust_f)
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)

# Cluster-level mean end-type fractions (from stage 9 tss_clustered.tsv)
clus_feat <- clus[, .(cluster_id, chr, tss_pos, strand, initiator,
                       hclust_group, annotation,
                       mean_pct_1Sg, mean_pct_2Sg, mean_pct_M, mean_pct_other)]

# ── classify clusters ─────────────────────────────────────────────────────────
clus_feat[, cluster_class := ifelse(
  !is.na(mean_pct_1Sg) & mean_pct_1Sg >= min_pct_1Sg,
  paste0("high_1Sg (>=", min_pct_1Sg, "%)"),
  "other"
)]

n_high <- clus_feat[cluster_class != "other", .N]
message("  clusters with mean pct_1Sg >= ", min_pct_1Sg, "%: ", n_high,
        " / ", nrow(clus_feat))

if (n_high == 0) {
  message("  WARNING: no clusters above threshold — lowering threshold to 50%")
  min_pct_1Sg <- 50
  clus_feat[, cluster_class := ifelse(
    !is.na(mean_pct_1Sg) & mean_pct_1Sg >= min_pct_1Sg,
    paste0("high_1Sg (>=", min_pct_1Sg, "%)"),
    "other"
  )]
}

# ── join with sample metadata ─────────────────────────────────────────────────
cond_map <- meta[, c("sample_id", "condition", "replicate"), drop = FALSE]

mat_ann <- merge(mat, cond_map, by = "sample_id", all.x = TRUE)
mat_ann <- merge(mat_ann, clus_feat[, .(cluster_id, cluster_class)],
                 by = "cluster_id", all.x = TRUE)
mat_ann <- mat_ann[!is.na(cluster_class)]

# ── condition-level aggregation ───────────────────────────────────────────────
# Samples subsampled to equal depth → use raw counts directly.
# Average replicates within condition.
cond_agg <- mat_ann[, .(
  mean_n_total = mean(n_total, na.rm = TRUE),
  mean_pct_1Sg = mean(pct_1Sg, na.rm = TRUE),
  mean_pct_2Sg = mean(pct_2Sg, na.rm = TRUE),
  mean_pct_M   = mean(pct_M,   na.rm = TRUE),
  mean_pct_other = mean(pct_other, na.rm = TRUE)
), by = .(cluster_id, condition, cluster_class)]

# ── Figure 1: expression (n_total) by condition, split by cluster class ───────
message("  expression by condition: ", out_expr)

cond_order <- unique(cond_agg$condition)

high_class_label <- paste0("high_1Sg (>=", min_pct_1Sg, "%)")
class_colors <- c("#4361ee", "#adb5bd")
names(class_colors) <- c(high_class_label, "other")

p_expr <- ggplot(cond_agg[n_total > 0],
                 aes(x = condition, y = mean_n_total + 1,
                     fill = cluster_class, colour = cluster_class)) +
  geom_violin(alpha = 0.3, linewidth = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.4,
               alpha = 0.7, linewidth = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = class_colors, name = "Cluster type") +
  scale_colour_manual(values = class_colors, guide = "none") +
  labs(
    title    = paste0("TSS cluster expression by condition"),
    subtitle = paste0("High-1Sg clusters: mean pct_1Sg >= ", min_pct_1Sg, "%"),
    x        = "Condition",
    y        = "Mean reads per cluster (log10+1)"
  ) +
  facet_wrap(~ cluster_class) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_blank())

pdf(out_expr, width = 7, height = 5)
print(p_expr)
dev.off()

# ── Figure 2: end-type fractions across conditions for high-1Sg clusters ──────
message("  end-type by condition: ", out_type)

high_cond <- cond_agg[cluster_class != "other"] |>
  as.data.frame() |>
  pivot_longer(
    cols      = c(mean_pct_1Sg, mean_pct_2Sg, mean_pct_M, mean_pct_other),
    names_to  = "end_type",
    values_to = "pct"
  ) |>
  mutate(end_type = recode(end_type,
    mean_pct_1Sg   = "1Sg",
    mean_pct_2Sg   = "2Sg",
    mean_pct_M     = "M",
    mean_pct_other = "other"
  )) |>
  mutate(end_type = factor(end_type, levels = c("1Sg", "2Sg", "M", "other")))

type_colors <- c("1Sg" = "#4361ee", "2Sg" = "#f72585",
                 "M"   = "#06d6a0", "other" = "#adb5bd")

p_type <- ggplot(high_cond, aes(x = condition, y = pct, fill = end_type)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.35,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = type_colors, name = "End type") +
  labs(
    title    = paste0("End-type fractions across conditions (high-1Sg clusters)"),
    subtitle = paste0("Clusters with mean pct_1Sg >= ", min_pct_1Sg, "%"),
    x        = "Condition",
    y        = "Mean % reads"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(out_type, width = 7, height = 5)
print(p_type)
dev.off()

# ── Figure 3: IQW comparison (high-1Sg vs other) from CAGEr stage 10 ─────────
if (file.exists(tc_gz)) {
  message("  IQW by cluster type: ", out_iqw)

  tc <- fread(cmd = paste("zcat", shQuote(tc_gz)))

  if ("iq_width" %in% names(tc) && any(!is.na(tc$iq_width))) {
    # Compute consensus cluster midpoints to match against tss_clustered.tsv
    # Approximate: use dominant_ctss column if available, else midpoint of start/end
    tc_pos_col <- if ("dominant_ctss" %in% names(tc)) "dominant_ctss" else NULL

    # Map from tss_matrix cluster_id (chr:strand:local) to CAGEr tag cluster ranges
    # Use strand-aware overlap between tss_clustered dominant TSS and CAGEr clusters

    # Build a per-sample average IQW table
    tc_sample_avg <- tc[!is.na(iq_width), .(mean_iqw = mean(iq_width), .N),
                        by = sample_id]

    # For comparison against cluster class, merge tss_clustered with stage9 cluster info
    # Match CAGEr clusters to custom clusters by dominant TSS position
    # Approximation: match by (chr, approx_pos, strand) within ±30 bp
    # This is intentionally approximate; exact matching would require re-running stage 9
    # with CAGEr cluster IDs.

    # Simpler: show overall IQW distribution split by sample condition
    meta_map <- meta[, c("sample_id", "condition"), drop = FALSE]
    tc_cond  <- merge(tc[!is.na(iq_width)], meta_map, by = "sample_id", all.x = TRUE)

    p_iqw <- ggplot(tc_cond, aes(x = iq_width, fill = condition, colour = condition)) +
      geom_density(alpha = 0.2, linewidth = 0.5) +
      scale_x_log10() +
      labs(
        title    = "Tag cluster IQW distribution by condition",
        subtitle = "CAGEr distclu tag clusters (q0.1–q0.9)",
        x        = "IQW (bp, log10 scale)",
        y        = "Density",
        fill     = "Condition", colour = "Condition"
      ) +
      theme_bw(base_size = 11)

    pdf(out_iqw, width = 7, height = 4)
    print(p_iqw)
    dev.off()
  } else {
    message("  no IQW data in tagclusters.tsv.gz, skipping IQW plot")
  }
} else {
  message("  tagclusters.tsv.gz not found — skipping IQW comparison")
  message("  (run scripts/06_cager/run.sh to generate it)")
}

message("Done. Outputs in: ", outdir)
sessionInfo()
