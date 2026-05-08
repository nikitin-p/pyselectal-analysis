#!/usr/bin/env Rscript
# Stage 9 step 2: visualize hierarchical TSS clusters.
#
# Outputs:
#   results/figures/09_dendrogram.pdf          — hclust dendrogram with cut line
#   results/figures/09_cluster_heatmap.pdf     — end-type fractions heatmap
#                                                 (rows ordered by hclust, annotated by group)
#   results/figures/09_cluster_endtypes.pdf    — stacked bar: mean end-type fractions per cluster
#   results/figures/09_cluster_annotation.pdf  — annotation category distribution per cluster
#                                                 (only if annotation column is populated)
#   results/figures/09_cluster_initiator.pdf   — initiator composition per hclust group
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/05_plot_clusters.R \
#       [--clustered  results/tss/tss_clustered.tsv] \
#       [--hclust_rds results/tss/tss_hclust.rds]    \
#       [--outdir     results/figures]               \
#       [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggdendro)
  library(readr)
})

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

clustered_f <- get_arg("--clustered",  "results/tss/tss_clustered.tsv")
hclust_rds  <- get_arg("--hclust_rds", "results/tss/tss_hclust.rds")
outdir      <- get_arg("--outdir",     "results/figures")
force       <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_dendro  <- file.path(outdir, "09_dendrogram.pdf")
out_heat    <- file.path(outdir, "09_cluster_heatmap.pdf")
out_bar     <- file.path(outdir, "09_cluster_endtypes.pdf")
out_annot   <- file.path(outdir, "09_cluster_annotation.pdf")
out_init    <- file.path(outdir, "09_cluster_initiator.pdf")

if (!force && all(file.exists(c(out_dendro, out_heat, out_bar)))) {
  message("Outputs exist, skipping (use --force to rerun).")
  quit(save = "no", status = 0)
}

# ── load data ─────────────────────────────────────────────────────────────────
message("Stage 9 step 2: plot clusters")
dt  <- fread(clustered_f)
hc  <- readRDS(hclust_rds)

k_clusters  <- uniqueN(dt$hclust_group)
type_colors <- c("1Sg" = "#4361ee", "2Sg" = "#f72585",
                 "M"   = "#06d6a0", "other" = "#adb5bd")

# Cluster group colors
group_levels <- paste0("C", sort(as.integer(sub("C", "", unique(dt$hclust_group)))))
group_pal    <- setNames(
  scales::hue_pal()(k_clusters),
  group_levels
)

# ── Figure 1: dendrogram ──────────────────────────────────────────────────────
message("  dendrogram: ", out_dendro)

# Height of cut (from cutree k): find the midpoint between the k-th and (k+1)-th merge height
merge_heights <- sort(hc$height, decreasing = TRUE)
cut_h <- if (k_clusters <= length(merge_heights)) {
  (merge_heights[k_clusters - 1] + merge_heights[k_clusters]) / 2
} else {
  merge_heights[length(merge_heights)]
}

ddata <- dendro_data(hc, type = "rectangle")
p_dendro <- ggplot() +
  geom_segment(data = segment(ddata),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3, colour = "grey40") +
  geom_hline(yintercept = cut_h, linetype = "dashed", colour = "#e63946", linewidth = 0.6) +
  annotate("text", x = 1, y = cut_h * 1.02,
           label = paste0("k = ", k_clusters), hjust = 0,
           colour = "#e63946", size = 3) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(title = "TSS hierarchical clustering (Euclidean + Ward.D2)",
       y = "Height", x = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank())

pdf(out_dendro, width = 8, height = 4)
print(p_dendro)
dev.off()

# ── Figure 2: heatmap — rows ordered by hclust leaf order ────────────────────
message("  heatmap: ", out_heat)

leaf_order  <- hc$order        # integer indices into hc$labels
leaf_labels <- hc$labels[leaf_order]   # cluster_id in hclust order

dt_ord <- dt[match(leaf_labels, dt$cluster_id)]
dt_ord[, row_rank := .I]
dt_ord$hclust_group <- factor(dt_ord$hclust_group, levels = group_levels)

heat_long <- dt_ord |>
  as.data.frame() |>
  select(row_rank, hclust_group,
         `1Sg` = mean_pct_1Sg, `2Sg` = mean_pct_2Sg,
         M = mean_pct_M, other = mean_pct_other) |>
  pivot_longer(cols = c(`1Sg`, `2Sg`, M, other),
               names_to = "end_type", values_to = "pct") |>
  mutate(end_type = factor(end_type, levels = c("1Sg", "2Sg", "M", "other")))

n_tss <- nrow(dt_ord)

p_heat <- ggplot(heat_long, aes(x = end_type, y = row_rank, fill = pct)) +
  geom_raster() +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      name = "Mean %", limits = c(0, 100), na.value = "grey90") +
  scale_y_continuous(expand = c(0, 0), labels = NULL, breaks = NULL) +
  labs(title = "TSS end-type composition (hclust order)",
       subtitle = paste0(n_tss, " clusters, k = ", k_clusters),
       x = "End type", y = "TSS (hclust order)") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

# Add group sidebar as a separate panel via patchwork if available, else standalone
p_sidebar <- ggplot(dt_ord, aes(x = 1, y = row_rank, fill = hclust_group)) +
  geom_raster() +
  scale_fill_manual(values = group_pal, name = "Group") +
  scale_y_continuous(expand = c(0, 0), labels = NULL, breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), labels = NULL, breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p_sidebar + p_heat + plot_layout(widths = c(0.07, 1))
  pdf(out_heat, width = 7, height = max(4, n_tss / 200 + 3))
  print(combined)
  dev.off()
} else {
  pdf(out_heat, width = 6, height = max(4, n_tss / 200 + 3))
  print(p_heat)
  dev.off()
}

# ── Figure 3: stacked bar — mean end-type fractions per hclust group ──────────
message("  cluster end-type bar: ", out_bar)

bar_dat <- dt |>
  as.data.frame() |>
  group_by(hclust_group) |>
  summarise(
    `1Sg`  = mean(mean_pct_1Sg,   na.rm = TRUE),
    `2Sg`  = mean(mean_pct_2Sg,   na.rm = TRUE),
    M      = mean(mean_pct_M,     na.rm = TRUE),
    other  = mean(mean_pct_other, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(`1Sg`, `2Sg`, M, other),
               names_to = "end_type", values_to = "pct") |>
  mutate(
    end_type    = factor(end_type, levels = c("1Sg", "2Sg", "M", "other")),
    hclust_group = factor(hclust_group, levels = group_levels),
    label       = paste0(hclust_group, "\n(n=", n, ")")
  )

p_bar <- ggplot(bar_dat, aes(x = hclust_group, y = pct, fill = end_type)) +
  geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.25) +
  geom_text(data = bar_dat |> distinct(hclust_group, n),
            aes(x = hclust_group, y = 102, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  scale_fill_manual(values = type_colors, name = "End type") +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(title = "End-type composition per TSS cluster",
       x = "Cluster", y = "Mean % reads") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

pdf(out_bar, width = max(5, k_clusters * 1.1), height = 5)
print(p_bar)
dev.off()

# ── Figure 4: annotation distribution per cluster ─────────────────────────────
if (!all(is.na(dt$annotation))) {
  message("  annotation bar: ", out_annot)

  annot_dat <- dt |>
    as.data.frame() |>
    filter(!is.na(annotation)) |>
    count(hclust_group, annotation) |>
    group_by(hclust_group) |>
    mutate(pct = 100 * n / sum(n),
           hclust_group = factor(hclust_group, levels = group_levels)) |>
    ungroup()

  annot_colors <- c(
    sense_genic     = "#2b9348",
    antisense_genic = "#e07a5f",
    intergenic      = "#adb5bd"
  )

  p_annot <- ggplot(annot_dat, aes(x = hclust_group, y = pct, fill = annotation)) +
    geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = annot_colors, name = "Annotation") +
    labs(title = "Gene annotation per TSS cluster",
         x = "Cluster", y = "% TSS") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  pdf(out_annot, width = max(5, k_clusters * 1.1), height = 5)
  print(p_annot)
  dev.off()
} else {
  message("  annotation column is empty — skipping annotation plot")
}

# ── Figure 5: initiator composition per hclust group ─────────────────────────
if (!all(is.na(dt$initiator))) {
  message("  initiator bar: ", out_init)

  init_dat <- dt |>
    as.data.frame() |>
    filter(!is.na(initiator)) |>
    count(hclust_group, initiator) |>
    group_by(hclust_group) |>
    mutate(pct = 100 * n / sum(n),
           hclust_group = factor(hclust_group, levels = group_levels)) |>
    ungroup()

  # Keep initiators present in ≥1% of any group
  keep_init <- init_dat |>
    group_by(initiator) |>
    filter(any(pct >= 1)) |>
    pull(initiator) |>
    unique()
  init_dat <- filter(init_dat, initiator %in% keep_init)

  p_init <- ggplot(init_dat, aes(x = hclust_group, y = pct, fill = initiator)) +
    geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.25) +
    labs(title = "Initiator dinucleotide composition per TSS cluster",
         x = "Cluster", y = "% TSS") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  pdf(out_init, width = max(5, k_clusters * 1.1), height = 5)
  print(p_init)
  dev.off()
}

message("Done. Outputs in: ", outdir)
sessionInfo()
