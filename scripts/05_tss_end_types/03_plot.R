#!/usr/bin/env Rscript
# Stage 8 step 3: visualize TSS end-type matrix.
#
# Outputs:
#   results/figures/08_tss_heatmap_pct1Sg.pdf  — heatmap rows=TSS sorted by mean pct_1Sg,
#                                                  cols=conditions, fill=pct_1Sg
#   results/figures/08_initiator_boxplot.pdf    — pct_1Sg / pct_2Sg by initiator dinucleotide
#   results/figures/08_initiator_stacked.pdf    — stacked end-type fractions per initiator
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/03_plot.R \
#       [--matrix  results/tss/tss_matrix.tsv] \
#       [--samples config/samples.tsv]          \
#       [--outdir  results/figures]             \
#       [--params  config/params.yaml]          \
#       [--min_cond_reads 5]                    \
#       [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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

matrix_f        <- get_arg("--matrix",  "results/tss/tss_matrix.tsv")
samples_f       <- get_arg("--samples", "config/samples.tsv")
outdir          <- get_arg("--outdir",  "results/figures")
params_f        <- get_arg("--params",  "config/params.yaml")
min_cond_reads  <- as.integer(get_arg("--min_cond_reads", "5"))
force           <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_heat   <- file.path(outdir, "08_tss_heatmap_pct1Sg.pdf")
out_box    <- file.path(outdir, "08_initiator_boxplot.pdf")
out_stack  <- file.path(outdir, "08_initiator_stacked.pdf")

if (!force && all(file.exists(c(out_heat, out_box, out_stack)))) {
  message("All outputs exist, skipping (use --force to rerun).")
  quit(save = "no", status = 0)
}

# ── load data ─────────────────────────────────────────────────────────────────
message("Stage 8 step 3: plotting")
mat  <- fread(matrix_f)
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)

# Join condition info
cond_map <- meta[, c("sample_id", "condition", "replicate"), drop = FALSE]
mat <- merge(mat, cond_map, by = "sample_id", all.x = TRUE)

# ── condition-level aggregation (sum replicates, recompute pct) ───────────────
cond_mat <- mat[, .(
  n_1Sg   = sum(n_1Sg),
  n_2Sg   = sum(n_2Sg),
  n_M     = sum(n_M),
  n_other = sum(n_other),
  n_total = sum(n_total)
), by = .(cluster_id, chr, tss_pos, strand, initiator, condition)]

cond_mat[n_total >= min_cond_reads, pct_1Sg   := 100 * n_1Sg   / n_total]
cond_mat[n_total >= min_cond_reads, pct_2Sg   := 100 * n_2Sg   / n_total]
cond_mat[n_total >= min_cond_reads, pct_M     := 100 * n_M     / n_total]
cond_mat[n_total >= min_cond_reads, pct_other := 100 * n_other / n_total]
cond_mat[n_total < min_cond_reads,
         c("pct_1Sg","pct_2Sg","pct_M","pct_other") := NA_real_]

# ── Figure 1: heatmap pct_1Sg sorted by mean across conditions ───────────────
message("  heatmap: ", out_heat)

# Compute mean pct_1Sg per cluster (ignoring NA)
cluster_order <- cond_mat[!is.na(pct_1Sg),
                           .(mean_pct1Sg = mean(pct_1Sg)), by = cluster_id]
setorder(cluster_order, mean_pct1Sg)
cluster_order[, row_rank := .I]

heat_dat <- merge(cond_mat[!is.na(pct_1Sg)], cluster_order, by = "cluster_id")

# Limit label density: only show tick for every Nth cluster
n_clusters <- uniqueN(heat_dat$cluster_id)
message("  clusters in heatmap: ", n_clusters)

p_heat <- ggplot(heat_dat, aes(x = condition, y = row_rank, fill = pct_1Sg)) +
  geom_raster() +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      name = "% 1Sg", limits = c(0, 100), na.value = "grey90") +
  scale_y_continuous(expand = c(0, 0),
                     labels = NULL, breaks = NULL) +
  labs(
    title = "TSS end-type: % 1Sg per cluster × condition",
    subtitle = paste0(n_clusters, " tag clusters (sorted by mean % 1Sg, low → high)"),
    x = "Condition", y = "TSS cluster"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

pdf(out_heat, width = 5, height = max(4, n_clusters / 300 + 3))
print(p_heat)
dev.off()

# ── Figure 2: boxplot pct_1Sg and pct_2Sg by initiator dinucleotide ──────────
message("  initiator boxplot: ", out_box)

# Focus on biologically meaningful initiators (≥5 clusters)
init_counts <- cond_mat[!is.na(initiator) & !is.na(pct_1Sg),
                         .N, by = initiator]
keep_init   <- init_counts[N >= 5, initiator]

box_dat <- cond_mat[initiator %in% keep_init & !is.na(pct_1Sg)] |>
  as.data.frame() |>
  pivot_longer(cols = c(pct_1Sg, pct_2Sg),
               names_to = "end_type", values_to = "pct") |>
  mutate(end_type = recode(end_type, pct_1Sg = "1Sg", pct_2Sg = "2Sg"))

# Order initiators by median pct_1Sg descending
init_order <- box_dat |>
  filter(end_type == "1Sg") |>
  group_by(initiator) |>
  summarise(med = median(pct, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(med))
box_dat$initiator <- factor(box_dat$initiator, levels = init_order$initiator)

p_box <- ggplot(box_dat, aes(x = initiator, y = pct, fill = end_type)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.35) +
  scale_fill_manual(values = c("1Sg" = "#4361ee", "2Sg" = "#f72585"),
                    name = "End type") +
  labs(
    title = "RT template switching by initiator dinucleotide",
    subtitle = "% reads per TSS cluster × condition; sorted by median % 1Sg",
    x = "Initiator (TSS−1 | TSS)", y = "% reads"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

pdf(out_box, width = 8, height = 5)
print(p_box)
dev.off()

# ── Figure 3: stacked bar — mean end-type fractions per initiator ─────────────
message("  stacked bar: ", out_stack)

stack_dat <- cond_mat[initiator %in% keep_init & !is.na(pct_1Sg),
                       .(mean_1Sg   = mean(pct_1Sg,   na.rm = TRUE),
                         mean_2Sg   = mean(pct_2Sg,   na.rm = TRUE),
                         mean_M     = mean(pct_M,     na.rm = TRUE),
                         mean_other = mean(pct_other, na.rm = TRUE)),
                       by = initiator] |>
  as.data.frame() |>
  pivot_longer(cols = starts_with("mean_"),
               names_to = "end_type", values_to = "pct") |>
  mutate(end_type = sub("mean_", "", end_type),
         end_type = factor(end_type, levels = c("1Sg", "2Sg", "M", "other")))

stack_dat$initiator <- factor(stack_dat$initiator, levels = init_order$initiator)

type_colors <- c("1Sg" = "#4361ee", "2Sg" = "#f72585",
                 "M"   = "#06d6a0", "other" = "#adb5bd")

p_stack <- ggplot(stack_dat, aes(x = initiator, y = pct, fill = end_type)) +
  geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = type_colors, name = "End type") +
  labs(
    title = "Mean 5′-end type composition by initiator",
    x = "Initiator (TSS−1 | TSS)", y = "Mean % reads"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

pdf(out_stack, width = 8, height = 5)
print(p_stack)
dev.off()

message("Done. Outputs in: ", outdir)
sessionInfo()
