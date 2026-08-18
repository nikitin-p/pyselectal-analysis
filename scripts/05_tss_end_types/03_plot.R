#!/usr/bin/env Rscript
# Stage 8 step 3: visualize TSS end-type matrix.
#
# Outputs:
#   results/figures/08_tss_heatmap_pct1Sg.pdf  — heatmap rows=TSS sorted by mean pct_1Sg,
#                                                  cols=conditions, fill=pct_1Sg
#   results/figures/08_initiator_boxplot.pdf    — pct_1Sg / pct_2Sgg / pct_3Sggg by initiator dinucleotide
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

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
# When this script is source()d (e.g. by a test) rather than run directly,
# commandArgs(FALSE) reflects the OUTER script's --file, not this file's own
# path, so the path above can resolve to the wrong directory depth. Fall back
# to the project-root-relative path in that case.
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

# fmt_n(): compact count labels drawn on top of every histogram bar
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

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
n_cols <- end_type_cols("n_", c(END_TYPES, "total"))
cond_mat <- mat[, lapply(.SD, sum), .SDcols = n_cols,
                by = .(cluster_id, chr, tss_pos, strand, initiator, condition)]

for (typ in END_TYPES) {
  n_col <- paste0("n_", typ)
  p_col <- paste0("pct_", typ)
  cond_mat[n_total >= min_cond_reads, (p_col) := 100 * get(n_col) / n_total]
  cond_mat[n_total < min_cond_reads, (p_col) := NA_real_]
}

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

# ── Figure 2: boxplot pct_1Sg and pct_2Sgg by initiator dinucleotide ──────────
message("  initiator boxplot: ", out_box)

# Focus on biologically meaningful initiators (≥5 clusters)
init_counts <- cond_mat[!is.na(initiator) & !is.na(pct_1Sg),
                         .N, by = initiator]
keep_init   <- init_counts[N >= 5, initiator]

box_dat <- cond_mat[initiator %in% keep_init & !is.na(pct_1Sg)] |>
  as.data.frame() |>
  pivot_longer(cols = all_of(end_type_cols("pct_", END_TYPES_REAL)),
               names_to = "end_type", values_to = "pct") %>%
  mutate(end_type = factor(sub("^pct_", "", end_type), levels = END_TYPES))

# Order initiators by median pct_1Sg descending
init_order <- box_dat |>
  filter(end_type == "1Sg") |>
  group_by(initiator) |>
  summarise(med = median(pct, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(med))
canonical <- c("CA", "CG", "TA", "TG")
box_dat$initiator <- factor(box_dat$initiator, levels = init_order$initiator)
init_label_colors <- ifelse(init_order$initiator %in% canonical, "#e63946", "black")

p_box <- ggplot(box_dat, aes(x = initiator, y = pct, fill = end_type)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.35) +
  scale_fill_manual(values = END_TYPE_COLORS_BRIGHT, name = "End type") +
  labs(
    title = "RT template switching by initiator dinucleotide",
    subtitle = "% alignments per TSS cluster × condition; sorted by median % 1Sg",
    x = "Initiator (TSS−1 | TSS)", y = "% alignments"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, colour = init_label_colors),
    panel.grid.minor = element_blank()
  )

pdf(out_box, width = 8, height = 5)
print(p_box)
dev.off()

# ── Figure 3: stacked bar — mean end-type fractions per initiator ─────────────
message("  stacked bar: ", out_stack)

pct_cols  <- end_type_cols("pct_")
stack_dat <- cond_mat[initiator %in% keep_init & !is.na(pct_1Sg),
                       lapply(.SD, mean, na.rm = TRUE), .SDcols = pct_cols,
                       by = initiator] |>
  as.data.frame() |>
  pivot_longer(cols = all_of(pct_cols),
               names_to = "end_type", values_to = "pct") |>
  mutate(end_type = sub("^pct_", "", end_type),
         end_type = factor(end_type, levels = END_TYPES))

stack_dat$initiator <- factor(stack_dat$initiator, levels = init_order$initiator)

# Each bar averages one percentage per cluster x condition, so the count that
# qualifies it is the number of those observations, not a read count.
stack_n <- cond_mat[initiator %in% keep_init & !is.na(pct_1Sg), .N, by = initiator] |>
  as.data.frame() |>
  mutate(initiator = factor(initiator, levels = init_order$initiator))

p_stack <- ggplot(stack_dat, aes(x = initiator, y = pct, fill = end_type)) +
  geom_bar(stat = "identity", position = "stack", colour = "white", linewidth = 0.2) +
  geom_text(data = stack_n, aes(x = initiator, y = 101, label = fmt_n(N)),
            inherit.aes = FALSE, vjust = 0, size = 2.5, colour = "grey25") +
  scale_fill_manual(values = END_TYPE_COLORS_BRIGHT, name = "End type") +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108)) +
  labs(
    title = "Mean 5′-end type composition by initiator",
    subtitle = "Labels give the number of TSS cluster × condition observations averaged",
    x = "Initiator (TSS−1 | TSS)", y = "Mean % alignments"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, colour = init_label_colors),
    panel.grid.minor = element_blank()
  )

pdf(out_stack, width = 8, height = 5)
print(p_stack)
dev.off()

message("Done. Outputs in: ", outdir)
sessionInfo()
