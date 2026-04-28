#!/usr/bin/env Rscript
# Stage 3, step 2: aggregate pyselectal --count TSVs, draw histograms,
# heatmap of replicate/condition reproducibility, chi-square tests.
#
# Usage (from project root):
#   Rscript scripts/03_five_prime/02_plot.R \
#       --counts_dir results/five_prime \
#       --samples    config/samples.tsv \
#       --outdir     results/figures \
#       --params     config/params.yaml

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(pheatmap)
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

counts_dir <- get_arg("--counts_dir", "results/five_prime")
samples_f  <- get_arg("--samples",    "config/samples.tsv")
outdir     <- get_arg("--outdir",     "results/figures")
params_f   <- get_arg("--params",     "config/params.yaml")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)
hclust_method  <- params$replicate_hclust_method %||% "complete"
chisq_min_count <- params$chisq_min_count %||% 5

# ── helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

classify_type <- function(type) {
  # 1Sg  → soft-clip length 1, base G
  # 2Sg  → soft-clip length 2, first base G (e.g. 2Sgg, 2Sgc, 2Sgt, 2Sga)
  # M    → mapped 5' end (no soft-clip), e.g. 36Mctcac
  # other → all remaining soft-clips
  dplyr::case_when(
    type == "1Sg"                       ~ "1Sg",
    grepl("^2Sg",  type, ignore.case = TRUE) ~ "2Sg",
    grepl("^[0-9]+M", type)             ~ "M",
    TRUE                                ~ "other"
  )
}

# ── load counts ───────────────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"))

tsv_files <- list.files(counts_dir, pattern = "_counts\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0) stop("No *_counts.tsv files found in: ", counts_dir)

counts_raw <- lapply(tsv_files, function(f) {
  sid <- sub("_counts\\.tsv$", "", basename(f))
  df  <- read_tsv(f, col_types = cols(type = "c", count = "d"))
  df$sample_id <- sid
  df
}) |> bind_rows()

# ── aggregate into broad categories ──────────────────────────────────────────
counts_agg <- counts_raw |>
  mutate(category = classify_type(type)) |>
  group_by(sample_id, category) |>
  summarise(count = sum(count), .groups = "drop") |>
  group_by(sample_id) |>
  mutate(freq = count / sum(count)) |>
  ungroup()

# join metadata
counts_agg <- left_join(counts_agg, meta, by = "sample_id")

cat_levels <- c("1Sg", "2Sg", "M", "other")
counts_agg$category <- factor(counts_agg$category, levels = cat_levels)

# ── per-sample histograms ─────────────────────────────────────────────────────
p_hist <- ggplot(counts_agg, aes(x = category, y = freq, fill = category)) +
  geom_col() +
  facet_wrap(~ sample_id, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("1Sg" = "#E64B35", "2Sg" = "#F39B7F",
                                "M"   = "#4DBBD5", "other" = "#B2B2B2")) +
  labs(title = "5' end type frequencies per sample",
       x = NULL, y = "Frequency", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "03_per_sample_histogram.pdf"),
       p_hist, width = 10, height = 6)
message("Saved: 03_per_sample_histogram.pdf")

# ── heatmap of frequencies: samples × categories ─────────────────────────────
freq_mat <- counts_agg |>
  select(sample_id, category, freq) |>
  pivot_wider(names_from = category, values_from = freq, values_fill = 0) |>
  column_to_rownames("sample_id") |>
  as.matrix()

# keep only categories with any variation
freq_mat <- freq_mat[, apply(freq_mat, 2, sd) > 0, drop = FALSE]

# annotation rows: condition + replicate (if columns exist)
ann_row <- NULL
if ("condition" %in% colnames(meta)) {
  ann_row <- meta |>
    select(sample_id, any_of(c("condition", "replicate"))) |>
    column_to_rownames("sample_id") |>
    as.data.frame()
  ann_row <- ann_row[rownames(freq_mat), , drop = FALSE]
}

pdf(file.path(outdir, "03_replicate_heatmap.pdf"), width = 6, height = max(4, nrow(freq_mat) * 0.5 + 2))
pheatmap(freq_mat,
         clustering_method = hclust_method,
         annotation_row    = ann_row,
         color             = colorRampPalette(c("#FFFFFF", "#E64B35"))(100),
         border_color      = NA,
         display_numbers   = TRUE,
         number_format     = "%.2f",
         fontsize          = 9,
         main              = "5' end type frequencies — replicate reproducibility")
dev.off()
message("Saved: 03_replicate_heatmap.pdf")

# ── chi-square tests between conditions ───────────────────────────────────────
if (!"condition" %in% colnames(counts_agg)) {
  message("No 'condition' column — skipping chi-square tests.")
} else {
  # Pool reads across replicates within each condition
  pooled <- counts_agg |>
    group_by(condition, category) |>
    summarise(count = sum(count), .groups = "drop")

  cond_levels <- unique(pooled$condition)

  # Build count matrix: conditions × categories
  count_mat <- pooled |>
    pivot_wider(names_from = category, values_from = count, values_fill = 0) |>
    column_to_rownames("condition") |>
    as.matrix()

  # Global chi-square across all conditions
  global_test <- chisq.test(count_mat)

  # Pairwise chi-square for every pair of conditions
  pairs <- combn(cond_levels, 2, simplify = FALSE)
  pairwise <- lapply(pairs, function(pair) {
    sub_mat <- count_mat[pair, , drop = FALSE]
    # Warn if expected counts are low
    exp_counts <- chisq.test(sub_mat)$expected
    if (any(exp_counts < chisq_min_count)) {
      warning(sprintf("Low expected counts in pair %s vs %s — chi-square may be unreliable",
                      pair[1], pair[2]))
    }
    res <- chisq.test(sub_mat)
    tibble(cond1 = pair[1], cond2 = pair[2],
           chi2 = res$statistic, df = res$parameter, p_value = res$p.value)
  }) |> bind_rows() |>
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  # Write results
  chisq_out <- file.path(outdir, "03_chisq_pairwise.tsv")
  write_tsv(pairwise, chisq_out)
  message("Saved: 03_chisq_pairwise.tsv")
  message(sprintf("Global chi-square: X²=%.2f  df=%d  p=%.3e",
                  global_test$statistic, global_test$parameter, global_test$p.value))
  print(pairwise)
}

# ── session info ──────────────────────────────────────────────────────────────
message("\n--- sessionInfo ---")
sessionInfo()
