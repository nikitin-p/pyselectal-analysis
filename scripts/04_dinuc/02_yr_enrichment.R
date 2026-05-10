#!/usr/bin/env Rscript
# Stage 7: YR vs YC enrichment and dynamics across conditions.
#
# Usage (from project root):
#   Rscript scripts/04_dinuc/02_yr_enrichment.R \
#       --dinuc_1sg  results/dinuc_1Sg/all_dinuc_proportions.tsv \
#       --dinuc_2sg  results/dinuc_2Sg/all_dinuc_proportions.tsv \
#       --dinuc_all  results/dinuc_all/all_dinuc_proportions.tsv \
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

dinuc_1sg <- get_arg("--dinuc_1sg", "results/dinuc_1Sg/all_dinuc_proportions.tsv")
dinuc_2sg <- get_arg("--dinuc_2sg", "results/dinuc_2Sg/all_dinuc_proportions.tsv")
dinuc_all <- get_arg("--dinuc_all", "results/dinuc_all/all_dinuc_proportions.tsv")
samples_f <- get_arg("--samples",   "config/samples.tsv")
outdir    <- get_arg("--outdir",    "results/figures")
params_f  <- get_arg("--params",    "config/params.yaml")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
params <- yaml.load_file(params_f)
chisq_min_count <- params$chisq_min_count %||% 5
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── dinucleotide classification ───────────────────────────────────────────────
# Positions are CTSS ±1 bp (strand-aware): first base = −1, second base = +1
# Y = pyrimidine (C/T), R = purine (A/G)
YR <- c("CA", "CG", "TA", "TG")
YC <- c("CC", "TC")

classify_dinuc <- function(dn) {
  dplyr::case_when(
    dn %in% YR ~ "YR",
    dn %in% YC ~ "YC",
    TRUE        ~ "other"
  )
}

# ── load data ─────────────────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c")) |>
  select(sample_id = sample_id, condition, replicate)

load_dinuc <- function(path, type_label) {
  read_tsv(path, col_types = cols(sample = "c", dinucleotide = "c",
                                   sum_score = "d", proportion = "d")) |>
    rename(sample_id = sample) |>
    mutate(end_type = type_label,
           category = classify_dinuc(dinucleotide))
}

d1 <- load_dinuc(dinuc_1sg, "1Sg")
d2 <- load_dinuc(dinuc_2sg, "2Sg")
dall <- load_dinuc(dinuc_all, "all")
dat <- bind_rows(d1, d2, dall) |>
  left_join(meta, by = "sample_id")

# ── aggregate to YR / YC / other per sample ───────────────────────────────────
sample_cat <- dat |>
  group_by(sample_id, end_type, condition, replicate, category) |>
  summarise(count = sum(sum_score), .groups = "drop") |>
  group_by(sample_id, end_type) |>
  mutate(proportion = count / sum(count)) |>
  ungroup()

cat_levels <- c("YR", "YC", "other")
sample_cat$category <- factor(sample_cat$category, levels = cat_levels)

# ── Figure 1: per-sample YR/YC/other stacked bar, faceted by end_type ─────────
p_bar <- ggplot(sample_cat,
                aes(x = sample_id, y = proportion, fill = category)) +
  geom_col() +
  facet_wrap(~ end_type, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("YR" = "#E64B35", "YC" = "#4DBBD5",
                                "other" = "#B2B2B2")) +
  labs(title = "Dinucleotide category proportions per sample",
       x = NULL, y = "Proportion", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(outdir, "07_yr_per_sample.pdf"), p_bar, width = 10, height = 7)
message("Saved: 07_yr_per_sample.pdf")

# ── Figure 2: YR proportion across conditions, 1Sg vs 2Sg ─────────────────────
yr_cond <- sample_cat |>
  filter(category == "YR") |>
  group_by(end_type, condition) |>
  summarise(mean_prop = mean(proportion),
            sd_prop   = sd(proportion),
            n         = n(),
            .groups = "drop")

p_cond <- ggplot(yr_cond,
                 aes(x = condition, y = mean_prop, fill = end_type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_prop - sd_prop,
                    ymax = mean_prop + sd_prop),
                position = position_dodge(0.7), width = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("1Sg" = "#E64B35", "2Sg" = "#F39B7F",
                                "all" = "#4DBBD5")) +
  labs(title = "YR proportion by condition",
       x = NULL, y = "YR proportion (mean ± SD across replicates)",
       fill = "End type") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "07_yr_by_condition.pdf"), p_cond, width = 7, height = 5)
message("Saved: 07_yr_by_condition.pdf")

# ── Figure 3: full dinucleotide profile heatmap per end_type ──────────────────
# average proportion across replicates within condition
dn_heat <- dat |>
  group_by(end_type, condition, dinucleotide) |>
  summarise(mean_prop = mean(proportion), .groups = "drop") |>
  mutate(label = paste0(end_type, "\n", condition))

for (et in c("1Sg", "2Sg", "all")) {
  df_et <- dn_heat |> filter(end_type == et)
  dinuc_lvls <- sort(unique(df_et$dinucleotide))
  df_et$dinucleotide <- factor(df_et$dinucleotide, levels = dinuc_lvls)
  heat_label_colors <- ifelse(dinuc_lvls %in% YR, "#e63946", "black")
  p_heat <- ggplot(df_et, aes(x = condition, y = dinucleotide, fill = mean_prop)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", mean_prop)), size = 2.5) +
    scale_fill_gradient(low = "white", high = "#E64B35",
                        labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = paste0("Dinucleotide proportions — ", et),
         x = NULL, y = "Dinucleotide", fill = "Proportion") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(colour = heat_label_colors))
  ggsave(file.path(outdir, paste0("07_dinuc_heatmap_", et, ".pdf")),
         p_heat, width = 6, height = 5)
  message("Saved: 07_dinuc_heatmap_", et, ".pdf")
}

# ── Figure 4: z-score dinucleotide heatmap with hierarchical clustering ───────
# Rows = initiator dinucleotides; columns = samples.
# Values = z-score of alignment counts (sum_score) within each sample so that
# sample-level depth differences are removed before clustering.
# Clustering: Euclidean distance, Ward.D2 linkage (params: dinuc_dist_method,
# dinuc_hclust_method).
hclust_method_dinuc <- params$dinuc_hclust_method %||% "ward.D2"
dist_method_dinuc   <- params$dinuc_dist_method   %||% "euclidean"

for (et in c("1Sg", "2Sg", "all")) {
  df_et <- dat |>
    filter(end_type == et) |>
    select(sample_id, dinucleotide, sum_score)

  count_mat <- df_et |>
    pivot_wider(names_from = sample_id, values_from = sum_score, values_fill = 0) |>
    column_to_rownames("dinucleotide") |>
    as.matrix()

  # Remove constant rows (no variation across samples — cannot z-score)
  row_sd <- apply(count_mat, 1, sd)
  count_mat <- count_mat[row_sd > 0, , drop = FALSE]

  if (nrow(count_mat) == 0) {
    message("No variable dinucleotides for ", et, " — skipping z-score heatmap")
    next
  }

  # Z-score within each sample (column-wise): removes inter-sample depth differences
  zscore_mat <- scale(count_mat)

  # Column annotation: condition
  ann_col <- meta |>
    filter(sample_id %in% colnames(zscore_mat)) |>
    select(sample_id, condition) |>
    column_to_rownames("sample_id") |>
    as.data.frame()
  ann_col <- ann_col[colnames(zscore_mat), , drop = FALSE]

  ann_row <- data.frame(
    canonical = ifelse(rownames(zscore_mat) %in% YR, "canonical", "other"),
    row.names = rownames(zscore_mat)
  )
  ann_row_colors <- list(canonical = c(canonical = "#e63946", other = "grey80"))

  out_pdf <- file.path(outdir, paste0("07_dinuc_zscore_heatmap_", et, ".pdf"))
  pdf(out_pdf,
      width  = max(6, ncol(zscore_mat) * 0.45 + 2.5),
      height = max(4, nrow(zscore_mat) * 0.35 + 2))
  pheatmap(zscore_mat,
           clustering_method        = hclust_method_dinuc,
           clustering_distance_rows = dist_method_dinuc,
           clustering_distance_cols = dist_method_dinuc,
           annotation_col           = ann_col,
           annotation_row           = ann_row,
           annotation_colors        = ann_row_colors,
           color                    = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
           border_color             = NA,
           fontsize                 = 9,
           main                     = paste0("Initiator dinucleotide z-scores — ", et,
                                             "\n(distance: ", dist_method_dinuc,
                                             ", linkage: ", hclust_method_dinuc, ")"))
  dev.off()
  message("Saved: ", basename(out_pdf))
}

# ── Chi-square: YR vs non-YR across conditions, per end_type ─────────────────
chisq_results <- lapply(c("1Sg", "2Sg", "all"), function(et) {
  # pool replicates per condition
  pooled <- sample_cat |>
    filter(end_type == et) |>
    group_by(condition, category) |>
    summarise(count = sum(count), .groups = "drop")

  count_mat <- pooled |>
    pivot_wider(names_from = category, values_from = count, values_fill = 0) |>
    column_to_rownames("condition") |>
    as.matrix()

  # ensure YR and non-YR columns exist
  if (!"YR" %in% colnames(count_mat)) {
    message("No YR counts for ", et, " — skipping chi-square")
    return(NULL)
  }
  mat2 <- cbind(YR = count_mat[, "YR"],
                nonYR = rowSums(count_mat[, colnames(count_mat) != "YR", drop = FALSE]))

  global <- chisq.test(mat2)
  message(sprintf("[%s] Global chi-square YR vs non-YR: X²=%.2f  df=%d  p=%.3e",
                  et, global$statistic, global$parameter, global$p.value))

  conds <- rownames(mat2)
  pairs <- combn(conds, 2, simplify = FALSE)
  pairwise <- lapply(pairs, function(pair) {
    sub <- mat2[pair, , drop = FALSE]
    exp_c <- chisq.test(sub)$expected
    if (any(exp_c < chisq_min_count))
      warning(sprintf("Low expected counts: %s vs %s", pair[1], pair[2]))
    res <- chisq.test(sub)
    tibble(end_type = et, cond1 = pair[1], cond2 = pair[2],
           chi2 = res$statistic, df = res$parameter, p_value = res$p.value)
  }) |> bind_rows()

  pairwise |> mutate(p_adj = p.adjust(p_value, method = "BH"))
}) |> bind_rows()

chisq_out <- file.path(outdir, "07_yr_chisq_pairwise.tsv")
write_tsv(chisq_results, chisq_out)
message("Saved: 07_yr_chisq_pairwise.tsv")
print(chisq_results)

# ── session info ──────────────────────────────────────────────────────────────
message("\n--- sessionInfo ---")
sessionInfo()
