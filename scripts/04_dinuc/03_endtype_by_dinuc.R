#!/usr/bin/env Rscript
# Analyze end-type distribution within each dinucleotide
# H0: G-appending frequency is independent of initiator dinucleotide

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Paths
results_dir <- "results"
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load all dinucleotide data
dinuc_all <- read.delim(file.path(results_dir, "dinuc_all", "all_dinuc_proportions.tsv"))
dinuc_1Sg <- read.delim(file.path(results_dir, "dinuc_1Sg", "all_dinuc_proportions.tsv"))
dinuc_2Sgg <- read.delim(file.path(results_dir, "dinuc_2Sgg", "all_dinuc_proportions.tsv"))
dinuc_3Sggg <- read.delim(file.path(results_dir, "dinuc_3Sggg", "all_dinuc_proportions.tsv"))

# Rename sum_score columns and join
dinuc_all <- dinuc_all %>% select(sample, dinucleotide, n_all = sum_score)
dinuc_1Sg <- dinuc_1Sg %>% select(sample, dinucleotide, n_1Sg = sum_score)
dinuc_2Sgg <- dinuc_2Sgg %>% select(sample, dinucleotide, n_2Sgg = sum_score)
dinuc_3Sggg <- dinuc_3Sggg %>% select(sample, dinucleotide, n_3Sggg = sum_score)

# Combine - use left join from dinuc_all to keep all dinucleotides
combined <- dinuc_all %>%
  left_join(dinuc_1Sg, by = c("sample", "dinucleotide")) %>%
  left_join(dinuc_2Sgg, by = c("sample", "dinucleotide")) %>%
  left_join(dinuc_3Sggg, by = c("sample", "dinucleotide")) %>%
  replace_na(list(n_1Sg = 0, n_2Sgg = 0, n_3Sggg = 0)) %>%
  mutate(
    n_M = n_all - n_1Sg - n_2Sgg - n_3Sggg,
    n_M = pmax(n_M, 0)  # safety
  )

# Calculate proportions within each dinucleotide
combined <- combined %>%
  mutate(
    pct_1Sg = n_1Sg / n_all * 100,
    pct_2Sgg = n_2Sgg / n_all * 100,
    pct_3Sggg = n_3Sggg / n_all * 100,
    pct_M = n_M / n_all * 100
  )

# Add first base category for testing
combined <- combined %>%
  mutate(first_base = substr(dinucleotide, 1, 1))

# Save combined data
write.table(combined, file.path(results_dir, "dinuc_all", "endtype_by_dinuc.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Data saved to results/dinuc_all/endtype_by_dinuc.tsv\n")

# ============================================================================
# PLOT 1: Stacked bar histogram ordered by expression (total reads per dinuc)
# ============================================================================

# Aggregate across samples, order by total expression
agg_by_dinuc <- combined %>%
  group_by(dinucleotide) %>%
  summarise(
    total = sum(n_all),
    n_1Sg = sum(n_1Sg),
    n_2Sgg = sum(n_2Sgg),
    n_3Sggg = sum(n_3Sggg),
    n_M = sum(n_M),
    .groups = "drop"
  ) %>%
  mutate(
    pct_1Sg = n_1Sg / total * 100,
    pct_2Sgg = n_2Sgg / total * 100,
    pct_3Sggg = n_3Sggg / total * 100,
    pct_M = n_M / total * 100
  ) %>%
  arrange(desc(total))

# Pivot for stacked bar
agg_long <- agg_by_dinuc %>%
  select(dinucleotide, total, pct_1Sg, pct_2Sgg, pct_3Sggg, pct_M) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "end_type", values_to = "percentage") %>%
  mutate(
    end_type = factor(gsub("pct_", "", end_type), levels = c("M", "3Sggg", "2Sgg", "1Sg")),
    dinucleotide = factor(dinucleotide, levels = agg_by_dinuc$dinucleotide)
  )

# Colors for end types
end_type_colors <- c("1Sg" = "#2166ac", "2Sgg" = "#67a9cf", "3Sggg" = "#d1e5f0", "M" = "#b2182b")

p1 <- ggplot(agg_long, aes(x = dinucleotide, y = percentage, fill = end_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = end_type_colors, name = "End type") +
  labs(
    title = "End-type distribution by initiator dinucleotide",
    subtitle = "Ordered by total expression (left = highest)",
    x = "Initiator dinucleotide",
    y = "Percentage of reads"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(fig_dir, "dinuc_endtype_histogram.pdf"), p1, width = 10, height = 6)
cat("Saved: figures/dinuc_endtype_histogram.pdf\n")

# ============================================================================
# PLOT 2: Heatmap - samples x dinucleotides, showing pct_2Sgg (or pct_1Sg)
# ============================================================================

# Heatmap of pct_1Sg by sample x dinucleotide
mat_1Sg <- combined %>%
  select(sample, dinucleotide, pct_1Sg) %>%
  pivot_wider(names_from = dinucleotide, values_from = pct_1Sg) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Order columns by total expression
col_order <- agg_by_dinuc$dinucleotide
mat_1Sg <- mat_1Sg[, col_order]

# YR annotation for columns
yr_dinucs <- c("CA", "CG", "TA", "TG")
col_annotation <- data.frame(
  row.names = col_order,
  Initiator = ifelse(col_order %in% yr_dinucs, "YR (canonical)", "Other")
)

ann_colors <- list(Initiator = c("YR (canonical)" = "#e63946", "Other" = "grey70"))

pdf(file.path(fig_dir, "dinuc_pct1Sg_heatmap.pdf"), width = 10, height = 6)
pheatmap(
  mat_1Sg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Percentage of 1Sg reads by dinucleotide",
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("white", "#2166ac"))(50)
)
dev.off()
cat("Saved: figures/dinuc_pct1Sg_heatmap.pdf\n")

# Heatmap of pct_2Sgg
mat_2Sgg <- combined %>%
  select(sample, dinucleotide, pct_2Sgg) %>%
  pivot_wider(names_from = dinucleotide, values_from = pct_2Sgg) %>%
  column_to_rownames("sample") %>%
  as.matrix()
mat_2Sgg <- mat_2Sgg[, col_order]

pdf(file.path(fig_dir, "dinuc_pct2Sgg_heatmap.pdf"), width = 10, height = 6)
pheatmap(
  mat_2Sgg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Percentage of 2Sgg reads by dinucleotide",
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("white", "#67a9cf"))(50)
)
dev.off()
cat("Saved: figures/dinuc_pct2Sgg_heatmap.pdf\n")

# ============================================================================
# STATISTICAL TEST: Is G-appending proportion dependent on first base?
# ============================================================================
# H0: The proportion of 2Sgg (among soft-clipped reads) is equal across
#     dinucleotides with different first bases (A, C, T).
# Dinucleotides starting with G are excluded (no soft-clipped G possible).

cat("\n========== STATISTICAL TESTS ==========\n\n")

# For the test, focus on dinucleotides NOT starting with G
# (those are the only ones that can have 1Sg/2Sgg/3Sggg)
test_data <- combined %>%
  filter(!grepl("^G", dinucleotide)) %>%
  group_by(first_base) %>%
  summarise(
    n_soft = sum(n_1Sg + n_2Sgg + n_3Sggg),
    n_2Sgg = sum(n_2Sgg),
    n_3Sggg = sum(n_3Sggg),
    .groups = "drop"
  ) %>%
  mutate(
    pct_2Sgg = n_2Sgg / n_soft * 100,
    pct_3Sggg = n_3Sggg / n_soft * 100
  )

cat("Summary by first base (pooled across samples):\n")
print(test_data)
cat("\n")

# Chi-square test for 2Sgg proportion across first bases
# Contingency table: rows = first_base, cols = (2Sgg, other_soft)
cont_2Sgg <- test_data %>%
  mutate(n_other = n_soft - n_2Sgg) %>%
  select(first_base, n_2Sgg, n_other) %>%
  column_to_rownames("first_base") %>%
  as.matrix()

cat("H0 (2Sgg): Proportion of 2Sgg among soft-clipped reads is independent of first base\n")
chi_2Sgg <- chisq.test(cont_2Sgg)
print(chi_2Sgg)
cat("\n")

# Chi-square test for 3Sggg proportion
cont_3Sggg <- test_data %>%
  mutate(n_other = n_soft - n_3Sggg) %>%
  select(first_base, n_3Sggg, n_other) %>%
  column_to_rownames("first_base") %>%
  as.matrix()

cat("H0 (3Sggg): Proportion of 3Sggg among soft-clipped reads is independent of first base\n")
chi_3Sggg <- chisq.test(cont_3Sggg)
print(chi_3Sggg)
cat("\n")

# ============================================================================
# PAIRWISE TESTS with correction
# ============================================================================

# Pairwise proportion tests for 2Sgg
first_bases <- c("A", "C", "T")
pairwise_results <- list()

for (i in 1:(length(first_bases) - 1)) {
  for (j in (i + 1):length(first_bases)) {
    b1 <- first_bases[i]
    b2 <- first_bases[j]

    d1 <- test_data %>% filter(first_base == b1)
    d2 <- test_data %>% filter(first_base == b2)

    # 2Sgg test
    test_2Sgg <- prop.test(
      x = c(d1$n_2Sgg, d2$n_2Sgg),
      n = c(d1$n_soft, d2$n_soft)
    )

    # 3Sggg test
    test_3Sggg <- prop.test(
      x = c(d1$n_3Sggg, d2$n_3Sggg),
      n = c(d1$n_soft, d2$n_soft)
    )

    pairwise_results[[paste0(b1, "_vs_", b2)]] <- data.frame(
      comparison = paste0(b1, " vs ", b2),
      type = c("2Sgg", "3Sggg"),
      pct_1 = c(d1$pct_2Sgg, d1$pct_3Sggg),
      pct_2 = c(d2$pct_2Sgg, d2$pct_3Sggg),
      p_value = c(test_2Sgg$p.value, test_3Sggg$p.value)
    )
  }
}

pairwise_df <- bind_rows(pairwise_results) %>%
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_fdr = p.adjust(p_value, method = "fdr"),
    significant_bonf = p_bonferroni < 0.05,
    significant_fdr = p_fdr < 0.05
  )

cat("Pairwise proportion tests (by first base):\n")
print(pairwise_df)
cat("\n")

# Save results
write.table(pairwise_df, file.path(results_dir, "dinuc_all", "first_base_pairwise_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: dinuc_all/first_base_pairwise_tests.tsv\n")

# ============================================================================
# ADDITIONAL: Test by individual dinucleotide (finer granularity)
# ============================================================================

cat("\n========== PER-DINUCLEOTIDE TESTS ==========\n\n")

# Aggregate by dinucleotide (pooled samples)
dinuc_agg <- combined %>%
  filter(!grepl("^G", dinucleotide)) %>%
  group_by(dinucleotide) %>%
  summarise(
    n_soft = sum(n_1Sg + n_2Sgg + n_3Sggg),
    n_2Sgg = sum(n_2Sgg),
    n_3Sggg = sum(n_3Sggg),
    .groups = "drop"
  ) %>%
  mutate(
    pct_2Sgg = n_2Sgg / n_soft * 100,
    pct_3Sggg = n_3Sggg / n_soft * 100,
    first_base = substr(dinucleotide, 1, 1)
  ) %>%
  arrange(desc(n_soft))

cat("Per-dinucleotide summary:\n")
print(dinuc_agg)
cat("\n")

# Chi-square: is 2Sgg proportion different across all 12 dinucleotides?
cont_dinuc <- dinuc_agg %>%
  mutate(n_other = n_soft - n_2Sgg) %>%
  select(dinucleotide, n_2Sgg, n_other) %>%
  column_to_rownames("dinucleotide") %>%
  as.matrix()

cat("H0: Proportion of 2Sgg is equal across all dinucleotides (not starting with G)\n")
chi_dinuc <- chisq.test(cont_dinuc)
print(chi_dinuc)
cat("\n")

# If significant, do post-hoc pairwise with FDR
if (chi_dinuc$p.value < 0.05) {
  cat("Global test significant. Running post-hoc pairwise tests...\n")

  dinucs <- dinuc_agg$dinucleotide
  pairwise_dinuc <- list()

  for (i in 1:(length(dinucs) - 1)) {
    for (j in (i + 1):length(dinucs)) {
      d1 <- dinuc_agg %>% filter(dinucleotide == dinucs[i])
      d2 <- dinuc_agg %>% filter(dinucleotide == dinucs[j])

      test <- prop.test(
        x = c(d1$n_2Sgg, d2$n_2Sgg),
        n = c(d1$n_soft, d2$n_soft)
      )

      pairwise_dinuc[[length(pairwise_dinuc) + 1]] <- data.frame(
        dinuc_1 = dinucs[i],
        dinuc_2 = dinucs[j],
        pct_1 = d1$pct_2Sgg,
        pct_2 = d2$pct_2Sgg,
        p_value = test$p.value
      )
    }
  }

  pairwise_dinuc_df <- bind_rows(pairwise_dinuc) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "fdr"),
      significant = p_fdr < 0.05
    ) %>%
    arrange(p_fdr)

  cat("\nTop 20 pairwise comparisons by FDR-adjusted p-value:\n")
  print(head(pairwise_dinuc_df, 20))

  write.table(pairwise_dinuc_df,
              file.path(results_dir, "dinuc_all", "dinuc_pairwise_2Sgg_tests.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nSaved: dinuc_all/dinuc_pairwise_2Sgg_tests.tsv\n")
}

cat("\n========== DONE ==========\n")
sessionInfo()
