#!/usr/bin/env Rscript
# Stage 3, step 2: aggregate pyselectal --count TSVs, draw histograms,
# heatmap of replicate/condition reproducibility, chi-square tests.
#
# Usage (from project root):
#   Rscript scripts/03_five_prime/02_plot.R \
#       --counts_dir results/five_prime \
#       --samples    config/samples.tsv \
#       --outdir     results/figures \
#       --params     config/params.yaml \
#       --types      coarse|fine
#
# --types selects the vocabulary the raw pyselectal --count types are folded into.
#   coarse (default) — the historical Stage 3 categories: the G-appending ladder, its two
#     two-base siblings and M, with every other clip swept into a residual. This is the set
#     END_TYPES_COARSE names in scripts/utils/end_types.R.
#   fine — every type of END_TYPES_REAL resolved by name, so 1Sa, 1St and 1Sc leave the
#     residual and appear as their own columns. Its outputs carry a _fine suffix, so a fine
#     run is a companion to the coarse one rather than a replacement.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(patchwork)
  library(ggdendro)
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
type_set   <- get_arg("--types",      "coarse")
samples_f  <- get_arg("--samples",    "config/samples.tsv")
outdir     <- get_arg("--outdir",     "results/figures")
params_f   <- get_arg("--params",     "config/params.yaml")

if (!type_set %in% c("coarse", "fine"))
  stop("--types must be one of: coarse, fine")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Output stems. The fine run writes beside the coarse one instead of over it.
tag <- if (type_set == "coarse") "" else paste0("_", type_set)
out_f <- function(stem, ext) file.path(outdir, paste0(stem, tag, ".", ext))

params <- yaml.load_file(params_f)
hclust_method  <- params$replicate_hclust_method %||% "complete"
chisq_min_count <- params$chisq_min_count %||% 5

# ── helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# fmt_n(): compact alignment-count labels drawn on top of every histogram bar
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
# When this script is source()d rather than run directly, commandArgs(FALSE)
# reflects the OUTER script's --file, not this file's own path, so the path
# above can resolve to the wrong directory depth. Fall back to the
# project-root-relative path in that case.
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

# Category labels follow the current pyselectal selector syntax, in which the
# soft-clip length is implied by the pattern itself (Sg = 1Sg, Sgg = 2Sgg,
# Sggg = 3Sggg, Stg = 2Stg, Sga = 2Sga). The `--count` output still prints the
# explicit length prefix, so the raw type strings are matched with the numeric
# form and relabelled here. The labels are literal clip sequences written
# 5' -> 3': "Sgg" means the two clipped bases are exactly "GG", not "first
# base is G", so the match against the raw type must be exact rather than a
# prefix match.
classify_type <- function(type) {
  # Sg                 → soft-clip length 1, base G
  # Sgg                → soft-clip length 2, both bases G
  # Sggg               → soft-clip length 3, all G
  # Stg                → soft-clip length 2, bases T then G
  # Sga                → soft-clip length 2, bases G then A
  # M                  → mapped 5' end (no soft-clip), e.g. 36Mctcac
  # other_soft_clipped → pyselectal collapsed rare soft-clipped types
  # other_mapped       → pyselectal collapsed rare mapped types
  dplyr::case_when(
    type == "1Sg"                      ~ "Sg",
    type == "2Sgg"                     ~ "Sgg",
    type == "3Sggg"                    ~ "Sggg",
    type == "2Stg"                     ~ "Stg",
    type == "2Sga"                     ~ "Sga",
    grepl("^[0-9]+M", type)            ~ "M",
    grepl("^other_soft_clipped", type) ~ "other_soft_clipped",
    grepl("^other_mapped", type)       ~ "other_mapped",
    TRUE                               ~ "other_soft_clipped"
  )
}

# The fine vocabulary needs no relabelling table at all: pyselectal's own raw type strings
# for the clipped types are exactly the END_TYPES_REAL labels, so the
# clipped arm is a membership test against the vocabulary rather than a list of cases. Only
# the mapped types need matching, since --count prints them with their length and their
# mapped prefix (e.g. "48Maagaa"). The residual keeps Stage 3's soft-clipped / mapped split,
# which the fine vocabulary's single "other" would throw away.
CLIPPED_TYPES <- setdiff(END_TYPES_REAL, "M")

classify_type_fine <- function(type) {
  dplyr::case_when(
    type %in% CLIPPED_TYPES            ~ type,
    grepl("^[0-9]+M", type)            ~ "M",
    grepl("^other_mapped", type)       ~ "other_mapped",
    TRUE                               ~ "other_soft_clipped"
  )
}

# ── load counts ───────────────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"))

tsv_files <- list.files(counts_dir, pattern = "_counts\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0) stop("No *_counts.tsv files found in: ", counts_dir)

# The counts directory is written once per library and is never pruned, so it
# outlives the sample sheet: a library commented out of the sheet keeps its
# count table on disk. Restricting to the sheet here is what makes such a
# removal take effect, since the join below would otherwise carry the orphaned
# library through under its raw identifier and plot it with the rest.
keep <- sub("_counts\\.tsv$", "", basename(tsv_files)) %in% meta$sample_id
if (!any(keep)) stop("No *_counts.tsv file in ", counts_dir,
                     " matches a sample in: ", samples_f)
if (any(!keep))
  message("Ignoring ", sum(!keep), " count table(s) absent from the sample sheet: ",
          paste(sub("_counts\\.tsv$", "", basename(tsv_files[!keep])), collapse = ", "))
tsv_files <- tsv_files[keep]

counts_raw <- lapply(tsv_files, function(f) {
  sid <- sub("_counts\\.tsv$", "", basename(f))
  df  <- read_tsv(f, col_types = cols(type = "c", count = "d"))
  df$sample_id <- sid
  df
}) |> bind_rows()

# ── aggregate into broad categories ──────────────────────────────────────────
counts_agg <- counts_raw |>
  mutate(category = if (type_set == "coarse") classify_type(type)
                    else classify_type_fine(type)) |>
  group_by(sample_id, category) |>
  summarise(count = sum(count), .groups = "drop") |>
  group_by(sample_id) |>
  mutate(freq = count / sum(count)) |>
  ungroup()

# join metadata
counts_agg <- left_join(counts_agg, meta, by = "sample_id")

# Readable sample labels, matching the condition names used in the stage 7
# figures ("glucose r1" rather than "yeast_cage_glu_r1"). Samples whose metadata
# is missing keep their raw identifier.
counts_agg$sample_label <- counts_agg$sample_id
if ("condition" %in% colnames(counts_agg)) {
  lab <- counts_agg$condition
  if ("replicate" %in% colnames(counts_agg)) {
    has_rep <- !is.na(lab) & !is.na(counts_agg$replicate)
    lab[has_rep] <- paste0(lab[has_rep], " r", counts_agg$replicate[has_rep])
  }
  counts_agg$sample_label <- ifelse(is.na(lab), counts_agg$sample_id, lab)
}

sample_labels <- setNames(counts_agg$sample_label, counts_agg$sample_id)
sample_labels <- sample_labels[!duplicated(names(sample_labels))]

# Fixed display order: the soft-clipped types first, then the residual soft clips, then the
# mapped ends. The fine vocabulary keeps that shape and simply spells out every clipped type,
# in the canonical order END_TYPES defines.
if (type_set == "coarse") {
  cat_levels <- c("Sg", "Sgg", "Sggg", "Stg", "Sga",
                  "other_soft_clipped", "M", "other_mapped")

  cat_colours <- c("Sg"                 = "#E64B35",
                   "Sgg"                = "#F39B7F",
                   "Sggg"               = "#FCBBA1",
                   "Stg"                = unname(END_TYPE_COLORS_BRIGHT["2Stg"]),
                   "Sga"                = unname(END_TYPE_COLORS_BRIGHT["2Sga"]),
                   "other_soft_clipped" = "#B2B2B2",
                   "M"                  = "#4DBBD5",
                   "other_mapped"       = "#7F7F7F")
} else {
  cat_levels <- c(CLIPPED_TYPES, "other_soft_clipped", "M", "other_mapped")

  cat_colours <- c(END_TYPE_COLORS_BRIGHT[CLIPPED_TYPES],
                   "other_soft_clipped" = "#B2B2B2",
                   "M"                  = "#4DBBD5",
                   "other_mapped"       = "#7F7F7F")
}
counts_agg$category <- factor(counts_agg$category, levels = cat_levels)

# ── abundance check: every displayed category should hold > 1% of the data ────
cat_overall <- counts_agg |>
  group_by(category) |>
  summarise(count = sum(count), .groups = "drop") |>
  mutate(pct = 100 * count / sum(count)) |>
  arrange(match(category, cat_levels))

write_tsv(cat_overall, out_f("03_category_overall", "tsv"))
message("Overall category shares (all samples pooled):")
print(as.data.frame(cat_overall))

low_cats <- cat_overall |> filter(pct <= 1)
if (nrow(low_cats) > 0) {
  warning(sprintf("Categories below 1%% of all alignments: %s",
                  paste(sprintf("%s (%.2f%%)", low_cats$category, low_cats$pct),
                        collapse = ", ")))
}

# ── per-sample histograms ─────────────────────────────────────────────────────
p_hist <- ggplot(counts_agg, aes(x = category, y = freq, fill = category)) +
  geom_col() +
  geom_text(aes(label = fmt_n(count)), vjust = -0.35, size = 2.3, colour = "grey25") +
  facet_wrap(~ sample_label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  scale_fill_manual(values = cat_colours) +
  labs(title = "5' end type frequencies per sample",
       subtitle = "Bar labels give the number of alignments in each category",
       x = NULL, y = "Frequency", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(out_f("03_per_sample_histogram", "pdf"),
       p_hist, width = 10 + 0.5 * max(0, length(cat_levels) - 8), height = 6)
message("Saved: ", basename(out_f("03_per_sample_histogram", "pdf")))

# ── heatmap of Z-scores: samples × categories ────────────────────────────────
# The end types are not clustered here either: the columns are ordered by their
# pooled abundance, from the most to the least frequent end type, and only the
# samples are clustered, which is what the replicate check needs.
count_mat_full <- counts_agg |>
  select(sample_id, category, count) |>
  pivot_wider(names_from = category, values_from = count, values_fill = 0) |>
  column_to_rownames("sample_id") |>
  as.matrix()

# render_replicate_heatmap(): builds and saves one replicate Z-score heatmap.
# `drop_cats` removes categories from the matrix before clustering and display
# (not just from the dendrogram distance), which is what the "no Sg" companion
# figure below uses: 1Sg/Sg dominates every profile by an order of magnitude,
# so leaving it in the panel while excluding it only from clustering would
# still let the printed Z-scores and column widths be dominated by it.
render_replicate_heatmap <- function(drop_cats = character(0), out_stem = "03_replicate_heatmap",
                                     title_suffix = "") {
  cat_overall_i <- cat_overall |> filter(!category %in% drop_cats)
  cat_desc_i    <- as.character(cat_overall_i$category[order(-cat_overall_i$count)])
  keep_cats_i   <- intersect(cat_desc_i, colnames(count_mat_full))
  count_mat_i   <- count_mat_full[, keep_cats_i, drop = FALSE]

  # Z-score of log-transformed counts (row-wise): each sample normalized independently
  zscore_mat <- t(scale(t(log10(count_mat_i + 1))))

  # sample order comes from hierarchical clustering of the Z-score profiles
  sample_hc    <- hclust(dist(zscore_mat), method = hclust_method)
  sample_order <- rownames(zscore_mat)[sample_hc$order]

  z_long <- as.data.frame(zscore_mat) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "category", values_to = "z") |>
    mutate(category  = factor(category, levels = keep_cats_i),
           sample_id = factor(sample_id, levels = sample_order))

  # top panel: total alignments per end type, pooled over all samples, on a log
  # scale — the same view as 03_per_sample_histogram but for the whole dataset
  top_df <- cat_overall_i |>
    filter(category %in% keep_cats_i) |>
    mutate(category = factor(as.character(category), levels = keep_cats_i))

  p_top <- ggplot(top_df, aes(x = category, y = count, fill = category)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", fmt_n(count), pct)),
              vjust = -0.2, size = 2.5, colour = "grey25", lineheight = 0.9) +
    # The axis sits on the right: anything drawn left of this panel would push the
    # aligned heatmap column away from the dendrogram.
    scale_y_log10(expand = expansion(mult = c(0, 0.35)), position = "right") +
    scale_x_discrete(limits = keep_cats_i, drop = FALSE) +
    scale_fill_manual(values = cat_colours, guide = "none") +
    labs(x = NULL, y = "Alignments (log10)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin  = margin(2, 2, 2, 2))

  # The fine run has too many columns for per-cell numbers to stay readable, so
  # only the coarse figure prints the Z-score inside each tile.
  p_heat <- ggplot(z_long, aes(x = category, y = sample_id, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    (if (type_set == "coarse")
       geom_text(aes(label = sprintf("%.2f", z)), size = 2.6, colour = "grey15")
     else NULL) +
    scale_x_discrete(limits = keep_cats_i, drop = FALSE) +
    scale_y_discrete(position = "right", labels = sample_labels) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFFF", high = "#E64B35",
                         midpoint = 0, name = "Z-score within row (sample)") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          panel.grid      = element_blank(),
          legend.position = "bottom",
          legend.key.height = unit(0.3, "cm"),
          plot.margin     = margin(2, 2, 2, 0))

  # left panel: sample dendrogram, drawn to match the heatmap row order. Its right
  # margin is zero so that the leaves sit directly against the first tile column.
  dd  <- ggdendro::dendro_data(as.dendrogram(sample_hc), type = "rectangle")
  p_dend <- ggplot(ggdendro::segment(dd)) +
    geom_segment(aes(x = y, y = x, xend = yend, yend = xend), linewidth = 0.3) +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0))) +
    scale_y_continuous(limits = c(0.5, nrow(zscore_mat) + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 2))

  subtitle <- paste0(
    "Z-scores are taken along rows: each sample is scaled across the ",
    length(keep_cats_i), " end types, z = (log10 count - row mean) / row SD.\n",
    "End types are ordered by their pooled abundance and are not clustered; samples are clustered on their Z-score profiles.")
  if (length(drop_cats) > 0) {
    subtitle <- paste0(subtitle, "\nExcluded from this figure (matrix, clustering and display): ",
                       paste(drop_cats, collapse = ", "), ".")
  }
  if (nrow(low_cats) > 0) {
    subtitle <- paste0(subtitle, "\nBelow 1% of all alignments: ",
                       paste(sprintf("%s (%.2f%%)", low_cats$category, low_cats$pct),
                             collapse = ", "))
  }

  # Figure geometry: the tile side is fixed, and the panel sizes, the relative
  # patchwork widths/heights and the final canvas are all derived from it, so the
  # heatmap cells stay square whatever the number of samples and end types.
  n_row <- nrow(zscore_mat)
  n_col <- ncol(zscore_mat)
  cell_in    <- 0.62                       # tile side, inches
  dend_in    <- 0.75                       # dendrogram panel width
  rowlab_in  <- 1.0 + 0.062 * max(nchar(sample_labels))    # row labels + right axis
  collab_in  <- 0.25 + 0.045 * max(nchar(keep_cats_i))     # rotated column labels
  top_in     <- 1.7                        # top histogram panel
  chrome_in  <- 1.5                        # title, subtitle and legend strip

  heat_w <- n_col * cell_in
  heat_h <- n_row * cell_in

  combined <- (plot_spacer() + p_top + p_dend + p_heat) +
    plot_layout(ncol    = 2,
                widths  = unit(c(dend_in, heat_w), "in"),
                heights = unit(c(top_in, heat_h), "in")) +
    plot_annotation(title = paste0("5' end type Z-scores per sample (row-wise) - replicate reproducibility",
                                   if (type_set == "coarse") "" else " - full end-type vocabulary",
                                   title_suffix),
                    subtitle = subtitle,
                    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                                  plot.subtitle = element_text(size = 8, colour = "grey30")))

  # The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
  for (ext in c("pdf", "png"))
    ggsave(out_f(out_stem, ext), combined,
           width  = dend_in + heat_w + rowlab_in + 0.2,
           height = top_in + heat_h + collab_in + chrome_in, dpi = 200, bg = "white")
  message("Saved: ", basename(out_f(out_stem, "pdf")))
}

# The fine figure draws only the named vocabulary. Its two residual rows are
# not end types of their own but whatever the named types leave over, so they
# are left out of the matrix; the bar labels still give each type's share of
# ALL alignments, because cat_overall is computed before the drop.
residual_cats <- if (type_set == "fine") {
  c("other_soft_clipped", "other_mapped")
} else character(0)

render_replicate_heatmap(drop_cats = residual_cats)

# ── companion figure: replicate heatmap with the dominant G-appending clip
# excluded from clustering ────────────────────────────────────────────────────
# 1Sg/Sg is the CAGE/oligo-capping signature and typically dwarfs every other
# category, so it alone can drive the sample dendrogram regardless of whether
# the remaining, much rarer end types carry any independent, condition-linked
# signal. This second figure drops it (Sg in the coarse vocabulary, 1Sg in the
# fine one) before scaling and clustering, so the question it answers is
# whether the *relative proportions of the minor end types* alone are enough
# to reproduce sensible replicate/condition grouping, i.e. whether they are a
# biological feature in their own right rather than noise riding on Sg's tail.
dominant_cat <- if (type_set == "coarse") "Sg" else END_TYPE_CANONICAL
if (dominant_cat %in% colnames(count_mat_full)) {
  render_replicate_heatmap(drop_cats    = c(dominant_cat, residual_cats),
                           out_stem     = "03_replicate_heatmap_no_sg",
                           title_suffix = paste0(" - ", dominant_cat, " excluded"))
} else {
  message("'", dominant_cat, "' not present among the categories — skipping the no-",
          dominant_cat, " companion heatmap.")
}

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
  chisq_out <- out_f("03_chisq_pairwise", "tsv")
  write_tsv(pairwise, chisq_out)
  message("Saved: ", basename(chisq_out))
  message(sprintf("Global chi-square: X²=%.2f  df=%d  p=%.3e",
                  global_test$statistic, global_test$parameter, global_test$p.value))
  print(pairwise)
}

# ── session info ──────────────────────────────────────────────────────────────
message("\n--- sessionInfo ---")
sessionInfo()
