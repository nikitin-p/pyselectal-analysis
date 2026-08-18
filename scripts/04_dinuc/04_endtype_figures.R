#!/usr/bin/env Rscript
# Three alternative visualisations of the dinucleotide x end-type composition.
#
# The three heatmaps produced by 03_endtype_by_dinuc.R split one composition
# across three files, so the 100% constraint is invisible. These figures each
# present the same table in a way that keeps the constraint (or deliberately
# drops it in favour of a single statistical statement).
#
#   A. dinuc_mosaic.pdf        - marimekko + rare-type zoom
#   B. dinuc_concordance.pdf   - cap-G vs YR-initiator odds ratios
#   C. dinuc_faceted_100pct.pdf - 100% stacked bars, one facet per sample
#
# Usage:
#   Rscript scripts/04_dinuc/04_endtype_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# fmt_n(): compact count labels drawn on top of every histogram bar
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
# When this script is source()d (e.g. by test/test_build_ctss_types.R) rather
# than run directly, commandArgs(FALSE) reflects the OUTER script's --file,
# not this file's own path, so the path above can resolve to the wrong
# directory depth. Fall back to the project-root-relative path in that case.
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

results_dir <- "results"
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Shared conventions ---------------------------------------------------------

YR_DINUCS <- c("CA", "CG", "TA", "TG")
YR_COLOR <- "#e63946"

END_TYPE_COLORS <- END_TYPE_COLORS_TAB10

# Stacking order, bottom to top: the cap-G signal at the base, mapped starts on
# top, so the "genuine initiation" fraction is read off the bottom of each bar.
END_TYPE_LEVELS <- c("1Sg", "2Sgg", "3Sggg", "other", "M")

rename_samples <- function(x) {
  x <- gsub("yeast_cage_", "", x)
  x <- gsub("_r", " r", x)
  x <- gsub("ana", "Anaerobic", x)
  x <- gsub("eth", "Ethanol", x)
  x <- gsub("glu", "Glucose", x)
  x <- gsub("nit", "Nitrogen", x)
  x
}

combined <- read.delim(file.path(results_dir, "dinuc_all", "endtype_by_dinuc.tsv"))
stopifnot(all(c("n_all", "n_1Sg", "n_2Sgg", "n_3Sggg", "n_other", "n_M") %in% names(combined)))

# Pooled counts per dinucleotide, ordered by abundance.
pooled <- combined %>%
  group_by(dinucleotide) %>%
  summarise(across(c(n_all, n_1Sg, n_2Sgg, n_3Sggg, n_other, n_M), sum), .groups = "drop") %>%
  arrange(desc(n_all)) %>%
  mutate(share = n_all / sum(n_all))

dinuc_order <- pooled$dinucleotide

# ---------------------------------------------------------------------------
# A. Mosaic (marimekko) + rare-type zoom
# ---------------------------------------------------------------------------
# Column width is proportional to read abundance, column height is the 100%
# composition. Both the constraint and the fact that four dinucleotides carry
# ~88% of the data are visible at once.

GAP <- 0.002  # gap between columns, in units of total width

mosaic_x <- pooled %>%
  mutate(
    w = share - GAP,
    xmax = cumsum(share),
    xmin = xmax - w,
    xmid = (xmin + xmax) / 2
  )

mosaic <- mosaic_x %>%
  select(dinucleotide, xmin, xmax, xmid, share, n_all,
         n_1Sg, n_2Sgg, n_3Sggg, n_other, n_M) %>%
  pivot_longer(starts_with("n_") & !any_of("n_all"),
               names_to = "end_type", values_to = "n") %>%
  mutate(end_type = factor(sub("^n_", "", end_type), levels = END_TYPE_LEVELS)) %>%
  arrange(dinucleotide, end_type) %>%
  group_by(dinucleotide) %>%
  mutate(
    frac = n / sum(n),
    ymax = cumsum(frac),
    ymin = ymax - frac
  ) %>%
  ungroup()

# Column widths span ~200x, so labels placed at column midpoints collide in the
# tail. Labels are instead spread evenly across the axis and tied back to their
# column with a leader line, which keeps all 16 readable and in order.
n_dinuc <- nrow(mosaic_x)
label_df <- mosaic_x %>%
  mutate(
    label = dinucleotide,
    label_x = seq(0.5 / n_dinuc, 1 - 0.5 / n_dinuc, length.out = n_dinuc),
    color = ifelse(dinucleotide %in% YR_DINUCS, YR_COLOR, "grey20"),
    # Share plus the read count behind it; only wide enough columns can carry it.
    share_label = ifelse(share >= 0.03,
                         paste0(percent(share, accuracy = 0.1), "\n", fmt_n(n_all)),
                         "")
  )

p_mosaic <- ggplot(mosaic) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = end_type),
            color = NA) +
  geom_segment(data = label_df,
               aes(x = xmid, xend = xmid, y = 0, yend = -0.03),
               color = "grey65", linewidth = 0.25) +
  geom_segment(data = label_df,
               aes(x = xmid, xend = label_x, y = -0.03, yend = -0.07),
               color = "grey65", linewidth = 0.25) +
  geom_text(data = label_df,
            aes(x = label_x, y = -0.095, label = label, color = color),
            hjust = 0.5, vjust = 1, size = 2.9, show.legend = FALSE) +
  geom_text(data = label_df,
            aes(x = xmid, y = 1.01, label = share_label),
            vjust = 0, lineheight = 0.9, size = 2.5, color = "grey30") +
  scale_color_identity() +
  scale_fill_manual(values = END_TYPE_COLORS, name = "End type",
                    breaks = END_TYPE_LEVELS) +
  scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
  coord_cartesian(ylim = c(-0.12, 1.12), clip = "off") +
  labs(
    title = "End-type composition by initiator dinucleotide",
    subtitle = paste0("Column width is proportional to alignment abundance; each column sums to 100%. ",
                      "YR initiators in red.\nColumn labels: share of all alignments and the alignment count."),
    y = "Share of alignments at that initiator", x = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(color = "grey35")
  )

# Zoom: the three rare types on a log axis, one point per sample so that
# replicate spread is visible.
rare <- combined %>%
  select(sample, dinucleotide, pct_2Sgg, pct_3Sggg, pct_other) %>%
  pivot_longer(starts_with("pct_"), names_to = "end_type", values_to = "pct") %>%
  mutate(
    end_type = factor(sub("^pct_", "", end_type), levels = c("2Sgg", "3Sggg", "other")),
    dinucleotide = factor(dinucleotide, levels = dinuc_order)
  ) %>%
  filter(pct > 0)

n_dropped <- sum(combined$n_2Sgg == 0) + sum(combined$n_3Sggg == 0) + sum(combined$n_other == 0)

rare_axis_colors <- ifelse(dinuc_order %in% YR_DINUCS, YR_COLOR, "grey20")

p_rare <- ggplot(rare, aes(x = dinucleotide, y = pct, color = end_type)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.2, alpha = 0.8) +
  scale_color_manual(values = END_TYPE_COLORS[c("2Sgg", "3Sggg", "other")],
                     name = "End type") +
  scale_y_log10(labels = function(x) paste0(x, "%"),
                breaks = c(0.01, 0.1, 1, 10, 100)) +
  annotation_logticks(sides = "l", size = 0.3) +
  labs(
    subtitle = "Rare end types on a log scale; one point per sample. Zero values are omitted (log scale).",
    x = "Initiator dinucleotide, ordered by abundance", y = "% of alignments"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(color = rare_axis_colors),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(color = "grey35")
  )

pdf(file.path(fig_dir, "dinuc_mosaic.pdf"), width = 10, height = 8)
print(p_mosaic / p_rare + plot_layout(heights = c(1.5, 1)))
dev.off()
cat("Saved: figures/dinuc_mosaic.pdf\n")

# ---------------------------------------------------------------------------
# B. Concordance of the two independent initiation markers
# ---------------------------------------------------------------------------
# Marker 1: a soft-clipped G (1Sg), interpreted as the cap-derived G.
# Marker 2: a canonical YR initiator dinucleotide.
# If both mark genuine initiation, they should agree.
#
# G-starting initiators are excluded: a local aligner extends the alignment
# through a genomic G, so a soft-clipped G can never precede one and 1Sg is
# mechanically zero there.

or_from_2x2 <- function(a, b, c, d) {
  # a = YR & 1Sg, b = YR & not-1Sg, c = non-YR & 1Sg, d = non-YR & not-1Sg
  # as.numeric guards against integer overflow: a * d exceeds .Machine$integer.max.
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  log_or <- log((a * d) / (b * c))
  se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  data.frame(
    or = exp(log_or),
    lo = exp(log_or - 1.96 * se),
    hi = exp(log_or + 1.96 * se),
    p = 2 * pnorm(-abs(log_or / se))
  )
}

conc_input <- combined %>%
  filter(substr(dinucleotide, 1, 1) != "G") %>%
  mutate(yr = dinucleotide %in% YR_DINUCS)

per_sample <- conc_input %>%
  group_by(sample, yr) %>%
  summarise(n_1Sg = sum(n_1Sg), n_all = sum(n_all), .groups = "drop") %>%
  mutate(n_not = n_all - n_1Sg) %>%
  select(sample, yr, n_1Sg, n_not) %>%
  pivot_wider(names_from = yr, values_from = c(n_1Sg, n_not),
              names_glue = "{.value}_{ifelse(yr, 'yr', 'nonyr')}")

conc <- per_sample %>%
  rowwise() %>%
  mutate(or_from_2x2(n_1Sg_yr, n_not_yr, n_1Sg_nonyr, n_not_nonyr)) %>%
  ungroup() %>%
  mutate(
    label = rename_samples(sample),
    pct_yr = n_1Sg_yr / (n_1Sg_yr + n_not_yr) * 100,
    pct_nonyr = n_1Sg_nonyr / (n_1Sg_nonyr + n_not_nonyr) * 100,
    kind = "sample"
  )

pooled_2x2 <- conc_input %>%
  group_by(yr) %>%
  summarise(n_1Sg = sum(n_1Sg), n_all = sum(n_all), .groups = "drop") %>%
  mutate(n_not = n_all - n_1Sg)

pooled_conc <- with(pooled_2x2, or_from_2x2(
  n_1Sg[yr], n_not[yr], n_1Sg[!yr], n_not[!yr]
)) %>%
  mutate(
    label = "Pooled",
    pct_yr = pooled_2x2$n_1Sg[pooled_2x2$yr] / pooled_2x2$n_all[pooled_2x2$yr] * 100,
    pct_nonyr = pooled_2x2$n_1Sg[!pooled_2x2$yr] / pooled_2x2$n_all[!pooled_2x2$yr] * 100,
    kind = "pooled"
  )

forest <- bind_rows(
  conc %>% select(label, or, lo, hi, p, pct_yr, pct_nonyr, kind),
  pooled_conc %>% select(label, or, lo, hi, p, pct_yr, pct_nonyr, kind)
) %>%
  mutate(label = factor(label, levels = rev(c(sort(conc$label), "Pooled"))))

write.table(forest, file.path(results_dir, "dinuc_all", "yr_cap_concordance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_forest <- ggplot(forest, aes(x = or, y = label, color = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(aes(size = kind)) +
  geom_text(aes(label = sprintf("%.2f  (%.1f%% vs %.1f%%)", or, pct_yr, pct_nonyr)),
            hjust = 0, nudge_x = 0.15, size = 3, color = "grey25") +
  scale_color_manual(values = c(sample = "grey35", pooled = YR_COLOR), guide = "none") +
  scale_size_manual(values = c(sample = 2, pooled = 3.2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(
    title = "A canonical initiator and a cap-derived G mark the same positions",
    subtitle = paste0(
      "Odds of a soft-clipped G (1Sg) at YR vs non-YR initiators. ",
      "G-starting initiators excluded.\nAnnotation: odds ratio (% 1Sg at YR vs at non-YR)."
    ),
    x = "Odds ratio", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(color = "grey35")
  )

pdf(file.path(fig_dir, "dinuc_concordance.pdf"), width = 9, height = 5)
print(p_forest)
dev.off()
cat("Saved: figures/dinuc_concordance.pdf\n")

# ---------------------------------------------------------------------------
# C. Faceted 100% stacked bars
# ---------------------------------------------------------------------------

faceted <- combined %>%
  select(sample, dinucleotide, starts_with("pct_")) %>%
  pivot_longer(starts_with("pct_"), names_to = "end_type", values_to = "pct") %>%
  mutate(
    end_type = factor(sub("^pct_", "", end_type), levels = rev(END_TYPE_LEVELS)),
    dinucleotide = factor(dinucleotide, levels = dinuc_order),
    sample = factor(rename_samples(sample))
  )

facet_axis_colors <- ifelse(dinuc_order %in% YR_DINUCS, YR_COLOR, "grey20")

# Read count behind each 100% bar. Sixteen labels per facet only stay legible
# rotated, so they are printed vertically above the bars.
facet_n <- combined %>%
  transmute(
    sample = factor(rename_samples(sample), levels = levels(faceted$sample)),
    dinucleotide = factor(dinucleotide, levels = dinuc_order),
    n_all
  )

p_facet <- ggplot(faceted, aes(x = dinucleotide, y = pct, fill = end_type)) +
  geom_col(width = 0.8) +
  geom_text(data = facet_n, aes(x = dinucleotide, y = 101, label = fmt_n(n_all)),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5,
            size = 1.9, color = "grey30") +
  facet_wrap(~ sample, ncol = 2) +
  scale_fill_manual(values = END_TYPE_COLORS, name = "End type",
                    breaks = END_TYPE_LEVELS) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 25),
                     limits = c(0, 118), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "End-type composition by initiator dinucleotide, per sample",
    subtitle = paste0("Every bar sums to 100%. Dinucleotides ordered by abundance; YR initiators in red.\n",
                      "Rotated labels give the number of alignments at each initiator."),
    x = "Initiator dinucleotide", y = "Share of alignments"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               color = facet_axis_colors, size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(color = "grey35")
  )

pdf(file.path(fig_dir, "dinuc_faceted_100pct.pdf"), width = 10, height = 11)
print(p_facet)
dev.off()
cat("Saved: figures/dinuc_faceted_100pct.pdf\n")

cat("\nZero-valued rare-type cells omitted from the log-scale panel:", n_dropped, "\n")
cat("Concordance table: results/dinuc_all/yr_cap_concordance.tsv\n")

sessionInfo()
