#!/usr/bin/env Rscript
# Stage 18 step 3: the bubble plots, the per-position profile and the
# sensitivity comparison.
#
# Bubble size encodes |log2FC| by AREA rather than by radius, because radius
# encoding exaggerates large values quadratically. The x-axis is on a sqrt scale
# and is not capped: 02_a8_test.R computes -log10(q) in log space so every value
# is finite, and the sqrt compresses the extreme tail enough for the moderate
# pairs to remain separable.
#
# Usage:
#   Rscript scripts/11_a8_saturation/03_plot.R \
#       [--tests   results/a8/a8_pair_tests.tsv]      \
#       [--profile results/a8/a8_position_profile.tsv]\
#       [--figdir  results/figures]                   \
#       [--params  config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag); args[i + 1]
}
tests_f   <- get_arg("--tests",   "results/a8/a8_pair_tests.tsv")
profile_f <- get_arg("--profile", "results/a8/a8_position_profile.tsv")
figdir    <- get_arg("--figdir",  "results/figures")
params_f  <- get_arg("--params",  "config/params.yaml")
force     <- "--force" %in% args

# Every number the figures state about the analysis is read from the same
# configuration the analysis itself reads, so a change to any of them cannot
# leave a PDF describing a run that no longer exists.
params    <- yaml.load_file(params_f)
fdr       <- as.numeric(params$a8_fdr)
bg_plat   <- as.integer(params$a8_bg_plateau)
min_reads <- as.integer(params$a8_min_type_reads)

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
f_q <- file.path(figdir, "19_a8_bubble_by_q.pdf")
if (!force && file.exists(f_q)) {
  message("Output exists, skipping (use --force to rerun): ", f_q)
  quit(save = "no", status = 0)
}

YR_RED  <- "#e63946"   # project convention for canonical initiators
CANON   <- c("CA", "CG", "TA", "TG")

res  <- fread(tests_f)
prim <- res[config == "primary"]
if (nrow(prim) == 0L) stop("No tested pairs in the primary configuration")

prim[, label := paste0(end_type, "  |  ", dinuc)]
prim[, abs_fc := abs(log2fc)]

bubble <- function(dt, order_col, decreasing, fname, subtitle) {
  d <- copy(dt)
  setorderv(d, order_col, if (decreasing) -1L else 1L)
  d[, label := factor(label, levels = rev(d$label))]
  # Colour the y-axis label red when the initiator is canonical.
  lv   <- levels(d$label)
  dn   <- sub("^.*\\|\\s*", "", lv)
  cols <- ifelse(dn %in% CANON, YR_RED, "grey20")

  p <- ggplot(d, aes(x = neg_log10_q, y = label)) +
    geom_vline(xintercept = -log10(fdr), linetype = "dotted", colour = "grey50") +
    geom_point(aes(size = abs_fc, fill = log2fc,
                   stroke = fifelse(boot_robust, 0.7, 0.2)),
               shape = 21, colour = "grey20") +
    scale_size_area(max_size = 9, name = expression(group("|", log[2]*FC, "|"))) +
    scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = YR_RED,
                         midpoint = 0, name = expression(log[2]*FC)) +
    scale_x_sqrt(name = expression(-log[10]~(adjusted~italic(p)))) +
    labs(y = NULL,
         title = "A-8 enrichment by <5' end type, initiator dinucleotide>",
         subtitle = subtitle) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(colour = cols),
          panel.grid.major.y = element_line(linewidth = 0.2))

  ggsave(file.path(figdir, fname), p, width = 8,
         height = max(4, 0.16 * nrow(d) + 2), limitsize = FALSE)
}

bubble(prim, "neg_log10_q", TRUE, "19_a8_bubble_by_q.pdf",
       paste0("position-weighted, background ", min(bg_plat), "..", max(bg_plat),
              ", membership >= ", min_reads, " alignments; ",
              "red = enriched, blue = depleted; ",
              "thick outline = bootstrap CI excludes 0"))
bubble(prim, "log2fc", TRUE, "19_a8_bubble_by_fc.pdf",
       paste("same data ordered by log2FC, which separates the enriched and",
             "depleted ends of the spectrum"))

# ── the per-position profile ──────────────────────────────────────────────────
prof <- fread(profile_f)
pp <- ggplot(prof[coord < 0], aes(coord, pct_A)) +
  geom_line(colour = "grey40") +
  geom_point(aes(colour = coord == -8L), size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = YR_RED), guide = "none") +
  scale_x_continuous(breaks = seq(-20, -2, 2)) +
  facet_wrap(~ set, scales = "free_y") +
  labs(x = "position relative to the CTSS (+1 = first transcribed base)",
       y = "A frequency (%), position-weighted",
       title = "Per-position A frequency around the CTSS",
       subtitle = "red = position -8; the flat -20..-11 stretch is the frozen background") +
  theme_bw(base_size = 9)
ggsave(file.path(figdir, "19_a8_profile.pdf"), pp, width = 10, height = 6)

# ── the sensitivity comparison ────────────────────────────────────────────────
sens <- merge(prim[, .(pair, primary = log2fc)],
              res[config != "primary", .(pair, config, log2fc)],
              by = "pair")
ps <- ggplot(sens, aes(primary, log2fc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
  geom_point(alpha = 0.7, size = 1.6, colour = "grey20") +
  facet_wrap(~ config) +
  labs(x = expression(log[2]*FC~", primary configuration"),
       y = expression(log[2]*FC~", sensitivity configuration"),
       title = "Sensitivity of the effect size to the three fixed choices",
       subtitle = "bg_proximal = background -2..-7 and -9; w_reads = alignment-weighted; member_1 = 1-alignment membership") +
  theme_bw(base_size = 9)
ggsave(file.path(figdir, "19_a8_sensitivity.pdf"), ps, width = 9, height = 3.6)

sessionInfo()
