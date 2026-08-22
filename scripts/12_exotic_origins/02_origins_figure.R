#!/usr/bin/env Rscript
# Stage 19 step 2: the composite exotic-origins figure.
#
# Panel A, the split bar, asks what kind of start site each 5' end type comes
# from, in four categories defined relative to the type of the bar itself:
#
#   below_floor  the position is too shallow (< exo_min_depth pooled reads) for
#                dominance to be assessed at all, so nothing is claimed about
#                it beyond the fact that it produces this type;
#   Sg_dom       the position is callable and 1Sg dominates it, so this type is
#                a minority end at a canonical start site;
#   self_dom     the position is callable and this type dominates it, so the
#                type is the site's own signature;
#   other_dom    the position is callable and neither 1Sg nor this type
#                dominates: a third type does, or nothing reaches threshold.
#
# 1Sg carries a bar of its own, drawn first, as the reference the exotics are
# read against. Its own positions fall into self_dom rather than Sg_dom, so
# that bar shows three classes and not four, which is the intended comparison:
# how much of the canonical class sits at sites it dominates, versus how much
# of each exotic class does.
#
# Each type carries two adjacent bars, one counting alignments and one counting
# the distinct CTSS positions that carry them, and each is stacked to 100 %.
# The two currencies are shown side by side because they answer different
# questions and routinely disagree: a type can draw most of its reads from a
# handful of deep 1Sg-dominated sites while occurring mostly at shallow ones.
# The bars are proportions rather than raw counts because 1Sg and M outweigh
# every other type by orders of magnitude, so a shared count axis would flatten
# all the others; the totals are printed above each bar instead.
#
# Panels B and C are rainclouds of the positions each type DOMINATES, that is
# of the self_dom class of panel A resolved to individual positions. B is that
# type's own depth at those positions and C is the pooled depth of the same
# positions over every type, which separates "this type dominates deep sites"
# from "this type is itself deeply expressed where it dominates". Every
# position contributes to exactly one row, the type that dominates it, and 1Sg
# is a row on the same footing as the exotics. Rows are ordered by descending
# total expression, so the reading order is the order of abundance.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/02_origins_figure.R \
#       [--dominance      results/exotic_origins/position_dominance.tsv.gz] \
#       [--ctss           results/dominance/fine_typed_ctss.tsv.gz]         \
#       [--cluster_counts results/dominance/cluster_type_counts.tsv]        \
#       [--outdir         results/exotic_origins]                           \
#       [--figdir         results/figures]                                  \
#       [--params         config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
dom_f    <- get_arg("--dominance", "results/exotic_origins/position_dominance.tsv.gz")
ctss_f   <- get_arg("--ctss",      "results/dominance/fine_typed_ctss.tsv.gz")
clus_f   <- get_arg("--cluster_counts", "results/dominance/cluster_type_counts.tsv")
outdir   <- get_arg("--outdir",    "results/exotic_origins")
figdir   <- get_arg("--figdir",    "results/figures")
params_f <- get_arg("--params",    "config/params.yaml")
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "exo_helpers.R"))
source(file.path(script_dir, "..", "utils", "end_types.R"))

params    <- yaml.load_file(params_f)
top_n     <- as.integer(params$exo_top_types)
point_cap <- as.integer(params$exo_point_cap)
seed      <- as.integer(params$exo_seed)
min_total <- as.integer(params$dom_min_total)
min_depth <- as.integer(if (is.null(params$exo_min_depth)) params$dom_min_depth
                        else params$exo_min_depth)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "exotic_tallies.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

# ── input ─────────────────────────────────────────────────────────────────────

ctss <- fread(ctss_f)
ctss[, position_id := paste(chr, strand, pos, sep = ":")]
pooled <- ctss[, .(n = sum(n)), by = .(position_id, type)]

dom  <- fread(dom_f)
cmap <- dom[, .(position_id, call = as.character(call))]
stopifnot(!anyNA(cmap$call))

usual <- END_TYPE_CANONICAL   # the canonical class; "exotic" means everything else
tal   <- exotic_tallies(pooled, cmap, usual = usual, top_n = top_n)

# The bars of the rainclouds are the bars of panel A: the reference type and
# every exotic type drawn individually, with the pooled tail kept pooled.
bar_types <- setdiff(as.character(levels(tal$type_label)), ALL_EXOTICS_LABEL)
expr <- dominated_expression(pooled, cmap,
                             types = setdiff(bar_types, OTHER_TYPES_LABEL))

# ── invariants ────────────────────────────────────────────────────────────────

# Every read must be drawn exactly once: the per-type bars, the pooled
# "All exotics" bar, the reference bar and the raw table must all agree.
n_exotic <- pooled[type != usual, sum(n)]
n_usual  <- pooled[type == usual, sum(n)]
per_type_bars <- setdiff(bar_types, usual)
stopifnot(tal[type_label %in% per_type_bars & measure == "alignments",
              sum(count)] == n_exotic)
stopifnot(tal[type_label == ALL_EXOTICS_LABEL & measure == "alignments",
              sum(count)] == n_exotic)
stopifnot(tal[type_label == usual & measure == "alignments",
              sum(count)] == n_usual)

# The reference bar can never carry the Sg_dom category, because a position
# that 1Sg dominates is that bar's own self_dom. A non-zero cell here would
# mean the classifier had stopped being relative to the bar's type.
stopifnot(tal[type_label == usual & category == "Sg_dom", sum(count)] == 0)

# The rainclouds must be exactly the self_dom class of the bars, resolved to
# positions: one row per dominated position, counted once.
stopifnot(nrow(expr) == uniqueN(expr$position_id))
stopifnot(nrow(expr) == cmap[!call %in% c("below_floor", "none"), .N])

# Cross-check against Stage 17. Its cluster-level table covers exactly the
# positions that reached the Stage 8 floor, so the two tallies must agree once
# this stage's positions are restricted the same way. A mismatch means the two
# stages are reading different data and every comparison between them is void.
if (file.exists(clus_f)) {
  ref <- fread(clus_f)[, .(ref_n = sum(n)), by = type]
  deep <- pooled[, tot := sum(n), by = position_id][tot >= min_total]
  own  <- deep[, .(own_n = sum(n)), by = type]
  cmp  <- merge(ref, own, by = "type", all = TRUE)
  cmp[is.na(ref_n), ref_n := 0][is.na(own_n), own_n := 0]
  if (!isTRUE(all.equal(cmp$ref_n, cmp$own_n))) {
    stop("Per-type totals disagree with Stage 17's cluster_type_counts.tsv ",
         "when restricted to positions with pooled depth >= ", min_total,
         ". Check that dom_min_total matches the value Stage 17 ran with.")
  }
  pooled[, tot := NULL]
  message("Cross-check against Stage 17 cluster_type_counts.tsv: OK (",
          nrow(cmp), " types)")
} else {
  message("Stage 17 cluster_type_counts.tsv not found, cross-check skipped")
}

# ── row order and summaries ───────────────────────────────────────────────────

# The rainclouds are ordered by descending total expression, that is by the sum
# of the type's own reads over the positions it dominates, which is the
# quantity panel B draws.
row_ord <- expr[, .(total = sum(own_depth)), by = type][order(-total)]$type
expr[, type := factor(as.character(type), levels = row_ord)]

expr_summ <- expr[, .(n_positions = .N,
                      total_own_depth = sum(own_depth),
                      median_own    = as.numeric(median(own_depth)),
                      q25_own       = as.numeric(quantile(own_depth, 0.25)),
                      q75_own       = as.numeric(quantile(own_depth, 0.75)),
                      median_total  = as.numeric(median(total_depth)),
                      q25_total     = as.numeric(quantile(total_depth, 0.25)),
                      q75_total     = as.numeric(quantile(total_depth, 0.75))),
                  by = type]
setorder(expr_summ, -total_own_depth)
fwrite(expr_summ, file.path(outdir, "dominated_expression_summary.tsv"),
       sep = "\t")
print(expr_summ)

# The per-position table is written out as well as summarised, because the
# Quarto report redraws the rainclouds at render time and a median with an
# interquartile range cannot reproduce a distribution. It is the only table of
# this stage that is per-position rather than per-class, and it is small enough
# to keep compressed alongside the summaries.
fwrite(expr[, .(position_id, type, own_depth, total_depth)],
       file.path(outdir, "dominated_expression.tsv.gz"), sep = "\t",
       compress = "gzip")

# ── panel A: the split bar ────────────────────────────────────────────────────

# Grey for the class that is a limit of the measurement rather than a property
# of the positions in it; blue for the canonical start site, green for the
# type's own, orange for everything else that was called.
CATEGORY_COLS <- c(Sg_dom      = "#1f78b4",
                   self_dom    = "#33a02c",
                   other_dom   = "#ff7f00",
                   below_floor = "#bdbdbd")

tal[, measure_lab := factor(measure, levels = c("alignments", "positions"),
                            labels = c("aln", "pos"))]
bar_tot <- tal[, .(total = sum(count)), by = .(type_label, measure_lab)]
fmt_n <- function(x) formatC(x, format = "d", big.mark = " ")

pA <- ggplot(tal, aes(measure_lab, count, fill = category)) +
  geom_col(position = "fill", width = 0.75) +
  geom_text(data = bar_tot, inherit.aes = FALSE,
            aes(measure_lab, 1.03, label = fmt_n(total)),
            size = 2.1, angle = 90, hjust = 0) +
  facet_grid(cols = vars(type_label), switch = "x") +
  # The headroom for the printed totals is set through coord_cartesian and not
  # through the scale's limits: position = "fill" is applied after the scale, so
  # a limit here would drop every raw count above it before the bar is stacked.
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0, 1.3)) +
  scale_fill_manual(values = CATEGORY_COLS, breaks = CATEGORY_LEVELS,
                    labels = CATEGORY_LABELS[CATEGORY_LEVELS]) +
  labs(x = NULL, y = "share of the bar",
       fill = "the start site this read comes from",
       title = "A  What kind of start site each 5' end type comes from",
       subtitle = paste0(
         "The categories are relative to the bar: \"this end type\" means the ",
         "type named below it. ", usual, " is drawn first as the\n",
         "reference, and it can only take three of them, because a position it ",
         "dominates is its own rather than\nanother type's. Exotic means every ",
         "fine type other than ", usual, ", so M is included. The left bar of ",
         "each pair\ncounts alignments and the right counts the distinct CTSS ",
         "positions; totals are printed above. A position\nis assessed for ",
         "dominance only from ", min_depth, " pooled reads upwards.")) +
  theme_classic(base_size = 9) +
  theme(panel.spacing.x = unit(1.5, "pt"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "top")

# ── panels B and C: the rainclouds ────────────────────────────────────────────

# Types outside the shared vocabulary (long or rare clips, and the pooled tail)
# have no colour of their own and are drawn grey.
row_cols <- END_TYPE_COLORS_BRIGHT[row_ord]
names(row_cols) <- row_ord
row_cols[is.na(row_cols)] <- "#adb5bd"

# ggdist is not installed in this project's R libraries, so each row is drawn
# from three layers by hand: the half violin as a polygon above the row's
# baseline, a narrow box just below it, and thinned jittered points below that.
raincloud <- function(dt, value_col, title, xlab) {
  d <- copy(dt)
  d[, val := log10(get(value_col))]
  d <- d[is.finite(val)]
  # The rows are drawn bottom-up, so the order is reversed here to put the most
  # expressed type at the top.
  ord <- rev(row_ord)
  d[, idx := match(as.character(type), ord)]
  hv  <- half_violin_data(d, "val", "type", ord, height = 0.42)
  pts <- thin_by_group(d, "type", cap = point_cap, seed = seed)

  # A group too small or too flat for a kernel density contributes no polygon,
  # which happens on the test fixture and could happen on a sparse real class.
  # The layer is then omitted rather than fed an empty frame, because an empty
  # frame carries no columns for the aesthetics to resolve against.
  violin_layer <- if (nrow(hv) > 0) {
    geom_polygon(data = hv, aes(x, y, group = group, fill = group),
                 alpha = 0.75, colour = NA)
  } else NULL

  ggplot() +
    violin_layer +
    geom_boxplot(data = d, aes(x = val, y = idx - 0.14, group = type),
                 orientation = "y", width = 0.16, linewidth = 0.3,
                 outlier.shape = NA, fill = "white") +
    geom_jitter(data = pts, aes(x = val, y = idx - 0.34),
                height = 0.09, size = 0.15, alpha = 0.25, colour = "grey25") +
    scale_y_continuous(breaks = seq_along(ord), labels = ord,
                       expand = expansion(add = c(0.5, 0.6))) +
    scale_fill_manual(values = row_cols, guide = "none") +
    labs(x = xlab, y = NULL, title = title) +
    theme_bw(base_size = 9)
}

pB <- raincloud(expr, "own_depth",
                "B  Positions each type dominates: that type's own reads",
                expression(log[10]~"alignments of this type at the position"))
pC <- raincloud(expr, "total_depth",
                "C  The same positions, every end type pooled",
                expression(log[10]~"pooled alignments at the position"))

p <- pA / pB / pC + plot_layout(heights = c(1.35, 1, 1))
ggsave(file.path(figdir, "20_exo_origins.pdf"), p, width = 11, height = 13)
message("Written: ", file.path(figdir, "20_exo_origins.pdf"))

# The tables are the sentinel, so the figure is written first: a failing ggsave
# must not leave a table behind that makes the next run skip this step.
fwrite(tal, out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()
