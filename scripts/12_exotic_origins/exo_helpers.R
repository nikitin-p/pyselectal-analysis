#!/usr/bin/env Rscript
# Stage 19 helpers: the category a CTSS position falls into for a given 5' end
# type, the tallies the split bar plot is drawn from, and the per-position
# expression of the positions each type dominates. Pure functions only, so that
# test/exotic_origins/*.R can source this file without touching the filesystem.
#
# The dominance cascade itself is NOT reimplemented here. Stage 19 calls
# scripts/10_dominant_endtype/dom_helpers.R::call_dominance() with the CTSS
# position substituted for the distclu cluster as the unit of analysis, which
# keeps the four filters, the background shares and the FDR correction bitwise
# identical to Stage 17. Only the depth floor differs, because a single position
# is far shallower than a distclu cluster.

# ── the four categories ───────────────────────────────────────────────────────

# The categories are defined RELATIVE TO THE TYPE whose bar is being drawn, not
# globally: the same position is "1Sg-dominated" in the bar of one exotic type
# and "this type dominates" in the bar of another. That is the point of the
# figure. For a type t at a position p:
#
#   below_floor  p is too shallow to be called at all, so its composition is
#                not interpreted. It produces t, and that is all that is
#                claimed about it;
#   Sg_dom       p is callable and 1Sg dominates it, so t is a minority end at
#                a canonical start site;
#   self_dom     p is callable and t itself dominates it, so t is the start
#                site's own signature rather than a contaminant of someone
#                else's;
#   other_dom    p is callable and neither 1Sg nor t dominates it: either a
#                third type does, or nothing reaches the threshold.
#
# The order below is the drawing order, and it is also the priority order of
# the classifier, which matters for the 1Sg bar: 1Sg's own positions are
# self_dom rather than Sg_dom, so that bar carries exactly three classes.
CATEGORY_LEVELS <- c("Sg_dom", "self_dom", "other_dom", "below_floor")

# Human-readable labels for the figure legend.
CATEGORY_LABELS <- c(
  Sg_dom      = "1Sg dominates the position",
  self_dom    = "this end type dominates the position",
  other_dom   = "another type dominates, or none does",
  below_floor = "below the depth floor, dominance not assessed"
)

# The label of the pooled bar, and of the tail of types pooled into one bar.
ALL_EXOTICS_LABEL <- "All exotics"
OTHER_TYPES_LABEL <- "other"

# Reduce a position-level dominance table to one call per position:
#   "below_floor"  the position failed the depth filter;
#   "none"         the position is deep enough but no type was called;
#   <type>         the fine type that dominates it.
# `dom` must carry position_id, dom_type, pass_depth and dominated, which is
# what call_dominance() returns once its cluster_id has been renamed. The depth
# floor is applied LAST, so a position that failed it can never be overwritten
# by a dominated row; `dominated` already implies pass_depth in the well-formed
# case, but this ordering guarantees the floor wins even if that is violated.
position_call <- function(dom) {
  stopifnot(all(c("position_id", "dom_type", "pass_depth", "dominated") %in%
                  names(dom)))
  cl <- rep("none", nrow(dom))
  cl[dom$dominated] <- as.character(dom$dom_type)[dom$dominated]
  cl[!dom$pass_depth] <- "below_floor"
  cl
}

# Assign the (type, position) pair to one of the four categories, given the
# position's call. `call` and `type` are character vectors of the same length;
# `usual` is the baseline type the Sg_dom category is named for, the same one
# that defines "exotic" for the tallies below.
classify_category <- function(call, type, usual = "1Sg") {
  stopifnot(length(call) == length(type))
  cls <- ifelse(call == "below_floor", "below_floor",
         ifelse(call == type,          "self_dom",
         ifelse(call == usual,         "Sg_dom", "other_dom")))
  factor(cls, levels = CATEGORY_LEVELS)
}

# ── the tallies ───────────────────────────────────────────────────────────────

# Count the signal of each end type by category, in two currencies at once.
#
# `pooled`: long data.table with position_id, type, n, already summed over
#   samples. `call_map`: data.table with position_id and call, covering every
#   position that appears in `pooled`.
# `usual`: the types that are not exotic. The default excludes 1Sg only, so M
#   is tallied as an exotic type.
#
# An alignment count is the number of reads of that type; a position count is
# the number of distinct CTSS carrying at least one read of that type. The two
# answer different questions — how much of this type's signal sits at
# 1Sg-dominated sites, and how much of its occurrence does — and neither is
# recoverable from the other. Position counts deliberately do not sum to the
# pooled bar, because one position usually carries several exotic types.
#
# Three kinds of bar are produced. `usual` gets a bar of its own, drawn first,
# as the reference the exotics are read against: 1Sg's own positions are
# self_dom, so that bar shows only three classes, which is the comparison the
# figure exists to make. The pooled bar sums every exotic type. The remaining
# bars are the top `top_n` exotic types by alignment count, with the tail pooled
# into `other_label` so the figure has a bounded number of bars while every read
# is still drawn somewhere. Categories are always computed from the row's own
# fine type, never from its bar label, so a pooled bar is a sum of correctly
# classified rows rather than a class of its own.
exotic_tallies <- function(pooled, call_map, usual = "1Sg", top_n = 12L,
                           all_label = ALL_EXOTICS_LABEL,
                           other_label = OTHER_TYPES_LABEL) {
  stopifnot(all(c("position_id", "type", "n") %in% names(pooled)))
  stopifnot(all(c("position_id", "call") %in% names(call_map)))

  n_before <- nrow(pooled)
  d <- merge(pooled, call_map, by = "position_id")
  if (nrow(d) != n_before) {
    stop("call_map does not cover every position: ", n_before,
         " rows before the join, ", nrow(d), " after. Every position must be ",
         "classified, including the ones below the depth floor.")
  }
  d[, category := classify_category(as.character(call), as.character(type),
                                    usual = usual[1])]

  ex <- d[!type %in% usual, .(position_id, type, n, category)]
  totals <- ex[, .(n = sum(n)), by = type][order(-n)]
  keep   <- utils::head(totals$type, top_n)
  ex[, type_label := data.table::fifelse(type %in% keep, as.character(type),
                                         other_label)]

  tally <- function(x) x[, .(alignments = sum(n),
                             positions  = data.table::uniqueN(position_id)),
                         by = .(type_label, category)]

  per_type <- tally(ex)
  overall  <- tally(ex[, .(type_label = all_label, position_id, category, n)])
  ref      <- tally(d[type %in% usual,
                      .(type_label = as.character(type), position_id,
                        category, n)])

  res <- data.table::rbindlist(list(ref, overall, per_type), use.names = TRUE)

  # Fill the absent (type_label, category) combinations with zeros, so that the
  # written table is rectangular and a downstream consumer cannot mistake an
  # absent combination for an unmeasured one.
  x_levels <- c(intersect(usual, unique(d$type)), all_label, keep,
                if (other_label %in% ex$type_label) other_label)
  grid <- data.table::CJ(type_label = x_levels, category = CATEGORY_LEVELS,
                         unique = TRUE)
  res <- merge(grid, res, by = c("type_label", "category"), all.x = TRUE)
  res[is.na(alignments), alignments := 0L]
  res[is.na(positions),  positions  := 0L]

  res <- data.table::melt(res, id.vars = c("type_label", "category"),
                          measure.vars = c("alignments", "positions"),
                          variable.name = "measure", value.name = "count")
  res[, type_label := factor(type_label, levels = x_levels)]
  res[, category := factor(category, levels = CATEGORY_LEVELS)]
  data.table::setorder(res, type_label, measure, category)
  res[]
}

# ── expression of the positions a type dominates ──────────────────────────────

# For every position that a type actually dominates, return that type's own
# depth there and the position's pooled depth over all types.
#
# This is the self_dom class of the tallies above, resolved to individual
# positions: the rainclouds ask how productive a start site is once it has been
# established that this end type is the site's own signature. Each position
# contributes to exactly one type, the one that dominates it, so the rows are a
# partition of the callable dominated positions rather than an overlapping set.
# 1Sg is included on the same footing as the exotics, which is what makes the
# comparison between the canonical class and the rest readable off one axis.
#
# `types` optionally restricts and relabels the output the way the bars are
# labelled, pooling everything outside it into `other_label`; pass NULL to keep
# every dominating type separate.
#
# The libraries were subsampled to equal depth in Stage 2, so pooled raw counts
# are comparable across positions without further normalisation.
dominated_expression <- function(pooled, call_map, types = NULL,
                                 other_label = OTHER_TYPES_LABEL) {
  stopifnot(all(c("position_id", "type", "n") %in% names(pooled)))
  stopifnot(all(c("position_id", "call") %in% names(call_map)))

  tot <- pooled[, .(total_depth = sum(n)), by = position_id]
  called <- call_map[!call %in% c("below_floor", "none")]
  d <- merge(pooled, called, by = "position_id")
  d <- d[as.character(type) == as.character(call)]
  d <- merge(d[, .(position_id, type, own_depth = n)], tot, by = "position_id")

  if (!is.null(types)) {
    d[, type := data.table::fifelse(type %in% types, as.character(type),
                                    other_label)]
  }
  d[]
}

# ── raincloud geometry ────────────────────────────────────────────────────────

# ggdist and gghalves are absent from this project's R libraries, so the
# half-violin is built by hand: a kernel density per group, scaled to a common
# height and returned as a polygon sitting above that group's baseline. The
# categorical axis is numeric (the group's index), which lets the violin, the
# box and the points be offset from one another within the same row.
#
# `dt` must carry the value column and the group column; `value_col` is
# evaluated on the log scale by the caller, not here.
half_violin_data <- function(dt, value_col, group_col, levels_order,
                             height = 0.42, n = 512L) {
  parts <- lapply(seq_along(levels_order), function(i) {
    v <- dt[[value_col]][as.character(dt[[group_col]]) == levels_order[i]]
    v <- v[is.finite(v)]
    if (length(v) < 3L || stats::var(v) == 0) return(NULL)
    d <- stats::density(v, n = n)
    data.table::data.table(
      group = levels_order[i],
      idx   = i,
      x     = c(d$x, rev(d$x)),
      y     = c(rep(i, length(d$x)), rev(i + height * d$y / max(d$y))))
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(data.table::data.table())
  out <- data.table::rbindlist(parts)
  out[, group := factor(group, levels = levels_order)]
  out[]
}

# Thin each group to at most `cap` rows, so that the jittered point layer stays
# a readable and reasonably sized PDF. Seeded, because the figure must redraw
# identically on a rerun.
thin_by_group <- function(dt, group_col, cap, seed = 42L) {
  if (cap <= 0L) return(dt[0])
  set.seed(seed)
  dt[, .SD[if (.N <= cap) seq_len(.N) else sample.int(.N, cap)],
     by = c(group_col)]
}
