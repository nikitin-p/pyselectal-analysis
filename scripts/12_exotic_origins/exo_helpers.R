#!/usr/bin/env Rscript
# Stage 19 helpers: origin classification of CTSS positions and the exotic
# tallies that the split bar plot is drawn from. Pure functions only, so that
# test/exotic_origins/*.R can source this file without touching the filesystem.
#
# The dominance cascade itself is NOT reimplemented here. Stage 19 calls
# scripts/10_dominant_endtype/dom_helpers.R::call_dominance() with the CTSS
# position substituted for the distclu cluster as the unit of analysis, which
# keeps the four filters, the background shares and the FDR correction bitwise
# identical to Stage 17.

# ── origin classes ────────────────────────────────────────────────────────────

# Drawing order, from the canonical class to the undecidable one.
ORIGIN_LEVELS <- c("Sg_dom", "Sgg_dom", "other_called", "not_callable")

# Human-readable labels for the figure legend.
ORIGIN_LABELS <- c(
  Sg_dom       = "1Sg-dominated",
  Sgg_dom      = "2Sgg-dominated",
  other_called = "callable, other or no dominant type",
  not_callable = "below the depth floor"
)

# Assign each row of a position-level dominance table to one of the four
# origins. `dom` must carry dom_type, pass_depth and dominated, which is what
# call_dominance() returns. A position that failed the depth filter is never
# callable, and `dominated` already implies pass_depth, so the three assignments
# below cannot overlap.
classify_origin <- function(dom, sg_type = "1Sg", sgg_type = "2Sgg") {
  stopifnot(all(c("dom_type", "pass_depth", "dominated") %in% names(dom)))
  cls <- rep("other_called", nrow(dom))
  cls[!dom$pass_depth] <- "not_callable"
  cls[dom$dominated & dom$dom_type == sg_type]  <- "Sg_dom"
  cls[dom$dominated & dom$dom_type == sgg_type] <- "Sgg_dom"
  factor(cls, levels = ORIGIN_LEVELS)
}

# ── exotic tallies ────────────────────────────────────────────────────────────

# Count exotic signal by type and by origin, in two currencies at once.
#
# `pooled`: long data.table with position_id, type, n, already summed over
#   samples. `origin_map`: data.table with position_id and origin, covering
#   every position that appears in `pooled`.
# `usual`: the types that are not exotic. The default excludes 1Sg only, so M
#   is tallied as an exotic type.
#
# An alignment count is the number of reads of that type; a position count is
# the number of distinct CTSS carrying at least one read of that type. The two
# answer different questions — how much of this exotic's signal comes from
# 1Sg-dominated sites, and how much of its occurrence does — and neither is
# recoverable from the other. Position counts deliberately do not sum to the
# "All exotics" row, because one position usually carries several exotic types.
#
# Types outside the top `top_n` by alignment count pool into `other_label`, so
# that the figure has a bounded number of bars while every read is still drawn
# somewhere.
exotic_tallies <- function(pooled, origin_map, usual = "1Sg", top_n = 12L,
                           all_label = "All exotics", other_label = "other") {
  stopifnot(all(c("position_id", "type", "n") %in% names(pooled)))
  stopifnot(all(c("position_id", "origin") %in% names(origin_map)))

  ex <- pooled[!type %in% usual]
  n_before <- nrow(ex)
  ex <- merge(ex, origin_map, by = "position_id")
  if (nrow(ex) != n_before) {
    stop("origin_map does not cover every exotic position: ", n_before,
         " rows before the join, ", nrow(ex), " after. Every position must be ",
         "classified, including the ones below the depth floor.")
  }

  totals <- ex[, .(n = sum(n)), by = type][order(-n)]
  keep   <- utils::head(totals$type, top_n)
  ex[, type_label := data.table::fifelse(type %in% keep, type, other_label)]

  per_type <- ex[, .(alignments = sum(n),
                     positions  = data.table::uniqueN(position_id)),
                 by = .(type_label, origin)]
  overall  <- ex[, .(alignments = sum(n),
                     positions  = data.table::uniqueN(position_id)),
                 by = .(origin)]
  overall[, type_label := all_label]

  res <- data.table::rbindlist(list(per_type, overall), use.names = TRUE)

  # Fill the absent (type_label, origin) combinations with zeros, so that the
  # written table is rectangular and a downstream consumer cannot mistake an
  # absent combination for an unmeasured one.
  x_levels <- c(all_label, keep, if (other_label %in% ex$type_label) other_label)
  grid <- data.table::CJ(type_label = x_levels, origin = ORIGIN_LEVELS,
                         unique = TRUE)
  res <- merge(grid, res, by = c("type_label", "origin"), all.x = TRUE)
  res[is.na(alignments), alignments := 0L]
  res[is.na(positions),  positions  := 0L]

  res <- data.table::melt(res, id.vars = c("type_label", "origin"),
                          measure.vars = c("alignments", "positions"),
                          variable.name = "measure", value.name = "count")
  res[, type_label := factor(type_label, levels = x_levels)]
  res[, origin := factor(origin, levels = ORIGIN_LEVELS)]
  data.table::setorder(res, type_label, measure, origin)
  res[]
}

# Per-position expression of the two quantities the rainclouds are drawn from,
# restricted to positions carrying at least one exotic read.
#   total_depth  — reads of every type at the position, 1Sg included;
#   exotic_depth — reads of the exotic types only.
# The libraries were subsampled to equal depth in Stage 2, so pooled raw counts
# are comparable across positions without further normalisation.
exotic_expression <- function(pooled, origin_map, usual = "1Sg") {
  stopifnot(all(c("position_id", "type", "n") %in% names(pooled)))
  expr <- pooled[, .(total_depth  = sum(n),
                     exotic_depth = sum(n[!type %in% usual])),
                 by = position_id]
  expr <- expr[exotic_depth > 0]
  expr <- merge(expr, origin_map, by = "position_id")
  expr[, origin := factor(as.character(origin), levels = ORIGIN_LEVELS)]
  expr[]
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
