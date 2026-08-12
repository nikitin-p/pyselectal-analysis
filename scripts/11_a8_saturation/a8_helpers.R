#!/usr/bin/env Rscript
# Stage 18 helpers: the pair universe, the A-8 proportion test, BH in log space
# and the cluster bootstrap. Sourced by 01_build_windows.R, 02_a8_test.R and
# 03_plot.R, and unit-tested by test/a8/test_helpers.R.

# ── the pair universe ─────────────────────────────────────────────────────────
A8_TYPES  <- c("1Sa", "1Sc", "1Sg", "1St", "2Sgg", "3Sggg", "M")
A8_BASES  <- c("A", "C", "G", "T")
A8_DINUCS <- as.vector(t(outer(A8_BASES, A8_BASES, paste0)))

# The 5' clip base each type implies. Under STAR --alignEndsType Local a read is
# only clipped when the genomic base at -1 differs from the clip base, so a
# clipped type can never occur on a dinucleotide starting with its own base.
# M carries no clip and is therefore unconstrained.
A8_CLIP_BASE <- c("1Sa" = "A", "1Sc" = "C", "1Sg" = "G", "1St" = "T",
                  "2Sgg" = "G", "3Sggg" = "G", "M" = NA_character_)

is_structurally_empty <- function(type, dinuc) {
  cb <- A8_CLIP_BASE[type]
  !is.na(cb) & substr(dinuc, 1L, 1L) == cb
}

# ── the window coordinate system ──────────────────────────────────────────────
# There is no coordinate zero: +1 is the first transcribed base and -1 the base
# immediately before it. The initiator dinucleotide is the pair (-1, +1).
A8_COORDS <- c(-20:-1, 1L)

a8_coord_index <- function(coord) {
  i <- match(as.integer(coord), A8_COORDS)
  if (anyNA(i)) stop("coordinate outside the window: ",
                     paste(coord[is.na(i)], collapse = ", "))
  i
}

# ── effect size ───────────────────────────────────────────────────────────────
# p_obs = n_A8 / n_obs and p_exp = n_A_bg / (n_bg_pos * n_obs), so n_obs cancels
# from the ratio and the fold change depends only on the two A tallies.
a8_log2fc <- function(n_A8, n_A_bg, n_bg_pos) {
  log2((n_bg_pos * n_A8) / n_A_bg)
}

# ── the two-sided proportion test ─────────────────────────────────────────────
# Returns -log10(p) rather than p. At the depths this stage reaches, binom.test
# underflows to exactly zero, which would make -log10(q) infinite for a large
# share of the 88 pairs and destroy both the x-axis and the row ordering.
#
# pbinom with log.p = TRUE evaluates the exact binomial tail in log space and
# never underflows, so this needs no approximate fallback and therefore has no
# branch point at which -log10(p) could jump. An earlier draft switched to a
# normal approximation past a fixed threshold and jumped by 4.4 units at the
# switch, on top of losing precision to denormals below 1e-308.
#
# Two-sided is by doubling the smaller tail. At a symmetric p_exp this coincides
# with binom.test's density-based rule; away from symmetry it is the more
# conservative of the two standard conventions.
a8_pvalue <- function(n_A8, n_obs, p_exp) {
  stopifnot(n_obs > 0, n_A8 >= 0, n_A8 <= n_obs, p_exp > 0, p_exp < 1)
  k <- round(n_A8); n <- round(n_obs)
  lower <- stats::pbinom(k,      n, p_exp, lower.tail = TRUE,  log.p = TRUE)
  upper <- stats::pbinom(k - 1L, n, p_exp, lower.tail = FALSE, log.p = TRUE)
  ln_p  <- min(log(2) + min(lower, upper), 0)   # capped at p = 1
  -ln_p / log(10)
}

# ── Benjamini-Hochberg on -log10 p ────────────────────────────────────────────
# p.adjust cannot be used, because its input would be a vector of zeros. Sorting
# p ascending is sorting -log10(p) descending, and BH's step-up minimum over
# n/j * p_(j) becomes a reverse cumulative maximum of L_(j) - log10(n/j).
bh_neg_log10 <- function(neg_log10_p) {
  n <- length(neg_log10_p)
  if (n == 0L) return(numeric(0))
  o   <- order(neg_log10_p, decreasing = TRUE)
  adj <- neg_log10_p[o] - log10(n / seq_len(n))
  q   <- pmax(rev(cummax(rev(adj))), 0)
  out <- numeric(n)
  out[o] <- q
  out
}

# ── bootstrap blocks for sub-floor positions ──────────────────────────────────
# Positions that did not reach Stage 8's min_total carry no cluster_id. They
# were excluded by depth rather than by isolation, so many of them sit inside a
# real cluster's footprint and share most of their 21 bp window with it. Each is
# therefore attached to the nearest cluster within max_dist on the same
# chromosome and strand; a position with no such neighbour becomes its own
# block, which for that position is accurate rather than assumed.
#
# `pos` must carry chr, pos, strand and a cluster_id that is NA exactly where no
# Stage 8 cluster was joined.
attach_unclustered <- function(pos, max_dist) {
  dt <- data.table::copy(pos)
  data.table::setorder(dt, chr, strand, pos)
  dt[, cluster_id := {
    cid <- cluster_id
    have <- which(!is.na(cid))
    miss <- which(is.na(cid))
    if (length(have) > 0L && length(miss) > 0L) {
      hp <- pos[have]
      j  <- findInterval(pos[miss], hp)
      # distance to the clustered position below and above each gap
      d_lo <- ifelse(j >= 1L,          pos[miss] - hp[pmax(j, 1L)],          Inf)
      d_hi <- ifelse(j + 1L <= length(hp), hp[pmin(j + 1L, length(hp))] - pos[miss], Inf)
      take <- ifelse(d_lo <= d_hi, pmax(j, 1L), pmin(j + 1L, length(hp)))
      cid[miss] <- ifelse(pmin(d_lo, d_hi) <= max_dist, cid[have[take]], NA_character_)
    }
    cid
  }, by = .(chr, strand)]
  dt[is.na(cluster_id),
     cluster_id := paste("unclustered", chr, strand, pos, sep = ":")]
  dt[]
}

# ── the cluster bootstrap ─────────────────────────────────────────────────────
# Resampling is done over pre-aggregated cluster totals rather than over rows,
# because a pair can carry hundreds of thousands of positions but only tens of
# thousands of clusters, and the statistic depends on the rows only through
# these three sums.
#   W   — per-cluster total weight (positions, or reads when read-weighted)
#   A8  — per-cluster weight carrying A at -8
#   BGA — per-cluster count of A across the background positions
bootstrap_log2fc <- function(W, A8, BGA, n_bg_pos, n_boot, seed) {
  stopifnot(length(W) == length(A8), length(W) == length(BGA))
  k <- length(W)
  if (k == 0L) return(rep(NA_real_, n_boot))
  set.seed(seed)
  out <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    i <- sample.int(k, k, replace = TRUE)
    a <- sum(A8[i]); g <- sum(BGA[i])
    out[b] <- if (a > 0 && g > 0) a8_log2fc(a, g, n_bg_pos) else NA_real_
  }
  out
}
