#!/usr/bin/env Rscript
# Stage 17 step 3: compare expression across dominance classes.
#
# The libraries were subsampled to equal depth in Stage 2, so pooled raw counts
# are directly comparable across clusters and no further normalisation is
# applied. Classes too small to support a test are kept in the descriptive
# table and excluded from testing, with the exclusion recorded explicitly.
#
# Usage (from project root):
#   Rscript scripts/10_dominant_endtype/03_expression.R \
#       [--calls  results/dominance/dominance_calls.tsv] \
#       [--outdir results/dominance]                     \
#       [--figdir results/figures]                       \
#       [--params config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
calls_f  <- get_arg("--calls",  "results/dominance/dominance_calls.tsv")
outdir   <- get_arg("--outdir", "results/dominance")
figdir   <- get_arg("--figdir", "results/figures")
params_f <- get_arg("--params", "config/params.yaml")
force    <- "--force" %in% args

params  <- yaml.load_file(params_f)
# Classes smaller than this are described but not tested. It is an analysis
# threshold, so it lives in config alongside the other dom_* keys.
min_set <- as.integer(params$dom_min_class_size)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "expression_by_class.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

calls <- fread(calls_f)

# Class assignment: dominated clusters take their dominant type as class;
# clusters that passed the depth filter but were never called dominated form
# the "no_dominant" reference class. Shallow clusters are excluded entirely,
# because their end-type composition is not estimable.
calls[, class := fifelse(dominated, dom_type, "no_dominant")]
calls <- calls[pass_depth == TRUE]
calls[, log10_expr := log10(n_total)]

summ <- calls[, .(n_clusters = .N,
                  median_expr = as.numeric(median(n_total)),
                  q25 = as.numeric(quantile(n_total, 0.25)),
                  q75 = as.numeric(quantile(n_total, 0.75)),
                  mean_log10 = mean(log10_expr)),
              by = .(threshold, class)][order(threshold, -n_clusters)]

# ── the same summary under weaker cascades ───────────────────────────────────
# The binomial filter's stringency is not the same for every class, because it
# compares each class against its own genome-wide background share: a 1Sg
# cluster is tested against a background of roughly three quarters of all reads
# and an M cluster against roughly a sixth, so a shallow 1Sg cluster needs an
# extreme fraction to survive while an M cluster of the same fraction and depth
# survives easily. Clusters the filter rejects are relabelled "no_dominant",
# which raises the surviving class median and lowers the reference. The
# quantity below makes that selection effect visible rather than leaving it to
# be inferred, by recomputing the same summary under two weaker cascades:
# "pre_binomial" applies only the depth and fraction filters, and
# "pre_binomial_reproducible" adds the reproducibility filter back.
class_summary <- function(dt, cls) {
  d <- data.table::copy(dt)
  d[, class := cls]
  d[, .(n_clusters = .N,
        median_expr = as.numeric(median(n_total)),
        q25 = as.numeric(quantile(n_total, 0.25)),
        q75 = as.numeric(quantile(n_total, 0.75)),
        mean_log10 = mean(log10_expr)),
    by = .(threshold, class)]
}
variants <- rbindlist(list(
  cbind(variant = "pre_binomial",
        class_summary(calls, fifelse(calls$pass_frac, calls$dom_type,
                                     "no_dominant"))),
  cbind(variant = "pre_binomial_reproducible",
        class_summary(calls, fifelse(calls$pass_frac & calls$pass_reprod,
                                     calls$dom_type, "no_dominant"))),
  cbind(variant = "full_cascade", class_summary(calls, calls$class))))
setorder(variants, threshold, variant, -n_clusters)
fwrite(variants, file.path(outdir, "expression_prebinom.tsv"), sep = "\t")

# ── tests ─────────────────────────────────────────────────────────────────────
# Rank-biserial correlation from the Wilcoxon statistic: r = 2W/(n_a n_b) - 1.
rank_biserial <- function(a, b) {
  w <- suppressWarnings(wilcox.test(a, b))$statistic
  as.numeric(2 * w / (length(a) * length(b)) - 1)
}

tests <- list()
for (thr in unique(calls$threshold)) {
  sub <- calls[threshold == thr]
  keep <- sub[, .N, by = class][N >= min_set, class]
  skipped <- setdiff(unique(sub$class), keep)
  if (length(skipped) > 0) {
    message("  threshold ", thr, ": classes excluded from testing (n < ",
            min_set, "): ", paste(skipped, collapse = ", "))
  }
  sub <- sub[class %in% keep]
  if (uniqueN(sub$class) < 2) next

  kw <- kruskal.test(log10_expr ~ factor(class), data = sub)
  message("  threshold ", thr, ": Kruskal-Wallis p = ",
          format.pval(kw$p.value))

  cmb <- combn(sort(keep), 2)
  rows <- lapply(seq_len(ncol(cmb)), function(j) {
    a <- sub[class == cmb[1, j], n_total]
    b <- sub[class == cmb[2, j], n_total]
    data.table(threshold = thr,
               comparison = paste(cmb[1, j], "vs", cmb[2, j]),
               n_a = length(a), n_b = length(b),
               p_wilcox = suppressWarnings(wilcox.test(a, b))$p.value,
               rank_biserial = rank_biserial(a, b),
               kruskal_p = kw$p.value)
  })
  t_dt <- rbindlist(rows)
  t_dt[, q_wilcox := p.adjust(p_wilcox, method = "BH")]
  tests[[as.character(thr)]] <- t_dt
}
if (length(tests) > 0) {
  fwrite(rbindlist(tests), file.path(outdir, "expression_tests.tsv"), sep = "\t")
}

# ── figures ───────────────────────────────────────────────────────────────────
ord <- summ[, .(m = median(median_expr)), by = class][order(-m), class]
calls[, class := factor(class, levels = ord)]

p1 <- ggplot(calls, aes(class, log10_expr, fill = class)) +
  geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ threshold, ncol = 1,
             labeller = labeller(threshold = function(x)
               paste0("dominance threshold ", x))) +
  labs(x = "dominance class", y = expression(log[10]~"pooled alignments"),
       title = "Expression of TSS clusters by dominant 5' end type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figdir, "18_dom_expression_box.pdf"), p1,
       width = 9, height = 8)

p2 <- ggplot(calls, aes(log10_expr, colour = class)) +
  stat_ecdf(linewidth = 0.7) +
  facet_wrap(~ threshold, ncol = 1) +
  labs(x = expression(log[10]~"pooled alignments"), y = "cumulative fraction",
       colour = "class",
       title = "Cumulative expression distribution by dominance class") +
  theme_bw()
ggsave(file.path(figdir, "18_dom_expression_ecdf.pdf"), p2,
       width = 9, height = 8)

# The idempotency sentinel is out_tsv, so it is written LAST: an interrupted run
# or a failing ggsave must not leave a sentinel behind that makes the next run
# skip this step while its figures are missing.
fwrite(summ, out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()
