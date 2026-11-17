#!/usr/bin/env Rscript
# NAD branch step 8: NAD captureSeq signal and YAAG at CAGE single-A-clip positions.
#
# Reverses step 6: starts from nuclear CAGE Y|A positions (-1 = C or T, +1 = A,
# pooled depth >= 25) rich in the single-A clip and asks whether they carry
# NAD captureSeq signal (S over N) and the YAAG motif more often than matched
# positions. Matching is by stratum = -1 base x A-run length from +1 (capped at
# 5) x floor(log2 CAGE depth), as in step 6, and the test is a stratified label
# permutation.
#
# The NAD signal comes from the step-4 raw-read matrix scaled by the step-5 size
# factors, which are recomputed here exactly (DESeq2 median of ratios over the
# rows that step 5 used: the positions of <strain>_S_vs_N.tsv.gz with biotype
# "none") and checked against nad/results/07_chrM_sites/sanity.tsv.
# Spec: docs/specs/2026-11-04-nad-step8-1Sa-positions.md
#
# Tests: 1 NAD signal at single-A-clip positions (detected share, mean log2 S);
#        2 dose response by A-run and clip-share bin;
#        3 G at +3 among A-run-2 positions (YAAG), clip-rich vs 1Sg-majority;
#        4 chrM positive control.
#
# Usage (from project root):
#   Rscript nad/scripts/08_1Sa_positions/08_1Sa_positions.R [--force]
#     [--ctss results/tss/typed_ctss.tsv.gz] [--matrix nad/results/04_ctss/matrix_reads.tsv.gz]
#     [--manifest nad/scripts/00_download/priority1_runs.tsv] [--step5-dir nad/results/05_deseq]
#     [--sanity nad/results/07_chrM_sites/sanity.tsv] [--fasta data/sacCer3.fa]
#     [--outdir nad/results/08_1Sa_positions] [--figdir nad/figures]
#     [--draws 1000] [--seed 1] [--min-depth 25]
#     [--null | --plant P]   synthetic checks: shuffle the NAD values within strata,
#                            or make a fraction P of clip-rich positions NAD-detected

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
ctss_f    <- get_arg("--ctss",      "results/tss/typed_ctss.tsv.gz")
matrix_f  <- get_arg("--matrix",    "nad/results/04_ctss/matrix_reads.tsv.gz")
manifest_f <- get_arg("--manifest", "nad/scripts/00_download/priority1_runs.tsv")
s5dir     <- get_arg("--step5-dir", "nad/results/05_deseq")
sanity_f  <- get_arg("--sanity",    "nad/results/07_chrM_sites/sanity.tsv")
fasta_f   <- get_arg("--fasta",     "data/sacCer3.fa")
outdir    <- get_arg("--outdir",    "nad/results/08_1Sa_positions")
figdir    <- get_arg("--figdir",    "nad/figures")
N_DRAWS   <- as.integer(get_arg("--draws", "1000"))
SEED      <- as.integer(get_arg("--seed", "1"))
MIN_DEPTH <- as.integer(get_arg("--min-depth", "25"))
PLANT     <- as.numeric(get_arg("--plant", "0"))
NULL_MODE <- "--null" %in% args
force     <- "--force" %in% args
THRESHOLDS <- c(0.05, 0.10, 0.20)
MAIN_T     <- 0.10
ARUN_CAP   <- 5L
SHARE_BREAKS <- c(-Inf, 0, 0.02, 0.05, 0.10, 0.20, Inf)
SHARE_LABELS <- c("0", "0-2", "2-5", "5-10", "10-20", ">=20")

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "..", "..", "scripts", "utils", "end_types.R"))

out_main <- file.path(outdir, "nad_signal.tsv")
if (file.exists(out_main) && !force) { message("[skip] ", out_main, " exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

TYPE_A1 <- END_TYPES_REAL[END_TYPE_CLIP_BASE[END_TYPES_REAL] %in% "A" & startsWith(END_TYPES_REAL, "1S")]
stopifnot(length(TYPE_A1) == 1)
COL_A1 <- paste0("n_", TYPE_A1)
COL_SG <- paste0("n_", END_TYPE_CANONICAL)

# ── CAGE candidates ──────────────────────────────────────────────────────────
ctss <- fread(cmd = sprintf("gzip -dc '%s'", ctss_f), select = c("chr", "pos", "strand", COL_A1, COL_SG, "n_total"))
cand <- ctss[, .(n_a1 = sum(get(COL_A1)), n_sg = sum(get(COL_SG)), n_total = sum(n_total)),
             by = .(chr, pos, strand)][n_total >= MIN_DEPTH]
rm(ctss)

genome <- readDNAStringSet(fasta_f)
names(genome) <- sub("\\s.*", "", names(genome))
CTX_LEN <- 1L + ARUN_CAP        # -1..+5
cand <- cand[chr %in% names(genome)]
cand[, ctx := NA_character_]
for (ch in unique(cand$chr)) {
  s <- as.character(genome[[ch]]); n <- nchar(s)
  i <- which(cand$chr == ch); plus <- cand$strand[i] == "+"
  from <- ifelse(plus, cand$pos[i] - 1L, cand$pos[i] - ARUN_CAP + 1L); to <- from + CTX_LEN - 1L
  ok <- from >= 1L & to <= n
  x <- rep(NA_character_, length(i)); x[ok] <- toupper(substring(s, from[ok], to[ok]))
  m <- ok & !plus
  if (any(m)) x[m] <- as.character(reverseComplement(DNAStringSet(x[m])))
  cand[i, ctx := x]
}
cand <- cand[!is.na(ctx) & substr(ctx, 1, 1) %in% c("C", "T") & substr(ctx, 2, 2) == "A"]
# A run from +1; regexpr rather than sub() with a back-reference, so the length
# does not depend on the regex engine's handling of greedy groups.
cand[, arun := pmin(attr(regexpr("^A+", substr(ctx, 2, CTX_LEN)), "match.length"), ARUN_CAP)]
cand[, `:=`(m1 = substr(ctx, 1, 1), g3 = substr(ctx, 4, 4) == "G",
            dbin = floor(log2(n_total)), share_a1 = n_a1 / n_total, share_sg = n_sg / n_total)]
cand[, stratum := paste(m1, arun, dbin, sep = "|")]
cand[, share_bin := cut(share_a1, SHARE_BREAKS, SHARE_LABELS)]
nuc <- cand[chr != "chrM"]; mito <- cand[chr == "chrM"]
msg("candidates: nuclear %d (single-A clip >= %.0f %%: %d), chrM %d", nrow(nuc), 100 * MAIN_T,
    nuc[share_a1 >= MAIN_T, .N], nrow(mito))

# ── NAD signal from the step-4 matrix ────────────────────────────────────────
man <- fread(manifest_f, select = c("run", "strain", "library", "replicate"))
mat <- fread(cmd = sprintf("gzip -dc '%s'", matrix_f))
stopifnot(all(man$run %in% names(mat)))
strains <- unique(man$strain)
sanity <- if (file.exists(sanity_f)) fread(sanity_f) else NULL

size_factors <- function(cnt) {                 # DESeq2 estimateSizeFactorsForMatrix
  lg <- rowMeans(log(cnt)); ok <- is.finite(lg)
  sub <- cnt[ok, , drop = FALSE]; lgo <- lg[ok]
  apply(sub, 2, function(c) exp(median((log(c) - lgo)[c > 0])))
}
checks <- list(); sf_tab <- list()
for (st in strains) {
  cs <- man[strain == st]; runs <- cs$run
  s5 <- fread(cmd = sprintf("gzip -dc '%s'", file.path(s5dir, sprintf("%s_S_vs_N.tsv.gz", st))),
              select = c("chr", "pos", "strand", "biotype", "meanS", "meanN"))
  ctrl <- mat[s5[biotype == "none"], on = .(chr, pos, strand), nomatch = NULL]
  sf <- size_factors(as.matrix(ctrl[, ..runs]))
  ref <- if (is.null(sanity)) rep(NA_real_, length(runs)) else
    sanity[strain == st][match(paste0("size_factor_", runs), item), step5]
  sf_tab[[st]] <- data.table(strain = st, run = runs, size_factor = sf, step5_log = ref)
  nc <- sweep(as.matrix(mat[, ..runs]), 2, sf, "/")
  isS <- cs$library == "S"
  sraw <- as.matrix(mat[, runs[isS], with = FALSE])
  mat[, paste0(c("meanS_", "meanN_", "nS_", "det_"), st) :=
        list(rowMeans(nc[, isS, drop = FALSE]), rowMeans(nc[, !isS, drop = FALSE]),
             rowSums(sraw >= 1), NA)]
  mat[, (paste0("det_", st)) := get(paste0("nS_", st)) >= 2 & get(paste0("meanS_", st)) > 2 * get(paste0("meanN_", st))]
  j <- merge(mat[, .(chr, pos, strand, s = get(paste0("meanS_", st)), n = get(paste0("meanN_", st)))],
             s5, by = c("chr", "pos", "strand"))
  rel <- max(abs(c(j$s - j$meanS, j$n - j$meanN)) / pmax(1e-9, abs(c(j$meanS, j$meanN))))
  checks[[length(checks) + 1]] <- c(sprintf("%s: normalised S/N reproduce step 5 on %d positions (max rel diff %.2g)",
                                            st, nrow(j), rel), rel < 1e-6)
  checks[[length(checks) + 1]] <- c(sprintf("%s: size factors equal sanity.tsv to 3 decimals", st),
                                    isTRUE(all(abs(round(sf, 3) - ref) < 1e-9)))
}
sf_tab <- rbindlist(sf_tab)

join_nad <- function(dt) {
  cols <- grep("^(meanS|meanN|nS|det)_", names(mat), value = TRUE)
  x <- mat[, c("chr", "pos", "strand", cols), with = FALSE][dt, on = .(chr, pos, strand)]
  for (st in strains) {
    x[is.na(get(paste0("meanS_", st))), c(paste0(c("meanS_", "meanN_", "nS_"), st), paste0("det_", st)) :=
        list(0, 0, 0L, FALSE)]
    x[, (paste0("log2S_", st)) := log2(get(paste0("meanS_", st)) + 1)]
  }
  x
}
nuc <- join_nad(nuc); mito <- join_nad(mito)

# ── synthetic modes ──────────────────────────────────────────────────────────
set.seed(SEED)
if (NULL_MODE) {                       # break any link between clip share and NAD
  for (st in strains) for (v in c("det_", "log2S_")) {
    col <- paste0(v, st)
    nuc[, (col) := get(col)[sample.int(.N)], by = stratum]
  }
  msg("synthetic: NAD values shuffled within strata")
}
if (PLANT > 0) {
  for (st in strains) {
    i <- which(nuc$share_a1 >= MAIN_T & runif(nrow(nuc)) < PLANT)
    nuc[i, (paste0("det_", st)) := TRUE]
    nuc[i, (paste0("log2S_", st)) := get(paste0("log2S_", st)) + 3]
  }
  msg("synthetic: %.0f %% of single-A-clip-rich positions planted as NAD-detected", 100 * PLANT)
}

# ── stratified label permutation ─────────────────────────────────────────────
# lab: TRUE for the test set; ctl: TRUE for eligible controls; value: numeric.
perm_test <- function(dt, lab, ctl, value) {
  d <- data.table(stratum = dt$stratum, lab = lab, ctl = ctl & !lab, v = as.numeric(value))[lab | ctl]
  st <- d[, .(k = sum(lab), n_ctl = sum(ctl)), by = stratum]
  dropped <- st[n_ctl == 0, sum(k)]
  keep <- st[k > 0 & n_ctl > 0, stratum]
  d <- d[stratum %in% keep]
  obs <- d[lab == TRUE, mean(v)]
  # control mean reweighted to the stratum composition of the test set
  w <- d[, .(k = sum(lab), cm = mean(v[ctl])), by = stratum]
  ctl_mean <- w[, sum(k * cm) / sum(k)]
  pools <- split(d$v, d$stratum); ks <- w[match(names(pools), stratum), k]
  ktot <- sum(ks)
  null <- vapply(seq_len(N_DRAWS), function(b)
    sum(unlist(mapply(function(p, k) p[sample.int(length(p), k)], pools, ks, SIMPLIFY = FALSE))) / ktot,
    numeric(1))
  med <- median(null)
  data.table(n_test = ktot, n_dropped = dropped, n_ctl = d[, sum(ctl)], observed = obs, control = ctl_mean,
             ratio = obs / ctl_mean, null_lo = quantile(null, 0.025), null_med = med,
             null_hi = quantile(null, 0.975),
             p = (1 + sum(abs(null - med) >= abs(obs - med))) / (N_DRAWS + 1))
}

# Test 1
res1 <- rbindlist(lapply(strains, function(st) rbindlist(lapply(THRESHOLDS, function(t) {
  lab <- nuc$share_a1 >= t; ctl <- nuc$share_a1 < t
  rbind(perm_test(nuc, lab, ctl, nuc[[paste0("det_", st)]])[, statistic := "detected_share"],
        perm_test(nuc, lab, ctl, nuc[[paste0("log2S_", st)]])[, statistic := "mean_log2S"])[
          , `:=`(strain = st, threshold = t)]
}))))
setcolorder(res1, c("strain", "threshold", "statistic"))

strata <- nuc[, .(clip_rich = sum(share_a1 >= MAIN_T), controls = sum(share_a1 < MAIN_T)), by = stratum][order(stratum)]

# Test 2. The clip share correlates with CAGE depth (a non-zero share of a few
# percent needs a deep position), and detection rises with depth, so each
# position also gets the detection rate expected from its stratum (all
# candidates of the same -1 base, A-run and depth bin), and the dose response is
# reported as observed / expected as well as raw.
dose <- rbindlist(lapply(strains, function(st) {
  x <- nuc[, .(arun, share_bin, stratum, det = get(paste0("det_", st)), l2 = get(paste0("log2S_", st)))]
  x[, expd := mean(det), by = stratum]
  x[, .(strain = st, n = .N, detected_share = mean(det), expected_share = mean(expd),
        obs_over_exp = mean(det) / mean(expd), mean_log2S = mean(l2)), by = .(arun, share_bin)]
}))
setorder(dose, strain, arun, share_bin)

# Test 3: A-run 2 only, stratum = -1 base x depth bin
a2 <- nuc[arun == 2]
a2[, stratum := paste(m1, dbin, sep = "|")]
lab <- a2$share_a1 >= MAIN_T
res3 <- rbind(perm_test(a2, lab, a2$share_sg > 0.5, a2$g3)[, controls := sprintf("%s-majority", END_TYPE_CANONICAL)],
              perm_test(a2, lab, !lab, a2$g3)[, controls := "all other A-run-2"])
setcolorder(res3, "controls")

# Test 4
chrm <- rbindlist(lapply(strains, function(st)
  rbind(mito[share_a1 >= MAIN_T, .(strain = st, set = "chrM, clip-rich", n = .N,
                                   detected_share = mean(get(paste0("det_", st))))],
        mito[share_a1 < MAIN_T, .(strain = st, set = "chrM, other", n = .N,
                                  detected_share = mean(get(paste0("det_", st))))],
        nuc[share_a1 >= MAIN_T, .(strain = st, set = "nuclear, clip-rich", n = .N,
                                  detected_share = mean(get(paste0("det_", st))))],
        nuc[share_a1 < MAIN_T, .(strain = st, set = "nuclear, other", n = .N,
                                 detected_share = mean(get(paste0("det_", st))))])))

# ── outputs ──────────────────────────────────────────────────────────────────
keep_cols <- c("chr", "pos", "strand", "ctx", "m1", "arun", "dbin", "stratum", "n_total", "n_a1", "n_sg",
               "share_a1", "share_sg", "share_bin", "g3", grep("^(meanS|meanN|nS|det|log2S)_", names(nuc), value = TRUE))
fwrite(nuc[, ..keep_cols], file.path(outdir, "positions.tsv.gz"), sep = "\t")
fwrite(res1,   out_main, sep = "\t")
fwrite(dose,   file.path(outdir, "dose_response.tsv"), sep = "\t")
fwrite(res3,   file.path(outdir, "motif.tsv"), sep = "\t")
fwrite(chrm,   file.path(outdir, "chrM_control.tsv"), sep = "\t")
fwrite(strata, file.path(outdir, "strata.tsv"), sep = "\t")
fwrite(sf_tab, file.path(outdir, "size_factors.tsv"), sep = "\t")

options(width = 200)
print(res1[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
print(res3[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
print(chrm)

# ── figure ───────────────────────────────────────────────────────────────────
# A: at the main threshold, how often NAD captureSeq sees something at nuclear
#    1Sa-rich positions, at matched ordinary nuclear positions (median and 95 %
#    range of the random matched sets) and at chrM 1Sa-rich positions, the
#    positive control. B: the same comparison at every threshold, as a ratio to
#    the matched median (the style of figure 06). C: YAAG. The dose response by
#    A run and clip-share bin stays in dose_response.tsv.
COL_HI <- "#e63946"
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.5))
strain_name <- function(s) factor(fifelse(s == "WT", "WT", "triple mutant"), c("WT", "triple mutant"))
rich_lab <- sprintf("%s-rich", TYPE_A1)

m1 <- res1[statistic == "detected_share" & threshold == MAIN_T]
pa <- rbind(
  m1[, .(strain, set = sprintf("nuclear,\n%s", rich_lab), share = observed, lo = NA_real_, hi = NA_real_,
         n = n_test, p)],
  m1[, .(strain, set = "nuclear,\nmatched ordinary", share = null_med, lo = null_lo, hi = null_hi,
         n = NA_integer_, p = NA_real_)],
  chrm[set == "chrM, clip-rich", .(strain, set = sprintf("chrM,\n%s", rich_lab), share = detected_share,
                                   lo = NA_real_, hi = NA_real_, n, p = NA_real_)])
pa[, set := factor(set, unique(set))]
pa[, strain := strain_name(strain)]
pa[, lab := fifelse(is.na(n), sprintf("%.0f %%", 100 * share), sprintf("%.0f %%\n(%d)", 100 * share, n))]
pa[, fill := c("rich", "matched", "chrM")[as.integer(set)]]
pa[, py := 100 * max(c(share, hi)[as.integer(set) < 3], na.rm = TRUE) + 14, by = strain]
pA <- ggplot(pa, aes(set, 100 * share)) +
  geom_col(aes(fill = fill), width = 0.6) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.2, na.rm = TRUE) +
  geom_text(aes(y = 100 * pmax(share, hi, na.rm = TRUE), label = lab), vjust = -0.3, size = 2.5, lineheight = 0.9) +
  geom_text(data = pa[!is.na(p)], aes(x = 1.5, y = py, label = sprintf("p = %.2g", p)), size = 2.5) +
  scale_fill_manual(values = c(rich = COL_HI, matched = "grey65", chrM = "grey25"), guide = "none") +
  scale_y_continuous(limits = c(0, 92), breaks = seq(0, 80, 20)) +
  facet_wrap(~ strain) +
  labs(x = NULL, y = "positions where NAD captureSeq\nsees a signal (%)",
       title = sprintf("A. Does NAD captureSeq see CAGE %s-rich positions? (%s >= %.0f %% of CAGE reads)",
                       TYPE_A1, TYPE_A1, 100 * MAIN_T),
       subtitle = paste0("Signal: S reads in >= 2 of 3 replicates and S > 2 x N. Matched ordinary: same -1 base, A-run length\n",
                         "and CAGE depth; bar = median, whisker = 95 % of 1000 random matched sets. chrM = positive control.")) +
  th

r1 <- res1[statistic == "detected_share"]
# a null range that reaches 0 is clipped at the axis minimum, which a log axis needs
r1[, `:=`(strain = strain_name(strain), rel = observed / null_med, lo = pmax(null_lo / null_med, 0.5),
          hi = null_hi / null_med,
          thr = factor(sprintf(">= %.0f %%", 100 * threshold), rev(sprintf(">= %.0f %%", 100 * THRESHOLDS))))]
pB <- ggplot(r1, aes(y = thr)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 3, colour = "grey80") +
  geom_point(aes(x = rel, colour = p < 0.05), size = 2.6) +
  geom_text(aes(x = 6, label = sprintf("p = %.2g", p)), hjust = 1, size = 2.5, colour = "grey30") +
  scale_colour_manual(values = c(`TRUE` = COL_HI, `FALSE` = "black"), guide = "none") +
  scale_x_log10(limits = c(0.5, 6), breaks = c(0.5, 1, 2, 4), labels = c("×0.5", "×1", "×2", "×4")) +
  facet_wrap(~ strain) +
  labs(x = "signal share, 1Sa-rich / matched ordinary (log)",
       y = sprintf("%s threshold", TYPE_A1),
       title = "B. Same comparison at every 1Sa threshold",
       subtitle = "Dot: 1Sa-rich positions (red: p < 0.05).\nGrey: 95 % of random matched sets. ×1 = no difference.") +
  th + theme(panel.grid.major.y = element_blank())

r3 <- rbind(res3[1, .(set = rich_lab, share = observed, n = n_test)],
            res3[, .(set = controls, share = control, n = n_ctl)])
r3[, set := factor(set, unique(set))]
pC <- ggplot(r3, aes(set, 100 * share)) +
  geom_col(fill = c(COL_HI, "grey65", "grey65"), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f %%\n(of %s)", 100 * share, prettyNum(n, big.mark = ","))),
            vjust = -0.2, size = 2.5, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 55)) +
  labs(x = "nuclear Y|A positions with an A run of exactly 2 (Y|AA?)", y = "G at +3 (%)",
       title = "C. The nuclear NAD motif YAAG: is there a G at +3?",
       subtitle = sprintf("p = %.2g, %s-rich against %s-majority positions.", res3$p[1], TYPE_A1,
                          END_TYPE_CANONICAL)) +
  th

ggsave(file.path(figdir, "08_1Sa_positions.pdf"),
       pA / (pB | pC + plot_layout()) + plot_layout(heights = c(1.1, 1)),
       width = 10, height = 8)

for (k in checks) msg("[%s] %s", if (k[2] == "TRUE") "PASS" else "FAIL", k[1])
msg("[INFO] nuclear candidates %d, single-A-clip-rich (>= %.0f %%) %d", nrow(nuc), 100 * MAIN_T,
    nuc[share_a1 >= MAIN_T, .N])
sessionInfo()
