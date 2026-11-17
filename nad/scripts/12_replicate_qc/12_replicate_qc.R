#!/usr/bin/env Rscript
# NAD branch step 12: replicate concordance of NAD captureSeq libraries.
#
# For each strain, on the fixed step-5 tested position set, the script compares
# every pair of runs on DESeq2 VST values (blind), relative to a Poisson ceiling:
# the correlation two runs of the same library would reach at their size factors
# if they differed only by counting noise. It adds S_i/N_i enrichment
# reproducibility, PCA, sample distances, Cook's distance, the positions that
# drive each S run's disagreement, and leave-one-out DESeq2 refits, and applies
# a three-part outlier rule to every run. Step 5 is not changed.
# Spec: docs/specs/2026-11-16-nad-step12-replicate-qc.md
#
# Outputs (in --outdir): read_classes.tsv pairs.tsv enrich_pairs.tsv pca.tsv
#   cooks.tsv drivers_<strain>.tsv loo.tsv verdict.tsv
#   and in --figdir: 12_pairs_<strain>.pdf 12_pca.pdf 12_distance.pdf 12_loo.pdf
# --self-test writes the same tables for one synthetic strain into
#   <outdir>/selftest/ and checks that the rule flags only the noised run.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--matrix",      type = "character", default = "nad/results/04_ctss/matrix_reads.tsv.gz"),
  make_option("--manifest",    type = "character", default = "nad/scripts/00_download/priority1_runs.tsv"),
  make_option("--step5-dir",   type = "character", default = "nad/results/05_deseq"),
  make_option("--outdir",      type = "character", default = "nad/results/12_replicate_qc"),
  make_option("--figdir",      type = "character", default = "nad/figures"),
  make_option("--min-reads",   type = "integer",   default = 10L),
  make_option("--min-samples", type = "integer",   default = 3L),
  make_option("--padj",        type = "double",    default = 0.05),
  make_option("--lfc",         type = "double",    default = 1),
  make_option("--draws",       type = "integer",   default = 200L),
  make_option("--margin",      type = "double",    default = 0.05),
  make_option("--seed",        type = "integer",   default = 1L),
  make_option("--self-test",   action = "store_true", default = FALSE)
)))
SELF <- opt$`self-test`
outdir <- if (SELF) file.path(opt$outdir, "selftest") else opt$outdir
figdir <- if (SELF) outdir else opt$figdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))
set.seed(opt$seed)
CLASSES <- c("nuclear", "chrM", "ncRNA")

# ── inputs ───────────────────────────────────────────────────────────────────
man <- fread(opt$manifest, select = c("run", "strain", "library", "replicate"))
# fread reads .gz itself only with R.utils, which r_nad lacks.
mat <- fread(cmd = sprintf("gzip -dcf '%s'", opt$matrix))
stopifnot(all(man$run %in% names(mat)))
read_s5 <- function(st) fread(cmd = sprintf("gzip -dcf '%s'",
                                            file.path(opt$`step5-dir`, sprintf("%s_S_vs_N.tsv.gz", st))))
pos_class <- function(chr, biotype) ifelse(biotype != "none", "ncRNA", ifelse(chr == "chrM", "chrM", "nuclear"))

prefilter <- function(cnt) rowSums(cnt >= opt$`min-reads`) >= opt$`min-samples`

# Step-5 fit: design ~ library, N reference, size factors over biotype "none".
fit <- function(cnt, lib, ctrl) {
  cd <- data.frame(library = factor(lib, levels = c("N", "S")), row.names = colnames(cnt))
  dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ library)
  dds <- estimateSizeFactors(dds, controlGenes = which(ctrl))
  DESeq(dds, quiet = TRUE)
}
calls <- function(dds, ctrl) {
  res <- results(dds, contrast = c("library", "S", "N"), alpha = opt$padj)
  data.table(id = rownames(dds), lfc = res$log2FoldChange,
             enriched = ctrl & !is.na(res$padj) & res$padj < opt$padj & res$log2FoldChange > opt$lfc)
}

# The VST of DESeq2 getVarianceStabilizedData, rebuilt as a fixed function of
# normalised counts so that real and simulated counts go through exactly the
# same transform. Parametric fit: closed form. Local fit (DESeq2 substitutes it
# when y = a/x + b does not fit): DESeq2's own numerical integration on its own
# grid (0 to the real maximum, 1000 points), continued linearly in asinh(q)
# beyond the grid end with the slope at the end, so that simulated counts above
# the real maximum are not extrapolated by the cubic spline. A finer grid was
# rejected because it departs from DESeq2's values by up to 0.08. Both branches
# are checked against DESeq2's output on the real data.
vst_fun <- function(df, nc, sf) {
  if (identical(attr(df, "fitType"), "parametric")) {
    cf <- attr(df, "coefficients"); a0 <- cf[["asymptDisp"]]; a1 <- cf[["extraPois"]]
    return(function(q) log((1 + a1 + 2 * a0 * q + 2 * sqrt(a0 * q * (1 + a1 + a0 * q))) / (4 * a0)) / log(2))
  }
  if (!identical(attr(df, "fitType"), "local")) stop("unsupported dispersion fitType: ", attr(df, "fitType"))
  xg <- sinh(seq(asinh(0), asinh(max(nc)), length.out = 1000))[-1]
  xim <- mean(1 / sf)
  ig <- 1 / sqrt(df(xg) * xg^2 + xim * xg)
  splf <- splinefun(asinh((xg[-1] + xg[-length(xg)]) / 2),
                    cumsum((xg[-1] - xg[-length(xg)]) * (ig[-1] + ig[-length(ig)]) / 2))
  rm_ <- rowMeans(nc); h1 <- quantile(rm_, 0.95); h2 <- quantile(rm_, 0.999)
  eta <- (log2(h2) - log2(h1)) / (splf(asinh(h2)) - splf(asinh(h1)))
  xi <- log2(h1) - eta * splf(asinh(h1))
  aend <- asinh(max(nc)); fend <- splf(aend); dend <- splf(aend, deriv = 1)
  function(q) {
    a <- asinh(q); hi <- a > aend
    r <- numeric(length(a)); r[!hi] <- splf(a[!hi]); r[hi] <- fend + dend * (a[hi] - aend)
    r <- eta * r + xi; if (is.matrix(q)) dim(r) <- dim(q); r
  }
}

# ── one strain ───────────────────────────────────────────────────────────────
analyse <- function(st, cnt, pos, cs, s5 = NULL) {
  runs <- cs$run; lab <- paste0(cs$library, "_", cs$replicate); names(lab) <- runs
  colnames(cnt) <- runs
  cls <- pos_class(pos$chr, pos$biotype); ctrl <- pos$biotype == "none"
  rownames(cnt) <- pos$key
  msg("%s: %d positions (%s)", st, nrow(cnt),
      paste(sprintf("%s %d", CLASSES, sapply(CLASSES, function(k) sum(cls == k))), collapse = ", "))

  dds <- fit(cnt, cs$library, ctrl)
  sf <- sizeFactors(dds); nc <- counts(dds, normalized = TRUE)
  isS <- cs$library == "S"
  if (!is.null(s5)) {
    pk <- pos$key; s5m <- s5[match(pk, s5$key)]
    rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
    dS <- rel(rowMeans(nc[, isS]), s5m$meanS); dN <- rel(rowMeans(nc[, !isS]), s5m$meanN)
    msg("%s: recomputed meanS / meanN vs step 5, max rel diff %.2g / %.2g", st, dS, dN)
    if (dS > 1e-6 || dN > 1e-6) stop(sprintf("%s: normalisation does not reproduce step 5", st))
  }

  # blind VST: dispersions re-estimated under ~ 1
  ddb <- dds; design(ddb) <- ~ 1
  ddb <- estimateDispersions(ddb, fitType = "parametric", quiet = TRUE)
  dfun <- dispersionFunction(ddb)
  ref <- unname(assay(varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")))
  V <- vst_fun(dfun, nc, sf)
  vs <- V(nc)
  dv <- max(abs(vs - ref))
  msg("%s: rebuilt VST (%s fit) vs DESeq2 vst, max abs diff %.2g", st, attr(dfun, "fitType"), dv)
  if (dv > 1e-6) stop(sprintf("%s: rebuilt VST does not match DESeq2", st))
  dimnames(vs) <- list(rownames(cnt), runs)

  # ceiling: Poisson draws from the pooled profile of each library
  sim_r <- list()
  for (lb in c("S", "N")) {
    rr <- runs[cs$library == lb]; mu <- rowMeans(nc[, rr, drop = FALSE])
    prs <- combn(rr, 2, simplify = FALSE)
    draws <- replicate(opt$draws, {
      sv <- sapply(rr, function(r) V(rpois(length(mu), mu * sf[[r]]) / sf[[r]]))
      unlist(lapply(prs, function(p) sapply(CLASSES, function(k) {
        i <- cls == k; if (sum(i) < 3) NA_real_ else cor(sv[i, p[1]], sv[i, p[2]]) })))
    })
    nm <- unlist(lapply(prs, function(p) paste(p[1], p[2], CLASSES, sep = "|")))
    sim_r[[lb]] <- setNames(apply(matrix(draws, ncol = opt$draws), 1, median, na.rm = TRUE), nm)
  }
  sim_r <- unlist(unname(sim_r))

  prs <- combn(runs, 2, simplify = FALSE)
  pairs <- rbindlist(lapply(prs, function(p) rbindlist(lapply(CLASSES, function(k) {
    i <- cls == k
    ce <- sim_r[paste(p[1], p[2], k, sep = "|")]
    data.table(strain = st, class = k, run_a = p[1], run_b = p[2], lab_a = lab[[p[1]]], lab_b = lab[[p[2]]],
               n = sum(i),
               pearson_vst   = if (sum(i) >= 3) cor(vs[i, p[1]], vs[i, p[2]]) else NA_real_,
               spearman_norm = if (sum(i) >= 3) cor(nc[i, p[1]], nc[i, p[2]], method = "spearman") else NA_real_,
               ceiling_r = if (is.na(ce)) NA_real_ else unname(ce))
  }))))
  pairs[, ratio := pearson_vst / ceiling_r]

  # S_i / N_i enrichment reproducibility (own denominator per replicate)
  nuc <- cls == "nuclear"
  reps <- intersect(cs[library == "S"]$replicate, cs[library == "N"]$replicate)
  lr <- sapply(reps, function(r) {
    s <- cs[library == "S" & replicate == r]$run; n <- cs[library == "N" & replicate == r]$run
    log2((nc[nuc, s] + 1) / (nc[nuc, n] + 1)) })
  enrich <- rbindlist(lapply(combn(reps, 2, simplify = FALSE), function(p)
    data.table(strain = st, rep_a = p[1], rep_b = p[2], n = sum(nuc), r = cor(lr[, p[1]], lr[, p[2]]))))

  # PCA and distances on nuclear VST
  pc <- prcomp(t(vs[nuc, ]))
  pve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pca <- data.table(strain = st, run = runs, lab = lab[runs], library = cs$library,
                    pc$x[, 1:min(4, ncol(pc$x)), drop = FALSE])
  pca[, pve := paste(sprintf("%.1f", pve[1:min(4, length(pve))]), collapse = ",")]
  D <- as.matrix(dist(t(vs[nuc, ])))
  dist_dt <- data.table(strain = st, a = lab[rownames(D)][row(D)], b = lab[colnames(D)][col(D)], d = as.vector(D))

  # Cook's distance: positions where a run carries the largest Cook's above the cutoff
  ck <- assays(dds)[["cooks"]]
  np <- ncol(model.matrix(design(dds), colData(dds)))
  cut <- qf(0.99, np, ncol(dds) - np)               # DESeq2's default cooksCutoff
  top <- max.col(ck, ties.method = "first"); big <- apply(ck, 1, max) > cut
  cooks <- data.table(strain = st, run = runs, lab = lab[runs], cutoff = cut,
                      n_max_cooks = sapply(seq_along(runs), function(j) sum(big & top == j, na.rm = TRUE)),
                      n_max_cooks_nuclear = sapply(seq_along(runs), function(j) sum(big & top == j & nuc, na.rm = TRUE)))

  # drivers: residual of each S run against the mean of the other two
  sr <- runs[isS]
  drivers <- rbindlist(lapply(sr, function(r) {
    o <- rowMeans(vs[, setdiff(sr, r), drop = FALSE]); res <- residuals(lm(vs[, r] ~ o))
    d <- data.table(run = r, lab = lab[[r]], pos[, .(chr, pos, strand, base_p1, ctx, biotype)], class = cls,
                    vst_run = vs[, r], vst_others = o, residual = res)
    d[order(-abs(residual))][1:min(.N, 200)]
  }))
  msg("%s: top-200 residual positions by class:", st)
  print(dcast(drivers[, .N, by = .(lab, class)], lab ~ class, value.var = "N", fill = 0))
  print(drivers[, .(share_positive = mean(residual > 0)), by = lab])

  # leave-one-out refits
  # The refits keep the full-fit position set. Re-applying ">= 10 reads in >= 3
  # runs" to 5 runs would drop every position seen only in S once an S run is
  # removed (N is near-empty at mRNA starts), which collapsed the calls ~15-fold
  # in the self-test.
  full <- calls(dds, ctrl)
  loo <- rbindlist(lapply(runs, function(r) {
    keep_r <- runs != r
    sub <- cnt[, keep_r, drop = FALSE]; k <- rowSums(sub) > 0
    d <- fit(sub[k, ], cs$library[keep_r], ctrl[k]); cl <- calls(d, ctrl[k])
    a <- full[enriched == TRUE]$id; b <- cl[enriched == TRUE]$id
    sh <- merge(full[, .(id, f = lfc)], cl[, .(id, l = lfc)], by = "id")
    data.table(strain = st, dropped_run = r, lab = lab[[r]], tested = sum(k), enriched = length(b),
               full_enriched = length(a), kept_of_full = length(intersect(a, b)) / length(a),
               gained = length(setdiff(b, a)), jaccard = length(intersect(a, b)) / length(union(a, b)),
               spearman_lfc_vs_full = cor(sh$f, sh$l, method = "spearman", use = "complete.obs"))
  }))
  msg("%s: leave-one-out done", st)

  # outlier rule
  pn <- pairs[class == "nuclear"]
  verdict <- rbindlist(lapply(runs, function(r) {
    lb <- cs[run == r]$library; mates <- setdiff(runs[cs$library == lb], r)
    pr <- function(a, b) pn[(run_a == a & run_b == b) | (run_a == b & run_b == a)]$ratio
    own <- sapply(mates, function(m) pr(r, m)); rest <- if (length(mates) == 2) pr(mates[1], mates[2]) else NA_real_
    c1 <- isTRUE(all(own <= rest - opt$margin))
    tot <- sapply(c(r, mates), function(x) sum(D[x, setdiff(c(r, mates), x)]))
    c2 <- names(which.max(tot)) == r
    le <- loo[match(c(r, mates), dropped_run)]$enriched
    c3 <- le[1] > max(le[-1])
    data.table(strain = st, run = r, lab = lab[[r]],
               ratio_with_mates = paste(sprintf("%.3f", own), collapse = ","), ratio_other_pair = rest,
               c1_low_ratio = c1, distance_sum = tot[[1]], c2_most_distant = c2,
               loo_enriched = le[1], loo_enriched_mates = paste(le[-1], collapse = ","), c3_drop_gains = c3,
               outlier = c1 && c2 && c3)
  }))

  read_classes <- rbindlist(lapply(runs, function(r) data.table(
    strain = st, run = r, lab = lab[[r]], class = CLASSES,
    reads = sapply(CLASSES, function(k) sum(cnt[cls == k, r])))))[, share := reads / sum(reads), by = run]

  # figures
  long <- rbindlist(lapply(prs, function(p) data.table(
    pair = paste(lab[[p[1]]], "vs", lab[[p[2]]]), x = vs[nuc, p[1]], y = vs[nuc, p[2]])))
  ann <- pn[, .(pair = paste(lab_a, "vs", lab_b),
                txt = ifelse(is.na(ceiling_r), sprintf("r = %.3f", pearson_vst),
                             sprintf("r = %.3f\nceiling %.3f", pearson_vst, ceiling_r)))]
  p <- ggplot(long, aes(x, y)) + geom_abline(linewidth = 0.3, colour = "grey50") +
    geom_point(size = 0.3, alpha = 0.3) + facet_wrap(~ pair, ncol = 5) +
    geom_text(data = ann, aes(-Inf, Inf, label = txt), hjust = -0.05, vjust = 1.1, size = 2.6) +
    labs(title = sprintf("%s: VST of nuclear unmasked positions (n = %d)", st, sum(nuc)),
         x = "VST, first run", y = "VST, second run") + theme_bw(base_size = 9)
  ggsave(file.path(figdir, sprintf("12_pairs_%s.pdf", st)), p, width = 12, height = 7.5)

  list(pairs = pairs, enrich = enrich, pca = pca, dist = dist_dt, cooks = cooks,
       drivers = drivers, loo = loo, verdict = verdict, read_classes = read_classes)
}

# ── run ──────────────────────────────────────────────────────────────────────
inputs <- list()
if (!SELF) {
  for (st in unique(man$strain)) {
    cs <- man[strain == st]; runs <- cs$run
    s5 <- read_s5(st); s5[, key := paste(chr, pos, strand, sep = ":")]
    m <- mat[, c("chr", "pos", "strand", "base_p1", runs), with = FALSE]
    keep <- prefilter(as.matrix(m[, ..runs]))
    m <- m[keep]; m[, key := paste(chr, pos, strand, sep = ":")]
    if (!setequal(m$key, s5$key)) stop(sprintf("%s: pre-filter set differs from step 5 (%d vs %d)", st, nrow(m), nrow(s5)))
    m <- m[match(s5$key, key)]
    pos <- cbind(m[, .(key, chr, pos, strand, base_p1)], s5[, .(ctx, biotype)])
    inputs[[st]] <- list(cnt = as.matrix(m[, ..runs]), pos = pos, cs = cs, s5 = s5)
  }
} else {
  # Synthetic strain: Poisson draws from the triple-mutant pooled profile at the
  # real size factors, with log-normal noise (sdlog 0.5) on S rep1.
  st <- "dxo1_rai1_npy1"; cs <- man[strain == st]; runs <- cs$run
  s5 <- read_s5(st); s5[, key := paste(chr, pos, strand, sep = ":")]
  m <- mat[, c("chr", "pos", "strand", "base_p1", runs), with = FALSE][, key := paste(chr, pos, strand, sep = ":")]
  m <- m[match(s5$key, key)]
  cnt <- as.matrix(m[, ..runs]); ctrl <- s5$biotype == "none"
  d0 <- fit(cnt, cs$library, ctrl); sf <- sizeFactors(d0); nc <- counts(d0, normalized = TRUE)
  syn <- sapply(seq_along(runs), function(j) {
    mu <- rowMeans(nc[, cs$library == cs$library[j], drop = FALSE]) * sf[[j]]
    if (cs$library[j] == "S" && cs$replicate[j] == "rep1") mu <- mu * rlnorm(length(mu), 0, 0.5)
    rpois(length(mu), mu) })
  keep <- prefilter(syn)
  cs_syn <- copy(cs)[, `:=`(run = paste0("SYN_", library, "_", replicate), strain = "synthetic")]
  colnames(syn) <- cs_syn$run
  pos <- cbind(m[, .(key, chr, pos, strand, base_p1)], s5[, .(ctx, biotype)])[keep]
  inputs[["synthetic"]] <- list(cnt = syn[keep, ], pos = pos, cs = cs_syn, s5 = NULL)
}

out <- lapply(names(inputs), function(st) with(inputs[[st]], analyse(st, cnt, pos, cs, s5)))
tab <- function(k) rbindlist(lapply(out, `[[`, k))
fw <- function(x, f) fwrite(x, file.path(outdir, f), sep = "\t")
fw(tab("read_classes"), "read_classes.tsv")
fw(tab("pairs"), "pairs.tsv")
fw(tab("enrich"), "enrich_pairs.tsv")
fw(tab("pca"), "pca.tsv")
fw(tab("cooks"), "cooks.tsv")
fw(tab("loo"), "loo.tsv")
for (i in seq_along(out)) fw(out[[i]]$drivers, sprintf("drivers_%s.tsv", names(inputs)[i]))

pca <- tab("pca")
p <- ggplot(pca, aes(PC1, PC2, colour = library, label = lab)) + geom_point(size = 2) +
  geom_text(size = 2.6, vjust = -0.8, show.legend = FALSE) + facet_wrap(~ strain, scales = "free") +
  scale_colour_manual(values = c(S = "#e63946", N = "#457b9d")) +
  labs(title = "PCA of VST, nuclear unmasked positions") + theme_bw(base_size = 10)
ggsave(file.path(figdir, "12_pca.pdf"), p, width = 9, height = 4.5)
p <- ggplot(tab("dist"), aes(a, b, fill = d)) + geom_tile() + geom_text(aes(label = round(d)), size = 2.4) +
  facet_wrap(~ strain, scales = "free") + scale_fill_gradient(low = "#ffffff", high = "#1d3557", name = "distance") +
  labs(x = NULL, y = NULL, title = "Euclidean distance of VST, nuclear unmasked positions") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figdir, "12_distance.pdf"), p, width = 9, height = 4.5)
loo <- tab("loo")
p <- ggplot(loo, aes(lab, enriched, fill = substr(lab, 1, 1))) + geom_col() +
  geom_hline(aes(yintercept = full_enriched), linetype = "dashed", linewidth = 0.3) +
  geom_text(aes(label = sprintf("J=%.2f", jaccard)), vjust = -0.3, size = 2.4) +
  facet_wrap(~ strain, scales = "free_y") + scale_fill_manual(values = c(S = "#e63946", N = "#457b9d"), guide = "none") +
  labs(x = "dropped run", y = "enriched unmasked positions",
       title = "Leave-one-out DESeq2 (dashed: full fit; J: Jaccard with full calls)") + theme_bw(base_size = 10)
ggsave(file.path(figdir, "12_loo.pdf"), p, width = 9, height = 4)

verdict <- tab("verdict"); fw(verdict, "verdict.tsv")
cat("\n== nuclear unmasked pairs\n"); print(tab("pairs")[class == "nuclear", .(strain, lab_a, lab_b, n,
  pearson_vst = round(pearson_vst, 3), spearman = round(spearman_norm, 3), ceiling = round(ceiling_r, 3), ratio = round(ratio, 3))])
cat("\n== enrichment pairs\n"); print(tab("enrich"))
cat("\n== Cook's\n"); print(tab("cooks"))
cat("\n== leave-one-out\n"); print(loo)
cat("\n== verdict\n"); print(verdict)

if (SELF) {
  pn <- tab("pairs")[class == "nuclear" & substr(lab_a, 1, 1) == substr(lab_b, 1, 1)]
  noisy <- "SYN_S_rep1"
  clean <- pn[run_a != noisy & run_b != noisy]
  bad <- pn[run_a == noisy | run_b == noisy]
  ok <- c(sprintf("clean same-library pairs have ratio in [0.97, 1.03]: %s",
                  paste(sprintf("%.3f", clean$ratio), collapse = " ")),
          sprintf("noised run has ratio <= min clean ratio - %.2f with both partners: %s (min clean %.3f)",
                  opt$margin, paste(sprintf("%.3f", bad$ratio), collapse = " "), min(clean$ratio)),
          sprintf("outlier rule flags only the noised run: %s", paste(verdict[outlier == TRUE]$run, collapse = ",")))
  pass <- c(all(clean$ratio >= 0.97 & clean$ratio <= 1.03), all(bad$ratio <= min(clean$ratio) - opt$margin),
            identical(verdict[outlier == TRUE]$run, noisy))
  for (k in seq_along(ok)) msg("[%s] %s", if (pass[k]) "PASS" else "FAIL", ok[k])
  if (!all(pass)) quit(save = "no", status = 1)
}

sessionInfo()
