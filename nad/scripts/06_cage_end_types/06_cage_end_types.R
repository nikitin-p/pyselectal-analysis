#!/usr/bin/env Rscript
# NAD branch step 6: CAGE end types at NAD captureSeq 5' ends.
#
# Asks whether the CAGE reads that start at NAD-enriched positions (step 5:
# unmasked, A at +1, padj < 0.05, log2FC > 1) carry a different end-type
# composition, above all a higher single-A clip share, than CAGE reads at
# matched control positions. Only Y|A starts (-1 = C or T, +1 = A) are tested,
# because a clipped A is impossible on a genomic A at -1, and controls are
# matched on the -1 base, the A-run length from +1 (capped at 5; step 1 showed
# slippage on A tracts drives the clip) and the log2 bin of the CAGE depth.
# The test is a matched permutation of read-pooled shares.
# Spec: docs/specs/2026-11-04-nad-step6-cage-end-types.md
#
# Runs: A  triple mutant, all CAGE libraries          (primary)
#       B  WT, all CAGE libraries
#       C  triple mutant, glucose library only
#       D  triple mutant, all libraries, controls only from NAD-tested,
#          non-enriched positions
#
# Usage (from project root):
#   Rscript nad/scripts/06_cage_end_types/06_cage_end_types.R \
#       [--ctss results/tss/typed_ctss.tsv.gz] [--nad-dir nad/results/05_deseq] \
#       [--fasta data/sacCer3.fa] [--outdir nad/results/06_cage_end_types] \
#       [--figdir nad/figures] [--draws 1000] [--seed 1] [--min-depth 5] \
#       [--glucose-sample yeast_cage_glu_r2] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
ctss_f    <- get_arg("--ctss",    "results/tss/typed_ctss.tsv.gz")
nad_dir   <- get_arg("--nad-dir", "nad/results/05_deseq")
fasta_f   <- get_arg("--fasta",   "data/sacCer3.fa")
outdir    <- get_arg("--outdir",  "nad/results/06_cage_end_types")
figdir    <- get_arg("--figdir",  "nad/figures")
N_DRAWS   <- as.integer(get_arg("--draws", "1000"))
SEED      <- as.integer(get_arg("--seed", "1"))
MIN_DEPTH <- as.integer(get_arg("--min-depth", "5"))
GLU       <- get_arg("--glucose-sample", "yeast_cage_glu_r2")
PADJ      <- 0.05
LFC       <- 1
ARUN_CAP  <- 5
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "..", "..", "scripts", "utils", "end_types.R"))

out_comp <- file.path(outdir, "composition.tsv")
if (file.exists(out_comp) && !force) {
  message("[skip] ", out_comp, " exists (use --force)"); quit(save = "no")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

# Types named from the vocabulary rather than as literals: the single-A clip is
# primary; the reported set adds every other A-abutting clip, M and the
# canonical type.
TYPE_A1 <- END_TYPES_REAL[END_TYPE_CLIP_BASE[END_TYPES_REAL] %in% "A" &
                            startsWith(END_TYPES_REAL, "1S")]
stopifnot(length(TYPE_A1) == 1)
TYPE_M <- END_TYPES_REAL[is.na(END_TYPE_CLIP_BASE[END_TYPES_REAL])]
stopifnot(length(TYPE_M) == 1)
REPORTED <- unique(c(TYPE_A1,
                     END_TYPES_REAL[END_TYPE_CLIP_BASE[END_TYPES_REAL] %in% "A"],
                     TYPE_M, END_TYPE_CANONICAL))
NCOLS <- end_type_cols("n_")

# ── NAD sets from step 5 ─────────────────────────────────────────────────────
nad_files <- list.files(nad_dir, pattern = "_S_vs_N\\.tsv\\.gz$", full.names = TRUE)
if (length(nad_files) == 0) stop("no *_S_vs_N.tsv.gz in ", nad_dir)
nad <- rbindlist(lapply(nad_files, function(f) {
  d <- fread(cmd = sprintf("gzip -dc '%s'", f),
             select = c("chr", "pos", "strand", "base_p1", "biotype", "log2FoldChange", "padj"))
  d[, strain := sub("_S_vs_N\\.tsv\\.gz$", "", basename(f))]
}))
nad[, enriched := !is.na(padj) & padj < PADJ & log2FoldChange > LFC]
strains <- unique(nad$strain)
MUT <- setdiff(strains, "WT")
if (length(MUT) != 1) stop("expected one non-WT strain, found: ", paste(MUT, collapse = ", "))
nad_pos <- nad[enriched & biotype == "none" & base_p1 == "A", .(chr, pos, strand, strain)]
msg("NAD positions (unmasked, A at +1, enriched): %s",
    paste(sprintf("%s %d", names(table(nad_pos$strain)), table(nad_pos$strain)), collapse = ", "))
# Every enriched position of either strain, masked or not, is kept out of the controls.
any_enriched <- unique(nad[enriched == TRUE, .(chr, pos, strand)])

# ── CAGE: pooled over all libraries and the glucose library alone ────────────
ctss <- fread(cmd = sprintf("gzip -dc '%s'", ctss_f))
stopifnot(all(NCOLS %in% names(ctss)), GLU %in% ctss$sample_id)
ctss <- ctss[chr != "chrM"]
pool <- function(d) d[, lapply(.SD, sum), by = .(chr, pos, strand), .SDcols = c(NCOLS, "n_total")]
cage <- list(all = pool(ctss), glu = pool(ctss[sample_id == GLU]))
rm(ctss)
msg("CAGE positions: all %d, %s %d", nrow(cage$all), GLU, nrow(cage$glu))

# ── context -1..+5 on the RNA strand ─────────────────────────────────────────
genome <- readDNAStringSet(fasta_f)
names(genome) <- sub("\\s.*", "", names(genome))
CTX_LEN <- 1L + ARUN_CAP
add_context <- function(dt) {
  dt[, ctx := NA_character_]
  for (ch in unique(dt$chr)) {
    s <- as.character(genome[[ch]]); n <- nchar(s)
    i <- which(dt$chr == ch)
    plus <- dt$strand[i] == "+"
    from <- ifelse(plus, dt$pos[i] - 1L, dt$pos[i] - ARUN_CAP + 1L)
    to <- from + CTX_LEN - 1L
    ok <- from >= 1L & to <= n
    x <- rep(NA_character_, length(i))
    x[ok] <- toupper(substring(s, from[ok], to[ok]))
    m <- ok & !plus
    if (any(m)) x[m] <- as.character(reverseComplement(DNAStringSet(x[m])))
    dt[i, ctx := x]
  }
  dt[, `:=`(m1 = substr(ctx, 1, 1), p1 = substr(ctx, 2, 2),
            arun = pmin(nchar(sub("^(A*).*$", "\\1", substr(ctx, 2, CTX_LEN))), ARUN_CAP))]
  dt[, ya := !is.na(ctx) & p1 == "A" & m1 %in% c("C", "T")]
  dt
}
for (k in names(cage)) {
  cage[[k]] <- add_context(cage[[k]][n_total >= MIN_DEPTH])
  cage[[k]][, dbin := floor(log2(n_total))]
  cage[[k]][, stratum := paste(m1, arun, dbin, sep = "|")]
}

# ── overlap of NAD positions with CAGE (all libraries) ───────────────────────
all_raw <- pool(fread(cmd = sprintf("gzip -dc '%s'", ctss_f))[chr != "chrM"])
overlap <- rbindlist(lapply(strains, function(st) {
  p <- nad_pos[strain == st]
  near <- function(w) {
    q <- p[, .(chr, strand, pos, start = pos - w, end = pos + w)]
    s <- all_raw[, .(chr, strand, start = pos, end = pos)]
    setkey(s, chr, strand, start, end)
    nrow(unique(foverlaps(q, s, nomatch = NULL), by = c("chr", "strand", "pos")))
  }
  ex <- merge(p, all_raw[, .(chr, pos, strand, n_total)], by = c("chr", "pos", "strand"))
  ex <- add_context(ex)
  data.table(strain = st, nad_positions = nrow(p), cage_exact = nrow(ex),
             cage_pm1 = near(1L), cage_pm5 = near(5L),
             exact_ya = sum(ex$ya), exact_ya_depth = sum(ex$ya & ex$n_total >= MIN_DEPTH))
}))
rm(all_raw)
fwrite(overlap, file.path(outdir, "overlap.tsv"), sep = "\t")
print(overlap)

# ── permutation test ─────────────────────────────────────────────────────────
tested_bg <- function(st) unique(nad[strain == st & biotype == "none" & !enriched, .(chr, pos, strand)])
runs <- list(
  A = list(strain = MUT,  cage = "all", bg = NULL),
  B = list(strain = "WT", cage = "all", bg = NULL),
  C = list(strain = MUT,  cage = "glu", bg = NULL),
  D = list(strain = MUT,  cage = "all", bg = "tested"))
runs <- runs[vapply(runs, function(r) r$strain %in% strains, logical(1))]

set.seed(SEED)
comp <- list(); strata_out <- list(); draws_out <- list(); pos_out <- list()
for (rn in names(runs)) {
  r <- runs[[rn]]; cg <- cage[[r$cage]]
  is_nad <- cg[nad_pos[strain == r$strain], on = .(chr, pos, strand), which = TRUE, nomatch = NULL]
  obs <- cg[is_nad][ya == TRUE]
  ctl <- cg[ya == TRUE][!any_enriched, on = .(chr, pos, strand)]
  if (identical(r$bg, "tested")) ctl <- ctl[tested_bg(r$strain), on = .(chr, pos, strand), nomatch = NULL]

  st_n <- obs[, .(nad = .N), by = stratum]
  st_c <- ctl[, .(controls = .N), by = stratum]
  st <- merge(st_n, st_c, by = "stratum", all.x = TRUE)[is.na(controls), controls := 0L]
  strata_out[[rn]] <- st[, run := rn]
  dropped <- st[controls == 0, sum(nad)]
  obs <- obs[stratum %in% st[controls > 0, stratum]]
  msg("run %s (%s, CAGE %s%s): %d NAD positions tested, %d dropped for lack of controls (%.1f %%), %d controls",
      rn, r$strain, r$cage, if (is.null(r$bg)) "" else ", NAD-tested background",
      nrow(obs), dropped, 100 * dropped / max(1, nrow(obs) + dropped), nrow(ctl))
  pos_out[[rn]] <- obs[, c("chr", "pos", "strand", "ctx", "m1", "arun", "dbin", "n_total", NCOLS), with = FALSE][, run := rn]

  # Stratified label permutation: within each stratum the NAD and control
  # positions are pooled and the NAD label is reassigned to as many positions
  # as the stratum has NAD positions, without replacement. Unlike drawing
  # controls with replacement, this stays exact when the control pool is
  # small (run D), because the pool's own sampling error is part of the null.
  both <- rbind(obs[, c("stratum", "n_total", NCOLS), with = FALSE],
                ctl[stratum %in% obs$stratum, c("stratum", "n_total", NCOLS), with = FALSE])
  cm <- as.matrix(both[, ..NCOLS]); ctot <- both$n_total
  pools <- split(seq_len(nrow(both)), both$stratum)[st[controls > 0, stratum]]
  need <- st[controls > 0, nad]
  draw_share <- matrix(NA_real_, N_DRAWS, length(NCOLS), dimnames = list(NULL, END_TYPES))
  draw_mean  <- draw_share
  for (d in seq_len(N_DRAWS)) {
    idx <- unlist(mapply(function(p, k) p[sample.int(length(p), k)],
                         pools, need, SIMPLIFY = FALSE), use.names = FALSE)
    draw_share[d, ] <- colSums(cm[idx, , drop = FALSE]) / sum(ctot[idx])
    draw_mean[d, ]  <- colMeans(cm[idx, , drop = FALSE] / ctot[idx])
  }
  om <- as.matrix(obs[, ..NCOLS])
  o_share <- colSums(om) / sum(obs$n_total)
  o_mean  <- colMeans(om / obs$n_total)
  med <- apply(draw_share, 2, median)
  pval <- vapply(seq_along(END_TYPES), function(j)
    (1 + sum(abs(draw_share[, j] - med[j]) >= abs(o_share[j] - med[j]))) / (N_DRAWS + 1), numeric(1))
  # The effect is measured against the controls alone, reweighted to the NAD
  # stratum composition: the permutation null mixes NAD and control positions
  # and would dilute the ratio when NAD positions are a large share of a pool.
  cs <- ctl[stratum %in% obs$stratum, c(list(n_ctl = .N, tot = mean(n_total)),
                                        lapply(.SD, mean)), by = stratum, .SDcols = NCOLS]
  cs <- merge(cs, st[controls > 0, .(stratum, k = nad)], by = "stratum")
  c_share <- colSums(as.matrix(cs[, ..NCOLS]) * cs$k) / sum(cs$tot * cs$k)
  comp[[rn]] <- data.table(run = rn, strain = r$strain, cage = r$cage,
                           background = if (is.null(r$bg)) "all" else r$bg,
                           type = END_TYPES, n_nad_pos = nrow(obs), n_nad_reads = sum(obs$n_total),
                           nad_share = o_share, ctrl_share = c_share,
                           ratio = o_share / c_share,
                           null_median = med,
                           null_lo = apply(draw_share, 2, quantile, 0.025),
                           null_hi = apply(draw_share, 2, quantile, 0.975),
                           p = pval,
                           nad_mean_pos_share = o_mean,
                           null_mean_pos_share = apply(draw_mean, 2, median))
  draws_out[[rn]] <- melt(as.data.table(draw_share)[, `:=`(run = rn, draw = .I)],
                          id.vars = c("run", "draw"), variable.name = "type", value.name = "share")
}
comp <- rbindlist(comp)
fwrite(comp, out_comp, sep = "\t")
fwrite(rbindlist(strata_out), file.path(outdir, "strata.tsv"), sep = "\t")
fwrite(rbindlist(pos_out), file.path(outdir, "nad_positions_cage.tsv.gz"), sep = "\t")

msg("reported types:")
print(comp[type %in% REPORTED, .(run, type, n_nad_pos, nad_share = signif(nad_share, 3),
                                 ctrl_share = signif(ctrl_share, 3), ratio = round(ratio, 2), p)])

# ── figure ───────────────────────────────────────────────────────────────────
dr <- rbindlist(draws_out)[type %in% REPORTED]
dr[, type := factor(type, levels = REPORTED)]
ob <- comp[type %in% REPORTED][, type := factor(type, levels = REPORTED)]
lab <- function(dt) dt[, run_lab := sprintf("%s: %s, CAGE %s, bg %s", run, strain, cage, background)]
rl <- unique(lab(copy(ob))[, .(run, run_lab)])
dr <- merge(dr, rl, by = "run"); ob <- lab(ob)
p <- ggplot(dr, aes(x = 100 * share)) +
  geom_histogram(bins = 40, fill = "grey70") +
  geom_vline(data = ob, aes(xintercept = 100 * nad_share), colour = "#e63946", linewidth = 0.8) +
  geom_text(data = ob, aes(x = 100 * nad_share, y = Inf, label = sprintf("x%.2f\np=%.3g", ratio, p)),
            vjust = 1.2, hjust = -0.05, size = 2.6, colour = "#e63946") +
  facet_grid(run_lab ~ type, scales = "free") +
  labs(x = "share of CAGE reads, % (grey: stratified label-permutation null, red: NAD positions)", y = "draws",
       title = "CAGE end types at NAD captureSeq 5' ends (Y|A starts)") +
  theme_bw(base_size = 9) + theme(strip.text.y = element_text(angle = 0, hjust = 0))
ggsave(file.path(figdir, "06_end_types_nad_vs_control.pdf"), p,
       width = 2.4 * length(REPORTED) + 2.5, height = 1.8 * length(runs) + 1)

sessionInfo()
