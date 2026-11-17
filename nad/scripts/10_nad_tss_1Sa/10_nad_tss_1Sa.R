#!/usr/bin/env Rscript
# NAD branch step 10: how many NAD TSSs produce at least one 1Sa read in CAGE,
# against how many are expected to.
#
# NAD TSSs are the step-5 positions enriched in S over N (padj < 0.05,
# log2FC > 1) with an A at +1; nuclear positions must be unmasked
# (biotype "none"), while chrM positions keep every biotype, because the chrM
# "rRNA" mask covers the real primary transcripts of 21S and 15S. Each set is
# also cut to its top quartile by mean normalised S (meanS >= Q3).
#
# CAGE is the pool of the seven yeast libraries in typed_ctss.tsv.gz. A clipped
# A is impossible when the genomic base at -1 is A, so the observed and expected
# counts are taken over Y|A starts (-1 = C or T) present in CAGE only.
#
# Expectation: for every NAD TSS, the share of background positions with >= 1
# 1Sa read in the same stratum, summed over the set. Background = nuclear Y|A
# CAGE positions not enriched in either strain; stratum = -1 base x A-run from
# +1 (capped at 5) x floor(log2 CAGE depth) (capped at 12); a stratum without
# background positions falls back to the same -1 base and A-run at depth bins
# >= 10 (column depth_fallback). Variance = sum f(1-f);
# p = Poisson upper tail at the observed count. For chrM sets a second
# expectation uses chrM Y|A background positions, stratified by depth only
# (capped at bin 8), because chrM has too few positions for the full stratum.
#
# Outputs (in --outdir):
#   summary.tsv          one row per strain x compartment x set
#   positions.tsv.gz     every NAD TSS with CAGE depth, 1Sa reads, context and
#                        its expected probability of >= 1 1Sa read
#   <figdir>/10_nad_tss_1Sa.{pdf,png}
#
# Usage (from project root):
#   Rscript nad/scripts/10_nad_tss_1Sa/10_nad_tss_1Sa.R [--force]
#     [--ctss results/tss/typed_ctss.tsv.gz] [--nad-dir nad/results/05_deseq]
#     [--fasta data/sacCer3.fa] [--outdir nad/results/10_nad_tss_1Sa]
#     [--figdir nad/figures]

suppressPackageStartupMessages({
  library(data.table); library(Biostrings); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
ctss_f  <- get_arg("--ctss",    "results/tss/typed_ctss.tsv.gz")
nad_dir <- get_arg("--nad-dir", "nad/results/05_deseq")
fasta_f <- get_arg("--fasta",   "data/sacCer3.fa")
outdir  <- get_arg("--outdir",  "nad/results/10_nad_tss_1Sa")
figdir  <- get_arg("--figdir",  "nad/figures")
force   <- "--force" %in% args
PADJ <- 0.05; LFC <- 1; ARUN_CAP <- 5; DBIN_CAP <- 12; DBIN_CAP_M <- 8; HI_BIN <- 10
STRAINS <- c(WT = "WT", triple = "dxo1_rai1_npy1")

out_sum <- file.path(outdir, "summary.tsv")
if (file.exists(out_sum) && !force) { message("[skip] ", out_sum, " exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

# ── CAGE, pooled, with genomic context ───────────────────────────────────────
ct <- fread(ctss_f, select = c("chr", "pos", "strand", "n_1Sa", "n_total"))
cage <- ct[, .(n = sum(n_total), sa = sum(n_1Sa)), by = .(chr, pos, strand)]
rm(ct)
g <- readDNAStringSet(fasta_f); names(g) <- sub(" .*", "", names(g))
ctx_of <- function(chr, pos, strand) {   # -1 and +1..+5 on the read strand
  L <- width(g)[match(chr, names(g))]
  s <- ifelse(strand == "+", pos - 1, pos - 5); e <- ifelse(strand == "+", pos + 5, pos + 1)
  ok <- !is.na(L) & s >= 1 & e <= L; out <- rep(NA_character_, length(chr))
  sq <- subseq(g[chr[ok]], s[ok], e[ok]); mi <- strand[ok] == "-"
  sq[mi] <- reverseComplement(sq[mi]); out[ok] <- toupper(as.character(sq)); out
}
add_ctx <- function(d) {
  d[, ctx := ctx_of(chr, pos, strand)]
  d[, `:=`(m1 = substr(ctx, 1, 1), p1 = substr(ctx, 2, 2))]
  d[, arun := pmin(nchar(sub("^(A*).*$", "\\1", substr(ctx, 2, 6))), ARUN_CAP)]
  d[, dbin := pmin(floor(log2(pmax(n, 1))), DBIN_CAP)]
  d[, has := sa > 0]
  d
}
cage <- add_ctx(cage)[!is.na(ctx)]
msg("CAGE positions with context: %d", nrow(cage))

# ── NAD sets and background ──────────────────────────────────────────────────
nad <- rbindlist(lapply(STRAINS, function(s)
  fread(file.path(nad_dir, sprintf("%s_S_vs_N.tsv.gz", s)))), idcol = "strain")
nad[, enr := !is.na(padj) & padj < PADJ & log2FoldChange > LFC]
excl <- unique(nad[enr == TRUE, paste(chr, pos, strand)])
isYA <- function(d) d$p1 == "A" & d$m1 %in% c("C", "T")
bg  <- cage[chr != "chrM" & isYA(cage) & !(paste(chr, pos, strand) %in% excl)]
bgM <- cage[chr == "chrM" & isYA(cage) & !(paste(chr, pos, strand) %in% excl)]
frac  <- bg[, .(f = mean(has)), by = .(m1, arun, dbin)]
fracM <- bgM[, .(fM = mean(has)), by = .(dbinM = pmin(dbin, DBIN_CAP_M))]
# A very deep position can fall in a stratum the background does not fill (for
# example chrM:58011, T|AAAAA at 13k reads); it then takes the background of
# its -1 base and A-run pooled over every depth bin >= HI_BIN.
frac_hi <- bg[dbin >= HI_BIN, .(f_hi = mean(has)), by = .(m1, arun)]
msg("background Y|A positions: nuclear %d, chrM %d", nrow(bg), nrow(bgM))

pos_all <- list(); res <- list()
for (st in names(STRAINS)) for (comp in c("nuclear", "chrM")) {
  d <- nad[strain == st & enr == TRUE & base_p1 == "A" &
           (if (comp == "nuclear") chr != "chrM" & biotype == "none" else chr == "chrM"),
           .(chr, pos, strand, biotype, meanS)]
  d[, topQ := meanS >= quantile(meanS, 0.75)]
  d <- merge(d, cage[, .(chr, pos, strand, n, sa)], by = c("chr", "pos", "strand"), all.x = TRUE)
  d[is.na(n), `:=`(n = 0L, sa = 0L)]
  d <- add_ctx(d)
  d[, dbinM := pmin(dbin, DBIN_CAP_M)]
  d <- merge(d, frac, by = c("m1", "arun", "dbin"), all.x = TRUE)
  d <- merge(d, frac_hi, by = c("m1", "arun"), all.x = TRUE)
  d[, fallback := is.na(f) & dbin >= HI_BIN & n > 0 & m1 %in% c("C", "T")]
  d[fallback == TRUE, f := f_hi]
  d[, f_hi := NULL]
  d <- merge(d, fracM, by = "dbinM", all.x = TRUE)
  d[!(n > 0 & m1 %in% c("C", "T")), `:=`(f = NA, fM = NA)]
  d[, `:=`(strain = st, comp = comp)]
  pos_all[[length(pos_all) + 1]] <- d
  for (sub in c("all", "topQ")) {
    x <- if (sub == "all") d else d[topQ == TRUE]
    ya <- x[n > 0 & m1 %in% c("C", "T")]
    e <- sum(ya$f, na.rm = TRUE)
    res[[length(res) + 1]] <- data.table(
      strain = st, comp = comp, set = sub, n_nad = nrow(x),
      in_cage = sum(x$n > 0), in_cage_ge10 = sum(x$n >= 10), median_cage = median(x$n),
      m1_A_in_cage = sum(x$n > 0 & x$m1 == "A"), ya_in_cage = nrow(ya),
      obs_ge1_1Sa = sum(ya$has), exp_nuclear = e,
      sd_nuclear = sqrt(sum(ya$f * (1 - ya$f), na.rm = TRUE)),
      p_nuclear = 1 - ppois(sum(ya$has) - 1, e),
      exp_chrM = if (comp == "chrM") sum(ya$fM, na.rm = TRUE) else NA_real_,
      reads_1Sa = sum(ya$sa), reads = sum(ya$n),
      depth_fallback = sum(ya$fallback), no_stratum = sum(is.na(ya$f)))
  }
}
r <- rbindlist(res)
stopifnot(all(r$no_stratum == 0))
pos_all <- rbindlist(pos_all, use.names = TRUE)
fwrite(r, out_sum, sep = "\t")
fwrite(pos_all[order(strain, comp, -meanS),
               .(strain, comp, chr, pos, strand, biotype, meanS, topQ, ctx, m1, arun,
                 cage_reads = n, cage_1Sa = sa, p_ge1_nuclear = f, depth_fallback = fallback,
                 p_ge1_chrM = fM)],
       file.path(outdir, "positions.tsv.gz"), sep = "\t")
print(r)

# ── figure ───────────────────────────────────────────────────────────────────
WHAT <- c("наблюдается", "ожидание: ядерный фон", "ожидание: фон chrM")
r[, lab := sprintf("%s · %s", ifelse(strain == "WT", "WT", "dxo1Δ rai1Δ npy1Δ"),
                   ifelse(comp == "nuclear", "ядро", "chrM"))]
r[, setl := ifelse(set == "all", "все NAD TSS", "верхняя квартиль")]
r[, panel := sprintf("%s · %s\nNAD TSS %d, в CAGE %d (%.0f%%), Y|A %d\nмедиана глубины CAGE %s; p = %s",
                     lab, setl, n_nad, in_cage, 100 * in_cage / n_nad, ya_in_cage,
                     format(round(median_cage), big.mark = " ", trim = TRUE), signif(p_nuclear, 2))]
r[, panel := factor(panel, levels = panel)]
long <- rbind(
  r[, .(panel, what = WHAT[1], v = obs_ge1_1Sa, lo = NA_real_, hi = NA_real_)],
  r[, .(panel, what = WHAT[2], v = exp_nuclear, lo = exp_nuclear - 1.96 * sd_nuclear,
        hi = exp_nuclear + 1.96 * sd_nuclear)],
  r[comp == "chrM", .(panel, what = WHAT[3], v = exp_chrM, lo = NA_real_, hi = NA_real_)])
long[, what := factor(what, levels = WHAT)]
long[, xl := pmax(v, hi, na.rm = TRUE)]
p <- ggplot(long, aes(x = v, y = what, fill = what)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, na.rm = TRUE) +
  geom_text(aes(x = xl, label = round(v, 1)), hjust = -0.25, size = 3.3) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_fill_manual(values = setNames(c("#e63946", "grey55", "#457b9d"), WHAT)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "число Y|A NAD TSS, у которых в CAGE есть ≥1 рид 1Sa", y = NULL,
       title = "NAD TSS с хотя бы одной 1Sa в CAGE: наблюдаемое против ожидаемого",
       subtitle = paste0("Ожидание = сумма долей фоновых позиций с ≥1 1Sa в той же страте (основание −1, A-run, глубина CAGE).\n",
                         "Усы = ±1.96 SD. p — превышение над ядерным фоном (Пуассон).")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(size = 8.5, hjust = 0),
        panel.grid.major.y = element_blank())
ggsave(file.path(figdir, "10_nad_tss_1Sa.pdf"), p, width = 9, height = 10, device = cairo_pdf)
ggsave(file.path(figdir, "10_nad_tss_1Sa.png"), p, width = 9, height = 10, dpi = 150)
msg("written: %s, %s", out_sum, file.path(figdir, "10_nad_tss_1Sa.{pdf,png}"))
sessionInfo()
