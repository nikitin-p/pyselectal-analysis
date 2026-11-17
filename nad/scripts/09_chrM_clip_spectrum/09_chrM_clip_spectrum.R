#!/usr/bin/env Rscript
# NAD branch step 9: is the extra 5' base on chrM CAGE reads NAD-specific?
#
# On chrM, CAGE reads carry a single extra soft-clipped base not only at the
# NAD-capped promoter starts but also at processing sites that NAD captureSeq
# never sees (16471, 13987, 25343). This step asks whether the extra base at
# NAD 5' ends differs from the one at non-NAD chrM 5' ends. A clipped base can
# never equal the genomic base at -1, so every comparison is made within the
# same -1 base.
#
# Classes of chrM 5' ends (CAGE depth >= --min-depth, both strands), first match wins:
#   promoter, NAD          step-7 promoter start (offset 0) called in WT or the triple mutant;
#   promoter, no NAD       step-7 promoter start not called (6470, 46079, and 6548);
#   group II intron 5′ end first base of a chrM intron that starts GTGYG (the
#                          group II consensus); these introns are excised as
#                          lariats, and reverse transcription across the 2′–5′
#                          branch can add the branch base as an extra 5′ base;
#   A, S-captured internal +1 = A, not a promoter, NAD-detected in WT or the triple
#                          mutant (all three S runs > 0, mean normalised S >= --min-s,
#                          log2((S + 0.5) / (N + 0.5)) >= --min-lfc); NAD can sit only on
#                          a primary 5′ end, so these are read starts inside NAD-RNAs,
#                          such as premature reverse-transcription stops in the S library;
#   A, no S                +1 = A and at most one raw S read over all six S runs;
#   A, weak S              every other A start (listed, not tested);
#   other non-A            +1 = C, G or T, which cannot carry NAD.
# Per position: share of reads with each single extra base (1Sa, 1Sc, 1Sg,
# 1St), with no clip (M), with a 2-3 nt clip and with a longer clip; and the
# A fraction among single-base clips. Positions are weighted equally, so that
# 58008 does not dominate. Tests: stratified label permutation (labels shuffled
# within the -1 base), statistic = difference of the mean per-position value.
#
# Inputs: nad/results/09_chrM_clip_spectrum/fine_types.tsv.gz (fine_clip_types.py),
#         nad/results/04_ctss/matrix_reads.tsv.gz, the step-8 size factors.
# Usage (from project root):
#   Rscript nad/scripts/09_chrM_clip_spectrum/09_chrM_clip_spectrum.R [--force]
#     [--fine ...] [--matrix ...] [--sf nad/results/08_1Sa_positions/size_factors.tsv]
#     [--manifest nad/scripts/00_download/priority1_runs.tsv] [--ctss results/tss/typed_ctss.tsv.gz]
#     [--fasta data/sacCer3.fa] [--outdir nad/results/09_chrM_clip_spectrum] [--figdir nad/figures]
#     [--min-depth 50] [--min-s 5] [--min-lfc 3] [--draws 5000] [--seed 1]

suppressPackageStartupMessages({
  library(data.table); library(Biostrings); library(ggplot2); library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
fine_f   <- get_arg("--fine",     "nad/results/09_chrM_clip_spectrum/fine_types.tsv.gz")
matrix_f <- get_arg("--matrix",   "nad/results/04_ctss/matrix_reads.tsv.gz")
sf_f     <- get_arg("--sf",       "nad/results/08_1Sa_positions/size_factors.tsv")
man_f    <- get_arg("--manifest", "nad/scripts/00_download/priority1_runs.tsv")
ctss_f   <- get_arg("--ctss",     "results/tss/typed_ctss.tsv.gz")
prom_f   <- get_arg("--sites",    "nad/results/07_chrM_sites/motif_sites.tsv")
gtf_f    <- get_arg("--gtf",      "data/sacCer3.ensGene.gtf")
POL3     <- strsplit(get_arg("--pol3", "SCR1,NME1,RPR1"), ",")[[1]]   # nuclear 5'-ppp controls
fasta_f  <- get_arg("--fasta",    "data/sacCer3.fa")
outdir   <- get_arg("--outdir",   "nad/results/09_chrM_clip_spectrum")
figdir   <- get_arg("--figdir",   "nad/figures")
MIN_DEPTH <- as.integer(get_arg("--min-depth", "50"))
MIN_S     <- as.numeric(get_arg("--min-s", "5"))
MIN_LFC   <- as.numeric(get_arg("--min-lfc", "3"))
N_DRAWS   <- as.integer(get_arg("--draws", "5000"))
SEED      <- as.integer(get_arg("--seed", "1"))
force     <- "--force" %in% args

out_main <- file.path(outdir, "positions.tsv")
if (file.exists(out_main) && !force) { message("[skip] ", out_main, " exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))
checks <- list()

# ── CAGE fine types on chrM ──────────────────────────────────────────────────
fine_all <- fread(cmd = sprintf("gzip -dc '%s'", fine_f))
fine <- fine_all[chr == "chrM"]
fine[, cls := fcase(type == "M", "M", type %in% c("1Sa", "1Sc", "1Sg", "1St"), type,
                    type == "4+S", "clip4plus", default = "clip2_3")]
pos <- dcast(fine, chr + pos + strand ~ cls, value.var = "n", fun.aggregate = sum, fill = 0)
for (k in c("M", "1Sa", "1Sc", "1Sg", "1St", "clip2_3", "clip4plus")) if (!k %in% names(pos)) pos[, (k) := 0]
pos[, n_total := M + `1Sa` + `1Sc` + `1Sg` + `1St` + clip2_3 + clip4plus]

# check: totals and 1Sa agree with typed_ctss (same BAMs, same 5' end rule)
tc <- fread(cmd = sprintf("gzip -dc '%s'", ctss_f), select = c("chr", "pos", "strand", "n_1Sa", "n_total"))[chr == "chrM"]
tc <- tc[, .(t_total = sum(n_total), t_1Sa = sum(n_1Sa)), by = .(chr, pos, strand)]
cmp <- merge(pos[, .(chr, pos, strand, n_total, `1Sa`)], tc, by = c("chr", "pos", "strand"), all = TRUE)
checks$ctss <- c("chrM totals and 1Sa equal typed_ctss at every position",
                 isTRUE(cmp[, all(!is.na(n_total) & !is.na(t_total) & n_total == t_total & `1Sa` == t_1Sa)]))
pos <- pos[n_total >= MIN_DEPTH]

# genome context -1, +1 (RNA sense)
g <- readDNAStringSet(fasta_f); names(g) <- sub("\\s.*", "", names(g))
mt <- as.character(g[["chrM"]])
base_at <- function(p) toupper(substring(mt, p, p))
cmpl <- c(A = "T", C = "G", G = "C", T = "A")
pos[, `:=`(m1 = ifelse(strand == "+", base_at(pos - 1), cmpl[base_at(pos + 1)]),
           p1 = ifelse(strand == "+", base_at(pos), cmpl[base_at(pos)]))]

# ── NAD captureSeq signal ────────────────────────────────────────────────────
man <- fread(man_f, select = c("run", "strain", "library"))
sf  <- fread(sf_f)[, .(run, size_factor)]
man <- merge(man, sf, by = "run")
mat <- fread(cmd = sprintf("gzip -dc '%s'", matrix_f))[chr == "chrM"]
nad <- mat[, .(chr, pos, strand)]
for (st in unique(man$strain)) {
  rs <- man[strain == st & library == "S"]; rn <- man[strain == st & library == "N"]
  S <- sweep(as.matrix(mat[, rs$run, with = FALSE]), 2, rs$size_factor, "/")
  N <- sweep(as.matrix(mat[, rn$run, with = FALSE]), 2, rn$size_factor, "/")
  tag <- if (st == "WT") "wt" else "tm"
  nad[, paste0("s_", tag) := rowMeans(S)]; nad[, paste0("n_", tag) := rowMeans(N)]
  nad[, paste0("det_", tag) := rowSums(S > 0) == ncol(S) & rowMeans(S) >= MIN_S &
                                 log2((rowMeans(S) + 0.5) / (rowMeans(N) + 0.5)) >= MIN_LFC]
}
nad[, s_raw := rowSums(as.matrix(mat[, man[library == "S", run], with = FALSE]))]
pos <- merge(pos, nad, by = c("chr", "pos", "strand"), all.x = TRUE)
for (k in c("s_wt", "n_wt", "s_tm", "n_tm", "s_raw")) pos[is.na(get(k)), (k) := 0]
for (k in c("det_wt", "det_tm")) pos[is.na(get(k)), (k) := FALSE]
# promoter starts: step-7 sites at offset 0 (ATATAAGTA final A, plus the AAGTA-core sites)
prom <- fread(prom_f)[offset == 0, .(called = any(called %in% TRUE)), by = .(chr, pos, strand)]
# group II intron 5' ends: first intron base after a chrM exon whose intron starts GTGYG
gtf <- fread(cmd = sprintf("awk -F'\\t' '$1==\"chrM\" && $3==\"exon\"' '%s'", gtf_f), header = FALSE, sep = "\t")
gtf[, tx := sub('.*transcript_id "([^"]+)".*', "\\1", V9)]
setorder(gtf, tx, V4)
intr <- gtf[, .(start = head(V5, -1) + 1L, end = tail(V4, -1) - 1L, strand = V7[1]), by = tx][strand == "+" & end > start]
intr <- unique(intr[, .(chr = "chrM", pos = start, strand)])
intr[, head5 := toupper(substring(mt, pos, pos + 4))]
intr2 <- intr[grepl("^GTG[CT]G", head5)]
checks$introns <- c(sprintf("group II intron 5' ends found: %s", paste(intr2$pos, collapse = ",")), nrow(intr2) >= 3)
pos[, is_prom := paste(chr, pos, strand) %in% prom[, paste(chr, pos, strand)]]
pos[, prom_nad := paste(chr, pos, strand) %in% prom[called == TRUE, paste(chr, pos, strand)]]
pos[, is_g2 := paste(chr, pos, strand) %in% intr2[, paste(chr, pos, strand)]]
pos[, class := fcase(is_prom & prom_nad, "promoter, NAD",
                     is_prom, "promoter, no NAD",
                     is_g2, "group II intron 5′ end",
                     p1 == "A" & (det_wt | det_tm), "A, S-captured internal",
                     p1 == "A" & s_raw <= 1, "A, no S",
                     p1 == "A", "A, weak S",
                     default = "other non-A")]
CLASSES <- c("promoter, NAD", "promoter, no NAD", "group II intron 5′ end",
             "A, S-captured internal", "A, no S", "other non-A")

# per-position shares
for (k in c("M", "1Sa", "1Sc", "1Sg", "1St", "clip2_3", "clip4plus")) pos[, paste0("p_", k) := get(k) / n_total]
pos[, n_1S := `1Sa` + `1Sc` + `1Sg` + `1St`]
pos[, f_a := ifelse(n_1S > 0, `1Sa` / n_1S, NA_real_)]
setorder(pos, class, -n_total)
fwrite(pos, out_main, sep = "\t")
msg("chrM 5' ends with CAGE >= %d: %d (%s)", MIN_DEPTH, nrow(pos),
    paste(pos[, .N, by = class][, sprintf("%s %d", class, N)], collapse = ", "))

# ── summaries ────────────────────────────────────────────────────────────────
VARS <- c("p_M", "p_1Sa", "p_1Sc", "p_1Sg", "p_1St", "p_clip2_3", "p_clip4plus")
summ <- pos[class %in% CLASSES, c(.(n_pos = .N, reads = sum(n_total)),
                                  lapply(.SD, mean), .(f_a = mean(f_a, na.rm = TRUE))),
            by = .(m1, class), .SDcols = VARS][order(m1, class)]
fwrite(summ, file.path(outdir, "spectrum_by_class.tsv"), sep = "\t")

# ── nuclear 5'-triphosphate controls (Pol III primary transcripts) ───────────
gtf_all <- fread(cmd = sprintf("awk -F'\\t' '$3==\"transcript\"' '%s'", gtf_f), header = FALSE, sep = "\t")
gtf_all[, gene := sub('.*gene_id "([^"]+)".*', "\\1", V9)]
p3 <- gtf_all[gene %in% POL3, .(gene, chr = V1, pos = ifelse(V7 == "+", V4, V5), strand = V7)]
p3[, ctx := mapply(function(ch, p, s) { x <- toupper(substring(as.character(g[[ch]]), p - 1, p + 1))
  if (s == "+") x else as.character(reverseComplement(DNAString(x))) }, chr, pos, strand)]
p3[, `:=`(m1 = substr(ctx, 1, 1), p1 = substr(ctx, 2, 2))]
fa <- fine_all[p3, on = .(chr, pos, strand), nomatch = 0]
fa[, cls := fcase(type == "M", "M", type %in% c("1Sa", "1Sc", "1Sg", "1St"), type,
                  type == "4+S", "clip4plus", default = "clip2_3")]
pol3 <- dcast(fa, gene + chr + pos + strand + m1 + p1 ~ cls, value.var = "n", fun.aggregate = sum, fill = 0)
for (k in c("M", "1Sa", "1Sc", "1Sg", "1St", "clip2_3", "clip4plus")) if (!k %in% names(pol3)) pol3[, (k) := 0]
pol3[, n_total := M + `1Sa` + `1Sc` + `1Sg` + `1St` + clip2_3 + clip4plus]
fwrite(pol3, file.path(outdir, "pol3_controls.tsv"), sep = "\t")

# ── stratified permutation tests ─────────────────────────────────────────────
set.seed(SEED)
perm_test <- function(d, a, b, v) {
  d <- d[class %in% c(a, b) & !is.na(get(v)) & m1 != "A"]     # -1 = A forbids 1Sa
  d <- d[m1 %in% d[, .(k = uniqueN(class)), by = m1][k == 2, m1]]  # strata holding both classes
  if (nrow(d) == 0) return(NULL)
  x <- d[[v]]; lab <- d$class == a; strat <- d$m1
  stat <- function(l) mean(x[l]) - mean(x[!l])
  obs <- stat(lab)
  idx <- split(seq_along(x), strat)
  null <- replicate(N_DRAWS, { l <- lab; for (i in idx) l[i] <- sample(lab[i]); stat(l) })
  data.table(contrast = paste(a, "vs", b), variable = v, strata = paste(sort(unique(strat)), collapse = ","),
             n_a = sum(lab), n_b = sum(!lab), mean_a = mean(x[lab]), mean_b = mean(x[!lab]), diff = obs,
             p = (1 + sum(abs(null) >= abs(obs))) / (N_DRAWS + 1))
}
CONTRASTS <- list(c("promoter, NAD", "A, S-captured internal"), c("promoter, NAD", "A, no S"),
                  c("promoter, NAD", "other non-A"), c("group II intron 5′ end", "other non-A"),
                  c("A, no S", "other non-A"))
tests <- rbindlist(lapply(c("p_1Sa", "f_a", "p_M"), function(v)
  rbindlist(lapply(CONTRASTS, function(cc) perm_test(pos, cc[1], cc[2], v)))))
fwrite(tests, file.path(outdir, "tests.tsv"), sep = "\t")

# null check: labels shuffled within -1 strata before testing
nullp <- sapply(1:20, function(i) { d <- copy(pos)[class %in% c("A, no S", "other non-A")]
  d[, class := sample(class), by = m1]; perm_test(d, "A, no S", "other non-A", "p_1Sa")$p })
checks$null <- c(sprintf("shuffled labels: %d of 20 p < 0.05 (A, no S vs other non-A, p_1Sa)", sum(nullp < 0.05)),
                 sum(nullp < 0.05) <= 3)

# ── figure ───────────────────────────────────────────────────────────────────
red <- "#e63946"
cls_cols <- setNames(c(red, "#f4a3a8", "#6a4c93", "#e9c46a", "grey55", "#457b9d", "black"),
                     c(CLASSES, "nuclear 5′-ppp (SCR1)"))
th <- theme_classic(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9),
                                           plot.subtitle = element_text(size = 7.5, colour = "grey30"))
# A: composition at -1 = T, the context of every chrM promoter start
VARS2 <- c(VARS, "n_total")
scr <- pol3[gene == "SCR1" & m1 == "T", c(.(class = "nuclear 5′-ppp (SCR1)", n_pos = 1L),
          lapply(.SD, function(x) x / n_total)), .SDcols = c("M", "1Sa", "1Sc", "1Sg", "1St", "clip2_3", "clip4plus")]
setnames(scr, c("M", "1Sa", "1Sc", "1Sg", "1St", "clip2_3", "clip4plus"), VARS)
tA <- rbind(summ[m1 == "T", c("class", "n_pos", VARS), with = FALSE], scr[, c("class", "n_pos", VARS), with = FALSE])
lp <- melt(tA, id.vars = c("class", "n_pos"), measure.vars = VARS, variable.name = "what", value.name = "share")
lp[, what := factor(sub("^p_", "", what), c("1Sa", "1St", "1Sc", "1Sg", "clip2_3", "clip4plus", "M"))]
lp[, class := factor(class, names(cls_cols))]
lp[, lab := factor(sprintf("%s\n(n=%d)", class, n_pos))]
lp[, lab := factor(lab, unique(lab[order(class)]))]
base_cols <- c("1Sa" = red, "1St" = "#f4a261", "1Sc" = "#2a9d8f", "1Sg" = "#264653",
               "clip2_3" = "#b5838d", "clip4plus" = "#6d6875", "M" = "grey88")
pA <- ggplot(lp, aes(lab, share, fill = what)) + geom_col(width = 0.75) +
  scale_fill_manual(values = base_cols, name = "5′ end of\nCAGE read") + scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "mean share per position",
       title = "A  CAGE 5′ ends with T at -1 (the context of all chrM promoter starts)",
       subtitle = sprintf("chrM 5′ ends with >= %d CAGE reads; a T at -1 makes 1St impossible, so the extra base can be A, C or G", MIN_DEPTH)) +
  th + theme(axis.text.x = element_text(size = 7))
# B: per-position extra A vs extra G at -1 = T
dB <- rbind(pos[m1 == "T" & class %in% CLASSES, .(class, pos = as.character(pos), n_total, a = p_1Sa, gg = p_1Sg)],
            pol3[m1 == "T", .(class = "nuclear 5′-ppp (SCR1)", pos = gene, n_total, a = `1Sa` / n_total, gg = `1Sg` / n_total)])
dB[, class := factor(class, names(cls_cols))]
lab_pos <- dB[class %in% c("promoter, NAD", "promoter, no NAD", "group II intron 5′ end", "nuclear 5′-ppp (SCR1)") |
              (a > 0.1 | gg > 0.3)]
pB <- ggplot(dB, aes(gg, a, colour = class)) + geom_point(aes(size = n_total), alpha = 0.8) +
  ggrepel::geom_text_repel(data = lab_pos, aes(label = pos), size = 2.3, min.segment.length = 0, box.padding = 0.4, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = cls_cols, name = NULL, drop = FALSE) +
  scale_size_area(max_size = 6, trans = "log10", name = "CAGE reads", breaks = c(100, 1000, 10000, 100000)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) + scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  labs(x = "extra G (1Sg)", y = "extra A (1Sa)", title = "B  Extra A versus extra G, one dot per 5′ end (-1 = T)",
       subtitle = "Promoter starts get A almost only; the lariat end 16471 gets A and G; the nuclear ppp start SCR1 gets G") +
  th + theme(legend.position = "bottom", legend.box = "vertical") + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 2))
tl <- tests[variable %in% c("p_1Sa", "f_a")]
tl[, lab := sprintf("%s\n%s", contrast, ifelse(variable == "p_1Sa", "1Sa share of reads", "A among 1-nt clips"))]
pC <- ggplot(tl, aes(diff, lab, colour = p < 0.05)) + geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 2.6) + geom_text(aes(label = sprintf("p = %.2g", p)), nudge_y = 0.35, size = 2.4, colour = "grey30") +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = "grey40"), guide = "none") +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "difference of means (first minus second)", y = NULL,
       title = "C  Permutation tests (labels shuffled within the -1 base)") + th + theme(axis.text.y = element_text(size = 6.5))
p <- pA / (pB | pC) + plot_layout(heights = c(1, 1.15))
ggsave(file.path(figdir, "09_chrM_clip_spectrum.pdf"), p, width = 13, height = 10, device = cairo_pdf)
ggsave(file.path(figdir, "09_chrM_clip_spectrum.png"), p, width = 13, height = 10, dpi = 150)

# ── report ───────────────────────────────────────────────────────────────────
for (ck in checks) message(sprintf("CHECK %s: %s", if (isTRUE(as.logical(ck[2]))) "PASS" else "FAIL", ck[1]))
print(summ[, .(m1, class, n_pos, reads, p_M = round(p_M, 3), p_1Sa = round(p_1Sa, 3), p_1St = round(p_1St, 3),
               p_1Sc = round(p_1Sc, 3), p_1Sg = round(p_1Sg, 3), f_a = round(f_a, 3))])
print(pol3[, .(gene, m1, p1, n_total, M, `1Sa`, `1St`, `1Sc`, `1Sg`)]); print(tests[, .(contrast, variable, strata, n_a, n_b, mean_a = round(mean_a, 3), mean_b = round(mean_b, 3), p = signif(p, 2))])
sessionInfo()
