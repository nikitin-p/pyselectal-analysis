#!/usr/bin/env Rscript
# NAD branch step 7 (Mac part): join the chrM ATATAAGTA sites to CAGE and plot.
#
# Reads the cluster output of deseq_chrM_sites.R, copied to the Mac:
#   scp -P 22 'pnikitin@85.89.112.24:~/pyselectal-analysis/nad/results/07_chrM_sites/*.tsv' \
#       nad/results/07_chrM_sites/
# and joins every motif-window position to the CAGE counts of
# results/tss/typed_ctss.tsv.gz, pooled over the libraries. The nuclear NAD
# reference is read from the step-5 tables in nad/results/05_deseq/.
# Spec: docs/specs/2026-11-04-nad-step7-chrM-promoter-sites.md
#
# Outputs:
#   nad/results/07_chrM_sites/motif_sites_cage.tsv
#   nad/figures/07_chrM_sites.pdf
#     A  log2FC S/N at +1: nuclear NAD positions, motif sites, processing
#        sites, chrM tRNA. +1 rather than the best window offset, because at
#        the 21S promoter the best offset is +2 (58010), the processed 5' end
#        of the mature rRNA (Q0158 starts at 58009), which cannot carry NAD;
#     B  CAGE single-A clip share against NAD log2FC per motif site at +1;
#        sites with CAGE reads but no NAD reads are drawn at x = 0;
#     C  mean normalised S reads per site across the -3..+3 window, with the
#        mature 21S end marked.
#
# Usage (from project root):
#   Rscript nad/scripts/07_chrM_sites/07_join_cage.R [--ctss results/tss/typed_ctss.tsv.gz]
#       [--indir nad/results/07_chrM_sites] [--step5-dir nad/results/05_deseq]
#       [--figdir nad/figures] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
ctss_f  <- get_arg("--ctss",      "results/tss/typed_ctss.tsv.gz")
indir   <- get_arg("--indir",     "nad/results/07_chrM_sites")
s5dir   <- get_arg("--step5-dir", "nad/results/05_deseq")
figdir  <- get_arg("--figdir",    "nad/figures")
force   <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "..", "..", "scripts", "utils", "end_types.R"))

out_tsv <- file.path(indir, "motif_sites_cage.tsv")
fig_pdf <- file.path(figdir, "07_chrM_sites.pdf")
if (file.exists(fig_pdf) && !force) {
  message("output exists, skipping (use --force): ", fig_pdf)
  quit(save = "no", status = 0)
}
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# The single-A clip, read from the vocabulary as in step 1.
TYPE_A1 <- END_TYPES_REAL[END_TYPE_CLIP_BASE[END_TYPES_REAL] %in% "A" &
                          startsWith(END_TYPES_REAL, "1S")]
stopifnot(length(TYPE_A1) == 1)
col_A1 <- paste0("n_", TYPE_A1)
COL_HI <- "#e63946"
# Gene context of the sites, for figure labels only (step-1 top-site table).
SITE_GENE <- c(`chrM+58008` = "21S", `chrM+73704` = "COX2", `chrM+6548` = "15S",
               `chrM+6470` = "15S prom.", `chrM+46079` = "ATP9", `chrM+13357` = "COX1")
MATURE_21S <- data.table(site_id = "chrM+58008", offset = 2L, note = "mature 21S 5' end")

ms <- fread(file.path(indir, "motif_sites.tsv"))
strain_lv <- unique(ms$strain)
ct <- fread(file.path(indir, "controls.tsv"))

# ── CAGE, chrM only, pooled over libraries ───────────────────────────────────
cage <- fread(cmd = sprintf("gzip -dc '%s' | awk -F'\\t' 'NR==1 || $1==\"chrM\"'", ctss_f),
              select = c("chr", "pos", "strand", col_A1, "n_total", "sample_id"))
message("CAGE libraries: ", uniqueN(cage$sample_id))
cage <- cage[, .(cage_n_total = sum(n_total), cage_n_A1 = sum(get(col_A1))), by = .(chr, pos, strand)]
ms <- merge(ms, cage, by = c("chr", "pos", "strand"), all.x = TRUE, sort = FALSE)
ms[is.na(cage_n_total), `:=`(cage_n_total = 0L, cage_n_A1 = 0L)]
ms[, cage_frac_A1 := fifelse(cage_n_total > 0, cage_n_A1 / cage_n_total, NA_real_)]
setnames(ms, c("cage_n_A1", "cage_frac_A1"), paste0("cage_", c("n_", "frac_"), TYPE_A1))
setorder(ms, strain, strand, site, offset)
fwrite(ms, out_tsv, sep = "\t")
fA1 <- paste0("cage_frac_", TYPE_A1)

# ── summaries printed to the log ─────────────────────────────────────────────
best <- ms[order(-meanS, abs(offset)), .SD[1], by = .(strain, site_id)]
best[meanS == 0, offset := 0L]                       # no S reads anywhere: report +1
print(best[order(strain, -meanS), .(strain, class, site_id, best_offset = offset, meanS = round(meanS, 1),
                                    meanN = round(meanN, 1), lfc = round(log2FoldChange, 2),
                                    padj = signif(padj, 2), called)])
print(merge(ms[offset == 0, .(sites = .N, with_S_reads = sum(meanS > 0), called_at_p1 = sum(called)), by = strain],
            best[, .(called_best_offset = sum(called)), by = strain], by = "strain"))

# ── reference: nuclear NAD positions of step 5 ───────────────────────────────
nuc <- rbindlist(lapply(unique(ms$strain), function(st) {
  d <- fread(file.path(s5dir, sprintf("%s_S_vs_N.tsv.gz", st)))
  d[chr != "chrM" & biotype == "none" & base_p1 == "A" & !is.na(padj) & padj < 0.05 &
      log2FoldChange > 1, .(strain = st, log2FoldChange)]
}))

site_lab <- function(id) ifelse(id %in% names(SITE_GENE),
                                sprintf("%s %s", sub("chrM", "", id), SITE_GENE[id]), sub("chrM", "", id))
grp_lv <- c("nuclear NAD", "ATATAAGTA (+1)", "AAGTA core, step 1 (+1)",
            "processing", "chrM tRNA")
p1 <- ms[offset == 0]
pa <- rbind(nuc[, .(strain, grp = grp_lv[1], log2FoldChange, called = TRUE, lab = NA_character_)],
            p1[, .(strain, grp = fifelse(class == "motif", grp_lv[2], grp_lv[3]), log2FoldChange,
                   called, lab = site_lab(site_id))],
            ct[class == "processing", .(strain, grp = grp_lv[4], log2FoldChange, called,
                                        lab = sub("chrM", "", site_id))],
            ct[class == "tRNA", .(strain, grp = grp_lv[5], log2FoldChange, called, lab = NA_character_)])
pa <- pa[is.finite(log2FoldChange)]
pa[, grp := factor(grp, grp_lv)]
pa[, strain := factor(strain, strain_lv)]

pA <- ggplot(pa, aes(grp, log2FoldChange)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_violin(data = pa[grp %in% grp_lv[c(1, 5)]], fill = "grey85", colour = NA, scale = "width") +
  geom_point(data = pa[!grp %in% grp_lv[c(1, 5)]], aes(colour = called), size = 1.6) +
  geom_text_repel(data = pa[grp %in% grp_lv[2:3]], aes(label = lab), size = 2,
                  direction = "y", nudge_x = 0.35, hjust = 0, segment.size = 0.2,
                  min.segment.length = 0, seed = 1) +
  scale_colour_manual("padj < 0.05 and\nlog2FC > 1", values = c(`TRUE` = COL_HI, `FALSE` = "grey40")) +
  scale_x_discrete(limits = grp_lv) +
  facet_wrap(~ strain) +
  labs(x = NULL, y = "log2 fold change S / N",
       title = "A. NAD captureSeq S vs N at +1 of the chrM promoter starts",
       caption = "Sites without any S or N reads have no log2FC and are not drawn.") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

ms[, strain := factor(strain, strain_lv)]
pb <- ms[offset == 0 & cage_n_total > 0]
pb[, `:=`(x = fifelse(is.finite(log2FoldChange), log2FoldChange, 0),
          nad = fifelse(is.finite(log2FoldChange), "NAD reads", "no NAD reads"),
          lab = site_lab(site_id))]
pB <- ggplot(pb, aes(x, 100 * get(fA1))) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = cage_n_total, colour = called, shape = nad), stroke = 0.6) +
  geom_text_repel(aes(label = lab), size = 2, segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.4, max.overlaps = Inf, seed = 1) +
  scale_size("CAGE reads", trans = "log10", range = c(1, 5), breaks = c(10, 100, 1e3, 1e4, 1e5),
             labels = c("10", "100", "1,000", "10,000", "100,000")) +
  scale_shape_manual(NULL, values = c(`NAD reads` = 16, `no NAD reads` = 1)) +
  scale_colour_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey40"), guide = "none") +
  facet_wrap(~ strain) +
  labs(x = "NAD log2 fold change S / N at +1", y = sprintf("CAGE %s share at +1 (%%)", TYPE_A1),
       title = sprintf("B. CAGE %s share against NAD enrichment", TYPE_A1)) +
  theme_bw(base_size = 9)

ord <- unique(ms[order(strand, site)]$site_id)
ms[, site_lab := factor(site_lab(site_id), site_lab(ord))]
mk <- merge(unique(ms[, .(strain, site_id, site_lab)]), MATURE_21S, by = "site_id")
pC <- ggplot(ms, aes(factor(offset), site_lab, fill = log10(meanS + 1))) +
  geom_tile() +
  geom_point(data = ms[called == TRUE], shape = 4, size = 0.8) +
  geom_tile(data = mk, aes(factor(offset), site_lab), inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.4) +
  scale_fill_gradient("log10 mean\nnormalised S + 1", low = "white", high = COL_HI) +
  facet_wrap(~ strain) +
  labs(x = "offset from the final A of ATATAAGTA (nt)", y = NULL,
       title = "C. S reads across the window (x: called; boxed: mature 21S end, not an NAD start)") +
  theme_bw(base_size = 9)

p <- pA / pB / pC + plot_layout(heights = c(1, 1, 1.2))
ggsave(fig_pdf, p, width = 8, height = 11)
message("wrote ", out_tsv, " and ", fig_pdf)

sessionInfo()
