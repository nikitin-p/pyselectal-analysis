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
SITE_GENE <- c(`chrM+58008` = "21S", `chrM+73704` = "COX2", `chrM+6548` = "15S mature end",
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

site_lab <- function(id) {
  s <- sub("^chrM([+-])([0-9]+)$", "\\2 (\\1)", id)
  ifelse(id %in% names(SITE_GENE), paste(s, SITE_GENE[id]), s)
}
strain_name <- function(s) fifelse(s == "WT", "WT", "triple mutant")

# ── A: one row per promoter start, CAGE and NAD side by side ─────────────────
# Rows are the sites with any CAGE or NAD read at +1, ordered by the WT NAD
# signal. Every value is read at +1 (offset 0).
p1 <- ms[offset == 0]
keep <- p1[, .(any = any(cage_n_total > 0 | meanS > 0 | meanN > 0)), by = site_id][any == TRUE, site_id]
p1 <- p1[site_id %in% keep]
ord_a <- p1[strain == "WT"][order(meanS)]$site_id
p1[, row := factor(site_lab(site_id), site_lab(ord_a))]
cg <- unique(p1[, .(row, cage_n_total, frac = 100 * get(fA1))])

th_row <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 9, face = "bold"))
no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

a1 <- ggplot(cg, aes(cage_n_total, row)) +
  geom_col(fill = "grey60", width = 0.6, orientation = "y") +
  geom_text(aes(label = format(cage_n_total, big.mark = ",")), hjust = -0.15, size = 2.5) +
  scale_x_log10(limits = c(1, 3e7), breaks = c(10, 1e3, 1e5), labels = c("10", "1,000", "100,000")) +
  labs(x = "reads (log)", y = NULL, title = "CAGE: reads at +1") + th_row
a2 <- ggplot(cg, aes(frac, row)) +
  geom_col(fill = COL_HI, alpha = 0.7, width = 0.6, orientation = "y") +
  geom_text(aes(label = sprintf("%.0f %%", frac)), hjust = -0.15, size = 2.5) +
  scale_x_continuous(limits = c(0, 135), breaks = c(0, 50, 100)) +
  labs(x = "% of reads", y = NULL, title = sprintf("CAGE: %s share", TYPE_A1)) + th_row + no_y
nad_panel <- function(st) {
  d <- p1[strain == st]
  dl <- melt(d[, .(row, called, S = meanS, N = meanN)], id.vars = c("row", "called"),
             variable.name = "lib", value.name = "reads")
  ggplot(d, aes(y = row)) +
    geom_segment(aes(x = meanN + 1, xend = meanS + 1, yend = row), colour = "grey75", linewidth = 0.6) +
    geom_point(data = dl, aes(x = reads + 1, shape = lib,
                              colour = fifelse(lib == "S" & called, "called", as.character(lib))),
               size = 2.2, stroke = 0.7) +
    scale_shape_manual(NULL, values = c(S = 16, N = 1)) +
    scale_colour_manual(NULL, values = c(called = COL_HI, S = "black", N = "grey40"),
                        breaks = "called", labels = "S, enriched (padj < 0.05)") +
    scale_x_log10(limits = c(1, 3000), breaks = c(1, 11, 101, 1001), labels = c("0", "10", "100", "1,000")) +
    labs(x = "normalised reads (log)", y = NULL, title = sprintf("NAD captureSeq, %s", strain_name(st))) +
    th_row + no_y
}
a_nad <- lapply(strain_lv, nad_panel)
pA <- wrap_plots(c(list(a1, a2), a_nad), nrow = 1, widths = c(1.1, 0.8, 1.1, 1.1), guides = "collect") +
  plot_annotation(title = "A. chrM promoter starts (+1 = final A of ATATAAGTA): CAGE against NAD captureSeq",
                  theme = theme(plot.title = element_text(size = 10, face = "bold"))) &
  theme(legend.position = "bottom")

# ── B: where the S reads start around +1 ─────────────────────────────────────
ref <- rbind(nuc[, .(what = "nuclear NAD positions", strain, log2FoldChange)],
             ct[class == "tRNA" & is.finite(log2FoldChange), .(what = "chrM tRNA (not NAD)", strain, log2FoldChange)])
ref <- ref[, .(med = median(log2FoldChange)), by = .(what, strain)]
ref_txt <- paste0("For reference, median log2(S/N): ",
                  paste(ref[, sprintf("%s, %s %.1f", what, strain_name(strain), med)], collapse = "; "),
                  ".\nPromoter starts not drawn have no read in CAGE, S or N.")

pc <- ms[site_id %in% keep]
pc[, row := factor(site_lab(site_id), site_lab(ord_a))]
pc[, strain := factor(strain_name(strain), strain_name(strain_lv))]
mk <- merge(unique(pc[, .(strain, site_id, row)]), MATURE_21S, by = "site_id")
pB <- ggplot(pc, aes(factor(offset), row, fill = log10(meanS + 1))) +
  geom_tile(colour = "grey92") +
  geom_tile(data = mk, aes(factor(offset), row), inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_gradient("S reads\n(normalised)", low = "white", high = COL_HI,
                      breaks = 0:3, labels = c("0", "10", "100", "1,000")) +
  facet_wrap(~ strain) +
  labs(x = "position relative to the promoter A (0 = +1)", y = NULL,
       title = "B. Where the S reads start around the promoter A (boxed: mature 21S rRNA 5' end, a T, cannot carry NAD)",
       caption = ref_txt) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 9, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7))

p <- wrap_elements(pA) / pB + plot_layout(heights = c(1, 0.9))
ggsave(fig_pdf, p, width = 10, height = 8)
message("wrote ", out_tsv, " and ", fig_pdf)

sessionInfo()
