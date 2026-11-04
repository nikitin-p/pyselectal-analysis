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
#     A  log2FC S/N: nuclear NAD positions, motif sites at their best window
#        offset (the offset with the most S reads), processing sites, chrM tRNA;
#     B  CAGE single-A clip share against NAD log2FC per motif site at +1;
#     C  mean normalised S reads per site across the -3..+3 window.
#
# Usage (from project root):
#   Rscript nad/scripts/07_chrM_sites/07_join_cage.R [--ctss results/tss/typed_ctss.tsv.gz]
#       [--indir nad/results/07_chrM_sites] [--step5-dir nad/results/05_deseq]
#       [--figdir nad/figures] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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

grp_lv <- c("nuclear NAD", "ATATAAGTA (best offset)", "AAGTA core, step 1 (best offset)",
            "processing", "chrM tRNA")
pa <- rbind(nuc[, .(strain, grp = grp_lv[1], log2FoldChange, called = TRUE, lab = NA_character_)],
            best[, .(strain, grp = fifelse(class == "motif", grp_lv[2], grp_lv[3]), log2FoldChange,
                     called, lab = sub("chrM", "", site_id))],
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
  geom_jitter(data = pa[!grp %in% grp_lv[c(1, 5)]], aes(colour = called), width = 0.15, height = 0, size = 1.6) +
  geom_text(data = pa[grp %in% grp_lv[2:3]], aes(label = lab), size = 2, hjust = -0.25, check_overlap = TRUE) +
  scale_colour_manual("padj < 0.05 and\nlog2FC > 1", values = c(`TRUE` = COL_HI, `FALSE` = "grey40")) +
  scale_x_discrete(limits = grp_lv) +
  facet_wrap(~ strain) +
  labs(x = NULL, y = "log2 fold change S / N",
       title = "A. NAD captureSeq S vs N at the chrM ATATAAGTA starts") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

ms[, strain := factor(strain, strain_lv)]
pb <- ms[offset == 0 & is.finite(log2FoldChange)]
pB <- ggplot(pb, aes(log2FoldChange, 100 * get(fA1))) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = cage_n_total, colour = called)) +
  geom_text(aes(label = sub("chrM", "", site_id)), size = 2, vjust = -1.1, check_overlap = TRUE) +
  scale_size_area("CAGE reads", max_size = 4) +
  scale_colour_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey40"), guide = "none") +
  facet_wrap(~ strain) +
  labs(x = "NAD log2 fold change S / N at +1", y = sprintf("CAGE %s share at +1 (%%)", TYPE_A1),
       title = sprintf("B. CAGE %s share against NAD enrichment", TYPE_A1)) +
  theme_bw(base_size = 9)

pC <- ggplot(ms, aes(factor(offset), site_id, fill = log10(meanS + 1))) +
  geom_tile() +
  geom_point(data = ms[called == TRUE], shape = 4, size = 0.8) +
  scale_fill_gradient("log10 mean\nnormalised S + 1", low = "white", high = COL_HI) +
  facet_wrap(~ strain) +
  labs(x = "offset from the final A of ATATAAGTA (nt)", y = NULL,
       title = "C. S reads across the window (x: called)") +
  theme_bw(base_size = 9)

p <- pA / pB / pC + plot_layout(heights = c(1, 1, 1.2))
ggsave(fig_pdf, p, width = 8, height = 11)
message("wrote ", out_tsv, " and ", fig_pdf)

sessionInfo()
