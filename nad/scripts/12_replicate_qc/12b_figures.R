#!/usr/bin/env Rscript
# NAD branch step 12b: PCA and leave-one-out figures from the step-12 tables.
#
# Draws from pca.tsv and loo.tsv only, so it runs on the Mac without DESeq2.
#   12_pca.pdf  PC1 x PC2 per strain, axes labelled with % variance
#   12_loo.pdf  calls lost (down) and gained (up) relative to the full fit for
#               each dropped run, so a difference of 287 vs 293 calls is visible
#
# Usage (from project root):
#   Rscript nad/scripts/12_replicate_qc/12b_figures.R
#     [--indir nad/results/12_replicate_qc] [--figdir nad/figures]

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork); library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir", default = "nad/results/12_replicate_qc"),
  make_option("--figdir", default = "nad/figures"))))
for (f in c("pca.tsv", "loo.tsv"))
  if (!file.exists(file.path(opt$indir, f))) stop("missing ", file.path(opt$indir, f))
dir.create(opt$figdir, showWarnings = FALSE, recursive = TRUE)

lib_col <- c(S = "#e63946", N = "#457b9d")
strain_lab <- c(WT = "WT", dxo1_rai1_npy1 = "dxo1 rai1 npy1")
strains <- c("WT", "dxo1_rai1_npy1")

# ── PCA ─────────────────────────────────────────────────────────────────────
pca <- fread(file.path(opt$indir, "pca.tsv"))
pca_panel <- function(st) {
  d <- pca[strain == st]
  pve <- as.numeric(strsplit(d$pve[1], ",")[[1]])
  ggplot(d, aes(PC1, PC2, colour = library, label = sub("_", " ", lab))) +
    geom_point(size = 2.5) +
    geom_text_repel(size = 3, min.segment.length = 0, box.padding = 0.4, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = lib_col, name = "library") +
    labs(title = strain_lab[[st]],
         x = sprintf("PC1 (%.0f %% of variance)", pve[1]), y = sprintf("PC2 (%.0f %% of variance)", pve[2])) +
    theme_bw(base_size = 11)
}
p <- wrap_plots(lapply(strains, pca_panel), nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Replicates on PCA of VST counts, nuclear unmasked positions",
    subtitle = "PC1 separates S from N in both strains; the spread inside a library is the replicate noise.\nAxes have different scales: compare distances using the % of variance on each axis.")
ggsave(file.path(opt$figdir, "12_pca.pdf"), p, width = 10, height = 4.8)

# ── leave-one-out ───────────────────────────────────────────────────────────
loo <- fread(file.path(opt$indir, "loo.tsv"))
loo[, kept := round(kept_of_full * full_enriched)]
loo[, lost := full_enriched - kept]
stopifnot(loo[, all(kept + gained == enriched)])
loo[, library := substr(lab, 1, 1)]
loo[, lab2 := sub("_", " ", lab)]
loo_panel <- function(st) {
  d <- loo[strain == st]
  d[, lab2 := factor(lab2, levels = d[order(library != "S", lab)]$lab2)]
  long <- rbind(d[, .(lab2, library, kind = "lost", n = -lost)], d[, .(lab2, library, kind = "gained", n = gained)])
  lim <- max(d$lost, d$gained) * 1.25
  ggplot(long, aes(lab2, n, fill = library, alpha = kind)) +
    geom_col(width = 0.7) + geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(data = d, aes(lab2, -lost, label = sprintf("%d calls\n(J = %.2f)", enriched, jaccard)),
              inherit.aes = FALSE, vjust = 1.2, size = 2.8, lineheight = 0.9) +
    geom_text(data = d[gained > 0], aes(lab2, gained, label = paste0("+", gained)),
              inherit.aes = FALSE, vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = lib_col, guide = "none") +
    scale_alpha_manual(values = c(lost = 1, gained = 0.45), guide = "none") +
    scale_y_continuous(limits = c(-lim * 1.35, lim * 0.6)) +
    labs(title = sprintf("%s: full fit %s calls", strain_lab[[st]], format(d$full_enriched[1], big.mark = ",")),
         x = "run dropped", y = "calls relative to the full fit") +
    theme_bw(base_size = 11)
}
p <- wrap_plots(lapply(strains, loo_panel), nrow = 1) +
  plot_annotation(
    title = "Leave-one-out DESeq2: calls lost (solid, below 0) and gained (pale, above 0) without each run",
    subtitle = "A noisy run would gain calls when dropped; every run here mostly loses calls, so each one adds signal.\nJ = Jaccard overlap with the full-fit calls.")
ggsave(file.path(opt$figdir, "12_loo.pdf"), p, width = 10, height = 4.8)

message("wrote ", file.path(opt$figdir, c("12_pca.pdf", "12_loo.pdf")), collapse = " ")
sessionInfo()
