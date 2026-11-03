#!/usr/bin/env Rscript
# Stage 22 step 6: figure 27_longclip_m_counterparts.pdf — the unclipped (M)
# counterparts of the long-clip reads, and why STAR clips some spliced reads.
#
#   A  one locus: long-clip reads as STAR aligned them, the same reads realigned
#      by step 4, and split M reads with the same 5' end;
#   B  share of long-clip reads with an M read starting d bases from their
#      inferred 5' end, per origin class (d = 0 is the counterpart, d != 0 the
#      local background);
#   C  share of spliced reads that STAR soft-clipped rather than split, by the
#      length of the exon-1 piece, for reads with and without an extra 5' base.
#
# Usage (from project root):
#   Rscript scripts/14_longclip/06_m_counterparts_plot.R \
#       [--indir results/longclip] [--genes data/sacCer3.genes.bed] \
#       [--figdir results/figures] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
indir  <- get_arg("--indir",  "results/longclip")
genes_f <- get_arg("--genes", "data/sacCer3.genes.bed")
figdir <- get_arg("--figdir", "results/figures")
force  <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "lc_helpers.R"))

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
fig_f <- file.path(figdir, "27_longclip_m_counterparts.pdf")
if (!force && file.exists(fig_f)) {
  message("Output exists, skipping (use --force to rerun): ", fig_f)
  quit(save = "no", status = 0)
}

prof <- fread(file.path(indir, "lc_m_profile.tsv"))
ex1  <- fread(file.path(indir, "lc_m_exon1.tsv"))
exm  <- fread(file.path(indir, "lc_m_example.tsv"))

th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.key.size = unit(0.35, "cm"))

PART_COLORS <- c(soft_clip = "#dc2626", untemplated = "#f59e0b",
                 exon1_piece = "#1d4ed8", aligned = "#6b7280")
PART_LABELS <- c(soft_clip = "soft-clipped by STAR", untemplated = "untemplated 5' base (cap G)",
                 exon1_piece = "exon-1 piece", aligned = "exon 2")

# ── A: the example locus ──────────────────────────────────────────────────────

t5 <- exm$tss5[1]; st <- exm$strand[1]; ex_chr <- exm$chr[1]
iv <- as.integer(strsplit(sub("^[^:]+:([0-9]+-[0-9]+):.*$", "\\1", exm$intron_id[1]), "-")[[1]])
genes <- fread(genes_f, header = FALSE, col.names = c("chr", "start", "end", "name", "score", "strand"))
gene <- genes[chr == ex_chr & strand == st & start < iv[1] & end > iv[2], name]
gene <- if (length(gene)) gene[1] else "gene"

# Read-direction distance from the 5' end, with the middle of the intron cut out
# so that both of its ends stay at full scale.
sgn <- if (st == "+") 1 else -1
d_of <- function(p) sgn * (p - t5)
i_lo <- min(d_of(iv)); i_hi <- max(d_of(iv))   # intron in distance units
KEEP <- 25; CUT <- 8
cut_lo <- i_lo + KEEP; cut_hi <- i_hi - KEEP
squeeze <- function(d) ifelse(d < cut_lo, d,
                       ifelse(d > cut_hi, d - (cut_hi - cut_lo) + CUT,
                              cut_lo + (d - cut_lo) / (cut_hi - cut_lo) * CUT))
exm[, `:=`(x0 = squeeze(pmin(d_of(start), d_of(end))) - 0.5,
           x1 = squeeze(pmax(d_of(start), d_of(end))) + 0.5)]

views <- c("STAR alignment", "realigned", "split M read")
view_lab <- c("STAR alignment" = "long-clip reads,\nas STAR aligned them",
              "realigned" = "the same reads,\nclip placed on exon 1",
              "split M read" = "M reads from the\nsame 5' end")
exm[, view := factor(view, views)]
exm[, row := paste(view, read)]
rows <- unique(exm[order(view, read), row])
exm[, y := match(row, rev(rows))]
gaps <- exm[, .(ymid = mean(y), ytop = max(y)), by = view]
exm[, y := y + 0.9 * (length(views) - as.integer(view))]   # space between groups
lab <- exm[, .(y = mean(y)), by = view]

x_e1 <- squeeze(i_lo) - 0.5; x_e2 <- squeeze(i_hi) + 0.5
x_cut <- c(squeeze(cut_lo), squeeze(cut_hi))
intr <- exm[part == "intron"]
blk  <- exm[part != "intron"]
ytop <- max(exm$y) + 1
n_ex <- if (nrow(blk[view == "realigned" & part == "exon1_piece"]))
  blk[view == "realigned" & part == "exon1_piece", round(x1[1] - x0[1])] else NA

pA <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = x_e1, ymin = -Inf, ymax = Inf, fill = "#dbeafe", alpha = 0.5) +
  annotate("rect", xmin = x_e2, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92", alpha = 0.6) +
  annotate("text", x = x_e1 / 2, y = ytop, label = "exon 1", size = 2.8, fontface = "bold") +
  annotate("text", x = (x_e1 + x_e2) / 2, y = ytop,
           label = sprintf("intron, %d nt (middle cut out)", iv[2] - iv[1] + 1), size = 2.8) +
  annotate("text", x = x_e2 + 20, y = ytop, label = "exon 2", size = 2.8, fontface = "bold") +
  annotate("segment", x = x_cut, xend = x_cut, y = 0.3, yend = ytop - 0.5,
           linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_segment(data = intr, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.3) +
  geom_rect(data = blk, aes(xmin = x0, xmax = x1, ymin = y - 0.32, ymax = y + 0.32, fill = part)) +
  geom_text(data = lab, aes(x = min(exm$x0) - 2, y = y, label = view_lab[as.character(view)]),
            hjust = 1, size = 2.6, lineheight = 0.9) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = PART_COLORS, labels = PART_LABELS, name = NULL,
                    breaks = names(PART_COLORS)) +
  scale_x_continuous(expand = expansion(mult = c(0.32, 0.02))) +
  scale_y_continuous(expand = expansion(add = 0.6)) +
  labs(title = "A  A long-clip read is a spliced mRNA whose short exon-1 piece STAR soft-clipped",
       subtitle = sprintf(paste0("%s, %s:%d (%s strand); reads 5'->3' left to right. ",
                                 "STAR placed the clipped reads a few bases inside the intron; ",
                                 "realigned, they carry the same %s-nt exon-1 piece as the split M reads, plus one untemplated G."),
                          gene, ex_chr, t5, st, n_ex),
       x = NULL, y = NULL) +
  th + theme(axis.line = element_blank(), axis.text = element_blank(),
             axis.ticks = element_blank(), legend.position = "bottom")

# ── B: M reads around the long-clip 5' end ────────────────────────────────────

prof[, class := factor(class, LC_CLASSES)]
bg <- prof[abs(offset) >= 5, .(bg = mean(pct)), by = .(class, kind)]
at0 <- prof[offset == 0, .(class, kind, pct)]
fold <- merge(at0, bg, by = c("class", "kind"))
fold[, label := sprintf("d = 0: %.0f%%;  flanks: %.0f%%", pct, bg)]
KIND_COL <- c(any_M = "grey20", spliced_M_same_intron = "#1d4ed8")
KIND_LAB <- c(any_M = "any M read", spliced_M_same_intron = "M read split across the same intron")
fold[, y := ifelse(kind == "any_M", 122, 112)]
note <- data.table(class = factor(c("splice_annotated", "splice_unannotated"), LC_CLASSES),
                   x = 11, y = 70,
                   label = "d > 0: exon-1 piece too\nshort for STAR to split (C)")

pB <- ggplot(prof, aes(offset, pct, colour = kind)) +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_line(linewidth = 0.4) + geom_point(size = 0.5) +
  geom_text(data = fold, aes(x = -20, y = y, label = label), hjust = 0, size = 2.3,
            show.legend = FALSE) +
  geom_text(data = note, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            size = 2.2, colour = "#1d4ed8", lineheight = 0.9) +
  facet_wrap(~ class, nrow = 1) +
  scale_colour_manual(values = KIND_COL, labels = KIND_LAB, name = NULL) +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 100, 25)) +
  labs(title = "B  Spliced long clips have split M reads at the same 5' end; templated_distal clips have none, adjacent and unplaced clips a modest peak",
       subtitle = paste0("Share of long-clip reads with at least one M read starting d bases from their inferred 5' end (the first templated clip base; the aligned 5' end for unplaced).\n",
                         "Flanks = mean over |d| >= 5; they sit inside the same TSS clusters, so they are not a neutral background."),
       x = "d, bases from the long-clip read's inferred 5' end (negative = upstream along the read strand)",
       y = "% of long-clip reads") +
  th + theme(legend.position = "bottom", strip.background = element_blank(),
             strip.text = element_text(face = "bold"))

# ── C: when does STAR clip rather than split? ─────────────────────────────────

ex1[, group := ifelse(kind %in% c("clipped_extra5", "split_short_clip"),
                      "extra 5' base (e.g. cap G)", "starts on the exon-1 base")]
ex1[, clipped := startsWith(kind, "clipped")]
sw <- ex1[exon1_len <= 25, .(n = sum(n_reads), n_clip = sum(n_reads[clipped])),
          by = .(group, exon1_len)][n >= 20]
sw[, pct := 100 * n_clip / n]
thr <- sw[pct < 50, .(thr = min(exon1_len)), by = group]
GROUP_COL <- c("starts on the exon-1 base" = "grey20", "extra 5' base (e.g. cap G)" = "#f59e0b")

pC <- ggplot(sw, aes(exon1_len, pct, colour = group)) +
  geom_hline(yintercept = 50, linetype = "dotted", colour = "grey70") +
  geom_line(linewidth = 0.5) + geom_point(aes(size = n)) +
  scale_size_area(max_size = 3.5, name = "reads", breaks = c(100, 1000, 10000, 50000),
                  labels = scales::comma) +
  scale_colour_manual(values = GROUP_COL, name = NULL) +
  scale_x_continuous(breaks = seq(0, 25, 5)) +
  labs(title = "C  STAR splits once the exon-1 piece is long enough, and an extra 5' G raises that length by 2-4 nt",
       subtitle = sprintf(paste0("Reads through the %d introns that splice-class long clips use. ",
                                 "Fewer than half are clipped from %s; lengths with < 20 reads omitted."),
                          uniqueN(ex1$intron_id),
                          paste(sprintf("%d nt when the read %s", thr$thr,
                                        ifelse(startsWith(thr$group, "extra"), "has an extra 5' base",
                                               "starts on the exon-1 base")), collapse = " and ")),
       x = "length of the exon-1 piece (templated bases upstream of the intron)",
       y = "% soft-clipped by STAR\n(rest split across the intron)") +
  th + theme(legend.position = "right")

fig <- pA / pB / pC + plot_layout(heights = c(1.1, 1, 1))
ggsave(fig_f, fig, width = 13, height = 12)
message("Written: ", fig_f)

fwrite(fold[, .(class, kind, pct_at_0 = pct, pct_background = round(bg, 2))],
       file.path(indir, "lc_m_counterpart_summary.tsv"), sep = "\t")
fwrite(sw, file.path(indir, "lc_m_clip_vs_split.tsv"), sep = "\t")
message("Written: ", file.path(indir, "lc_m_counterpart_summary.tsv"), ", ",
        file.path(indir, "lc_m_clip_vs_split.tsv"))
sessionInfo()
