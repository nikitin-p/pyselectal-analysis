#!/usr/bin/env Rscript
# Stage 22 step 10: why spliced long-clip reads align one or more bases inside
# the intron.
#
# For every spliced long-clip read (class splice_*), the aligned 5' end is
# expressed as an offset from the first base of exon 2 (0 = exon start,
# negative = inside the intron), and k is the length of the suffix that the
# last K = 12 bases of exon 1 and of the intron share. STAR moves into the
# intron exactly the bases that read the same on both sides of the junction.
#   A  offset by the last base of exon 1 (the intron always ends in AG, so a
#      -1 shift occurs only when exon 1 ends in G);
#   B  offset by k;
#   C  among unshifted reads, the last clip base against the last base of exon 1.
#
# Outputs: <outdir>/lc_splice_shift_reads.tsv (one row per alignment);
#          <figdir>/30_longclip_splice_shift.{pdf,png}
#
# Usage (from project root):
#   Rscript scripts/14_longclip/10_splice_shift.R \
#       [--reads results/longclip/longclip_reads.tsv.gz] [--fasta data/sacCer3.fa] \
#       [--outdir results/longclip] [--figdir results/figures] [--force]
suppressPackageStartupMessages({library(data.table);library(Rsamtools);library(Biostrings)
  library(ggplot2);library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
reads_f <- get_arg("--reads",  "results/longclip/longclip_reads.tsv.gz")
fasta_f <- get_arg("--fasta",  "data/sacCer3.fa")
outdir  <- get_arg("--outdir", "results/longclip")
figdir  <- get_arg("--figdir", "results/figures")
force   <- "--force" %in% args
fig_f <- file.path(figdir, "30_longclip_splice_shift.pdf")
if (file.exists(fig_f) && !force) { message("[skip] ", fig_f, " exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

r <- fread(reads_f)[grepl("^splice", class)]
r[, c("ic","is","ie","ist") := tstrsplit(intron_id, "[:-]")][, `:=`(is=as.integer(is), ie=as.integer(ie))]
r[, ex2 := ifelse(strand=="+", ie+1L, is-1L)]
r[, off := ifelse(strand=="+", pos5-ex2, ex2-pos5)]
fa <- FaFile(fasta_f)
gs <- function(chr, s, e, st) { x <- getSeq(fa, GRanges(chr, IRanges(s, e), strand=st)); as.character(x) }
K <- 12L
# exon-1 tail (K bases ending at donor) and intron tail (K bases ending at last intron base), read strand
r[, ex1 := gs(chr, ifelse(strand=="+", donor-K+1L, donor), ifelse(strand=="+", donor, donor+K-1L), strand)]
r[, itl := gs(chr, ifelse(strand=="+", ie-K+1L, is), ifelse(strand=="+", ie, is+K-1L), strand)]
sufmatch <- function(a, b) mapply(function(x, y) { x <- rev(strsplit(x,"")[[1]]); y <- rev(strsplit(y,"")[[1]]); m <- which(x != y); if (length(m)) m[1]-1L else length(x) }, a, b)
r[, k := sufmatch(ex1, itl)]
r[, ex1_last := substr(ex1, K, K)][, intr_last2 := substr(itl, K-1, K)]
r[, clip_last := substr(clip, clen, clen)]
cat("splice reads:", sum(r$n_reads), "\n")
cat("\nacceptor (last 2 intron bases):\n"); print(r[, .(n=sum(n_reads)), by=intr_last2][order(-n)])
cat("\noffset by last base of exon 1 (read-weighted %):\n")
t1 <- r[, .(n=sum(n_reads)), by=.(ex1_last, off)][, pct := round(100*n/sum(n),1), by=ex1_last]
print(dcast(t1, ex1_last ~ off, value.var="pct", fill=0))
print(r[, .(n_reads=sum(n_reads), n_introns=uniqueN(intron_id)), by=ex1_last])
cat("\noffset == -k (exon1/intron shared suffix):", round(100*r[off == -k, sum(n_reads)]/sum(r$n_reads),1), "%\n")
cat("offset == -1 :", round(100*r[off == -1, sum(n_reads)]/sum(r$n_reads),1), "%\n")
print(r[, .(n=sum(n_reads)), by=.(k=pmin(k,6), off)][order(k, off)][, pct_k := round(100*n/sum(n),1), by=k][])
r <- r[, .(intron_id, class, ex1_last, intr_last2, k, off, clip_last, clip, n_reads, n_untemplated, identity)]
fwrite(r, file.path(outdir, "lc_splice_shift_reads.tsv"), sep="\t")

th <- theme_classic(base_size=9) + theme(plot.title=element_text(face="bold", size=10),
  plot.subtitle=element_text(size=8, colour="grey30"), axis.line=element_blank(), axis.ticks=element_blank())
tile <- function(d, x, y, xl, yl, title, sub, rowlab) {
  d <- copy(d)[, pct := 100*n/sum(n), by=y]
  nl <- d[, .(n=sum(n)), by=y]
  ggplot(d, aes(factor(.data[[x]]), factor(.data[[y]], levels=rev(sort(unique(as.character(d[[y]]))))), fill=pct)) +
    geom_tile(colour="white", linewidth=0.8) +
    geom_text(data=d[pct >= 1], aes(label=sprintf("%.0f", pct), colour=pct > 50), size=2.6) +
    geom_text(data=nl, aes(x=Inf, y=factor(.data[[y]], levels=rev(sort(unique(as.character(d[[y]]))))), label=sprintf("  n = %s", format(n, big.mark=",", trim=TRUE))),
              inherit.aes=FALSE, hjust=0, size=2.5, colour="grey30") +
    scale_colour_manual(values=c(`TRUE`="white", `FALSE`="grey15"), guide="none") +
    scale_fill_gradient(low="#eff6ff", high="#1d4ed8", limits=c(0,100), name="% of row") +
    coord_cartesian(clip="off") + labs(x=xl, y=yl, title=title, subtitle=sub) + th +
    theme(plot.margin=margin(5, 55, 5, 5))
}
OFFX <- "aligned 5' end relative to the first base of exon 2\n(0 = exon start, negative = inside the intron)"
pA <- tile(r[, .(n=sum(n_reads)), by=.(ex1_last, off)], "off", "ex1_last", OFFX, "last base of exon 1",
  "A  A shift of -1 occurs only when exon 1 ends in G",
  "the intron always ends in AG; a cap G next to exon 2 would match it at every intron")
r[, k6 := ifelse(k >= 4, "4+", as.character(k))]
pB <- tile(r[, .(n=sum(n_reads)), by=.(k6, off)], "off", "k6", OFFX, "bases shared by the tails\nof exon 1 and the intron",
  "B  The shift equals the shared suffix of exon 1 and the intron",
  "STAR moves into the intron exactly the bases that read the same on both sides of the junction")
pC <- tile(r[off == 0, .(n=sum(n_reads)), by=.(ex1_last, clip_last)], "clip_last", "ex1_last",
  "last clip base (next to exon 2)", "last base of exon 1",
  "C  Unshifted reads: the base next to exon 2 is exon 1's, not a cap G",
  "reads with offset 0; G-ending exons have only 8 such reads")
fig <- (pA / pB / pC) + plot_layout(heights=c(4, 5, 4), guides="collect") & theme(legend.position="bottom", legend.key.height=unit(0.3,"cm"))
ggsave(fig_f, fig, width=7, height=9)
ggsave(sub("\\.pdf$", ".png", fig_f), fig, width=7, height=9, dpi=150)
message("Written: ", fig_f)
sessionInfo()
