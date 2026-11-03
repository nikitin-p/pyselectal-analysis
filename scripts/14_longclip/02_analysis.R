#!/usr/bin/env Rscript
# Stage 22 step 2: the tables and the figure that answer the long-clip
# questions, from the per-read table of step 1.
#
#   lc_class_summary.tsv       reads per origin class, beside the reversed-clip
#                              control's reads in the same class
#   lc_acceptor_offsets.tsv    aligned 5' end relative to the first base of the
#                              downstream exon, overall and per intron
#   lc_acceptor_composition.tsv  per-position base composition, -20..+5, of the
#                              genome at all annotated acceptors, at the aligned
#                              5' ends of the splice reads and of the unplaced
#                              reads, and of the unplaced clips themselves
#   lc_tss_class.tsv           TSS class of the inferred 5' end, per origin class
#   lc_coproduction.tsv        1Sg reads within +-1 bp of the inferred 5' end,
#                              per origin class, against unclipped (M) reads as
#                              the baseline
#   lc_donor_distance.tsv      inferred 5' end of the splice reads measured from
#                              the donor, negative toward the promoter, beside
#                              the 1Sg 5' ends upstream of the same donors
#   <figdir>/26_longclip_origins.pdf   panels A-D
#
# In the composition table, position -1 is the base immediately upstream of the
# reference point and +1 the reference base itself (the exon's first base, or
# the first aligned base), so the clip of a read occupies -clen..-1.
#
# Usage (from project root):
#   Rscript scripts/14_longclip/02_analysis.R \
#       [--reads results/longclip/longclip_reads.tsv.gz] \
#       [--fine  results/dominance/fine_typed_ctss.tsv.gz] \
#       [--promoters results/tss/promoter_clusters.tsv] \
#       [--fasta data/sacCer3.fa] [--gtf data/sacCer3.ensGene.gtf] \
#       [--outdir results/longclip] [--figdir results/figures] \
#       [--params config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
reads_f  <- get_arg("--reads",     "results/longclip/longclip_reads.tsv.gz")
fine_f   <- get_arg("--fine",      "results/dominance/fine_typed_ctss.tsv.gz")
prom_f   <- get_arg("--promoters", "results/tss/promoter_clusters.tsv")
fasta_f  <- get_arg("--fasta",     "data/sacCer3.fa")
gtf_f    <- get_arg("--gtf",       "data/sacCer3.ensGene.gtf")
outdir   <- get_arg("--outdir",    "results/longclip")
figdir   <- get_arg("--figdir",    "results/figures")
params_f <- get_arg("--params",    "config/params.yaml")
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "lc_helpers.R"))

params <- yaml.load_file(params_f)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
fig_f <- file.path(figdir, "26_longclip_origins.pdf")
if (!force && file.exists(fig_f)) {
  message("Output exists, skipping (use --force to rerun): ", fig_f)
  quit(save = "no", status = 0)
}

reads <- fread(reads_f)
reads[, class := factor(class, levels = LC_CLASSES)]
reads[, tss_class := factor(tss_class, levels = LC_TSS_CLASSES)]
fa <- FaFile(fasta_f)
seqlens <- seqlengths(fa)
UNCLIPPED <- END_TYPES_REAL[!grepl("S", END_TYPES_REAL)]   # the unclipped type, M
stopifnot(length(UNCLIPPED) == 1)

# ── class summary ─────────────────────────────────────────────────────────────

ctrl <- reads[, .(control_reads = sum(n_reads)), by = .(class = control_class)]
summ <- reads[, .(n_reads = sum(n_reads), n_alignments = .N,
                  n_short = sum(n_reads[unplaced_reason %in% "short"]),
                  n_not_found = sum(n_reads[unplaced_reason %in% "not_found"])),
              by = class]
summ <- merge(summ, ctrl[, class := factor(class, levels = LC_CLASSES)],
              by = "class", all.x = TRUE)
summ[is.na(control_reads), control_reads := 0L]
summ[, `:=`(pct_reads = 100 * n_reads / sum(n_reads),
            pct_control = 100 * control_reads / sum(control_reads))]
setorder(summ, class)
fwrite(summ, file.path(outdir, "lc_class_summary.tsv"), sep = "\t")
print(summ)

# ── offsets at the acceptor ───────────────────────────────────────────────────

spl <- reads[class == "splice_annotated" & !is.na(acceptor_offset)]
offs <- rbind(
  spl[, .(n_reads = sum(n_reads)), by = acceptor_offset][, intron_id := "all"],
  spl[, .(n_reads = sum(n_reads)), by = .(intron_id, acceptor_offset)])
offs[, pct_of_intron := 100 * n_reads / sum(n_reads), by = intron_id]
setorder(offs, intron_id, acceptor_offset)
fwrite(offs, file.path(outdir, "lc_acceptor_offsets.tsv"), sep = "\t")

# How consistent the offset is per intron: the share of each intron's reads at
# its own modal offset, for introns with at least 20 reads.
per_intron <- offs[intron_id != "all"][, .(n = sum(n_reads), modal_share = max(pct_of_intron)),
                                       by = intron_id][n >= 20]
message("Introns with >= 20 splice reads: ", nrow(per_intron),
        "; median share at the intron's modal offset: ",
        round(median(per_intron$modal_share), 1), " %")

# ── base composition around the reference points ──────────────────────────────

FLANK_UP <- 20L; FLANK_DN <- 5L
composition <- function(chr, ref, strand, weight, series) {
  # ref is the +1 base; the window runs -FLANK_UP..+FLANK_DN in read direction.
  st <- ifelse(strand == "+", ref - FLANK_UP, ref - FLANK_DN + 1L)
  en <- ifelse(strand == "+", ref + FLANK_DN - 1L, ref + FLANK_UP)
  ok <- st >= 1 & en <= seqlens[chr]
  s  <- as.character(getSeq(fa, GRanges(chr[ok], IRanges(st[ok], en[ok]), strand = strand[ok])))
  m  <- do.call(rbind, strsplit(s, ""))
  pos <- c(-FLANK_UP:-1, 1:FLANK_DN)
  rbindlist(lapply(seq_along(pos), function(j)
    data.table(base = m[, j], w = weight[ok])[, .(w = sum(w)), by = base][
      , .(series = series, position = pos[j], base, frac = w / sum(w))]))
}

gtf <- import(gtf_f)
ex  <- as.data.table(gtf[gtf$type == "exon"])[, .(chr = as.character(seqnames), start, end,
                                                   strand = as.character(strand),
                                                   transcript_id, gene_id)]
setorder(ex, transcript_id, start)
intr <- ex[, if (.N > 1) .(i_start = end[-.N] + 1L, i_end = start[-1] - 1L),
           by = .(transcript_id, gene_id, chr, strand)]
intr <- unique(intr[, .(gene_id, chr, strand, i_start, i_end)])
intr[, intron_id := paste0(chr, ":", i_start, "-", i_end, ":", strand)]
intr[, ex2 := ifelse(strand == "+", i_end + 1L, i_start - 1L)]
intr[, don := ifelse(strand == "+", i_start - 1L, i_end + 1L)]

unp <- reads[class == "unplaced"]
comp <- rbind(
  composition(intr$chr, intr$ex2, intr$strand, rep(1, nrow(intr)), "genome at all annotated acceptors"),
  composition(spl$chr, spl$pos5, spl$strand, spl$n_reads, "genome at splice_annotated aligned starts"),
  composition(unp$chr, unp$pos5, unp$strand, unp$n_reads, "genome at unplaced aligned starts"))
# The unplaced clips themselves, placed so that their last base sits at -1.
clipm <- unp[, .(pos = -(clen:1), base = strsplit(clip, "")[[1]], w = n_reads),
             by = seq_len(nrow(unp))][pos >= -FLANK_UP]
comp <- rbind(comp, clipm[, .(w = sum(w)), by = .(position = pos, base)][
  , .(series = "unplaced clip bases", position, base, frac = w / sum(w)), by = position][
  , .(series, position, base, frac)])
fwrite(comp, file.path(outdir, "lc_acceptor_composition.tsv"), sep = "\t")

g_at <- function(s, p) comp[series == s & position %in% p & base == "G", round(100 * mean(frac), 1)]
message("G content at -20..-1: all acceptors ", g_at("genome at all annotated acceptors", -20:-1),
        " %; unplaced clip bases ", g_at("unplaced clip bases", -20:-1), " %")

# ── TSS class ────────────────────────────────────────────────────────────────

tss <- reads[, .(n_reads = sum(n_reads), n_positions = uniqueN(paste(chr, strand, tss5))),
             by = .(class, tss_class)]
tss[, pct_reads := 100 * n_reads / sum(n_reads), by = class]
setorder(tss, class, tss_class)
fwrite(tss, file.path(outdir, "lc_tss_class.tsv"), sep = "\t")

# ── co-production with 1Sg ────────────────────────────────────────────────────

fine <- fread(fine_f)
sg <- fine[type == END_TYPE_CANONICAL, .(n_sg = sum(n)), by = .(chr, strand, pos)]
co <- reads[, .(n_reads = sum(n_reads),
                pct_with_sg_pm1 = 100 * sum(n_reads[sg_reads_pm1 > 0]) / sum(n_reads),
                pct_sg_dominated_pm1 = 100 * sum(n_reads[sg_dominated_pm1]) / sum(n_reads),
                median_sg_reads_pm1 = as.numeric(median(rep(sg_reads_pm1, n_reads)))),
            by = class]
co[, class := as.character(class)]
# Baseline: unclipped reads, scored with the same +-1 bp rule.
mm <- fine[type == UNCLIPPED, .(n_m = sum(n)), by = .(chr, strand, pos)]
nb <- mm[, .(pos2 = pos + (-1:1)), by = .(chr, strand, pos, n_m)]
nb <- merge(nb, sg, by.x = c("chr", "strand", "pos2"), by.y = c("chr", "strand", "pos"),
            all.x = TRUE)
mb <- nb[, .(sg = sum(n_sg, na.rm = TRUE)), by = .(chr, strand, pos, n_m)]
co <- rbind(co, data.table(class = paste0("baseline: ", UNCLIPPED, " reads"),
                           n_reads = sum(mb$n_m),
                           pct_with_sg_pm1 = 100 * sum(mb$n_m[mb$sg > 0]) / sum(mb$n_m),
                           pct_sg_dominated_pm1 = NA_real_,
                           median_sg_reads_pm1 = as.numeric(median(rep(mb$sg, mb$n_m)))))
fwrite(co, file.path(outdir, "lc_coproduction.tsv"), sep = "\t")
print(co)

# ── the histogram: distance from the donor toward the promoter ────────────────

# The window reaches well past the longest clip (43 nt for a 48 nt read with a
# 5 nt aligned part), so the 1Sg reference shows where ordinary starts sit
# beyond the range a long-clip read can report.
REACH <- 150L
sp <- reads[class %in% c("splice_annotated", "splice_unannotated") & !is.na(donor)]
sp[, dist := ifelse(strand == "+", tss5 - donor, donor - tss5)]
# Reference: every 1Sg 5' end within REACH bases upstream of the same donors,
# on the same strand, so the long-clip ends can be read against where ordinary
# capped starts sit.
dons <- unique(sp[, .(chr, strand, donor)])
ref <- dons[, .(pos = if (strand == "+") (donor - REACH):donor else donor:(donor + REACH)),
            by = .(chr, strand, donor)]
ref <- merge(ref, sg, by = c("chr", "strand", "pos"))
ref[, dist := ifelse(strand == "+", pos - donor, donor - pos)]
# The canonical TSS of each annotated intron's gene, as a rug.
prom  <- fread(prom_f)
canon <- prom[prom[, .I[which.max(reads)], by = gene]$V1, .(gene, tss_pos)]
cr <- merge(unique(sp[class == "splice_annotated", .(intron_id)]), intr, by = "intron_id")
cr <- merge(cr, canon, by.x = "gene_id", by.y = "gene")
cr[, dist := ifelse(strand == "+", tss_pos - don, don - tss_pos)]

dd <- rbind(sp[, .(series = "long-clip inferred 5' end", class, tss_class, dist, n = n_reads)],
            ref[, .(series = "1Sg 5' ends, same donors", class = NA, tss_class = NA,
                    dist, n = n_sg)], fill = TRUE)
fwrite(dd, file.path(outdir, "lc_donor_distance.tsv"), sep = "\t")

# ── figure ────────────────────────────────────────────────────────────────────

th <- theme_classic(base_size = 9) + theme(legend.key.size = unit(0.35, "cm"))

pA <- ggplot(offs[intron_id == "all"], aes(acceptor_offset, n_reads)) +
  geom_col(fill = LC_CLASS_COLORS[["splice_annotated"]], width = 0.8) +
  scale_x_continuous(breaks = seq(-10, 2, 1)) +
  labs(title = "A  Aligned 5' end at annotated acceptors",
       x = "offset from the first base of the downstream exon (negative = inside the intron)",
       y = "reads") + th

lc_sp <- sp[dist >= -REACH]
sg_ref <- ref[dist >= -REACH, .(n = sum(n_sg)), by = dist]
scale_k <- sum(lc_sp$n_reads) / sum(sg_ref$n)
pB <- ggplot() +
  geom_col(data = lc_sp[, .(n = sum(n_reads)), by = .(dist, tss_class)],
           aes(dist, n, fill = tss_class), width = 1) +
  geom_line(data = sg_ref[order(dist)], aes(dist, n * scale_k, colour = "1Sg 5' ends (scaled to the same total)"),
            linewidth = 0.3) +
  geom_rug(data = cr[dist >= -REACH], aes(x = dist), sides = "b", colour = "black",
           length = unit(0.04, "npc")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = LC_TSS_COLORS, drop = FALSE, name = "TSS class") +
  scale_colour_manual(values = c("1Sg 5' ends (scaled to the same total)" = "#dc2626"), name = NULL) +
  labs(title = "B  Inferred 5' end of spliced long-clip reads, from the donor toward the promoter",
       subtitle = paste("0 = last base of exon 1; rug = canonical TSS of each annotated intron's gene;",
                        "a 48 nt read can report a 5' end at most ~43 nt from the donor"),
       x = "position relative to the 5' splice site (nt, read strand)", y = "reads") + th

pC <- ggplot(comp[base %in% c("A", "C", "G", "T")], aes(position, 100 * frac, colour = series)) +
  geom_line(linewidth = 0.4) + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
  facet_wrap(~ base, nrow = 1) +
  scale_colour_manual(values = c("genome at all annotated acceptors" = "grey50",
                                 "genome at splice_annotated aligned starts" = LC_CLASS_COLORS[["splice_annotated"]],
                                 "genome at unplaced aligned starts" = "#f87171",
                                 "unplaced clip bases" = LC_CLASS_COLORS[["unplaced"]]), name = NULL) +
  labs(title = "C  Base composition around the reference point (-1 = base upstream of it)",
       x = "position", y = "% of reads") + th + theme(legend.position = "bottom")

d1 <- summ[, .(class, reads = n_reads, control = control_reads)]
d1 <- melt(d1, id.vars = "class", variable.name = "set", value.name = "n")
pD1 <- ggplot(d1, aes(class, n, fill = class, alpha = set)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = LC_CLASS_COLORS, guide = "none") +
  scale_alpha_manual(values = c(reads = 1, control = 0.35), name = NULL,
                     labels = c("reads", "reversed-clip control")) +
  labs(title = "D  Origin classes", x = NULL, y = "reads") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
pD2 <- ggplot(tss, aes(class, pct_reads, fill = tss_class)) + geom_col(width = 0.75) +
  scale_fill_manual(values = LC_TSS_COLORS, name = "TSS class") +
  labs(title = "E  TSS class of the inferred 5' end", x = NULL, y = "% of reads") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
cod <- co[, .(class, pct_with_sg_pm1)]
cod[, class := factor(class, levels = c(LC_CLASSES, cod$class[!cod$class %in% LC_CLASSES]))]
pD3 <- ggplot(cod, aes(class, pct_with_sg_pm1)) + geom_col(fill = "grey30", width = 0.75) +
  labs(title = "F  Reads with any 1Sg read within 1 bp of their 5' end", x = NULL, y = "% of reads") +
  th + theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig <- (pA | pB) / pC / (pD1 | pD2 | pD3) + plot_layout(heights = c(1, 0.9, 1))
ggsave(fig_f, fig, width = 15, height = 12)
message("Written: ", fig_f)
sessionInfo()
