#!/usr/bin/env Rscript
# NAD branch step 1: the soft-clipped-A end type on the mitochondrial genome.
#
# Mitochondrial transcripts carry no m7G cap, yet yeast mtRNAP initiates with
# NAD efficiently enough that up to 60 % of some mitochondrial mRNAs are
# NAD-capped (Bird et al. 2018). If an NAD cap is read by CAGE as a single
# untemplated A, the single-A clip type should therefore be strongly enriched on
# chrM and sit on mtRNAP start sites. This script tests that from the existing
# typed CTSS table alone, without any NAD captureSeq data.
#
# A soft-clipped A can never sit on a genomic A at -1 (see end_types.R), so a
# fair nuclear baseline must be restricted to the same start context. Every
# comparison of the single-A clip is therefore made on Y|A starts (+1 = A,
# -1 = C or T). It is further stratified by the length of the A run beginning at
# +1, because the nuclear single-A clip rises steeply with that run, which is
# the signature of slippage on an A tract rather than of a cap.
#
# Panels:
#   A  end-type composition of all reads, nuclear versus chrM;
#   B  single-A clip share at Y|A starts by A-run length, nuclear versus chrM,
#      as the pooled value with the range over libraries;
#   C  the chrM positions with the most single-A clip reads, their clip share,
#      depth, nearest gene and sequence around +1;
#   D  sequence logos of chrM positions rich in the single-A clip against
#      nuclear Y|A positions where that clip exceeds 20 % of reads.
#
# Usage (from project root):
#   Rscript nad/scripts/01_chrM_end_types/01_chrM_end_types.R \
#       [--ctss   results/tss/typed_ctss.tsv.gz] \
#       [--params config/params.yaml]            \
#       [--outdir nad/results/01_chrM_end_types] \
#       [--figdir nad/figures] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggseqlogo)
  library(yaml)
  library(GenomicRanges)
  library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
ctss_f   <- get_arg("--ctss",   "results/tss/typed_ctss.tsv.gz")
params_f <- get_arg("--params", "config/params.yaml")
outdir   <- get_arg("--outdir", "nad/results/01_chrM_end_types")
figdir   <- get_arg("--figdir", "nad/figures")
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "..", "..", "scripts", "utils", "end_types.R"))
source(file.path(script_dir, "..", "..", "..", "scripts", "utils", "plot_helpers.R"))

fig_pdf <- file.path(figdir, "01_chrM_single_A_clip.pdf")
if (file.exists(fig_pdf) && !force) {
  message("output exists, skipping (use --force): ", fig_pdf)
  quit(save = "no", status = 0)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)
fa     <- FaFile(params$genome_fasta)
gtf_f  <- params$annotation_gtf

# The single-A clip is read from the vocabulary rather than written as a
# literal: the one-base clip whose genome-adjacent base is A.
TYPE_A1 <- END_TYPES_REAL[END_TYPE_CLIP_BASE[END_TYPES_REAL] %in% "A" &
                          startsWith(END_TYPES_REAL, "1S")]
stopifnot(length(TYPE_A1) == 1)
col_A1  <- paste0("n_", TYPE_A1)

COL_NUC  <- "#8a8f98"
COL_MITO <- "#d9681a"
CTX_UP   <- 10   # bases upstream of +1 kept for context and logos
CTX_DN   <- 5    # bases from +1 onward
ARUN_CAP <- 5    # A runs of this length or longer share one stratum

# ── load and annotate positions ───────────────────────────────────────────────
ctss  <- fread(ctss_f)
ncols <- end_type_cols("n_")
libs  <- sort(unique(ctss$sample_id))
message("libraries: ", length(libs), " (", paste(libs, collapse = ", "), ")")

ctss[, compartment := fifelse(chr == "chrM", "chrM", "nuclear")]

pos <- ctss[, c(lapply(.SD, sum), list(n_total = sum(n_total))),
            by = .(chr, pos, strand, compartment), .SDcols = ncols]
sl  <- seqlengths(fa)
pos <- pos[pos > CTX_UP & pos <= sl[chr] - CTX_UP]

# Sequence from -CTX_UP to +CTX_DN in transcript orientation, so character
# CTX_UP is -1 and character CTX_UP + 1 is +1 on either strand.
ctx_gr <- GRanges(pos$chr,
                  IRanges(ifelse(pos$strand == "+", pos$pos - CTX_UP, pos$pos - CTX_DN + 1),
                          width = CTX_UP + CTX_DN),
                  strand = pos$strand)
pos[, ctx := as.character(getSeq(fa, ctx_gr))]
pos[, `:=`(m1 = substr(ctx, CTX_UP, CTX_UP), p1 = substr(ctx, CTX_UP + 1, CTX_UP + 1))]
pos[, arun := nchar(sub("^(A*).*$", "\\1", substr(ctx, CTX_UP + 1, CTX_UP + CTX_DN)))]
pos[, ya := p1 == "A" & m1 %in% c("C", "T")]
key <- pos[, .(chr, pos, strand, ctx, m1, p1, arun, ya)]

# ── A: composition of all reads ───────────────────────────────────────────────
comp <- melt(pos[, lapply(.SD, sum), by = compartment, .SDcols = c(ncols, "n_total")],
             id.vars = c("compartment", "n_total"), variable.name = "type", value.name = "n")
comp[, type := factor(sub("^n_", "", type), levels = END_TYPES)]
comp[, pct := 100 * n / n_total]
fwrite(comp, file.path(outdir, "composition.tsv"), sep = "\t")

comp_w <- dcast(comp, type ~ compartment, value.var = "pct")
comp_w[, ratio := chrM / nuclear]

# ── B: single-A clip at Y|A starts by A-run length ────────────────────────────
lib_ya <- merge(ctss[, .(chr, pos, strand, compartment, sample_id, nA1 = get(col_A1), n_total)],
                key[ya == TRUE], by = c("chr", "pos", "strand"))
lib_ya[, stratum := pmin(arun, ARUN_CAP)]

strata_lib <- lib_ya[, .(reads = sum(n_total), nA1 = sum(nA1)),
                     by = .(compartment, stratum, sample_id)][, pct := 100 * nA1 / reads]
strata <- lib_ya[, .(sites = uniqueN(paste(chr, pos, strand)), reads = sum(n_total),
                     nA1 = sum(nA1)), by = .(compartment, stratum)]
strata[, pct := 100 * nA1 / reads]
strata <- merge(strata, strata_lib[, .(lo = min(pct), hi = max(pct)), by = .(compartment, stratum)],
                by = c("compartment", "stratum"))
setorder(strata, compartment, stratum)
fwrite(strata, file.path(outdir, "single_A_clip_by_A_run.tsv"), sep = "\t")
fwrite(strata_lib, file.path(outdir, "single_A_clip_by_A_run_per_library.tsv"), sep = "\t")

nuc_ya_pct <- pos[compartment == "nuclear" & ya == TRUE, 100 * sum(get(col_A1)) / sum(n_total)]
mito_ya_pct <- pos[compartment == "chrM" & ya == TRUE, 100 * sum(get(col_A1)) / sum(n_total)]

# ── C: the chrM positions richest in the single-A clip ────────────────────────
# Each position is labelled by the chrM transcript it lies in, or else by the
# nearest transcript starting downstream on the same strand. The ensGene GTF
# names mitochondrial genes by systematic ID only, so the common names of the
# genes that can appear here are supplied below.
MITO_NAMES <- c(Q0020 = "15S rRNA", Q0045 = "COX1", Q0050 = "COX1", Q0055 = "COX1",
                Q0060 = "COX1", Q0065 = "COX1", Q0070 = "COX1", Q0080 = "ATP8",
                Q0085 = "ATP6", Q0105 = "COB", Q0110 = "COB", Q0115 = "COB",
                Q0120 = "COB", Q0130 = "ATP9", Q0140 = "VAR1", Q0158 = "21S rRNA",
                Q0250 = "COX2", Q0275 = "COX3")
gtf <- fread(cmd = paste("awk -F'\\t' '$1==\"chrM\" && $3==\"transcript\"'", shQuote(gtf_f)),
             header = FALSE, sep = "\t")
tx  <- gtf[, .(start = V4, end = V5, strand = V7,
               id = sub('.*gene_id "([^"]+)".*', "\\1", V9))]
tx[, name := fifelse(id %in% names(MITO_NAMES), MITO_NAMES[id], id)]
tx  <- unique(tx[, .(start, end, strand, name)])

label_site <- function(p, s) {
  t  <- tx[strand == s]
  tss <- if (s == "+") t$start else t$end
  d   <- if (s == "+") tss - p else p - tss            # > 0: gene starts downstream
  at  <- which(abs(d) <= 3)
  if (length(at)) return(paste(t$name[at[1]], "5' end"))
  inside <- which(t$start <= p & t$end >= p)
  if (length(inside)) return(paste(t$name[inside[1]], "internal"))
  up <- which(d > 0 & d <= 1500)
  if (length(up)) { i <- up[which.min(d[up])]; return(sprintf("%s, %d nt upstream", t$name[i], d[i])) }
  "intergenic"
}

top <- pos[compartment == "chrM" & get(col_A1) >= 30]
top[, `:=`(nA1 = get(col_A1), frac = 100 * get(col_A1) / n_total)]
top[, gene := mapply(label_site, pos, strand)]
# The mtRNAP promoter nonanucleotide ATATAAGTA initiates at its final A; its
# conserved core AAGTA ending at +1 marks a promoter start.
top[, promoter_core := substr(ctx, CTX_UP - 3, CTX_UP + 1) == "AAGTA"]
top[, seq_lab := paste0(substr(ctx, CTX_UP - 5, CTX_UP), "|", substr(ctx, CTX_UP + 1, CTX_UP + 3))]
top[, site := sprintf("%s %s  %s", format(pos, big.mark = ","), strand, gene)]
setorder(top, -nA1)
fwrite(top[, .(chr, pos, strand, nA1, n_total, frac, gene, promoter_core, ctx)],
       file.path(outdir, "chrM_top_single_A_clip_sites.tsv"), sep = "\t")

# ── D: logo sets ──────────────────────────────────────────────────────────────
logo_mito <- pos[compartment == "chrM" & get(col_A1) >= 20, ctx]
logo_nuc  <- pos[compartment == "nuclear" & ya == TRUE & n_total >= 20 &
                 get(col_A1) / n_total > 0.2, ctx]

# ── figure ────────────────────────────────────────────────────────────────────
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major.y = element_blank(),
                plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(colour = "grey35", size = 8.5),
                legend.position = "top", legend.title = element_blank()))
comp_cols <- c(nuclear = COL_NUC, chrM = COL_MITO)

pa_d <- copy(comp)[, type := factor(type, levels = rev(END_TYPES))]
pa_seg <- dcast(pa_d, type ~ compartment, value.var = "pct")
pa_seg <- merge(pa_seg, comp_w[, .(type, ratio)], by = "type")
p_a <- ggplot(pa_d, aes(pct, type)) +
  geom_segment(data = pa_seg, aes(x = nuclear, xend = chrM, y = type, yend = type),
               colour = "grey80", linewidth = 0.8) +
  geom_point(aes(colour = compartment), size = 2.6) +
  geom_text(data = pa_seg, aes(x = 250, y = type, label = sprintf("×%.2g", ratio)),
            hjust = 1, size = 2.8,
            colour = ifelse(pa_seg$type == TYPE_A1, "black", "grey45"),
            fontface = ifelse(pa_seg$type == TYPE_A1, "bold", "plain")) +
  scale_x_log10(labels = function(x) paste0(x, "%"), breaks = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  coord_cartesian(xlim = c(0.0008, 250)) +
  scale_colour_manual(values = comp_cols,
                      labels = c(nuclear = sprintf("nuclear (%s reads)", fmt_n(comp[compartment == "nuclear", n_total[1]])),
                                 chrM = sprintf("chrM (%s reads)", fmt_n(comp[compartment == "chrM", n_total[1]])))) +
  labs(x = "share of reads (log scale)", y = NULL,
       title = "A  End-type composition, nuclear vs mitochondrial",
       subtitle = sprintf("Pooled over %d CAGE libraries; right column is chrM / nuclear.", length(libs)))

pb_d <- copy(strata)[, stratum_lab := factor(ifelse(stratum >= ARUN_CAP, paste0("≥", ARUN_CAP), stratum),
                                              levels = c(1:(ARUN_CAP - 1), paste0("≥", ARUN_CAP)))]
p_b <- ggplot(pb_d, aes(stratum_lab, pct, colour = compartment, group = compartment)) +
  geom_line(linewidth = 0.7, alpha = 0.6) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 2, alpha = 0.35) +
  geom_point(size = 2.6) +
  geom_text(aes(label = paste0(fmt_n(sites), " sites")), size = 2.4, show.legend = FALSE,
            vjust = ifelse(pb_d$compartment == "chrM", -1.1, 2), colour = "grey30") +
  scale_y_log10(labels = function(x) paste0(x, "%"), expand = expansion(mult = 0.15)) +
  scale_colour_manual(values = comp_cols) +
  labs(x = "length of the A run starting at +1", y = paste(TYPE_A1, "share of reads (log scale)"),
       title = sprintf("B  %s at Y|A starts, stratified by A-run length", TYPE_A1),
       subtitle = "Dot: pooled share. Bar: range over libraries. Nuclear signal tracks A runs; chrM does not.") +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92"))

pc_d <- copy(top)[, site := factor(site, levels = rev(site))]
p_c <- ggplot(pc_d, aes(frac, site)) +
  geom_vline(xintercept = nuc_ya_pct, linetype = "dashed", colour = COL_NUC) +
  geom_segment(aes(x = 0.1, xend = frac, yend = site), colour = "grey80") +
  geom_point(aes(shape = promoter_core), colour = COL_MITO, size = 2.8, stroke = 1) +
  geom_text(aes(x = 150, label = sprintf("%s   %s", seq_lab, fmt_n(n_total))),
            hjust = 0, family = "mono", size = 2.5, colour = "grey25") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "mtRNAP promoter core AAGT|A", `FALSE` = "other context")) +
  scale_x_log10(labels = function(x) paste0(x, "%"), breaks = c(0.1, 1, 10, 100)) +
  coord_cartesian(xlim = c(0.1, 100), clip = "off") +
  labs(x = paste(TYPE_A1, "share of reads at the position (log scale)"), y = NULL,
       title = sprintf("C  chrM positions with ≥ 30 %s reads", TYPE_A1),
       subtitle = sprintf("Dashed: nuclear Y|A baseline (%.2f%%). Right: sequence -6..-1 | +1..+3 and total depth.",
                          nuc_ya_pct)) +
  theme(plot.margin = margin(5, 110, 5, 5))

pos_lab <- c(paste0("-", CTX_UP:1), paste0("+", 1:CTX_DN))
logo <- function(seqs, title) {
  ggseqlogo(seqs, method = "bits") +
    scale_x_continuous(breaks = seq_along(pos_lab), labels = pos_lab, expand = c(0, 0)) +
    labs(title = title, x = NULL) +
    theme(plot.title = element_text(size = 9, face = "plain"), axis.text.x = element_text(size = 7),
          legend.position = "none")
}
p_d <- logo(logo_mito, sprintf("chrM, positions with ≥ 20 %s reads (n = %d)", TYPE_A1, length(logo_mito))) /
       logo(logo_nuc,  sprintf("nuclear Y|A, %s > 20%% of ≥ 20 reads (n = %d)", TYPE_A1, length(logo_nuc))) +
  plot_annotation(title = "D  Sequence around +1 of single-A-clip-rich positions")
p_d <- wrap_elements(full = p_d)

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = sprintf("The soft-clipped-A end type (%s) is enriched on the uncapped mitochondrial genome", TYPE_A1),
    subtitle = sprintf(paste0("At Y|A starts, %s is %.2f%% of chrM reads against %.2f%% of nuclear reads (×%.0f). ",
                              "Candidate explanation: NAD caps, which yeast mtRNAP adds at its promoter starts."),
                       TYPE_A1, mito_ya_pct, nuc_ya_pct, mito_ya_pct / nuc_ya_pct),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9.5, colour = "grey30")))

ggsave(fig_pdf, fig, width = 15, height = 10.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", fig_pdf), fig, width = 15, height = 10.5, dpi = 200, bg = "white")
message("wrote ", fig_pdf)

fwrite(data.table(metric = c("nuclear_YA_pct", "chrM_YA_pct", "fold"),
                  value = c(nuc_ya_pct, mito_ya_pct, mito_ya_pct / nuc_ya_pct)),
       file.path(outdir, "headline.tsv"), sep = "\t")
sessionInfo()
