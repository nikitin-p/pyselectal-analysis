#!/usr/bin/env Rscript
# Stage 22 step 7: where the long-clip reads sit in gene structure, how long
# their clips are there, and whether the structure of the 5' end (including
# the untemplated "junk" sequence) depends on the location.
#
# Each read is located twice with the shared Stage 10 vocabulary
# (scripts/utils/localization.R, fast mode, cluster-defined promoters under
# promoter_source: clusters):
#   observed  the aligned 5' end (pos5), which is where a 5'-end count lands;
#   inferred  the 5' end once the clip is put back at its origin (tss5), which
#             equals pos5 for an unplaced read.
# sacCer3 carries no real UTRs (see below), so a gene body is CDS or Intron, and a position on
# the antisense strand of a gene, or past its 3' end, is Intergenic.
#
# The 5' structure of a read is one of:
#   G + exon 1          spliced, one untemplated G before the exon-1 piece;
#   G + upstream genome templated_adjacent, one untemplated G;
#   placed, other prefix  spliced or templated_adjacent with 0, 2 or 3
#                         untemplated bases, or one that is not G;
#   junk: low-complexity  unplaced (or templated_distal, which is at the
#                         chance level of the reversed-clip control) and at
#                         least LOWCOMP_FRAC of one base;
#   junk: G-initial, short / not found  the remaining unplaced G-initial clips,
#                         split by whether the clip was long enough to search;
#   junk: not G-initial   every other unplaced clip.
#
# Outputs (in --outdir):
#   lc_location_summary.tsv     reads per origin class x coordinate x location
#   lc_location_cliplen.tsv     clip length n / mean / sd / median per location,
#                               overall and per origin class (observed location)
#   lc_location_structure.tsv   reads per location x 5' structure, with the
#                               chi-square of structure against location
#   lc_location_composition.tsv base composition of clip positions 1..8 by
#                               location
#   lc_location_junk_recurrence.tsv  per unplaced clip sequence: reads, distinct
#                               5' end sites, samples, locations
#   <figdir>/28_longclip_location.pdf  panels A-E
#
# All shares are of reads (n_reads), not of distinct alignments.
#
# Usage (from project root):
#   Rscript scripts/14_longclip/07_location.R \
#       [--reads results/longclip/longclip_reads.tsv.gz] \
#       [--gtf data/sacCer3.ensGene.gtf] \
#       [--promoters results/tss/promoter_clusters.tsv] \
#       [--outdir results/longclip] [--figdir results/figures] \
#       [--params config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
reads_f  <- get_arg("--reads",     "results/longclip/longclip_reads.tsv.gz")
gtf_f    <- get_arg("--gtf",       "data/sacCer3.ensGene.gtf")
prom_f   <- get_arg("--promoters", "results/tss/promoter_clusters.tsv")
outdir   <- get_arg("--outdir",    "results/longclip")
figdir   <- get_arg("--figdir",    "results/figures")
params_f <- get_arg("--params",    "config/params.yaml")
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "lc_helpers.R"))
source(file.path(script_dir, "..", "utils", "localization.R"))
source(file.path(script_dir, "..", "utils", "txdb.R"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
fig_f <- file.path(figdir, "28_longclip_location.pdf")
if (!force && file.exists(fig_f)) {
  message("Output exists, skipping (use --force to rerun): ", fig_f)
  quit(save = "no", status = 0)
}

LOWCOMP_FRAC <- 0.8   # a clip this rich in one base counts as low-complexity
COMP_POS     <- 8L    # clip positions shown in the composition panel
set.seed(1L)

params <- yaml.load_file(params_f)
set_promoter_window(params)
prom_cl <- load_promoter_clusters(params, prom_f)

# ── annotation layers (as Stage 10 fast mode) ─────────────────────────────────
txdb <- suppressWarnings(make_txdb_from_gff(gtf_f, format = "gtf"))
tx <- transcripts(txdb)
promoter_gr <- if (is.null(prom_cl))
  promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN) else prom_cl
# The sacCer3 GTF yields only four 5' UTRs, each two bases long, where an exon
# starts two bases before its CDS. Kept, they turn one 2,066-read hotspot into
# a "5' UTR" category of its own, so the UTR layers are left empty here.
utr5_gr   <- GRanges()
utr3_gr   <- GRanges()
cds_gr    <- unlist(cdsBy(txdb, by = "tx"))
intron_gr <- unlist(intronsByTranscript(txdb))

reads <- fread(reads_f)
stopifnot(all(reads$class %in% LC_CLASSES))
locate <- function(pos) annotate_positions_fast(
  GRanges(reads$chr, IRanges(pos, width = 1L), strand = reads$strand),
  promoter_gr, utr5_gr, utr3_gr, cds_gr, intron_gr)
reads[, loc_observed := locate(pos5)]
reads[, loc_inferred := locate(fifelse(is.na(tss5), pos5, tss5))]
LOCS <- intersect(LOC_CATEGORY_ORDER, unique(c(reads$loc_observed, reads$loc_inferred)))
LOC_COLS <- LOC_CATEGORY_COLORS[LOCS]

# ── 5' structure ──────────────────────────────────────────────────────────────
STRUCT <- c("G + exon 1", "G + upstream genome", "placed, other prefix",
            "junk: low-complexity", "junk: G-initial, short",
            "junk: G-initial, not found", "junk: not G-initial")
STRUCT_COLS <- setNames(c("#1d4ed8", "#059669", "#93c5fd",
                          "#7f1d1d", "#f87171", "#dc2626", "#fbbf24"), STRUCT)

max_base_frac <- function(s) vapply(strsplit(s, ""), function(x)
  max(tabulate(match(x, c("A", "C", "G", "T")), 4L)) / length(x), numeric(1))
reads[, prefix := fifelse(!is.na(n_untemplated) & n_untemplated > 0L,
                          substr(clip, 1L, pmax(n_untemplated, 0L)), "")]
reads[, lowcomp := max_base_frac(clip) >= LOWCOMP_FRAC]
reads[, placed := class %in% c("splice_annotated", "splice_unannotated", "templated_adjacent")]
reads[, structure := fcase(
  placed & grepl("^splice", class) & prefix == "G", STRUCT[1],
  placed & prefix == "G",                            STRUCT[2],
  placed,                                            STRUCT[3],
  lowcomp,                                           STRUCT[4],
  substr(clip, 1, 1) == "G" & unplaced_reason == "short", STRUCT[5],
  substr(clip, 1, 1) == "G",                         STRUCT[6],
  default = STRUCT[7])]
reads[, structure := factor(structure, levels = STRUCT)]

# ── tables ────────────────────────────────────────────────────────────────────
long <- rbind(reads[, .(class, coordinate = "observed", location = loc_observed, n_reads)],
              reads[, .(class, coordinate = "inferred", location = loc_inferred, n_reads)])
loc_sum <- long[, .(n_reads = sum(n_reads)), by = .(class, coordinate, location)]
loc_sum <- rbind(loc_sum,
                 long[, .(class = "all", n_reads = sum(n_reads)), by = .(coordinate, location)])
loc_sum[, pct_reads := 100 * n_reads / sum(n_reads), by = .(class, coordinate)]
setorder(loc_sum, class, coordinate, location)
fwrite(loc_sum, file.path(outdir, "lc_location_summary.tsv"), sep = "\t")

wstats <- function(x, w) {
  m <- sum(w * x) / sum(w)
  list(n_reads = sum(w), mean = m, sd = sqrt(sum(w * (x - m)^2) / (sum(w) - 1)),
       median = as.numeric(median(rep(x, w))))
}
len <- rbind(
  reads[, c(list(class = "all", location = "all"), wstats(clen, n_reads))],
  reads[, c(list(class = "all"), wstats(clen, n_reads)), by = .(location = loc_observed)],
  reads[, c(list(location = "all"), wstats(clen, n_reads)), by = class],
  reads[, wstats(clen, n_reads), by = .(class, location = loc_observed)], use.names = TRUE)
len[, sd := fifelse(n_reads > 1, sd, NA_real_)]
fwrite(len, file.path(outdir, "lc_location_cliplen.tsv"), sep = "\t")

st <- reads[, .(n_reads = sum(n_reads)), by = .(location = loc_observed, structure)]
st[, pct_reads := 100 * n_reads / sum(n_reads), by = location]
ct <- as.matrix(dcast(st, location ~ structure, value.var = "n_reads", fill = 0)[, -1])
chi <- suppressWarnings(chisq.test(ct))
cramer <- sqrt(unname(chi$statistic) / (sum(ct) * (min(dim(ct)) - 1)))
st[, `:=`(chisq = unname(chi$statistic), df = unname(chi$parameter),
          p_value = chi$p.value, cramers_v = cramer)]
setorder(st, location, structure)
fwrite(st, file.path(outdir, "lc_location_structure.tsv"), sep = "\t")

comp <- reads[, .(pos = seq_len(min(clen, COMP_POS)),
                  base = strsplit(substr(clip, 1, COMP_POS), "")[[1]], n_reads),
              by = .(row = seq_len(nrow(reads)), location = loc_observed)]
comp <- comp[, .(n_reads = sum(n_reads)), by = .(location, pos, base)]
comp[, frac := n_reads / sum(n_reads), by = .(location, pos)]
fwrite(comp, file.path(outdir, "lc_location_composition.tsv"), sep = "\t")

junk <- reads[!(placed), .(n_reads = sum(n_reads),
                           n_sites = uniqueN(paste(chr, strand, pos5)),
                           max_samples = max(n_samples), lowcomp = lowcomp[1],
                           locations = paste(sort(unique(loc_observed)), collapse = ",")),
              by = clip][order(-n_reads)]
fwrite(junk, file.path(outdir, "lc_location_junk_recurrence.tsv"), sep = "\t")

# ── figure ────────────────────────────────────────────────────────────────────
th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.key.size = unit(0.35, "cm"))
CLASS_LEVELS <- c("all", LC_CLASSES)
n_lab <- reads[, .(n = sum(n_reads)), by = class]
n_lab <- rbind(data.table(class = "all", n = sum(reads$n_reads)), n_lab)
class_lab <- setNames(sprintf("%s\n(n = %s)", n_lab$class, format(n_lab$n, big.mark = ",")),
                      n_lab$class)[CLASS_LEVELS]

# A: location of the 5' end, observed and inferred, per origin class
pa_d <- copy(loc_sum)[, `:=`(class = factor(class, CLASS_LEVELS, labels = class_lab),
                             location = factor(location, rev(LOCS)),
                             coordinate = factor(coordinate, c("observed", "inferred"),
                               c("observed (aligned 5' end)", "inferred (clip at its origin)")))]
pA <- ggplot(pa_d, aes(class, pct_reads, fill = location)) +
  geom_col(width = 0.8) + facet_wrap(~coordinate) +
  scale_fill_manual(values = LOC_COLS, breaks = LOCS, name = "location") +
  labs(x = NULL, y = "% of reads",
       title = "A  Where long-clip 5' ends lie",
       subtitle = "spliced reads align inside introns but start in promoters once the clip is restored") +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

# B: clip length per observed location, with mean +- SD
pb_len <- reads[, .(n_reads = sum(n_reads)), by = .(location = loc_observed, clen)]
pb_len[, frac := n_reads / sum(n_reads), by = location]
pb_len[, location := factor(location, LOCS)]
pb_st <- len[class == "all" & location != "all"][, location := factor(location, LOCS)]
pb_st[, lab := sprintf("%.1f ± %.1f nt\n(n = %s)", mean, sd, format(n_reads, big.mark = ","))]
ov <- len[class == "all" & location == "all"]
pB <- ggplot(pb_len, aes(clen, frac, fill = location)) +
  geom_col(width = 0.85) +
  geom_vline(data = pb_st, aes(xintercept = mean), linetype = 2, linewidth = 0.4) +
  geom_errorbarh(data = pb_st, aes(y = max(pb_len$frac) * 1.05, xmin = mean - sd, xmax = mean + sd),
                 height = 0, inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = pb_st, aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.2,
            size = 2.6, inherit.aes = FALSE) +
  facet_wrap(~location, nrow = 1) +
  scale_fill_manual(values = LOC_COLS, guide = "none") +
  labs(x = "clip length (nt)", y = "fraction of reads",
       title = "B  Clip length by location (mean ± SD, read-weighted)",
       subtitle = sprintf("all long clips: %.2f ± %.2f nt, median %g (clips > %s nt)",
                          ov$mean, ov$sd, ov$median, params$dom_max_clip)) + th

# C: 5' structure by observed location
pc_d <- copy(st)[, location := factor(location, LOCS)]
pc_n <- st[, .(n = sum(n_reads)), by = location][, location := factor(location, LOCS)]
pC <- ggplot(pc_d, aes(location, pct_reads, fill = structure)) +
  geom_col(width = 0.75) +
  geom_text(data = pc_n, aes(location, 102, label = format(n, big.mark = ",")),
            size = 2.5, inherit.aes = FALSE) +
  scale_fill_manual(values = STRUCT_COLS, name = "5' structure") +
  labs(x = "observed location", y = "% of reads",
       title = "C  5' structure by location",
       subtitle = sprintf("chi-square = %.0f, df = %d, Cramér's V = %.2f",
                          unname(chi$statistic), unname(chi$parameter), cramer)) + th

# D: base composition of the first clip positions by location
pd_d <- copy(comp)[, `:=`(location = factor(location, LOCS),
                          base = factor(base, c("A", "C", "G", "T")))]
pD <- ggplot(pd_d, aes(pos, frac, fill = base)) +
  geom_col(width = 0.9) + facet_wrap(~location, nrow = 1) +
  scale_fill_manual(values = c(A = "#16a34a", C = "#2563eb", G = "#f59e0b", T = "#dc2626")) +
  scale_x_continuous(breaks = seq_len(COMP_POS)) +
  labs(x = "clip position from the read 5' end", y = "base fraction",
       title = "D  Clip composition by location") + th

# E: recurrence of the unplaced clip sequences
lab_d <- junk[seq_len(min(10L, .N))]
pE <- ggplot(junk, aes(n_sites, n_reads, colour = lowcomp)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_text_repel(data = lab_d, aes(label = clip), size = 2.2, show.legend = FALSE,
                  min.segment.length = 0, seed = 1L, max.overlaps = Inf) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = "#7f1d1d", `FALSE` = "#9ca3af"),
                      labels = c(`TRUE` = "low-complexity", `FALSE` = "other"), name = NULL) +
  labs(x = "distinct 5' end sites", y = "reads",
       title = "E  Recurrence of unplaced clips",
       subtitle = "few sites: locus-specific; many sites: low-complexity junk") + th

fig <-pA / pB / (pC | pE) / pD + plot_layout(heights = c(1.1, 0.8, 1, 0.7))
ggsave(fig_f, fig, width = 15, height = 17)
message("Written: ", fig_f)
sessionInfo()
