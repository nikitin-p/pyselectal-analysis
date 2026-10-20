#!/usr/bin/env Rscript
# Stage 8 step 9: dinucleotide class (YR / YY / RR / RY) by position, across 5' end types
#
# Motivation
# ----------
# The canonical initiator dinucleotides (CA, CG, TA, TG — the YR class) should be a property
# of the transcription start site itself, not of the sequence around it. This script tests
# that directly: for every 5' end type it tabulates the YR/YY/RR/RY composition of the
# dinucleotide ENDING at each position in a window around the CTSS, with all samples
# combined, and compares each cell against the genomic expectation.
#
# The dinucleotide AT position q is the pair (prev(q), q), where prev(+1) = -1. The +1 column
# is therefore the (-1,+1) initiator and the -4 column is (-5,-4). The eight output columns
# are the eight consecutive dinucleotides of one nine-base window spanning -5..+4.
#
# The alignment constraint at -1
# ------------------------------
# STAR ran with --alignEndsType Local. A read is only classified 1Sg when the genomic base at
# -1 is not G, since otherwise the aligner would have extended the match and the read would
# be M. The same holds for 2Sgg and 3Sggg (a G at -1 turns a 2Sgg read into a 1Sg read), and
# for 1Sa / 1Sc / 1St with their own clip base. So for a G-clip type the -1 base is drawn from
# {A,C,T}, which puts Y = {C,T} at about two thirds of the -1 pool by construction rather than
# about one half — inflating YR at exactly the cell under test.
#
# The genome expectation is therefore CONDITIONED for the two columns whose dinucleotide
# touches -1: at column +1 the tally is restricted to dinucleotides whose first base differs
# from the clip base, and at column -1 to those whose second base differs from it. M reads
# carry no soft clip and so impose no constraint anywhere.
#
# Note on the unit of analysis: reads at one genomic position are not independent, so with
# millions of 1Sg reads the read-based test clears any BH cutoff regardless of effect size —
# most of its p-values underflow to zero outright. The figure is therefore built entirely on
# distinct TSS positions, each counted once: the fill is log2_ratio_pos and the asterisk is
# q_bh_positions < Q_CUT. The read-weighted obs_pct, log2_ratio, p and q_bh stay in the TSV
# as a sensitivity view, but they must not drive fill and asterisk separately — the two
# ratios disagree by up to 4 log2 units, so mixing them would star near-white cells.
#
# Usage:
#   Rscript scripts/05_tss_end_types/09_dinuc_class_by_position.R \
#       [--ctss    results/tss/typed_ctss.tsv.gz] \
#       [--params  config/params.yaml]            \
#       [--samples config/samples.tsv]            \
#       [--genome  data/sacCer3.fa]               \
#       [--outdir  results/dinuc_position]        \
#       [--figdir  results/figures]               \
#       [--flank   4]                             \
#       [--types   all|coarse]                    \
#       [--force]
#
# --types coarse restricts the figure to END_TYPES_COARSE, the end-type set drawn in
# results/figures/03_replicate_heatmap.pdf, and writes its own outputs (14b_*, *_coarse.tsv)
# so that the full-vocabulary figure is left in place beside it.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(yaml)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

.val_flags  <- character(0)   # flags that consume a following value
.bool_flags <- character(0)   # standalone switches

get_arg <- function(flag, default = NULL) {
  .val_flags <<- union(.val_flags, flag)
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

has_flag <- function(flag) {
  .bool_flags <<- union(.bool_flags, flag)
  flag %in% args
}

validate_args <- function() {
  unknown <- character(0)
  i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% .val_flags)  { i <- i + 2L; next }
    if (args[i] %in% .bool_flags) { i <- i + 1L; next }
    unknown <- c(unknown, args[i])
    i <- i + 1L
  }
  if (length(unknown) > 0) {
    stop("Unrecognized argument(s): ", paste(unknown, collapse = " "),
         "\n  Known flags: ", paste(sort(c(.val_flags, .bool_flags)), collapse = " "),
         call. = FALSE)
  }
}

ctss_f    <- get_arg("--ctss",    "results/tss/typed_ctss.tsv.gz")
params_f  <- get_arg("--params",  "config/params.yaml")
samples_f <- get_arg("--samples", "config/samples.tsv")
genome_f  <- get_arg("--genome",  NULL)
outdir    <- get_arg("--outdir",  "results/dinuc_position")
figdir    <- get_arg("--figdir",  "results/figures")
flank     <- as.integer(get_arg("--flank", "4"))
type_set  <- get_arg("--types",   "all")
log2_cap  <- as.numeric(get_arg("--log2_cap", "3"))
force     <- has_flag("--force")

validate_args()

if (!type_set %in% c("all", "coarse"))
  stop("--types must be one of: all, coarse")
if (is.na(flank) || flank < 1L) stop("--flank must be a positive integer")
if (is.na(log2_cap) || log2_cap <= 0) stop("--log2_cap must be a positive number")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# The restricted run is a companion to the full one, not a replacement, so it takes its own
# output stems and its own figure number rather than overwriting the full-vocabulary files.
tag        <- if (type_set == "all") "" else paste0("_", type_set)
fig_number <- if (type_set == "all") "14" else "14b"

out_main <- file.path(outdir, paste0("dinuc_class_by_position", tag, ".tsv"))
out_bg   <- file.path(outdir, paste0("genome_background", tag, ".tsv"))
out_fig  <- file.path(figdir, paste0(fig_number, "_dinuc_class_by_position", tag, ".pdf"))

if (!force && file.exists(out_main)) {
  message("Output exists, skipping (use --force to rerun): ", out_main)
  quit(save = "no", status = 0)
}

YR_RED   <- "#e63946"   # canonical initiator colour (project convention)
OTHER_GY <- "#6c757d"

# 'other' is excluded on purpose: it mixes clip lengths and clip first bases, so it has no
# single well-defined constraint at -1 and therefore no valid conditional expectation.
# END_TYPE_CLIP_BASE (scripts/utils/end_types.R) is keyed over END_TYPES_REAL, which
# includes "other"'s siblings but not "other" itself, so no further filtering is needed here.
CLASSES   <- c("YR", "YY", "RR", "RY")

# ── genome resolution and loading ─────────────────────────────────────────────
if (is.null(genome_f)) {
  meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                   show_col_types = FALSE)
  if (!"genome" %in% colnames(meta)) stop("No 'genome' column in ", samples_f)
  g <- unique(meta$genome[!is.na(meta$genome)])
  if (length(g) != 1) {
    stop("Expected exactly one genome in ", samples_f, ", found: ",
         paste(g, collapse = ", "),
         "\n  One combined analysis cannot span two genomes; pass --genome explicitly.")
  }
  genome_f <- g
}
message("Genome: ", genome_f)

load_genome <- function(spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", spec, ignore.case = TRUE)) {
    if (!file.exists(spec)) stop("FASTA not found: ", spec)
    if (!file.exists(paste0(spec, ".fai"))) Rsamtools::indexFa(spec)
    fa <- Rsamtools::FaFile(spec)
    open(fa)
    return(fa)
  }
  pkg <- switch(spec,
    hg38    = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias: ", spec)
  )
  if (!requireNamespace(pkg, quietly = TRUE)) stop("BSgenome not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

genome <- load_genome(genome_f)

# ── genome-wide dinucleotide background, both strands ─────────────────────────
# Reads arrive on both strands, so the null must count both. A reverse-strand dinucleotide
# corresponds exactly to its reverse complement on the forward strand, so the reverse tally
# is obtained by relabelling rather than by reverse-complementing the whole genome.
message("Tallying genome dinucleotides (both strands)...")
gseqs <- if (inherits(genome, "FaFile")) scanFa(genome) else getSeq(genome)
dn_fwd <- colSums(dinucleotideFrequency(gseqs))
rc_map <- as.character(reverseComplement(DNAStringSet(names(dn_fwd))))
dn     <- dn_fwd + dn_fwd[rc_map]
rm(gseqs)

DINUCS <- names(dn)
di_b1  <- substr(DINUCS, 1L, 1L)
di_b2  <- substr(DINUCS, 2L, 2L)
di_cls <- paste0(ifelse(di_b1 %in% c("C", "T"), "Y", "R"),
                 ifelse(di_b2 %in% c("C", "T"), "Y", "R"))
names(di_cls) <- DINUCS

# keep_mask -> expected proportion of each class within that restricted dinucleotide pool
exp_props <- function(keep) {
  tot <- sum(dn[keep])
  if (tot == 0) return(setNames(rep(NA_real_, length(CLASSES)), CLASSES))
  vapply(CLASSES, function(cl) sum(dn[keep & di_cls == cl]) / tot, numeric(1))
}

# The restriction that applies to end type `type` at output position `q`.
bg_mask <- function(type, q) {
  b <- END_TYPE_CLIP_BASE[[type]]
  if (is.na(b))    return(rep(TRUE, length(DINUCS)))   # M: unconstrained everywhere
  if (q ==  1L)    return(di_b1 != b)                  # dinuc (-1,+1): -1 is the first base
  if (q == -1L)    return(di_b2 != b)                  # dinuc (-2,-1): -1 is the second base
  rep(TRUE, length(DINUCS))
}

# ── load typed CTSS and combine samples ───────────────────────────────────────
message("Loading typed CTSS: ", ctss_f)
ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
message("  rows: ", nrow(ctss))

sel_types <- names(END_TYPE_CLIP_BASE)
if (type_set == "coarse") {
  sel_types <- intersect(END_TYPES_COARSE, sel_types)
  message("Restricted to the Stage 3 end-type set: ", paste(sel_types, collapse = ", "))
}
tcols <- paste0("n_", sel_types)
missing_cols <- setdiff(tcols, colnames(ctss))
if (length(missing_cols) > 0) {
  stop("Missing end-type columns in ", ctss_f, ": ", paste(missing_cols, collapse = ", "))
}

message("Combining samples...")
agg <- ctss[, lapply(.SD, sum), by = .(chr, pos, strand), .SDcols = tcols]
rm(ctss)
message("  unique positions: ", nrow(agg))

# Facet order is data-driven: most abundant end type first.
totals <- vapply(tcols, function(cc) sum(agg[[cc]]), numeric(1))
TYPES  <- sub("^n_", "", names(sort(totals, decreasing = TRUE)))
message("  end types by abundance: ", paste(TYPES, collapse = ", "))

# ── window extraction ─────────────────────────────────────────────────────────
# +1 is the first templated base; -1 is immediately upstream; there is no zero position.
# The window spans -(flank+1)..-1 and +1..+flank, so its width is 2*flank + 1.
w <- 2L * flank + 1L
is_plus <- agg$strand != "-"
start_v <- ifelse(is_plus, agg$pos - (flank + 1L), agg$pos - flank + 1L)

gr <- GRanges(seqnames = agg$chr,
              ranges   = IRanges(start = start_v, width = w),
              strand   = ifelse(is_plus, "+", "-"))

si     <- seqinfo(genome)
common <- intersect(seqlevels(gr), seqlevels(si))
if (length(common) == 0) stop("No common seqlevels between CTSS and genome")

keep <- as.character(seqnames(gr)) %in% common
if (!all(keep)) {
  message("  dropping ", sum(!keep), " positions on seqlevels absent from the genome")
  gr <- gr[keep]; agg <- agg[keep]
}
gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
seqinfo(gr) <- si[seqlevels(gr)]

# Out-of-bound windows at chromosome edges are dropped rather than trimmed, so that every
# column keeps the same denominator and the four class rows sum to exactly 100 %.
lens <- seqlengths(gr)[as.character(seqnames(gr))]
inb  <- start(gr) >= 1L & end(gr) <= lens
if (!all(inb)) {
  message("  dropping ", sum(!inb), " positions whose window runs off a chromosome end")
  gr <- gr[inb]; agg <- agg[inb]
}

message("Extracting ", w, "-base windows for ", length(gr), " positions...")
s <- as.character(getSeq(genome, gr))   # getSeq is strand-aware and returns 5'->3'

has_n <- grepl("[^ACGT]", s)
if (any(has_n)) {
  message("  dropping ", sum(has_n), " windows containing non-ACGT bases")
  s <- s[!has_n]; agg <- agg[!has_n]
}
if (length(s) == 0) stop("No usable windows remain after filtering")

# Position labels: the window covers pos_full, and the dinucleotide at column j ends at
# pos_full[j + 1].
pos_full <- c(seq.int(-(flank + 1L), -1L), seq.int(1L, flank))
ends     <- pos_full[-1]
stopifnot(length(pos_full) == w, length(ends) == 2L * flank)

# ── observed counts per (end type, position, class) ───────────────────────────
message("Counting dinucleotide classes at ", length(ends), " positions...")
obs_list <- vector("list", length(ends))

for (j in seq_along(ends)) {
  cl  <- factor(di_cls[substr(s, j, j + 1L)], levels = CLASSES)
  tmp <- cbind(data.table(cl = cl), agg[, ..tcols])

  reads <- tmp[, lapply(.SD, sum),              by = cl, .SDcols = tcols]
  npos  <- tmp[, lapply(.SD, function(x) sum(x > 0L)), by = cl, .SDcols = tcols]

  reads <- melt(reads, id.vars = "cl", variable.name = "col", value.name = "obs_reads")
  npos  <- melt(npos,  id.vars = "cl", variable.name = "col", value.name = "obs_positions")

  d <- merge(reads, npos, by = c("cl", "col"))
  d[, `:=`(end_type = sub("^n_", "", as.character(col)),
           position = ends[j],
           col      = NULL)]
  setnames(d, "cl", "class")
  obs_list[[j]] <- d
}

res <- rbindlist(obs_list)
res[is.na(obs_reads),     obs_reads     := 0]
res[is.na(obs_positions), obs_positions := 0]

# Denominators: constant across positions by construction, but computed per column anyway so
# that any future change to the filtering cannot silently break the 100 % property.
res[, `:=`(n_reads     = sum(obs_reads),
           n_positions = sum(obs_positions)), by = .(end_type, position)]
res[, obs_pct     := 100 * obs_reads / n_reads]
res[, obs_pct_pos := 100 * obs_positions / n_positions]

# ── expected proportions ──────────────────────────────────────────────────────
message("Computing conditional genome expectations...")
exp_tbl <- rbindlist(lapply(TYPES, function(t) {
  rbindlist(lapply(ends, function(q) {
    e <- exp_props(bg_mask(t, q))
    data.table(end_type = t, position = q, class = CLASSES, exp_prop = as.numeric(e))
  }))
}))

res <- merge(res, exp_tbl, by = c("end_type", "position", "class"))
res[, exp_pct    := 100 * exp_prop]
res[, log2_ratio     := log2((obs_pct     / 100) / exp_prop)]
res[, log2_ratio_pos := log2((obs_pct_pos / 100) / exp_prop)]

# ── proportion tests ──────────────────────────────────────────────────────────
message("Running proportion tests...")
prop_p <- function(x, n, p0) {
  if (is.na(p0) || n <= 0 || p0 <= 0 || p0 >= 1) return(NA_real_)
  suppressWarnings(stats::prop.test(x = x, n = n, p = p0)$p.value)
}

res[, p           := mapply(prop_p, obs_reads,     n_reads,     exp_prop)]
res[, p_positions := mapply(prop_p, obs_positions, n_positions, exp_prop)]

# One BH family across the whole figure: 4 classes x 2*flank positions x n end types.
res[, q_bh           := p.adjust(p,           method = "BH")]
res[, q_bh_positions := p.adjust(p_positions, method = "BH")]
Q_CUT <- 0.001
# Significance is judged on distinct positions, not reads: reads within a position are not
# independent, and the read-based p-values underflow to zero for most cells.
res[, sig := !is.na(q_bh_positions) & q_bh_positions < Q_CUT]
message("  tests: ", nrow(res), "; significant at BH q_positions < ", Q_CUT, ": ", sum(res$sig),
        " (read-based q would give ", sum(!is.na(res$q_bh) & res$q_bh < Q_CUT), ")")

# ── write tables ──────────────────────────────────────────────────────────────
res[, end_type := factor(end_type, levels = TYPES)]
setorder(res, end_type, position, class)

out_cols <- c("end_type", "position", "class", "exp_pct",
              "obs_positions", "n_positions", "obs_pct_pos", "log2_ratio_pos",
              "p_positions", "q_bh_positions",
              "obs_reads", "n_reads", "obs_pct", "log2_ratio", "p", "q_bh")
fwrite(res[, ..out_cols], out_main, sep = "\t")
message("Written: ", out_main)

bg <- data.table(dinucleotide = DINUCS, class = di_cls, genome_count = as.numeric(dn))
bg[, genome_freq := genome_count / sum(genome_count)]
for (b in setdiff(unique(END_TYPE_CLIP_BASE), NA)) {
  bg[[paste0("used_at_plus1_clip", b)]]  <- di_b1 != b
  bg[[paste0("used_at_minus1_clip", b)]] <- di_b2 != b
}
fwrite(bg, out_bg, sep = "\t")
message("Written: ", out_bg)

# ── figure ────────────────────────────────────────────────────────────────────
message("Drawing figure...")
res[, position_lab := factor(sprintf("%+d", position), levels = sprintf("%+d", ends))]
# Levels run bottom to top, so YR is listed last to place it on the top row.
res[, class := factor(as.character(class), levels = c("RY", "RR", "YY", "YR"))]

ylab_cols <- ifelse(levels(res$class) == "YR", YR_RED, OTHER_GY)

# The initiator column is the point of the figure, so it is marked twice: its axis label is
# drawn in the canonical-initiator colour and a box is stroked around the column itself.
x_levels  <- levels(res$position_lab)
plus1_i   <- which(x_levels == "+1")
xlab_cols <- ifelse(x_levels == "+1", YR_RED, OTHER_GY)
xlab_face <- ifelse(x_levels == "+1", "bold", "plain")

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid    = element_blank(),
    axis.text.y   = element_text(colour = ylab_cols, face = "bold"),
    axis.text.x   = element_text(size = 7, colour = xlab_cols, face = xlab_face),
    axis.title    = element_blank(),
    strip.text    = element_text(face = "bold", size = 9),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )

# The scale is capped rather than stretched to the data range. A near-absent class produces
# log2 values around -9 (e.g. 1Sg at +1, RY), and letting those set the limits washes every
# informative cell to near-white. Values beyond the cap are squished onto the end colour.
obs_max <- max(abs(res$log2_ratio_pos[is.finite(res$log2_ratio_pos)]), na.rm = TRUE)
lim     <- min(log2_cap, obs_max)
message("  log2 range in data: +/-", round(obs_max, 2), "; colour scale capped at +/-", lim)

n_cls <- nlevels(res$class)

# The end types wrap into a near-square grid rather than one row: sixteen panels of eight
# columns side by side would make a canvas about forty inches wide.
n_types   <- length(unique(res$end_type))
facet_col <- ceiling(sqrt(n_types))
facet_row <- ceiling(n_types / facet_col)

fig <- ggplot(res, aes(position_lab, class, fill = log2_ratio_pos)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 3.6, vjust = 0.75) +
  annotate("rect", xmin = plus1_i - 0.5, xmax = plus1_i + 0.5,
           ymin = 0.5, ymax = n_cls + 0.5,
           fill = NA, colour = YR_RED, linewidth = 0.7) +
  facet_wrap(~ end_type, ncol = facet_col) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       oob = scales::squish,
                       name = paste0("log2(obs/exp)\nper TSS position\ncapped at +/-", lim)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  labs(title = paste0("Enrichment over the conditional genome expectation",
                      if (type_set == "all") ""
                      else " - end-type set of 03_replicate_heatmap")) +
  base_theme

cap <- paste0(
  "Position +1 is the first templated base; there is no zero position. The dinucleotide at ",
  "column q is (prev(q), q), so the boxed +1 column is the (-1,+1) initiator. All ",
  length(unique(res$end_type)), " end types use all samples combined,\nalignment-weighted. ",
  if (type_set == "all") "" else paste0(
    "The end types shown are restricted to the set named in 03_replicate_heatmap.pdf, which ",
    "resolves the G-appending ladder, its two two-base siblings and M and sweeps every ",
    "single-base non-G clip into a residual category;\nthat residual is not plotted here ",
    "because it mixes clip bases and so has no single conditional expectation. "),
  "Expected frequencies come from the ", basename(genome_f), " dinucleotide composition on ",
  "both strands. At columns -1 and +1 the expectation is conditioned on the STAR local-",
  "alignment constraint\n(the genomic base at -1 cannot equal the clip base), which is why ",
  "the +1 expectation differs between clip types and M. Fill and asterisk both count each ",
  "distinct TSS position once rather than each read;\nasterisk: BH q < ", Q_CUT,
  " across all ", nrow(res), " cells. Reads at one position are not independent, so a ",
  "read-weighted test calls every cell significant regardless of effect size. The ",
  "read-weighted percentages and q-values are\nkept in ", basename(out_main),
  " as a sensitivity view, alongside the per-position ones plotted here."
)

fig <- fig +
  labs(caption = cap) +
  theme(plot.caption = element_text(hjust = 0, size = 6.5, colour = OTHER_GY))

# coord_fixed() makes the cells square, so the canvas is sized from a target cell edge rather
# than the other way round: whichever of width or height binds first would otherwise be padded
# out with whitespace.
CELL_IN <- 0.30
STRIP_IN <- 0.45                          # facet strip plus the x-axis labels of each row
fig_w <- 1.9 + CELL_IN * facet_col * length(ends) + 0.12 * facet_col
fig_h <- 1.9 + (CELL_IN * n_cls + STRIP_IN) * facet_row
# The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
for (f in c(out_fig, sub("\\.pdf$", ".png", out_fig)))
  ggsave(f, fig, width = fig_w, height = fig_h, limitsize = FALSE, dpi = 200, bg = "white")
message("Written: ", out_fig)

sessionInfo()
