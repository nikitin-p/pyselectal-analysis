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
# Note on p-values: reads at one genomic position are not independent, so with millions of
# 1Sg reads nearly every cell clears q < 0.05 regardless of effect size. The asterisks certify
# direction, not importance — panel B carries the effect size. The same test computed on
# distinct positions instead of reads is written to the TSV as a sensitivity column.
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
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(yaml)
})

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
log2_cap  <- as.numeric(get_arg("--log2_cap", "3"))
force     <- has_flag("--force")

validate_args()

if (is.na(flank) || flank < 1L) stop("--flank must be a positive integer")
if (is.na(log2_cap) || log2_cap <= 0) stop("--log2_cap must be a positive number")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

out_main <- file.path(outdir, "dinuc_class_by_position.tsv")
out_bg   <- file.path(outdir, "genome_background.tsv")
out_fig  <- file.path(figdir, "14_dinuc_class_by_position.pdf")

if (!force && file.exists(out_main)) {
  message("Output exists, skipping (use --force to rerun): ", out_main)
  quit(save = "no", status = 0)
}

YR_RED   <- "#e63946"   # canonical initiator colour (project convention)
OTHER_GY <- "#6c757d"

# 'other' is excluded on purpose: it mixes clip lengths and clip first bases, so it has no
# single well-defined constraint at -1 and therefore no valid conditional expectation.
CLIP_BASE <- c("1Sg" = "G", "2Sgg" = "G", "3Sggg" = "G",
               "1Sa" = "A", "1Sc" = "C", "1St" = "T", "M" = NA_character_)
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
  b <- CLIP_BASE[[type]]
  if (is.na(b))    return(rep(TRUE, length(DINUCS)))   # M: unconstrained everywhere
  if (q ==  1L)    return(di_b1 != b)                  # dinuc (-1,+1): -1 is the first base
  if (q == -1L)    return(di_b2 != b)                  # dinuc (-2,-1): -1 is the second base
  rep(TRUE, length(DINUCS))
}

# ── load typed CTSS and combine samples ───────────────────────────────────────
message("Loading typed CTSS: ", ctss_f)
ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
message("  rows: ", nrow(ctss))

tcols <- paste0("n_", names(CLIP_BASE))
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
res[, obs_pct := 100 * obs_reads / n_reads]

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
res[, log2_ratio := log2((obs_pct / 100) / exp_prop)]

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
res[, sig := !is.na(q_bh) & q_bh < 0.05]
message("  tests: ", nrow(res), "; significant at BH 5%: ", sum(res$sig))

# ── write tables ──────────────────────────────────────────────────────────────
res[, end_type := factor(end_type, levels = TYPES)]
setorder(res, end_type, position, class)

out_cols <- c("end_type", "position", "class", "obs_reads", "n_reads", "obs_pct",
              "exp_pct", "log2_ratio", "p", "q_bh",
              "obs_positions", "n_positions", "p_positions", "q_bh_positions")
fwrite(res[, ..out_cols], out_main, sep = "\t")
message("Written: ", out_main)

bg <- data.table(dinucleotide = DINUCS, class = di_cls, genome_count = as.numeric(dn))
bg[, genome_freq := genome_count / sum(genome_count)]
for (b in setdiff(unique(CLIP_BASE), NA)) {
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

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid    = element_blank(),
    axis.text.y   = element_text(colour = ylab_cols, face = "bold"),
    axis.text.x   = element_text(size = 7),
    axis.title    = element_blank(),
    strip.text    = element_text(face = "bold", size = 9),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )

cell_label <- function(pct, sig) paste0(sprintf("%.1f", pct), ifelse(sig, "*", ""))

pA <- ggplot(res, aes(position_lab, class, fill = obs_pct)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label(obs_pct, sig)), size = 1.9) +
  facet_grid(. ~ end_type) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", name = "observed %") +
  labs(title = "A. Observed composition (each position column sums to 100 %)") +
  base_theme

# The scale is capped rather than stretched to the data range. A near-absent class produces
# log2 values around -9 (e.g. 1Sg at +1, RY), and letting those set the limits washes every
# informative cell to near-white. Values beyond the cap are squished onto the end colour.
obs_max <- max(abs(res$log2_ratio[is.finite(res$log2_ratio)]), na.rm = TRUE)
lim     <- min(log2_cap, obs_max)
message("  log2 range in data: +/-", round(obs_max, 2), "; colour scale capped at +/-", lim)

pB <- ggplot(res, aes(position_lab, class, fill = log2_ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label(obs_pct, sig)), size = 1.9) +
  facet_grid(. ~ end_type) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       oob = scales::squish,
                       name = paste0("log2(obs/exp)\ncapped at +/-", lim)) +
  labs(title = "B. Enrichment over the conditional genome expectation") +
  base_theme

cap <- paste0(
  "Position +1 is the first templated base; there is no zero position. The dinucleotide at ",
  "column q is (prev(q), q), so +1 is the (-1,+1) initiator. All ",
  length(unique(res$end_type)), " end types use all samples combined, read-weighted.\n",
  "Expected frequencies come from the ", basename(genome_f), " dinucleotide composition on ",
  "both strands. At columns -1 and +1 the expectation is conditioned on the STAR local-",
  "alignment constraint\n(the genomic base at -1 cannot equal the clip base), which is why ",
  "the +1 expectation differs between clip types and M. Asterisk: BH q < 0.05 across all ",
  nrow(res), " cells.\n",
  "Reads at one position are not independent, so with millions of reads the asterisks ",
  "certify direction rather than importance; panel B carries the effect size."
)

fig <- (pA / pB) +
  plot_annotation(caption = cap,
                  theme = theme(plot.caption = element_text(hjust = 0, size = 6.5,
                                                            colour = OTHER_GY)))

fig_w <- 3 + 2.4 * length(TYPES)
ggsave(out_fig, fig, width = fig_w, height = 7, limitsize = FALSE)
message("Written: ", out_fig)

sessionInfo()
