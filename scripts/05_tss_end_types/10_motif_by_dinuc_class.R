#!/usr/bin/env Rscript
# Stage 8 step 10: sequence motifs per 5' end type, split by initiator dinucleotide class
#
# Motivation
# ----------
# A single logo per 5' end type mixes start sites that differ in the one feature the whole
# project is about: the initiator dinucleotide at (-1, +1). Because the YR class alone accounts
# for 95 % of 1Sg reads and for well under half of the M reads, a plain per-type comparison of
# logos partly reports a difference in class composition rather than a difference in sequence
# context. This script therefore emits five motifs per end type: one over all reads of that
# type, and one for each of YR, YY, RR and RY. Reading a row across gives the context of that
# end type at a canonical start versus a non-canonical one; reading a column down gives the
# end types compared at matched initiator class.
#
# Which end types are included
# ----------------------------
# The default set (--types) omits 3Sggg. Splitting its 16,369 reads four ways leaves the
# non-YR classes with too few reads for a readable logo. Every other type carries at least
# 44,000 reads. See results/dinuc_position/dinuc_class_by_position.tsv for the per-type depth
# that this choice is based on.
#
# What one site is: --level
# -------------------------
# The default, --level ctss, treats every single-base 5' position as its own site, so a promoter
# that fires from a spread of nearby bases contributes one window per base. This is the level at
# which the end-type question is actually posed: the clip is a property of a start base, not of
# a promoter.
#
# --level tss instead treats a distclu tag cluster from results/tss/tss_matrix.tsv as one site,
# anchoring the window on the cluster's dominant position and weighting it by the whole
# cluster's typed counts. Because the dominant position of a cluster is overwhelmingly a
# canonical initiator, that level collapses the non-YR panels to a few hundred sites and makes
# most of the grid unreadable; it answers "what does a promoter look like" rather than "what
# does a start base look like". Its outputs carry a _tss tag so the two levels never collide.
#
# How the matrices are built
# --------------------------
# Windows are extracted once for the union of site positions, and each (end type, class)
# matrix is then a read-weighted base tally over the relevant subset of rows. No positions are
# sampled, so the matrices are exact rather than Monte Carlo estimates, and a position that
# carries a thousand reads counts a thousand times as much as a position carrying one. This
# differs from 07_motif_by_occurrence.R, which subsamples positions and then weights them
# equally.
#
# Read weighting has a cost: a panel whose reads pile onto a handful of loci produces a logo
# that describes those loci rather than the end type. Each panel therefore reports Kish's
# effective sample size, (sum w)^2 / sum(w^2) — the number of equally weighted positions the
# panel is really worth. CAGE depth is always skewed, so even 1Sg over 327,535 positions has
# an ESS near 1,900; that is expected and not a defect. The pathological case is 1St | YR,
# which has 71,993 reads over 4,572 positions but an ESS of 1.4, so its near-fixed poly-T
# logo is one locus and nothing else.
#
# Because no single ESS cutoff separates "skewed as usual" from "one locus", both weightings
# are emitted for every panel rather than one being chosen:
#
#   reads     — each position contributes in proportion to its read count. This is what the
#               sequencer saw, and it is the right view when asking what the bulk of the
#               material of a given end type looks like.
#   positions — each site contributes once. This is immune to locus domination, and it is the
#               right view when asking what sequence context permits an end type at all.
#
# Where the two agree the motif is a property of the end type. Where they disagree the read-
# weighted panel is being carried by a few deep loci, and the ESS says how few.
#
# The alignment constraint at -1
# ------------------------------
# STAR ran with --alignEndsType Local, so a read is classified 1Sg only when the genomic base
# at -1 is not G (otherwise the aligner extends the match and the read becomes M). The same
# holds for 2Sgg with G, and for 1Sa / 1Sc / 1St with their own clip base. The -1 letter of
# every clipped-type logo is therefore depleted of the clip base by construction, and the
# class split inherits that depletion. M reads carry no clip and no constraint. This affects
# only the -1 column; it is annotated on the figure rather than corrected, since the class
# split itself already conditions on that base.
#
# Usage:
#   Rscript scripts/05_tss_end_types/10_motif_by_dinuc_class.R \
#       [--level      tss]                            \
#       [--tss        results/tss/tss_matrix.tsv]     \
#       [--ctss       results/tss/typed_ctss.tsv.gz]  \
#       [--samples    config/samples.tsv]             \
#       [--genome     data/sacCer3.fa]                \
#       [--outdir     results/motifs]                 \
#       [--figdir     results/figures]                \
#       [--types      1Sg,M,2Sgg,1Sa,1Sc,1St]         \
#       [--upstream   15]                             \
#       [--downstream 15]                             \
#       [--min_reads  200]                            \
#       [--min_ess    100]                            \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(readr)
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

level      <- get_arg("--level",      "ctss")
ctss_f     <- get_arg("--ctss",       "results/tss/typed_ctss.tsv.gz")
tss_f      <- get_arg("--tss",        "results/tss/tss_matrix.tsv")
samples_f  <- get_arg("--samples",    "config/samples.tsv")
genome_f   <- get_arg("--genome",     NULL)
outdir     <- get_arg("--outdir",     "results/motifs")
figdir     <- get_arg("--figdir",     "results/figures")
types_arg  <- get_arg("--types",      "1Sg,M,2Sgg,1Sa,1Sc,1St")
upstream   <- as.integer(get_arg("--upstream",   "15"))
downstream <- as.integer(get_arg("--downstream", "15"))
min_reads  <- as.integer(get_arg("--min_reads",  "200"))
min_ess    <- as.numeric(get_arg("--min_ess",    "100"))
force      <- has_flag("--force")

validate_args()

if (is.na(upstream)   || upstream   < 1L) stop("--upstream must be a positive integer")
if (is.na(downstream) || downstream < 1L) stop("--downstream must be a positive integer")
if (is.na(min_reads)  || min_reads  < 0L) stop("--min_reads must be a non-negative integer")
if (is.na(min_ess)    || min_ess    <= 0) stop("--min_ess must be a positive number")
if (!level %in% c("tss", "ctss")) stop("--level must be 'tss' or 'ctss', got: ", level)

TYPES <- trimws(strsplit(types_arg, ",", fixed = TRUE)[[1]])
TYPES <- TYPES[nzchar(TYPES)]
if (length(TYPES) == 0) stop("--types is empty")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# The CTSS level is the default, so it keeps the plain output names; a TSS-level run is tagged
# so that the two levels cannot overwrite each other's tables and figures.
sfx        <- if (level == "tss") "_tss" else ""

out_pwm    <- file.path(outdir, paste0("motif_by_dinuc_class_pwm", sfx, ".tsv"))
out_counts <- file.path(outdir, paste0("motif_by_dinuc_class_counts", sfx, ".tsv"))

if (!force && file.exists(out_pwm)) {
  message("Output exists, skipping (use --force to rerun): ", out_pwm)
  quit(save = "no", status = 0)
}

if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  stop("ggseqlogo not installed. Run: install.packages('ggseqlogo')")
}
library(ggseqlogo)

YR_RED   <- "#e63946"   # canonical initiator colour (project convention)
OTHER_GY <- "#6c757d"

BASES      <- c("A", "C", "G", "T")
CLASSES    <- c("YR", "YY", "RR", "RY")
PANELS     <- c("all", CLASSES)
WEIGHTINGS <- c("reads", "positions")

# What one row of the site table is, in words, for titles and captions.
SITE_PL <- if (level == "tss") "TSSs" else "CTSS positions"
SITE_SG <- if (level == "tss") "TSS"  else "CTSS position"

# The clip base each type implies at -1 under local alignment; M is unconstrained.
# END_TYPE_CLIP_BASE is defined once in scripts/utils/end_types.R over the full
# END_TYPES_REAL vocabulary (see the note there on why it's the label's LAST character).

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

# ── load the site table and combine samples ───────────────────────────────────
# Under --level tss each site is a distclu tag cluster: the window is anchored on the cluster's
# dominant position and weighted by the typed counts of the whole cluster, so a promoter with a
# spread of start sites counts once rather than once per base it fires from.
# Under --level ctss each site is a single-base CTSS, which is the older behaviour.
tcols <- paste0("n_", TYPES)

check_types <- function(dt, src) {
  missing <- setdiff(tcols, colnames(dt))
  if (length(missing) > 0) {
    stop("Missing end-type columns in ", src, ": ", paste(missing, collapse = ", "),
         "\n  Available: ", paste(grep("^n_", colnames(dt), value = TRUE), collapse = ", "))
  }
}

if (level == "tss") {
  message("Loading TSS matrix: ", tss_f)
  if (!file.exists(tss_f)) stop("TSS matrix not found: ", tss_f,
                                "\n  Run scripts/05_tss_end_types/run.sh first, or pass --level ctss.")
  tss <- fread(tss_f)
  message("  rows: ", nrow(tss))
  check_types(tss, tss_f)
  message("Combining samples over tag clusters...")
  # Tag clusters are disjoint within a (chr, strand), so their dominant positions are unique and
  # this grouping only collapses the per-sample rows of one cluster.
  agg <- tss[, lapply(.SD, sum), by = .(chr, pos = tss_pos, strand), .SDcols = tcols]
  rm(tss); invisible(gc(verbose = FALSE))
} else {
  message("Loading typed CTSS: ", ctss_f)
  ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
  message("  rows: ", nrow(ctss))
  check_types(ctss, ctss_f)
  message("Combining samples...")
  agg <- ctss[, lapply(.SD, sum), by = .(chr, pos, strand), .SDcols = tcols]
  rm(ctss); invisible(gc(verbose = FALSE))
}
message("  unique ", SITE_PL, ": ", nrow(agg))

# Drop sites carrying no read of any requested type; they cost window extraction time
# and contribute nothing to any matrix.
keep_any <- rowSums(as.matrix(agg[, ..tcols])) > 0
agg <- agg[keep_any]
message("  ", SITE_PL, " with at least one read of a requested type: ", nrow(agg))

# ── window extraction ─────────────────────────────────────────────────────────
# +1 is the first templated base and sits at `pos`; -1 is the base immediately upstream.
# There is no zero position, so the window covers -upstream..-1 and +1..+downstream.
w        <- upstream + downstream
pos_lab  <- c(seq.int(-upstream, -1L), seq.int(1L, downstream))
stopifnot(length(pos_lab) == w)

is_plus <- agg$strand != "-"
start_v <- ifelse(is_plus, agg$pos - upstream, agg$pos - downstream + 1L)

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

# Windows running off a chromosome end are dropped rather than trimmed, so that every column
# of a matrix shares one denominator and each column sums to exactly 1.
lens <- seqlengths(gr)[as.character(seqnames(gr))]
inb  <- start(gr) >= 1L & end(gr) <= lens
if (!all(inb)) {
  message("  dropping ", sum(!inb), " positions whose window runs off a chromosome end")
  gr <- gr[inb]; agg <- agg[inb]
}

message("Extracting ", w, "-base windows for ", length(gr), " positions...")
seqs <- getSeq(genome, gr)   # getSeq is strand-aware and returns 5'->3'
rm(gr); invisible(gc(verbose = FALSE))

# Windows with an ambiguity code are dropped whole, for the same denominator reason.
has_n <- alphabetFrequency(seqs, baseOnly = TRUE)[, "other"] > 0
if (any(has_n)) {
  message("  dropping ", sum(has_n), " windows containing non-ACGT bases")
  seqs <- seqs[!has_n]; agg <- agg[!has_n]
}
if (length(seqs) == 0) stop("No usable windows remain after filtering")

# One integer code matrix (1 = A, 2 = C, 3 = G, 4 = T) reused by every weighted tally below.
# It is filled column by column so that peak memory stays at one column-worth of characters
# rather than the whole window set expanded at once.
message("Encoding bases...")
code <- matrix(NA_integer_, nrow = length(seqs), ncol = w)
for (j in seq_len(w)) {
  code[, j] <- match(as.character(subseq(seqs, j, j)), BASES)
}
rm(seqs); invisible(gc(verbose = FALSE))
stopifnot(!anyNA(code))

# ── initiator dinucleotide class per position ─────────────────────────────────
# The initiator is (-1, +1), i.e. the last upstream column and the first downstream column of
# the extracted window, so it needs no second genome query.
i_minus1 <- upstream          # column of -1
i_plus1  <- upstream + 1L     # column of +1

is_pyr   <- function(b) b %in% c(2L, 4L)   # C or T
di_class <- factor(paste0(ifelse(is_pyr(code[, i_minus1]), "Y", "R"),
                          ifelse(is_pyr(code[, i_plus1]),  "Y", "R")),
                   levels = CLASSES)
agg[, dinuc_class := di_class]

message("Initiator class composition over all retained positions:")
print(agg[, .N, by = dinuc_class][order(-N)])

# ── read-weighted base matrices ───────────────────────────────────────────────
# For one group of rows the matrix entry (base b, column j) is the total number of reads whose
# window carries base b at column j, divided by the group's read total.
weighted_pfm <- function(rows, weight) {
  m <- matrix(0, nrow = 4L, ncol = w, dimnames = list(BASES, NULL))
  if (length(rows) == 0) return(m)
  wv <- weight[rows]
  for (j in seq_len(w)) {
    s <- rowsum(wv, code[rows, j], reorder = FALSE)
    m[as.integer(rownames(s)), j] <- m[as.integer(rownames(s)), j] + as.numeric(s)
  }
  m
}

message("\nBuilding matrices under both weightings...")
pfms     <- list()   # keyed "end_type|panel|weighting"
depth_dt <- list()   # keyed "end_type|panel"

for (typ in TYPES) {
  wv <- as.numeric(agg[[paste0("n_", typ)]])
  for (pn in PANELS) {
    rows <- if (pn == "all") which(wv > 0) else which(wv > 0 & agg$dinuc_class == pn)
    n_rd <- sum(wv[rows])
    key  <- paste(typ, pn, sep = "|")

    pfms[[paste(key, "reads",     sep = "|")]] <- weighted_pfm(rows, wv)
    pfms[[paste(key, "positions", sep = "|")]] <-
      weighted_pfm(rows, rep(1, length(wv)))

    # Kish's effective sample size: the number of equally weighted positions the read-weighted
    # panel is really worth. It equals n_positions when depth is even and collapses towards 1
    # when one locus dominates.
    ess <- if (length(rows) == 0) 0 else sum(wv[rows])^2 / sum(wv[rows]^2)

    depth_dt[[key]] <- data.table(end_type = typ, panel = pn,
                                  n_reads = n_rd, n_positions = length(rows),
                                  ess     = ess,
                                  usable  = n_rd >= min_reads & length(rows) >= min_reads,
                                  read_weighting_ok = ess >= min_ess)
    message(sprintf("  %-5s %-3s  reads = %12s   positions = %8s   ESS = %9.1f%s",
                    typ, pn, format(n_rd, big.mark = ","),
                    format(length(rows), big.mark = ","), ess,
                    if (n_rd < min_reads || length(rows) < min_reads)
                      "   [below --min_reads: too little data either way]"
                    else if (ess < min_ess)
                      "   [ESS below --min_ess: read a few loci; use the position-weighted panel]"
                    else ""))
  }
}

depth <- rbindlist(depth_dt)
depth[, level := level]
fwrite(depth, out_counts, sep = "\t")
message("Written: ", out_counts)

# ── long PWM table ────────────────────────────────────────────────────────────
pwm_rows <- lapply(names(pfms), function(key) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  m     <- pfms[[key]]
  tot   <- colSums(m)
  freq  <- if (all(tot > 0)) sweep(m, 2, tot, "/") else m * NA_real_
  data.table(end_type   = parts[1],
             panel      = parts[2],
             weighting  = parts[3],
             position   = rep(pos_lab, each = 4L),
             nucleotide = rep(BASES, times = w),
             count      = as.numeric(m),
             freq       = as.numeric(freq))
})
pwm_long <- rbindlist(pwm_rows)
pwm_long <- merge(pwm_long, depth, by = c("end_type", "panel"), sort = FALSE)
setorderv(pwm_long, c("end_type", "panel", "weighting", "position", "nucleotide"))
fwrite(pwm_long, out_pwm, sep = "\t")
message("Written: ", out_pwm)

# ── figures ───────────────────────────────────────────────────────────────────
# ggseqlogo takes a probability matrix directly, so the logos are drawn from the exact
# read-weighted frequencies without expanding anything back into sequences.
to_ppm <- function(m) {
  tot <- colSums(m)
  if (any(tot == 0)) return(NULL)
  mm <- sweep(m, 2, tot, "/")
  rownames(mm) <- BASES
  mm
}

# The title carries the two counts needed to judge the panel: how many alignments it rests on,
# and how many distinct CTSSs those alignments come from. The effective sample size behind the
# CONCENTRATED flag stays in the counts table rather than on the figure, because the flag is the
# part a reader acts on and the number itself needs a paragraph to interpret.
unit_label <- function(w) if (w == "reads") "alignments" else SITE_PL
wt_label   <- function(w) if (w == "reads") "alignment" else SITE_SG

panel_title <- function(d, weighting) {
  lab <- if (d$panel == "all") "all initiators" else d$panel
  n   <- if (weighting == "reads") d$n_reads else d$n_positions
  unit <- unit_label(weighting)
  flag <- if (!d$usable) ", LOW"
          else if (weighting == "reads" && !d$read_weighting_ok) ", CONCENTRATED"
          else ""
  # The alignment-weighted panel also names the number of sites those alignments come from,
  # because an alignment count alone does not say whether the logo rests on many loci or a few.
  extra <- if (weighting == "reads")
             paste0(" over ", format(d$n_positions, big.mark = ","), " ", SITE_PL) else ""
  paste0(d$end_type, " | ", lab, "  (", format(n, big.mark = ","), " ", unit, extra,
         flag, ")")
}

x_breaks <- which(pos_lab %% 5 == 0 | pos_lab %in% c(-1L, 1L))

logo_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid   = element_blank(),
        strip.text   = element_text(face = "bold", size = 8),
        axis.text.x  = element_text(size = 6),
        plot.title   = element_text(face = "bold", size = 11),
        plot.caption = element_text(hjust = 0, size = 6.5, colour = OTHER_GY))

weight_blurb <- function(weighting) {
  if (weighting == "reads") {
    paste0("Each position contributes in proportion to its alignment count, so this is what the ",
           "bulk of the sequenced material of that end type looks like.")
  } else {
    paste0("Each ", SITE_SG, " contributes once regardless of depth, so this is what sequence ",
           "context permits that end type at all, immune to a few deep loci.")
  }
}

# One sentence saying what a site is, since the whole figure hangs on it.
site_blurb <- if (level == "tss") {
  paste0("Each site is a distclu tag cluster from ", basename(tss_f), ": the window is anchored ",
         "on the cluster's dominant position and carries the typed counts of the whole cluster, ",
         "so a promoter counts once rather than once per base it fires from.")
} else {
  paste0("Each site is a single-base CTSS from ", basename(ctss_f),
         ", so a promoter that fires from several bases contributes one window per base.")
}

cap_for <- function(weighting) paste0(
  "Position +1 is the first templated base; there is no zero position, and the initiator ",
  "dinucleotide is (-1, +1). Weighting: ", toupper(unit_label(weighting)), ". ", weight_blurb(weighting),
  "\nEvery ", SITE_SG, " at which the end type occurs is used; nothing is sampled, so the ",
  "matrices are ",
  "exact. The 'all initiators' column is the union of the four class columns.\n",
  "CONCENTRATED marks a panel whose alignments pile onto so few loci that the alignment-weighted ",
  "logo describes those loci rather than the end type; read its ", SITE_SG, "-weighted companion ",
  "instead.\nThe effective sample size that sets the mark is in results/motifs/*_counts.tsv. ",
  "LOW marks fewer than ", format(min_reads, big.mark = ","), " alignments or ", SITE_PL, ".\n",
  "STAR ran with --alignEndsType Local, so for every soft-clipped type the genomic base at ",
  "-1 cannot equal the clip base (G for 1Sg and 2Sgg, A/C/T for 1Sa/1Sc/1St); the -1 letter ",
  "of those logos is depleted of that base by construction.\nM alignments carry no clip and no ",
  "such constraint. Genome: ", basename(genome_f), ". All samples combined.\n", site_blurb
)

build_list <- function(keys, weighting) {
  lst <- list()
  for (key in keys) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    d     <- depth[end_type == parts[1] & panel == parts[2]]
    ppm   <- to_ppm(pfms[[paste(key, weighting, sep = "|")]])
    if (is.null(ppm)) next
    lst[[panel_title(d, weighting)]] <- ppm
  }
  lst
}

# One grid per weighting: rows are end types, the five columns are 'all' then YR, YY, RR, RY.
grid_keys <- unlist(lapply(TYPES, function(t) paste(t, PANELS, sep = "|")))

for (weighting in WEIGHTINGS) {
  grid_list <- build_list(grid_keys, weighting)

  p_grid <- ggseqlogo(grid_list, method = "bits", ncol = length(PANELS)) +
    scale_x_continuous(breaks = x_breaks, labels = pos_lab[x_breaks]) +
    labs(title = paste0("Sequence context by 5' end type and initiator dinucleotide class (",
                        wt_label(weighting), "-weighted)"),
         x = "Position relative to the 5' end", y = "Bits", caption = cap_for(weighting)) +
    logo_theme

  out_g <- file.path(figdir, paste0("15_motif_by_dinuc_class_", weighting, sfx, ".pdf"))
  ggsave(out_g, p_grid,
         width = 3.6 * length(PANELS), height = 1.9 * length(TYPES) + 1.8,
         limitsize = FALSE)
  message("Written: ", out_g)
}

# One page per end type, both weightings interleaved, so a single type can be read at full
# width and the two views of the same panel sit next to each other.
for (typ in TYPES) {
  lst <- list()
  for (pn in PANELS) {
    for (weighting in WEIGHTINGS) {
      d   <- depth[end_type == typ & panel == pn]
      ppm <- to_ppm(pfms[[paste(typ, pn, weighting, sep = "|")]])
      if (is.null(ppm)) next
      lst[[paste0(panel_title(d, weighting), "  [", unit_label(weighting), "]")]] <- ppm
    }
  }
  if (length(lst) == 0) next

  p <- ggseqlogo(lst, method = "bits", ncol = length(WEIGHTINGS)) +
    scale_x_continuous(breaks = x_breaks, labels = pos_lab[x_breaks]) +
    labs(title = paste0(typ, ": sequence context by initiator dinucleotide class"),
         subtitle = paste0("Left column alignment-weighted, right column ",
                           SITE_SG, "-weighted"),
         x = "Position relative to the 5' end", y = "Bits",
         caption = paste0(cap_for("reads"), "\nThe right-hand column of this page is ",
                          SITE_SG, "-weighted: ", weight_blurb("positions"))) +
    logo_theme

  out_t <- file.path(figdir, paste0("15_motif_", typ, "_by_dinuc_class", sfx, ".pdf"))
  ggsave(out_t, p, width = 7.5 * length(WEIGHTINGS),
         height = 1.9 * length(PANELS) + 2.0, limitsize = FALSE)
  message("Written: ", out_t)
}

# Difference heatmaps: within one class, every other end type against 1Sg. Holding the class
# fixed removes the class-composition confound from the comparison.
if ("1Sg" %in% TYPES && length(TYPES) > 1) {
  for (wt in WEIGHTINGS) {
    freq_dt <- pwm_long[weighting == wt,
                        .(end_type, panel, position, nucleotide, freq)]
    ref     <- freq_dt[end_type == "1Sg", .(panel, position, nucleotide, ref_freq = freq)]
    cmp     <- merge(freq_dt[end_type != "1Sg"], ref,
                     by = c("panel", "position", "nucleotide"))
    cmp[, diff := freq - ref_freq]
    cmp[, panel    := factor(panel, levels = PANELS)]
    cmp[, end_type := factor(end_type, levels = setdiff(TYPES, "1Sg"))]

    p_diff <- ggplot(cmp, aes(position, nucleotide, fill = diff)) +
      geom_tile() +
      facet_grid(end_type ~ panel) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                           midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish,
                           name = "frequency\nminus 1Sg") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
      theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(),
            strip.text = element_text(face = "bold", size = 8),
            plot.caption = element_text(hjust = 0, size = 6.5, colour = OTHER_GY)) +
      labs(title = paste0("Each end type minus 1Sg, at matched initiator dinucleotide class (",
                          wt_label(wt), "-weighted)"),
           x = "Position relative to the 5' end", y = "",
           caption = paste0("Columns share the initiator class, so a residual signal is a ",
                            "difference in context rather than in class composition. ",
                            "The scale is capped at +/-0.3.\nThe -1 column carries the local-",
                            "alignment constraint described on the logo figure and is not ",
                            "comparable between types. Weighting: ", toupper(unit_label(wt)), ". ",
                            weight_blurb(wt)))

    out_diff <- file.path(figdir, paste0("15_motif_diff_vs_1Sg_by_class_", wt, sfx, ".pdf"))
    ggsave(out_diff, p_diff, width = 15, height = 1.5 * length(unique(cmp$end_type)) + 2.2,
           limitsize = FALSE)
    message("Written: ", out_diff)
  }
}

message("\nDone.")
sessionInfo()
