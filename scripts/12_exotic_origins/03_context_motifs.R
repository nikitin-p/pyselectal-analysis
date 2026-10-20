#!/usr/bin/env Rscript
# Stage 19 step 3: genomic context motifs around the positions each exotic
# 5' end type occupies -- a pair of alignment-weighted logos per type, and a
# per-type difference panel against the canonical type.
#
# Steps 1 and 2 establish which CTSS positions each fine type dominates, and
# how productive those positions are. They say nothing about what the genome
# looks like there. This step asks that question directly: if a type is a start
# site's own signature rather than a contaminant of a 1Sg site, the sequence
# around that site should carry the type's own signal, and the shape of that
# signal is what separates a genuine untemplated 5' extension from a mapping
# artefact.
#
# The unit is the DOMINATED position, taken from step 1's per-position call, so
# each position contributes to exactly one type — the one that dominates it.
# That makes the motifs a partition of the callable dominated positions rather
# than an overlapping set, and it is the same self_dom class that panels B and
# C of the step 2 figure draw. Positions below the depth floor, and callable
# positions where nothing was called, carry no type and are excluded here:
# the question is what a type's OWN sites look like, and those positions have
# no owner.
#
# The canonical type is drawn on the same axis as the exotics, as the reference
# every exotic motif is read against. Its panel is what an ordinary CAGE start
# site looks like in this data, so a difference panel against it is the actual
# claim; the second figure draws exactly that difference.
#
# Coordinates are written TSS-relative and strand-aware, with +1 the first
# templated base (the TSS itself) and -1 the base immediately upstream. There
# is no zero, which is the usual promoter convention and is the same frame the
# initiator dinucleotide [-1, +1] is quoted in everywhere else in this
# pipeline. Under STAR --alignEndsType Local a clipped base mismatches the
# genome, so a type whose clip is templated shows its clip sequence sitting in
# the upstream positions -n..-1 — that is the artefact signature, and it is
# read straight off the logo.
#
# Three figures come out of this. The two logo figures split every type's
# alignments into groups and draw one logo per group, so that the contrast a
# reader is asked to make sits inside the figure: 21_exo_context_logos.pdf
# splits in two, on whether the canonical type dominates the position, and
# 21_exo_context_logos_self.pdf splits in three, on WHICH type dominates it --
# the type itself, the canonical type, or a third type (or none). The third
# column is what separates the contaminant reading of a minority end from the
# case where the site belongs to neither the type nor the canonical type. Their windows
# are weighted by alignment count and are narrower (exo_logo_flank) than the
# tables', because every feature they are read for lies within ten bases of the
# TSS. The third figure, 21_exo_context_delta.pdf, keeps the per-type dominated
# partition and the wider window, and is the difference against the canonical
# type's own motif.
#
# A type is drawn only if it dominates at least exo_motif_min_positions
# positions. Below that a per-base frequency is noise, and a logo drawn from a
# handful of sites invites a claim the data cannot carry. Types below the floor
# are still written to the table, flagged rather than dropped, so that nothing
# silently disappears between step 1's call summary and this step's output.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/03_context_motifs.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
#       [--ctss      results/dominance/fine_typed_ctss.tsv.gz]          \
#       [--samples   config/samples.tsv]                               \
#       [--outdir    results/exotic_origins]                           \
#       [--figdir    results/figures]                                  \
#       [--params    config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(ggseqlogo)
  library(readr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
dom_f     <- get_arg("--dominance", "results/exotic_origins/position_dominance.tsv.gz")
ctss_f    <- get_arg("--ctss",      "results/dominance/fine_typed_ctss.tsv.gz")
samples_f <- get_arg("--samples",   "config/samples.tsv")
outdir    <- get_arg("--outdir",    "results/exotic_origins")
figdir    <- get_arg("--figdir",    "results/figures")
params_f  <- get_arg("--params",    "config/params.yaml")
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
# LOGO_SCHEMES, group_labels() and assign_group(): the two-group split the
# alignment-weighted logos below are drawn by, kept beside the four-category
# classifier the step 2 bar plot uses so that the two figures partition the
# same positions in the same vocabulary.
source(file.path(script_dir, "exo_helpers.R"))

params <- yaml.load_file(params_f)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "context_motifs.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

flank <- if (is.null(params$exo_motif_flank)) 15L else as.integer(params$exo_motif_flank)
min_n <- if (is.null(params$exo_motif_min_positions)) 20L else
           as.integer(params$exo_motif_min_positions)
message("Window: ", flank, " bases each side of the dominated position; ",
        "a type needs ", min_n, " positions for its own panel")

# ── genome ────────────────────────────────────────────────────────────────────
# Loading and seqlevel harmonisation are carried over from
# scripts/10_dominant_endtype/04_context.R, which took them in turn from
# scripts/05_tss_end_types/02_cluster_annotate.R. Setting seqlengths from the
# FaFile is what lets trim() clip ranges that run off a contig end.
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec)) stop("FASTA not found: ", genome_spec)
    if (!file.exists(paste0(genome_spec, ".fai"))) indexFa(genome_spec)
    fa <- FaFile(genome_spec); open(fa); return(fa)
  }
  pkg <- switch(genome_spec,
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias: ", genome_spec))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("BSgenome not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

harmonize <- function(gr, genome) {
  fa_si <- seqinfo(genome)
  keep  <- intersect(seqlevels(gr), as.character(seqnames(fa_si)))
  gr    <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
  seqlengths(gr) <- seqlengths(fa_si)[seqlevels(gr)]
  gr
}

meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)
genome <- load_genome(unique(na.omit(meta$genome))[1])
on.exit(if (inherits(genome, "FaFile")) close(genome), add = TRUE)

# ── the dominated positions ───────────────────────────────────────────────────

if (!file.exists(dom_f)) {
  stop("Position-level dominance table not found: ", dom_f,
       "\nRun scripts/12_exotic_origins/01_position_dominance.R first.")
}
dom <- fread(dom_f, select = c("position_id", "chr", "pos", "strand", "call"))
dom <- dom[!call %in% c("below_floor", "none")]
if (nrow(dom) == 0) stop("No dominated positions to build motifs from.")
message("Dominated positions: ", nrow(dom), " across ", uniqueN(dom$call),
        " types")

gr <- GRanges(dom$chr, IRanges(dom$pos, width = 1), strand = dom$strand)
gr <- harmonize(gr, genome)
# harmonize() prunes ranges on seqlevels the FASTA does not carry, after which
# every dom[i, ...] assignment below would index dom by a position in the
# SHORTENED gr and attach a sequence to the wrong position. That never happens
# on sacCer3, where the FASTA carries every contig, so a mismatch is an error
# rather than a case to handle.
stopifnot(length(gr) == nrow(dom))

# promoters(upstream = flank, downstream = flank) spans [-flank .. -1, +1 ..
# +flank] in TSS-relative, strand-aware coordinates, and getSeq reverse
# complements the minus strand for us.
win <- trim(promoters(gr, upstream = flank, downstream = flank))
full <- width(win) == 2L * flank
n_trimmed <- sum(!full)
if (n_trimmed > 0) {
  message("Dropping ", n_trimmed, " position(s) whose window runs off a ",
          "contig end; a truncated window cannot be stacked with the rest")
}
dom <- dom[full]
seqs <- getSeq(genome, win[full])

# Windows containing an N would bias a per-base frequency toward the three
# bases actually present at that offset, so they are dropped as a unit and
# counted, rather than being let through with an unbalanced denominator.
has_n <- vcountPattern("N", seqs) > 0
if (any(has_n)) {
  message("Dropping ", sum(has_n), " window(s) containing N")
  dom  <- dom[!has_n]
  seqs <- seqs[!has_n]
}
stopifnot(length(seqs) == nrow(dom))

# ── per-type base frequencies ─────────────────────────────────────────────────

offsets <- c(seq(-flank, -1L), seq_len(flank))
BASES   <- c("A", "C", "G", "T")

freq_of <- function(idx) {
  cm <- consensusMatrix(seqs[idx], baseOnly = TRUE)[BASES, , drop = FALSE]
  data.table(offset = rep(offsets, each = length(BASES)),
             base   = rep(BASES, times = length(offsets)),
             count  = as.vector(cm))
}

types <- dom[, .(n_positions = .N), by = .(type = call)][order(-n_positions)]
tab <- rbindlist(lapply(seq_len(nrow(types)), function(i) {
  tp <- types$type[i]
  f  <- freq_of(which(dom$call == tp))
  f[, `:=`(type = tp, n_positions = types$n_positions[i])]
  f[, freq := count / n_positions]
  f
}))
tab[, drawn := n_positions >= min_n]
setcolorder(tab, c("type", "n_positions", "drawn", "offset", "base", "count",
                   "freq"))
setorder(tab, -n_positions, type, offset, base)

# Shannon information content per offset, in bits, on the same 2-bit maximum a
# DNA sequence logo's height is drawn against. This is the number the logos are
# read by, written out so a downstream consumer does not have to re-derive it
# from the frequencies.
ic <- tab[, {
  p <- freq[freq > 0]
  .(bits = 2 + sum(p * log2(p)))
}, by = .(type, n_positions, drawn, offset)]
setorder(ic, -n_positions, type, offset)

ref <- END_TYPE_CANONICAL
drawn_types <- types[n_positions >= min_n]$type
if (!ref %in% drawn_types) {
  stop("The canonical type ", ref, " dominates fewer than ", min_n,
       " positions, so there is no reference to read the exotics against.")
}
message("Types with their own panel: ", paste(drawn_types, collapse = ", "))
skipped <- setdiff(types$type, drawn_types)
if (length(skipped) > 0) {
  message("Below the ", min_n, "-position floor, tabulated but not drawn: ",
          paste(skipped, collapse = ", "))
}

# ── panel naming and the x axis shared by the figures below ───────────────────
# The offset frame the tables are written in has no zero, so a tick position is
# looked up through the offsets vector rather than computed. Ticks are every
# fifth base plus the two flanking the TSS, which is where the initiator sits.
tick_at <- sort(unique(c(1L, flank, flank + 1L, 2L * flank,
                         which(offsets %% 5 == 0))))
panel_name <- function(tp) sprintf("%s  (n = %s)", tp,
                                   format(types[type == tp]$n_positions,
                                          big.mark = ","))

# ── figures 1a/1b: alignment-weighted context logos, two groups per type ──────
#
# The logo above the fold is not the per-type motif this script used to draw.
# That motif answered "what does a type's own site look like", and answered it
# one position at a time -- a position dominated by a type contributed one
# window whether it carried three of that type's alignments or three hundred.
# Two changes make the figure say something the per-type motif could not.
#
# First, the window is split in two, by one of the two schemes exo_helpers.R
# defines. Each type is drawn as a PAIR of logos, so the contrast is inside the
# figure rather than in the reader's memory of another panel. Second, a window
# is weighted by the number of alignments OF THAT TYPE at its position, so the
# logo describes where a type's reads actually are rather than where its
# positions are. Those two currencies differ by more than an order of magnitude
# at the canonical type and are the same two the step 2 bar plot draws.
#
# The display window is exo_logo_flank rather than exo_motif_flank. Everything
# these logos are read for -- the clip sitting upstream of +1, the initiator,
# the A at -8 -- lies within ten bases of the TSS, and the wider window the
# tables and the difference panel use spends half its width on the flat AT-rich
# background that carries no per-type signal.
#
# Both figures are drawn over ALL callable positions, not only the dominated
# ones: group B exists precisely to hold the positions a type does not own, and
# those positions are outside the partition the tables above are built on.

logo_flank <- if (is.null(params$exo_logo_flank)) 10L else
                as.integer(params$exo_logo_flank)
if (logo_flank > flank) {
  stop("exo_logo_flank (", logo_flank, ") exceeds exo_motif_flank (", flank,
       "); the logos would be wider than the window the tables are built on.")
}

if (!file.exists(ctss_f)) {
  stop("Fine-typed CTSS table not found: ", ctss_f,
       "\nRun scripts/10_dominant_endtype/run.sh first.")
}
ctss <- fread(ctss_f, select = c("chr", "pos", "strand", "type", "n"))
ctss[, position_id := paste(chr, strand, pos, sep = ":")]
pooled <- ctss[, .(n = sum(n)), by = .(position_id, type)]
rm(ctss)

# The logo rows run from the most to the least abundant type, which is the
# order the columns of 03_replicate_heatmap_fine take, so that the two figures
# read in the same sequence. The abundance is pooled over every position, not
# only the callable ones, because that is the total Stage 3 counts.
type_abund <- pooled[, .(n = sum(n)), by = type]
by_abundance <- function(types) {
  n <- type_abund$n[match(types, type_abund$type)]
  types[order(-ifelse(is.na(n), 0, n))]
}

# Every callable position, which is the union of the two groups. below_floor is
# excluded here rather than in assign_group(), which refuses it.
call_all <- fread(dom_f, select = c("position_id", "chr", "pos", "strand",
                                    "call"))
call_all <- call_all[call != "below_floor"]
message("Callable positions: ", nrow(call_all), " (", uniqueN(pooled$type),
        " fine types carry alignments somewhere)")

gr_all <- GRanges(call_all$chr, IRanges(call_all$pos, width = 1),
                  strand = call_all$strand)
gr_all <- harmonize(gr_all, genome)
stopifnot(length(gr_all) == nrow(call_all))

win_all  <- trim(promoters(gr_all, upstream = logo_flank,
                           downstream = logo_flank))
full_all <- width(win_all) == 2L * logo_flank
if (sum(!full_all) > 0) {
  message("Dropping ", sum(!full_all), " callable position(s) whose logo ",
          "window runs off a contig end")
}
call_all <- call_all[full_all]
seqs_all <- getSeq(genome, win_all[full_all])

has_n_all <- vcountPattern("N", seqs_all) > 0
if (any(has_n_all)) {
  message("Dropping ", sum(has_n_all), " callable window(s) containing N")
  call_all <- call_all[!has_n_all]
  seqs_all <- seqs_all[!has_n_all]
}
stopifnot(length(seqs_all) == nrow(call_all))
call_all[, row_idx := .I]

# One character matrix, built once: the per-group tally is a column sum over a
# subset of its rows, which is cheaper and clearer than re-extracting sequences
# per group.
base_mat  <- as.matrix(seqs_all)
logo_offs <- c(seq(-logo_flank, -1L), seq_len(logo_flank))

# The type's own alignment count at every callable position it appears at.
# Positions where a type has no read contribute nothing to either of its
# groups: a logo of windows with zero weight would be a logo of nothing.
wt <- merge(pooled, call_all[, .(position_id, row_idx, call)],
            by = "position_id")
wt <- wt[n > 0]

# A weighted position frequency matrix, bases in rows, offsets in columns, in
# the counts ggseqlogo normalises per column before taking bits.
weighted_pfm <- function(rows, w) {
  m <- base_mat[rows, , drop = FALSE]
  out <- vapply(BASES, function(b) colSums(w * (m == b)), numeric(ncol(m)))
  t(out)
}

# The name of a matrix in the list handed to ggseqlogo. It is only a key: the
# figure labels its rows and its columns once each, on the outside of the grid,
# rather than repeating both on every cell, so this string is split back into
# its two fields and never drawn as written. PANEL_SEP is therefore chosen to
# be a sequence no type name and no group label can contain.
PANEL_SEP <- "@@"

panel_of <- function(tp, lab) paste0(tp, PANEL_SEP, lab)

# What the figure title says the split is, per scheme. The group labels are
# phrased as column captions and do not read as the object of "split by".
scheme_title <- function(scheme, usual) {
  switch(scheme,
    sg    = paste0("whether ", usual, " dominates the position"),
    self  = "whether this end type dominates the position",
    self3 = "which type dominates the position",
    stop("unknown scheme: ", scheme))
}

# ggseqlogo fills its grid row-wise from the order of the matrix list, so a
# type-major, group-minor list puts one type per row and one group per column
# only as long as EVERY type contributes a panel to every column. A type with
# an empty group would slide every later panel one cell to the left and put the
# figure's columns out of register with its labels, which under "self3" is not
# hypothetical: the canonical type's middle column is empty by construction.
# Such a cell is drawn from a uniform frequency matrix, which carries zero bits
# and therefore renders as a blank panel of the right size, labelled for what
# it is. It contributes no row to the written table -- there is no measurement
# behind it.
blank_pfm <- function() {
  m <- matrix(0.25, nrow = length(BASES), ncol = 2L * logo_flank,
              dimnames = list(BASES, NULL))
  # ggseqlogo replaces the information content of an ENTIRELY flat matrix with
  # the 2-bit maximum, on the assumption that a flat logo means too few input
  # sequences rather than a genuinely uniform one, which would draw this cell
  # as a full-height wall of every base -- the loudest panel in the figure.
  # One negligible perturbation puts the matrix outside that special case, and
  # 1e-6 of a base is far below the resolution of any glyph, so the cell draws
  # as the blank it is meant to be.
  m[1L, 1L] <- m[1L, 1L] + 1e-6
  m
}

# The page width a scheme's grid needs. A cell only has to be wide enough for
# its 2 * exo_logo_flank positions, and the per-cell titles that used to force
# the columns apart are gone, so a column is narrower than it was: the
# three-column self figure fits a page instead of running well past one. The
# constant covers the y axis on the left and the type strip on the right.
fig_width <- function(ncol) 0.9 + 3.2 * ncol

# The project's canonical red, the same one the initiator plots highlight the
# YR dinucleotides in. Here it marks offset -8 rather than an initiator.
A8_COL <- "#e63946"

draw_scheme <- function(scheme, types) {
  types   <- by_abundance(types)
  labs_ab <- group_labels(scheme, ref)
  groups  <- names(labs_ab)
  wt[, grp := assign_group(as.character(call), as.character(type),
                           scheme = scheme, usual = ref)]

  mats  <- list()
  rows  <- list()
  for (tp in types) {
    for (g in groups) {
      s <- wt[type == tp & grp == g]
      if (nrow(s) == 0L) {
        message("  ", tp, " / ", labs_ab[[g]], ": no alignments, panel blank")
        mats[[panel_of(tp, labs_ab[[g]])]] <- blank_pfm()
        next
      }
      nm <- panel_of(tp, labs_ab[[g]])
      mats[[nm]] <- weighted_pfm(s$row_idx, s$n)
      rows[[length(rows) + 1L]] <-
        data.table(scheme = scheme, type = tp, group = g,
                   group_label = labs_ab[[g]], n_positions = nrow(s),
                   n_alignments = sum(s$n),
                   offset = rep(logo_offs, each = length(BASES)),
                   base = rep(BASES, times = length(logo_offs)),
                   weight = as.vector(mats[[nm]]))
    }
  }
  # -8 is a tick of its own, and its label alone is drawn in the canonical red,
  # because it is the offset step 5's A-8 saturation is about: the base an A
  # there is claimed to act from is the one a reader comes to these logos to
  # check, and hunting for it between the -10 and -5 ticks is exactly the
  # reading error the mark removes. The mark stays on the axis -- nothing is
  # drawn inside the cells, where a rule would only compete with the glyphs.
  #
  # -1 and +1 are one base apart, which at this label size would set the two
  # numbers touching and read as a single "-11". The TSS boundary is therefore
  # marked once, by a tick halfway between the two columns labelled "-1|+1",
  # which is where the boundary actually falls.
  i_a8  <- which(logo_offs == -8L)
  i_num <- sort(unique(c(1L, i_a8, which(logo_offs %% 5 == 0),
                         2L * logo_flank)))
  minus <- function(x) sub("^-", "−", as.character(x))
  tick     <- c(as.numeric(i_num), logo_flank + 0.5)
  tick_lab <- c(minus(logo_offs[i_num]), "−1|+1")
  ord      <- order(tick)
  tick     <- tick[ord]
  tick_lab <- tick_lab[ord]
  tick_col <- ifelse(tick %in% i_a8, A8_COL, "grey20")
  # The figure is now only as wide as its cells need, so the title and the
  # subtitle no longer run off the page on their own: both are wrapped to the
  # page width, at roughly the characters per inch each of their point sizes
  # gives. The column captions the subtitle used to spell out have moved onto
  # the columns themselves, as the grid's top strips.
  page_w <- fig_width(length(groups))
  wrap_to <- function(txt, cpi) paste(strwrap(txt, width = cpi * page_w),
                                      collapse = "\n")
  p <- ggseqlogo(mats, method = "bits", ncol = length(groups)) +
    scale_x_continuous(breaks = tick, labels = tick_lab,
                       expand = expansion(add = 0.5)) +
    # Bits run 0-2 for a four-letter alphabet, so whole bits are the only
    # gradations worth a tick; intermediate labels only add ink to a row that
    # is now deliberately short.
    scale_y_continuous(breaks = 0:2, expand = expansion(mult = c(0, 0.05))) +
    labs(x = "position relative to the CTSS (+1 = first templated base)",
         y = "information content (bits)",
         title = wrap_to(
           paste0("Genomic context of each 5' end type's alignments, ",
                  "split by ", scheme_title(scheme, ref)), 12),
         subtitle = wrap_to(paste0(
           "Each window is weighted by that type's alignment count ",
           "at its position. Callable positions only; ", ref,
           " is the reference the exotics are read against. ",
           "A clip that appears upstream of +1 is templated. ",
           "Per-cell counts are in context_motifs_by_group.tsv."), 19)) +
    theme(axis.text.x = element_text(size = 12, colour = tick_col),
          plot.subtitle = element_text(size = 7))

  # ggseqlogo lays the list out with facet_wrap over the panel key, which can
  # only label a cell with the whole key and so repeats the type and the group
  # on every cell. The grid is a true crossing of type by group, so it is
  # relabelled as one: the key is split back into its two fields, facet_grid
  # draws the type once per row on the right and the group once per column on
  # top, and no cell carries a title of its own. The row and column order is
  # carried by the factor levels rather than by the list order, which also
  # makes the blank placeholder cells land where they belong.
  ld <- p$layers[[1]]$data
  keys <- do.call(rbind, strsplit(as.character(ld$seq_group), PANEL_SEP,
                                  fixed = TRUE))
  ld$panel_type  <- factor(keys[, 1], levels = types)
  ld$panel_group <- factor(keys[, 2], levels = unname(labs_ab))
  p$layers[[1]]$data <- ld
  p <- p +
    # axes = "all_x" repeats the coordinate axis under every cell instead of
    # only under the bottom row, so a logo in the middle of the grid can be
    # read against its offsets without counting rows down to the page's foot.
    facet_grid(panel_type ~ panel_group, scales = "free_x", axes = "all_x",
               labeller = labeller(panel_group = label_wrap_gen(width = 26))) +
    theme(strip.text.x.top = element_text(size = 11, face = "bold",
                                          margin = margin(b = 4)),
          strip.text.y.right = element_text(size = 15, face = "bold",
                                            angle = 0,
                                            margin = margin(l = 5)))
  list(plot = p, rows = rbindlist(rows), n_panels = length(mats),
       ncol = length(groups))
}

# The self figure is drawn under "self3" rather than "self": its right-hand
# column is split into the positions the canonical type dominates and the
# positions a third type or nothing does, because a minority end AT A CANONICAL
# START SITE and a minority end at a site some third type owns are the two
# readings this stage exists to tell apart, and they need not share a motif.
# Step 5's A-8 bubbles keep the two-group "self", where a third column would
# thin each panel past readability.
fig_schemes <- c("sg", "self3")

group_tab <- list()
for (scheme in fig_schemes) {
  message("Scheme \"", scheme, "\": ",
          paste(group_labels(scheme, ref), collapse = "  vs  "))
  out <- draw_scheme(scheme, drawn_types)
  group_tab[[scheme]] <- out$rows
  # The "sg" scheme keeps the historical figure name, because it is the figure
  # the report already points at; the self split is written beside it under the
  # name it has always had.
  fig <- file.path(figdir, if (scheme == "sg") "21_exo_context_logos.pdf" else
                             "21_exo_context_logos_self.pdf")
  # The header takes a fixed allowance on top of the rows of panels. A row is
  # kept short: a bits logo is read across, so the cells only need enough height
  # for the glyphs, and at 1.75 in a nine-type grid ran to eighteen inches of
  # page that no one can see at once.
  h <- 1.8 + 1.1 * ceiling(out$n_panels / out$ncol)
  ggsave(fig, out$plot, width = fig_width(out$ncol), height = min(h, 40),
         limitsize = FALSE)
  message("Written: ", fig)
}
fwrite(rbindlist(group_tab), file.path(outdir, "context_motifs_by_group.tsv"),
       sep = "\t")
message("Written: ", file.path(outdir, "context_motifs_by_group.tsv"))

# ── same two schemes, every fine type rather than only drawn_types ────────────
# "sg"'s split is whether the CANONICAL type dominates a position, which has
# nothing to do with whether the row's own type dominates anything; "self3"'s B
# and C columns are the same. The exo_motif_min_positions floor only ever
# controlled whether a type's OWN column (self3's A, and hence drawn_types
# membership) had enough positions to trust -- it was never a floor on the
# type's occurrence at Sg- or third-dominated sites. These figures draw every
# END_TYPES_REAL type; a type below the floor simply gets a blank self-dominated
# cell via blank_pfm(), the same mechanism that already blanks a zero-alignment
# cell.
group_tab_all <- list()
for (scheme in fig_schemes) {
  out <- draw_scheme(scheme, END_TYPES_REAL)
  group_tab_all[[scheme]] <- out$rows
  fig <- file.path(figdir,
                   if (scheme == "sg") "21_exo_context_logos_all_types.pdf" else
                     "21_exo_context_logos_self_all_types.pdf")
  h <- 1.8 + 1.1 * ceiling(out$n_panels / out$ncol)
  # The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
  for (f in c(fig, sub("\\.pdf$", ".png", fig)))
    ggsave(f, out$plot, width = fig_width(out$ncol), height = min(h, 40),
           limitsize = FALSE, dpi = 200, bg = "white")
  message("Written: ", fig)
}
fwrite(rbindlist(group_tab_all),
       file.path(outdir, "context_motifs_by_group_all_types.tsv"), sep = "\t")
message("Written: ", file.path(outdir, "context_motifs_by_group_all_types.tsv"))

# ── figure 2: difference against the canonical type ───────────────────────────

# A logo is dominated by whatever the genome does everywhere — the AT-richness
# of a yeast promoter swamps the per-type signal — so the panel that carries
# the actual claim is the difference against the canonical type's own motif at
# the same offsets. Anything at zero is the background; anything away from zero
# is specific to this type's sites.
#
# draw_delta_figure() takes the set of types to include as an argument so the
# same panel can be drawn twice: once restricted to drawn_types (n_positions
# >= min_n, i.e. enough of the type's own sites to trust the per-base
# frequency), and once over every type that self-dominates at least one
# position at all, min_n or not. A type absent from `types` entirely (zero
# self-dominated positions) still has no row in either version, because there
# is no sequence context to take a delta of -- that boundary is "no data", not
# "below the floor", and no relaxation of the floor changes it.
ref_freq <- tab[type == ref, .(offset, base, ref_freq = freq)]
draw_delta_figure <- function(type_set, fname) {
  ts <- intersect(setdiff(type_set, ref), unique(tab$type))
  delta <- merge(tab[type %in% ts], ref_freq, by = c("offset", "base"))
  delta[, d := freq - ref_freq]
  delta[, type_label := factor(vapply(type, panel_name, character(1)),
                               levels = vapply(ts, panel_name, character(1)))]
  lim <- max(abs(delta$d))

  # geom_tile rather than geom_raster: the offset frame has no zero, so the
  # columns are unevenly spaced by design and a raster would shift them to
  # close the gap the TSS boundary sits in.
  p2 <- ggplot(delta, aes(offset, base, fill = d)) +
    geom_tile() +
    facet_wrap(~ type_label, ncol = 2) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey20") +
    scale_fill_gradient2(low = "#457b9d", mid = "white", high = "#e63946",
                         midpoint = 0, limits = c(-lim, lim),
                         name = paste0("frequency\nminus ", ref)) +
    scale_x_continuous(breaks = offsets[tick_at], expand = c(0, 0)) +
    scale_y_discrete(limits = rev(BASES), expand = c(0, 0)) +
    labs(x = "position relative to the dominated CTSS (+1 = first templated base)",
         y = NULL,
         title = paste0("Context of each exotic type's own sites, minus the ",
                        ref, " reference"),
         subtitle = paste0("Red is enriched over ", ref, " sites, blue depleted; ",
                           "the line marks the TSS. A templated clip shows as a ",
                           "red run immediately upstream.")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 6),
          strip.text = element_text(face = "bold"))
  h2 <- 1.6 + 1.2 * ceiling(length(ts) / 2)
  ggsave(file.path(figdir, fname), p2, width = 13, height = min(h2, 40),
         limitsize = FALSE)
  message("Written: ", file.path(figdir, fname))
}
draw_delta_figure(drawn_types, "21_exo_context_delta.pdf")
draw_delta_figure(END_TYPES_REAL, "21_exo_context_delta_all_types.pdf")

# The main table is the idempotency sentinel, so it is written LAST: an
# interrupted run or a failing ggsave must not leave a sentinel behind that
# makes the next run skip this step while its figures are missing.
fwrite(ic, file.path(outdir, "context_motif_ic.tsv"), sep = "\t")
message("Written: ", file.path(outdir, "context_motif_ic.tsv"))
fwrite(tab, out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()
