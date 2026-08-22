#!/usr/bin/env Rscript
# Stage 19 step 3: genomic context motifs around the positions each exotic
# 5' end type dominates, one motif per type.
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
# A type is drawn only if it dominates at least exo_motif_min_positions
# positions. Below that a per-base frequency is noise, and a logo drawn from a
# handful of sites invites a claim the data cannot carry. Types below the floor
# are still written to the table, flagged rather than dropped, so that nothing
# silently disappears between step 1's call summary and this step's output.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/03_context_motifs.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
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
samples_f <- get_arg("--samples",   "config/samples.tsv")
outdir    <- get_arg("--outdir",    "results/exotic_origins")
figdir    <- get_arg("--figdir",    "results/figures")
params_f  <- get_arg("--params",    "config/params.yaml")
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))

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

# ── figure 1: one sequence logo per type ──────────────────────────────────────

# ggseqlogo indexes the x axis 1..2*flank, so the axis is relabelled to the
# TSS-relative frame the table is written in. Ticks are every fifth base plus
# the two flanking the TSS, which is where the initiator sits.
tick_at <- sort(unique(c(1L, flank, flank + 1L, 2L * flank,
                         which(offsets %% 5 == 0))))
relabel_x <- scale_x_continuous(breaks = tick_at, labels = offsets[tick_at],
                                expand = expansion(add = 0.5))

panel_name <- function(tp) sprintf("%s  (n = %s)", tp,
                                   format(types[type == tp]$n_positions,
                                          big.mark = ","))
seq_list <- setNames(lapply(drawn_types,
                            function(tp) as.character(seqs[dom$call == tp])),
                     vapply(drawn_types, panel_name, character(1)))

p1 <- ggseqlogo(seq_list, method = "bits", ncol = 2) +
  relabel_x +
  labs(x = "position relative to the dominated CTSS (+1 = first templated base)",
       y = "information content (bits)",
       title = "Genomic context of the positions each 5' end type dominates",
       subtitle = paste0("One motif per dominating type; ", ref,
                         " is the reference the exotics are read against. ",
                         "A clip that appears upstream of +1 is templated.")) +
  theme(axis.text.x = element_text(size = 6),
        strip.text = element_text(size = 9, face = "bold"))
h1 <- 1.6 + 1.5 * ceiling(length(drawn_types) / 2)
ggsave(file.path(figdir, "21_exo_context_logos.pdf"), p1,
       width = 13, height = min(h1, 40), limitsize = FALSE)
message("Written: ", file.path(figdir, "21_exo_context_logos.pdf"))

# ── figure 2: difference against the canonical type ───────────────────────────

# A logo is dominated by whatever the genome does everywhere — the AT-richness
# of a yeast promoter swamps the per-type signal — so the panel that carries
# the actual claim is the difference against the canonical type's own motif at
# the same offsets. Anything at zero is the background; anything away from zero
# is specific to this type's sites.
ref_freq <- tab[type == ref, .(offset, base, ref_freq = freq)]
delta <- merge(tab[drawn == TRUE & type != ref], ref_freq,
               by = c("offset", "base"))
delta[, d := freq - ref_freq]
delta[, type_label := factor(vapply(type, panel_name, character(1)),
                             levels = vapply(setdiff(drawn_types, ref),
                                             panel_name, character(1)))]
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
h2 <- 1.6 + 1.2 * ceiling(length(setdiff(drawn_types, ref)) / 2)
ggsave(file.path(figdir, "21_exo_context_delta.pdf"), p2,
       width = 13, height = min(h2, 40), limitsize = FALSE)
message("Written: ", file.path(figdir, "21_exo_context_delta.pdf"))

# The main table is the idempotency sentinel, so it is written LAST: an
# interrupted run or a failing ggsave must not leave a sentinel behind that
# makes the next run skip this step while its figures are missing.
fwrite(ic, file.path(outdir, "context_motif_ic.tsv"), sep = "\t")
message("Written: ", file.path(outdir, "context_motif_ic.tsv"))
fwrite(tab, out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()
