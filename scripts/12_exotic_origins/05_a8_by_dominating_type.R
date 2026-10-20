#!/usr/bin/env Rscript
# Stage 19 step 5: the A-8 saturation test applied to the set step 3 draws its
# context logos from, one bubble row per <5' end type, dominance group> pair.
#
# Stage 18 asks whether A at -8 is enriched at <end type, initiator
# dinucleotide> pairs, where a position joins a pair as soon as it carries
# a8_min_type_reads alignments of that type. A type's enrichment there is
# measured over every site its reads land on, whoever owns the site. This step
# asks the same question of the split step 3 draws: a type's callable positions
# are partitioned by who dominates them, and each part is tested on its own. A
# type whose A-8 signal is carried by the sites it owns, and one whose signal
# is inherited from the canonical sites it contaminates, look identical in
# Stage 18 and separate here.
#
# Two of exo_helpers.R's schemes are drawn, into two figures. Under "sg" the
# splitter is the canonical type, which asks whether a type's reads at
# canonical start sites differ from its reads elsewhere; under "self3" the
# splitter is the type itself, which asks whether its own sites differ from the
# canonical sites it contaminates and from the sites a third type owns, those
# last two kept apart rather than merged. The row unit is the pair, not Stage 18's
# <type, dinucleotide>: splitting these sets across sixteen dinucleotides as
# well would leave almost every cell below any defensible testing floor, so the
# dinucleotide is not conditioned on.
#
# Membership in a row is OCCURRENCE -- the position carries at least one
# alignment of the type -- while the group encodes dominance. Using dominance
# for membership as well would leave one group of every row empty. A type gets
# rows at all only once it dominates exo_motif_min_positions positions, which is
# the same floor that decides whether step 3 draws its logo, and each group is
# tested only if it reaches that floor itself.
#
# Three of Stage 18's four fixed choices carry over unchanged: the position
# under test (a8_obs_coord), the proximal background (a8_bg_proximal), and
# position weighting, where each position counts once. The fourth, the
# a8_min_type_reads membership rule, is not applied: a single alignment of a
# type at a callable position is what puts that position in the type's row, and
# the dominance cascade of 01_position_dominance.R has already decided the
# group. There is correspondingly no sensitivity configuration -- the point
# estimate has no membership knob left to vary, and the smaller sets are too
# small for a second background to say anything the bootstrap does not.
#
# Note that the logos weight a window by the type's alignment count there and
# these tests do not. That difference is deliberate: a frequency should follow
# the reads, whereas a binomial tail counts a position's genomic context once
# however many reads start on it, and weighting it by depth would inflate every
# count and shrink every p value with no added evidence. The alignment count is
# written to the table so the two figures stay reconcilable.
#
# The position set is taken from the same window extraction step 3 performs, at
# exo_motif_flank and with the same two exclusions (a window running off a
# contig end, a window containing an N), so that the bubbles and the logos are
# drawn from a bitwise identical set of positions.
#
# The bootstrap resamples BLOCKS rather than positions. Positions inside one
# promoter share most of their upstream window and are not independent, so a
# position bootstrap would report a CI that is too narrow -- the same reason
# Stage 18 resamples clusters. Blocks are built by chaining every retained
# position within a8_max_dist on the same chromosome and strand, which is the
# distclu rule Stage 8 clusters by, and they are built over ALL retained
# positions before any split by type or group.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/05_a8_by_dominating_type.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
#       [--ctss      results/dominance/fine_typed_ctss.tsv.gz]         \
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
  library(scales)
  library(grid)          # unit(), for the legend key height
  library(readr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag)
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
# a8_log2fc(), a8_pvalue(), bh_neg_log10() and bootstrap_log2fc() are Stage 18's
# own, sourced rather than reimplemented so that the two bubble figures are
# comparable by construction: the same effect size, the same exact binomial tail
# in log space, and the same BH step-up on -log10 p.
source(file.path(script_dir, "..", "11_a8_saturation", "a8_helpers.R"))
# group_labels() and assign_group(): the same splits the step 3 context logos
# are drawn by, so that a bubble row and a logo panel describe the same set of
# positions.
source(file.path(script_dir, "exo_helpers.R"))

params <- yaml.load_file(params_f)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "a8_by_dominating_type.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

flank    <- if (is.null(params$exo_motif_flank)) 15L else as.integer(params$exo_motif_flank)
min_n    <- if (is.null(params$exo_motif_min_positions)) 20L else
              as.integer(params$exo_motif_min_positions)
obs_c    <- as.integer(params$a8_obs_coord)
bg_c     <- as.integer(params$a8_bg_proximal)
fdr      <- as.numeric(params$a8_fdr)
n_boot   <- as.integer(params$a8_n_boot)
seed     <- as.integer(params$a8_seed)
max_dist <- as.integer(params$a8_max_dist)

# The window step 3 extracts must actually contain every coordinate the test
# reads. It does at the configured values, but the two are set independently in
# the params file, so the containment is asserted rather than assumed.
need <- c(obs_c, bg_c)
if (any(need > 0L) || any(abs(need) > flank)) {
  stop("The context window (exo_motif_flank = ", flank, ") does not contain ",
       "every tested coordinate: ", paste(sort(need), collapse = ", "),
       ".\nEither widen exo_motif_flank or narrow a8_bg_proximal.")
}
if (obs_c %in% bg_c) {
  stop("a8_obs_coord (", obs_c, ") is inside a8_bg_proximal, so the position ",
       "under test would be part of its own background.")
}
message("Testing coordinate ", obs_c, " against the background ",
        paste(sort(bg_c), collapse = ", "), "; a type needs ", min_n,
        " dominated positions to be tested")

# ── genome ────────────────────────────────────────────────────────────────────
# Carried over verbatim from 03_context_motifs.R, because this step must select
# and trim exactly the same positions that step does.
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

# ── the callable positions ────────────────────────────────────────────────────
# The set is every CALLABLE position, not only the dominated ones. Group B of
# each scheme is by definition made of positions the splitter does not
# dominate, and those lie outside the dominated partition step 3's tables are
# built on. below_floor positions are excluded as a set: they are not callable,
# so they can be assigned to neither group.

if (!file.exists(dom_f)) {
  stop("Position-level dominance table not found: ", dom_f,
       "\nRun scripts/12_exotic_origins/01_position_dominance.R first.")
}
dom <- fread(dom_f, select = c("position_id", "chr", "pos", "strand", "call"))
dom <- dom[call != "below_floor"]
if (nrow(dom) == 0) stop("No callable positions to test.")
dom_counts <- dom[call != "none", .N, by = .(type = call)]
setorder(dom_counts, -N)
message("Callable positions: ", nrow(dom), "; dominated by ",
        nrow(dom_counts), " types")

gr <- GRanges(dom$chr, IRanges(dom$pos, width = 1), strand = dom$strand)
gr <- harmonize(gr, genome)
stopifnot(length(gr) == nrow(dom))

win  <- trim(promoters(gr, upstream = flank, downstream = flank))
full <- width(win) == 2L * flank
if (sum(!full) > 0) {
  message("Dropping ", sum(!full), " position(s) whose window runs off a ",
          "contig end")
}
dom  <- dom[full]
seqs <- getSeq(genome, win[full])

has_n <- vcountPattern("N", seqs) > 0
if (any(has_n)) {
  message("Dropping ", sum(has_n), " window(s) containing N")
  dom  <- dom[!has_n]
  seqs <- seqs[!has_n]
}
stopifnot(length(seqs) == nrow(dom))

# ── the A indicators ──────────────────────────────────────────────────────────
# There is no coordinate zero: the window spans -flank..-1 then +1..+flank, so
# offsets and string indices differ by one across the TSS boundary and are
# mapped through the offsets vector rather than by arithmetic.
offsets <- c(seq(-flank, -1L), seq_len(flank))
idx_of  <- function(coord) {
  i <- match(as.integer(coord), offsets)
  if (anyNA(i)) stop("coordinate outside the window: ", coord[is.na(i)])
  i
}

mat <- as.matrix(seqs)          # n_positions x 2*flank, one base per cell
is_a <- mat == "A"
dom[, a_obs := as.integer(is_a[, idx_of(obs_c)])]
dom[, a_bg  := as.integer(rowSums(is_a[, idx_of(bg_c), drop = FALSE]))]
rm(mat, is_a)
n_bg_pos <- length(bg_c)

# ── bootstrap blocks ──────────────────────────────────────────────────────────
# Single-linkage chaining: a new block starts wherever the gap to the previous
# position on the same chromosome and strand exceeds max_dist. This is the same
# rule distclu applies in Stage 8, so a block is a promoter-scale footprint
# rather than an arbitrary bin, and the blocks are built over ALL retained
# callable positions before any split, so that two rows drawing on neighbouring
# positions in one promoter resample the same block.
setorder(dom, chr, strand, pos)
dom[, block := {
  cumsum(c(TRUE, diff(pos) > max_dist))
}, by = .(chr, strand)]
dom[, block_id := paste(chr, strand, block, sep = ":")]
message("Bootstrap blocks: ", uniqueN(dom$block_id), " over ", nrow(dom),
        " positions (max_dist = ", max_dist, ")")

# ── which positions carry which type ──────────────────────────────────────────
# A row of the figure is a <type, group> pair, and a position joins a type's row
# only where that type actually has an alignment. Membership is therefore
# occurrence, not dominance: dominance is what the GROUP encodes, and using it
# for membership as well would leave group B empty by construction.

if (!file.exists(ctss_f)) {
  stop("Fine-typed CTSS table not found: ", ctss_f,
       "\nRun scripts/10_dominant_endtype/run.sh first.")
}
ctss <- fread(ctss_f, select = c("chr", "pos", "strand", "type", "n"))
ctss[, position_id := paste(chr, strand, pos, sep = ":")]
pooled <- ctss[, .(n = sum(n)), by = .(position_id, type)][n > 0]
rm(ctss)

wt <- merge(pooled, dom[, .(position_id, call, a_obs, a_bg, block_id)],
            by = "position_id")
rm(pooled)

ref <- END_TYPE_CANONICAL
drawn_types <- dom_counts[N >= min_n]$type
if (!ref %in% drawn_types) {
  stop("The canonical type ", ref, " dominates fewer than ", min_n,
       " positions, so there is no reference row to read the others against.")
}
message("Types with a row: ", paste(drawn_types, collapse = ", "))
# wt itself is left unrestricted -- test_scheme() filters to its own `types`
# argument, which is drawn_types for the historical figures and END_TYPES_REAL
# for the all-types figures below.

# Short y-axis tags. The full sentence each one abbreviates is in the subtitle,
# because a row label has to stay short enough to read down a column.
short_tags <- function(scheme) {
  switch(scheme,
    sg    = c(A = paste0(ref, "-dominated"), B = paste0("not ", ref, "-dominated")),
    self  = c(A = "self-dominated", B = "not self-dominated"),
    self3 = c(A = "self-dominated", B = paste0(ref, "-dominated"),
              C = "third-type- or un-dominated"),
    stop("unknown scheme: ", scheme))
}

# The point shape carries the group, so that the halves of a type can be told
# apart without reading the axis text. Open shapes with a settable fill only:
# the fill encodes the sign of the effect and the outline encodes the bootstrap.
GROUP_SHAPES <- c(21, 22, 24)

# ── the per-row test ──────────────────────────────────────────────────────────
# Each position counts ONCE in its row, exactly as Stage 18 and the previous
# version of this figure count it. The logos of step 3 weight a window by the
# type's alignment count there, which is the right currency for a frequency;
# it is the wrong one for a binomial tail, where a position's genomic context
# is one observation however many reads start on it, and weighting by depth
# would inflate every count and shrink every p value with no added evidence.
# The alignment count is carried in the table for reference instead.

test_scheme <- function(scheme, types) {
  labs_ab <- group_labels(scheme, ref)
  tags    <- short_tags(scheme)
  shape_names <- c("circles", "squares", "triangles")
  w <- copy(wt)[type %in% types]
  w[, grp := assign_group(as.character(call), as.character(type),
                          scheme = scheme, usual = ref)]

  sizes <- w[, .(n_positions  = .N,
                 n_alignments = sum(n),
                 n_blocks     = uniqueN(block_id),
                 n_a_obs      = sum(a_obs),
                 n_a_bg       = sum(a_bg)), by = .(type, grp)]
  sizes[, tested := n_positions >= min_n]
  setorder(sizes, -n_positions)

  bad <- sizes[tested == TRUE & (n_a_bg == 0L | n_a_obs > n_positions)]
  if (nrow(bad) > 0) {
    stop("Degenerate tally for ", paste(bad$type, bad$grp, collapse = ", "),
         ": a background with no A, or more A at ", obs_c, " than positions.")
  }

  res <- sizes[tested == TRUE]
  if (nrow(res) == 0) {
    stop("No <type, group> pair reaches ", min_n, " positions under scheme ",
         scheme, ".")
  }
  res[, p_exp  := n_a_bg / (n_bg_pos * n_positions)]
  res[, p_obs  := n_a_obs / n_positions]
  res[, log2fc := a8_log2fc(n_a_obs, n_a_bg, n_bg_pos)]
  res[, neg_log10_p := mapply(a8_pvalue, n_a_obs, n_positions, p_exp)]
  # The BH step-up runs over the rows of ONE scheme. The two schemes are two
  # readings of the same data and pooling them would correct each row for a
  # test it is not being read beside.
  res[, neg_log10_q := bh_neg_log10(neg_log10_p)]
  res[, significant := neg_log10_q >= -log10(fdr)]

  # The bootstrap is per row over that row's own blocks; the seed is offset by
  # the row index so that two rows do not share a resampling pattern.
  boot <- rbindlist(lapply(seq_len(nrow(res)), function(i) {
    bl <- w[type == res$type[i] & grp == res$grp[i],
            .(W = .N, A8 = sum(a_obs), BGA = sum(a_bg)), by = block_id]
    b  <- bootstrap_log2fc(bl$W, bl$A8, bl$BGA, n_bg_pos, n_boot, seed + i)
    b  <- b[is.finite(b)]
    if (length(b) == 0L) {
      return(data.table(type = res$type[i], grp = res$grp[i],
                        boot_lo = NA_real_, boot_hi = NA_real_, boot_n = 0L))
    }
    q <- unname(quantile(b, c(0.025, 0.975)))
    data.table(type = res$type[i], grp = res$grp[i], boot_lo = q[1],
               boot_hi = q[2], boot_n = length(b))
  }))
  res <- merge(res, boot, by = c("type", "grp"), sort = FALSE)
  # Robust means the percentile CI lies wholly on one side of zero. A resample
  # that produced no finite estimate cannot support the claim, so it is FALSE
  # rather than NA: the outline is a visual cue, not a test result.
  res[, boot_robust := !is.na(boot_lo) & !is.na(boot_hi) &
                       (boot_lo > 0 | boot_hi < 0)]

  tab <- merge(sizes, res[, .(type, grp, p_exp, p_obs, log2fc, neg_log10_p,
                              neg_log10_q, significant, boot_lo, boot_hi,
                              boot_n, boot_robust)],
               by = c("type", "grp"), all.x = TRUE)
  tab[, `:=`(scheme = scheme, group_label = labs_ab[as.character(grp)],
             obs_coord = obs_c, n_bg_positions = n_bg_pos)]
  setcolorder(tab, c("scheme", "type", "grp", "group_label", "tested",
                     "n_positions", "n_alignments", "n_blocks", "obs_coord",
                     "n_a_obs", "n_a_bg", "n_bg_positions", "p_obs", "p_exp",
                     "log2fc", "neg_log10_p", "neg_log10_q", "significant",
                     "boot_lo", "boot_hi", "boot_n", "boot_robust"))
  setorder(tab, -n_positions)

  skipped <- sizes[tested == FALSE]
  if (nrow(skipped) > 0) {
    message("  below the ", min_n, "-position floor, tabulated but not tested: ",
            paste(skipped$type, skipped$grp, sep = "/", collapse = ", "))
  }

  # ── the bubble figure ───────────────────────────────────────────────────────
  # Every encoding is Stage 18's, so that the figures can be read side by side:
  # area encodes |log2FC|, the fill its sign, the dotted line the FDR threshold,
  # and a thick outline a bootstrap CI clear of zero. The group is the one new
  # dimension and it is carried by the point shape as well as by the row label,
  # so the two halves of a type can be found without reading the axis text.
  #
  # Stage 18 puts -log10(q) on a sqrt scale, which works there because its 74
  # pairs are spread across the range. This set is bimodal instead -- the
  # canonical rows run into the hundreds and several exotic rows sit below 1.6 --
  # so under sqrt the small rows collapse onto the axis and the one that clears
  # the FDR threshold cannot be told from those that do not, which is precisely
  # the comparison this figure exists to make. A pseudo-log transform separates
  # them while staying monotone and uncapped, so no value is hidden.
  d <- copy(res)
  d[, label := sprintf("%s  |  %s   (n = %s)", type, tags[as.character(grp)],
                       format(n_positions, big.mark = ","))]
  d[, abs_fc := abs(log2fc)]
  d[, group_label := labs_ab[as.character(grp)]]
  setorder(d, -neg_log10_q)
  d[, label := factor(label, levels = rev(d$label))]

  p <- ggplot(d, aes(x = neg_log10_q, y = label)) +
    geom_vline(xintercept = -log10(fdr), linetype = "dotted", colour = "grey50") +
    geom_point(aes(size = abs_fc, fill = log2fc, shape = group_label,
                   stroke = fifelse(boot_robust, 0.7, 0.2)),
               colour = "grey20") +
    scale_shape_manual(values = GROUP_SHAPES[seq_along(labs_ab)],
                       breaks = unname(labs_ab),
                       name = NULL, guide = guide_legend(order = 3)) +
    scale_size_area(max_size = 9, name = expression(group("|", log[2]*FC, "|")),
                    guide = guide_legend(order = 1)) +
    scale_fill_gradient2(low = BLUE, mid = "grey95", high = YR_RED,
                         midpoint = 0, name = expression(log[2]*FC),
                         guide = guide_colourbar(order = 2)) +
    scale_x_continuous(name = expression(-log[10]~(adjusted~italic(p))),
                       trans  = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 10, 100, 1000)) +
    labs(y = NULL,
         # The scheme is named in the subtitle, not the title: spelled out in
         # full up here it runs past the panel on the "self" figure.
         title = paste0("A", obs_c, " enrichment by 5' end type and ",
                        "dominance group"),
         subtitle = paste0(
           "Same set as the step 3 context logos: the callable positions that\n",
           "carry the type, split into ",
           # One clause per group, so the two-group and three-group schemes are
           # described by the same code rather than by two hand-written strings.
           paste0(labs_ab, " (", shape_names[seq_along(labs_ab)], ")",
                  collapse = ",\n"), ".\n",
           min_n, "-position floor, position-weighted. Background ",
           min(bg_c), "..", max(bg_c), " excluding ", obs_c, ".\n",
           "Fill: red = enriched, blue = depleted. Thick outline = ",
           "block-bootstrap CI excludes 0.")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_line(linewidth = 0.2),
          plot.title = element_text(size = 9.5),
          plot.subtitle = element_text(size = 6.2, lineheight = 1.15),
          legend.key.height = unit(0.62, "lines"),
          legend.key.width = unit(0.62, "lines"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(3, "pt"),
          legend.spacing.y = unit(2, "pt"),
          axis.text = element_text(size = 6.5),
          plot.margin = margin(4, 4, 2, 2))

  list(tab = tab, plot = p, n_rows = nrow(d))
}

YR_RED <- "#e63946"   # project convention for canonical initiators
BLUE   <- "#2166ac"

# The width is deliberately tight. Most of the panel is empty by construction:
# the canonical rows sit far to the right of every other row, and that stretch
# is real and is not compressed away -- an axis break would hide exactly how far
# apart the groups are, which is the finding -- so the panel is made narrow
# instead, and the legends are shrunk so that they, rather than the data, stop
# setting the width.
# The self figure is drawn under "self3" rather than "self", to match the step 3
# context logos it is read beside: its right-hand group is split into the
# positions the canonical type dominates and the positions a third type or
# nothing does, because a minority end AT A CANONICAL START SITE and a minority
# end at a site some third type owns are the two readings this stage exists to
# tell apart. Unlike the logo grid, which draws a blank placeholder cell to keep
# its columns in register, a bubble row exists only where positions do, so the
# canonical type's empty middle group simply produces no row.
fig_schemes <- c("sg", "self3")

tabs <- list()
for (scheme in fig_schemes) {
  message("Scheme \"", scheme, "\": ",
          paste(group_labels(scheme, ref), collapse = "  vs  "))
  out <- test_scheme(scheme, drawn_types)
  tabs[[scheme]] <- out$tab
  # The "sg" scheme keeps the historical figure name, because it is the figure
  # the report already points at; "self" is written beside it.
  fig <- file.path(figdir, if (scheme == "sg") "23_exo_a8_bubble_by_q.pdf" else
                             "23_exo_a8_bubble_by_q_self.pdf")
  # The five-line subtitle, not the data, sets the height allowance here.
  ggsave(fig, out$plot, width = 6.6,
         height = max(3, 0.26 * out$n_rows + 2.2), limitsize = FALSE)
  message("Written: ", fig)
}

# ── same two schemes, every fine type rather than only drawn_types ────────────
# "sg"'s split is whether the CANONICAL type dominates a position, which has
# nothing to do with whether the row's own type self-dominates anything, and
# "self3"'s B/C groups are the same. Membership in a row was always
# OCCURRENCE (a_obs/a_bg from any callable position the type has a read at),
# never dominance; drawn_types only ever gated which types got a row at all.
# Dropping that whole-type gate in favour of END_TYPES_REAL lets more types
# reach test_scheme()'s existing per-<type, group> min_n floor (`tested`),
# which is untouched and still drops a thin cell rather than testing it.
tabs_all <- list()
for (scheme in fig_schemes) {
  out <- test_scheme(scheme, END_TYPES_REAL)
  tabs_all[[scheme]] <- out$tab
  fig <- file.path(figdir,
                   if (scheme == "sg") "23_exo_a8_bubble_by_q_all_types.pdf" else
                     "23_exo_a8_bubble_by_q_self_all_types.pdf")
  ggsave(fig, out$plot, width = 6.6,
         height = max(3, 0.26 * out$n_rows + 2.2), limitsize = FALSE)
  message("Written: ", fig)
}
fwrite(rbindlist(tabs_all),
       file.path(outdir, "a8_by_dominating_type_all_types.tsv"), sep = "\t")
message("Written: ", file.path(outdir, "a8_by_dominating_type_all_types.tsv"))

# The table is the idempotency sentinel, so it is written LAST: an interrupted
# run must not leave a sentinel behind that makes the next run skip this step
# while its figures are missing.
fwrite(rbindlist(tabs), out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()