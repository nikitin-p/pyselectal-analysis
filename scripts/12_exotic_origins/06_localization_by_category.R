#!/usr/bin/env Rscript
# Stage 19 step 6: genomic localization of EVERY position an end type occurs
# at, split by the dominance category of the position.
#
# Step 4 (04_dominated_localization.R) annotates only the positions a type
# dominates, which is the self_dom class of step 2. That answers what a type's
# OWN sites look like, and it leaves the rest of the type's signal unlocated:
# its reads at canonical start sites, its reads at sites a third type owns, and
# its reads at positions too shallow to be called at all. Step 3 already makes
# that contrast for the sequence context, by drawing one logo per dominance
# group, and this step makes it for the genome: one set of localization bars
# per category, so that a type's own sites can be read against the sites where
# it is a minority end.
#
# The split is the four categories of step 2 (CATEGORY_LEVELS in
# exo_helpers.R), assigned by classify_category() relative to the type whose
# bar is being drawn. That is the three-group "self3" split of the step 3
# logos with below_floor kept as a fourth column rather than dropped: a logo of
# unassessed positions invites a claim about dominance that cannot be made, but
# a localization bar only says where the positions are, and below the floor is
# where most of an exotic type's occurrence lies. The self_dom column of this
# figure is step 4's figure, redrawn as one column among four.
#
# The unit is the <end type, position> pair with at least one alignment of the
# type, so a position is drawn once in the bar of every type it carries, in
# the column of the category it falls into FOR THAT TYPE. The bars therefore
# overlap between types, which step 4's deliberately do not. Both currencies
# of step 2 are kept: an alignment count is the type's own depth at the
# position, and a position count is the number of distinct positions.
#
# The category set, the priority order, the promoter window and the palette
# come from scripts/utils/localization.R, as in step 4 and Stage 10. Enhancers
# are not detected here; Stage 10's enhancer BED is reused when it exists.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/06_localization_by_category.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
#       [--ctss      results/dominance/fine_typed_ctss.tsv.gz]          \
#       [--gtf       data/sacCer3.ensGene.gtf]                         \
#       [--enhancers results/localization/enhancer_regions.bed]        \
#       [--outdir    results/exotic_origins]                           \
#       [--figdir    results/figures]                                  \
#       [--params    config/params.yaml] [--fast] [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
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
gtf_f     <- get_arg("--gtf",       "data/sacCer3.ensGene.gtf")
enh_f     <- get_arg("--enhancers", "results/localization/enhancer_regions.bed")
outdir    <- get_arg("--outdir",    "results/exotic_origins")
figdir    <- get_arg("--figdir",    "results/figures")
params_f  <- get_arg("--params",    "config/params.yaml")
fast_mode <- "--fast" %in% args
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "exo_helpers.R"))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "..", "utils", "localization.R"))
source(file.path(script_dir, "..", "utils", "txdb.R"))
source(file.path(script_dir, "..", "utils", "plot_helpers.R"))

params <- yaml.load_file(params_f)
set_promoter_window(params)
promoter_cl <- load_promoter_clusters(params, get_arg("--promoters"))
top_n  <-as.integer(params$exo_top_types)
# The bar set follows the same rule as steps 2 and 4, so the figures line up
# type for type.
exo_type_set <- if (is.null(params$exo_type_set)) "top" else params$exo_type_set
if (!exo_type_set %in% c("vocabulary", "top"))
  stop("exo_type_set must be 'vocabulary' or 'top', not '", exo_type_set, "'")
exo_types <- if (exo_type_set == "vocabulary")
               c(setdiff(END_TYPES_REAL, END_TYPE_CANONICAL), "longclip") else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "localization_by_category.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

if (!fast_mode) {
  suppressPackageStartupMessages(library(ChIPseeker))
}

usual <- END_TYPE_CANONICAL

# ── every position, with its call ─────────────────────────────────────────────

# Unlike step 4, nothing is dropped: below_floor and "none" positions carry
# alignments of some type, and locating them is the point of this step.
pos <- fread(dom_f, select = c("position_id", "chr", "pos", "strand", "call"))
stopifnot(nrow(pos) == uniqueN(pos$position_id))

ctss <- fread(ctss_f, select = c("chr", "pos", "strand", "type", "n"))
ctss[, position_id := paste(chr, strand, pos, sep = ":")]
pooled <- ctss[, .(n = sum(n)), by = .(position_id, type)][n > 0]
rm(ctss)
message("Positions: ", nrow(pos), "; <end type, position> pairs: ", nrow(pooled))

# ── annotation ────────────────────────────────────────────────────────────────

message("Loading gene annotation: ", gtf_f)
txdb <- make_txdb_from_gff(gtf_f, format = "gtf")

pos_gr <- GRanges(seqnames = pos$chr,
                  ranges   = IRanges(start = pos$pos, width = 1),
                  strand   = pos$strand)

enhancer_gr <- GRanges()
if (!fast_mode && !is.null(enh_f) && nzchar(enh_f) && file.exists(enh_f)) {
  enhancer_gr <- import(enh_f, format = "BED")
  message("Reusing Stage 10 enhancer regions: ", length(enhancer_gr), " from ", enh_f)
} else if (!fast_mode) {
  message("No enhancer BED at ", enh_f,
          "; the Enhancer category will be absent from this figure")
}

# Each position is annotated once, whatever number of types it carries.
if (fast_mode) {
  message("Using fast mode (TxDb-based annotation)")
  tx <- transcripts(txdb)
  promoter_gr <- if (is.null(promoter_cl))
    promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN) else promoter_cl
  utr5_gr <- unlist(fiveUTRsByTranscript(txdb))
  utr3_gr <- unlist(threeUTRsByTranscript(txdb))
  cds_gr <- unlist(cdsBy(txdb, by = "tx"))
  intron_gr <- unlist(intronsByTranscript(txdb))
  pos[, category := annotate_positions_fast(pos_gr, promoter_gr, utr5_gr,
                                            utr3_gr, cds_gr, intron_gr)]
} else {
  message("Annotating ", nrow(pos), " positions with ChIPseeker...")
  pos[, category := annotate_positions_chipseeker(pos_gr, txdb, enhancer_gr,
                                                    promoter_gr = promoter_cl)]
}
stopifnot(!anyNA(pos$category))

# ── the <end type, position> pairs ────────────────────────────────────────────

n_before <- nrow(pooled)
d <- merge(pooled, pos[, .(position_id, call, category)], by = "position_id")
if (nrow(d) != n_before) {
  stop("The dominance table does not cover every position: ", n_before,
       " pairs before the join, ", nrow(d), " after.")
}
# `category` is the genomic category throughout this script, so the dominance
# category of exo_helpers.R is held as dom_class.
d[, dom_class := classify_category(as.character(call), as.character(type),
                                   usual = usual)]

# ── bars ──────────────────────────────────────────────────────────────────────

# The canonical reference first, then the pooled exotic bar, then the exotic
# types drawn individually with the tail pooled into one bar, as in panel A of
# step 2. A position is counted once per bar, so the pooled bars are smaller
# in positions than the sum of their members.
totals <- d[type != usual, .(n = sum(n)), by = type][order(-n)]
keep   <- exotic_bar_types(totals$type, top_n, exo_types)
d[, type_label := fifelse(type == usual, usual,
                  fifelse(type %in% keep, as.character(type), OTHER_TYPES_LABEL))]

tally <- function(x) x[, .(alignments = sum(n),
                           positions  = uniqueN(position_id)),
                       by = .(type_label, dom_class, category)]

per_type <- tally(d)
overall  <- tally(d[type != usual,
                    .(type_label = ALL_EXOTICS_LABEL, position_id, n,
                      dom_class, category)])
res <- rbindlist(list(per_type, overall), use.names = TRUE)

x_levels <- c(if (usual %in% d$type) usual, ALL_EXOTICS_LABEL, keep,
              if (OTHER_TYPES_LABEL %in% d$type_label) OTHER_TYPES_LABEL)
cat_order <- intersect(LOC_CATEGORY_ORDER, unique(d$category))

# The table is rectangular, so an absent combination reads as a measured zero.
grid <- CJ(type_label = x_levels, dom_class = CATEGORY_LEVELS,
           category = cat_order, unique = TRUE)
res[, dom_class := as.character(dom_class)]
res <- merge(grid, res, by = c("type_label", "dom_class", "category"),
             all.x = TRUE)
res[is.na(alignments), alignments := 0L]
res[is.na(positions),  positions  := 0L]
# An empty bar has no composition: the canonical type has no Sg_dom class by
# construction, and its percentages are NA rather than a division by zero.
pct <- function(x) if (sum(x) > 0) 100 * x / sum(x) else NA_real_
res[, pct_alignments := pct(alignments), by = .(type_label, dom_class)]
res[, pct_positions  := pct(positions),  by = .(type_label, dom_class)]
res[, dom_class_label := CATEGORY_LABELS[dom_class]]
res[, type_label := factor(type_label, levels = x_levels)]
res[, dom_class  := factor(dom_class, levels = CATEGORY_LEVELS)]
res[, category   := factor(category, levels = cat_order)]
setorder(res, type_label, dom_class, category)

# Every alignment is drawn exactly once in the per-type bars, the pooled bar
# repeats the exotics without adding to them, and the canonical type has no
# Sg_dom class.
per_bars <- setdiff(x_levels, ALL_EXOTICS_LABEL)
stopifnot(res[type_label %in% per_bars, sum(alignments)] == d[, sum(n)])
stopifnot(res[type_label == ALL_EXOTICS_LABEL, sum(alignments)] ==
            d[type != usual, sum(n)])
stopifnot(res[type_label == usual & dom_class == "Sg_dom", sum(alignments)] == 0)

# The self_dom column must reproduce step 4, which annotates the same
# positions under the same vocabulary. The comparison is skipped when step 4
# has not been run, or was run in the other annotation mode.
s4_f <- file.path(outdir, "dominated_localization.tsv")
if (file.exists(s4_f)) {
  s4 <- fread(s4_f)
  if (setequal(s4$category, as.character(cat_order))) {
    cmp <- merge(s4[type_label == usual, .(category, s4 = positions)],
                 res[type_label == usual & dom_class == "self_dom",
                     .(category = as.character(category), positions)],
                 by = "category")
    stopifnot(nrow(cmp) > 0, all(cmp$s4 == cmp$positions))
    message("self_dom column agrees with step 4 for ", usual)
  }
}

# ── figure ────────────────────────────────────────────────────────────────────

MEASURE_LABELS <- c(pct_alignments = "Alignments of the type",
                    pct_positions  = "Positions carrying the type")
plot_dt <- melt(res, id.vars = c("type_label", "dom_class", "category"),
                measure.vars = names(MEASURE_LABELS),
                variable.name = "measure", value.name = "pct")
plot_dt <- plot_dt[!is.na(pct)]

bar_n <- rbind(
  res[, .(measure = "pct_alignments", n = sum(alignments)),
      by = .(type_label, dom_class)],
  res[, .(measure = "pct_positions", n = sum(positions)),
      by = .(type_label, dom_class)])
bar_n <- bar_n[n > 0]
for (x in list(plot_dt, bar_n))
  x[, measure := factor(measure, levels = names(MEASURE_LABELS),
                        labels = MEASURE_LABELS)]

p_bar <- ggplot(plot_dt, aes(x = type_label, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(data = bar_n, aes(x = type_label, y = 101, label = fmt_n(n)),
            inherit.aes = FALSE, vjust = 0, size = 2.4, colour = "grey25") +
  # One row per dominance category and one column per currency, so a type's
  # four bars are read down a column.
  facet_grid(dom_class ~ measure,
             labeller = labeller(dom_class = as_labeller(CATEGORY_LABELS,
                                                         label_wrap_gen(24)))) +
  scale_fill_manual(values = LOC_CATEGORY_COLORS, name = "Category",
                    drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108)) +
  labs(title = paste0("Genomic localization of every position a 5' end type ",
                      "occurs at,\nby the dominance category of the position"),
       subtitle = paste0("Every <end type, position> pair with at least one ",
                         "alignment of the type",
                         if (fast_mode) " (fast mode)" else " (ChIPseeker)",
                         ";\nlabels give the total behind each bar, and a ",
                         "missing bar is an empty class"),
       x = "End type", y = "Percentage") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_f <- file.path(figdir, "22_exo_localization_by_category.pdf")
# The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
# The table is the sentinel, so the figure is written first.
for (f in c(fig_f, sub("\\.pdf$", ".png", fig_f)))
  ggsave(f, p_bar, width = 5 + 0.5 * length(x_levels),
         height = 2 + 2.6 * length(CATEGORY_LEVELS), dpi = 200, bg = "white")
message("Written: ", fig_f)

fwrite(res, out_tsv, sep = "\t")
message("Written: ", out_tsv)
print(dcast(res[category == "Promoter"], type_label ~ dom_class,
            value.var = "pct_positions"), digits = 3)

sessionInfo()
