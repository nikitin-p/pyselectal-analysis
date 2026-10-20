#!/usr/bin/env Rscript
# Stage 19 step 4: genomic localization of the CTSS positions each fine 5' end
# type dominates.
#
# Stage 10 (scripts/06_localization/01_localization.R) annotates the 5' end of
# every alignment of every end type, so its bars answer "where in the genome
# does this type's signal land". That is not the question the rest of this
# stage asks. Steps 1 to 3 restrict every claim to the positions a type
# DOMINATES, because a type's reads at a 1Sg-dominated site say nothing about
# the type — they are a minority end at somebody else's start site. This step
# puts the localization bars on that same restricted footing: the unit is the
# dominated position from step 1's per-position call, so each position
# contributes to exactly one type, the one that owns it, and the bars are a
# partition of the callable dominated positions rather than an overlapping set.
# It is the same self_dom class that panels B and C of step 2 draw and that
# step 3 builds its motifs from.
#
# The canonical type is drawn on the same axis as the exotics, as the reference
# every exotic bar is read against, and a pooled "All exotics" bar sums every
# non-canonical type. Both currencies of step 2 are kept: an alignment count is
# the dominating type's own depth at the position, and a position count is the
# number of dominated positions, which answers a different question and is not
# recoverable from the other.
#
# The category set, the priority order that resolves overlapping features, the
# promoter window and the palette all come from scripts/utils/localization.R,
# which Stage 10 sources as well — otherwise the two figures could not be read
# against each other. Enhancers are not detected here, because this stage opens
# no BAM; if Stage 10 has already written its enhancer BED, it is reused, and
# otherwise the Enhancer category is simply absent.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/04_dominated_localization.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
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
# exo_type_set "vocabulary" gives every non-canonical type of END_TYPES_REAL a
# bar, plus longclip, which the vocabulary cannot name (its clips vary in
# length) but which is common enough to keep apart from the pooled tail; "top"
# is the historical rule, the exo_top_types most abundant exotics.
exo_type_set <- if (is.null(params$exo_type_set)) "top" else params$exo_type_set
if (!exo_type_set %in% c("vocabulary", "top"))
  stop("exo_type_set must be 'vocabulary' or 'top', not '", exo_type_set, "'")
exo_types <- if (exo_type_set == "vocabulary")
               c(setdiff(END_TYPES_REAL, END_TYPE_CANONICAL), "longclip") else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "dominated_localization.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

if (!fast_mode) {
  suppressPackageStartupMessages(library(ChIPseeker))
}

# ── the dominated positions ───────────────────────────────────────────────────

# step 1 wrote one row per CTSS position, with `call` already reduced to the
# type that dominates it, "below_floor" or "none". Only the owned positions are
# kept here: a position nothing was called at has no type to attribute its
# location to, and the question of this step is what a type's OWN sites look
# like.
dom <- fread(dom_f)
stopifnot(all(c("position_id", "chr", "pos", "strand", "call", "n_dom") %in%
                names(dom)))

owned <- dom[!call %in% c("below_floor", "none"),
             .(position_id, chr, pos, strand,
               type = as.character(call), n_dom)]
stopifnot(nrow(owned) == uniqueN(owned$position_id))
message("Dominated positions: ", nrow(owned), " of ", nrow(dom),
        " called positions, carrying ", sum(owned$n_dom),
        " alignments of their own type")

usual <- END_TYPE_CANONICAL

# ── annotation ────────────────────────────────────────────────────────────────

message("Loading gene annotation: ", gtf_f)
txdb <- make_txdb_from_gff(gtf_f, format = "gtf")

pos_gr <- GRanges(seqnames = owned$chr,
                  ranges   = IRanges(start = owned$pos, width = 1),
                  strand   = owned$strand)

enhancer_gr <- GRanges()
if (!fast_mode && !is.null(enh_f) && nzchar(enh_f) && file.exists(enh_f)) {
  enhancer_gr <- import(enh_f, format = "BED")
  message("Reusing Stage 10 enhancer regions: ", length(enhancer_gr), " from ", enh_f)
} else if (!fast_mode) {
  message("No enhancer BED at ", enh_f,
          "; the Enhancer category will be absent from this figure")
}

if (fast_mode) {
  message("Using fast mode (TxDb-based annotation)")
  tx <- transcripts(txdb)
  promoter_gr <- if (is.null(promoter_cl))
    promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN) else promoter_cl
  utr5_gr <- unlist(fiveUTRsByTranscript(txdb))
  utr3_gr <- unlist(threeUTRsByTranscript(txdb))
  cds_gr <- unlist(cdsBy(txdb, by = "tx"))
  intron_gr <- unlist(intronsByTranscript(txdb))
  owned[, category := annotate_positions_fast(pos_gr, promoter_gr, utr5_gr,
                                              utr3_gr, cds_gr, intron_gr)]
} else {
  message("Annotating with ChIPseeker...")
  owned[, category := annotate_positions_chipseeker(pos_gr, txdb, enhancer_gr,
                                                    promoter_gr = promoter_cl)]
}
stopifnot(!anyNA(owned$category))

# ── bars ──────────────────────────────────────────────────────────────────────

# The bar set mirrors panel A of step 2, so the two figures line up type for
# type: the canonical reference first, then the pooled exotic bar, then the top
# `exo_top_types` exotic types by alignment count with the tail pooled into one
# bar. A type's bar label never enters the annotation — every row is annotated
# from its own position — so a pooled bar is a sum of correctly located rows
# rather than a class of its own.
totals <- owned[type != usual, .(n = sum(n_dom)), by = type][order(-n)]
keep   <- exotic_bar_types(totals$type, top_n, exo_types)
owned[, type_label := fifelse(type == usual, usual,
                       fifelse(type %in% keep, type, OTHER_TYPES_LABEL))]

tally <- function(x) x[, .(alignments = sum(n_dom), positions = .N),
                       by = .(type_label, category)]

per_type <- tally(owned)
overall  <- owned[type != usual, .(alignments = sum(n_dom), positions = .N),
                  by = category][, type_label := ALL_EXOTICS_LABEL][]

res <- rbindlist(list(per_type, overall), use.names = TRUE)

x_levels <- c(if (usual %in% owned$type) usual, ALL_EXOTICS_LABEL, keep,
              if (OTHER_TYPES_LABEL %in% owned$type_label) OTHER_TYPES_LABEL)
cat_order <- intersect(LOC_CATEGORY_ORDER, unique(owned$category))

# Fill the absent (type_label, category) combinations with zeros, so that the
# written table is rectangular and a downstream consumer cannot mistake an
# absent combination for an unmeasured one.
grid <- CJ(type_label = x_levels, category = cat_order, unique = TRUE)
res <- merge(grid, res, by = c("type_label", "category"), all.x = TRUE)
res[is.na(alignments), alignments := 0L]
res[is.na(positions),  positions  := 0L]
res[, pct_alignments := 100 * alignments / sum(alignments), by = type_label]
res[, pct_positions  := 100 * positions  / sum(positions),  by = type_label]
res[, type_label := factor(type_label, levels = x_levels)]
res[, category   := factor(category, levels = cat_order)]
setorder(res, type_label, category)

# Every dominated position and every alignment of a dominating type must be
# drawn exactly once in the per-type bars, and the pooled bar must repeat the
# exotics without adding to them.
per_bars <- setdiff(x_levels, c(usual, ALL_EXOTICS_LABEL))
stopifnot(res[type_label %in% per_bars, sum(positions)] ==
            owned[type != usual, .N])
stopifnot(res[type_label == ALL_EXOTICS_LABEL, sum(alignments)] ==
            owned[type != usual, sum(n_dom)])
stopifnot(res[type_label == usual, sum(positions)] == owned[type == usual, .N])

fwrite(res, out_tsv, sep = "\t")
message("Written: ", out_tsv)
print(res[category == "Promoter", .(type_label, positions, pct_positions)])

# ── figure ────────────────────────────────────────────────────────────────────

plot_dt <- melt(res, id.vars = c("type_label", "category"),
                measure.vars = c("pct_alignments", "pct_positions"),
                variable.name = "measure", value.name = "pct")
plot_dt[, measure := factor(measure,
                            levels = c("pct_alignments", "pct_positions"),
                            labels = c("Alignments of the dominating type",
                                       "Dominated positions"))]

bar_n <- rbind(
  res[, .(measure = "Alignments of the dominating type",
          n = sum(alignments)), by = type_label],
  res[, .(measure = "Dominated positions", n = sum(positions)), by = type_label])
bar_n[, measure := factor(measure, levels = levels(plot_dt$measure))]

p_bar <- ggplot(plot_dt, aes(x = type_label, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(data = bar_n, aes(x = type_label, y = 101, label = fmt_n(n)),
            inherit.aes = FALSE, vjust = 0, size = 2.4, colour = "grey25") +
  facet_wrap(~measure) +
  scale_fill_manual(values = LOC_CATEGORY_COLORS, name = "Category",
                    drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108)) +
  labs(title = "Genomic localization of the positions each 5' end type dominates",
       subtitle = paste0("Dominated CTSS positions only, one owner per position",
                         if (fast_mode) " (fast mode)" else " (ChIPseeker)",
                         "; labels give the total behind each bar"),
       x = "End type", y = "Percentage") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig_f <- file.path(figdir, "22_exo_localization.pdf")
# The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
for (f in c(fig_f, sub("\\.pdf$", ".png", fig_f)))
  ggsave(f, p_bar, width = 10, height = 6, dpi = 200, bg = "white")
message("Written: ", fig_f)

sessionInfo()
