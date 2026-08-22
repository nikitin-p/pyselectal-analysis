#!/usr/bin/env Rscript
# Stage 19 step 1: run the Stage 17 dominance cascade with the CTSS position,
# rather than the distclu cluster, as the unit of analysis.
#
# Stage 17 asks which tag clusters are dominated by one 5' end type. This stage
# asks the same question of individual start sites, because the origin of an
# exotic alignment is a property of the base it starts on and not of the twenty
# base window around it. The four filters, the background shares and the FDR
# correction are reused unchanged from scripts/10_dominant_endtype/dom_helpers.R
# by renaming position_id to cluster_id across the call, so nothing about the
# statistics is re-specified here.
#
# Two consequences of the finer unit are deliberate and are recorded rather
# than hidden. First, Stage 8's per-position floor (dom_min_total) is NOT
# applied: every position that carries a read is classified, because a position
# too shallow to have joined a cluster still has exotic alignments that have to
# be accounted for somewhere. Second, the depth floor is this stage's own
# exo_min_depth rather than Stage 17's dom_min_depth, because dom_min_depth is
# calibrated for a distclu cluster and is far above the depth of a typical
# single position; at that floor virtually every position would be uncallable
# and the figure would describe nothing. Positions still below exo_min_depth
# are carried as an explicit "below_floor" call, which the figure draws instead
# of dropping.
#
# Usage (from project root):
#   Rscript scripts/12_exotic_origins/01_position_dominance.R \
#       [--ctss       results/dominance/fine_typed_ctss.tsv.gz] \
#       [--background results/dominance/background_shares.tsv]  \
#       [--outdir     results/exotic_origins]                   \
#       [--params     config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
ctss_f   <- get_arg("--ctss",       "results/dominance/fine_typed_ctss.tsv.gz")
bg_f     <- get_arg("--background", "results/dominance/background_shares.tsv")
outdir   <- get_arg("--outdir",     "results/exotic_origins")
params_f <- get_arg("--params",     "config/params.yaml")
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))
source(file.path(script_dir, "exo_helpers.R"))

params <- yaml.load_file(params_f)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "position_dominance.tsv.gz")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

# ── input ─────────────────────────────────────────────────────────────────────

if (!file.exists(ctss_f)) {
  stop("Fine-typed CTSS table not found: ", ctss_f,
       "\nRun scripts/10_dominant_endtype/run.sh first.")
}
if (!file.exists(bg_f)) {
  stop("Background shares not found: ", bg_f,
       "\nThis stage reuses Stage 17's genome-wide shares rather than ",
       "recomputing them, so that the binomial test is the same test.")
}

ctss <- fread(ctss_f)
stopifnot(all(c("chr", "pos", "strand", "sample_id", "type", "n") %in%
              names(ctss)))
ctss[, cluster_id := paste(chr, strand, pos, sep = ":")]
message("Loaded ", nrow(ctss), " rows, ", uniqueN(ctss$cluster_id),
        " CTSS positions, ", uniqueN(ctss$sample_id), " samples")

bg_dt <- fread(bg_f)
bg    <- setNames(bg_dt$background, bg_dt$type)

# ── the cascade, at the position ──────────────────────────────────────────────

counts <- ctss[, .(n = sum(n)), by = .(cluster_id, sample_id, type)]

# The dominance threshold is the primary one, dom_thresholds[1]. Stage 17
# evaluates both of its thresholds because its dominance calls are the result;
# here they are only a way to label the origin of an exotic alignment, so the
# secondary threshold would double the figure for no added claim.
threshold <- as.numeric(params$dom_thresholds)[1]

# The depth floor is this stage's own. A params file that predates it (the
# human one, which does not run this stage) falls back to Stage 17's value
# rather than failing, and says so.
min_depth <- if (is.null(params$exo_min_depth)) {
  message("exo_min_depth absent from the params file, falling back to ",
          "dom_min_depth = ", params$dom_min_depth)
  as.integer(params$dom_min_depth)
} else {
  as.integer(params$exo_min_depth)
}
message("Cascade at threshold ", threshold, ", depth floor ", min_depth)

res <- call_dominance(counts, bg,
                      min_depth            = min_depth,
                      threshold            = threshold,
                      min_samples          = as.integer(params$dom_min_samples),
                      min_reads_per_sample = as.integer(params$dom_min_reads_per_sample),
                      fdr                  = as.numeric(params$dom_fdr))
setnames(res, "cluster_id", "position_id")
res[, threshold := threshold]
res[, min_depth := min_depth]

# ── the per-position call ─────────────────────────────────────────────────────

# One call per position: "below_floor", "none", or the dominating fine type.
# The figure turns that call into a category relative to the type of the bar it
# is drawing, which is why no category is assigned here.
res[, call := position_call(res)]

# The call must cover the positions exactly once each: no position
# unclassified, none counted twice. This is asserted rather than trusted,
# because every proportion in the figure is a proportion of this partition.
stopifnot(!anyNA(res$call))
stopifnot(nrow(res) == uniqueN(res$position_id))

# Positions carry chr, strand and pos in their id; splitting them back out
# keeps the table joinable to the rest of the pipeline without re-parsing.
res[, c("chr", "strand", "pos") := tstrsplit(position_id, ":", fixed = TRUE)]
res[, pos := as.integer(pos)]
setcolorder(res, c("position_id", "chr", "pos", "strand", "call", "dom_type",
                   "n_dom", "n_total", "frac"))

summ <- res[, .(n_positions = .N,
                n_alignments = sum(n_total),
                median_depth = as.numeric(median(n_total))), by = call]
setorder(summ, -n_positions)
summ[, pct_positions := 100 * n_positions / sum(n_positions)]
fwrite(summ, file.path(outdir, "call_summary.tsv"), sep = "\t")
message("Written: ", file.path(outdir, "call_summary.tsv"))
print(summ)

# The sentinel is written last, so that an interrupted run leaves no file that
# would make the next run skip this step.
fwrite(res, out_tsv, sep = "\t", compress = "gzip")
message("Written: ", out_tsv)
sessionInfo()
