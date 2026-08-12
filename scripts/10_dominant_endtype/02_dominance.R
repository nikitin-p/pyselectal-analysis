#!/usr/bin/env Rscript
# Stage 17 step 2: re-derive the Stage 8 clusters from the fine-typed CTSS,
# assert they are identical, and call per-cluster dominance.
#
# Re-typing the reads changes no position and no per-position total, so the
# clustering MUST reproduce Stage 8's cluster set. The assertion below is what
# licenses every comparison in this stage; if it fails, the run aborts rather
# than proceeding on a silently divergent cluster set.
#
# Usage (from project root):
#   Rscript scripts/10_dominant_endtype/02_dominance.R \
#       [--ctss   results/dominance/fine_typed_ctss.tsv.gz] \
#       [--matrix results/tss/tss_matrix.tsv]               \
#       [--outdir results/dominance]                        \
#       [--params config/params.yaml]                       \
#       [--max_dist N] [--skip_assert] [--force]
#
# --max_dist overrides the dom_max_dist config key, which must equal Stage 8's
# clustering distance for the cluster-identity assertion below to hold.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

ctss_f   <- get_arg("--ctss",   "results/dominance/fine_typed_ctss.tsv.gz")
matrix_f <- get_arg("--matrix", "results/tss/tss_matrix.tsv")
outdir   <- get_arg("--outdir", "results/dominance")
params_f <- get_arg("--params", "config/params.yaml")
skip_assert <- "--skip_assert" %in% args   # for the synthetic fixture only
force    <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "dom_helpers.R"))

params <- yaml.load_file(params_f)
# The clustering distance is as load-bearing as the per-position floor, so it
# lives in config; the CLI flag remains available as an override.
max_dist <- as.integer(get_arg("--max_dist", as.character(params$dom_max_dist)))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_calls <- file.path(outdir, "dominance_calls.tsv")
if (!force && file.exists(out_calls)) {
  message("Output exists, skipping (use --force to rerun): ", out_calls)
  quit(save = "no", status = 0)
}

message("Stage 17 step 2: clustering + dominance cascade")
ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
message("  fine CTSS rows: ", nrow(ctss),
        "   distinct types: ", uniqueN(ctss$type))

# ── re-derive clusters ────────────────────────────────────────────────────────
min_total  <- as.integer(params$dom_min_total)
total_ctss <- ctss[, .(n_total = sum(n)), by = .(chr, pos, strand)]
total_ctss <- total_ctss[n_total >= min_total]
message("  CTSS positions (n_total >= ", min_total, "): ", nrow(total_ctss))

clustered <- derive_clusters(total_ctss, max_dist)
message("  tag clusters: ", uniqueN(clustered$cluster_id))

# ── assert identity with Stage 8 ──────────────────────────────────────────────
if (!skip_assert) {
  if (!file.exists(matrix_f)) stop("Stage 8 matrix not found: ", matrix_f)
  ref <- unique(fread(matrix_f, select = c("cluster_id", "chr", "tss_pos",
                                           "strand")))
  dom <- clustered[clustered[, .I[which.max(n_total)][1L],
                             by = cluster_id]$V1,
                   .(cluster_id, chr, tss_pos = pos, strand)]
  setkey(ref, cluster_id); setkey(dom, cluster_id)
  ok <- nrow(ref) == nrow(dom) &&
        setequal(ref$cluster_id, dom$cluster_id) &&
        identical(ref[dom$cluster_id, tss_pos], dom$tss_pos)
  if (!ok) {
    stop("Re-derived clusters differ from Stage 8 (", nrow(dom), " vs ",
         nrow(ref), " clusters). Every comparison in this stage assumes the ",
         "two cluster sets are identical. Check that dom_min_total (",
         min_total, ") equals the --min_total used by ",
         "scripts/05_tss_end_types/02_cluster_annotate.R.")
  }
  message("  cluster identity with Stage 8: OK (", nrow(ref), " clusters)")
}

# ── per-cluster per-type counts ───────────────────────────────────────────────
pos_map <- clustered[, .(chr, pos, strand, cluster_id)]
joined  <- merge(ctss, pos_map, by = c("chr", "pos", "strand"))
counts  <- joined[, .(n = sum(n)), by = .(cluster_id, sample_id, type)]

pooled <- counts[, .(n = sum(n)), by = .(cluster_id, type)]
setorder(pooled, cluster_id, -n)
fwrite(pooled, file.path(outdir, "cluster_type_counts.tsv"), sep = "\t")

# ── genome-wide background shares ─────────────────────────────────────────────
bg_dt <- ctss[, .(n = sum(n)), by = type]
bg    <- setNames(bg_dt$n / sum(bg_dt$n), bg_dt$type)
fwrite(data.table(type = names(bg), background = as.numeric(bg))[order(-background)],
       file.path(outdir, "background_shares.tsv"), sep = "\t")

# ── cascade at every threshold ────────────────────────────────────────────────
thresholds <- as.numeric(params$dom_thresholds)
all_calls  <- list()
all_aggr   <- list()
cascade    <- list()

for (thr in thresholds) {
  res <- call_dominance(counts, bg,
                        min_depth            = as.integer(params$dom_min_depth),
                        threshold            = thr,
                        min_samples          = as.integer(params$dom_min_samples),
                        min_reads_per_sample = as.integer(params$dom_min_reads_per_sample),
                        fdr                  = as.numeric(params$dom_fdr))
  res[, threshold := thr]
  all_calls[[as.character(thr)]] <- res

  # Supplementary aggregate test, run at the same thresholds and with the same
  # filters. It is written to its own table and never merged into the per-type
  # calls, which are left exactly as they were.
  agg <- call_aggregate_unusual(counts, bg,
                        min_depth            = as.integer(params$dom_min_depth),
                        threshold            = thr,
                        min_samples          = as.integer(params$dom_min_samples),
                        min_reads_per_sample = as.integer(params$dom_min_reads_per_sample),
                        fdr                  = as.numeric(params$dom_fdr))
  agg[, threshold := thr]
  all_aggr[[as.character(thr)]] <- agg
  message("  threshold ", thr, ": ",
          agg[pass_depth & pass_frac, .N],
          " clusters exceed the threshold on the pooled unusual share, ",
          agg[(aggregate_unusual), .N], " pass the full aggregate test")

  cascade[[as.character(thr)]] <- data.table(
    threshold  = thr,
    step       = c("all_clusters", "pass_depth", "+pass_frac", "+pass_binom",
                   "+pass_reprod"),
    n_clusters = c(nrow(res),
                   res[(pass_depth), .N],
                   res[pass_depth & pass_frac, .N],
                   res[pass_depth & pass_frac & pass_binom, .N],
                   res[(dominated), .N]))
  message("  threshold ", thr, ": ", res[(dominated), .N], " dominated clusters")
}

calls <- rbindlist(all_calls)
aggr  <- rbindlist(all_aggr)
setcolorder(aggr, c("cluster_id", "threshold"))

# The idempotency sentinel is out_calls, so it is written LAST: an interrupted
# or failing run must not leave a sentinel behind that makes the next run skip a
# step whose remaining outputs were never produced.
fwrite(rbindlist(cascade), file.path(outdir, "filter_cascade.tsv"), sep = "\t")
fwrite(aggr, file.path(outdir, "aggregate_unusual.tsv"), sep = "\t")
fwrite(calls, out_calls, sep = "\t")

message("Written: ", out_calls)
message("Written: ", file.path(outdir, "aggregate_unusual.tsv"))
for (thr in thresholds) {
  message("  dominant types at ", thr, ":")
  print(calls[threshold == thr & dominated, .N, by = dom_type][order(-N)])
}

sessionInfo()
