#!/usr/bin/env Rscript
# Stage 18 step 2: form <end type, dinucleotide> pairs, test A-8 against each
# pair's own upstream background, and bootstrap the effect size over clusters.
#
# Four configurations are emitted. The primary is position-weighted with the
# proximal -10..-2 background (minus -8 itself) and a 2-read membership rule;
# the other three vary exactly one of those three choices, and each is
# BH-corrected on its own.
#
# Usage:
#   Rscript scripts/11_a8_saturation/02_a8_test.R \
#       [--windows results/a8/pair_windows.tsv.gz] \
#       [--outdir  results/a8]                     \
#       [--params  config/params.yaml] [--force]

suppressPackageStartupMessages({ library(data.table); library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag); args[i + 1]
}
has_flag <- function(flag) flag %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "a8_helpers.R"))

win_f    <- get_arg("--windows", "results/a8/pair_windows.tsv.gz")
outdir   <- get_arg("--outdir",  "results/a8")
params_f <- get_arg("--params",  "config/params.yaml")
force    <- has_flag("--force")

params    <- yaml.load_file(params_f)
min_reads <- as.integer(params$a8_min_type_reads)
min_pos   <- as.integer(params$a8_min_positions)
n_boot    <- as.integer(params$a8_n_boot)
seed      <- as.integer(params$a8_seed)
fdr       <- as.numeric(params$a8_fdr)
n_plat    <- length(as.integer(params$a8_bg_plateau))
n_prox    <- length(as.integer(params$a8_bg_proximal))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_tests <- file.path(outdir, "a8_pair_tests.tsv")
out_sizes <- file.path(outdir, "a8_pair_sizes.tsv")
if (!force && file.exists(out_tests)) {
  message("Output exists, skipping (use --force to rerun): ", out_tests)
  quit(save = "no", status = 0)
}

win <- fread(win_f)

# config, background column, number of background positions, weighting,
# membership threshold, and whether to bootstrap.
CONFIGS <- list(
  primary     = list(bg = "n_a_proximal", k = n_prox, w = "positions", m = min_reads, boot = TRUE),
  bg_plateau  = list(bg = "n_a_plateau",  k = n_plat, w = "positions", m = min_reads, boot = FALSE),
  w_reads     = list(bg = "n_a_proximal", k = n_prox, w = "reads",     m = min_reads, boot = FALSE),
  member_1    = list(bg = "n_a_proximal", k = n_prox, w = "positions", m = 1L,        boot = FALSE)
)

grid <- as.data.table(expand.grid(end_type = A8_TYPES, dinuc = A8_DINUCS,
                                  stringsAsFactors = FALSE))
grid[, pair := paste0("<", end_type, ", ", dinuc, ">")]
grid[, structurally_empty := is_structurally_empty(end_type, dinuc)]

# The structural-emptiness assertion. A populated cell here means the window
# extraction got a strand or an offset wrong, so this is a hard failure.
for (i in which(grid$structurally_empty)) {
  tc <- paste0("n_", grid$end_type[i])
  n  <- win[dinuc == grid$dinuc[i] & get(tc) > 0, .N]
  if (n > 0L) {
    stop("Structurally empty pair ", grid$pair[i], " carries ", n, " positions.\n",
         "  Under --alignEndsType Local a ", grid$end_type[i], " read cannot occur ",
         "on a dinucleotide starting with ", A8_CLIP_BASE[grid$end_type[i]], ".\n",
         "  This indicates a strand or offset error in 01_build_windows.R.")
  }
}
message("  structural-emptiness assertion passed for ",
        sum(grid$structurally_empty), " pairs")

all_tests <- list(); all_sizes <- list()

for (cfg_name in names(CONFIGS)) {
  cfg <- CONFIGS[[cfg_name]]
  message("  config: ", cfg_name)
  rows <- vector("list", nrow(grid)); szs <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    et <- grid$end_type[i]; dn <- grid$dinuc[i]; tc <- paste0("n_", et)
    sub <- win[dinuc == dn & get(tc) >= cfg$m]
    wgt <- if (cfg$w == "reads") sub[[tc]] else rep(1L, nrow(sub))

    n_positions <- nrow(sub)
    n_reads     <- if (n_positions > 0L) sum(sub[[tc]]) else 0L
    tested <- !grid$structurally_empty[i] && n_positions >= min_pos
    reason <- if (grid$structurally_empty[i]) "structurally empty"
              else if (n_positions < min_pos) paste0("below the ", min_pos, "-position floor")
              else "tested"

    szs[[i]] <- data.table(config = cfg_name, end_type = et, dinuc = dn,
                           pair = grid$pair[i],
                           structurally_empty = grid$structurally_empty[i],
                           n_positions = n_positions, n_reads = n_reads,
                           tested = tested, reason = reason)
    if (!tested) next

    n_obs  <- sum(wgt)
    n_A8   <- sum(wgt[sub$a8 == 1L])
    n_A_bg <- sum(wgt * sub[[cfg$bg]])
    n_bg   <- cfg$k * n_obs
    p_obs  <- n_A8 / n_obs
    p_exp  <- n_A_bg / n_bg

    # Cluster-level aggregation, which is what the bootstrap resamples.
    cl <- data.table(cluster_id = sub$cluster_id, w = wgt,
                     a8 = wgt * sub$a8, bga = wgt * sub[[cfg$bg]]
                     )[, .(W = sum(w), A8 = sum(a8), BGA = sum(bga)), by = cluster_id]

    if (isTRUE(cfg$boot)) {
      b  <- bootstrap_log2fc(cl$W, cl$A8, cl$BGA, cfg$k, n_boot, seed)
      ci <- unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))
    } else {
      ci <- c(NA_real_, NA_real_)
    }

    rows[[i]] <- data.table(
      config = cfg_name, end_type = et, dinuc = dn, pair = grid$pair[i],
      n_positions = n_positions, n_clusters = nrow(cl),
      n_obs = n_obs, n_A8 = n_A8, p_obs = p_obs,
      n_bg = n_bg, n_A_bg = n_A_bg, p_exp = p_exp,
      log2fc = a8_log2fc(n_A8, n_A_bg, cfg$k),
      neg_log10_p = a8_pvalue(n_A8, n_obs, p_exp),
      boot_lo = ci[1], boot_hi = ci[2])
  }

  tt <- rbindlist(rows)
  if (nrow(tt) > 0L) {
    tt[, neg_log10_q := bh_neg_log10(neg_log10_p)]
    tt[, significant := neg_log10_q > -log10(fdr)]
    tt[, boot_robust := !is.na(boot_lo) & (boot_lo > 0 | boot_hi < 0)]
    message("    tested ", nrow(tt), " pairs, ", sum(tt$significant), " significant")
  }
  all_tests[[cfg_name]] <- tt
  all_sizes[[cfg_name]] <- rbindlist(szs)
}

fwrite(rbindlist(all_tests), out_tests, sep = "\t")
fwrite(rbindlist(all_sizes), out_sizes, sep = "\t")
message("  wrote ", out_tests, " and ", out_sizes)

sessionInfo()
