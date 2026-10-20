#!/usr/bin/env Rscript
# Stage 18 step 4: the same A-8 test as step 2, but with the initiator
# dinucleotides pooled into their four purine/pyrimidine classes.
#
# Step 2 tests one cell per <end type, dinucleotide> pair, which is the finest
# resolution the data support and which leaves several pairs below the
# 200-position floor. Pooling the 16 dinucleotides into YR, RY, RR and YY tests
# the same question at the resolution the biology is usually stated in: the
# canonical initiators are exactly the YR class, so a G-clipped type enriched
# for A-8 on YR is the claim the pair figure makes one dinucleotide at a time.
#
# The pooling is done over the WINDOW ROWS, not over the pair results: each cell
# re-runs the whole test on the positions of all its member dinucleotides, so
# its background, its bootstrap and its p-value are computed from the pooled
# tallies rather than averaged from the per-pair ones. Pooling the results would
# have been wrong, because the pair-level log2FCs carry different denominators.
#
# Structural emptiness is partial rather than total at this resolution. Under
# STAR --alignEndsType Local a clipped type cannot occur on a dinucleotide
# starting with its own clip base, which removes some MEMBERS of a class but
# never a whole class, so every cell here is testable and the table records
# which dinucleotides actually contributed.
#
# Only the primary configuration is run. The three sensitivity configurations
# vary choices whose effect is already shown, at the pair level, by
# results/figures/19_a8_sensitivity.pdf; pooling does not interact with any of
# them.
#
# Usage:
#   Rscript scripts/11_a8_saturation/04_class_test.R \
#       [--windows results/a8/pair_windows.tsv.gz] \
#       [--outdir  results/a8]                     \
#       [--figdir  results/figures]                \
#       [--params  config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(yaml)
})

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
figdir   <- get_arg("--figdir",  "results/figures")
params_f <- get_arg("--params",  "config/params.yaml")
force    <- has_flag("--force")

params    <- yaml.load_file(params_f)
min_reads <- as.integer(params$a8_min_type_reads)
min_pos   <- as.integer(params$a8_min_positions)
n_boot    <- as.integer(params$a8_n_boot)
seed      <- as.integer(params$a8_seed)
fdr       <- as.numeric(params$a8_fdr)
bg_prox   <- as.integer(params$a8_bg_proximal)
n_prox    <- length(bg_prox)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tests <- file.path(outdir, "a8_class_tests.tsv")
f_q       <- file.path(figdir, "19_a8_bubble_by_q_class.pdf")
if (!force && file.exists(out_tests) && file.exists(f_q)) {
  message("Output exists, skipping (use --force to rerun): ", out_tests)
  quit(save = "no", status = 0)
}

# ── the four classes ──────────────────────────────────────────────────────────
# Y = C or T, R = A or G, and the class of a dinucleotide is the class of its
# first base followed by that of its second. This is the same rule Stage 8
# step 9 uses; YR is the canonical initiator class (CA, CG, TA, TG).
CLASSES <- c("YR", "RY", "RR", "YY")
dinuc_class <- function(d)
  paste0(ifelse(substr(d, 1L, 1L) %in% c("C", "T"), "Y", "R"),
         ifelse(substr(d, 2L, 2L) %in% c("C", "T"), "Y", "R"))

win <- fread(win_f)
win[, dclass := dinuc_class(dinuc)]

grid <- as.data.table(expand.grid(end_type = A8_TYPES, dclass = CLASSES,
                                  stringsAsFactors = FALSE))
grid[, cell := paste0("<", end_type, ", ", dclass, ">")]

rows <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  et <- grid$end_type[i]; cl <- grid$dclass[i]; tc <- paste0("n_", et)
  sub <- win[dclass == cl & get(tc) >= min_reads]

  # The class members that are structurally reachable for this end type, and
  # the ones that actually carried a qualifying position.
  members  <- A8_DINUCS[dinuc_class(A8_DINUCS) == cl]
  reachable <- members[!is_structurally_empty(et, members)]
  observed  <- sort(unique(sub$dinuc))

  n_positions <- nrow(sub)
  if (n_positions < min_pos) {
    message("  skipping ", grid$cell[i], ": ", n_positions,
            " positions, below the ", min_pos, "-position floor")
    next
  }

  wgt    <- rep(1L, n_positions)          # position-weighted, the primary rule
  n_obs  <- sum(wgt)
  n_A8   <- sum(wgt[sub$a8 == 1L])
  n_A_bg <- sum(wgt * sub$n_a_proximal)
  n_bg   <- n_prox * n_obs

  cl_tot <- data.table(cluster_id = sub$cluster_id, w = wgt,
                       a8 = wgt * sub$a8, bga = wgt * sub$n_a_proximal
                       )[, .(W = sum(w), A8 = sum(a8), BGA = sum(bga)), by = cluster_id]

  b  <- bootstrap_log2fc(cl_tot$W, cl_tot$A8, cl_tot$BGA, n_prox, n_boot, seed)
  ci <- unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))

  rows[[i]] <- data.table(
    config = "primary", end_type = et, dinuc_class = cl, cell = grid$cell[i],
    n_dinucs_reachable = length(reachable), n_dinucs_observed = length(observed),
    dinucs_observed = paste(observed, collapse = ","),
    n_positions = n_positions, n_clusters = nrow(cl_tot),
    n_obs = n_obs, n_A8 = n_A8, p_obs = n_A8 / n_obs,
    n_bg = n_bg, n_A_bg = n_A_bg, p_exp = n_A_bg / n_bg,
    log2fc = a8_log2fc(n_A8, n_A_bg, n_prox),
    neg_log10_p = a8_pvalue(n_A8, n_obs, n_A_bg / n_bg),
    boot_lo = ci[1], boot_hi = ci[2])
}

res <- rbindlist(rows)
if (nrow(res) == 0L) stop("No testable <end type, dinucleotide class> cells")

# BH over this family alone. The class cells are not a subset of the pair
# family: they are a different set of hypotheses on the same reads, so mixing
# the two into one correction would count the same evidence twice.
res[, neg_log10_q := bh_neg_log10(neg_log10_p)]
res[, significant := neg_log10_q > -log10(fdr)]
res[, boot_robust := !is.na(boot_lo) & (boot_lo > 0 | boot_hi < 0)]
message("  tested ", nrow(res), " cells, ", sum(res$significant), " significant")

fwrite(res, out_tests, sep = "\t")
message("  wrote ", out_tests)

# ── the bubble plots ──────────────────────────────────────────────────────────
# Same encoding as 03_plot.R: area for |log2FC|, a diverging fill for its sign,
# a sqrt x-axis for -log10(q) and a thick outline where the bootstrap CI
# excludes zero.
YR_RED   <- "#e63946"   # project convention for canonical initiators
BLUE     <- "#2166ac"
OTHER_GY <- "grey20"

res[, label := paste0(end_type, "  |  ", dinuc_class)]
res[, abs_fc := abs(log2fc)]

bubble <- function(dt, order_col, fname, subtitle) {
  d <- copy(dt)
  setorderv(d, order_col, -1L)
  d[, label := factor(label, levels = rev(d$label))]
  lv   <- levels(d$label)
  cols <- ifelse(d$dinuc_class[match(lv, d$label)] == "YR", YR_RED, OTHER_GY)

  p <- ggplot(d, aes(x = neg_log10_q, y = label)) +
    geom_vline(xintercept = -log10(fdr), linetype = "dotted", colour = "grey50") +
    geom_point(aes(size = abs_fc, fill = log2fc,
                   stroke = fifelse(boot_robust, 0.7, 0.2)),
               shape = 21, colour = "grey20") +
    scale_size_area(max_size = 9, name = expression(group("|", log[2]*FC, "|"))) +
    scale_fill_gradient2(low = BLUE, mid = "grey95", high = YR_RED,
                         midpoint = 0, name = expression(log[2]*FC)) +
    scale_x_sqrt(name = expression(-log[10]~(adjusted~italic(p)))) +
    labs(y = NULL,
         title = "A-8 enrichment by <5' end type, initiator dinucleotide class>",
         subtitle = subtitle) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(colour = cols),
          panel.grid.major.y = element_line(linewidth = 0.2))

  # The PNG twin is what the Quarto report embeds, since HTML cannot show a PDF.
  for (f in c(fname, sub("\\.pdf$", ".png", fname)))
    ggsave(file.path(figdir, f), p, width = 8,
           height = max(4, 0.16 * nrow(d) + 2), limitsize = FALSE, dpi = 200, bg = "white")
}

bubble(res, "neg_log10_q", "19_a8_bubble_by_q_class.pdf",
       paste0("dinucleotides pooled into YR / RY / RR / YY; position-weighted, ",
              "background ", min(bg_prox), "..", max(bg_prox), " excluding -8,\n",
              "membership >= ", min_reads, " alignments; fill: red = enriched, ",
              "blue = depleted; thick outline = bootstrap CI excludes 0\n",
              "labels: red = YR, the canonical initiator class; a clipped end ",
              "type reaches only the class\nmembers that do not start with its ",
              "clip base"))
bubble(res, "log2fc", "19_a8_bubble_by_fc_class.pdf",
       paste("same data ordered by log2FC, which separates the enriched and",
             "depleted ends of the spectrum"))
message("  wrote ", f_q)

sessionInfo()
