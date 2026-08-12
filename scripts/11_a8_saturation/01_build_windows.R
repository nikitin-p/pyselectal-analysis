#!/usr/bin/env Rscript
# Stage 18 step 1: build one compact row per CTSS carrying the initiator
# dinucleotide, the A-8 indicator, the two background A tallies and a distclu
# cluster label.
#
# Sequence is read once and reduced immediately to these few columns, because
# keeping 21 bases for three quarters of a million positions is far more memory
# than the downstream tests need.
#
# Usage:
#   Rscript scripts/11_a8_saturation/01_build_windows.R \
#       [--ctss    results/tss/typed_ctss.tsv.gz] \
#       [--samples config/samples.tsv]            \
#       [--matrix  results/tss/tss_matrix.tsv]    \
#       [--genome  data/sacCer3.fa]               \
#       [--outdir  results/a8]                    \
#       [--params  config/params.yaml]            \
#       [--skip_assert] [--force]

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Rsamtools)
  library(Biostrings); library(yaml); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag)
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "a8_helpers.R"))

ctss_f    <- get_arg("--ctss",    "results/tss/typed_ctss.tsv.gz")
samples_f <- get_arg("--samples", "config/samples.tsv")
matrix_f  <- get_arg("--matrix",  "results/tss/tss_matrix.tsv")
genome_f  <- get_arg("--genome",  NULL)
outdir    <- get_arg("--outdir",  "results/a8")
params_f  <- get_arg("--params",  "config/params.yaml")
skip_assert <- has_flag("--skip_assert")
force       <- has_flag("--force")

params    <- yaml.load_file(params_f)
upstream  <- as.integer(params$a8_upstream)
obs_coord <- as.integer(params$a8_obs_coord)
bg_plat   <- as.integer(params$a8_bg_plateau)
bg_prox   <- as.integer(params$a8_bg_proximal)
min_total <- as.integer(params$a8_min_total)
max_dist  <- as.integer(params$a8_max_dist)

stopifnot(upstream == 20L)   # A8_COORDS is fixed at -20..+1
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_win  <- file.path(outdir, "pair_windows.tsv.gz")
out_prof <- file.path(outdir, "a8_position_profile.tsv")

if (!force && file.exists(out_win)) {
  message("Output exists, skipping (use --force to rerun): ", out_win)
  quit(save = "no", status = 0)
}

TYPE_COLS <- paste0("n_", A8_TYPES)

# ── genome ────────────────────────────────────────────────────────────────────
if (is.null(genome_f)) {
  meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                   show_col_types = FALSE)
  g <- unique(meta$genome[!is.na(meta$genome)])
  if (length(g) != 1) stop("Expected exactly one genome in ", samples_f,
                           ", found: ", paste(g, collapse = ", "))
  genome_f <- g
}
message("Genome: ", genome_f)
if (!file.exists(paste0(genome_f, ".fai"))) indexFa(genome_f)
fa <- FaFile(genome_f); open(fa)
on.exit(close(fa), add = TRUE)
sl <- seqlengths(fa)

# ── pool the per-sample table ─────────────────────────────────────────────────
ctss <- fread(ctss_f)
pos  <- ctss[, c(lapply(.SD, sum), .(n_total = sum(n_total))),
             by = .(chr, pos, strand), .SDcols = TYPE_COLS]
message("  unique CTSS: ", nrow(pos))

# ── clusters, reproducing Stage 8 exactly ─────────────────────────────────────
# Stage 8 applies the min_total floor BEFORE clustering, so this must too.
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))
clustered <- derive_clusters(pos[n_total >= min_total, .(chr, pos, strand, n_total)],
                             max_dist)

if (!skip_assert) {
  stage8 <- unique(fread(matrix_f)$cluster_id)
  mine   <- unique(clustered$cluster_id)
  if (!setequal(stage8, mine)) {
    stop("Cluster set differs from ", matrix_f, ": ",
         length(mine), " here vs ", length(stage8), " there.\n",
         "  Check that a8_min_total (", min_total, ") and a8_max_dist (",
         max_dist, ") equal the --min_total and --max_dist used by ",
         "scripts/05_tss_end_types/02_cluster_annotate.R.")
  }
  message("  cluster set matches Stage 8: ", length(mine), " clusters")
}

pos <- merge(pos, clustered[, .(chr, pos, strand, cluster_id)],
             by = c("chr", "pos", "strand"), all.x = TRUE)
# Positions below the Stage 8 floor were excluded by depth, not by isolation:
# on the real data about 40 % of them lie within max_dist of a clustered
# position and therefore share most of their window with it. Treating each as
# its own bootstrap block would count neighbours as independent and narrow the
# very interval the bootstrap exists to widen, so each is attached to the
# nearest cluster within max_dist and only genuinely isolated positions remain
# singletons.
pos <- attach_unclustered(pos, max_dist)

# ── window extraction ─────────────────────────────────────────────────────────
# Window spans -20..+1 in transcription order. On the plus strand that is
# [pos - 20, pos]; on the minus strand it is [pos, pos + 20]. getSeq on a FaFile
# is strand-aware and always returns 5'->3'.
start_g <- ifelse(pos$strand == "+", pos$pos - upstream, pos$pos)
end_g   <- ifelse(pos$strand == "+", pos$pos,            pos$pos + upstream)
gr <- GRanges(pos$chr, IRanges(start_g, end_g), strand = pos$strand)

inb <- start(gr) >= 1L & end(gr) <= sl[as.character(seqnames(gr))]
if (any(!inb)) message("  dropped out-of-bound windows: ", sum(!inb))
gr <- gr[inb]; pos <- pos[inb]

sq <- getSeq(fa, gr)
stopifnot(all(width(sq) == length(A8_COORDS)))

bases <- as.matrix(sq)                  # n x 21 character
is_A  <- bases == "A"                   # n x 21 logical
i_m1  <- a8_coord_index(-1L); i_p1 <- a8_coord_index(1L)

pos[, dinuc        := paste0(bases[, i_m1], bases[, i_p1])]
pos[, a8           := as.integer(is_A[, a8_coord_index(obs_coord)])]
pos[, n_a_plateau  := as.integer(rowSums(is_A[, a8_coord_index(bg_plat), drop = FALSE]))]
pos[, n_a_proximal := as.integer(rowSums(is_A[, a8_coord_index(bg_prox), drop = FALSE]))]

# Positions containing an ambiguous base at -1 or +1 cannot be assigned a pair.
bad <- !pos$dinuc %in% A8_DINUCS
if (any(bad)) { message("  dropped windows with an ambiguous initiator: ", sum(bad))
                pos <- pos[!bad]; is_A <- is_A[!bad, , drop = FALSE] }

setcolorder(pos, c("chr", "pos", "strand", "dinuc", "cluster_id", "a8",
                   "n_a_plateau", "n_a_proximal", "n_total", TYPE_COLS))
fwrite(pos, out_win, sep = "\t")
message("  wrote ", nrow(pos), " windows -> ", out_win)

# ── the per-position profile ──────────────────────────────────────────────────
prof_for <- function(idx, label) {
  data.table(set = label, coord = A8_COORDS,
             n_A = colSums(is_A[idx, , drop = FALSE]), n_tot = length(idx))
}
prof <- list(prof_for(seq_len(nrow(pos)), "all CTSS"))
for (tc in TYPE_COLS) {
  idx <- which(pos[[tc]] > 0)
  if (length(idx) > 0) prof[[length(prof) + 1L]] <- prof_for(idx, sub("^n_", "", tc))
}
prof <- rbindlist(prof)
prof[, pct_A := 100 * n_A / n_tot]
fwrite(prof, out_prof, sep = "\t")
message("  wrote profile -> ", out_prof)

sessionInfo()
