#!/usr/bin/env Rscript
# Stage 8 step 2: cluster typed CTSS into tag clusters (distclu, maxDist=20 as per
# CAGEr default), identify dominant TSS per cluster, extract initiator dinucleotide
# [TSS-1, TSS] strand-aware, and build a per-cluster per-sample count matrix.
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/02_cluster_annotate.R \
#       [--ctss    results/tss/typed_ctss.tsv.gz] \
#       [--samples config/samples.tsv]            \
#       [--outdir  results/tss]                   \
#       [--params  config/params.yaml]            \
#       [--max_dist 20]                            \
#       [--min_total 5]                            \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(readr)
  library(yaml)
})

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

ctss_f    <- get_arg("--ctss",      "results/tss/typed_ctss.tsv.gz")
samples_f <- get_arg("--samples",   "config/samples.tsv")
outdir    <- get_arg("--outdir",    "results/tss")
params_f  <- get_arg("--params",    "config/params.yaml")
max_dist  <- as.integer(get_arg("--max_dist", "20"))  # CAGEr default
min_total <- as.integer(get_arg("--min_total", "5"))  # min reads across all samples
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

out_tsv <- file.path(outdir, "tss_matrix.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

# ── load data ─────────────────────────────────────────────────────────────────
message("Stage 8 step 2: cluster + annotate")
message("  loading typed CTSS: ", ctss_f)

ctss <- fread(ctss_f)
message("  rows: ", nrow(ctss))

# ── genome loading (same helper as Stage 6) ───────────────────────────────────
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec)) stop("FASTA not found: ", genome_spec)
    fai <- paste0(genome_spec, ".fai")
    if (!file.exists(fai)) {
      message("  indexing FASTA: ", genome_spec)
      Rsamtools::indexFa(genome_spec)
    }
    fa <- Rsamtools::FaFile(genome_spec)
    open(fa)
    return(fa)
  }
  pkg <- switch(genome_spec,
    hg38    = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias: ", genome_spec)
  )
  if (!requireNamespace(pkg, quietly = TRUE)) stop("BSgenome not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

harmonize_seqlevels_to_genome <- function(gr, genome) {
  if (inherits(genome, "FaFile")) {
    fa_si   <- seqinfo(genome)
    fa_seqs <- as.character(seqnames(fa_si))
    keep    <- intersect(seqlevels(gr), fa_seqs)
    gr      <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
    seqlengths(gr) <- seqlengths(fa_si)[seqlevels(gr)]
    return(gr)
  }
  try({ seqlevelsStyle(gr) <- seqlevelsStyle(genome) }, silent = TRUE)
  gseqs  <- seqlevels(genome)
  grseqs <- seqlevels(gr)
  if (all(grepl("^chr", gseqs)) && any(!grepl("^chr", grseqs))) {
    new <- ifelse(grseqs %in% c("M", "MT"), grseqs, paste0("chr", grseqs))
    new[new == "chrMT"] <- "chrM"
    suppressWarnings(seqlevels(gr) <- new)
  }
  keep <- intersect(seqlevels(gr), gseqs)
  keepSeqlevels(gr, keep, pruning.mode = "coarse")
}

# Detect genome from samples.tsv
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)
genome_spec <- unique(na.omit(meta$genome))
if (length(genome_spec) > 1) {
  warning("Multiple genome specs found; using the first: ", genome_spec[1])
  genome_spec <- genome_spec[1]
}
genome <- load_genome(genome_spec)
on.exit(if (inherits(genome, "FaFile")) close(genome), add = TRUE)

# ── sum total counts across samples ──────────────────────────────────────────
total_ctss <- ctss[, .(n_total = sum(n_total)), by = .(chr, pos, strand)]

# Apply minimum total count filter
total_ctss <- total_ctss[n_total >= min_total]
message("  CTSS positions (n_total >= ", min_total, "): ", nrow(total_ctss))

# ── distclu clustering ────────────────────────────────────────────────────────
# CAGEr default: merge adjacent positions within maxDist bp on same chr+strand.
# Reference: CAGEr::clusterCTSS() method="distclu", maxDist=20
setorder(total_ctss, chr, strand, pos)

total_ctss[, cluster_local := {
  if (.N == 1L) return(1L)
  gap <- c(Inf, diff(pos))
  cumsum(gap > max_dist) + 1L
}, by = .(chr, strand)]

total_ctss[, cluster_id := paste(chr, strand, cluster_local, sep = ":")]

# Dominant TSS = position with highest total coverage in each cluster
dominant <- total_ctss[, .SD[which.max(n_total)][1L], by = cluster_id]
dominant <- dominant[, .(cluster_id, chr, pos, strand)]
setnames(dominant, "pos", "tss_pos")

message("  tag clusters: ", nrow(dominant))

# ── extract initiator dinucleotide ────────────────────────────────────────────
# [TSS-1, TSS]: promoters(gr, upstream=1, downstream=1) strand-aware,
# then getSeq gives the sequence in the TSS's strand orientation.
gr_dom <- GRanges(
  seqnames = dominant$chr,
  ranges   = IRanges(start = dominant$tss_pos, width = 1),
  strand   = dominant$strand
)
gr_dom <- harmonize_seqlevels_to_genome(gr_dom, genome)

win  <- promoters(gr_dom, upstream = 1, downstream = 1)
win  <- trim(win)
keep <- which(width(win) == 2)

initiator <- rep(NA_character_, length(gr_dom))
if (length(keep) > 0) {
  initiator[keep] <- as.character(getSeq(genome, win[keep]))
}

dominant$initiator <- initiator

# ── build per-position cluster membership map ─────────────────────────────────
# Map every (chr, pos, strand) in total_ctss to its cluster_id
pos_to_cluster <- total_ctss[, .(chr, pos, strand, cluster_id)]

# ── sum typed counts per cluster per sample ───────────────────────────────────
# Keep only CTSS positions that passed the min_total filter
ctss_filtered <- merge(ctss, pos_to_cluster, by = c("chr", "pos", "strand"))

tss_matrix <- ctss_filtered[, .(
  n_1Sg   = sum(n_1Sg),
  n_2Sg   = sum(n_2Sg),
  n_M     = sum(n_M),
  n_other = sum(n_other),
  n_total = sum(n_total)
), by = .(cluster_id, sample_id)]

# Compute percentages
tss_matrix[n_total > 0, pct_1Sg   := 100 * n_1Sg   / n_total]
tss_matrix[n_total > 0, pct_2Sg   := 100 * n_2Sg   / n_total]
tss_matrix[n_total > 0, pct_M     := 100 * n_M     / n_total]
tss_matrix[n_total > 0, pct_other := 100 * n_other / n_total]
tss_matrix[n_total == 0, c("pct_1Sg","pct_2Sg","pct_M","pct_other") := NA_real_]

# Attach dominant TSS annotation
tss_matrix <- merge(tss_matrix, dominant, by = "cluster_id")

setcolorder(tss_matrix, c(
  "cluster_id", "chr", "tss_pos", "strand", "initiator",
  "sample_id",
  "n_1Sg", "n_2Sg", "n_M", "n_other", "n_total",
  "pct_1Sg", "pct_2Sg", "pct_M", "pct_other"
))
setorder(tss_matrix, cluster_id, sample_id)

# ── write output ──────────────────────────────────────────────────────────────
fwrite(tss_matrix, out_tsv, sep = "\t")
message("Written: ", out_tsv, "  (", uniqueN(tss_matrix$cluster_id), " clusters, ",
        uniqueN(tss_matrix$sample_id), " samples)")

sessionInfo()
