#!/usr/bin/env Rscript
# Stage 6: per-sample initiator dinucleotide proportion (score-weighted, CTSS ±1 bp).
# Reads subsampled BAMs, extracts 5' CTSS positions, queries the genome for the
# 2-bp window centred on each CTSS (strand-aware), computes raw-count-weighted
# proportions, and saves one PDF + one TSV row per sample.
# Combined TSV is written to {outdir}/all_dinuc_proportions.tsv.
#
# Usage (from project root):
#   Rscript scripts/04_dinuc/01_dinuc_proportion.R \
#       [--bam_dir  results/bam_subsampled] \
#       [--samples  config/samples.tsv]     \
#       [--outdir   results/dinuc]          \
#       [--params   config/params.yaml]     \
#       [--force]
#
# The 'genome' column in samples.tsv accepts either:
#   - a BSgenome alias: hg38 | sacCer3
#   - a path to a FASTA file (.fa / .fasta / .fa.gz) for custom / test genomes

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(BSgenome)
  library(Rsamtools)
  library(data.table)
  library(dplyr)
  library(ggplot2)
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

bam_dir   <- get_arg("--bam_dir", "results/bam_subsampled")
samples_f <- get_arg("--samples",  "config/samples.tsv")
outdir    <- get_arg("--outdir",   "results/dinuc")
params_f  <- get_arg("--params",   "config/params.yaml")
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

# ── genome loading ────────────────────────────────────────────────────────────
# Returns either a BSgenome object or an open FaFile (for custom / test FASTA).
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec))
      stop("FASTA not found: ", genome_spec)
    fai <- paste0(genome_spec, ".fai")
    if (!file.exists(fai)) {
      message("    indexing FASTA: ", genome_spec)
      Rsamtools::indexFa(genome_spec)
    }
    fa <- Rsamtools::FaFile(genome_spec)
    open(fa)
    return(fa)
  }
  pkg <- switch(genome_spec,
    hg38    = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias '", genome_spec,
         "'. Use hg38, sacCer3, or a path to a FASTA file.")
  )
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("BSgenome package not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

# ── seqlevel harmonisation ────────────────────────────────────────────────────
harmonize_seqlevels_to_genome <- function(gr, genome) {
  if (inherits(genome, "FaFile")) {
    fa_seqs <- as.character(seqnames(seqinfo(genome)))
    keep <- intersect(seqlevels(gr), fa_seqs)
    return(keepSeqlevels(gr, keep, pruning.mode = "coarse"))
  }
  # BSgenome path
  try({ seqlevelsStyle(gr) <- seqlevelsStyle(genome) }, silent = TRUE)
  gseqs  <- seqlevels(genome)
  grseqs <- seqlevels(gr)
  if (all(grepl("^chr", gseqs)) && any(!grepl("^chr", grseqs))) {
    new <- ifelse(grseqs %in% c("M", "MT"), grseqs, paste0("chr", grseqs))
    new[new == "chrMT"] <- "chrM"
    suppressWarnings(seqlevels(gr) <- new)
  }
  if ("chrM" %in% gseqs && "chrMT" %in% seqlevels(gr)) {
    sl <- seqlevels(gr); sl[sl == "chrMT"] <- "chrM"
    suppressWarnings(seqlevels(gr) <- sl)
  }
  keep <- intersect(seqlevels(gr), gseqs)
  keepSeqlevels(gr, keep, pruning.mode = "coarse")
}

# ── BAM → CTSS GRanges with raw counts ───────────────────────────────────────
bam_to_ctss_gr <- function(bam_path) {
  param <- ScanBamParam(
    flag = scanBamFlag(
      isSecondaryAlignment     = FALSE,
      isSupplementaryAlignment = FALSE,
      isUnmappedQuery          = FALSE
    )
  )
  ga  <- readGAlignments(bam_path, param = param)
  strand_chr <- as.character(strand(ga))
  # treat unstranded reads (flag 0) as "+" — handles single-end CAGE where
  # all forward-strand reads have flag 0
  pos <- strand_chr %in% c("+", "*")
  neg <- strand_chr == "-"

  gr_list <- list()
  if (any(pos)) {
    gr_list$plus <- GRanges(
      seqnames = seqnames(ga)[pos],
      ranges   = IRanges(start = start(ga)[pos], width = 1),
      strand   = "+"
    )
  }
  if (any(neg)) {
    gr_list$minus <- GRanges(
      seqnames = seqnames(ga)[neg],
      ranges   = IRanges(start = end(ga)[neg], width = 1),
      strand   = "-"
    )
  }
  if (length(gr_list) == 0) return(GRanges())
  gr <- do.call(c, unname(gr_list))

  dt <- data.table(
    chr    = as.character(seqnames(gr)),
    pos    = start(gr),
    strand = as.character(strand(gr))
  )
  agg <- dt[, .(score = .N), by = .(chr, pos, strand)]

  GRanges(
    seqnames = agg$chr,
    ranges   = IRanges(start = agg$pos, width = 1),
    strand   = agg$strand,
    score    = agg$score
  )
}

# ── per-sample processing ─────────────────────────────────────────────────────
process_sample <- function(sample_id, bam_path, genome_spec, outdir, force) {
  out_pdf <- file.path(outdir, paste0(sample_id, "_dinuc_proportion.pdf"))
  out_tsv <- file.path(outdir, paste0(sample_id, "_dinuc_proportion.tsv"))

  if (!force && file.exists(out_pdf) && file.exists(out_tsv)) {
    message("  skip (outputs exist): ", sample_id)
    return(read_tsv(out_tsv, col_types = cols(.default = "c", sum_score = "d", proportion = "d"),
                   show_col_types = FALSE))
  }

  if (!file.exists(bam_path)) {
    warning("BAM not found, skipping: ", bam_path)
    return(NULL)
  }

  message("  processing: ", sample_id)
  genome <- load_genome(genome_spec)
  on.exit(if (inherits(genome, "FaFile")) close(genome), add = TRUE)

  gr <- bam_to_ctss_gr(bam_path)
  if (length(gr) == 0) {
    warning("No alignments extracted from: ", bam_path)
    return(NULL)
  }

  gr  <- harmonize_seqlevels_to_genome(gr, genome)
  if (length(gr) == 0) {
    warning("No CTSS left after seqlevel filtering for: ", sample_id)
    return(NULL)
  }

  win <- promoters(gr, upstream = 1, downstream = 1)
  win <- trim(win)
  keep2 <- width(win) == 2
  win   <- win[keep2]
  mcols(win)$score        <- mcols(gr)$score[keep2]
  mcols(win)$dinucleotide <- as.character(getSeq(genome, win))

  df <- data.frame(
    dinucleotide = mcols(win)$dinucleotide,
    score        = mcols(win)$score,
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(dinucleotide) & nchar(dinucleotide) == 2)

  if (nrow(df) == 0) {
    warning("No valid dinucleotides for: ", sample_id)
    return(NULL)
  }

  props <- df |>
    group_by(dinucleotide) |>
    summarise(sum_score = sum(score, na.rm = TRUE), .groups = "drop") |>
    mutate(
      sample     = sample_id,
      proportion = sum_score / sum(sum_score, na.rm = TRUE)
    ) |>
    arrange(proportion)

  props$dinucleotide <- factor(props$dinucleotide, levels = props$dinucleotide)

  p <- ggplot(props, aes(x = dinucleotide, y = proportion * 100)) +
    geom_bar(stat = "identity", fill = "#4cc9f0", colour = "grey30", linewidth = 0.25) +
    ggtitle(paste0("Initiator dinucleotide proportion\n", sample_id)) +
    xlab("Dinucleotide (CTSS ±1 bp)") +
    ylab("% of total score") +
    coord_flip() +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  pdf(out_pdf, width = 5, height = max(3, nrow(props) * 0.25 + 1.5))
  print(p)
  dev.off()
  message("    saved: ", out_pdf)

  write_tsv(props, out_tsv)

  props
}

# ── load sample manifest ──────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)

if (!all(c("sample_id", "genome") %in% colnames(meta)))
  stop("samples.tsv must have columns: sample_id, genome")
if (nrow(meta) == 0)
  stop("No samples in manifest (are all rows commented out?)")

# Support explicit 'bam' column (used in test manifests) in addition to bam_dir
has_bam_col <- "bam" %in% colnames(meta)

# ── process samples ───────────────────────────────────────────────────────────
message("Stage 6: dinucleotide proportions")
message("  BAM dir : ", bam_dir)
message("  Samples : ", nrow(meta))

results <- lapply(seq_len(nrow(meta)), function(i) {
  sid  <- meta$sample_id[i]
  gnm  <- meta$genome[i]
  bam_path <- if (has_bam_col && !is.na(meta$bam[i])) {
    meta$bam[i]
  } else {
    file.path(bam_dir, paste0(sid, ".bam"))
  }
  process_sample(sid, bam_path, gnm, outdir, force)
})

results <- Filter(Negate(is.null), results)

# ── combined TSV ──────────────────────────────────────────────────────────────
if (length(results) > 0) {
  combined <- bind_rows(results) |>
    select(sample, dinucleotide, sum_score, proportion)
  out_all <- file.path(outdir, "all_dinuc_proportions.tsv")
  write_tsv(combined, out_all)
  message("Combined TSV: ", out_all)
} else {
  message("No samples processed (check BAM paths and genome column in samples.tsv).")
}

sessionInfo()
