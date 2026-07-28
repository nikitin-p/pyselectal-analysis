#!/usr/bin/env Rscript
# Stage 8 step 6: Compare sequence motifs around TSS by dominant end type
#
# For each end type (1Sg, 2Sgg, 3Sggg, 1Sa, 1St, 1Sc, M), selects TSS where that
# type is enriched (top percentile), filters to promoter-localized + YR initiators,
# and generates sequence logos and nucleotide frequency comparisons.
#
# Usage:
#   Rscript scripts/05_tss_end_types/06_motif_analysis.R \
#       [--tss_matrix results/tss/tss_matrix.tsv] \
#       [--params     config/params.yaml]         \
#       [--gtf        data/sacCer3.ensGene.gtf]   \
#       [--outdir     results/motifs]             \
#       [--figdir     results/figures]            \
#       [--min_reads  10]                         \
#       [--top_pct    10]                         \
#       [--upstream   50]                         \
#       [--downstream 10]                         \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)

# Flags are registered as they are read, so validate_args() below can reject
# anything unrecognized instead of silently ignoring it.
.val_flags  <- character(0)   # flags that consume a following value
.bool_flags <- character(0)   # standalone switches

get_arg <- function(flag, default = NULL) {
  .val_flags <<- union(.val_flags, flag)
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

has_flag <- function(flag) {
  .bool_flags <<- union(.bool_flags, flag)
  flag %in% args
}

validate_args <- function() {
  unknown <- character(0)
  i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% .val_flags) { i <- i + 2L; next }
    if (args[i] %in% .bool_flags) { i <- i + 1L; next }
    unknown <- c(unknown, args[i])
    i <- i + 1L
  }
  if (length(unknown) > 0) {
    stop("Unrecognized argument(s): ", paste(unknown, collapse = " "),
         "\n  Known flags: ", paste(sort(c(.val_flags, .bool_flags)), collapse = " "),
         call. = FALSE)
  }
}

tss_f      <- get_arg("--tss_matrix", "results/tss/tss_matrix.tsv")
params_f   <- get_arg("--params",     "config/params.yaml")
gtf_f      <- get_arg("--gtf",        "data/sacCer3.ensGene.gtf")
outdir     <- get_arg("--outdir",     "results/motifs")
figdir     <- get_arg("--figdir",     "results/figures")
min_reads  <- as.integer(get_arg("--min_reads", "10"))
top_pct    <- as.numeric(get_arg("--top_pct", "10"))
upstream   <- as.integer(get_arg("--upstream", "50"))
downstream <- as.integer(get_arg("--downstream", "10"))
force      <- has_flag("--force")

validate_args()

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

out_logo <- file.path(figdir, "11_motif_logo_comparison.pdf")
if (!force && file.exists(out_logo)) {
  message("Output exists, skipping (use --force to rerun): ", out_logo)
  quit(save = "no", status = 0)
}

if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  stop("ggseqlogo not installed. Run: install.packages('ggseqlogo')")
}
library(ggseqlogo)

# ── genome loading ────────────────────────────────────────────────────────────
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec)) stop("FASTA not found: ", genome_spec)
    fai <- paste0(genome_spec, ".fai")
    if (!file.exists(fai)) Rsamtools::indexFa(genome_spec)
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

harmonize_seqlevels <- function(gr, genome) {
  if (inherits(genome, "FaFile")) {
    genome_seqinfo <- seqinfo(genome)
  } else {
    genome_seqinfo <- seqinfo(genome)
  }
  common <- intersect(seqlevels(gr), seqlevels(genome_seqinfo))
  if (length(common) == 0) stop("No common seqlevels between GRanges and genome")
  gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
  seqinfo(gr) <- genome_seqinfo[seqlevels(gr)]
  gr
}

extract_sequences <- function(dt, genome, up, down) {
  if (nrow(dt) == 0) return(DNAStringSet())

  gr <- GRanges(
    seqnames = dt$chr,
    ranges   = IRanges(start = dt$tss_pos, width = 1),
    strand   = dt$strand
  )
  gr <- harmonize_seqlevels(gr, genome)

  gr_ext <- gr
  start(gr_ext) <- ifelse(strand(gr) == "+", start(gr) - up, start(gr) - down)
  end(gr_ext) <- ifelse(strand(gr) == "+", end(gr) + down, end(gr) + up)
  gr_ext <- trim(gr_ext)

  expected_width <- up + 1 + down
  keep <- width(gr_ext) == expected_width
  gr_ext <- gr_ext[keep]

  if (length(gr_ext) == 0) return(DNAStringSet())

  # getSeq() already returns the reverse complement for minus-strand ranges,
  # so the sequences come back 5'->3' along the transcript. Do NOT revcomp again.
  seqs <- getSeq(genome, gr_ext)
  seqs
}

compute_pwm <- function(seqs) {
  if (length(seqs) == 0) return(NULL)
  mat <- consensusMatrix(seqs, as.prob = TRUE)[1:4, ]
  rownames(mat) <- c("A", "C", "G", "T")
  mat
}

# ── load TSS matrix ───────────────────────────────────────────────────────────
message("Loading TSS matrix: ", tss_f)
tss <- fread(tss_f)
message("  rows: ", nrow(tss), " (", uniqueN(tss$cluster_id), " clusters)")

# Detect available end types from column names
all_types <- c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M", "other")
avail_n_cols <- paste0("n_", all_types)
avail_types <- all_types[avail_n_cols %in% colnames(tss)]
message("  available end types: ", paste(avail_types, collapse = ", "))

# ── aggregate across samples ──────────────────────────────────────────────────
message("Aggregating counts across samples...")

agg_cols <- c("chr", "tss_pos", "strand", "initiator")
sum_cols <- intersect(paste0("n_", avail_types), colnames(tss))

tss_agg <- tss[, c(
  lapply(.SD[, ..agg_cols], first),
  lapply(.SD[, ..sum_cols, with = FALSE], sum),
  list(n_total = sum(n_total))
), by = cluster_id]

# Compute percentages
for (typ in avail_types) {
  n_col <- paste0("n_", typ)
  p_col <- paste0("pct_", typ)
  if (n_col %in% colnames(tss_agg)) {
    tss_agg[, (p_col) := 100 * get(n_col) / n_total]
  }
}

message("  aggregated clusters: ", nrow(tss_agg))

# ── filter: YR initiators ─────────────────────────────────────────────────────
yr_dinucs <- c("CA", "CG", "TA", "TG")
tss_yr <- tss_agg[initiator %in% yr_dinucs]
message("  after YR filter: ", nrow(tss_yr), " (", round(100 * nrow(tss_yr) / nrow(tss_agg), 1), "%)")

# ── filter: minimum reads ─────────────────────────────────────────────────────
tss_yr <- tss_yr[n_total >= min_reads]
message("  after min_reads filter (>=", min_reads, "): ", nrow(tss_yr))

# ── load TxDb for promoter annotation ─────────────────────────────────────────
message("Loading gene annotation: ", gtf_f)
txdb <- makeTxDbFromGFF(gtf_f, format = "gtf")

PROMOTER_UP <- 1000L
PROMOTER_DOWN <- 100L
tx <- transcripts(txdb)
promoter_gr <- reduce(promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN))
message("  promoter regions: ", length(promoter_gr))

# ── filter: promoter-localized ────────────────────────────────────────────────
tss_gr <- GRanges(
  seqnames = tss_yr$chr,
  ranges   = IRanges(start = tss_yr$tss_pos, width = 1),
  strand   = tss_yr$strand
)
in_promoter <- overlapsAny(tss_gr, promoter_gr, ignore.strand = TRUE)
tss_prom <- tss_yr[in_promoter]
message("  promoter-localized: ", nrow(tss_prom), " (", round(100 * nrow(tss_prom) / nrow(tss_yr), 1), "%)")

# ── load genome ───────────────────────────────────────────────────────────────
genome_spec <- params$genome_fasta_yeast
if (is.null(genome_spec) || genome_spec == "") {
  if (file.exists("data/sacCer3.fa")) {
    genome_spec <- "data/sacCer3.fa"
  } else {
    genome_spec <- "sacCer3"
  }
}
message("\nLoading genome: ", genome_spec)
genome <- load_genome(genome_spec)

# ── select high-enrichment TSS for each type ──────────────────────────────────
message("\nSelecting top ", top_pct, "% TSS for each end type...")

type_seqs <- list()
type_summary <- data.frame()

for (typ in avail_types) {
  pct_col <- paste0("pct_", typ)
  if (!pct_col %in% colnames(tss_prom)) next

  pct_vals <- tss_prom[[pct_col]]
  threshold <- quantile(pct_vals, 1 - top_pct/100, na.rm = TRUE)

  if (threshold < 5) {
    message("  ", typ, ": threshold too low (", round(threshold, 1), "%), skipping")
    next
  }

  high_tss <- tss_prom[pct_vals >= threshold]
  message("  ", typ, ": threshold=", round(threshold, 1), "%, n=", nrow(high_tss))

  seqs <- extract_sequences(high_tss, genome, upstream, downstream)

  if (length(seqs) >= 50) {
    type_seqs[[typ]] <- seqs
    type_summary <- rbind(type_summary, data.frame(
      type = typ,
      threshold_pct = round(threshold, 1),
      n_tss = length(seqs)
    ))
  } else {
    message("    too few sequences (", length(seqs), "), skipping")
  }
}

if (length(type_seqs) < 2) {
  stop("Not enough end types with sufficient sequences for comparison")
}

# ── write FASTA files ─────────────────────────────────────────────────────────
message("\nWriting FASTA files...")
for (typ in names(type_seqs)) {
  fasta_out <- file.path(outdir, paste0("tss_high_", typ, ".fa"))
  names(type_seqs[[typ]]) <- paste0("tss_", typ, "_", seq_along(type_seqs[[typ]]))
  writeXStringSet(type_seqs[[typ]], fasta_out)
  message("  ", fasta_out)
}

# ── compute PWMs ──────────────────────────────────────────────────────────────
pwms <- lapply(type_seqs, compute_pwm)

# ── generate sequence logos ───────────────────────────────────────────────────
message("\nGenerating sequence logos...")

logo_data <- lapply(type_seqs, function(s) as.character(s))
pos_labels <- seq(-upstream, downstream)

p_logo <- ggseqlogo(logo_data, method = "bits", ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  scale_x_continuous(
    breaks = seq(1, upstream + 1 + downstream, by = 10),
    labels = pos_labels[seq(1, length(pos_labels), by = 10)]
  ) +
  labs(x = "Position relative to TSS", y = "Bits")

ggsave(out_logo, p_logo, width = 12, height = 2.5 * length(type_seqs))
message("  saved: ", out_logo)

# ── difference heatmaps vs 1Sg ────────────────────────────────────────────────
if ("1Sg" %in% names(pwms)) {
  ref_pwm <- pwms[["1Sg"]]

  for (typ in setdiff(names(pwms), "1Sg")) {
    nuc_diff <- pwms[[typ]] - ref_pwm

    diff_df <- as.data.frame(t(nuc_diff)) %>%
      mutate(position = pos_labels) %>%
      pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "diff")

    p_diff <- ggplot(diff_df, aes(x = position, y = nucleotide, fill = diff)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                           midpoint = 0, limits = c(-0.3, 0.3),
                           name = "Freq diff") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      theme_minimal(base_size = 10) +
      labs(x = "Position relative to TSS", y = "Nucleotide",
           title = paste0("Nucleotide frequency difference: ", typ, " vs 1Sg")) +
      theme(panel.grid = element_blank())

    diff_out <- file.path(figdir, paste0("11_motif_diff_", typ, "_vs_1Sg.pdf"))
    ggsave(diff_out, p_diff, width = 12, height = 3)
    message("  saved: ", diff_out)
  }
}

# ── positional frequency line plots ──────────────────────────────────────────
freq_list <- list()
for (typ in names(pwms)) {
  freq_list[[typ]] <- as.data.frame(t(pwms[[typ]])) %>%
    mutate(position = pos_labels, type = typ) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "freq")
}
freq_all <- bind_rows(freq_list)
freq_all$type <- factor(freq_all$type, levels = names(pwms))

type_colors <- c(
  "1Sg" = "#1f77b4", "2Sgg" = "#ff7f0e", "3Sggg" = "#2ca02c",
  "1Sa" = "#d62728", "1St" = "#9467bd", "1Sc" = "#8c564b", "M" = "#7f7f7f"
)

p_freq <- ggplot(freq_all, aes(x = position, y = freq, color = type)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~nucleotide, ncol = 2, scales = "free_y") +
  scale_color_manual(values = type_colors, name = "End type") +
  theme_minimal(base_size = 10) +
  labs(x = "Position relative to TSS", y = "Frequency",
       title = "Positional nucleotide frequency by end type") +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

freq_out <- file.path(figdir, "11_motif_freq_comparison.pdf")
ggsave(freq_out, p_freq, width = 10, height = 7)
message("  saved: ", freq_out)

# ── summary table ─────────────────────────────────────────────────────────────
type_summary$window <- paste0("-", upstream, " to +", downstream)
summary_out <- file.path(outdir, "motif_summary.tsv")
write_tsv(type_summary, summary_out)
message("\n  saved summary: ", summary_out)

# ── key region analysis ───────────────────────────────────────────────────────
key_range <- (upstream - 4):(upstream + 6)  # positions -5 to +5
cat("\n=== Key region analysis (positions -5 to +5) ===\n")

for (typ in names(pwms)) {
  cat("\n", typ, " nucleotide frequencies:\n", sep = "")
  print(round(pwms[[typ]][, key_range], 3))
}

if ("1Sg" %in% names(pwms)) {
  ref_pwm <- pwms[["1Sg"]]
  for (typ in setdiff(names(pwms), "1Sg")) {
    cat("\nDifference (", typ, " minus 1Sg):\n", sep = "")
    print(round(pwms[[typ]][, key_range] - ref_pwm[, key_range], 3))
  }
}

message("\nDone. Outputs in ", outdir, " and ", figdir)
sessionInfo()
