#!/usr/bin/env Rscript
# Stage 8 step 7: Sequence motifs at positions where each end type OCCURS
#
# Unlike 06_motif_analysis.R (which requires TSS dominance), this script extracts
# sequences at every position where a given type has reads, weighted by read count.
# This allows comparison of all end types including rare ones (2Sgg, 3Sggg, 1Sa, etc.)
#
# Usage:
#   Rscript scripts/05_tss_end_types/07_motif_by_occurrence.R \
#       [--ctss     results/tss/typed_ctss.tsv.gz] \
#       [--params   config/params.yaml]            \
#       [--gtf      data/sacCer3.ensGene.gtf]      \
#       [--outdir   results/motifs]                \
#       [--figdir   results/figures]               \
#       [--max_pos  5000]                          \
#       [--upstream 20]                            \
#       [--downstream 20]                          \
#       [--promoter_only]                          \
#       [--yr_only]                                \
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

ctss_f       <- get_arg("--ctss",       "results/tss/typed_ctss.tsv.gz")
params_f     <- get_arg("--params",     "config/params.yaml")
gtf_f        <- get_arg("--gtf",        "data/sacCer3.ensGene.gtf")
outdir       <- get_arg("--outdir",     "results/motifs")
figdir       <- get_arg("--figdir",     "results/figures")
max_pos      <- as.integer(get_arg("--max_pos", "5000"))
upstream     <- as.integer(get_arg("--upstream", "20"))
downstream   <- as.integer(get_arg("--downstream", "20"))
promoter_only <- has_flag("--promoter_only")
yr_only       <- has_flag("--yr_only")
force         <- has_flag("--force")

validate_args()

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

out_logo <- file.path(figdir, "12_motif_by_occurrence.pdf")
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

# ── load typed CTSS ───────────────────────────────────────────────────────────
message("Loading typed CTSS: ", ctss_f)
ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
message("  rows: ", nrow(ctss))

# Aggregate across samples
message("Aggregating across samples...")
type_cols <- grep("^n_", colnames(ctss), value = TRUE)
type_cols <- setdiff(type_cols, "n_total")

ctss_agg <- ctss[, c(
  lapply(.SD[, ..type_cols], sum)
), by = .(chr, pos, strand)]

message("  unique positions: ", nrow(ctss_agg))

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

# ── optionally filter to promoters ────────────────────────────────────────────
if (promoter_only) {
  message("\nFiltering to promoter regions...")
  txdb <- makeTxDbFromGFF(gtf_f, format = "gtf")
  tx <- transcripts(txdb)
  promoter_gr <- reduce(promoters(tx, upstream = 1000, downstream = 100))

  pos_gr <- GRanges(
    seqnames = ctss_agg$chr,
    ranges = IRanges(start = ctss_agg$pos, width = 1),
    strand = ctss_agg$strand
  )
  in_promoter <- overlapsAny(pos_gr, promoter_gr, ignore.strand = TRUE)
  ctss_agg <- ctss_agg[in_promoter]
  message("  positions in promoters: ", nrow(ctss_agg))
}

# ── optionally filter to YR initiators ────────────────────────────────────────
if (yr_only) {
  message("\nFiltering to YR initiator positions...")
  pos_gr <- GRanges(
    seqnames = ctss_agg$chr,
    ranges = IRanges(start = ctss_agg$pos, width = 1),
    strand = ctss_agg$strand
  )
  pos_gr <- harmonize_seqlevels(pos_gr, genome)

  win <- promoters(pos_gr, upstream = 1, downstream = 1)
  win <- trim(win)
  keep_width <- width(win) == 2

  initiator <- rep(NA_character_, length(pos_gr))
  initiator[keep_width] <- as.character(getSeq(genome, win[keep_width]))

  yr_dinucs <- c("CA", "CG", "TA", "TG")
  is_yr <- initiator %in% yr_dinucs
  ctss_agg <- ctss_agg[is_yr]
  message("  positions with YR initiator: ", nrow(ctss_agg))
}

# ── extract sequences for each type ───────────────────────────────────────────
message("\nExtracting sequences for each end type...")

end_types <- gsub("^n_", "", type_cols)
end_types <- setdiff(end_types, "other")  # exclude "other" for cleaner comparison

extract_weighted_sequences <- function(dt, type_col, genome, up, down, max_n) {
  n_col <- paste0("n_", type_col)

  # Get positions where this type has reads
  has_type <- dt[[n_col]] > 0
  dt_type <- dt[has_type]

  if (nrow(dt_type) == 0) return(NULL)

  # Weight by read count, sample up to max_n positions
  weights <- dt_type[[n_col]]

  if (nrow(dt_type) > max_n) {
    # Weighted sampling
    set.seed(42)
    idx <- sample(nrow(dt_type), max_n, prob = weights, replace = FALSE)
    dt_type <- dt_type[idx]
  }

  gr <- GRanges(
    seqnames = dt_type$chr,
    ranges = IRanges(start = dt_type$pos, width = 1),
    strand = dt_type$strand
  )
  gr <- harmonize_seqlevels(gr, genome)

  # Extend around TSS
  gr_ext <- gr
  start(gr_ext) <- ifelse(strand(gr) == "+", start(gr) - up, start(gr) - down)
  end(gr_ext) <- ifelse(strand(gr) == "+", end(gr) + down, end(gr) + up)
  gr_ext <- trim(gr_ext)

  expected_width <- up + 1 + down
  keep <- width(gr_ext) == expected_width
  gr_ext <- gr_ext[keep]

  if (length(gr_ext) == 0) return(NULL)

  # getSeq() already returns the reverse complement for minus-strand ranges,
  # so the sequences come back 5'->3' along the transcript. Do NOT revcomp again.
  seqs <- getSeq(genome, gr_ext)

  seqs
}

type_seqs <- list()
type_counts <- data.frame()

for (typ in end_types) {
  seqs <- extract_weighted_sequences(ctss_agg, typ, genome, upstream, downstream, max_pos)

  if (!is.null(seqs) && length(seqs) >= 100) {
    type_seqs[[typ]] <- seqs
    n_col <- paste0("n_", typ)
    total_reads <- sum(ctss_agg[[n_col]])
    type_counts <- rbind(type_counts, data.frame(
      type = typ,
      n_positions = length(seqs),
      total_reads = total_reads
    ))
    message("  ", typ, ": ", length(seqs), " positions (", format(total_reads, big.mark=","), " total reads)")
  } else {
    message("  ", typ, ": insufficient data, skipping")
  }
}

if (length(type_seqs) < 2) {
  stop("Not enough end types with sufficient data")
}

# ── compute PWMs ──────────────────────────────────────────────────────────────
compute_pwm <- function(seqs) {
  mat <- consensusMatrix(seqs, as.prob = TRUE)[1:4, ]
  rownames(mat) <- c("A", "C", "G", "T")
  mat
}

pwms <- lapply(type_seqs, compute_pwm)

# ── generate sequence logos ───────────────────────────────────────────────────
message("\nGenerating sequence logos...")

logo_data <- lapply(type_seqs, function(s) as.character(s))
pos_labels <- seq(-upstream, downstream)

# Order types logically
type_order <- c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M")
type_order <- intersect(type_order, names(logo_data))
logo_data <- logo_data[type_order]
pwms <- pwms[type_order]

p_logo <- ggseqlogo(logo_data, method = "bits", ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 11)) +
  scale_x_continuous(
    breaks = seq(1, upstream + 1 + downstream, by = 5),
    labels = pos_labels[seq(1, length(pos_labels), by = 5)]
  ) +
  labs(x = "Position relative to TSS", y = "Bits")

ggsave(out_logo, p_logo, width = 12, height = 2.2 * length(logo_data))
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
           title = paste0(typ, " vs 1Sg")) +
      theme(panel.grid = element_blank())

    diff_out <- file.path(figdir, paste0("12_motif_diff_", typ, "_vs_1Sg.pdf"))
    ggsave(diff_out, p_diff, width = 12, height = 3)
    message("  saved: ", diff_out)
  }
}

# ── combined difference heatmap ───────────────────────────────────────────────
if ("1Sg" %in% names(pwms) && length(pwms) > 1) {
  ref_pwm <- pwms[["1Sg"]]

  diff_list <- list()
  for (typ in setdiff(names(pwms), "1Sg")) {
    nuc_diff <- pwms[[typ]] - ref_pwm
    diff_df <- as.data.frame(t(nuc_diff)) %>%
      mutate(position = pos_labels, type = typ) %>%
      pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "diff")
    diff_list[[typ]] <- diff_df
  }

  diff_all <- bind_rows(diff_list)
  diff_all$type <- factor(diff_all$type, levels = setdiff(type_order, "1Sg"))

  p_diff_all <- ggplot(diff_all, aes(x = position, y = nucleotide, fill = diff)) +
    geom_tile() +
    facet_wrap(~type, ncol = 1) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-0.3, 0.3),
                         name = "Freq diff\nvs 1Sg") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    theme_minimal(base_size = 10) +
    labs(x = "Position relative to TSS", y = "") +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold"))

  diff_all_out <- file.path(figdir, "12_motif_diff_all_vs_1Sg.pdf")
  ggsave(diff_all_out, p_diff_all, width = 14, height = 2 * (length(pwms) - 1))
  message("  saved: ", diff_all_out)
}

# ── frequency line plots ──────────────────────────────────────────────────────
freq_list <- list()
for (typ in names(pwms)) {
  freq_list[[typ]] <- as.data.frame(t(pwms[[typ]])) %>%
    mutate(position = pos_labels, type = typ) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "freq")
}
freq_all <- bind_rows(freq_list)
freq_all$type <- factor(freq_all$type, levels = type_order)

type_colors <- c(
  "1Sg" = "#1f77b4", "2Sgg" = "#ff7f0e", "3Sggg" = "#2ca02c",
  "1Sa" = "#d62728", "1St" = "#9467bd", "1Sc" = "#8c564b", "M" = "#7f7f7f"
)

p_freq <- ggplot(freq_all, aes(x = position, y = freq, color = type)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~nucleotide, ncol = 2, scales = "free_y") +
  scale_color_manual(values = type_colors, name = "End type") +
  theme_minimal(base_size = 10) +
  labs(x = "Position relative to TSS", y = "Frequency",
       title = "Nucleotide frequency by end type (occurrence-based)") +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

freq_out <- file.path(figdir, "12_motif_freq_by_occurrence.pdf")
ggsave(freq_out, p_freq, width = 10, height = 7)
message("  saved: ", freq_out)

# ── summary ───────────────────────────────────────────────────────────────────
type_counts$window <- paste0("-", upstream, " to +", downstream)
summary_out <- file.path(outdir, "motif_occurrence_summary.tsv")
write_tsv(type_counts, summary_out)
message("\n  saved summary: ", summary_out)

# ── key region analysis ───────────────────────────────────────────────────────
key_range <- (upstream - 4):(upstream + 6)
cat("\n=== Key region analysis (positions -5 to +5) ===\n")

for (typ in names(pwms)) {
  cat("\n", typ, ":\n", sep = "")
  print(round(pwms[[typ]][, key_range], 3))
}

if ("1Sg" %in% names(pwms)) {
  ref_pwm <- pwms[["1Sg"]]
  cat("\n=== Differences vs 1Sg ===\n")
  for (typ in setdiff(names(pwms), "1Sg")) {
    cat("\n", typ, " - 1Sg:\n", sep = "")
    print(round(pwms[[typ]][, key_range] - ref_pwm[, key_range], 3))
  }
}

message("\nDone. Outputs in ", outdir, " and ", figdir)
sessionInfo()
