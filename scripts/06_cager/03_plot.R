#!/usr/bin/env Rscript
# Stage 10 step 3: diagnostic plots from the CAGEexp object.
#
# Outputs:
#   results/figures/10_correlation.pdf          — sample-sample Pearson correlation
#   results/figures/10_reverse_cumulatives.pdf  — reverse cumulative CTSS distributions
#   results/figures/10_iqw.pdf                  — interquantile width distribution
#   results/figures/10_annotation_enrichment.pdf — sense/antisense/intergenic breakdown
#                                                   (only if annotation_enrichment.tsv exists)
#
# Usage (from project root):
#   Rscript scripts/06_cager/03_plot.R \
#       [--rds      results/cager/cageexp.rds]              \
#       [--annot    results/cager/annotation_enrichment.tsv] \
#       [--tagclusters results/cager/tagclusters_annotated.tsv] \
#       [--outdir   results/figures]                         \
#       [--params   config/params.yaml]                      \
#       [--force]

suppressPackageStartupMessages({
  library(CAGEr)
  library(data.table)
  library(ggplot2)
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

rds_f    <- get_arg("--rds",    "results/cager/cageexp.rds")
annot_f  <- get_arg("--annot",  "results/cager/annotation_enrichment.tsv")
tc_f     <- get_arg("--tagclusters", "results/cager/tagclusters_annotated.tsv")
outdir   <- get_arg("--outdir", "results/figures")
params_f <- get_arg("--params", "config/params.yaml")
force    <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_corr   <- file.path(outdir, "10_correlation.pdf")
out_rcum   <- file.path(outdir, "10_reverse_cumulatives.pdf")
out_iqw    <- file.path(outdir, "10_iqw.pdf")
out_annot  <- file.path(outdir, "10_annotation_enrichment.pdf")

if (!force && all(file.exists(c(out_corr, out_rcum, out_iqw)))) {
  message("Outputs exist, skipping (use --force to rerun).")
  quit(save = "no", status = 0)
}

if (!file.exists(rds_f)) {
  stop("Input not found: ", rds_f,
       "\nRun step 1 first (scripts/06_cager/01_cageexp.R).")
}

params <- yaml.load_file(params_f)
iq_tpm_threshold <- as.numeric(params$cager_iqw_tpm_threshold %||% 3)

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x[[1]]) && x[[1]] != "") x else y

message("Stage 10 step 3: plotting")
message("  loading RDS: ", rds_f)
ce <- readRDS(rds_f)

# ── 1. Correlation heatmap ────────────────────────────────────────────────────
message("  correlation: ", out_corr)
pdf(out_corr, width = 6, height = 5)
tryCatch(
  {
    p <- plotCorrelation2(
      ce,
      samples            = "all",
      tagCountThreshold  = 1,
      applyThresholdBoth = FALSE,
      method             = "pearson"
    )
    if (!is.null(p)) print(p)
  },
  error = function(e) message("  plotCorrelation2 failed: ", e$message)
)
dev.off()

# ── 2. Reverse cumulatives ────────────────────────────────────────────────────
message("  reverse cumulatives: ", out_rcum)
pdf(out_rcum, width = 7, height = 5)
tryCatch(
  plotReverseCumulatives(ce, fitInRange = c(5, 1000)),
  error = function(e) message("  plotReverseCumulatives failed: ", e$message)
)
dev.off()

# ── 3. Interquantile width ────────────────────────────────────────────────────
message("  IQW: ", out_iqw)
pdf(out_iqw, width = 7, height = 5)
tryCatch(
  {
    p <- plotInterquantileWidth(
      ce,
      clusters     = "tagClusters",
      tpmThreshold = iq_tpm_threshold,
      qLow         = 0.1,
      qUp          = 0.9
    )
    if (!is.null(p)) print(p)
  },
  error = function(e) message("  plotInterquantileWidth failed: ", e$message)
)
dev.off()

# ── 4. IQW distribution from tagclusters.tsv.gz (supplementary) ──────────────
tc_gz <- sub("_annotated\\.tsv$", ".tsv.gz",
             get_arg("--tagclusters", "results/cager/tagclusters_annotated.tsv"))
if (file.exists(tc_gz)) {
  tc_raw <- fread(cmd = paste("zcat", shQuote(tc_gz)))
  if ("iq_width" %in% names(tc_raw) && any(!is.na(tc_raw$iq_width))) {
    out_iqw2 <- file.path(outdir, "10_iqw_distribution.pdf")
    message("  IQW distribution (all samples): ", out_iqw2)

    tc_valid <- tc_raw[!is.na(iq_width) & iq_width > 0]

    p_iq <- ggplot(tc_valid, aes(x = iq_width, colour = sample_id, fill = sample_id)) +
      geom_density(alpha = 0.15, linewidth = 0.6) +
      scale_x_log10() +
      labs(
        title = "Tag cluster interquantile width (q0.1–q0.9)",
        x     = "IQW (bp, log10 scale)",
        y     = "Density",
        colour = "Sample", fill = "Sample"
      ) +
      theme_bw(base_size = 11)

    pdf(out_iqw2, width = 7, height = 4)
    print(p_iq)
    dev.off()
  }
}

# ── 5. Annotation enrichment bar ──────────────────────────────────────────────
if (file.exists(tc_f)) {
  cons_dt <- fread(tc_f)
  if ("annotation" %in% names(cons_dt) && any(!is.na(cons_dt$annotation))) {
    message("  annotation enrichment: ", out_annot)

    enr <- cons_dt[!is.na(annotation), .N, by = annotation]
    enr[, pct := 100 * N / sum(N)]
    enr$annotation <- factor(enr$annotation,
                              levels = c("sense_genic", "antisense_genic", "intergenic"))

    p_enr <- ggplot(enr, aes(x = annotation, y = pct, fill = annotation)) +
      geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("n=%d\n%.1f%%", N, pct)),
                vjust = -0.3, size = 3) +
      scale_fill_manual(
        values = c(sense_genic     = "#4361ee",
                   antisense_genic = "#f72585",
                   intergenic      = "#adb5bd"),
        guide = "none"
      ) +
      labs(
        title = "CAGEr tag cluster genomic distribution",
        x     = NULL, y     = "% clusters"
      ) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    pdf(out_annot, width = 4, height = 4)
    print(p_enr)
    dev.off()
  }
} else if (file.exists(annot_f)) {
  enr <- fread(annot_f)
  if (nrow(enr) > 0) {
    message("  annotation enrichment (from TSV): ", out_annot)
    enr$annotation <- factor(enr$annotation,
                              levels = c("sense_genic", "antisense_genic", "intergenic"))

    p_enr <- ggplot(enr, aes(x = annotation, y = pct, fill = annotation)) +
      geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("n=%d\n%.1f%%", N, pct)),
                vjust = -0.3, size = 3) +
      scale_fill_manual(
        values = c(sense_genic     = "#4361ee",
                   antisense_genic = "#f72585",
                   intergenic      = "#adb5bd"),
        guide = "none"
      ) +
      labs(title = "CAGEr tag cluster genomic distribution",
           x = NULL, y = "% clusters") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    pdf(out_annot, width = 4, height = 4)
    print(p_enr)
    dev.off()
  }
}

message("Done. Outputs in: ", outdir)
sessionInfo()
