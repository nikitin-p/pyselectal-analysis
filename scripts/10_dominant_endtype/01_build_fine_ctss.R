#!/usr/bin/env Rscript
# Stage 17 step 1: re-classify every read in each subsampled BAM by the literal
# sequence of its 5' soft clip, and aggregate to a long per-position per-sample
# typed CTSS table.
#
# This exists because scripts/05_tss_end_types/01_build_ctss.R collapses every
# clip that is not exactly g/gg/ggg into a single "other" category, which is
# 6.6 % of reads and contains the only clusters dominated by an unusual end.
#
# Usage (from project root):
#   Rscript scripts/10_dominant_endtype/01_build_fine_ctss.R \
#       [--samples config/samples.tsv] \
#       [--outdir  results/dominance]  \
#       [--params  config/params.yaml] \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(Rsamtools)
  library(data.table)
  library(readr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

samples_f <- get_arg("--samples", "config/samples.tsv")
outdir    <- get_arg("--outdir",  "results/dominance")
params_f  <- get_arg("--params",  "config/params.yaml")
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "dom_helpers.R"))

params   <- yaml.load_file(params_f)
max_clip <- as.integer(params$dom_max_clip)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_gz <- file.path(outdir, "fine_typed_ctss.tsv.gz")

if (!force && file.exists(out_gz)) {
  message("Output exists, skipping (use --force to rerun): ", out_gz)
  quit(save = "no", status = 0)
}

bam_to_fine_ctss <- function(bam_path, sample_id, max_clip) {
  message("  reading: ", sample_id)
  param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = FALSE,
                       isSecondaryAlignment = FALSE,
                       isSupplementaryAlignment = FALSE),
    what = "seq")
  ga <- readGAlignments(bam_path, param = param)
  if (length(ga) == 0) {
    warning("No alignments in: ", bam_path)
    return(NULL)
  }
  strand_chr <- as.character(strand(ga))
  is_plus    <- strand_chr %in% c("+", "*")   # treat unstranded as +

  pos <- integer(length(ga))
  pos[ is_plus] <- start(ga)[ is_plus]
  pos[!is_plus] <- end(ga)[!is_plus]

  seq_str <- as.character(mcols(ga)$seq)
  clen    <- clip5_len(cigar(ga), is_plus)
  type    <- fine_type(clen, clip5_seq(seq_str, clen, is_plus), max_clip)

  dt <- data.table(chr = as.character(seqnames(ga)), pos = pos,
                   strand = strand_chr, type = type)
  dt[strand == "*", strand := "+"]
  agg <- dt[, .(n = .N), by = .(chr, pos, strand, type)]
  agg[, sample_id := sample_id]
  agg[]
}

meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)
stopifnot("sample_id" %in% colnames(meta), nrow(meta) > 0)
has_bam_col <- "bam" %in% colnames(meta)

message("Stage 17 step 1: fine-typed CTSS")
message("  samples: ", nrow(meta), "   max_clip: ", max_clip)

results <- lapply(seq_len(nrow(meta)), function(i) {
  sid <- meta$sample_id[i]
  bam_path <- if (has_bam_col && !is.na(meta$bam[i])) meta$bam[i] else
    file.path(params$results_bam_sub, paste0(sid, ".bam"))
  if (!file.exists(bam_path)) {
    warning("BAM not found, skipping: ", bam_path)
    return(NULL)
  }
  bam_to_fine_ctss(bam_path, sid, max_clip)
})

results <- Filter(Negate(is.null), results)
if (length(results) == 0) stop("No samples processed.")

combined <- rbindlist(results)
setcolorder(combined, c("chr", "pos", "strand", "sample_id", "type", "n"))
setorder(combined, chr, strand, pos, sample_id, -n)

fwrite(combined, out_gz, sep = "\t", compress = "gzip")
message("Written: ", out_gz, "  (", nrow(combined), " rows, ",
        uniqueN(combined$type), " distinct types)")

sessionInfo()
