#!/usr/bin/env Rscript
# Stage 8 step 1: classify every read in each subsampled BAM by 5'-end type
# (the full END_TYPES vocabulary from scripts/utils/end_types.R: 1Sg / 2Sgg /
# 3Sggg / 2Stg / 1Sa / 2Sga / 1St / 1Sc / M / other) and aggregate to a
# per-position per-sample typed CTSS table.
#
# Classification is strand-aware using CIGAR + sequence. The 5' soft clip is
# extracted in the read's original 5'->3' orientation (matching clip5_seq() in
# scripts/10_dominant_endtype/dom_helpers.R, so that Stage 8 and Stage 17
# label the same read identically) and then mapped to a first-class type by
# classify_clip(), or "other" if it isn't in END_TYPES_REAL:
#   + strand: 5' soft clip = start of CIGAR, stored sequence read as-is
#   - strand: 5' soft clip = end of CIGAR, stored sequence reverse-complemented
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/01_build_ctss.R \
#       [--samples  config/samples.tsv] \
#       [--outdir   results/tss]        \
#       [--params   config/params.yaml] \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(Rsamtools)
  library(data.table)
  library(readr)
  library(yaml)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
# When this script is source()d (e.g. by test/test_build_ctss_types.R) rather
# than run directly, commandArgs(FALSE) reflects the OUTER script's --file,
# not this file's own path, so the path above can resolve to the wrong
# directory depth. Fall back to the project-root-relative path in that case.
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

get_arg <- function(flag, default = NULL, args = commandArgs(trailingOnly = TRUE)) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

# ── CIGAR helpers ─────────────────────────────────────────────────────────────

# Length of 5'-end soft clip (strand-aware).
# + strand: leading soft clip (e.g. "2S48M" → 2)
# - strand: trailing soft clip (e.g. "48M2S" → 2)
clip5_len <- function(cigar, is_plus) {
  len <- integer(length(cigar))
  len[ is_plus] <- as.integer(ifelse(
    grepl("^\\d+S", cigar[is_plus]),
    sub("^(\\d+)S.*", "\\1", cigar[is_plus]),
    "0"
  ))
  len[!is_plus] <- as.integer(ifelse(
    grepl("\\d+S$", cigar[!is_plus]),
    sub(".*?(\\d+)S$", "\\1", cigar[!is_plus]),
    "0"
  ))
  len
}

# The 5' soft clip in the read's ORIGINAL 5'->3' orientation, matching the
# convention of clip5_seq() in scripts/10_dominant_endtype/dom_helpers.R, so
# that Stage 8 and Stage 17 label the same read identically.
.comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
.revcomp <- function(s) {
  vapply(strsplit(s, "", fixed = TRUE),
         function(ch) paste(rev(.comp[ch]), collapse = ""),
         character(1))
}

clip5_seq_plus_minus <- function(seq_str, clen, is_plus) {
  out <- rep("", length(seq_str))
  p <- is_plus & clen > 0
  m <- !is_plus & clen > 0
  if (any(p)) out[p] <- substr(seq_str[p], 1, clen[p])
  if (any(m)) {
    n <- nchar(seq_str[m])
    out[m] <- .revcomp(substr(seq_str[m], n - clen[m] + 1, n))
  }
  out
}

# Map a clip length and clip sequence to a first-class type, or "other".
# The label is <length>S<lower-cased clip>, so membership in END_TYPES_REAL is
# the whole test and adding a type needs no change here.
classify_clip <- function(clen, clipseq) {
  lbl <- ifelse(clen == 0L, "M", paste0(clen, "S", tolower(clipseq)))
  ifelse(lbl %in% END_TYPES_REAL, lbl, "other")
}

# ── per-BAM read classification → typed CTSS ─────────────────────────────────
bam_to_typed_ctss <- function(bam_path, sample_id) {
  message("  reading: ", sample_id)

  param <- ScanBamParam(
    flag = scanBamFlag(
      isUnmappedQuery          = FALSE,
      isSecondaryAlignment     = FALSE,
      isSupplementaryAlignment = FALSE
    ),
    what = "seq"
  )

  ga <- readGAlignments(bam_path, param = param)
  if (length(ga) == 0) {
    warning("No alignments in: ", bam_path)
    return(NULL)
  }

  strand_chr <- as.character(strand(ga))
  is_plus    <- strand_chr %in% c("+", "*")  # treat unstranded as +

  # TSS position: leftmost mapped base for +, rightmost for -
  tss_pos          <- integer(length(ga))
  tss_pos[ is_plus]  <- start(ga)[ is_plus]
  tss_pos[!is_plus]  <- end(ga)[!is_plus]

  # 5'-end soft clip length
  cig  <- cigar(ga)
  clen <- clip5_len(cig, is_plus)

  seq_str <- as.character(mcols(ga)$seq)

  clipseq <- clip5_seq_plus_minus(seq_str, clen, is_plus)
  type    <- classify_clip(clen, clipseq)

  # Aggregate: (chr, pos, strand, type) → count
  # type is factored on END_TYPES so dcast()'s implicit sort follows the
  # canonical mechanistic order rather than falling back to alphabetical
  # order for a plain character column.
  dt <- data.table(
    chr    = as.character(seqnames(ga)),
    pos    = tss_pos,
    strand = strand_chr,
    type   = factor(type, levels = END_TYPES)
  )
  dt[strand == "*", strand := "+"]

  agg <- dt[, .(n = .N), by = .(chr, pos, strand, type)]

  # Pivot to wide
  wide <- dcast(agg, chr + pos + strand ~ type, value.var = "n", fill = 0L)

  # Ensure all type columns present. dcast() only emits columns for factor
  # levels actually observed in this BAM, so a type with zero reads here is
  # added back at the end of colnames(wide) and must be reordered onto
  # END_TYPES's canonical position rather than left trailing.
  for (col in END_TYPES) {
    if (!col %in% colnames(wide)) wide[[col]] <- 0L
  }
  setcolorder(wide, c("chr", "pos", "strand", END_TYPES))

  wide[, n_total := rowSums(.SD), .SDcols = END_TYPES]
  setnames(wide, END_TYPES, end_type_cols("n_", END_TYPES))
  wide[, sample_id := sample_id]

  wide
}

main <- function() {
  # ── argument parsing ──────────────────────────────────────────────────────
  args <- commandArgs(trailingOnly = TRUE)

  samples_f <- get_arg("--samples", "config/samples.tsv", args = args)
  outdir    <- get_arg("--outdir",  "results/tss",        args = args)
  params_f  <- get_arg("--params",  "config/params.yaml", args = args)
  force     <- "--force" %in% args

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  params <- yaml.load_file(params_f)

  # ── load sample manifest ────────────────────────────────────────────────
  meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                   show_col_types = FALSE)
  stopifnot("sample_id" %in% colnames(meta))
  stopifnot(nrow(meta) > 0)

  has_bam_col <- "bam" %in% colnames(meta)

  # ── output path ──────────────────────────────────────────────────────────
  out_gz <- file.path(outdir, "typed_ctss.tsv.gz")

  if (!force && file.exists(out_gz)) {
    message("Output exists, skipping (use --force to rerun): ", out_gz)
    quit(save = "no", status = 0)
  }

  message("Stage 8 step 1: typed CTSS table")
  message("  samples: ", nrow(meta))

  # ── process all samples ─────────────────────────────────────────────────
  results <- lapply(seq_len(nrow(meta)), function(i) {
    sid <- meta$sample_id[i]
    bam_path <- if (has_bam_col && !is.na(meta$bam[i])) {
      meta$bam[i]
    } else {
      file.path(params$results_bam_sub, paste0(sid, ".bam"))
    }

    if (!file.exists(bam_path)) {
      warning("BAM not found, skipping: ", bam_path)
      return(NULL)
    }
    bam_to_typed_ctss(bam_path, sid)
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) stop("No samples processed.")

  combined <- rbindlist(results)

  # ── write output ────────────────────────────────────────────────────────
  fwrite(combined, out_gz, sep = "\t", compress = "gzip")
  message("Written: ", out_gz, "  (", nrow(combined), " rows)")

  sessionInfo()
}

if (!nzchar(Sys.getenv("BUILD_CTSS_TESTING"))) {
  main()
}
