#!/usr/bin/env Rscript
# Stage 8 step 1: classify every read in each subsampled BAM by 5'-end type
# (1Sg / 2Sgg / 3Sggg / M / other) and aggregate to a per-position per-sample typed CTSS table.
#
# Classification is strand-aware using CIGAR + sequence:
#   + strand: 5' soft clip = start of CIGAR; first stored base checked for G
#   - strand: 5' soft clip = end of CIGAR; last stored base checked for C (RC of G)
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

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

samples_f <- get_arg("--samples", "config/samples.tsv")
outdir    <- get_arg("--outdir",  "results/tss")
params_f  <- get_arg("--params",  "config/params.yaml")
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

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

  # First original 5' base (strand-aware):
  #   + strand: first char of stored seq
  #   - strand: complement of last char of stored seq (BAM stores RC)
  seq_str <- as.character(mcols(ga)$seq)
  comp    <- c(A = "T", T = "A", G = "C", C = "G", N = "N")

  has_clip <- clen > 0
  plus_clip  <- is_plus & has_clip
  minus_clip <- !is_plus & has_clip

  # Get first 5' base for single soft-clips
  first_base <- character(length(ga))
  if (any(plus_clip)) {
    first_base[plus_clip] <- substr(seq_str[plus_clip], 1, 1)
  }
  if (any(minus_clip)) {
    last_chars <- substr(seq_str[minus_clip],
                         nchar(seq_str[minus_clip]),
                         nchar(seq_str[minus_clip]))
    first_base[minus_clip] <- comp[last_chars]  # reverse complement
  }

  # Check if first 2 bases are GG for 2Sgg
  first_two_gg <- logical(length(ga))
  len2_clip <- clen == 2 & has_clip
  if (any(len2_clip & is_plus)) {
    idx <- which(len2_clip & is_plus)
    first_two_gg[idx] <- substr(seq_str[idx], 1, 2) == "GG"
  }
  if (any(len2_clip & !is_plus)) {
    idx <- which(len2_clip & !is_plus)
    last2 <- substr(seq_str[idx], nchar(seq_str[idx]) - 1, nchar(seq_str[idx]))
    first_two_gg[idx] <- last2 == "CC"  # original 5' GG stored as CC after RC
  }

  # Check if first 3 bases are GGG for 3Sggg
  first_three_ggg <- logical(length(ga))
  len3_clip <- clen == 3 & has_clip
  if (any(len3_clip & is_plus)) {
    idx <- which(len3_clip & is_plus)
    first_three_ggg[idx] <- substr(seq_str[idx], 1, 3) == "GGG"
  }
  if (any(len3_clip & !is_plus)) {
    idx <- which(len3_clip & !is_plus)
    last3 <- substr(seq_str[idx], nchar(seq_str[idx]) - 2, nchar(seq_str[idx]))
    first_three_ggg[idx] <- last3 == "CCC"  # original 5' GGG stored as CCC after RC
  }

  # Assign type with granular categories
  type <- character(length(ga))
  type[clen == 0]                        <- "M"
  type[clen == 1 & first_base == "G"]    <- "1Sg"
  type[clen == 1 & first_base == "A"]    <- "1Sa"
  type[clen == 1 & first_base == "T"]    <- "1St"
  type[clen == 1 & first_base == "C"]    <- "1Sc"
  type[clen == 2 & first_two_gg]         <- "2Sgg"
  type[clen == 3 & first_three_ggg]      <- "3Sggg"
  type[type == ""]                       <- "other"

  # Aggregate: (chr, pos, strand, type) → count
  dt <- data.table(
    chr    = as.character(seqnames(ga)),
    pos    = tss_pos,
    strand = strand_chr,
    type   = type
  )
  dt[strand == "*", strand := "+"]

  agg <- dt[, .(n = .N), by = .(chr, pos, strand, type)]

  # Pivot to wide
  wide <- dcast(agg, chr + pos + strand ~ type, value.var = "n", fill = 0L)

  # Ensure all type columns present
  all_types <- c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M", "other")
  for (col in all_types) {
    if (!col %in% colnames(wide)) wide[[col]] <- 0L
  }

  wide[, n_total := `1Sg` + `2Sgg` + `3Sggg` + `1Sa` + `1St` + `1Sc` + M + other]
  setnames(wide, all_types,
           paste0("n_", all_types))
  wide[, sample_id := sample_id]

  wide
}

# ── load sample manifest ──────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)
stopifnot("sample_id" %in% colnames(meta))
stopifnot(nrow(meta) > 0)

has_bam_col <- "bam" %in% colnames(meta)

# ── output path ───────────────────────────────────────────────────────────────
out_gz <- file.path(outdir, "typed_ctss.tsv.gz")

if (!force && file.exists(out_gz)) {
  message("Output exists, skipping (use --force to rerun): ", out_gz)
  quit(save = "no", status = 0)
}

message("Stage 8 step 1: typed CTSS table")
message("  samples: ", nrow(meta))

# ── process all samples ───────────────────────────────────────────────────────
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

# ── write output ──────────────────────────────────────────────────────────────
fwrite(combined, out_gz, sep = "\t", compress = "gzip")
message("Written: ", out_gz, "  (", nrow(combined), " rows)")

sessionInfo()
