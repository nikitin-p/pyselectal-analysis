#!/usr/bin/env Rscript
# Stage 08 step 1: extract fixed-width, strand-oriented promoter windows around
# dominant TSSs for seqArchR.
#
# Input : results/tss/tss_clustered.tsv (Stage 09 output)
# Output: results/seqarchr/promoter_windows.fa  (5'->3', fixed width, N-free)
#         results/seqarchr/window_metadata.tsv  (per-sequence end-type labels)

suppressMessages({
    library(optparse)
    library(Biostrings)
    library(GenomicRanges)
    library(Rsamtools)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--tss",        type = "character"),
    make_option("--genome",     type = "character"),
    make_option("--outdir",     type = "character"),
    make_option("--upstream",   type = "integer", default = 150L),
    make_option("--downstream", type = "integer", default = 50L),
    make_option("--min_total",  type = "double",  default = 2),
    make_option("--force",      action = "store_true", default = FALSE)
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
fa_out  <- file.path(opt$outdir, "promoter_windows.fa")
meta_out <- file.path(opt$outdir, "window_metadata.tsv")

if (file.exists(fa_out) && file.exists(meta_out) && !opt$force) {
    message("[skip] outputs exist; use --force to regenerate")
    quit(save = "no")
}

tss <- read.delim(opt$tss, stringsAsFactors = FALSE)
message(sprintf("Loaded %d TSS clusters", nrow(tss)))

# ---- filter on expression -------------------------------------------------
tss <- tss[tss$mean_total >= opt$min_total, ]
message(sprintf("After mean_total >= %s: %d", opt$min_total, nrow(tss)))

# ---- dominant end type per TSS --------------------------------------------
pct_cols <- end_type_cols("mean_pct_", END_TYPES_REAL)
pct_mat  <- as.matrix(tss[, pct_cols])
tss$end_type <- sub("^mean_pct_", "", pct_cols[max.col(pct_mat, ties.method = "first")])
tss$end_type_purity <- round(apply(pct_mat, 1, max), 2)

# ---- build strand-oriented windows ----------------------------------------
fa <- FaFile(opt$genome)
open(fa)
sl <- seqlengths(scanFaIndex(fa))

gr <- GRanges(
    seqnames = tss$chr,
    ranges   = IRanges(start = tss$tss_pos, width = 1),
    strand   = tss$strand
)
mcols(gr) <- tss[, c("cluster_id", "initiator", "hclust_group", "end_type",
                     "end_type_purity", "mean_total", pct_cols)]

# keep only chromosomes present in the genome, then set lengths so that
# out-of-bound windows can be detected rather than silently truncated
keep <- as.character(seqnames(gr)) %in% names(sl)
if (any(!keep)) message(sprintf("Dropping %d TSSs on chromosomes absent from genome", sum(!keep)))
gr <- gr[keep]
seqlevels(gr) <- names(sl)
seqlengths(gr) <- sl

win <- promoters(gr, upstream = opt$upstream, downstream = opt$downstream)
width_expected <- opt$upstream + opt$downstream

# drop windows running off a chromosome end — seqArchR requires equal widths
inb <- start(win) >= 1 & end(win) <= seqlengths(win)[as.character(seqnames(win))]
message(sprintf("Dropping %d out-of-bound windows", sum(!inb)))
win <- win[inb]

seqs <- getSeq(fa, win)          # strand-aware: minus strand returned 5'->3'
close(fa)

# ---- drop windows containing N --------------------------------------------
has_n <- vcountPattern("N", seqs) > 0
if (any(has_n)) message(sprintf("Dropping %d windows containing N", sum(has_n)))
seqs <- seqs[!has_n]
win  <- win[!has_n]

stopifnot(all(width(seqs) == width_expected))

names(seqs) <- mcols(win)$cluster_id
writeXStringSet(seqs, fa_out)

meta <- as.data.frame(mcols(win), stringsAsFactors = FALSE)
meta$window_chr    <- as.character(seqnames(win))
meta$window_start  <- start(win)
meta$window_end    <- end(win)
meta$window_strand <- as.character(strand(win))
write.table(meta, meta_out, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Wrote %d sequences of width %d", length(seqs), width_expected))
message("End-type composition:")
print(table(meta$end_type))

sessionInfo()
