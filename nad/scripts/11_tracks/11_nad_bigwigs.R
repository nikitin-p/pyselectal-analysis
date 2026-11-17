#!/usr/bin/env Rscript
# NAD branch step 11: stranded 5'-end bigWigs of NAD captureSeq, for browsing
# beside the CAGE end-type tracks of Stage 21.
#
# Per strain (WT, dxo1 rai1 npy1) and strand, three tracks are written:
#   S     mean over the three S (+ADPRC) replicates of reads / size factor;
#   N     the same for the three N (no ADPRC) replicates;
#   S-N   S minus N at every position, negative where N exceeds S.
# The size factors are those of step 5 (DESeq2 median of ratios over unmasked
# tested positions), recomputed here exactly as step 8 does, so that the S and
# N tracks equal the meanS and meanN columns of <strain>_S_vs_N.tsv.gz at every
# tested position and the difference track shows what DESeq2 compared. Library
# depth (CPM) is deliberately not used: N libraries hold mostly background, and
# depth scaling would make S and N incomparable at mRNA starts. A log2 ratio
# track is not written either, because single-nucleotide counts are mostly 0-3
# and a ratio of such counts is noise.
#
# A bigBed of the step-5 calls (padj < 0.05, log2FC > 1) is written per strain
# beside the bigWigs, coloured by the base at +1 (A = red, other = grey), so
# that called positions can be found without reading the tables.
#
# Usage (from project root):
#   Rscript nad/scripts/11_tracks/11_nad_bigwigs.R \
#       [--matrix   nad/results/04_ctss/matrix_reads.tsv.gz]        \
#       [--manifest nad/scripts/00_download/priority1_runs.tsv]      \
#       [--s5dir    nad/results/05_deseq]                            \
#       [--fai      data/sacCer3.fa.fai]                             \
#       [--outdir   nad/results/11_tracks] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
matrix_f   <- get_arg("--matrix",   "nad/results/04_ctss/matrix_reads.tsv.gz")
manifest_f <- get_arg("--manifest", "nad/scripts/00_download/priority1_runs.tsv")
s5dir      <- get_arg("--s5dir",    "nad/results/05_deseq")
fai_f      <- get_arg("--fai",      "data/sacCer3.fa.fai")
outdir     <- get_arg("--outdir",   "nad/results/11_tracks")
force      <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_manifest <- file.path(outdir, "track_manifest.tsv")
if (!force && file.exists(out_manifest)) {
  message("Output exists, skipping (use --force to rerun): ", out_manifest)
  quit(save = "no", status = 0)
}
for (f in c(matrix_f, manifest_f, fai_f)) if (!file.exists(f)) stop("Input not found: ", f)

PADJ <- 0.05; LFC <- 1
STRANDS <- c(plus = "+", minus = "-")
COLORS  <- c(S = "0,90,160", N = "140,140,140", diff = "0,0,0")
DESCR   <- c(S = "NAD captureSeq S (+ADPRC), mean of 3 replicates, reads / step-5 size factor",
             N = "NAD captureSeq N (no ADPRC), mean of 3 replicates, reads / step-5 size factor",
             diff = "NAD captureSeq S minus N, normalised as S and N")

# ── input ─────────────────────────────────────────────────────────────────────

man <- fread(manifest_f, select = c("run", "strain", "library", "replicate"))
mat <- fread(cmd = sprintf("gzip -dc '%s'", matrix_f))
stopifnot(all(man$run %in% names(mat)))
fai <- fread(fai_f, header = FALSE, select = 1:2, col.names = c("chr", "len"))
seqlens <- setNames(fai$len, fai$chr)
stopifnot(all(unique(mat$chr) %in% names(seqlens)))
si <- Seqinfo(names(seqlens), seqlens)

size_factors <- function(cnt) {                 # DESeq2 estimateSizeFactorsForMatrix
  lg <- rowMeans(log(cnt)); ok <- is.finite(lg)
  sub <- cnt[ok, , drop = FALSE]; lgo <- lg[ok]
  apply(sub, 2, function(c) exp(median((log(c) - lgo)[c > 0])))
}

write_bw <- function(dt, path) {
  gr <- GRanges(dt$chr, IRanges(dt$pos, width = 1), score = as.numeric(dt$score), seqinfo = si)
  # Fixed zoom levels, as in Stage 21: data-derived ones make bbiWrite fail on
  # sparse tracks.
  export.bw(sort(gr), path, fixedSummaries = TRUE)
}

unlink(list.files(outdir, pattern = "\\.(bw|bb)$", full.names = TRUE))
manifest <- list(); checks <- list(); sf_tab <- list()

for (st in unique(man$strain)) {
  cs <- man[strain == st]; runs <- cs$run
  s5_f <- file.path(s5dir, sprintf("%s_S_vs_N.tsv.gz", st))
  if (!file.exists(s5_f)) stop("Input not found: ", s5_f)
  s5 <- fread(cmd = sprintf("gzip -dc '%s'", s5_f))

  ctrl <- mat[s5[biotype == "none", .(chr, pos, strand)], on = .(chr, pos, strand), nomatch = NULL]
  sf <- size_factors(as.matrix(ctrl[, ..runs]))
  sf_tab[[st]] <- data.table(strain = st, run = runs, library = cs$library, size_factor = sf)
  message(sprintf("%s size factors: %s", st, paste(sprintf("%s=%.3f", runs, sf), collapse = " ")))

  nc <- sweep(as.matrix(mat[, ..runs]), 2, sf, "/")
  sig <- data.table(mat[, .(chr, pos, strand)],
                    S = rowMeans(nc[, cs$library == "S", drop = FALSE]),
                    N = rowMeans(nc[, cs$library == "N", drop = FALSE]))
  sig[, diff := S - N]

  # S and N reproduce the step-5 table at every tested position.
  cmp <- sig[s5, on = .(chr, pos, strand), nomatch = NULL]
  checks[[length(checks) + 1]] <- c(
    sprintf("%s: S and N equal step-5 meanS/meanN at %d tested positions", st, nrow(cmp)),
    nrow(cmp) == nrow(s5) && max(abs(cmp$S - cmp$meanS), abs(cmp$N - cmp$meanN)) < 1e-8)

  for (tr in c("S", "N", "diff")) for (sn in names(STRANDS)) {
    sel <- sig[strand == STRANDS[[sn]] & get(tr) != 0, .(chr, pos, score = get(tr))]
    file <- sprintf("nad_%s_%s_%s.bw", st, tr, sn)
    write_bw(sel, file.path(outdir, file))
    back <- import.bw(file.path(outdir, file))
    checks[[length(checks) + 1]] <- c(sprintf("%s read back: sum and position count match", file),
                                      abs(sum(score(back) * width(back)) - sum(sel$score)) < 1e-3 * max(1, sum(abs(sel$score))) &&
                                        sum(width(back)) == nrow(sel))
    manifest[[length(manifest) + 1]] <- data.table(
      file = file, kind = "signal", strain = st, track = tr, strand = STRANDS[[sn]],
      color = COLORS[[tr]], n_positions = nrow(sel), total = sum(sel$score),
      description = paste0(DESCR[[tr]], ", ", st))
  }

  # Step-5 calls as a bigBed, one feature per position.
  calls <- s5[!is.na(padj) & padj < PADJ & log2FoldChange > LFC]
  if (nrow(calls)) {
    gr <- GRanges(calls$chr, IRanges(calls$pos, width = 1), strand = calls$strand, seqinfo = si,
                  name = sprintf("%s_lfc%.1f", calls$base_p1, calls$log2FoldChange),
                  score = pmin(1000L, as.integer(round(100 * calls$log2FoldChange))))
    # rtracklayer writes BED9 columns in order, so thick must precede itemRgb.
    mcols(gr)$thick   <- ranges(gr)
    mcols(gr)$itemRgb <- ifelse(calls$base_p1 == "A", "#e63946", "#8c8c8c")
    file <- sprintf("nad_%s_calls.bb", st)
    ok <- tryCatch({ export(sort(gr), file.path(outdir, file), format = "bigBed"); TRUE },
                   error = function(e) { unlink(file.path(outdir, file)); message("bigBed export failed (", conditionMessage(e), "); calls track left out"); FALSE })
    if (ok) manifest[[length(manifest) + 1]] <- data.table(
      file = file, kind = "calls", strain = st, track = "calls", strand = ".",
      color = "230,57,70", n_positions = nrow(calls), total = NA_real_,
      description = sprintf("NAD captureSeq S vs N calls (padj < %.2g, log2FC > %g), red = A at +1, %s", PADJ, LFC, st))
  }
}
manifest <- rbindlist(manifest)
sf_tab <- rbindlist(sf_tab)

# ── checks ────────────────────────────────────────────────────────────────────
# The size factors must equal step 8's table, which matched the step-5 log.
sf8_f <- "nad/results/08_1Sa_positions/size_factors.tsv"
if (file.exists(sf8_f)) {
  sf8 <- fread(sf8_f)[sf_tab, on = .(strain, run)]
  checks[[length(checks) + 1]] <- c("size factors equal step 8's size_factors.tsv",
                                    isTRUE(all(abs(sf8$size_factor - sf8$i.size_factor) < 1e-9)))
}
for (ck in checks) message(sprintf("[%s] %s", if (ck[[2]] == "TRUE") "PASS" else "FAIL", ck[[1]]))
if (any(vapply(checks, function(ck) ck[[2]] != "TRUE", logical(1)))) stop("at least one check failed")

print(manifest[, .(file, n_positions, total)], nrows = 100)
fwrite(sf_tab, file.path(outdir, "size_factors.tsv"), sep = "\t")
fwrite(manifest, out_manifest, sep = "\t")   # last, so an interrupted run is not skipped
message("Written: ", out_manifest)
sessionInfo()
