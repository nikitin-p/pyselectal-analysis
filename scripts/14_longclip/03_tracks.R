#!/usr/bin/env Rscript
# Stage 22 step 3: browser tracks of the long clips, coloured by origin class.
#
#   longclip_reads.bed    BED12, one item per distinct alignment: the clip as one
#                         block at its origin and the aligned part as the other
#                         block(s), joined by the browser's intron line when they
#                         are apart. The clip block is the thick part, so it is
#                         drawn taller than the aligned part; the item colour is
#                         the origin class; the name is clip|class|n_reads; the
#                         score is n_reads capped at 1000. An unplaced clip is
#                         drawn on the bases immediately upstream of the aligned
#                         5' end, where an untemplated clip would sit. A placed
#                         clip's untemplated leading bases are not drawn, because
#                         they have no genomic position; the name keeps them.
#   longclip_reads.bb     the same, when --bedtobigbed names UCSC's bedToBigBed
#   clip_<class>_{plus,minus}.bw    coverage of the drawn clip bases, read-weighted
#   start_<class>_{plus,minus}.bw   aligned 5' ends of the class's reads
#   track_manifest.tsv
#
# The BED is written by hand rather than through rtracklayer::export.bb, which
# mis-encodes BED12 block offsets; export.bw was checked separately and is
# correct. Every item's blocks are asserted to tile chromStart..chromEnd in
# order, which is the check bedToBigBed would otherwise make.
#
# Usage (from project root):
#   Rscript scripts/14_longclip/03_tracks.R \
#       [--reads  results/longclip/longclip_reads.tsv.gz] \
#       [--fai    data/sacCer3.fa.fai] \
#       [--outdir results/longclip/tracks] [--bedtobigbed PATH] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicAlignments)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
reads_f <- get_arg("--reads",       "results/longclip/longclip_reads.tsv.gz")
fai_f   <- get_arg("--fai",         "data/sacCer3.fa.fai")
outdir  <- get_arg("--outdir",      "results/longclip/tracks")
b2bb    <- get_arg("--bedtobigbed", "")
force   <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "lc_helpers.R"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_f <- file.path(outdir, "track_manifest.tsv")
if (!force && file.exists(manifest_f)) {
  message("Output exists, skipping (use --force to rerun): ", manifest_f)
  quit(save = "no", status = 0)
}

reads <- fread(reads_f)
fai <- fread(fai_f, header = FALSE, select = 1:2, col.names = c("chr", "len"))
seqlens <- setNames(fai$len, fai$chr)
si <- Seqinfo(names(seqlens), seqlens)

# ── the clip block ────────────────────────────────────────────────────────────

placed <- reads$class != "unplaced"
reads[, `:=`(
  clip_start = ifelse(placed, tmpl_start, ifelse(strand == "+", aln_start - clen, aln_end + 1L)),
  clip_end   = ifelse(placed, tmpl_end,   ifelse(strand == "+", aln_start - 1L,   aln_end + clen)))]
off_chrom <- reads[clip_start < 1 | clip_end > seqlens[chr]]
if (nrow(off_chrom)) message("Dropped ", nrow(off_chrom),
                             " alignments whose clip block runs off the chromosome")
reads <- reads[clip_start >= 1 & clip_end <= seqlens[chr]]

# ── BED12 ─────────────────────────────────────────────────────────────────────

ga <- GAlignments(seqnames = Rle(reads$chr), pos = reads$aln_start,
                  cigar = reads$cigar, strand = Rle(strand(reads$strand)))
aln_blocks <- as.data.table(grglist(ga))[, .(item = group, start, end)]
blocks <- rbind(aln_blocks,
                reads[, .(item = .I, start = clip_start, end = clip_end)])
setorder(blocks, item, start)
# Abutting blocks (a clip templated with no gap) are merged, because BED12
# blocks must not touch; the thick part still marks the clip.
blocks[, grp := cumsum(c(TRUE, start[-1] > end[-.N] + 1L)), by = item]
blocks <- blocks[, .(start = min(start), end = max(end)), by = .(item, grp)]

bed <- blocks[, .(chromStart = min(start) - 1L, chromEnd = max(end),
                  blockCount = .N,
                  blockSizes  = paste0(paste(end - start + 1L, collapse = ","), ","),
                  blockStarts = paste0(paste(start - min(start), collapse = ","), ","),
                  tiles = all(diff(start) > 0) && min(start) == start[1] &&
                          max(end) == end[.N] && all(start[-1] > end[-.N])),
              by = item]
stopifnot(all(bed$tiles), nrow(bed) == nrow(reads))

hex_rgb <- function(h) apply(col2rgb(h), 2, paste, collapse = ",")
bed[, `:=`(chrom = reads$chr[item],
           name = paste(reads$clip[item], reads$class[item], reads$n_reads[item], sep = "|"),
           score = pmin(reads$n_reads[item], 1000L),
           strand = reads$strand[item],
           thickStart = reads$clip_start[item] - 1L,
           thickEnd = reads$clip_end[item],
           itemRgb = hex_rgb(LC_CLASS_COLORS[reads$class[item]]))]
stopifnot(all(bed$thickStart >= bed$chromStart & bed$thickEnd <= bed$chromEnd))
setorder(bed, chrom, chromStart)
bed_f <- file.path(outdir, "longclip_reads.bed")
fwrite(bed[, .(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd,
               itemRgb, blockCount, blockSizes, blockStarts)],
       bed_f, sep = "\t", col.names = FALSE, quote = FALSE)
message("Written: ", bed_f, " (", nrow(bed), " items)")

if (nzchar(b2bb)) {
  sizes <- tempfile(fileext = ".sizes")
  fwrite(fai, sizes, sep = "\t", col.names = FALSE)
  bb_f <- file.path(outdir, "longclip_reads.bb")
  status <- system2(b2bb, c("-type=bed12", bed_f, sizes, bb_f))
  if (status != 0) stop("bedToBigBed failed with status ", status)
  message("Written: ", bb_f)
} else {
  message("No --bedtobigbed given; the bigBed was not written")
}

# ── bigWigs ───────────────────────────────────────────────────────────────────

unlink(list.files(outdir, pattern = "\\.bw$", full.names = TRUE))
write_cov <- function(gr, weight, path) {
  if (!length(gr)) return(FALSE)
  cv <- coverage(gr, weight = weight)
  out <- as(cv, "GRanges")
  out <- out[out$score > 0]
  seqinfo(out) <- si[seqlevels(out)]
  # Fixed zoom levels, as in Stage 21: UCSC's bbiWrite fails on data-derived
  # ones when a track holds a handful of positions on different chromosomes.
  export.bw(out, path, fixedSummaries = TRUE)
  TRUE
}

strands <- c(plus = "+", minus = "-")
manifest <- rbindlist(lapply(LC_CLASSES, function(cl) rbindlist(lapply(names(strands), function(sn) {
  r <- reads[class == cl & strand == strands[[sn]]]
  rows <- list()
  for (kind in c("clip", "start")) {
    gr <- if (kind == "clip")
      GRanges(r$chr, IRanges(r$clip_start, r$clip_end), seqinfo = si) else
      GRanges(r$chr, IRanges(r$pos5, width = 1), seqinfo = si)
    file <- paste0(kind, "_", cl, "_", sn, ".bw")
    ok <- write_cov(gr, r$n_reads, file.path(outdir, file))
    rows[[kind]] <- data.table(
      file = if (ok) file else NA_character_, kind = kind, class = cl,
      strand = strands[[sn]], color = hex_rgb(LC_CLASS_COLORS[[cl]]),
      n_reads = sum(r$n_reads),
      signal = if (kind == "clip") sum(r$n_reads * (r$clip_end - r$clip_start + 1L)) else sum(r$n_reads))
  }
  rbindlist(rows)
}))))

written  <- manifest[!is.na(file)]
readback <- vapply(written$file, function(f) {
  g <- import.bw(file.path(outdir, f)); sum(score(g) * width(g))
}, numeric(1))
bad <- written$file[abs(readback - written$signal) > 0.5]
if (length(bad)) stop("bigWig sums disagree with the manifest: ", paste(bad, collapse = ", "))
message("All ", nrow(written), " bigWigs read back with sums matching the manifest")
print(manifest[, .(file, n_reads, signal)])

fwrite(manifest, manifest_f, sep = "\t")
message("Written: ", manifest_f)
sessionInfo()
