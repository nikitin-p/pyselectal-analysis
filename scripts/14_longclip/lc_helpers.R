#!/usr/bin/env Rscript
# Stage 22 helpers: the origin classes of a long 5' soft clip, their palette,
# and the pure functions that place a clip in the genome. Nothing here touches
# a BAM or the filesystem.

# ── origin classes ────────────────────────────────────────────────────────────
# Ordered from the most to the least constrained placement. Every long-clip read
# takes exactly one of these.
LC_CLASSES <- c("splice_annotated", "splice_unannotated", "templated_adjacent",
                "templated_distal", "unplaced")

LC_CLASS_COLORS <- c(
  splice_annotated   = "#1d4ed8",
  splice_unannotated = "#7c3aed",
  templated_adjacent = "#059669",
  templated_distal   = "#d97706",
  unplaced           = "#dc2626"
)

# The TSS context of a read's inferred 5' end.
LC_TSS_CLASSES <- c("canonical", "alternative_promoter", "other_cluster", "no_cluster")
LC_TSS_COLORS  <- c(canonical = "#1f2937", alternative_promoter = "#6b7280",
                    other_cluster = "#d1d5db", no_cluster = "#fde68a")

# ── sequence helpers ──────────────────────────────────────────────────────────

# Fraction of identical positions between equal-length strings, vectorised.
identity_frac <- function(a, b) {
  stopifnot(length(a) == length(b))
  vapply(seq_along(a), function(i) {
    x <- strsplit(a[i], "")[[1]]; y <- strsplit(b[i], "")[[1]]
    if (length(x) == 0 || length(x) != length(y)) return(NA_real_)
    mean(x == y)
  }, numeric(1))
}

# Genome sequence of [start, end] read 5'->3' on `strand`, vectorised over rows.
# Out-of-range intervals return NA rather than failing the whole batch.
strand_seq <- function(fa, chr, start, end, strand, seqlens) {
  out <- rep(NA_character_, length(chr))
  ok  <- start >= 1 & end <= seqlens[chr] & end >= start
  if (any(ok)) {
    gr <- GenomicRanges::GRanges(chr[ok], IRanges::IRanges(start[ok], end[ok]),
                                 strand = strand[ok])
    out[ok] <- as.character(Rsamtools::getSeq(fa, gr))
  }
  out
}

# Convert a 1-based index into a read-strand window string back to a genomic
# coordinate. On + the window runs left to right from `wstart`; on - it is the
# reverse complement of [wstart, wend], so index 1 is wend.
win_to_genome <- function(idx, strand, wstart, wend) {
  ifelse(strand == "+", wstart + idx - 1L, wend - idx + 1L)
}

# First shift t in 0..max_shift at which the gap, slid t bases toward the read,
# reads GT...AG on the read strand. `ext` is the read-strand genome string from
# the first gap base through max_shift bases past the aligned 5' end, and `gap`
# is the gap length. Returns NA when no shift gives the motif. The slide is what
# recovers a junction the aligner shifted by extending the aligned part over
# intron bases that happen to equal the exon-1 tail.
splice_shift <- function(ext, gap, max_shift) {
  vapply(seq_along(ext), function(i) {
    if (is.na(ext[i]) || gap[i] < 4L) return(NA_integer_)
    for (t in 0:max_shift) {
      s <- 1L + t; e <- gap[i] + t
      if (e > nchar(ext[i])) break
      if (substr(ext[i], s, s + 1L) == "GT" && substr(ext[i], e - 1L, e) == "AG")
        return(as.integer(t))
    }
    NA_integer_
  }, integer(1))
}
