#!/usr/bin/env Rscript
# Shared genomic-localization vocabulary and annotation functions.
#
# Stage 10 (scripts/06_localization/) annotates the 5' end of every alignment of
# every end type, and Stage 19 step 4 (scripts/12_exotic_origins/) annotates the
# CTSS positions each fine type dominates. The two stages ask different
# questions of the same genome, so the category set, the priority order that
# resolves overlapping features, the promoter window and the palette must be
# identical between them or the two figures cannot be read against each other.
# This file is the single source of truth for all four, exactly as
# scripts/utils/end_types.R is for the end-type vocabulary.
#
# Two annotation modes are provided. The full mode delegates to ChIPseeker,
# which is what the pipeline uses by default; the fast mode resolves the same
# categories directly against a TxDb, which needs far less memory but knows
# nothing about enhancers and calls coding exons "CDS" rather than "Exon".
#
# The caller is responsible for having loaded GenomicRanges (both modes),
# ChIPseeker (full mode only) and any TxDb it passes in.

# The promoter window used everywhere localization is computed, in bases
# upstream and downstream of an annotated TSS.
PROMOTER_UP <- 1000L
PROMOTER_DOWN <- 100L
TSS_REGION <- c(-PROMOTER_UP, PROMOTER_DOWN)

# ChIPseeker holds the whole annotation in memory at once, so long position
# vectors are annotated in batches of this size.
BATCH_SIZE <- 200000L

# Drawing order of the categories, which is also the legend order. Both modes
# are covered by one vector; a caller intersects it with the categories its own
# mode actually produced.
LOC_CATEGORY_ORDER <- c("Promoter", "Enhancer", "5' UTR", "3' UTR",
                        "CDS", "Exon", "Intron", "Downstream", "Intergenic")

# One colour per category. CDS and Exon share a colour because they are the
# same feature under the two modes' different names.
LOC_CATEGORY_COLORS <- c(
  Promoter = "#2166ac", Enhancer = "#b2182b",
  `5' UTR` = "#92c5de", `3' UTR` = "#f4a582",
  CDS = "#4daf4a", Exon = "#4daf4a", Intron = "#984ea3",
  Downstream = "#ff7f00", Intergenic = "#999999"
)

# ── Annotation functions ──────────────────────────────────────────────────────

# Fast mode: direct TxDb-based annotation
annotate_positions_fast <- function(pos_gr, promoter_gr, utr5_gr, utr3_gr, cds_gr, intron_gr) {
  n <- length(pos_gr)
  if (n == 0) return(character(0))

  category <- rep("Intergenic", n)

  if (length(intron_gr) > 0) {
    hits <- findOverlaps(pos_gr, intron_gr, ignore.strand = FALSE)
    category[queryHits(hits)] <- "Intron"
  }
  if (length(cds_gr) > 0) {
    hits <- findOverlaps(pos_gr, cds_gr, ignore.strand = FALSE)
    category[queryHits(hits)] <- "CDS"
  }
  if (length(utr3_gr) > 0) {
    hits <- findOverlaps(pos_gr, utr3_gr, ignore.strand = FALSE)
    category[queryHits(hits)] <- "3' UTR"
  }
  if (length(utr5_gr) > 0) {
    hits <- findOverlaps(pos_gr, utr5_gr, ignore.strand = FALSE)
    category[queryHits(hits)] <- "5' UTR"
  }
  hits <- findOverlaps(pos_gr, promoter_gr, ignore.strand = FALSE)
  category[queryHits(hits)] <- "Promoter"

  category
}

# Full mode: ChIPseeker annotation (batched for memory)
annotate_positions_chipseeker <- function(pos_gr, txdb, enhancer_gr = GRanges(),
                                          tssRegion = TSS_REGION, batch_size = BATCH_SIZE) {
  n <- length(pos_gr)
  if (n == 0) return(character(0))

  if (n <= batch_size) {
    anno <- annotatePeak(pos_gr, TxDb = txdb, tssRegion = tssRegion,
                         level = "transcript", annoDb = NULL, verbose = FALSE)
    category <- as.data.frame(anno)$annotation
  } else {
    n_batches <- ceiling(n / batch_size)
    category <- character(n)
    for (b in seq_len(n_batches)) {
      idx_start <- (b - 1) * batch_size + 1
      idx_end <- min(b * batch_size, n)
      anno <- annotatePeak(pos_gr[idx_start:idx_end], TxDb = txdb,
                           tssRegion = tssRegion, level = "transcript",
                           annoDb = NULL, verbose = FALSE)
      category[idx_start:idx_end] <- as.data.frame(anno)$annotation
      rm(anno); gc(verbose = FALSE)
    }
  }

  category <- gsub(" \\(.*", "", category)
  category <- gsub("Distal Intergenic", "Intergenic", category)
  category <- gsub("Downstream.*", "Downstream", category)

  if (length(enhancer_gr) > 0) {
    enh_hits <- findOverlaps(pos_gr, enhancer_gr, ignore.strand = TRUE)
    category[queryHits(enh_hits)] <- "Enhancer"
  }

  category
}

