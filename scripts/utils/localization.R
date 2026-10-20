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
#
# "Promoter" itself has two definitions, chosen by promoter_source in the
# params file. Under "window" it is the fixed window around every annotated
# transcript start. Under "clusters" it is the set of TSS clusters the data
# show upstream of a gene on the gene's own strand, written by
# scripts/05_tss_end_types/12_promoter_clusters.R, and the window plays no part
# in the labelling.

# The promoter window used everywhere localization is computed, in bases
# upstream and downstream of an annotated TSS. The values below are only the
# fallback: the window is a property of the genome, not of the pipeline, so each
# params file sets promoter_upstream / promoter_downstream and every caller
# applies them with set_promoter_window() right after loading its params. Yeast
# needs a far narrower upstream arm than human, because the median gap between
# a yeast start codon and the gene upstream of it is under 400 bp and a 1 kb arm
# would label the 3' end of that neighbouring gene as a promoter.
PROMOTER_UP <- 1000L
PROMOTER_DOWN <- 100L
TSS_REGION <- c(-PROMOTER_UP, PROMOTER_DOWN)

# Overwrites the three globals above from a loaded params list. A key that is
# absent leaves the fallback in place.
set_promoter_window <- function(params) {
  if (!is.null(params$promoter_upstream))
    PROMOTER_UP <<- as.integer(params$promoter_upstream)
  if (!is.null(params$promoter_downstream))
    PROMOTER_DOWN <<- as.integer(params$promoter_downstream)
  if (is.na(PROMOTER_UP) || is.na(PROMOTER_DOWN) || PROMOTER_UP < 0L || PROMOTER_DOWN < 0L)
    stop("promoter_upstream and promoter_downstream must be non-negative integers")
  TSS_REGION <<- c(-PROMOTER_UP, PROMOTER_DOWN)
  message("Promoter window: -", PROMOTER_UP, " .. +", PROMOTER_DOWN)
  invisible(TSS_REGION)
}

# Returns the cluster-defined promoters as a stranded GRanges when
# promoter_source is "clusters", and NULL under "window" or when the key is
# absent, which is how a caller tells the two definitions apart. A missing table
# is an error rather than a silent fallback to the window, because the two
# definitions give different figures under the same file names.
load_promoter_clusters <- function(params, path = NULL) {
  src <- if (is.null(params$promoter_source)) "window" else params$promoter_source
  if (!src %in% c("window", "clusters"))
    stop("promoter_source must be 'window' or 'clusters', not '", src, "'")
  if (src == "window") return(NULL)
  if (is.null(path)) {
    root <- if (is.null(params$results_root)) "results" else params$results_root
    path <- file.path(root, "tss", "promoter_clusters.tsv")
  }
  if (!file.exists(path))
    stop("promoter_source is 'clusters' but ", path, " does not exist; run ",
         "scripts/05_tss_end_types/12_promoter_clusters.R first")
  tab <- unique(utils::read.delim(path, stringsAsFactors = FALSE)[
    , c("chr", "start", "end", "strand")])
  message("Promoters: ", nrow(tab), " TSS clusters from ", path)
  GRanges(tab$chr, IRanges(tab$start, tab$end), strand = tab$strand)
}

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
#
# With promoter_gr supplied (see load_promoter_clusters), ChIPseeker is run with
# an empty promoter window, so every position first takes the category it would
# have with no promoter at all, and a position inside a promoter cluster on its
# own strand is then relabelled "Promoter". ChIPseeker still calls the annotated
# start base itself "Promoter" under an empty window; that base is the first
# transcribed one, so it is given "5' UTR" where the annotation has one over it
# and "Exon" otherwise.
annotate_positions_chipseeker <- function(pos_gr, txdb, enhancer_gr = GRanges(),
                                          tssRegion = TSS_REGION, batch_size = BATCH_SIZE,
                                          promoter_gr = NULL) {
  n <- length(pos_gr)
  if (n == 0) return(character(0))
  if (!is.null(promoter_gr)) tssRegion <- c(0L, 0L)

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

  if (!is.null(promoter_gr)) {
    at_start <- which(category == "Promoter")
    if (length(at_start) > 0) {
      utr5_gr <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb))
      category[at_start] <- ifelse(
        overlapsAny(pos_gr[at_start], utr5_gr, ignore.strand = TRUE), "5' UTR", "Exon")
    }
    category[overlapsAny(pos_gr, promoter_gr, ignore.strand = FALSE)] <- "Promoter"
  }

  if (length(enhancer_gr) > 0) {
    enh_hits <- findOverlaps(pos_gr, enhancer_gr, ignore.strand = TRUE)
    category[queryHits(enh_hits)] <- "Enhancer"
  }

  category
}

