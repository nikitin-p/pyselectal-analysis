#!/usr/bin/env Rscript
# Stage 06: Genomic localization by 5' end type
#
# Uses ChIPseeker for genomic annotation (default) or TxDb-based fast mode.
# CAGEfightR for enhancer detection (requires pooled signal, memory-intensive).
#
# ChIPseeker priority: Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic
# Fast mode priority:  Promoter > 5'UTR > 3'UTR > CDS > Intron > Intergenic
#
# Usage:
#   Rscript scripts/06_localization/01_localization.R \
#       --samples config/samples.tsv \
#       --params  config/params.yaml \
#       --bam_sel results/bam_selected \
#       --bam_sub results/bam_subsampled \
#       --gtf     data/sacCer3.ensGene.gtf \
#       --outdir  results/localization \
#       --figdir  results/figures \
#       [--fast]   # Use TxDb-based annotation (low memory)
#       [--force]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(Rsamtools)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
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
params_f  <- get_arg("--params",  "config/params.yaml")
bam_sel   <- get_arg("--bam_sel", "results/bam_selected")
bam_sub   <- get_arg("--bam_sub", "results/bam_subsampled")
gtf_f     <- get_arg("--gtf",     "data/sacCer3.ensGene.gtf")
outdir    <- get_arg("--outdir",  "results/localization")
figdir    <- get_arg("--figdir",  "results/figures")
force     <- "--force" %in% args
fast_mode <- "--fast" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# fmt_n(): compact count labels drawn on top of every histogram bar
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

params  <- yaml.load_file(params_f)
samples <- read_tsv(samples_f, show_col_types = FALSE) %>%
  filter(!grepl("^#", sample_id))

PROMOTER_UP <- 1000L
PROMOTER_DOWN <- 100L
TSS_REGION <- c(-PROMOTER_UP, PROMOTER_DOWN)
BALANCE_THRESHOLD <- 0.95
BATCH_SIZE <- 200000L

# ── Load annotation-related packages ──────────────────────────────────────────
if (!fast_mode) {
  message("Loading ChIPseeker and CAGEfightR for full annotation mode...")
  suppressPackageStartupMessages({
    library(CAGEfightR)
    library(ChIPseeker)
    library(SummarizedExperiment)
  })
} else {
  message("Using fast mode (TxDb-based annotation)")
}

# ── Load TxDb ─────────────────────────────────────────────────────────────────
message("Loading gene annotation: ", gtf_f)
txdb <- makeTxDbFromGFF(gtf_f, format = "gtf")

if (fast_mode) {
  message("Building annotation layers for fast mode...")
  tx <- transcripts(txdb)
  promoter_gr <- promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN)
  utr5_gr <- unlist(fiveUTRsByTranscript(txdb))
  utr3_gr <- unlist(threeUTRsByTranscript(txdb))
  cds_gr <- unlist(cdsBy(txdb, by = "tx"))
  intron_gr <- unlist(intronsByTranscript(txdb))
  message("  promoters: ", length(promoter_gr), ", 5'UTR: ", length(utr5_gr),
          ", 3'UTR: ", length(utr3_gr), ", CDS: ", length(cds_gr), ", introns: ", length(intron_gr))
}

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

# ── CAGEfightR enhancer detection ─────────────────────────────────────────────
detect_enhancers_cagefightr <- function(ctss_list, promoter_gr, balance_threshold = BALANCE_THRESHOLD) {
  message("Pooling CTSS for enhancer detection...")

  pooled_ctss <- do.call(c, unname(ctss_list))
  pooled_dt <- data.table(
    chr    = as.character(seqnames(pooled_ctss)),
    pos    = start(pooled_ctss),
    strand = as.character(strand(pooled_ctss)),
    score  = pooled_ctss$score
  )
  pooled_agg <- pooled_dt[, .(score = sum(score)), by = .(chr, pos, strand)]

  pooled_gr <- GRanges(
    seqnames = pooled_agg$chr,
    ranges   = IRanges(start = pooled_agg$pos, width = 1),
    strand   = pooled_agg$strand,
    score    = pooled_agg$score
  )
  message("  pooled CTSS: ", length(pooled_gr))

  message("Detecting bidirectional clusters (balance threshold = ", balance_threshold, ")...")
  bidirectional_clusters <- clusterBidirectionally(pooled_gr, balanceThreshold = balance_threshold)
  message("  bidirectional clusters: ", length(bidirectional_clusters))

  if (length(bidirectional_clusters) == 0) {
    return(GRanges())
  }

  promoter_hits <- findOverlaps(bidirectional_clusters, promoter_gr)
  if (length(promoter_hits) > 0) {
    enhancer_gr <- bidirectional_clusters[-unique(queryHits(promoter_hits))]
  } else {
    enhancer_gr <- bidirectional_clusters
  }
  message("  enhancers (excluding promoters): ", length(enhancer_gr))

  enhancer_gr
}

# ── BAM reading functions ─────────────────────────────────────────────────────
bam_to_5prime <- function(bam_path) {
  param <- ScanBamParam(
    flag = scanBamFlag(
      isSecondaryAlignment     = FALSE,
      isSupplementaryAlignment = FALSE,
      isUnmappedQuery          = FALSE
    )
  )
  ga <- readGAlignments(bam_path, param = param)
  if (length(ga) == 0) return(GRanges())

  strand_chr <- as.character(strand(ga))
  pos_idx <- strand_chr %in% c("+", "*")
  neg_idx <- strand_chr == "-"

  pos_5p <- integer(length(ga))
  pos_5p[pos_idx] <- start(ga)[pos_idx]
  pos_5p[neg_idx] <- end(ga)[neg_idx]

  strand_out <- strand_chr
  strand_out[strand_out == "*"] <- "+"

  GRanges(
    seqnames = seqnames(ga),
    ranges   = IRanges(start = pos_5p, width = 1),
    strand   = strand_out
  )
}

bam_to_5prime_M <- function(bam_path) {
  param <- ScanBamParam(
    flag = scanBamFlag(
      isSecondaryAlignment     = FALSE,
      isSupplementaryAlignment = FALSE,
      isUnmappedQuery          = FALSE
    ),
    what = "cigar"
  )
  ga <- readGAlignments(bam_path, param = param, use.names = FALSE)
  if (length(ga) == 0) return(GRanges())

  cigar <- mcols(ga)$cigar
  strand_chr <- as.character(strand(ga))

  is_plus  <- strand_chr %in% c("+", "*")
  is_minus <- strand_chr == "-"

  is_M <- logical(length(ga))
  is_M[is_plus]  <- !grepl("^\\d+S", cigar[is_plus])
  is_M[is_minus] <- !grepl("\\d+S$", cigar[is_minus])

  ga_M <- ga[is_M]
  if (length(ga_M) == 0) return(GRanges())

  strand_M <- as.character(strand(ga_M))
  pos_idx <- strand_M %in% c("+", "*")
  neg_idx <- strand_M == "-"

  pos_5p <- integer(length(ga_M))
  pos_5p[pos_idx] <- start(ga_M)[pos_idx]
  pos_5p[neg_idx] <- end(ga_M)[neg_idx]

  strand_out <- strand_M
  strand_out[strand_out == "*"] <- "+"

  GRanges(
    seqnames = seqnames(ga_M),
    ranges   = IRanges(start = pos_5p, width = 1),
    strand   = strand_out
  )
}

bam_to_ctss <- function(bam_path) {
  pos_gr <- bam_to_5prime(bam_path)
  if (length(pos_gr) == 0) return(GRanges())

  dt <- data.table(
    chr    = as.character(seqnames(pos_gr)),
    pos    = start(pos_gr),
    strand = as.character(strand(pos_gr))
  )
  agg <- dt[, .(score = .N), by = .(chr, pos, strand)]

  GRanges(
    seqnames = agg$chr,
    ranges   = IRanges(start = agg$pos, width = 1),
    strand   = agg$strand,
    score    = agg$score
  )
}

# ── Enhancer detection (full mode only) ───────────────────────────────────────
enhancer_gr <- GRanges()

if (!fast_mode) {
  message("\nBuilding CTSS for enhancer detection...")
  ctss_list <- list()
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    bam_path <- file.path(bam_sub, paste0(sid, "_sub.bam"))
    if (!file.exists(bam_path)) next
    message("  ", sid)
    ctss_list[[sid]] <- bam_to_ctss(bam_path)
  }

  if (length(ctss_list) > 0) {
    tx <- transcripts(txdb)
    prom_gr <- promoters(tx, upstream = PROMOTER_UP, downstream = PROMOTER_DOWN)

    tryCatch({
      enhancer_gr <- detect_enhancers_cagefightr(ctss_list, prom_gr, BALANCE_THRESHOLD)
      if (length(enhancer_gr) > 0) {
        enh_bed <- file.path(outdir, "enhancer_regions.bed")
        export(enhancer_gr, enh_bed, format = "BED")
        message("  written: ", enh_bed)
      }
    }, error = function(e) {
      message("  CAGEfightR enhancer detection failed: ", e$message)
      enhancer_gr <<- GRanges()
    })
  }

  rm(ctss_list); gc(verbose = FALSE)
}

# ── Process all samples ───────────────────────────────────────────────────────
message("\nProcessing samples...")

end_types <- c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M")
results_list <- list()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]
  message("  ", sample_id)

  for (et in end_types) {
    # All types now have pre-selected BAM files from 00_select_bam.sh
    bam_path <- file.path(bam_sel, paste0(sample_id, "_", et, ".bam"))

    if (!file.exists(bam_path)) {
      message("    [", et, "] BAM not found: ", bam_path)
      next
    }

    pos_gr <- bam_to_5prime(bam_path)

    if (length(pos_gr) == 0) next

    if (fast_mode) {
      category <- annotate_positions_fast(pos_gr, promoter_gr, utr5_gr, utr3_gr, cds_gr, intron_gr)
    } else {
      category <- annotate_positions_chipseeker(pos_gr, txdb, enhancer_gr)
    }

    # Aggregate counts
    dt <- data.table(category = category)
    agg <- dt[, .(n_reads = .N), by = category]

    pos_dt <- data.table(
      chr      = as.character(seqnames(pos_gr)),
      pos      = start(pos_gr),
      strand   = as.character(strand(pos_gr)),
      category = category
    )
    pos_uniq <- unique(pos_dt, by = c("chr", "pos", "strand", "category"))
    agg_pos <- pos_uniq[, .(n_positions = .N), by = category]

    merged <- merge(agg, agg_pos, by = "category", all = TRUE)
    merged[is.na(n_reads), n_reads := 0L]
    merged[is.na(n_positions), n_positions := 0L]
    merged$sample_id <- sample_id
    merged$end_type  <- et

    results_list[[paste(sample_id, et, sep = "_")]] <- merged

    rm(pos_gr, category, dt, agg, pos_dt, pos_uniq, agg_pos, merged)
    gc(verbose = FALSE)
  }
}

results_dt <- rbindlist(results_list)
setcolorder(results_dt, c("sample_id", "end_type", "category", "n_reads", "n_positions"))

# ── Write outputs ─────────────────────────────────────────────────────────────
counts_f <- file.path(outdir, "category_counts.tsv")
fwrite(results_dt, counts_f, sep = "\t")
message("\nWritten: ", counts_f)

summary_dt <- results_dt[, .(
  total_reads     = sum(n_reads),
  total_positions = sum(n_positions)
), by = .(end_type, category)]

summary_dt[, pct_reads := 100 * total_reads / sum(total_reads), by = end_type]
summary_dt[, pct_positions := 100 * total_positions / sum(total_positions), by = end_type]

summary_f <- file.path(outdir, "category_summary.tsv")
fwrite(summary_dt, summary_f, sep = "\t")
message("Written: ", summary_f)

# ── Chi-square tests ──────────────────────────────────────────────────────────
categories_present <- unique(summary_dt$category)
message("  categories: ", paste(categories_present, collapse = ", "))

ref_reads <- summary_dt[end_type == "1Sg", .(category, total_reads)]
ref_prop <- setNames(ref_reads$total_reads / sum(ref_reads$total_reads), ref_reads$category)

chisq_results <- list()
for (et in intersect(c("2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M"), unique(summary_dt$end_type))) {
  obs <- summary_dt[end_type == et, .(category, total_reads)]
  common <- intersect(names(ref_prop), obs$category)
  if (length(common) < 2) next

  obs_vec <- setNames(obs$total_reads, obs$category)[common]
  ref_vec <- ref_prop[common]
  ref_vec <- ref_vec / sum(ref_vec)

  if (all(ref_vec * sum(obs_vec) > 5)) {
    chi <- chisq.test(obs_vec, p = ref_vec)
    chisq_results[[et]] <- data.table(
      end_type = et, chi_sq = chi$statistic, df = chi$parameter, p_value = chi$p.value
    )
  }
}
chisq_dt <- rbindlist(chisq_results)
chisq_f <- file.path(outdir, "chisq_results.tsv")
fwrite(chisq_dt, chisq_f, sep = "\t")
message("Written: ", chisq_f)

# ── Visualization ─────────────────────────────────────────────────────────────
message("\nGenerating figures...")

if (fast_mode) {
  cat_order <- c("Promoter", "5' UTR", "3' UTR", "CDS", "Intron", "Intergenic")
} else {
  cat_order <- c("Promoter", "Enhancer", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Intergenic")
}
cat_order <- intersect(cat_order, unique(summary_dt$category))

cat_colors <- c(
  Promoter = "#2166ac", Enhancer = "#b2182b",
  `5' UTR` = "#92c5de", `3' UTR` = "#f4a582",
  CDS = "#4daf4a", Exon = "#4daf4a", Intron = "#984ea3",
  Downstream = "#ff7f00", Intergenic = "#999999"
)

plot_dt <- melt(summary_dt, id.vars = c("end_type", "category"),
                measure.vars = c("pct_reads", "pct_positions"),
                variable.name = "level", value.name = "pct")
plot_dt[, level := gsub("pct_", "", level)]
plot_dt[, category := factor(category, levels = cat_order)]
plot_dt[, end_type := factor(end_type, levels = c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M"))]

# Totals behind each 100% bar: reads in the left facet, unique positions in the right.
bar_n <- rbind(
  summary_dt[!is.na(category), .(level = "reads", n = sum(total_reads)), by = end_type],
  summary_dt[!is.na(category), .(level = "positions", n = sum(total_positions)), by = end_type]
)
bar_n[, end_type := factor(end_type, levels = levels(plot_dt$end_type))]

p_bar <- ggplot(plot_dt[!is.na(category)], aes(x = end_type, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(data = bar_n, aes(x = end_type, y = 101, label = fmt_n(n)),
            inherit.aes = FALSE, vjust = 0, size = 2.6, colour = "grey25") +
  facet_wrap(~level, labeller = labeller(level = c(reads = "Alignment counts", positions = "Unique positions"))) +
  scale_fill_manual(values = cat_colors, name = "Category", drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 108)) +
  labs(title = paste0("Genomic localization by 5' end type", if(fast_mode) " (fast mode)" else " (ChIPseeker)"),
       subtitle = "Labels give the number of unique positions or alignments behind each bar",
       x = "End type", y = "Percentage") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

ggsave(file.path(figdir, "10_localization_barplot.pdf"), p_bar, width = 10, height = 6)
message("  written: 10_localization_barplot.pdf")

# Enrichment heatmap
enrich_dt <- summary_dt[, .(end_type, category, pct_reads)]
ref_pct <- enrich_dt[end_type == "1Sg", .(category, ref = pct_reads)]
enrich_dt <- merge(enrich_dt, ref_pct, by = "category")
enrich_dt[, log2_enrich := log2((pct_reads + 0.01) / (ref + 0.01))]

enrich_wide <- dcast(enrich_dt, end_type ~ category, value.var = "log2_enrich")
enrich_mat <- as.matrix(enrich_wide[, -1, with = FALSE])
rownames(enrich_mat) <- enrich_wide$end_type

avail_types <- intersect(c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M"), rownames(enrich_mat))
avail_cats <- intersect(cat_order, colnames(enrich_mat))
enrich_mat <- enrich_mat[avail_types, avail_cats, drop = FALSE]

pdf(file.path(figdir, "10_localization_enrichment.pdf"), width = 8, height = 4)
pheatmap(enrich_mat, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
         breaks = seq(-2, 2, length.out = 51),
         main = "Log2 enrichment vs 1Sg", display_numbers = TRUE, number_format = "%.2f")
dev.off()
message("  written: 10_localization_enrichment.pdf")

# Per-sample plot
# Percentages are taken within a sample x end_type bar, so category must stay in j:
# grouping by it as well makes every bar segment 100%.
sample_dt <- results_dt[, .(category, pct = 100 * n_reads / sum(n_reads)),
                        by = .(sample_id, end_type)]
sample_dt[, category := factor(category, levels = cat_order)]
sample_dt[, end_type := factor(end_type, levels = c("1Sg", "2Sgg", "3Sggg", "1Sa", "1St", "1Sc", "M"))]

sample_n <- results_dt[, .(n = sum(n_reads)), by = .(sample_id, end_type)]
sample_n[, end_type := factor(end_type, levels = levels(sample_dt$end_type))]

p_sample <- ggplot(sample_dt[!is.na(category)], aes(x = end_type, y = pct, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(data = sample_n, aes(x = end_type, y = 101, label = fmt_n(n)),
            inherit.aes = FALSE, vjust = 0, size = 2, colour = "grey25") +
  facet_wrap(~sample_id, ncol = 4) +
  scale_fill_manual(values = cat_colors, name = "Category", drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 110)) +
  labs(title = "Localization by sample",
       subtitle = "Labels give the number of alignments per end type",
       x = "End type", y = "Percentage") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figdir, "10_localization_by_sample.pdf"), p_sample, width = 12, height = 8)
message("  written: 10_localization_by_sample.pdf")

sessionInfo()
