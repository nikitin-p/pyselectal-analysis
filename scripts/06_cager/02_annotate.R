#!/usr/bin/env Rscript
# Stage 10 step 2: annotate CAGEr tag clusters with the yeast gene annotation BED.
# Strand-aware classification: sense_genic / antisense_genic / intergenic.
# Also writes per-cluster IQW statistics merged across samples.
#
# If annotation_bed_yeast is empty in params.yaml, writes unannotated clusters.
#
# Usage (from project root):
#   Rscript scripts/06_cager/02_annotate.R \
#       [--tagclusters results/cager/tagclusters.tsv.gz] \
#       [--samples     config/samples.tsv]                \
#       [--outdir      results/cager]                     \
#       [--params      config/params.yaml]                \
#       [--annot_bed   /path/to/annotation.bed]           \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
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

tc_f      <- get_arg("--tagclusters", "results/cager/tagclusters.tsv.gz")
samples_f <- get_arg("--samples",     "config/samples.tsv")
outdir    <- get_arg("--outdir",      "results/cager")
params_f  <- get_arg("--params",      "config/params.yaml")
annot_bed <- get_arg("--annot_bed",   NULL)
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_tsv  <- file.path(outdir, "tagclusters_annotated.tsv")
out_enr  <- file.path(outdir, "annotation_enrichment.tsv")

if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

if (!file.exists(tc_f)) {
  stop("Input not found: ", tc_f,
       "\nRun step 1 first (scripts/06_cager/01_cageexp.R).")
}

params <- yaml.load_file(params_f)

if (is.null(annot_bed) || annot_bed == "") {
  annot_bed <- params$annotation_bed_yeast
}

# ── load tag clusters ─────────────────────────────────────────────────────────
message("Stage 10 step 2: annotate tag clusters")
tc <- fread(cmd = paste("zcat", shQuote(tc_f)))
message("  tag clusters (all samples): ", nrow(tc))

# ── consensus cluster set (reduce across samples) ─────────────────────────────
# Union of all tag cluster ranges → strand-aware reduce
gr_all <- GRanges(
  seqnames = tc$chr,
  ranges   = IRanges(start = tc$start, end = tc$end),
  strand   = tc$strand
)
consensus <- reduce(gr_all)
message("  consensus clusters after reduce: ", length(consensus))

# Annotate consensus with mean IQW per cluster range
# Map each tag cluster to its consensus cluster by overlap
hits <- findOverlaps(gr_all, consensus, ignore.strand = FALSE)
tc[, consensus_idx := NA_integer_]
tc[queryHits(hits), consensus_idx := subjectHits(hits)]

iqw_summary <- tc[!is.na(consensus_idx) & !is.na(iq_width),
                  .(mean_iq_width = mean(iq_width, na.rm = TRUE),
                    n_samples     = .N),
                  by = consensus_idx]

cons_dt <- data.table(
  consensus_idx = seq_along(consensus),
  chr           = as.character(seqnames(consensus)),
  start         = start(consensus),
  end           = end(consensus),
  strand        = as.character(strand(consensus))
)
cons_dt <- merge(cons_dt, iqw_summary, by = "consensus_idx", all.x = TRUE)
cons_dt[, consensus_idx := NULL]

# ── BED annotation ────────────────────────────────────────────────────────────
if (!is.null(annot_bed) && annot_bed != "" && file.exists(annot_bed)) {
  message("  annotating against: ", annot_bed)

  genes <- tryCatch(
    import(annot_bed, format = "BED"),
    error = function(e) stop("Failed to import BED: ", e$message)
  )

  cons_gr <- GRanges(
    seqnames = cons_dt$chr,
    ranges   = IRanges(start = cons_dt$start, end = cons_dt$end),
    strand   = cons_dt$strand
  )

  hits_sense <- findOverlaps(cons_gr, genes, ignore.strand = FALSE)
  hits_any   <- findOverlaps(cons_gr, genes, ignore.strand = TRUE)
  hits_anti  <- hits_any[!(queryHits(hits_any) %in% queryHits(hits_sense))]

  annotation <- rep("intergenic", length(cons_gr))
  annotation[queryHits(hits_anti)]  <- "antisense_genic"
  annotation[queryHits(hits_sense)] <- "sense_genic"

  dist_nearest <- rep(NA_integer_, length(cons_gr))
  nearest_gene <- rep(NA_character_, length(cons_gr))

  intergenic_idx <- which(annotation == "intergenic")
  if (length(intergenic_idx) > 0) {
    nn <- distanceToNearest(cons_gr[intergenic_idx], genes, ignore.strand = TRUE)
    dist_nearest[intergenic_idx[queryHits(nn)]] <- as.integer(mcols(nn)$distance)
    if (!is.null(genes$name)) {
      nearest_gene[intergenic_idx[queryHits(nn)]] <- genes$name[subjectHits(nn)]
    }
  }

  cons_dt$annotation   <- annotation
  cons_dt$dist_nearest <- dist_nearest
  cons_dt$nearest_gene <- nearest_gene

  # Enrichment table
  enr <- cons_dt[, .N, by = annotation]
  enr[, pct := round(100 * N / sum(N), 2)]
  message("  annotation counts:")
  print(enr)
  fwrite(enr, out_enr, sep = "\t")
  message("  written: ", out_enr)

} else {
  if (!is.null(annot_bed) && annot_bed != "") {
    warning("Annotation BED not found: ", annot_bed, " — skipping annotation")
  } else {
    message("  no annotation BED provided (set annotation_bed_yeast in params.yaml)")
  }
  cons_dt$annotation   <- NA_character_
  cons_dt$dist_nearest <- NA_integer_
  cons_dt$nearest_gene <- NA_character_
}

# ── write output ──────────────────────────────────────────────────────────────
fwrite(cons_dt, out_tsv, sep = "\t")
message("Written: ", out_tsv, "  (", nrow(cons_dt), " consensus clusters)")

sessionInfo()
