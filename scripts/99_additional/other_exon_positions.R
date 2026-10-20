#!/usr/bin/env Rscript
# List the exonic CTSS positions that carry "other" 5' end types above an
# expression threshold.
#
# "Other" is the pooled tail bar of the Stage 19 figures: every fine type that
# is outside the vocabulary of scripts/utils/end_types.R and is not longclip.
# The expression of a position is the number of alignments of those types at
# it, pooled over the libraries, which were subsampled to equal depth in
# Stage 2. A position is listed when that count is ABOVE --min_reads, and the
# above_10 column marks the stricter of the two thresholds that were asked for.
#
# The genomic category comes from the same ChIPseeker call, promoter window and
# priority order as Stage 19 steps 4 and 6 (scripts/utils/localization.R), so a
# position listed here is one that those figures draw as "Exon".
#
# Usage (from project root):
#   Rscript scripts/99_additional/other_exon_positions.R \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
#       [--ctss      results/dominance/fine_typed_ctss.tsv.gz]          \
#       [--gtf       data/sacCer3.ensGene.gtf]                         \
#       [--out       results/exotic_origins/other_exon_positions.tsv]  \
#       [--min_reads 5]                                               \
#       [--promoters results/tss/promoter_clusters.tsv]                \
#       [--params    config/params.yaml]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
dom_f     <- get_arg("--dominance", "results/exotic_origins/position_dominance.tsv.gz")
ctss_f    <- get_arg("--ctss",      "results/dominance/fine_typed_ctss.tsv.gz")
gtf_f     <- get_arg("--gtf",       "data/sacCer3.ensGene.gtf")
out_f     <- get_arg("--out",       "results/exotic_origins/other_exon_positions.tsv")
min_reads <- as.integer(get_arg("--min_reads", "5"))
params_f  <- get_arg("--params",    "config/params.yaml")

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "..", "utils", "localization.R"))
params <- yaml::yaml.load_file(params_f)
set_promoter_window(params)
promoter_cl <- load_promoter_clusters(params, get_arg("--promoters"))
source(file.path(script_dir, "..", "utils", "txdb.R"))

# ── the "other" signal per position ───────────────────────────────────────────

ctss <- fread(ctss_f)
ctss[, position_id := paste(chr, strand, pos, sep = ":")]
ctss[, is_other := !type %in% c(END_TYPES_REAL, "longclip")]

per_type <- ctss[is_other == TRUE, .(n = sum(n)), by = .(position_id, type)]
setorder(per_type, position_id, -n, type)
other <- per_type[, .(n_other        = sum(n),
                      n_other_types  = .N,
                      top_other_type = type[1],
                      top_other_n    = n[1],
                      other_types    = paste0(type, ":", n, collapse = ";")),
                  by = position_id]
other <- merge(other,
               ctss[is_other == TRUE, .(n_samples_other = uniqueN(sample_id)),
                    by = position_id], by = "position_id")
other <- merge(other, ctss[, .(n_total = sum(n)), by = position_id],
               by = "position_id")
other <- other[n_other > min_reads]
message("Positions with more than ", min_reads, " 'other' alignments: ",
        nrow(other))

# The position's coordinates and its dominance call come from step 1's table.
dom <- fread(dom_f, select = c("position_id", "chr", "pos", "strand", "call"))
n_before <- nrow(other)
other <- merge(other, dom, by = "position_id")
stopifnot(nrow(other) == n_before)

# ── annotation ────────────────────────────────────────────────────────────────

txdb <- make_txdb_from_gff(gtf_f, format = "gtf")
gr <- GRanges(other$chr, IRanges(other$pos, width = 1), strand = other$strand)
anno <- as.data.frame(annotatePeak(gr, TxDb = txdb, tssRegion = TSS_REGION,
                                   level = "transcript", annoDb = NULL,
                                   verbose = FALSE))
stopifnot(nrow(anno) == nrow(other),
          all(anno$start == other$pos),
          all(as.character(anno$seqnames) == other$chr))
# The category goes through the shared function, so that it follows the
# promoter definition the params file selects; the gene and the distance below
# do not depend on that definition.
other[, category := annotate_positions_chipseeker(gr, txdb, promoter_gr = promoter_cl)]
other[, gene_id  := anno$geneId]
other[, distance_to_tss := anno$distanceToTSS]
# ChIPseeker assigns a feature without regard to strand, so an exonic position
# can lie antisense to its gene. geneStrand is coded 1 for + and 2 for -.
other[, same_strand := strand == c("+", "-")[anno$geneStrand]]
print(other[, .N, by = category][order(-N)])

# The GTF carries only systematic names, so the standard gene name and the SGD
# description are taken from org.Sc.sgd.db, which is keyed by the systematic
# name. A gene SGD has given no standard name keeps an empty gene_name, and the
# two columns are empty throughout when the package is absent or the gene
# identifiers are not yeast ORFs.
other[, c("gene_name", "gene_description") := ""]
if (requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
  sgd <- suppressMessages(AnnotationDbi::select(
    org.Sc.sgd.db::org.Sc.sgd.db, keys = unique(other$gene_id),
    keytype = "ORF", columns = c("GENENAME", "DESCRIPTION")))
  sgd <- sgd[!duplicated(sgd$ORF), ]
  i <- match(other$gene_id, sgd$ORF)
  other[, gene_name := fifelse(is.na(sgd$GENENAME[i]), "", sgd$GENENAME[i])]
  other[, gene_description := fifelse(is.na(sgd$DESCRIPTION[i]), "",
                                      sgd$DESCRIPTION[i])]
} else {
  message("org.Sc.sgd.db is not installed; gene_name is left empty")
}

res <- other[category == "Exon",
             .(chr, pos, strand, gene_id, gene_name, same_strand,
               distance_to_tss,
               n_other, above_10 = n_other > 10L, n_total,
               frac_other = round(n_other / n_total, 4),
               n_samples_other, n_other_types, top_other_type, top_other_n,
               other_types, position_call = call, gene_description)]
setorder(res, -n_other)

dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_f, sep = "\t")
message("Written: ", out_f, " (", nrow(res), " exonic positions, ",
        res[above_10 == TRUE, .N], " of them above 10)")
print(head(res[, !c("other_types", "gene_description")], 15))

sessionInfo()
