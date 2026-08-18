#!/usr/bin/env Rscript
# Stage 16 step 1: pool per-type read counts across samples, assign each TSS
# cluster to the gene it plausibly drives, and define the background universe.
#
# The assignment rule is strand-aware and windowed. A TSS is assigned to a gene
# when it lies at most --upstream bp before that gene's annotated 5' end on the
# same strand, or at most --downstream bp inside the gene from that 5' end. A
# TSS matching several genes goes to the nearest; a TSS matching none is
# dropped. Antisense TSSs are dropped, because a TSS on the opposite strand
# does not drive the gene.
#
# Percentages are pooled: each type's counts are summed over all samples and
# divided once at the end. A mean of per-sample percentages was rejected because
# it lets a sample contributing two reads outvote one contributing two thousand.
#
# Input : results/tss/tss_matrix.tsv, results/tss/tss_clustered.tsv,
#         data/sacCer3.genes.bed
# Output: results/go/tss_gene_assignment.tsv, results/go/universe.txt

suppressMessages({
    library(optparse)
    library(GenomicRanges)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--tss_matrix",    type = "character"),
    make_option("--tss_clustered", type = "character"),
    make_option("--genes_bed",     type = "character"),
    make_option("--outdir",        type = "character"),
    make_option("--upstream",      type = "integer", default = 500L),
    make_option("--downstream",    type = "integer", default = 100L),
    make_option("--min_total",     type = "double",  default = 2),
    make_option("--force",         action = "store_true", default = FALSE)
)))

TYPES <- END_TYPES

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
assign_out <- file.path(opt$outdir, "tss_gene_assignment.tsv")
univ_out   <- file.path(opt$outdir, "universe.txt")

if (file.exists(assign_out) && file.exists(univ_out) && !opt$force) {
    message("[skip] outputs exist; use --force to regenerate")
    quit(save = "no")
}

# ---- pooled per-type counts per TSS cluster -------------------------------
mat <- read.delim(opt$tss_matrix, stringsAsFactors = FALSE)
message(sprintf("Loaded %d rows from the TSS matrix", nrow(mat)))

n_cols <- paste0("n_", TYPES)
missing <- setdiff(n_cols, colnames(mat))
if (length(missing)) {
    stop("TSS matrix is missing count columns: ", paste(missing, collapse = ", "))
}

pooled <- aggregate(mat[, n_cols], by = list(cluster_id = mat$cluster_id), FUN = sum)
pooled$n_total <- rowSums(pooled[, n_cols])
for (ty in TYPES) {
    pooled[[paste0("pct_", ty)]] <-
        ifelse(pooled$n_total > 0, 100 * pooled[[paste0("n_", ty)]] / pooled$n_total, 0)
}
message(sprintf("Pooled to %d TSS clusters", nrow(pooled)))

# ---- TSS coordinates, strand, initiator and expression --------------------
clu <- read.delim(opt$tss_clustered, stringsAsFactors = FALSE)
keep <- c("cluster_id", "chr", "tss_pos", "strand", "initiator", "mean_total")
clu <- clu[, keep]

clu <- clu[clu$mean_total >= opt$min_total, ]
message(sprintf("After mean_total >= %s: %d TSS clusters", opt$min_total, nrow(clu)))

tss <- merge(clu, pooled, by = "cluster_id")
stopifnot(nrow(tss) > 0)

# ---- gene models ----------------------------------------------------------
genes <- read.delim(opt$genes_bed, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chr", "start", "end", "gene", "score", "strand"))
message(sprintf("Loaded %d genes", nrow(genes)))

shared <- intersect(unique(genes$chr), unique(tss$chr))
if (!length(shared)) {
    stop("No chromosome names are shared between the gene BED (",
         paste(head(unique(genes$chr)), collapse = ", "),
         ") and the TSS table (",
         paste(head(unique(tss$chr)), collapse = ", "),
         "). This is a reference mismatch, not a filtering artefact.")
}

# A gene's 5' end, 1-based: the low coordinate on +, the high coordinate on -.
gene_5p <- ifelse(genes$strand == "+", genes$start + 1L, genes$end)

# Acceptance zone around that 5' end, oriented by strand.
zone_start <- ifelse(genes$strand == "+",
                     gene_5p - opt$upstream, gene_5p - opt$downstream)
zone_end   <- ifelse(genes$strand == "+",
                     gene_5p + opt$downstream, gene_5p + opt$upstream)
zone_start <- pmax(zone_start, 1L)

zone_gr <- GRanges(seqnames = genes$chr,
                   ranges   = IRanges(start = zone_start, end = zone_end),
                   strand   = genes$strand)
mcols(zone_gr)$gene    <- genes$gene
mcols(zone_gr)$gene_5p <- gene_5p

tss_gr <- GRanges(seqnames = tss$chr,
                  ranges   = IRanges(start = tss$tss_pos, width = 1),
                  strand   = tss$strand)

# findOverlaps respects strand by default, so antisense pairs never appear.
hits <- findOverlaps(tss_gr, zone_gr)
message(sprintf("%d TSS-gene candidate pairs within the window",
                length(hits)))
if (!length(hits)) stop("No TSS fell within any gene's acceptance zone.")

qi <- queryHits(hits)
si <- subjectHits(hits)

# Signed distance: positive when the TSS is upstream of the gene's 5' end.
tp   <- tss$tss_pos[qi]
g5   <- mcols(zone_gr)$gene_5p[si]
strd <- tss$strand[qi]
dist <- ifelse(strd == "+", g5 - tp, tp - g5)

# Nearest gene wins; ties break on the first, which is deterministic because
# the BED is sorted.
ord  <- order(qi, abs(dist))
best <- ord[!duplicated(qi[ord])]

out <- tss[qi[best], ]
out$gene               <- mcols(zone_gr)$gene[si[best]]
out$dist_to_gene_start <- dist[best]

col_order <- c("cluster_id", "chr", "tss_pos", "strand", "initiator",
               "mean_total", "gene", "dist_to_gene_start",
               paste0("n_", TYPES), paste0("pct_", TYPES), "n_total")
out <- out[, col_order]
out <- out[order(out$chr, out$tss_pos), ]

write.table(out, assign_out, sep = "\t", quote = FALSE, row.names = FALSE)

universe <- sort(unique(out$gene))
writeLines(universe, univ_out)

message(sprintf("Assigned %d of %d TSS clusters (%.1f %%) to %d genes",
                nrow(out), nrow(tss), 100 * nrow(out) / nrow(tss),
                length(universe)))
message("Distance to gene 5' end (positive = upstream):")
print(summary(out$dist_to_gene_start))

sessionInfo()
