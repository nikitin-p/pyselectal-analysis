#!/usr/bin/env Rscript
# Stage 17 step 4: genomic context of dominated clusters.
#
# This is the control that decides whether a dominated cluster is interesting.
# A cluster dominated by, say, 2Sga whose genome reads "GA" immediately
# upstream of the start site is a mapping artefact: the clipped bases are
# templated and the aligner merely chose to clip them. A cluster whose upstream
# genome differs carries a genuine untemplated 5' extension.
#
# The script also reports the frequency of a genomic G at the first templated
# base per class, which is the masking confound Stage 15 identified for M and
# which must be reported now that M carries equal standing among the classes.
#
# Usage (from project root):
#   Rscript scripts/10_dominant_endtype/04_context.R \
#       [--calls   results/dominance/dominance_calls.tsv] \
#       [--matrix  results/tss/tss_matrix.tsv]            \
#       [--samples config/samples.tsv]                    \
#       [--outdir  results/dominance] [--figdir results/figures] \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
calls_f   <- get_arg("--calls",   "results/dominance/dominance_calls.tsv")
matrix_f  <- get_arg("--matrix",  "results/tss/tss_matrix.tsv")
samples_f <- get_arg("--samples", "config/samples.tsv")
outdir    <- get_arg("--outdir",  "results/dominance")
figdir    <- get_arg("--figdir",  "results/figures")
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(outdir, "genomic_context.tsv")
if (!force && file.exists(out_tsv)) {
  message("Output exists, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}

# Genome loading and seqlevel harmonisation, carried over verbatim from
# scripts/05_tss_end_types/02_cluster_annotate.R lines 62-102. Setting
# seqlengths from the FaFile is required for trim() to clip out-of-bound ranges.
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec)) stop("FASTA not found: ", genome_spec)
    if (!file.exists(paste0(genome_spec, ".fai"))) indexFa(genome_spec)
    fa <- FaFile(genome_spec); open(fa); return(fa)
  }
  pkg <- switch(genome_spec,
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias: ", genome_spec))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("BSgenome not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

harmonize <- function(gr, genome) {
  fa_si <- seqinfo(genome)
  keep  <- intersect(seqlevels(gr), as.character(seqnames(fa_si)))
  gr    <- keepSeqlevels(gr, keep, pruning.mode = "coarse")
  seqlengths(gr) <- seqlengths(fa_si)[seqlevels(gr)]
  gr
}

meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"),
                 show_col_types = FALSE)
genome <- load_genome(unique(na.omit(meta$genome))[1])
on.exit(if (inherits(genome, "FaFile")) close(genome), add = TRUE)

calls <- fread(calls_f)[dominated == TRUE]
coord <- unique(fread(matrix_f,
                      select = c("cluster_id", "chr", "tss_pos", "strand")))
dt <- merge(calls, coord, by = "cluster_id")
if (nrow(dt) == 0) stop("No dominated clusters to annotate.")
message("Stage 17 step 4: context for ", nrow(dt), " dominated cluster calls")

gr <- GRanges(dt$chr, IRanges(dt$tss_pos, width = 1), strand = dt$strand)
gr <- harmonize(gr, genome)
# harmonize() prunes ranges on seqlevels the FASTA does not carry, after which
# every dt[i, ...] assignment below would index dt by a position in the
# SHORTENED gr and silently attach a sequence to the wrong cluster. That never
# happens on sacCer3, where the FASTA carries every contig, so the mismatch is
# an error rather than a case to handle.
stopifnot(length(gr) == nrow(dt))

# Five bases immediately upstream of the TSS, strand-aware and 5'->3'.
up <- trim(promoters(gr, upstream = 5, downstream = 0))
dt[, upstream5 := NA_character_]
ok <- which(width(up) == 5)
if (length(ok) > 0) dt[ok, upstream5 := as.character(getSeq(genome, up[ok]))]

# The initiator dinucleotide [TSS-1, TSS], as used everywhere else in the
# pipeline, and the first templated base, which is the TSS base itself.
ini <- trim(promoters(gr, upstream = 1, downstream = 1))
dt[, initiator := NA_character_]
ok2 <- which(width(ini) == 2)
if (length(ok2) > 0) dt[ok2, initiator := as.character(getSeq(genome, ini[ok2]))]
dt[, first_templated := substr(initiator, 2, 2)]

# Templated test: does the genome immediately upstream reproduce the clip?
# For a dom_type of "<len>S<seq>", compare <seq> against the last <len> bases of
# upstream5. M, longclip and Nclip have no clip to compare and are NA.
dt[, clip_seq := toupper(sub("^\\d+S", "", dom_type))]
dt[dom_type %in% c("M", "longclip", "Nclip"), clip_seq := NA_character_]
dt[, templated := NA]
dt[!is.na(clip_seq) & !is.na(upstream5),
   templated := substr(upstream5, 6 - nchar(clip_seq), 5) == clip_seq]

frac <- dt[, .(n_clusters = .N,
               n_templated = sum(templated, na.rm = TRUE),
               frac_templated = mean(templated, na.rm = TRUE),
               frac_first_templated_G = mean(first_templated == "G",
                                             na.rm = TRUE)),
           by = .(threshold, class = dom_type)][order(threshold, -n_clusters)]
fwrite(frac, file.path(outdir, "templated_fraction.tsv"), sep = "\t")
print(frac)

# ── figures ───────────────────────────────────────────────────────────────────
pf <- melt(frac, id.vars = c("threshold", "class", "n_clusters"),
           measure.vars = c("frac_templated", "frac_first_templated_G"),
           variable.name = "measure", value.name = "fraction")
p1 <- ggplot(pf[!is.na(fraction)],
             aes(reorder(class, -n_clusters), fraction, fill = measure)) +
  geom_col(position = "dodge") +
  facet_wrap(~ threshold, ncol = 1) +
  scale_fill_manual(values = c(frac_templated = "#457b9d",
                               frac_first_templated_G = "#e63946"),
                    labels = c("clip matches upstream genome",
                               "genomic G at first templated base")) +
  labs(x = "dominance class", y = "fraction of clusters", fill = NULL,
       title = "Templated artefact and G-masking controls by dominance class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
ggsave(file.path(figdir, "18_dom_context_templated.pdf"), p1,
       width = 9, height = 8)

# Initiator composition, with canonical YR initiators in the project red.
CANON <- c("CA", "CG", "TA", "TG")
ic <- dt[!is.na(initiator), .N, by = .(threshold, class = dom_type, initiator)]
ic[, frac := N / sum(N), by = .(threshold, class)]
p2 <- ggplot(ic, aes(initiator, frac, fill = initiator %in% CANON)) +
  geom_col(show.legend = FALSE) +
  facet_grid(class ~ threshold, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#e63946", `FALSE` = "grey60")) +
  labs(x = "initiator dinucleotide", y = "fraction of clusters",
       title = "Initiator composition by dominance class (YR in red)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(figdir, "18_dom_context_initiator.pdf"), p2,
       width = 9, height = 11)

# The idempotency sentinel is out_tsv, so it is written LAST: an interrupted run
# or a failing ggsave must not leave a sentinel behind that makes the next run
# skip this step while its figures are missing.
fwrite(dt[, .(cluster_id, threshold, dom_type, chr, tss_pos, strand,
              upstream5, initiator, first_templated, templated)],
       out_tsv, sep = "\t")
message("Written: ", out_tsv)
sessionInfo()
