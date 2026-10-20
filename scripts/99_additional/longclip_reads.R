#!/usr/bin/env Rscript
# Select the longclip reads of one library into BAM files for the track hub,
# and list example reads by the genomic category of their first mapped base.
#
# A longclip read is one whose 5' soft clip is longer than dom_max_clip bases
# (4 for yeast), the same rule under which scripts/10_dominant_endtype/
# dom_helpers.R collapses such clips into the "longclip" fine type. Every read
# is 48 nt long, so "long" refers to the clip and not to the read. The 5' clip
# is the leading S of the CIGAR on the + strand and the trailing S on the -
# strand.
#
# Three BAMs are written, each sorted and indexed, because a UCSC hub "bam"
# track needs nothing more than the file and its .bai next to it:
#   <stem>_longclip.bam     both strands;
#   <stem>_longclip_F.bam   + strand only, for the hub's F composite;
#   <stem>_longclip_R.bam   - strand only, for the hub's R composite.
#
# The first mapped base is POS on the + strand and the alignment end on the -
# strand. It is annotated with the ChIPseeker call, promoter definition and
# priority order of scripts/utils/localization.R, so a read listed as "Exon"
# here is one Stages 10 and 19 would also call "Exon". The examples are drawn
# as evenly as possible across --categories. Within a category, 5' positions
# are ranked by the number of longclip reads that start there, so the most
# visible ones come first, and a position is skipped when it lies within
# --min_gap bp of one already chosen, so that the examples spread over
# different loci. One read represents each chosen position.
#
# Usage (from project root):
#   Rscript scripts/99_additional/longclip_reads.R \
#       [--bam        results/bam_subsampled/yeast_cage_ana_r2_sub.bam] \
#       [--outdir     results/tracks]                                   \
#       [--n          20]                                               \
#       [--categories "Promoter,Intron,Exon"]                           \
#       [--min_gap    500]                                              \
#       [--gtf        data/sacCer3.ensGene.gtf]                         \
#       [--enhancers  results/localization/enhancer_regions.bed]        \
#       [--params     config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(rtracklayer)
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
bam_f      <- get_arg("--bam",        "results/bam_subsampled/yeast_cage_ana_r2_sub.bam")
outdir     <- get_arg("--outdir",     "results/tracks")
n_examples <- as.integer(get_arg("--n", "20"))
categories <- strsplit(get_arg("--categories", "Promoter,Intron,Exon"), ",")[[1]]
min_gap    <- as.integer(get_arg("--min_gap", "500"))
gtf_f      <- get_arg("--gtf",        "data/sacCer3.ensGene.gtf")
enh_f      <- get_arg("--enhancers",  "results/localization/enhancer_regions.bed")
params_f   <- get_arg("--params",     "config/params.yaml")
force      <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "localization.R"))
source(file.path(script_dir, "..", "utils", "txdb.R"))

params <- yaml.load_file(params_f)
max_clip <- as.integer(params$dom_max_clip)
set_promoter_window(params)
promoter_cl <- load_promoter_clusters(params, get_arg("--promoters"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stem <- sub("_sub$", "", sub("\\.bam$", "", basename(bam_f)))
out_all <- file.path(outdir, paste0(stem, "_longclip.bam"))
out_f   <- file.path(outdir, paste0(stem, "_longclip_F.bam"))
out_r   <- file.path(outdir, paste0(stem, "_longclip_R.bam"))
tsv_f   <- file.path(outdir, paste0(stem, "_longclip_examples.tsv"))

# Length of the 5' soft clip, given the CIGAR and the reverse-strand flag.
five_prime_clip <- function(cigar, is_minus) {
  out <- integer(length(cigar))
  p <- !is_minus & grepl("^[0-9]+S", cigar)
  m <-  is_minus & grepl("[0-9]+S$", cigar)
  out[p] <- as.integer(sub("^([0-9]+)S.*$", "\\1", cigar[p]))
  out[m] <- as.integer(sub("^.*[^0-9]([0-9]+)S$", "\\1", cigar[m]))
  out
}

# ── BAMs ──────────────────────────────────────────────────────────────────────

if (!force && all(file.exists(c(out_all, out_f, out_r)))) {
  message("Longclip BAMs exist; reusing them (pass --force to rewrite)")
} else {
  message("Selecting reads with a 5' clip longer than ", max_clip, " from ", bam_f)
  rule <- FilterRules(list(longclip = function(x) {
    five_prime_clip(x$cigar, bamFlagTest(x$flag, "isMinusStrand")) > max_clip
  }))
  what <- c("cigar", "flag")
  # filterBam indexes each output itself unless indexDestination = FALSE.
  filterBam(bam_f, out_all, filter = rule, param = ScanBamParam(what = what),
            indexDestination = TRUE)
  filterBam(out_all, out_f, param = ScanBamParam(
    flag = scanBamFlag(isMinusStrand = FALSE)), indexDestination = TRUE)
  filterBam(out_all, out_r, param = ScanBamParam(
    flag = scanBamFlag(isMinusStrand = TRUE)), indexDestination = TRUE)
}

aln <- readGAlignments(out_all, use.names = TRUE)
message("Longclip alignments: ", length(aln), " (F ", sum(strand(aln) == "+"),
        ", R ", sum(strand(aln) == "-"), ")")

# ── examples ──────────────────────────────────────────────────────────────────

is_minus <- as.character(strand(aln)) == "-"
cig <- cigar(aln)
reads <- data.table(
  read_name = names(aln),
  chr       = as.character(seqnames(aln)),
  strand    = as.character(strand(aln)),
  aln_start = start(aln),
  aln_end   = end(aln),
  cigar     = cig,
  clip_len  = five_prime_clip(cig, is_minus)
)
reads[, five_prime := ifelse(strand == "+", aln_start, aln_end)]
reads[, n_at_pos := .N, by = .(chr, strand, five_prime)]

message("Loading gene annotation: ", gtf_f)
txdb <- make_txdb_from_gff(gtf_f, format = "gtf")
enhancer_gr <- if (file.exists(enh_f)) import(enh_f, format = "BED") else GRanges()

pos <- unique(reads[, .(chr, strand, five_prime)])
pos_gr <- GRanges(pos$chr, IRanges(pos$five_prime, width = 1), strand = pos$strand)
pos[, category := annotate_positions_chipseeker(pos_gr, txdb, enhancer_gr,
                                                promoter_gr = promoter_cl)]
reads <- merge(reads, pos, by = c("chr", "strand", "five_prime"))
message("Reads per category of the first mapped base:")
print(reads[, .N, by = category][order(-N)])

# Even allocation across categories, with any shortfall of a sparse category
# passed on to the others.
avail <- sapply(categories, function(k) uniqueN(reads[category == k],
                                                by = c("chr", "strand", "five_prime")))
quota <- setNames(integer(length(categories)), categories)
while (sum(quota) < n_examples && any(quota < avail)) {
  open <- names(quota)[quota < avail]
  k <- open[which.min(quota[open])]
  quota[k] <- quota[k] + 1L
}

pick <- function(k, q) {
  cand <- reads[category == k]
  setorder(cand, -n_at_pos, chr, five_prime, read_name)
  cand <- cand[!duplicated(cand[, .(chr, strand, five_prime)])]
  keep <- integer(0)
  for (i in seq_len(nrow(cand))) {
    if (length(keep) == q) break
    near <- cand$chr[keep] == cand$chr[i] &
            abs(cand$five_prime[keep] - cand$five_prime[i]) < min_gap
    if (!any(near)) keep <- c(keep, i)
  }
  cand[keep]
}
ex <- rbindlist(lapply(names(quota), function(k) pick(k, quota[[k]])))
if (nrow(ex) < n_examples)
  warning("Only ", nrow(ex), " examples found across ", paste(categories, collapse = ", "))

seq_of <- readGAlignments(out_all, use.names = TRUE,
                          param = ScanBamParam(what = "seq"))
clip_seq <- mcols(seq_of)$seq[match(ex$read_name, names(seq_of))]
clip_seq <- ifelse(ex$strand == "+",
                   as.character(subseq(clip_seq, 1, ex$clip_len)),
                   as.character(reverseComplement(
                     subseq(clip_seq, width(clip_seq) - ex$clip_len + 1))))

ex[, `:=`(clip_seq = clip_seq,
          browser  = sprintf("%s:%d-%d", chr, aln_start, aln_end))]
ex[, category := factor(category, levels = categories)]
setorder(ex, category, -n_at_pos)
out <- ex[, .(category, chr, strand, five_prime, aln_start, aln_end, browser,
              clip_len, clip_seq, n_longclip_at_pos = n_at_pos, cigar, read_name)]
fwrite(out, tsv_f, sep = "\t")
message("Wrote ", nrow(out), " examples to ", tsv_f)
print(out[, .(category, browser, strand, five_prime, clip_len, clip_seq,
              n_longclip_at_pos)])

sessionInfo()
