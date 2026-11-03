#!/usr/bin/env Rscript
# Stage 22 step 1: isolate the 5' soft clip of every long-clip read, find where
# in the genome it comes from, and annotate the read's inferred 5' end.
#
# A long-clip read is one whose 5' clip is longer than dom_max_clip and carries
# no N, which is exactly the "longclip" fine type of scripts/10_dominant_endtype/
# dom_helpers.R; the read count is asserted against Stage 17's table.
#
# Up to lc_max_untemplated leading clip bases may be untemplated (a cap G, for
# instance). The annotated-intron model takes the number giving the best
# identity, then the smallest; the exact search takes the smallest that hits. Each read takes
# the first origin class that applies (lc_helpers.R):
#   splice_annotated    the aligned 5' end lies at most lc_acceptor_window bases
#                       inside an annotated intron's 3' end, or on the first base
#                       of the downstream exon, and the clip is at least
#                       lc_min_identity identical to the exon-1 sequence that ends
#                       at the donor, shifted by the intron bases the aligned part
#                       covers; or route B below finds a GT..AG gap that equals an
#                       annotated intron;
#   splice_unannotated  route B: the clip occurs exactly within lc_search_window
#                       bases upstream on the read strand, the gap is at least
#                       lc_min_intron, and the gap, slid by up to
#                       lc_acceptor_window bases, reads GT...AG;
#   templated_adjacent  route B with a gap of at most lc_adjacent_gap bases;
#   templated_distal    route B with any other gap;
#   unplaced            none of the above ("short" when the searchable clip is
#                       below lc_min_search, "not_found" otherwise).
# The same procedure is run on the reversed clip as a chance-match control,
# whose class is recorded beside the real one.
#
# The inferred 5' end is the genomic position of the first templated clip base,
# or the aligned 5' end for an unplaced read. It is placed in a TSS class
# (canonical / alternative_promoter / other_cluster / no_cluster) against the
# Stage 8 distclu clusters and the promoter clusters of
# scripts/05_tss_end_types/12_promoter_clusters.R, and the 1Sg reads within
# +-1 bp of it and the Stage 19 call at it are attached.
#
# Usage (from project root):
#   Rscript scripts/14_longclip/01_classify.R \
#       [--samples   config/samples.tsv]                            \
#       [--fasta     data/sacCer3.fa] [--gtf data/sacCer3.ensGene.gtf] \
#       [--typed     results/tss/typed_ctss.tsv.gz]                 \
#       [--matrix    results/tss/tss_matrix.tsv]                    \
#       [--promoters results/tss/promoter_clusters.tsv]             \
#       [--fine      results/dominance/fine_typed_ctss.tsv.gz]      \
#       [--dominance results/exotic_origins/position_dominance.tsv.gz] \
#       [--outdir    results/longclip] [--params config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(Rsamtools)
  library(data.table)
  library(stringi)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
samples_f <- get_arg("--samples",   "config/samples.tsv")
fasta_f   <- get_arg("--fasta",     "data/sacCer3.fa")
gtf_f     <- get_arg("--gtf",       "data/sacCer3.ensGene.gtf")
typed_f   <- get_arg("--typed",     "results/tss/typed_ctss.tsv.gz")
mat_f     <- get_arg("--matrix",    "results/tss/tss_matrix.tsv")
prom_f    <- get_arg("--promoters", "results/tss/promoter_clusters.tsv")
fine_f    <- get_arg("--fine",      "results/dominance/fine_typed_ctss.tsv.gz")
pdom_f    <- get_arg("--dominance", "results/exotic_origins/position_dominance.tsv.gz")
outdir    <- get_arg("--outdir",    "results/longclip")
params_f  <- get_arg("--params",    "config/params.yaml")
force     <- "--force" %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))
source(file.path(script_dir, "lc_helpers.R"))

params    <- yaml.load_file(params_f)
max_clip  <- as.integer(params$dom_max_clip)
U         <- as.integer(params$lc_max_untemplated)
acc_win   <- as.integer(params$lc_acceptor_window)
min_srch  <- as.integer(params$lc_min_search)
W         <- as.integer(params$lc_search_window)
min_intr  <- as.integer(params$lc_min_intron)
adj_gap   <- as.integer(params$lc_adjacent_gap)
min_id    <- as.numeric(params$lc_min_identity)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out_gz <- file.path(outdir, "longclip_reads.tsv.gz")
if (!force && file.exists(out_gz)) {
  message("Output exists, skipping (use --force to rerun): ", out_gz)
  quit(save = "no", status = 0)
}

# ── 1. the long-clip reads ────────────────────────────────────────────────────

read_longclip <- function(bam_path, sample_id) {
  message("  reading: ", sample_id)
  param <- ScanBamParam(
    flag = scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                       isSupplementaryAlignment = FALSE),
    what = "seq")
  ga <- readGAlignments(bam_path, param = param)
  is_plus <- as.character(strand(ga)) %in% c("+", "*")
  clen <- clip5_len(cigar(ga), is_plus)
  keep <- clen > max_clip
  ga <- ga[keep]; clen <- clen[keep]; is_plus <- is_plus[keep]
  clip <- clip5_seq(as.character(mcols(ga)$seq), clen, is_plus)
  dt <- data.table(chr = as.character(seqnames(ga)),
                   strand = ifelse(is_plus, "+", "-"),
                   aln_start = start(ga), aln_end = end(ga),
                   cigar = cigar(ga), clen = clen, clip = clip)
  dt <- dt[!grepl("N", clip, fixed = TRUE)]
  dt[, .(n = .N), by = .(chr, strand, aln_start, aln_end, cigar, clen, clip)][
    , sample_id := sample_id][]
}

meta <- readr::read_tsv(samples_f, comment = "#", col_types = readr::cols(.default = "c"),
                        show_col_types = FALSE)
stopifnot("sample_id" %in% colnames(meta), nrow(meta) > 0)
if (!"bam" %in% colnames(meta)) meta$bam <- NA_character_
bams <- ifelse(!is.na(meta$bam) & meta$bam != "", meta$bam,
               file.path(params$results_bam_sub, paste0(meta$sample_id, ".bam")))
missing <- !file.exists(bams)
if (any(missing)) warning("BAM not found, skipping: ", paste(bams[missing], collapse = ", "))
message("Stage 22 step 1: long-clip origins (clip > ", max_clip, " nt)")
per_sample <- rbindlist(Map(read_longclip, bams[!missing], meta$sample_id[!missing]))

reads <- per_sample[, .(n_reads = sum(n), n_samples = uniqueN(sample_id)),
                    by = .(chr, strand, aln_start, aln_end, cigar, clen, clip)]
reads[, pos5 := ifelse(strand == "+", aln_start, aln_end)]
reads[, read_key := .I]
message("Long-clip reads: ", sum(reads$n_reads), " in ", nrow(reads),
        " distinct alignments")

fine <- fread(fine_f)
n_stage17 <- fine[type == "longclip", sum(n)]
if (sum(reads$n_reads) != n_stage17) {
  stop("Long-clip reads here (", sum(reads$n_reads), ") differ from Stage 17's ",
       "longclip count (", n_stage17, "); the BAMs or dom_max_clip have changed.")
}

fa      <- FaFile(fasta_f)
seqlens <- seqlengths(fa)

# ── 2. the annotated introns ──────────────────────────────────────────────────

txdb  <- suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGFF(gtf_f)))
intr  <- unique(unlist(GenomicFeatures::intronsByTranscript(txdb)))
introns <- unique(data.table(chr = as.character(seqnames(intr)),
                             strand = as.character(strand(intr)),
                             i_start = start(intr), i_end = end(intr)))
# ex2: first base of the downstream exon; don: last base of the upstream exon,
# both in the transcript's own orientation.
introns[, ex2 := ifelse(strand == "+", i_end + 1L, i_start - 1L)]
introns[, don := ifelse(strand == "+", i_start - 1L, i_end + 1L)]
introns[, intron_id := paste0(chr, ":", i_start, "-", i_end, ":", strand)]
message("Annotated introns: ", nrow(introns))

# ── 3. route A: the annotated acceptor ────────────────────────────────────────

route_a <- function(reads, query_fun) {
  cand <- introns[reads[, .(read_key, chr, strand, pos5, clip, clen)],
                  on = .(chr, strand), allow.cartesian = TRUE, nomatch = 0]
  cand[, off := ifelse(strand == "+", pos5 - ex2, ex2 - pos5)]
  cand <- cand[off >= -acc_win & off <= 0]
  if (!nrow(cand)) return(NULL)
  cand[, k := -off]
  hits <- rbindlist(lapply(0:U, function(u) {
    x <- copy(cand)[clen - u >= 1]
    x[, u := u]
    x[, q := query_fun(substr(clip, u + 1L, clen))]
    x[, nq := nchar(q)]
    # The template ends k bases before the donor, read-strand.
    x[, t3 := ifelse(strand == "+", don - k, don + k)]
    x[, t5 := ifelse(strand == "+", t3 - nq + 1L, t3 + nq - 1L)]
    x[, tmpl := strand_seq(fa, chr, pmin(t5, t3), pmax(t5, t3), strand, seqlens)]
    x[, ident := identity_frac(q, ifelse(is.na(tmpl), strrep("N", nq), tmpl))]
    x[ident >= min_id]
  }))
  if (!nrow(hits)) return(NULL)
  # Best identity first, then the fewest untemplated bases: a cap G that happens
  # to mismatch the template must not be drawn onto the genome just because the
  # remaining bases clear the identity floor with it included.
  setorder(hits, read_key, -ident, u)
  hits <- hits[, .SD[1], by = read_key]
  hits[, .(read_key, class = "splice_annotated", n_untemplated = u,
           tss5 = t5, tmpl_start = pmin(t5, t3), tmpl_end = pmax(t5, t3),
           gap = NA_integer_, donor = don, intron_id, acceptor_offset = off,
           identity = ident)]
}

# ── 4. route B: the local upstream search ─────────────────────────────────────

route_b <- function(reads, query_fun) {
  x <- copy(reads[, .(read_key, chr, strand, pos5, clip, clen)])
  x[, wstart := ifelse(strand == "+", pmax(1L, pos5 - W), pos5 + 1L)]
  x[, wend   := ifelse(strand == "+", pos5 - 1L, pmin(seqlens[chr], pos5 + W))]
  x[, win := strand_seq(fa, chr, wstart, wend, strand, seqlens)]
  x[, `:=`(u = NA_integer_, loc_s = NA_integer_, loc_e = NA_integer_)]
  # The loop variable is not named u, because inside x[...] that name would
  # resolve to the column rather than to the loop.
  for (n_u in 0:U) {
    todo <- which(is.na(x$u) & x$clen - n_u >= min_srch & !is.na(x$win))
    if (!length(todo)) next
    q   <- query_fun(substr(x$clip[todo], n_u + 1L, x$clen[todo]))
    loc <- stri_locate_last_fixed(x$win[todo], q)
    hit <- !is.na(loc[, "end"])
    set(x, todo[hit], "u",     n_u)
    set(x, todo[hit], "loc_s", as.integer(loc[hit, "start"]))
    set(x, todo[hit], "loc_e", as.integer(loc[hit, "end"]))
  }
  found <- x[!is.na(u)]
  found[, gap := nchar(win) - loc_e]
  found[, tss5 := win_to_genome(loc_s, strand, wstart, wend)]
  found[, t3   := win_to_genome(loc_e, strand, wstart, wend)]
  # The splice-motif test, on the gap slid by up to acc_win bases toward the read.
  found[, ext := NA_character_]
  sp <- found[gap >= min_intr, which = TRUE]
  if (length(sp)) {
    f <- found[sp]
    es <- ifelse(f$strand == "+", f$t3 + 1L, f$pos5 + 1L - acc_win)
    ee <- ifelse(f$strand == "+", f$pos5 - 1L + acc_win, f$t3 - 1L)
    found[sp, ext := strand_seq(fa, chr, es, ee, strand, seqlens)]
  }
  found[, shift := splice_shift(ext, gap, acc_win)]
  found[, `:=`(
    i_start = ifelse(strand == "+", t3 + 1L + shift, pos5 + 1L - shift),
    i_end   = ifelse(strand == "+", pos5 - 1L + shift, t3 - 1L - shift),
    donor   = ifelse(strand == "+", t3 + shift, t3 - shift))]
  found[, intron_id := ifelse(is.na(shift), NA_character_,
                              paste0(chr, ":", i_start, "-", i_end, ":", strand))]
  found[, class := fcase(
    !is.na(shift) & intron_id %in% introns$intron_id, "splice_annotated",
    !is.na(shift),                                    "splice_unannotated",
    gap <= adj_gap,                                   "templated_adjacent",
    default =                                         "templated_distal")]
  found[, acceptor_offset := NA_integer_]
  hit_a <- found[class == "splice_annotated"]
  if (nrow(hit_a)) {
    ex2 <- introns[match(hit_a$intron_id, introns$intron_id), ex2]
    found[class == "splice_annotated",
          acceptor_offset := ifelse(strand == "+", pos5 - ex2, ex2 - pos5)]
  }
  found[, .(read_key, class, n_untemplated = u, tss5,
            tmpl_start = pmin(tss5, t3), tmpl_end = pmax(tss5, t3),
            gap, donor = ifelse(is.na(shift), NA_integer_, donor),
            intron_id, acceptor_offset, identity = 1)]
}

classify <- function(reads, query_fun) {
  a <- route_a(reads, query_fun)
  rest <- if (is.null(a)) reads else reads[!read_key %in% a$read_key]
  placed <- rbind(a, route_b(rest, query_fun), fill = TRUE)
  un <- reads[!read_key %in% placed$read_key, .(read_key, pos5, clen)]
  un <- un[, .(read_key, class = "unplaced", n_untemplated = NA_integer_,
               tss5 = pos5, unplaced_reason = ifelse(clen >= min_srch, "not_found", "short"))]
  out <- rbind(placed, un, fill = TRUE)
  stopifnot(nrow(out) == nrow(reads), !anyDuplicated(out$read_key))
  out
}

message("Classifying the clips")
real <- classify(reads, identity)
message("Classifying the reversed-clip control")
ctrl <- classify(reads, stri_reverse)[, .(read_key, control_class = class)]

reads <- merge(reads, real, by = "read_key")
reads <- merge(reads, ctrl, by = "read_key")
reads[, class := factor(class, levels = LC_CLASSES)]

# ── 5. the TSS context of the inferred 5' end ─────────────────────────────────

typed <- fread(typed_f, select = c("chr", "pos", "strand", "n_total"))
pooled <- typed[, .(n_total = sum(n_total)), by = .(chr, pos, strand)]
members <- derive_clusters(pooled[n_total >= as.integer(params$a8_min_total)],
                           as.integer(params$a8_max_dist))
stage8 <- unique(fread(mat_f, select = "cluster_id"))$cluster_id
if (!setequal(unique(members$cluster_id), stage8)) {
  stop("Re-derived cluster set differs from ", mat_f,
       "; check that a8_min_total and a8_max_dist equal Stage 8's settings.")
}
extent <- members[, .(c_start = min(pos), c_end = max(pos)),
                  by = .(cluster_id, chr, strand)]

prom <- fread(prom_f)
canon <- prom[prom[, .I[which.max(reads)], by = gene]$V1, unique(cluster_id)]
extent[, tss_class := fcase(cluster_id %in% canon,          "canonical",
                            cluster_id %in% prom$cluster_id, "alternative_promoter",
                            default =                        "other_cluster")]

q5 <- reads[, .(read_key, chr, strand, s = tss5, e = tss5)]
setkey(extent, chr, strand, c_start, c_end)
ov <- foverlaps(q5, extent, by.x = c("chr", "strand", "s", "e"),
                by.y = c("chr", "strand", "c_start", "c_end"), nomatch = 0)
reads <- merge(reads, ov[, .(read_key, cluster_id, tss_class)], by = "read_key",
               all.x = TRUE)
reads[is.na(tss_class), tss_class := "no_cluster"]

# ── 6. co-production with 1Sg ─────────────────────────────────────────────────

sg <- fine[type == END_TYPE_CANONICAL, .(n_sg = sum(n)), by = .(chr, strand, pos)]
nb <- reads[, .(read_key, chr, strand, pos = tss5)][
  , .(pos = pos + (-1:1)), by = .(read_key, chr, strand)]
nb <- merge(nb, sg, by = c("chr", "strand", "pos"), all.x = TRUE)
pdom <- fread(pdom_f, select = c("chr", "pos", "strand", "call"))
nb <- merge(nb, pdom, by = c("chr", "strand", "pos"), all.x = TRUE)
co <- nb[, .(sg_reads_pm1 = sum(n_sg, na.rm = TRUE),
             sg_dominated_pm1 = any(call == END_TYPE_CANONICAL, na.rm = TRUE)),
         by = read_key]
reads <- merge(reads, co, by = "read_key")
reads <- merge(reads, pdom[, .(chr, strand, tss5 = pos, call_at_tss5 = call)],
               by = c("chr", "strand", "tss5"), all.x = TRUE)

setorder(reads, chr, strand, pos5)
setcolorder(reads, c("read_key", "chr", "strand", "pos5", "aln_start", "aln_end",
                     "cigar", "clen", "clip", "n_reads", "n_samples", "class",
                     "control_class", "n_untemplated", "tss5", "tmpl_start",
                     "tmpl_end", "gap", "donor", "intron_id", "acceptor_offset",
                     "identity", "unplaced_reason", "cluster_id", "tss_class",
                     "sg_reads_pm1", "sg_dominated_pm1", "call_at_tss5"))

summ <- reads[, .(n_reads = sum(n_reads), n_alignments = .N), by = class][order(class)]
summ[, pct_reads := round(100 * n_reads / sum(n_reads), 2)]
print(summ)

fwrite(reads, out_gz, sep = "\t", compress = "gzip")
message("Written: ", out_gz)
sessionInfo()
