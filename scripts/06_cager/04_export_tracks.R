#!/usr/bin/env Rscript
# Stage 10 step 4: browser tracks from the CAGEr object of step 1, which is built
# from the subsampled libraries.
#
#   <sample>_tpm_{plus,minus}.bw       power-law normalized TPM per CTSS, restricted
#                                      to the CTSSs that pass filterLowExpCTSS
#                                      (cager_tpm_threshold TPM in at least
#                                      cager_tpm_pass_samples samples); the filter
#                                      is across samples, as in CAGEr, so a kept
#                                      CTSS shows its value in every library even
#                                      where that value is below the threshold
#   <sample>_tagClusters.bed           the sample's own distclu tag clusters, BED12:
#                                      the dominant CTSS as the thick part and the
#                                      interquantile range as the block
#   consensusClusters.bed              consensus clusters across samples, BED12
#   *.bb                               the BEDs as bigBed, when --bedtobigbed is given
#   track_manifest.tsv
#
# Minus-strand values are written positive, as in Stage 21. The BEDs are written
# with rtracklayer::export.bed, which encodes BED12 blocks correctly (its
# export.bb does not); export.bw was verified against a UCSC-written bigWig.
#
# Usage (from project root):
#   Rscript scripts/06_cager/04_export_tracks.R \
#       [--rds results/cager/cageexp.rds] [--fai data/sacCer3.fa.fai] \
#       [--outdir results/cager/tracks] [--params config/params.yaml] \
#       [--bedtobigbed PATH] [--force]

suppressPackageStartupMessages({
  library(CAGEr)
  library(data.table)
  library(rtracklayer)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
rds_f    <- get_arg("--rds",         "results/cager/cageexp.rds")
fai_f    <- get_arg("--fai",         "data/sacCer3.fa.fai")
outdir   <- get_arg("--outdir",      "results/cager/tracks")
params_f <- get_arg("--params",      "config/params.yaml")
b2bb     <- get_arg("--bedtobigbed", "")
force    <- "--force" %in% args

params <- yaml.load_file(params_f)
iq_low <- as.numeric(params$cager_iqw_qlow %||% 0.1)
iq_up  <- as.numeric(params$cager_iqw_qup  %||% 0.9)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_f <- file.path(outdir, "track_manifest.tsv")
if (!force && file.exists(manifest_f)) {
  message("Output exists, skipping (use --force to rerun): ", manifest_f)
  quit(save = "no", status = 0)
}
if (!file.exists(rds_f)) stop("CAGEexp not found: ", rds_f, "\nRun step 1 first.")

# exportToTrack calls sampleLabels<- on a GRanges, which some CAGEr versions do
# not define; the same workaround as in 01_cageexp.R.
if (is.null(suppressWarnings(tryCatch(methods::getMethod("sampleLabels<-", "GRanges"),
                                      error = function(e) NULL)))) {
  setMethod("sampleLabels<-", "GRanges", function(object, value) {
    S4Vectors::metadata(object)$sampleLabels <- value
    object
  })
}

ce  <- readRDS(rds_f)
fai <- fread(fai_f, header = FALSE, select = 1:2, col.names = c("chr", "len"))
si  <- Seqinfo(fai$chr, fai$len)
unlink(list.files(outdir, pattern = "\\.(bw|bed|bb)$", full.names = TRUE))

# ── normalized, filtered CTSS signal ──────────────────────────────────────────

keep <- as.logical(filteredCTSSidx(ce))
message("CTSSs passing the filter: ", sum(keep), " of ", length(keep))

write_bw <- function(gr, path) {
  if (!length(gr)) return(FALSE)
  gr <- GRanges(seqnames(gr), IRanges(start(gr), width = 1), score = as.numeric(score(gr)),
                seqinfo = si)
  # Fixed zoom levels, as in Stage 21: UCSC's bbiWrite fails on data-derived
  # ones for very sparse tracks.
  export.bw(sort(gr), path, fixedSummaries = TRUE)
  TRUE
}

# The normalized matrix is read from the filtered experiment rather than through
# CTSSnormalizedTpmGR(), which returns only the sample's non-zero positions and
# so cannot be indexed by the filter.
se <- CTSStagCountSE(ce)[keep]
tpm_mat <- assay(se, "normalizedTpmMatrix")
manifest <- rbindlist(lapply(sampleLabels(ce), function(s) {
  gr <- rowRanges(se)
  score(gr) <- as.numeric(tpm_mat[[s]])
  gr <- gr[score(gr) > 0]
  rbindlist(lapply(c(plus = "+", minus = "-"), function(st) {
    g <- gr[as.character(strand(gr)) == st]
    sn <- if (st == "+") "plus" else "minus"
    file <- paste0(s, "_tpm_", sn, ".bw")
    ok <- write_bw(g, file.path(outdir, file))
    data.table(file = if (ok) file else NA_character_, kind = "tpm", sample = s,
               strand = st, n_features = length(g), total = sum(score(g)))
  }))
}))

# ── per-sample tag clusters and consensus clusters ────────────────────────────

# UCSC requires an integer BED score in 0..1000, whereas CAGEr writes the TPM
# there, so the score is rounded and capped and the exact TPM moves into the
# name, as scripts/99_additional/ucsc_tracks.sh does for the older CAGEr BEDs.
write_cluster_bed <- function(gr, path) {
  gr <- methods::as(gr, "GRanges")
  tpm <- as.numeric(score(gr))
  gr$name  <- sprintf("tpm=%.2f", tpm)
  score(gr) <- pmin(1000L, as.integer(round(tpm)))
  export.bed(gr, path)
  list(n = length(gr), total = sum(tpm))
}

tc <- exportToTrack(ce, what = "tagClusters", qLow = iq_low, qUp = iq_up, oneTrack = FALSE)
for (s in names(tc)) {
  file <- paste0(s, "_tagClusters.bed")
  w <- write_cluster_bed(tc[[s]], file.path(outdir, file))
  manifest <- rbind(manifest, data.table(file = file, kind = "tagClusters", sample = s,
                                         strand = "both", n_features = w$n, total = w$total))
}
cc <- exportToTrack(ce, what = "consensusClusters", qLow = iq_low, qUp = iq_up,
                    colorByExpressionProfile = FALSE)
w <- write_cluster_bed(cc, file.path(outdir, "consensusClusters.bed"))
manifest <- rbind(manifest, data.table(file = "consensusClusters.bed", kind = "consensusClusters",
                                       sample = "all", strand = "both",
                                       n_features = w$n, total = w$total))

if (nzchar(b2bb)) {
  sizes <- tempfile(fileext = ".sizes")
  fwrite(fai, sizes, sep = "\t", col.names = FALSE)
  for (f in manifest[grepl("\\.bed$", file), file]) {
    body <- tempfile(fileext = ".bed")
    x <- fread(file.path(outdir, f), header = FALSE, sep = "\t")
    setorder(x, V1, V2)
    fwrite(x, body, sep = "\t", col.names = FALSE, quote = FALSE)
    status <- system2(b2bb, c(paste0("-type=bed", ncol(x)), body,
                              sizes, file.path(outdir, sub("\\.bed$", ".bb", f))))
    if (status != 0) stop("bedToBigBed failed on ", f)
  }
  message("bigBeds written")
} else {
  message("No --bedtobigbed given; only the BED files were written")
}

# ── checks ────────────────────────────────────────────────────────────────────

# Each library's bigWigs, read back, must hold that library's filtered TPM.
bw <- manifest[kind == "tpm" & !is.na(file)]
readback <- vapply(bw$file, function(f) {
  g <- import.bw(file.path(outdir, f)); sum(score(g) * width(g))
}, numeric(1))
bad <- bw$file[abs(readback - bw$total) > 1e-3 * pmax(1, bw$total)]
if (length(bad)) stop("bigWig sums disagree with the manifest: ", paste(bad, collapse = ", "))
message("All ", nrow(bw), " bigWigs read back with sums matching the manifest")
print(manifest[, .(file, n_features, total = round(total, 1))], nrows = 100)

fwrite(manifest, manifest_f, sep = "\t")
message("Written: ", manifest_f)
sessionInfo()
