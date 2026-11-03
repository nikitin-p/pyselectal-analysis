#!/usr/bin/env Rscript
# Stage 21: stranded CTSS bigWigs from the subsampled libraries, one black track
# for the pooled total and one coloured track per 5' end type.
#
# Every value is a raw count of 5' ends pooled over the libraries in Stage 17's
# fine-typed CTSS table, which was read from the subsampled BAMs; those were
# brought to equal depth, so the pooled sum needs no further scaling. A type's
# track holds, at every position, the reads of that type there, whatever the
# other types at the same position do; whether a position is dominated by one
# type or mixed is left to be read off the hub by overlaying the tracks. Each
# read has exactly one type, so the type tracks add up to the black one. Types
# outside the vocabulary of scripts/utils/end_types.R are pooled into its last
# class ("other").
#
# Usage (from project root):
#   Rscript scripts/13_tracks/01_endtype_bigwigs.R \
#       [--ctss   results/dominance/fine_typed_ctss.tsv.gz] \
#       [--fai    data/sacCer3.fa.fai]                       \
#       [--outdir results/tracks_subsampled] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
ctss_f <- get_arg("--ctss",   "results/dominance/fine_typed_ctss.tsv.gz")
fai_f  <- get_arg("--fai",    "data/sacCer3.fa.fai")
outdir <- get_arg("--outdir", "results/tracks_subsampled")
force  <- "--force" %in% args

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
       else "scripts/utils/end_types.R")
OTHER <- setdiff(END_TYPES, END_TYPES_REAL)
stopifnot(length(OTHER) == 1)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_f <- file.path(outdir, "track_manifest.tsv")
if (!force && file.exists(manifest_f)) {
  message("Output exists, skipping (use --force to rerun): ", manifest_f)
  quit(save = "no", status = 0)
}

for (f in c(ctss_f, fai_f)) if (!file.exists(f)) stop("Input not found: ", f)

# ── input ─────────────────────────────────────────────────────────────────────

ctss <- fread(ctss_f)
stopifnot(all(c("chr", "pos", "strand", "sample_id", "type", "n") %in% names(ctss)))
message("Loaded ", nrow(ctss), " CTSS rows from ", uniqueN(ctss$sample_id),
        " libraries: ", paste(sort(unique(ctss$sample_id)), collapse = ", "))

ctss[, vtype := ifelse(type %in% END_TYPES_REAL, type, OTHER)]
by_type <- ctss[, .(depth = sum(n)), by = .(chr, strand, pos, vtype)]
total   <- by_type[, .(depth = sum(depth)), by = .(chr, strand, pos)]

# ── track order, labels and colours ───────────────────────────────────────────

rgb_str <- function(h) paste(as.integer(col2rgb(h)), collapse = ",")
track_meta <- rbind(
  data.table(track = "total", color = "0,0,0",
             description = "all 5' ends pooled over the subsampled libraries"),
  data.table(track = END_TYPES,
             color = vapply(END_TYPES, function(t) rgb_str(END_TYPE_COLORS_BRIGHT[[t]]), ""),
             description = ifelse(END_TYPES == OTHER,
                                  "5' ends of every type outside the end-type vocabulary",
                                  paste0("5' ends of type ", END_TYPES))))

# ── bigWigs ───────────────────────────────────────────────────────────────────

fai <- fread(fai_f, header = FALSE, select = 1:2, col.names = c("chr", "len"))
seqlens <- setNames(fai$len, fai$chr)
stopifnot(all(unique(total$chr) %in% names(seqlens)))

strands <- c(plus = "+", minus = "-")

write_bw <- function(dt, path) {
  setorder(dt, chr, pos)
  gr <- GRanges(dt$chr, IRanges(dt$pos, width = 1), score = as.numeric(dt$depth),
                seqinfo = Seqinfo(names(seqlens), seqlens))
  # The default zoom levels are derived from the data, and UCSC's bbiWrite fails
  # on them ("doesn't intersect ... Internal error bbiWrite.c 414") when a track
  # holds a handful of positions on different chromosomes; fixed zoom levels
  # avoid that and are used for every track so that all files share one layout.
  export.bw(sort(gr), path, fixedSummaries = TRUE)
}

# A failed earlier run can leave empty or partial files, and the tracks of the
# earlier dominance-based design would otherwise linger beside the new ones.
unlink(list.files(outdir, pattern = "\\.bw$", full.names = TRUE))

manifest <- rbindlist(lapply(seq_len(nrow(track_meta)), function(i) {
  tr <- track_meta$track[i]
  rbindlist(lapply(names(strands), function(sn) {
    sel <- if (tr == "total") total[strand == strands[[sn]]] else
      by_type[vtype == tr & strand == strands[[sn]]]
    # A rare type can be absent from one strand, and an empty bigWig cannot be
    # written; that row stays in the manifest with no file.
    file <- if (nrow(sel)) paste0(tr, "_", sn, ".bw") else NA_character_
    if (nrow(sel)) write_bw(copy(sel), file.path(outdir, file))
    data.table(file = file, track = tr, strand = strands[[sn]],
               color = track_meta$color[i], n_positions = nrow(sel),
               n_reads = sum(sel$depth), description = track_meta$description[i])
  }))
}))

# ── checks ────────────────────────────────────────────────────────────────────

# The type tracks sum, on each strand, to the black one and to the CTSS table.
for (s in strands) {
  tot  <- manifest[track == "total" & strand == s, n_reads]
  cols <- manifest[track != "total" & strand == s, sum(n_reads)]
  raw  <- ctss[strand == s, sum(n)]
  message(sprintf("strand %s: total %d, type tracks %d, CTSS table %d", s, tot, cols, raw))
  stopifnot(tot == cols, tot == raw)
}

# The files themselves, read back.
written  <- manifest[!is.na(file)]
readback <- vapply(written$file, function(f) {
  gr <- import.bw(file.path(outdir, f))
  sum(score(gr) * width(gr))
}, numeric(1))
bad <- written$file[abs(readback - written$n_reads) > 0.5]
if (length(bad)) stop("bigWig sums disagree with the manifest: ", paste(bad, collapse = ", "))
message("All ", nrow(written), " bigWigs read back with sums matching the manifest; ",
        sum(is.na(manifest$file)), " empty track-strand combinations were not written")

print(manifest[, .(file, n_positions, n_reads)], nrows = 200)

# The manifest is written last, so an interrupted run is not skipped next time.
fwrite(manifest, manifest_f, sep = "\t")
message("Written: ", manifest_f)
sessionInfo()
