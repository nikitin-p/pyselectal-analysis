#!/usr/bin/env Rscript
# Stage 21 step 2: a UCSC hub.txt for the tracks built from the subsampled
# libraries, generated from the track manifests so that it can never name a
# file that does not exist or miss one that does.
#
# Three sets are drawn, each as one composite per strand (F = +, R = -), so
# that a whole set is switched on or off from one settings page:
#   endtype   Stage 21: the black pooled total first, then one coloured
#             track per 5' end type, in manifest order;
#   cager     Stage 10 step 4: power-law TPM per library, filtered CTSSs only;
#   longclip  Stage 22: clipped-base coverage and aligned starts per origin class.
#   nad       NAD branch step 11: NAD captureSeq S, N and S-N per strain, on
#             one shared scale within each strain and strand, plus the S-vs-N
#             calls as a bigBed.
# Two more sets come from files rather than manifests: the per-base sequence
# letters (sequence_F.bb, sequence_R.bb from make_hub.sh) on top, and the
# Stage 22 step 4 long-clip BAMs, coloured by their YC tag, at the bottom.
# Every bigDataUrl is the bare file name, because Hub Space stores the files of
# one hub side by side; the files of the three folders must therefore be
# uploaded into the same hub as this hub.txt. A set whose manifest is absent is
# left out with a message rather than failing.
#
# Usage (from project root):
#   Rscript scripts/13_tracks/02_make_hub.R \
#       [--endtype  results/tracks_subsampled/track_manifest.tsv] \
#       [--cager    results/cager/tracks/track_manifest.tsv]      \
#       [--longclip results/longclip/tracks/track_manifest.tsv]   \
#       [--longclip-bam-dir results/longclip/bam]                \
#       [--nad      nad/results/11_tracks/track_manifest.tsv]     \
#       [--sequence-dir results/tracks]                          \
#       [--out      results/tracks_subsampled/hub.txt]            \
#       [--genome sacCer3] [--hub yeastCageSubsampled] [--email ADDRESS]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}
endtype_f  <- get_arg("--endtype",  "results/tracks_subsampled/track_manifest.tsv")
cager_f    <- get_arg("--cager",    "results/cager/tracks/track_manifest.tsv")
longclip_f <- get_arg("--longclip", "results/longclip/tracks/track_manifest.tsv")
lc_bam_dir <- get_arg("--longclip-bam-dir", "results/longclip/bam")
nad_f      <- get_arg("--nad",      "nad/results/11_tracks/track_manifest.tsv")
seq_dir    <- get_arg("--sequence-dir", "results/tracks")
out_f      <- get_arg("--out",      "results/tracks_subsampled/hub.txt")
genome     <- get_arg("--genome",   "sacCer3")
hub_name   <- get_arg("--hub",      "yeastCageSubsampled")
email      <- get_arg("--email",    "nikipavel14@gmail.com")

# UCSC truncates a shortLabel beyond 17 characters in the track list.
short <- function(x) substr(x, 1, 17)
safe  <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
strand_name  <- c("+" = "F", "-" = "R")
strand_long  <- c("+" = "forward (+) strand", "-" = "reverse (-) strand")
# Strand colours of the original hub (make_hub.sh), for sets without their own.
strand_color <- c("+" = "33,102,172", "-" = "178,24,43")

read_manifest <- function(f, what) {
  if (!file.exists(f)) { message("No manifest for ", what, " (", f, "); set left out"); return(NULL) }
  m <- fread(f)
  m[!is.na(file) & nzchar(file)]
}

composite <- function(id, short_label, long_label, prio) {
  c(sprintf("track %s", id), "compositeTrack on", sprintf("shortLabel %s", short(short_label)),
    sprintf("longLabel %s", long_label), "type bigWig", "autoScale on",
    "dragAndDrop subTracks", "allButtonPair on", "visibility full",
    sprintf("priority %d", prio), "")
}
subtrack <- function(id, parent, file, short_label, long_label, color, prio, on = TRUE) {
  c(sprintf("    track %s", id), sprintf("    parent %s%s", parent, if (on) " on" else " off"),
    sprintf("    bigDataUrl %s", file), sprintf("    shortLabel %s", short(short_label)),
    sprintf("    longLabel %s", long_label), "    type bigWig",
    sprintf("    color %s", color), sprintf("    priority %d", prio), "")
}

lines <- c(sprintf("hub %s", hub_name),
           "shortLabel Yeast CAGE subsampled",
           "longLabel Yeast CAGE 5' end tracks built from the depth-equalised (subsampled) libraries",
           "useOneFile on", sprintf("email %s", email), "",
           sprintf("genome %s", genome), "")
prio <- 1L
n_tracks <- 0L

# ── per-base sequence letters (make_hub.sh / make_sequence_bed.py) ───────────
# sequence_R.bb is the complement in genome order, read right to left for the
# reverse strand; the forward one is off by default because UCSC already draws
# the + strand bases. Both come before every signal track.
seq_files <- c(F = "sequence_F.bb", R = "sequence_R.bb")
seq_have  <- file.exists(file.path(seq_dir, seq_files))
if (!any(seq_have)) {
  message("No sequence bigBeds in ", seq_dir, "; set left out")
} else {
  lines <- c(lines, "track sequence", "compositeTrack on", "shortLabel Sequence",
             "longLabel Per-base genome letters: forward strand and its complement",
             "type bigBed 4", "color 255,255,255", "labelOnFeature on",
             "visibility dense", "priority 0", "")
  seq_long <- c(F = "Genome sequence, forward strand (one letter per base)",
                R = "Complement of the genome sequence, in genome order; read right to left for the reverse strand 5' to 3'")
  for (g in names(seq_files)[seq_have]) {
    lines <- c(lines, sprintf("    track sequence_%s", g),
               sprintf("    parent sequence %s", if (g == "R") "on" else "off"),
               sprintf("    bigDataUrl %s", seq_files[[g]]),
               sprintf("    shortLabel %s", if (g == "F") "Sequence +" else "Complement -"),
               sprintf("    longLabel %s", seq_long[[g]]),
               "    type bigBed 4", "    color 255,255,255", "    labelOnFeature on",
               sprintf("    priority %d", if (g == "F") 1L else 2L), "")
    n_tracks <- n_tracks + 1L
  }
}

# ── end-type partition (Stage 21) ─────────────────────────────────────────────
et <- read_manifest(endtype_f, "end-type tracks")
if (!is.null(et)) for (s in c("+", "-")) {
  g <- strand_name[[s]]; id <- paste0("endtype_", g)
  lines <- c(lines, composite(id, paste("End types", g),
                              paste("Subsampled CTSS counts: black = pooled total, colour = reads of one 5' end type,",
                                    strand_long[[s]]), prio))
  m <- et[strand == s]
  for (i in seq_len(nrow(m))) {
    lab <- if (m$track[i] == "total") "total" else m$track[i]
    lines <- c(lines, subtrack(paste0(id, "_", safe(sub("\\.bw$", "", m$file[i]))), id, m$file[i],
                               paste(lab, g), paste0(m$description[i], ", ", strand_long[[s]]),
                               m$color[i], i))
    n_tracks <- n_tracks + 1L
  }
  prio <- prio + 1L
}

# ── CAGEr TPM per library (Stage 10) ─────────────────────────────────────────
cg <- read_manifest(cager_f, "CAGEr TPM tracks")
if (!is.null(cg)) for (s in c("+", "-")) {
  g <- strand_name[[s]]; id <- paste0("cager_tpm_", g)
  lines <- c(lines, composite(id, paste("CAGEr TPM", g),
                              paste("Power-law TPM per subsampled library, CTSSs passing the 1 TPM filter,",
                                    strand_long[[s]]), prio))
  m <- cg[kind == "tpm" & strand == s]
  for (i in seq_len(nrow(m))) {
    lab <- sub("^yeast_cage_", "", m$sample[i])
    lines <- c(lines, subtrack(paste0(id, "_", safe(m$sample[i])), id, m$file[i],
                               paste(lab, "TPM", g), paste(m$sample[i], "power-law TPM,", strand_long[[s]]),
                               strand_color[[s]], i))
    n_tracks <- n_tracks + 1L
  }
  prio <- prio + 1L
}

# ── long-clip origins (Stage 22) ──────────────────────────────────────────────
lc <- read_manifest(longclip_f, "long-clip tracks")
if (!is.null(lc)) for (s in c("+", "-")) {
  g <- strand_name[[s]]; id <- paste0("longclip_", g)
  lines <- c(lines, composite(id, paste("Long clips", g),
                              paste("Long 5' soft clips by origin class: clip = clipped bases at their origin,",
                                    "start = aligned 5' ends,", strand_long[[s]]), prio))
  m <- lc[strand == s]
  for (i in seq_len(nrow(m))) {
    lines <- c(lines, subtrack(paste0(id, "_", safe(sub("\\.bw$", "", m$file[i]))), id, m$file[i],
                               paste(m$kind[i], m$class[i]),
                               paste0(m$kind[i], " track, ", m$class[i], ", ", strand_long[[s]]),
                               m$color[i], i, on = m$kind[i] == "clip"))
    n_tracks <- n_tracks + 1L
  }
  prio <- prio + 1L
}

# ── NAD captureSeq (NAD branch step 11) ──────────────────────────────────────
# S, N and S-N of one strain share one y scale (autoScale group), so that the
# S and N heights are directly comparable; N and S-N are off by default.
nd <- read_manifest(nad_f, "NAD captureSeq tracks")
if (!is.null(nd)) {
  for (st in unique(nd$strain)) for (s in c("+", "-")) {
    g <- strand_name[[s]]; id <- paste0("nad_", safe(st), "_", g)
    lines <- c(lines, sprintf("track %s", id), "compositeTrack on",
               sprintf("shortLabel %s", short(paste("NAD", sub("dxo1_rai1_npy1", "dxo1rai1npy1", st), g))),
               sprintf("longLabel NAD captureSeq 5' ends, %s: S (+ADPRC), N (no ADPRC) and S-N, reads / step-5 size factor, %s",
                       st, strand_long[[s]]),
               "type bigWig", "autoScale group", "alwaysZero on", "dragAndDrop subTracks",
               "allButtonPair on", "visibility full", sprintf("priority %d", prio), "")
    m <- nd[kind == "signal" & strain == st & strand == s]
    for (i in seq_len(nrow(m))) {
      lines <- c(lines, subtrack(paste0(id, "_", m$track[i]), id, m$file[i],
                                 paste("NAD", if (m$track[i] == "diff") "S-N" else m$track[i], g),
                                 paste0(m$description[i], ", ", strand_long[[s]]),
                                 m$color[i], i, on = m$track[i] == "S"))
      n_tracks <- n_tracks + 1L
    }
    prio <- prio + 1L
  }
  for (i in which(nd$kind == "calls")) {
    lines <- c(lines, sprintf("track nad_calls_%s", safe(nd$strain[i])),
               sprintf("bigDataUrl %s", nd$file[i]),
               sprintf("shortLabel %s", short(paste("NAD calls", sub("dxo1_rai1_npy1", "triple", nd$strain[i])))),
               sprintf("longLabel %s", nd$description[i]),
               "type bigBed 9", "itemRgb on", "visibility pack", sprintf("priority %d", prio), "")
    n_tracks <- n_tracks + 1L
    prio <- prio + 1L
  }
}

# ── original STAR CIGARs with the soft clips drawn (Stage 22 step 8) ─────────
# UCSC never draws the S of a BAM record, so the clip is a BED12 block of its
# own: thin = soft clip where STAR left it, thick = aligned, line = N; the name
# is the original CIGAR and the colour the origin class.
cigar_bb <- "longclip_cigar.bb"
if (!file.exists(file.path(lc_bam_dir, cigar_bb))) {
  message("No ", cigar_bb, " in ", lc_bam_dir, "; CIGAR track left out")
} else {
  lines <- c(lines, "track longclip_cigar", sprintf("bigDataUrl %s", cigar_bb),
             "shortLabel Long-clip CIGAR",
             "longLabel Long-clip reads as STAR aligned them: thin = soft clip (S), thick = aligned (M), line = intron (N); label = CIGAR, colour = origin class",
             "type bigBed 12", "itemRgb on", "maxItems 20000",
             "visibility pack", sprintf("priority %d", prio), "")
  n_tracks <- n_tracks + 1L
  prio <- prio + 1L
}

# ── long-clip reads as =/X BAMs (Stage 22 step 4) ─────────────────────────────
# Both BAMs carry YC:Z: the origin-class colour per read in the realigned BAM,
# green / red / grey per match / mismatch / untemplated run in the segments
# BAM. UCSC uses it only under bamColorMode tag. The segments BAM holds several
# records per read and overflows the default item limit at a busy TSS, so the
# limit is raised; it is off by default because the realigned BAM is the one to
# read first. UCSC finds each .bai next to its BAM by name.
lc_bams <- data.table(
  file  = c("longclip_realigned.bam", "longclip_segments.bam"),
  short = c("LC realigned", "LC segments"),
  long  = c("Long-clip reads with the clip aligned at its origin; =/X CIGAR, coloured by origin class (YC)",
            "Long-clip reads split into match (green), mismatch (red) and untemplated (grey) runs; BX:Z = read name"),
  on    = c(TRUE, FALSE))
lc_bams <- lc_bams[file.exists(file.path(lc_bam_dir, file)) &
                   file.exists(file.path(lc_bam_dir, paste0(file, ".bai")))]
if (nrow(lc_bams) == 0) {
  message("No indexed long-clip BAMs in ", lc_bam_dir, "; set left out")
} else {
  lines <- c(lines, "track longclip_reads", "compositeTrack on",
             "shortLabel Long-clip reads",
             "longLabel Long-clip reads (Stage 22 step 4): the whole read compared with the genome, both strands",
             "type bam", "bamColorMode tag", "bamColorTag YC",
             "baseColorDefault diffBases", "baseColorUseSequence lfExtra",
             "showNames off", "maxItems 20000", "dragAndDrop subTracks",
             "allButtonPair on", "visibility squish", sprintf("priority %d", prio), "")
  for (i in seq_len(nrow(lc_bams))) {
    lines <- c(lines, sprintf("    track longclip_reads_%s", safe(sub("\\.bam$", "", lc_bams$file[i]))),
               sprintf("    parent longclip_reads %s", if (lc_bams$on[i]) "on" else "off"),
               sprintf("    bigDataUrl %s", lc_bams$file[i]),
               sprintf("    shortLabel %s", short(lc_bams$short[i])),
               sprintf("    longLabel %s", lc_bams$long[i]),
               "    type bam", "    bamColorMode tag", "    bamColorTag YC",
               "    maxItems 20000", sprintf("    priority %d", i), "")
    n_tracks <- n_tracks + 1L
  }
  prio <- prio + 1L
}

# Track names must be unique across the whole hub.
ids <- sub("^ *track ", "", grep("^ *track ", lines, value = TRUE))
stopifnot(!anyDuplicated(ids), all(grepl("^[A-Za-z0-9_]+$", ids)))
files <- sub("^ *bigDataUrl ", "", grep("bigDataUrl", lines, value = TRUE))
stopifnot(!anyDuplicated(files), all(grepl("^[A-Za-z0-9._]+$", files)))

writeLines(lines, out_f)
message("Written: ", out_f, " (", prio - 1L, " composites, ", n_tracks, " tracks)")
