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

# Track names must be unique across the whole hub.
ids <- sub("^ *track ", "", grep("^ *track ", lines, value = TRUE))
stopifnot(!anyDuplicated(ids), all(grepl("^[A-Za-z0-9_]+$", ids)))
files <- sub("^ *bigDataUrl ", "", grep("bigDataUrl", lines, value = TRUE))
stopifnot(!anyDuplicated(files), all(grepl("^[A-Za-z0-9._]+$", files)))

writeLines(lines, out_f)
message("Written: ", out_f, " (", prio - 1L, " composites, ", n_tracks, " tracks)")
