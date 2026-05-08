#!/usr/bin/env Rscript
# Stage 10 step 1: build CAGEexp from typed_ctss.tsv.gz, normalize (powerLaw),
# filter low-expression CTSSs, cluster with distclu, compute IQW.
# Saves a CAGEexp RDS and a per-sample tag cluster TSV for downstream use.
#
# Usage (from project root):
#   Rscript scripts/06_cager/01_cageexp.R \
#       [--ctss    results/tss/typed_ctss.tsv.gz] \
#       [--samples config/samples.tsv]             \
#       [--outdir  results/cager]                  \
#       [--params  config/params.yaml]             \
#       [--force]

suppressPackageStartupMessages({
  library(CAGEr)
  library(data.table)
  library(readr)
  library(yaml)
  library(methods)
  library(GenomicRanges)
})

# ── helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x[[1]]) && x[[1]] != "") x else y

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

ctss_f    <- get_arg("--ctss",    "results/tss/typed_ctss.tsv.gz")
samples_f <- get_arg("--samples", "config/samples.tsv")
outdir    <- get_arg("--outdir",  "results/cager")
params_f  <- get_arg("--params",  "config/params.yaml")
force     <- "--force" %in% args

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(outdir, "cageexp.rds")
out_tc  <- file.path(outdir, "tagclusters.tsv.gz")

if (!force && file.exists(out_rds) && file.exists(out_tc)) {
  message("Outputs exist, skipping (use --force to rerun): ", out_rds)
  quit(save = "no", status = 0)
}

if (!file.exists(ctss_f)) {
  stop("Input not found: ", ctss_f,
       "\nRun stage 8 first (scripts/05_tss_end_types/run.sh).")
}

params <- yaml.load_file(params_f)

tpm_threshold   <- as.numeric(params$cager_tpm_threshold   %||% 1)
max_dist        <- as.integer(params$cager_max_dist        %||% 20)
keep_singletons <- as.integer(params$cager_keep_singletons_above %||% 5)
iq_low          <- as.numeric(params$cager_iqw_qlow        %||% 0.1)
iq_up           <- as.numeric(params$cager_iqw_qup         %||% 0.9)

message("Stage 10 step 1: build CAGEexp")
message("  distclu maxDist=", max_dist,
        "  keepSingletonsAbove=", keep_singletons,
        "  TPM filter=", tpm_threshold)

# ── load + aggregate CTSS ─────────────────────────────────────────────────────
message("  loading: ", ctss_f)
ctss <- fread(cmd = paste("zcat", shQuote(ctss_f)))
message("  rows: ", nrow(ctss))

# Sum typed counts to total per (chr, pos, strand, sample)
total <- ctss[, .(count = sum(n_total)), by = .(chr, pos, strand, sample_id)]

# ── build wide data.frame for CAGEexp ────────────────────────────────────────
wide <- dcast(total, chr + pos + strand ~ sample_id,
              value.var = "count", fun.aggregate = sum, fill = 0L)

sample_cols <- setdiff(names(wide), c("chr", "pos", "strand"))
for (j in sample_cols) wide[[j]] <- as.integer(wide[[j]])
wide$chr    <- as.character(wide$chr)
wide$pos    <- as.integer(wide$pos)
wide$strand <- as.character(wide$strand)

df_all <- as.data.frame(wide, stringsAsFactors = FALSE)

# CAGEr requires: integer pos > 0, strand %in% {+,-}, integer sample columns
df_all <- df_all[!is.na(df_all$pos) & df_all$pos > 0 &
                   df_all$strand %in% c("+", "-"), ]
stopifnot(
  is.integer(df_all$pos),
  all(df_all$strand %in% c("+", "-")),
  all(vapply(df_all[sample_cols], is.integer, TRUE))
)

message("  CTSS positions: ", nrow(df_all),
        "  samples: ", length(sample_cols))

# ── build CAGEexp ─────────────────────────────────────────────────────────────
ce <- as(df_all, "CAGEexp")

message("  library sizes:")
print(colData(ce)[, "librarySizes", drop = FALSE])

# ── normalize + filter ────────────────────────────────────────────────────────
message("  normalizing (powerLaw, T=1e6)...")
ce <- normalizeTagCount(ce, method = "powerLaw", T = 1e6)

message("  filtering low-expression CTSSs (TPM >= ", tpm_threshold, ")...")
ce <- filterLowExpCTSS(ce,
                       thresholdIsTpm    = TRUE,
                       nrPassThreshold   = 1,
                       threshold         = tpm_threshold)

# ── cluster ───────────────────────────────────────────────────────────────────
message("  distclu (maxDist=", max_dist, ")...")
ce <- distclu(ce, maxDist = max_dist, keepSingletonsAbove = keep_singletons)

# sampleLabels<- GRanges workaround (needed in some CAGEr versions)
has_sl_gr <- !is.null(
  suppressWarnings(tryCatch(
    methods::getMethod("sampleLabels<-", "GRanges"),
    error = function(e) NULL
  ))
)
if (!has_sl_gr) {
  if (!methods::isGeneric("sampleLabels<-"))
    setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))
  setMethod("sampleLabels<-", "GRanges", function(object, value) {
    S4Vectors::metadata(object)$sampleLabels <- value
    object
  })
}

# ── IQW ───────────────────────────────────────────────────────────────────────
message("  cumulative CTSS distribution...")
ce <- cumulativeCTSSdistribution(ce, clusters = "tagClusters",
                                 useMulticore = FALSE)

message("  quantile positions (qLow=", iq_low, " qUp=", iq_up, ")...")
ce <- quantilePositions(ce, clusters = "tagClusters",
                        qLow = iq_low, qUp = iq_up)

# ── save RDS ──────────────────────────────────────────────────────────────────
message("  saving RDS: ", out_rds)
saveRDS(ce, out_rds)

# ── export tag clusters TSV ───────────────────────────────────────────────────
message("  exporting tag clusters: ", out_tc)

tc_list <- lapply(sampleLabels(ce), function(sid) {
  gr <- tryCatch(tagClustersGR(ce, sample = sid), error = function(e) NULL)
  if (is.null(gr) || length(gr) == 0) return(NULL)

  dt <- data.table(
    sample_id      = sid,
    chr            = as.character(seqnames(gr)),
    start          = start(gr),
    end            = end(gr),
    strand         = as.character(strand(gr)),
    dominant_ctss  = if (!is.null(mcols(gr)$dominant_ctss)) mcols(gr)$dominant_ctss else NA_integer_,
    tpm            = if (!is.null(mcols(gr)$tpm.dominant))  mcols(gr)$tpm.dominant  else NA_real_
  )

  # IQW: try common column name variants across CAGEr versions
  qlo_col <- grep("^q_?0?\\.?1$|^q_low|^qLow", names(mcols(gr)), value = TRUE)[1]
  qup_col <- grep("^q_?0?\\.?9$|^q_up|^qUp",  names(mcols(gr)), value = TRUE)[1]
  if (!is.na(qlo_col) && !is.na(qup_col)) {
    dt$iq_width <- as.integer(mcols(gr)[[qup_col]]) -
                   as.integer(mcols(gr)[[qlo_col]]) + 1L
  } else {
    dt$iq_width <- NA_integer_
  }

  dt
})

tc_all <- rbindlist(Filter(Negate(is.null), tc_list))
fwrite(tc_all, out_tc, sep = "\t", compress = "gzip")
message("  tag clusters written: ", nrow(tc_all), " rows")

sessionInfo()
