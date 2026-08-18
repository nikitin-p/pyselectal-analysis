#!/usr/bin/env Rscript
# Stage 8 step 8: within-position analysis of the added-G ladder (1Sg / 2Sgg / 3Sggg)
#
# Motivation
# ----------
# 2Sgg and 3Sggg almost never dominate a TSS cluster, so dominance-based selection
# (06_motif_analysis.R) cannot study them. They do, however, overwhelmingly occur AT
# THE SAME genomic positions as 1Sg: 2,651 positions carry all three types and hold
# 79% of all 3Sggg reads. At such a position the genomic sequence is identical for
# every read, so a motif comparison between types is uninformative by construction.
#
# This script therefore treats the position as its own control and models the number
# of non-templated G's appended to a read as a per-position rate:
#
#     k = 1  ->  n_1Sg      k = 2  ->  n_2Sgg      k = 3  ->  n_3Sggg
#
# The hypothesis under test is that the extra G's are template-free additions by
# reverse transcriptase (the terminal-transferase activity that template switching
# exploits). That predicts (a) a geometric ladder, p(k+1)/p(k) constant, and (b) a
# per-step probability that does not depend on the surrounding genomic sequence.
# Systematic sequence or condition dependence would argue against a pure RT slip.
#
# Note on the alignment constraint: every position in the ladder set carries at least
# one soft-clipped-G read, which is only possible when the genomic base at -1 is not G
# (otherwise the aligner would have extended the match). This constraint is IDENTICAL
# for 1Sg, 2Sgg and 3Sggg, so it cancels within the ladder and does not confound the
# comparison. It is the reason base_m1 below only ever takes levels A/C/T.
#
# M-type reads are reported alongside for reference but are excluded from the ladder:
# they carry no soft clip, so they are subject to a different alignment constraint and
# are additionally confounded by uncapped and degraded RNA.
#
# Usage:
#   Rscript scripts/05_tss_end_types/08_g_ladder.R \
#       [--ctss      results/tss/typed_ctss.tsv.gz] \
#       [--params    config/params.yaml]            \
#       [--outdir    results/g_ladder]              \
#       [--figdir    results/figures]               \
#       [--min_depth 50]                            \
#       [--force]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(data.table)
  library(ggplot2)
  library(readr)
  library(yaml)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

.val_flags  <- character(0)
.bool_flags <- character(0)

get_arg <- function(flag, default = NULL) {
  .val_flags <<- union(.val_flags, flag)
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

has_flag <- function(flag) {
  .bool_flags <<- union(.bool_flags, flag)
  flag %in% args
}

validate_args <- function() {
  unknown <- character(0)
  i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% .val_flags)  { i <- i + 2L; next }
    if (args[i] %in% .bool_flags) { i <- i + 1L; next }
    unknown <- c(unknown, args[i])
    i <- i + 1L
  }
  if (length(unknown) > 0) {
    stop("Unrecognized argument(s): ", paste(unknown, collapse = " "),
         "\n  Known flags: ", paste(sort(c(.val_flags, .bool_flags)), collapse = " "),
         call. = FALSE)
  }
}

ctss_f       <- get_arg("--ctss",      "results/tss/typed_ctss.tsv.gz")
params_f     <- get_arg("--params",    "config/params.yaml")
outdir       <- get_arg("--outdir",    "results/g_ladder")
figdir       <- get_arg("--figdir",    "results/figures")
min_depth    <- as.integer(get_arg("--min_depth", "50"))
# Positions above this extra-G rate are a qualitatively different population (see
# section B'): they carry almost no 1Sg reads at all, so they are not 1Sg TSS that
# occasionally slip. They are reported separately rather than averaged in.
outlier_rate <- as.numeric(get_arg("--outlier_rate", "0.25"))
force        <- has_flag("--force")

validate_args()

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)

out_main <- file.path(outdir, "ladder_global.tsv")
if (!force && file.exists(out_main)) {
  message("Output exists, skipping (use --force to rerun): ", out_main)
  quit(save = "no", status = 0)
}

YR_RED   <- "#e63946"   # canonical initiator colour (project convention)
OTHER_GY <- "#6c757d"

# fmt_n(): compact count labels drawn on top of every histogram bar
.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

# ── genome loading ────────────────────────────────────────────────────────────
load_genome <- function(genome_spec) {
  if (grepl("\\.(fa|fasta)(\\.gz)?$", genome_spec, ignore.case = TRUE)) {
    if (!file.exists(genome_spec)) stop("FASTA not found: ", genome_spec)
    if (!file.exists(paste0(genome_spec, ".fai"))) Rsamtools::indexFa(genome_spec)
    fa <- Rsamtools::FaFile(genome_spec)
    open(fa)
    return(fa)
  }
  pkg <- switch(genome_spec,
    hg38    = "BSgenome.Hsapiens.UCSC.hg38",
    sacCer3 = "BSgenome.Scerevisiae.UCSC.sacCer3",
    stop("Unsupported genome alias: ", genome_spec)
  )
  if (!requireNamespace(pkg, quietly = TRUE)) stop("BSgenome not installed: ", pkg)
  get(pkg, envir = asNamespace(pkg))
}

harmonize_seqlevels <- function(gr, genome) {
  gi <- seqinfo(genome)
  common <- intersect(seqlevels(gr), seqlevels(gi))
  if (length(common) == 0) stop("No common seqlevels between GRanges and genome")
  gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
  seqinfo(gr) <- gi[seqlevels(gr)]
  gr
}

# Extract a window around each TSS, oriented 5'->3' along the transcript.
# getSeq() already reverse-complements minus-strand ranges; do NOT revcomp again.
extract_context <- function(dt, genome, up, down) {
  gr <- GRanges(dt$chr, IRanges(dt$pos, width = 1), strand = dt$strand)
  gr <- harmonize_seqlevels(gr, genome)
  start(gr) <- ifelse(strand(gr) == "+", start(gr) - up,   start(gr) - down)
  end(gr)   <- ifelse(strand(gr) == "+", end(gr)   + down, end(gr)   + up)
  gr <- trim(gr)
  ok <- width(gr) == up + 1 + down
  list(seqs = getSeq(genome, gr[ok]), keep = ok)
}

# ── load and aggregate ────────────────────────────────────────────────────────
message("Loading typed CTSS: ", ctss_f)
ctss <- fread(cmd = paste("gzip -cd", shQuote(ctss_f)))
message("  rows: ", nrow(ctss), "   samples: ", uniqueN(ctss$sample_id))

ladder_cols <- c("n_1Sg", "n_2Sgg", "n_3Sggg")
stopifnot(all(ladder_cols %in% names(ctss)))
stopifnot("n_2Stg" %in% names(ctss))

agg <- ctss[, c(lapply(.SD, sum)),
            .SDcols = c(ladder_cols, "n_2Stg", "n_M", "n_total"),
            by = .(chr, pos, strand)]
agg[, depth := n_1Sg + n_2Sgg + n_3Sggg]
message("  positions with any ladder read: ", format(sum(agg$depth > 0), big.mark = ","))

lad <- agg[depth >= min_depth]
message("  positions with ladder depth >= ", min_depth, ": ", format(nrow(lad), big.mark = ","),
        "  (", format(sum(lad$depth), big.mark = ","), " reads)")
if (nrow(lad) < 100) stop("Too few positions at min_depth=", min_depth)

# extra-G rate: fraction of ladder reads carrying MORE than one non-templated G
lad[, x_extra := n_2Sgg + n_3Sggg]
lad[, rate_extra := x_extra / depth]
lad[, rate_2 := n_2Sgg / depth]
lad[, rate_3 := n_3Sggg / depth]

# Split off the outlier population before anything is averaged.
lad[, is_outlier := rate_extra > outlier_rate]
clean <- lad[!(is_outlier)]

message("  outlier positions (extra-G rate > ", outlier_rate, "): ", sum(lad$is_outlier),
        "  holding ", round(100 * sum(lad[(is_outlier)]$x_extra) / sum(lad$x_extra), 1),
        "% of all extra-G reads")

# ── A. global ladder and geometric fit ────────────────────────────────────────
ladder_of <- function(d, label) {
  tot <- c(sum(d$n_1Sg), sum(d$n_2Sgg), sum(d$n_3Sggg))
  p <- tot / sum(tot)
  data.frame(
    set        = label,
    n_G        = 1:3,
    # The ladder proper is the pure G run. 2Stg abuts the genome with a G but
    # carries a 5' T, so it is drawn alongside the ladder rather than inside it:
    # folding it in would inflate the second rung with a different species.
    type       = c("1Sg", "2Sgg", "3Sggg"),
    reads      = as.integer(tot),
    proportion = round(p, 6),
    step_ratio = c(NA, round(p[2] / p[1], 5), round(p[3] / p[2], 5))
  )
}
global_tab <- rbind(ladder_of(lad, "all"), ladder_of(clean, "outliers removed"))
write_tsv(global_tab, out_main)

# 2Stg reference rate, reported alongside the ladder but never mixed into it:
# it is a separate row rather than a fourth "n_G" rung (see the note in
# ladder_of() above). Expressed as a share of the same ladder depth so it is
# directly comparable to the ladder's own proportions.
stg_of <- function(d, label) {
  data.frame(
    set        = label,
    type       = "2Stg",
    reads      = as.integer(sum(d$n_2Stg)),
    proportion = round(sum(d$n_2Stg) / sum(d$depth), 6)
  )
}
stg_tab <- rbind(stg_of(lad, "all"), stg_of(clean, "outliers removed"))
write_tsv(stg_tab, file.path(outdir, "ladder_2Stg_reference.tsv"))

geo_of <- function(d) {
  q1 <- sum(d$n_2Sgg) / sum(d$n_1Sg)
  q2 <- sum(d$n_3Sggg) / sum(d$n_2Sgg)
  c(q1 = q1, q2 = q2, ratio = q2 / q1)
}
g_all <- geo_of(lad); g_cln <- geo_of(clean)

cat("\n=== A. Global added-G ladder (positions with depth >=", min_depth, ") ===\n")
print(global_tab, row.names = FALSE)
cat("\n  Step ratios (see the WARNING below before interpreting these):\n")
cat(sprintf("    all positions      p2/p1=%.5f  p3/p2=%.5f  ratio=%.2f\n",
            g_all[["q1"]], g_all[["q2"]], g_all[["ratio"]]))
cat(sprintf("    outliers removed   p2/p1=%.5f  p3/p2=%.5f  ratio=%.2f\n",
            g_cln[["q1"]], g_cln[["q2"]], g_cln[["ratio"]]))
cat("\n  WARNING: do NOT read ratio != 1 as evidence that the additions are\n")
cat("  non-independent. Two reasons:\n")
cat("   (1) The first G is cap-derived while the extra G's come from terminal\n")
cat("       transferase, so p2/p1 and p3/p2 measure different mechanisms and\n")
cat("       there is no reason for them to be equal in the first place.\n")
cat("   (2) Pooling positions with different rates inflates p3/p2 relative to\n")
cat("       p2/p1 all by itself. Section E stratifies to remove this.\n")

# ── B. is the per-position rate constant? (overdispersion) ────────────────────
# Under a single global rate p_extra, x_extra ~ Binomial(depth, p_extra).
disp_of <- function(d, label) {
  pe <- sum(d$x_extra) / sum(d$depth)
  ev <- d$depth * pe * (1 - pe)
  chi <- sum((d$x_extra - d$depth * pe)^2 / ev)
  data.frame(
    set                 = label,
    global_extra_G_rate = round(pe, 6),
    n_positions         = nrow(d),
    pearson_chisq       = round(chi, 1),
    dispersion          = round(chi / (nrow(d) - 1), 3),
    p_value             = format.pval(pchisq(chi, nrow(d) - 1, lower.tail = FALSE), digits = 3)
  )
}
disp_tab <- rbind(disp_of(lad, "all"), disp_of(clean, "outliers removed"))
write_tsv(disp_tab, file.path(outdir, "dispersion.tsv"))

p_extra <- sum(clean$x_extra) / sum(clean$depth)
set.seed(42)
sim_rate <- rbinom(nrow(clean), clean$depth, p_extra) / clean$depth

cat("\n=== B. Is the extra-G rate the same at every position? ===\n")
print(disp_tab, row.names = FALSE)
cat("  dispersion ~1 => a single RT rate explains everything;",
    "\n  dispersion >>1 => positions genuinely differ.\n")

# ── B'. the outlier population ────────────────────────────────────────────────
out_pos <- lad[(is_outlier)][order(-x_extra)]
write_tsv(out_pos, file.path(outdir, "outlier_positions.tsv"))
cat("\n=== B'. Outlier positions (extra-G rate >", outlier_rate, ") ===\n")
cat(sprintf("  %d positions (%.2f%% of all) hold %.1f%% of all extra-G reads\n",
            nrow(out_pos), 100 * nrow(out_pos) / nrow(lad),
            100 * sum(out_pos$x_extra) / sum(lad$x_extra)))
cat(sprintf("  median 1Sg reads there = %.0f  (vs %.0f at typical positions)\n",
            median(out_pos$n_1Sg), median(clean$n_1Sg)))
cat("  These are NOT 1Sg TSS that occasionally slip - they carry almost no 1Sg at all,\n")
cat("  so they are a separate phenomenon and are excluded from the model below.\n")
cat("  Full list written to outlier_positions.tsv for follow-up.\n")

# ── C. sequence dependence ────────────────────────────────────────────────────
genome_spec <- params$genome_fasta_yeast
if (is.null(genome_spec) || genome_spec == "") {
  genome_spec <- if (file.exists("data/sacCer3.fa")) "data/sacCer3.fa" else "sacCer3"
}
message("\nLoading genome: ", genome_spec)
genome <- load_genome(genome_spec)

UP <- 5L; DOWN <- 5L
ctx <- extract_context(clean, genome, UP, DOWN)   # outliers excluded, see B'
lad_ctx <- clean[ctx$keep]
seqs <- ctx$seqs
message("  context extracted for ", length(seqs), " positions")

offsets <- setdiff(-UP:DOWN, 0L)   # position 0 is the TSS base itself
for (o in c(offsets, 0L)) {
  col <- paste0("base_", ifelse(o < 0, paste0("m", abs(o)), paste0("p", o)))
  lad_ctx[[col]] <- substr(as.character(seqs), o + UP + 1L, o + UP + 1L)
}
lad_ctx[, initiator := paste0(base_m1, base_p0)]
lad_ctx[, is_YR := initiator %in% c("CA", "CG", "TA", "TG")]

cat("\n=== C. Sequence dependence of the extra-G rate ===\n")
cat("  base at -1 takes levels:", paste(sort(unique(lad_ctx$base_m1)), collapse = "/"),
    "(no G, as forced by the soft-clip constraint)\n")

base_cols <- grep("^base_", names(lad_ctx), value = TRUE)
for (bc in base_cols) lad_ctx[[bc]] <- factor(lad_ctx[[bc]], levels = c("A", "C", "G", "T"))

fml <- as.formula(paste("cbind(x_extra, depth - x_extra) ~", paste(base_cols, collapse = " + ")))
# Quasibinomial, not binomial: section B showed the counts are heavily overdispersed
# (positions differ far more than a single rate allows). A plain binomial would
# understate every standard error by ~sqrt(dispersion) and report p = 0 throughout.
fit <- glm(fml, family = quasibinomial(), data = lad_ctx)

co <- summary(fit)$coefficients
glm_tab <- data.frame(
  term      = rownames(co),
  estimate  = round(co[, 1], 4),
  odds_ratio = round(exp(co[, 1]), 3),
  std_error = round(co[, 2], 4),
  p_value   = co[, 4],
  row.names = NULL
)
glm_tab <- glm_tab[order(glm_tab$p_value), ]
write_tsv(glm_tab, file.path(outdir, "glm_sequence_coefficients.tsv"))

cat("\n  Binomial GLM, extra-G vs single-G, on flanking bases (top 12 by p):\n")
print(head(glm_tab, 12), row.names = FALSE)

null_dev <- fit$null.deviance; res_dev <- fit$deviance
cat(sprintf("\n  deviance explained by sequence: %.2f%%  (null %.0f -> residual %.0f)\n",
            100 * (null_dev - res_dev) / null_dev, null_dev, res_dev))

# per-initiator summary
init_tab <- lad_ctx[, .(n_positions = .N,
                        reads       = sum(depth),
                        extra_G_rate = sum(x_extra) / sum(depth)),
                    by = .(initiator, is_YR)][order(-reads)]
write_tsv(init_tab, file.path(outdir, "rate_by_initiator.tsv"))
cat("\n  Extra-G rate by initiator dinucleotide (top 10 by depth):\n")
print(head(init_tab, 10), row.names = FALSE)

# ── D. condition / replicate consistency ──────────────────────────────────────
per_sample <- ctss[, .(n_1Sg = sum(n_1Sg), n_2Sgg = sum(n_2Sgg), n_3Sggg = sum(n_3Sggg)),
                   by = sample_id]
per_sample[, depth := n_1Sg + n_2Sgg + n_3Sggg]
per_sample[, extra_G_rate := (n_2Sgg + n_3Sggg) / depth]
per_sample[, step_p2_p1 := n_2Sgg / n_1Sg]
per_sample[, step_p3_p2 := n_3Sggg / n_2Sgg]
per_sample[, condition := sub("^yeast_cage_([a-z]+)_r[0-9]+$", "\\1", sample_id)]
setorder(per_sample, sample_id)
write_tsv(per_sample, file.path(outdir, "rate_by_sample.tsv"))

cat("\n=== D. Per-sample consistency ===\n")
print(per_sample[, .(sample_id, condition, depth, extra_G_rate = round(extra_G_rate, 5),
                     step_p2_p1 = round(step_p2_p1, 5), step_p3_p2 = round(step_p3_p2, 5))],
      row.names = FALSE)

# ── E. does position heterogeneity explain the step ratios? ───────────────────
# Strata are defined on replicate 1 and the ratios measured on replicate 2, so
# positions are not selected on the same counts used to score them.
rep_dt <- copy(ctss)
rep_dt[, rep := sub(".*_(r[0-9]+)$", "\\1", sample_id)]
strat_tab <- NULL
if (uniqueN(rep_dt$rep) == 2) {
  w <- dcast(rep_dt, chr + pos + strand ~ rep,
             value.var = ladder_cols, fun.aggregate = sum)
  r <- sort(unique(rep_dt$rep))
  w[, d1 := rowSums(.SD), .SDcols = paste0(ladder_cols, "_", r[1])]
  w[, d2 := rowSums(.SD), .SDcols = paste0(ladder_cols, "_", r[2])]
  w <- w[d1 >= min_depth & d2 >= min_depth]
  w[, rate1 := (get(paste0("n_2Sgg_", r[1])) + get(paste0("n_3Sggg_", r[1]))) / d1]
  w <- w[rate1 <= outlier_rate]
  w[, stratum := cut(rate1, c(-1e-9, .005, .01, .02, .04, .08, outlier_rate),
                     labels = c("0-0.5%", "0.5-1%", "1-2%", "2-4%", "4-8%",
                                paste0("8-", 100 * outlier_rate, "%")))]
  strat_tab <- w[, .(positions = .N,
                     p2_p1 = sum(get(paste0("n_2Sgg_",  r[2]))) / sum(get(paste0("n_1Sg_",  r[2]))),
                     p3_p2 = sum(get(paste0("n_3Sggg_", r[2]))) / sum(get(paste0("n_2Sgg_", r[2])))),
                 by = stratum][order(stratum)]
  strat_tab[, ratio := p3_p2 / p2_p1]
  write_tsv(strat_tab, file.path(outdir, "step_ratio_by_stratum.tsv"))

  cat("\n=== E. Step ratios within homogeneous strata ===\n")
  cat("  (strata from replicate 1, ratios measured in replicate 2)\n")
  print(strat_tab[, .(stratum, positions, p2_p1 = round(p2_p1, 4),
                      p3_p2 = round(p3_p2, 4), ratio = round(ratio, 2))], row.names = FALSE)
  cat("\n  If 'ratio' varies strongly across strata and crosses 1.0, the pooled\n")
  cat("  ratio in section A is a mixture artifact of position heterogeneity and\n")
  cat("  carries no information about dependence between successive additions.\n")
  cat("  Note that p2/p1 varies far more across strata than p3/p2 does: positions\n")
  cat("  differ mainly in WHETHER extra-G addition happens at all, much less in\n")
  cat("  whether it continues once started.\n")
} else {
  message("  (section E skipped: expected exactly 2 replicates per condition)")
}

# ── figures ───────────────────────────────────────────────────────────────────
message("\nWriting figures...")

# 13a: observed vs binomial-null spread of per-position rates
dens_df <- rbind(
  data.frame(rate = clean$rate_extra, src = "observed"),
  data.frame(rate = sim_rate,         src = "binomial null (single RT rate)")
)
p1 <- ggplot(dens_df, aes(x = rate, fill = src, colour = src)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  scale_fill_manual(values = c("observed" = YR_RED, "binomial null (single RT rate)" = OTHER_GY)) +
  scale_colour_manual(values = c("observed" = YR_RED, "binomial null (single RT rate)" = OTHER_GY)) +
  theme_minimal(base_size = 10) +
  labs(x = "extra-G rate at a position  ((2Sgg + 3Sggg) / ladder depth)", y = "Density",
       fill = NULL, colour = NULL,
       title = "Per-position extra-G rate vs a single-rate binomial null",
       subtitle = sprintf("%s positions, ladder depth >= %d, outliers removed;  dispersion = %.1f",
                          format(nrow(clean), big.mark = ","), min_depth,
                          disp_tab$dispersion[disp_tab$set == "outliers removed"])) +
  theme(legend.position = "bottom")
ggsave(file.path(figdir, "13_gladder_dispersion.pdf"), p1, width = 8, height = 4.5)

# 13b: rate vs depth (detection-bias control)
p2 <- ggplot(clean, aes(x = depth, y = rate_extra)) +
  geom_point(alpha = 0.15, size = 0.5, colour = OTHER_GY) +
  geom_smooth(method = "gam", formula = y ~ s(log10(x)), colour = YR_RED, linewidth = 0.7) +
  geom_hline(yintercept = p_extra, linetype = "dashed") +
  scale_x_log10() +
  theme_minimal(base_size = 10) +
  labs(x = "Ladder depth at position (log10)", y = "extra-G rate",
       title = "Extra-G rate vs sequencing depth",
       subtitle = "A flat trend means the rate estimate is not driven by detection bias")
ggsave(file.path(figdir, "13_gladder_rate_vs_depth.pdf"), p2, width = 8, height = 4.5)

# 13c: rate by initiator, YR highlighted per project convention
init_plot <- init_tab[reads >= 1000]
init_plot[, initiator := factor(initiator, levels = initiator[order(-extra_G_rate)])]
p3 <- ggplot(init_plot, aes(x = initiator, y = extra_G_rate, fill = is_YR)) +
  geom_col() +
  geom_hline(yintercept = p_extra, linetype = "dashed") +
  # after the reference line, so a label that lands on it stays readable
  geom_text(aes(label = fmt_n(reads)), vjust = -0.4, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = c("TRUE" = YR_RED, "FALSE" = OTHER_GY),
                    labels = c("TRUE" = "YR (canonical)", "FALSE" = "other"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 10) +
  labs(x = "Initiator dinucleotide (-1, +1)", y = "extra-G rate",
       title = "Extra-G rate by initiator dinucleotide",
       subtitle = paste("Dashed line = global rate; initiators with >= 1000 ladder alignments.",
                        "Labels give the ladder alignments per initiator.")) +
  theme(legend.position = "bottom")
ggsave(file.path(figdir, "13_gladder_by_initiator.pdf"), p3, width = 8, height = 4.5)

# 13d: per-sample ladder step ratios
ps_long <- melt(per_sample, id.vars = c("sample_id", "condition"),
                measure.vars = c("step_p2_p1", "step_p3_p2"),
                variable.name = "step", value.name = "ratio")
ps_n <- per_sample[, .(sample_id, condition, depth,
                       ytop = pmax(step_p2_p1, step_p3_p2, na.rm = TRUE))]

p4 <- ggplot(ps_long, aes(x = sample_id, y = ratio, fill = step)) +
  geom_col(position = "dodge") +
  geom_text(data = ps_n, aes(x = sample_id, y = ytop, label = fmt_n(depth)),
            inherit.aes = FALSE, vjust = -0.6, size = 2.3, colour = "grey25") +
  facet_wrap(~condition, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c(step_p2_p1 = YR_RED, step_p3_p2 = OTHER_GY),
                    labels = c("2Sgg / 1Sg", "3Sggg / 2Sgg"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal(base_size = 10) +
  labs(x = NULL, y = "Step ratio",
       title = "Added-G ladder step ratios across samples and conditions",
       subtitle = paste("Equal bars within a sample support a geometric (constant-probability) process.",
                        "Labels give the ladder depth per sample.")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
# 13e: the global ladder itself (1Sg -> 2Sgg -> 3Sggg), with the 2Stg reference
# rate drawn as its own dashed series rather than folded into the ladder line.
ladder_plot_df <- subset(global_tab, set == "outliers removed")
ladder_plot_df$type <- factor(ladder_plot_df$type, levels = c("1Sg", "2Sgg", "3Sggg"))
stg_rate <- stg_tab$proportion[stg_tab$set == "outliers removed"]

p5 <- ggplot(ladder_plot_df, aes(x = n_G, y = proportion)) +
  geom_line(aes(colour = "ladder (1Sg / 2Sgg / 3Sggg)"), linewidth = 0.8) +
  geom_point(aes(colour = type), size = 2.5) +
  geom_hline(aes(yintercept = stg_rate, colour = "2Stg (not part of the ladder)"),
             linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(
    name   = NULL,
    values = c(
      setNames(unname(END_TYPE_COLORS_TAB10[levels(ladder_plot_df$type)]),
                levels(ladder_plot_df$type)),
      "ladder (1Sg / 2Sgg / 3Sggg)"    = "grey40",
      "2Stg (not part of the ladder)"  = unname(END_TYPE_COLORS_TAB10["2Stg"])
    )
  ) +
  scale_x_continuous(breaks = 1:3, labels = levels(ladder_plot_df$type)) +
  theme_minimal(base_size = 10) +
  labs(x = NULL, y = "Proportion of ladder depth",
       title = "Added-G ladder vs the 2Stg reference rate",
       subtitle = "2Stg abuts the genome with G but carries a 5' T, so it is not a fourth rung") +
  theme(legend.position = "bottom")
ggsave(file.path(figdir, "13_gladder_global.pdf"), p5, width = 7, height = 5)

# ── per-position table ────────────────────────────────────────────────────────
keep_cols <- c("chr", "pos", "strand", "n_1Sg", "n_2Sgg", "n_3Sggg", "n_M", "depth",
               "x_extra", "rate_extra", "rate_2", "rate_3", "initiator", "is_YR", base_cols)
pos_out <- file.path(outdir, "position_rates.tsv.gz")
write_tsv(lad_ctx[, ..keep_cols], pos_out)
message("  saved: ", pos_out)

message("\nDone. Outputs in ", outdir, " and ", figdir)
sessionInfo()
