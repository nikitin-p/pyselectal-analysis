#!/usr/bin/env Rscript
# NAD branch step 5: S-versus-N enrichment of NAD captureSeq 5' ends (DESeq2).
#
# One DESeq2 fit per strain, design ~ library with N as reference, on raw read
# counts per single-nucleotide 5' end. Positions inside tRNA/rRNA/snoRNA/snRNA/
# ncRNA genes on the same strand are tested and reported with their biotype,
# but are excluded from size-factor estimation (controlGenes), so that a few
# very abundant ncRNA fragments cannot drive normalisation.
# Spec: docs/specs/2026-11-03-nad-step5-deseq2.md
#
# Outputs (in --outdir unless noted):
#   ncrna_mask.bed, <strain>_S_vs_N.tsv.gz, <strain>_enriched.tsv,
#   summary.tsv, context_enriched.tsv,
#   <figdir>/05_<strain>_ma.pdf, <figdir>/05_first_base.pdf

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(Biostrings)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--matrix",      type = "character"),
  make_option("--manifest",    type = "character"),
  make_option("--gtf",         type = "character", help = "Ensembl R64-1-1 GTF (.gz ok)"),
  make_option("--fasta",       type = "character"),
  make_option("--outdir",      type = "character"),
  make_option("--figdir",      type = "character"),
  make_option("--min-reads",   type = "integer", default = 10L),
  make_option("--min-samples", type = "integer", default = 3L),
  make_option("--padj",        type = "double",  default = 0.05),
  make_option("--lfc",         type = "double",  default = 1)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$figdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))

MASK_BIOTYPES <- c("tRNA", "rRNA", "snoRNA", "snRNA", "ncRNA")
BASES <- c("A", "C", "G", "T")

# ── non-coding RNA mask from the Ensembl GTF ─────────────────────────────────
gtf <- fread(cmd = sprintf("gzip -dcf '%s' | grep -v '^#'", opt$gtf),
             sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "src", "feature", "start", "end",
                           "score", "strand", "frame", "attr"))
gtf <- gtf[feature == "gene"]
gtf[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
gtf[, biotype := sub('.*gene_biotype "([^"]+)".*', "\\1", attr)]
mask <- gtf[biotype %in% MASK_BIOTYPES, .(chr, start, end, strand, gene_id, biotype)]
# R64-1-1 is the sacCer3 assembly: only the chromosome names differ.
mask[, chr := ifelse(chr == "Mito", "chrM", paste0("chr", chr))]
fwrite(mask[, .(chr, start = start - 1L, end, gene_id, biotype, strand)],
       file.path(opt$outdir, "ncrna_mask.bed"), sep = "\t", col.names = FALSE)
msg("mask: %d genes (%s)", nrow(mask),
    paste(sprintf("%s %d", names(table(mask$biotype)), table(mask$biotype)), collapse = ", "))

# ── inputs ───────────────────────────────────────────────────────────────────
man <- fread(opt$manifest)[, .(run, strain, library, replicate)]
mat <- fread(opt$matrix)
stopifnot(all(man$run %in% names(mat)))
msg("matrix: %d positions x %d runs", nrow(mat), nrow(man))

# Strand-aware overlap with the mask (same strand only); a position inside
# several masked genes takes the first one.
mat[, `:=`(start = pos, end = pos)]
setkey(mask, chr, strand, start, end)
ov <- foverlaps(mat[, .(chr, strand, pos, start, end)], mask,
                by.x = c("chr", "strand", "start", "end"), nomatch = NULL)
ov <- unique(ov, by = c("chr", "strand", "pos"))
mat[ov, on = .(chr, strand, pos), `:=`(biotype = i.biotype, gene_id = i.gene_id)]
mat[is.na(biotype), `:=`(biotype = "none", gene_id = NA_character_)]
mat[, c("start", "end") := NULL]

# ── sequence context -1..+4 on the RNA strand ────────────────────────────────
genome <- readDNAStringSet(opt$fasta)
names(genome) <- sub("\\s.*", "", names(genome))
add_context <- function(dt) {
  dt[, ctx := NA_character_]
  for (ch in unique(dt$chr)) {
    if (!ch %in% names(genome)) next
    s <- as.character(genome[[ch]]); n <- nchar(s)
    i <- which(dt$chr == ch)
    plus <- dt$strand[i] == "+"
    from <- ifelse(plus, dt$pos[i] - 1L, dt$pos[i] - 4L)
    to   <- from + 5L
    ok <- from >= 1L & to <= n
    x <- rep(NA_character_, length(i))
    x[ok] <- toupper(substring(s, from[ok], to[ok]))
    m <- ok & !plus
    if (any(m)) x[m] <- as.character(reverseComplement(DNAStringSet(x[m])))
    dt[i, ctx := x]
  }
  dt
}

# ── one fit per strain ───────────────────────────────────────────────────────
is_yaag <- function(ctx) grepl("^[CT]AAG", ctx)
base_freq <- function(ctx) {
  ctx <- ctx[!is.na(ctx)]
  rbindlist(lapply(1:5, function(k) {
    b <- substr(ctx, k, k)
    data.table(offset = c(-1L, 1:4)[k], base = BASES,
               frac = vapply(BASES, function(x) mean(b == x), numeric(1)))
  }))
}

summ <- list(); ctxs <- list(); fb <- list()
for (st in unique(man$strain)) {
  cs <- man[strain == st]
  runs <- cs$run
  cnt <- as.matrix(mat[, ..runs])
  keep <- rowSums(cnt >= opt$`min-reads`) >= opt$`min-samples`
  pos <- mat[keep, .(chr, pos, strand, base_p1, biotype, gene_id)]
  cnt <- cnt[keep, , drop = FALSE]
  rownames(cnt) <- paste(pos$chr, pos$pos, pos$strand, sep = ":")
  msg("%s: %d positions pass >= %d reads in >= %d of %d runs (%d masked)",
      st, nrow(pos), opt$`min-reads`, opt$`min-samples`, length(runs), sum(pos$biotype != "none"))

  cd <- data.frame(library = factor(cs$library, levels = c("N", "S")), row.names = runs)
  dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ library)
  dds <- estimateSizeFactors(dds, controlGenes = which(pos$biotype == "none"))
  msg("%s size factors: %s", st,
      paste(sprintf("%s=%.3f", runs, sizeFactors(dds)), collapse = " "))
  dds <- DESeq(dds, quiet = TRUE)   # keeps the size factors set above
  res <- results(dds, contrast = c("library", "S", "N"), alpha = opt$padj)

  nc <- counts(dds, normalized = TRUE)
  out <- add_context(copy(pos))
  out[, `:=`(baseMean = res$baseMean,
             meanS = rowMeans(nc[, cs$library == "S", drop = FALSE]),
             meanN = rowMeans(nc[, cs$library == "N", drop = FALSE]),
             log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
             pvalue = res$pvalue, padj = res$padj)]
  setcolorder(out, c("chr", "pos", "strand", "base_p1", "ctx", "biotype", "gene_id"))
  n_bad <- out[!is.na(ctx) & substr(ctx, 2, 2) != base_p1, .N]
  if (n_bad > 0) stop(sprintf("%s: base_p1 disagrees with the FASTA at %d positions", st, n_bad))
  out[, enriched := !is.na(padj) & padj < opt$padj & log2FoldChange > opt$lfc]

  fwrite(out[, !"enriched"], file.path(opt$outdir, sprintf("%s_S_vs_N.tsv.gz", st)), sep = "\t")
  fwrite(out[enriched == TRUE][order(padj), !"enriched"],
         file.path(opt$outdir, sprintf("%s_enriched.tsv", st)), sep = "\t")

  out[, class := ifelse(biotype == "none", "unmasked", "masked")]
  summ[[st]] <- out[, .(tested = .N, enriched = sum(enriched),
                        fracA_tested   = mean(base_p1 == "A"),
                        fracA_enriched = mean(base_p1[enriched] == "A"),
                        fracYAAG_tested_A   = mean(is_yaag(ctx[base_p1 == "A" & !is.na(ctx)])),
                        fracYAAG_enriched_A = mean(is_yaag(ctx[enriched & base_p1 == "A" & !is.na(ctx)]))),
                    by = class][, strain := st]

  u <- out[class == "unmasked"]
  ctxs[[st]] <- rbind(base_freq(u$ctx)[, set := "tested"],
                      base_freq(u[enriched == TRUE]$ctx)[, set := "enriched"])[, strain := st]
  fb[[st]] <- ctxs[[st]][offset == 1L]

  p <- ggplot(out[!is.na(log2FoldChange)],
              aes(log10(baseMean), log2FoldChange,
                  colour = ifelse(enriched, "enriched", ifelse(class == "masked", "ncRNA", "tested")))) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_hline(yintercept = c(0, opt$lfc), linetype = c("solid", "dashed"), linewidth = 0.3) +
    scale_colour_manual(NULL, values = c(enriched = "#e63946", ncRNA = "#457b9d", tested = "grey60")) +
    labs(title = sprintf("%s: S vs N", st), x = "log10 mean normalised reads",
         y = "log2 fold change S / N") +
    theme_bw(base_size = 10) + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  ggsave(file.path(opt$figdir, sprintf("05_%s_ma.pdf", st)), p, width = 6, height = 4.5)
}

summ <- rbindlist(summ); setcolorder(summ, "strain")
ctxs <- rbindlist(ctxs); setcolorder(ctxs, c("strain", "set", "offset", "base"))
fwrite(summ, file.path(opt$outdir, "summary.tsv"), sep = "\t")
fwrite(ctxs, file.path(opt$outdir, "context_enriched.tsv"), sep = "\t")

p <- ggplot(rbindlist(fb), aes(set, frac, fill = base)) +
  geom_col() + facet_wrap(~ strain) +
  scale_fill_manual(values = c(A = "#e63946", C = "#457b9d", G = "#f4a261", T = "#2a9d8f")) +
  labs(x = NULL, y = "share of positions", title = "+1 base, unmasked positions") +
  theme_bw(base_size = 10)
ggsave(file.path(opt$figdir, "05_first_base.pdf"), p, width = 5, height = 3.5)

# ── acceptance checks (spec): triple mutant, unmasked positions ──────────────
print(summ)
mut <- setdiff(unique(man$strain), "WT")
u <- summ[class == "unmasked"]
for (st in mut) {
  r <- u[strain == st]; w <- u[strain == "WT"]
  checks <- c(
    sprintf("A share at +1 enriched > tested: %.3f vs %.3f", r$fracA_enriched, r$fracA_tested),
    sprintf("YAAG among A-starts enriched > tested: %.3f vs %.3f", r$fracYAAG_enriched_A, r$fracYAAG_tested_A),
    sprintf("enriched unmasked positions %s > WT: %d vs %d", st, r$enriched, w$enriched))
  pass <- c(isTRUE(r$fracA_enriched > r$fracA_tested),
            isTRUE(r$fracYAAG_enriched_A > r$fracYAAG_tested_A),
            isTRUE(r$enriched > w$enriched))
  for (k in seq_along(checks)) msg("[%s] %s", if (pass[k]) "PASS" else "FAIL", checks[k])
}

sessionInfo()
