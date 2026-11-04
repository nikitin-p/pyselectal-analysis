#!/usr/bin/env Rscript
# NAD branch step 7: NAD captureSeq S-versus-N signal at the chrM ATATAAGTA starts.
#
# Refits the step-5 DESeq2 model per strain (design ~ library, N as reference,
# size factors by median of ratios over unmasked positions that pass the
# step-5 pre-filter), with one change: the -3..+3 windows around the final A of
# every ATATAAGTA on chrM, and the chrM processing-site controls, are added to
# the tested set even when they fail the pre-filter (zero rows if absent from
# the matrix). The size-factor rows are exactly those of step 5, so the
# normalisation is unchanged. padj is computed without independent filtering so
# that every added row gets one.
# Spec: docs/specs/2026-11-04-nad-step7-chrM-promoter-sites.md
#
# Outputs (--outdir): motif_sites.tsv, controls.tsv, sanity.tsv
# --dry-run stops before DESeq2 (used for the Mac test, where DESeq2 is absent).

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Biostrings)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--matrix",      type = "character"),
  make_option("--manifest",    type = "character"),
  make_option("--mask",        type = "character", help = "ncrna_mask.bed from step 5"),
  make_option("--fasta",       type = "character"),
  make_option("--step5-dir",   type = "character", default = NA_character_,
              help = "step-5 results dir, for the log2FC comparison"),
  make_option("--step5-log",   type = "character", default = NA_character_,
              help = "step-5 log, for the size-factor comparison"),
  make_option("--outdir",      type = "character"),
  make_option("--motif",       type = "character", default = "ATATAAGTA"),
  make_option("--window",      type = "integer", default = 3L),
  make_option("--core-sites",  type = "character", default = "73704:+,6548:+",
              help = "step-1 promoter starts with the AAGTA core but not the full motif, pos:strand"),
  make_option("--processing",  type = "character", default = "16471:+,13987:+,25343:+",
              help = "chrM processing-site controls, pos:strand, comma-separated"),
  make_option("--min-reads",   type = "integer", default = 10L),
  make_option("--min-samples", type = "integer", default = 3L),
  make_option("--padj",        type = "double",  default = 0.05),
  make_option("--lfc",         type = "double",  default = 1),
  make_option("--dry-run",     action = "store_true", default = FALSE)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...))
CHRM <- "chrM"

# ── genome, motif scan ───────────────────────────────────────────────────────
genome <- readDNAStringSet(opt$fasta)
names(genome) <- sub("\\s.*", "", names(genome))
stopifnot(CHRM %in% names(genome))
chrm <- genome[[CHRM]]
L <- nchar(opt$motif)
hp <- matchPattern(DNAString(opt$motif), chrm)
hm <- matchPattern(reverseComplement(DNAString(opt$motif)), chrm)
# + strand: the final A is the last base of the match; - strand: the final A of
# the motif on the RNA strand is the first base of the reverse-complement match.
sites <- rbind(data.table(site = end(hp),   strand = "+"),
               data.table(site = start(hm), strand = "-"))
sites[, class := "motif"]
msg("%s: %d sites (%d +, %d -)", opt$motif, nrow(sites), length(hp), length(hm))
stopifnot(nrow(sites) > 0)
parse_sites <- function(x) rbindlist(lapply(strsplit(strsplit(x, ",")[[1]], ":"), function(v)
  data.table(site = as.integer(v[1]), strand = v[2])))
# Step-1 promoter starts that carry only the AAGTA core (73704 at COX2, 6548 at
# 15S rRNA): reported as their own class, since they were chosen by CAGE.
core <- parse_sites(opt$`core-sites`)[, class := "step1_core"]
core <- core[!sites, on = .(site, strand)]
sites <- rbind(sites, core)
sites[, site_id := sprintf("%s%s%d", CHRM, strand, site)]

offs <- seq(-opt$window, opt$window)
targets <- sites[, .(offset = offs), by = .(site_id, site, strand, class)]
targets[, pos := ifelse(strand == "+", site + offset, site - offset)]
proc <- parse_sites(opt$processing)
proc[, `:=`(site_id = sprintf("%s%s%d", CHRM, strand, site), offset = 0L, pos = site,
            class = "processing")]
targets <- rbind(targets, proc[, names(targets), with = FALSE])
targets[, chr := CHRM]

# Sequence context -1..+4 on the RNA strand, as in step 5, and the +1 base.
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
targets <- add_context(targets)
targets[, base_p1 := substr(ctx, 2, 2)]
stopifnot(all(targets[class != "processing" & offset == 0, base_p1] == "A"))
# On the RNA strand the motif ends ...GTA, so the context at +1 must start "TA".
# step1_core sites end in AAGTA too, so the same holds for them.
stopifnot(all(startsWith(targets[class != "processing" & offset == 0, ctx], "TA")))

# ── matrix, mask, augmentation ───────────────────────────────────────────────
man <- fread(opt$manifest, select = c("run", "strain", "library", "replicate"))
mat <- fread(cmd = sprintf("gzip -dcf '%s'", opt$matrix))
stopifnot(all(man$run %in% names(mat)))
msg("matrix: %d positions x %d runs", nrow(mat), nrow(man))

add <- targets[!mat, on = .(chr, pos, strand), .(chr, pos, strand, base_p1)]
add <- unique(add)
if (nrow(add)) {
  z <- as.data.table(matrix(0L, nrow(add), nrow(man), dimnames = list(NULL, man$run)))
  mat <- rbind(mat, cbind(add, z), fill = TRUE)
  msg("added %d target positions absent from the matrix as zero rows", nrow(add))
}

mask <- fread(opt$mask, header = FALSE,
              col.names = c("chr", "start", "end", "gene_id", "biotype", "strand"))
mask[, start := start + 1L]                     # BED -> 1-based closed
mat[, `:=`(start = pos, end = pos)]
setkey(mask, chr, strand, start, end)
ov <- foverlaps(mat[, .(chr, strand, pos, start, end)], mask,
                by.x = c("chr", "strand", "start", "end"), nomatch = NULL)
ov <- unique(ov, by = c("chr", "strand", "pos"))
mat[ov, on = .(chr, strand, pos), `:=`(biotype = i.biotype, gene_id = i.gene_id)]
mat[is.na(biotype), `:=`(biotype = "none", gene_id = NA_character_)]
mat[, c("start", "end") := NULL]
mat[, is_target := FALSE][targets, on = .(chr, pos, strand), is_target := TRUE]

# ── per strain ───────────────────────────────────────────────────────────────
s5_sf <- function(st) {                          # size factors printed by step 5
  if (is.na(opt$`step5-log`) || !file.exists(opt$`step5-log`)) return(NULL)
  ln <- grep(sprintf("%s size factors:", st), readLines(opt$`step5-log`), value = TRUE, fixed = TRUE)
  if (!length(ln)) return(NULL)
  kv <- regmatches(ln[length(ln)], gregexpr("[A-Za-z0-9_]+=[0-9.]+", ln[length(ln)]))[[1]]
  setNames(as.numeric(sub(".*=", "", kv)), sub("=.*", "", kv))
}

motif_out <- list(); ctrl_out <- list(); sanity <- list()
for (st in unique(man$strain)) {
  cs <- man[strain == st]
  runs <- cs$run
  cnt <- as.matrix(mat[, ..runs])
  pass <- rowSums(cnt >= opt$`min-reads`) >= opt$`min-samples`
  keep <- pass | mat$is_target
  pos <- mat[keep, .(chr, pos, strand, base_p1, biotype, gene_id)]
  pos[, prefilter := pass[keep]]
  cnt <- cnt[keep, , drop = FALSE]
  rownames(cnt) <- paste(pos$chr, pos$pos, pos$strand, sep = ":")
  sf_rows <- which(pos$prefilter & pos$biotype == "none")   # = step 5 controlGenes
  msg("%s: %d tested (%d pass the pre-filter, %d added), %d size-factor rows",
      st, nrow(pos), sum(pos$prefilter), sum(!pos$prefilter), length(sf_rows))
  if (opt$`dry-run`) next

  suppressPackageStartupMessages(library(DESeq2))
  cd <- data.frame(library = factor(cs$library, levels = c("N", "S")), row.names = runs)
  dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ library)
  dds <- estimateSizeFactors(dds, controlGenes = sf_rows)
  msg("%s size factors: %s", st, paste(sprintf("%s=%.3f", runs, sizeFactors(dds)), collapse = " "))
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("library", "S", "N"), alpha = opt$padj,
                 independentFiltering = FALSE)

  nc <- counts(dds, normalized = TRUE)
  lab <- paste0(cs$library, "_", cs$replicate)
  out <- copy(pos)
  out[, (lab) := as.data.table(nc)]
  out[, `:=`(meanS = rowMeans(nc[, cs$library == "S", drop = FALSE]),
             meanN = rowMeans(nc[, cs$library == "N", drop = FALSE]),
             baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
             lfcSE = res$lfcSE, pvalue = res$pvalue, padj = res$padj)]
  out[, called := !is.na(padj) & padj < opt$padj & log2FoldChange > opt$lfc]
  out[, strain := st]

  cols <- c("strain", "class", "site_id", "site", "offset", "chr", "pos", "strand",
            "ctx", "base_p1", "biotype", "prefilter", lab, "meanS", "meanN",
            "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "called")
  tg <- merge(targets[, !"base_p1"], out, by = c("chr", "pos", "strand"), sort = FALSE)
  motif_out[[st]] <- tg[class != "processing", ..cols]
  trna <- out[chr == CHRM & biotype == "tRNA" & prefilter]
  trna[, `:=`(class = "tRNA", site_id = sprintf("%s%s%d", chr, strand, pos), site = pos, offset = 0L)]
  trna <- add_context(trna)
  ctrl_out[[st]] <- rbind(tg[class == "processing", ..cols], trna[, ..cols])

  # Sanity: size factors against the step-5 log, log2FC on shared rows.
  ref <- s5_sf(st)
  sf <- sizeFactors(dds)
  sanity[[st]] <- data.table(strain = st, item = paste0("size_factor_", runs),
                             step7 = round(sf, 3),
                             step5 = if (is.null(ref)) NA_real_ else round(ref[runs], 3))
  f5 <- if (is.na(opt$`step5-dir`)) "" else file.path(opt$`step5-dir`, sprintf("%s_S_vs_N.tsv.gz", st))
  if (file.exists(f5)) {
    s5 <- fread(cmd = sprintf("gzip -dcf '%s'", f5), select = c("chr", "pos", "strand", "log2FoldChange"))
    j <- merge(out[, .(chr, pos, strand, l7 = log2FoldChange)], s5, by = c("chr", "pos", "strand"))
    j <- j[is.finite(l7) & is.finite(log2FoldChange)]
    sanity[[st]] <- rbind(sanity[[st]],
      data.table(strain = st, item = c("lfc_shared_rows", "lfc_pearson_r", "lfc_max_abs_diff"),
                 step7 = c(nrow(j), round(cor(j$l7, j$log2FoldChange), 6),
                           round(max(abs(j$l7 - j$log2FoldChange)), 4)),
                 step5 = NA_real_))
  }
}

if (opt$`dry-run`) {
  fwrite(targets, file.path(opt$outdir, "targets_dry_run.tsv"), sep = "\t")
  msg("dry run: targets written, DESeq2 skipped")
  quit(save = "no", status = 0)
}

# Per-strain replicate columns differ only by name pattern (S_rep1 ...), so the
# tables bind by name.
motif_out <- rbindlist(motif_out, use.names = TRUE, fill = TRUE)
ctrl_out  <- rbindlist(ctrl_out,  use.names = TRUE, fill = TRUE)
sanity    <- rbindlist(sanity)
setorder(motif_out, strain, strand, site, offset)
fwrite(motif_out, file.path(opt$outdir, "motif_sites.tsv"), sep = "\t")
fwrite(ctrl_out,  file.path(opt$outdir, "controls.tsv"),    sep = "\t")
fwrite(sanity,    file.path(opt$outdir, "sanity.tsv"),      sep = "\t")

# ── acceptance checks (spec) ─────────────────────────────────────────────────
print(sanity)
print(motif_out[offset == 0, .(strain, site_id, prefilter, meanS = round(meanS, 1),
                               meanN = round(meanN, 1), lfc = round(log2FoldChange, 2),
                               padj = signif(padj, 2), called)])
n_exp <- nrow(sites) * length(offs) * uniqueN(man$strain)
sfr <- sanity[startsWith(item, "size_factor_")]
chk <- list(
  c(sprintf("size factors equal step 5 to 3 decimals (%d compared)", sfr[!is.na(step5), .N]),
    isTRUE(sfr[!is.na(step5), .N] > 0 && sfr[!is.na(step5), all(abs(step7 - step5) < 1e-3)])),
  c(sprintf("log2FC r on shared rows > 0.999 (min r %s)",
            paste(sanity[item == "lfc_pearson_r", step7], collapse = ", ")),
    isTRUE(sanity[item == "lfc_pearson_r", .N] > 0 && all(sanity[item == "lfc_pearson_r", step7] > 0.999))),
  c(sprintf("motif_sites rows = sites x offsets x strains: %d vs %d", nrow(motif_out), n_exp),
    nrow(motif_out) == n_exp),
  c(sprintf("chrM tRNA controls not enriched, median log2FC: %s",
            paste(ctrl_out[class == "tRNA", sprintf("%s %.2f", strain[1], median(log2FoldChange, na.rm = TRUE)),
                           by = strain]$V1, collapse = "; ")),
    all(ctrl_out[class == "tRNA", median(log2FoldChange, na.rm = TRUE), by = strain]$V1 < 0)))
for (k in chk) msg("[%s] %s", if (k[2] == "TRUE") "PASS" else "FAIL", k[1])

sessionInfo()
