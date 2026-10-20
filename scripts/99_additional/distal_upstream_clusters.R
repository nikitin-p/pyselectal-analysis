#!/usr/bin/env Rscript
# Draw the TSS clusters that lie upstream of a gene but beyond the longest
# promoter window, which the localization therefore labels "Intergenic".
#
# Under promoter_source: clusters, a promoter is found by the growing windows of
# promoter_cluster_windows (scripts/05_tss_end_types/promoter_helpers.R), so a
# same-strand cluster further upstream than the longest window is not a
# promoter. Such a cluster is one of two kinds:
#   distal_only    a cluster of a gene that no window reaches, so the gene has
#                  no promoter at all;
#   distal_behind  a further cluster of a gene that already has a nearer
#                  promoter cluster.
# Both are counted out to --reach bases upstream (300, the distance of the first
# cluster definition). A cluster is assigned to the first class that applies in
# the order promoter, distal_only, distal_behind, and every remaining cluster is
# "other". The figure compares the classes on depth, 1Sg share, M share and the
# YR initiator at the dominant TSS, to show whether the distal clusters look
# like start sites, and draws an example locus of each kind.
#
# Outputs:
#   <outdir>/distal_upstream_clusters.tsv       one row per cluster, with its class
#   <figdir>/25_distal_upstream_clusters.pdf
#
# Usage (from project root):
#   Rscript scripts/99_additional/distal_upstream_clusters.R \
#       [--ctss   results/tss/typed_ctss.tsv.gz] \
#       [--matrix results/tss/tss_matrix.tsv]    \
#       [--genes  data/sacCer3.genes.bed]        \
#       [--genome data/sacCer3.fa]               \
#       [--reach  300]                           \
#       [--outdir results/tss] [--figdir results/figures] \
#       [--params config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(yaml)
  library(Rsamtools); library(Biostrings)
  library(ggplot2); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag)
  args[i + 1]
}
has_flag <- function(flag) flag %in% args

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "10_dominant_endtype", "dom_helpers.R"))
source(file.path(script_dir, "..", "05_tss_end_types", "promoter_helpers.R"))

params  <- yaml.load_file(get_arg("--params", "config/params.yaml"))
ctss_f  <- get_arg("--ctss",   "results/tss/typed_ctss.tsv.gz")
mat_f   <- get_arg("--matrix", "results/tss/tss_matrix.tsv")
genes_f <- get_arg("--genes",  params$annotation_bed)
fasta_f <- get_arg("--genome", params$genome_fasta)
reach   <- as.integer(get_arg("--reach", "300"))
outdir  <- get_arg("--outdir", "results/tss")
figdir  <- get_arg("--figdir", "results/figures")
force   <- has_flag("--force")

windows <- sort(as.integer(unlist(params$promoter_cluster_windows)))
x_max   <- max(windows)
out_tsv <- file.path(outdir, "distal_upstream_clusters.tsv")
out_pdf <- file.path(figdir, "25_distal_upstream_clusters.pdf")
if (!force && all(file.exists(c(out_tsv, out_pdf)))) {
  message("Outputs exist, skipping (use --force to rerun): ", out_pdf)
  quit(save = "no", status = 0)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

CLASS_LAB <- c(promoter      = sprintf("Promoter (within %d bp)", x_max),
               distal_only   = sprintf("Distal, gene has no promoter (%d-%d bp)", x_max + 1L, reach),
               distal_behind = sprintf("Distal, behind a promoter (%d-%d bp)", x_max + 1L, reach),
               other         = "Other clusters")
CLASS_COL <- c(promoter = "#d62828", distal_only = "#f4a261",
               distal_behind = "#457b9d", other = "grey65")

# ── classify the clusters ─────────────────────────────────────────────────────
ext   <- load_cluster_extents(ctss_f, mat_f, as.integer(params$a8_min_total),
                              as.integer(params$a8_max_dist))
cl    <- ext$clusters
genes <- load_gene_starts(genes_f)
pairs <- candidate_pairs(cl, genes, reach)
prom  <- define_promoters(pairs, windows)

distal <- pairs[gap > x_max]
distal[, kind := fifelse(gene %in% prom$gene, "distal_behind", "distal_only")]

cl[, class := fifelse(cluster_id %in% prom$cluster_id, "promoter",
              fifelse(cluster_id %in% distal[kind == "distal_only"]$cluster_id, "distal_only",
              fifelse(cluster_id %in% distal[kind == "distal_behind"]$cluster_id,
                      "distal_behind", "other")))]

# Per-cluster composition from the pooled per-type counts of the members.
typed <- fread(ctss_f, select = c("chr", "pos", "strand", "n_1Sg", "n_M", "n_total"))
typed <- typed[, lapply(.SD, sum), by = .(chr, pos, strand)]
comp  <- merge(ext$members[, .(chr, pos, strand, cluster_id)], typed,
               by = c("chr", "pos", "strand"))[
  , .(n_1Sg = sum(n_1Sg), n_M = sum(n_M)), by = cluster_id]
cl <- merge(cl, comp, by = "cluster_id")
cl[, `:=`(pct_1Sg = 100 * n_1Sg / reads, pct_M = 100 * n_M / reads)]

# Dinucleotide -1/+1 at the dominant TSS, read in transcript orientation.
fa <- FaFile(fasta_f)
din_gr <- GRanges(cl$chr, IRanges(ifelse(cl$strand == "+", cl$tss_pos - 1L, cl$tss_pos),
                                  width = 2L), strand = cl$strand)
cl[, dinuc := as.character(getSeq(fa, din_gr))]
cl[, YR := substr(dinuc, 1, 1) %in% c("C", "T") & substr(dinuc, 2, 2) %in% c("A", "G")]

fwrite(cl[, .(cluster_id, chr, start, end, strand, tss_pos, reads, class,
              pct_1Sg, pct_M, dinuc, YR)], out_tsv, sep = "\t")
message("Written: ", out_tsv)

summ <- cl[, .(clusters = .N, reads = sum(reads)), by = class]
summ[, pct_reads := 100 * reads / sum(reads)]
summ[, class := factor(class, levels = names(CLASS_LAB))]
setorder(summ, class)
genes_only <- uniqueN(distal[kind == "distal_only"]$gene)
genes_behind <- uniqueN(distal[kind == "distal_behind"]$gene)
print(summ)
message("Genes with no promoter but a cluster within ", reach, " bp: ", genes_only,
        "; genes with a promoter and a further distal cluster: ", genes_behind)

# ── panel A: what the two cases look like ─────────────────────────────────────
sch_case <- function(title, clusters, ctss) {
  ggplot() +
    annotate("rect", xmin = -x_max, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.6) +
    annotate("text", x = -x_max / 2, y = 1.08, size = 2.8,
             label = sprintf("%d bp windows", x_max)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_segment(data = ctss, aes(x = x, xend = x, y = 0, yend = h), linewidth = 0.6) +
    geom_segment(data = clusters, aes(x = from, xend = to, y = -0.15, yend = -0.15,
                                      colour = class), linewidth = 2.2) +
    geom_text(data = clusters, aes(x = (from + to) / 2, y = -0.32, label = label,
                                   colour = class), size = 2.8) +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 0.6, linewidth = 0.6) +
    annotate("segment", x = 0, xend = 60, y = 0.6, yend = 0.6, linewidth = 0.6,
             arrow = arrow(length = unit(0.12, "cm"))) +
    annotate("text", x = 3, y = 0.7, label = "gene", hjust = 0, size = 3) +
    scale_colour_manual(values = CLASS_COL, guide = "none") +
    coord_cartesian(xlim = c(-reach, 70), ylim = c(-0.45, 1.15)) +
    labs(x = "Position relative to the gene start (bp)", y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
}
p_a1 <- sch_case(sprintf("A1. Distal, gene has no promoter: %d genes", genes_only),
                 data.table(from = -150, to = -110, class = "distal_only",
                            label = "labelled Intergenic"),
                 data.table(x = c(-148, -140, -133, -126, -118, -112),
                            h = c(0.2, 0.5, 0.9, 0.4, 0.3, 0.15)))
p_a2 <- sch_case(sprintf("A2. Distal, behind a promoter: %d genes", genes_behind),
                 data.table(from = c(-200, -35), to = c(-160, -12),
                            class = c("distal_behind", "promoter"),
                            label = c("labelled Intergenic", "Promoter")),
                 data.table(x = c(-195, -184, -176, -168, -30, -24, -18, -14),
                            h = c(0.3, 0.8, 0.5, 0.2, 0.4, 0.9, 0.6, 0.2)))

# ── panel B: how much signal is affected ──────────────────────────────────────
p_b <- ggplot(summ, aes(pct_reads, class, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f %% of reads, %s clusters", pct_reads,
                                format(clusters, big.mark = ","))),
            hjust = -0.05, size = 3) +
  scale_y_discrete(labels = CLASS_LAB, limits = rev) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.85)), breaks = seq(0, 75, 25)) +
  labs(x = "Clustered reads (%)", y = NULL,
       title = "B. Share of the clustered signal in each class") +
  theme_minimal(base_size = 10)

# ── panel C: how far upstream the distal clusters lie ────────────────────────
p_c <- ggplot(distal, aes(gap, fill = kind)) +
  geom_histogram(binwidth = 10, boundary = x_max, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = CLASS_COL, labels = CLASS_LAB, name = NULL) +
  labs(x = "Gap from the gene start to the cluster (bp)", y = "Gene-cluster pairs",
       title = "C. Distance of the distal clusters from the gene start") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", legend.direction = "vertical")

# ── panel D: do they look like start sites? ───────────────────────────────────
props <- rbind(
  cl[, .(class, value = log10(reads), metric = "Cluster depth (log10 reads)")],
  cl[, .(class, value = pct_1Sg, metric = "1Sg share of reads (%)")],
  cl[, .(class, value = pct_M, metric = "M share of reads (%)")])
props[, class := factor(class, levels = names(CLASS_LAB))]
p_d1 <- ggplot(props, aes(class, value, fill = class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, width = 0.65) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  scale_x_discrete(labels = c(promoter = "Promoter", distal_only = "Distal,\nno promoter",
                              distal_behind = "Distal,\nbehind", other = "Other")) +
  labs(x = NULL, y = NULL, title = "D. Properties of the clusters in each class") +
  theme_minimal(base_size = 10)
yr <- cl[, .(pct = 100 * mean(YR), n = .N), by = class]
yr[, class := factor(class, levels = names(CLASS_LAB))]
p_d2 <- ggplot(yr, aes(class, pct, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f %%", pct)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  scale_x_discrete(labels = c(promoter = "Promoter", distal_only = "Distal,\nno promoter",
                              distal_behind = "Distal,\nbehind", other = "Other")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Clusters (%)", title = "YR initiator at the dominant TSS") +
  theme_minimal(base_size = 10)

# ── panel E: one real locus of each kind ──────────────────────────────────────
real_case <- function(kind_sel, title) {
  # The deepest example whose distal cluster carries the class itself (not the
  # promoter of another gene), lies wholly 41..reach bp upstream and is at most
  # 80 bp wide, so the panel shows one clear case.
  cand <- merge(distal[kind == kind_sel], cl[, .(cluster_id, class)], by = "cluster_id")
  cand <- cand[class == kind_sel & gap > x_max &
               cluster_end - cluster_start < 80L &
               fifelse(strand == "+", gene_start - cluster_start, cluster_end - gene_start) < reach - 10L]
  setorder(cand, -reads)
  g <- cand$gene[1]
  gp <- pairs[gene == g]; gs <- gp$gene_start[1]; st <- gp$strand[1]
  rel <- function(p) if (st == "+") p - gs else gs - p
  cls <- merge(gp[, .(cluster_id)], cl[, .(cluster_id, start, end, class)], by = "cluster_id")
  cls[, `:=`(a = rel(start), b = rel(end))][, `:=`(from = pmin(a, b), to = pmax(a, b))]
  m <- ext$members[cluster_id %in% cls$cluster_id]
  m[, x := rel(pos)]
  m[, h := 0.9 * sqrt(n_total / max(n_total))]
  ggplot() +
    annotate("rect", xmin = -x_max, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_segment(data = m, aes(x = x, xend = x, y = 0, yend = h), linewidth = 0.5) +
    geom_segment(data = cls, aes(x = from, xend = to, y = -0.15, yend = -0.15, colour = class),
                 linewidth = 2.2) +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 1, linewidth = 0.6) +
    annotate("segment", x = 0, xend = 60, y = 1, yend = 1, linewidth = 0.6,
             arrow = arrow(length = unit(0.12, "cm"))) +
    annotate("text", x = 3, y = 1.08, label = g, hjust = 0, size = 3) +
    scale_colour_manual(values = CLASS_COL, labels = CLASS_LAB, name = NULL) +
    coord_cartesian(xlim = c(-reach, 70), ylim = c(-0.3, 1.15)) +
    labs(x = "Position relative to the gene start (bp)", y = NULL,
         title = sprintf("%s %s (%s strand), distal cluster %s reads", title, g, st,
                         format(cand$reads[1], big.mark = ","))) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
}
p_e1 <- real_case("distal_only",   "E1. No promoter:")
p_e2 <- real_case("distal_behind", "E2. Behind a promoter:")

# ── panel F: the options for discussion ──────────────────────────────────────
# The numbers are computed here so the text follows the params file.
targets   <- as.integer(unlist(params$promoter_target_genes))
tfrac     <- if (is.null(params$promoter_target_fraction)) 1 else
  as.numeric(params$promoter_target_fraction)
gaps      <- sort(pairs[, .(gap = min(gap)), by = gene]$gap)
x_full    <- if (length(targets) && max(targets) <= length(gaps)) gaps[max(targets)] else NA
pct_dist  <- sum(summ[class %in% c("distal_only", "distal_behind")]$pct_reads)
target_txt <- if (length(targets))
  sprintf(" (at least %g %% of the %s expressed genes)", 100 * tfrac,
          paste(format(range(targets), big.mark = ","), collapse = "-")) else ""
options_txt <- c(
  sprintf("1. Keep the %d bp windows and accept the mislabel. %s genes get a promoter%s, but %.1f %% of the clustered reads, which come from real start sites, are labelled Intergenic.",
          x_max, format(uniqueN(prom$gene), big.mark = ","), target_txt, pct_dist),
  sprintf("2. Add a \"Distal upstream\" localization category for clusters %d-%d bp upstream of a same-strand gene. The promoters keep their strict definition, and this signal is no longer counted as Intergenic.",
          x_max + 1L, reach),
  sprintf("3. Let the windows grow further, up to %s bp, which reaches all %s genes. This fixes part of the no-promoter case (A1, E1) but not the behind-a-promoter case (A2, E2), because the growing rule stops at the first cluster it reaches.",
          ifelse(is.na(x_full), paste0(">", reach), x_full),
          format(max(targets), big.mark = ",")),
  "4. Pick the strongest cluster instead of the nearest within a longer reach. This handles the behind-a-promoter case (E2), but it changes the definition, under which the nearest cluster decides.")
opt <- data.table(y = rev(seq_along(options_txt)),
                  label = vapply(options_txt, function(t) paste(strwrap(t, 150), collapse = "\n"),
                                 character(1)))
p_f <- ggplot(opt, aes(0, y, label = label)) +
  geom_text(hjust = 0, vjust = 0.5, size = 3.4, lineheight = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, length(options_txt) + 0.5), clip = "off") +
  labs(title = "F. Options to discuss") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 13, margin = margin(b = 6)))

pdf(out_pdf, width = 13, height = 17.5)
d_row <- (p_d1 | p_d2) + plot_layout(widths = c(3, 1))
print((p_a1 | p_a2) / (p_b | p_c) / d_row / (p_e1 | p_e2) / p_f +
        plot_layout(heights = c(0.8, 1, 1, 0.9, 0.55)) +
        plot_annotation(
          title = sprintf("TSS clusters upstream of a gene but beyond the %d bp promoter windows are labelled Intergenic", x_max),
          subtitle = sprintf("Windows of %s bp define the promoter; clusters %d-%d bp upstream of a same-strand gene start are counted here",
                             paste(windows, collapse = ", "), x_max + 1L, reach)))
invisible(dev.off())
message("Written: ", out_pdf)

sessionInfo()
