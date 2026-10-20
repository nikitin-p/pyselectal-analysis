#!/usr/bin/env Rscript
# Stage 8 step 13: compare window settings of the cluster-defined promoter and
# draw the definition.
#
# The promoter definition of promoter_helpers.R is applied under a grid of
# window settings: a single fixed window of each length in --grid, and the
# growing rule over every grid length from --grow_from up to that length. A
# fixed window and a growing series with the same longest window reach exactly
# the same genes, and differ only in that the growing rule drops a further
# cluster once a nearer one has decided the gene, so the comparison isolates
# that effect. Each setting is summarised by the genes it gives a promoter, the
# clustered reads it puts in promoter clusters, the clusters per gene, the
# length of a gene's promoter, and the share of promoter bases that fall inside
# the body of another gene, which is the cost of reaching further.
#
# Two further checks frame the choice of windows. The window length at which
# the genes with a promoter first reach each target count is reported, the
# targets being an external estimate of the expressed genes
# (promoter_target_genes, scaled by promoter_target_fraction). And the genes are
# split into those reached by the shortest window, those gained by the longer
# ones, and those without a promoter, and the window of each gene under the
# longest setting is checked for other genes' bodies, opposite-strand TSS
# clusters, a head-to-head partner cluster and enhancers.
#
# Outputs:
#   <outdir>/promoter_window_sweep.tsv        one row per setting
#   <outdir>/promoter_window_targets.tsv      window length reaching each target
#   <outdir>/promoter_window_features.tsv     window features by gene group
#   <figdir>/24_promoter_window_sweep.pdf     gap distribution and the sweep
#   <figdir>/24_promoter_definition.pdf       schematic and three real loci,
#                                             drawn under promoter_cluster_windows
#   <figdir>/24_promoter_overhang.pdf         how far promoter clusters run past
#                                             the gene start, with real loci
#
# Usage (from project root):
#   Rscript scripts/05_tss_end_types/13_promoter_window_sweep.R \
#       [--ctss      results/tss/typed_ctss.tsv.gz] \
#       [--matrix    results/tss/tss_matrix.tsv]    \
#       [--genes     data/sacCer3.genes.bed]        \
#       [--outdir    results/tss]                   \
#       [--figdir    results/figures]               \
#       [--grid      10,15,20,25,30,40,50,75,100,150,200,300] \
#       [--grow_from 15]                            \
#       [--enhancers results/localization/enhancer_regions.bed] \
#       [--divergent 300]                           \
#       [--params    config/params.yaml] [--force]

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(yaml)
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
source(file.path(script_dir, "promoter_helpers.R"))

params    <- yaml.load_file(get_arg("--params", "config/params.yaml"))
ctss_f    <- get_arg("--ctss",   "results/tss/typed_ctss.tsv.gz")
mat_f     <- get_arg("--matrix", "results/tss/tss_matrix.tsv")
genes_f   <- get_arg("--genes",  params$annotation_bed)
outdir    <- get_arg("--outdir", "results/tss")
figdir    <- get_arg("--figdir", "results/figures")
grid      <- sort(as.integer(strsplit(get_arg("--grid", "10,15,20,25,30,40,50,75,100,150,200,300"), ",")[[1]]))
grow_from <- as.integer(get_arg("--grow_from", "15"))
divergent <- as.integer(get_arg("--divergent", "300"))
force     <- has_flag("--force")

min_total <- as.integer(params$a8_min_total)
max_dist  <- as.integer(params$a8_max_dist)
chosen    <- sort(as.integer(unlist(params$promoter_cluster_windows)))
if (length(chosen) == 0L || anyNA(chosen))
  stop("promoter_cluster_windows must be a list of positive integers in the params file")

out_tsv   <- file.path(outdir, "promoter_window_sweep.tsv")
out_tgt   <- file.path(outdir, "promoter_window_targets.tsv")
out_feat  <- file.path(outdir, "promoter_window_features.tsv")
out_sweep <- file.path(figdir, "24_promoter_window_sweep.pdf")
out_def   <- file.path(figdir, "24_promoter_definition.pdf")
out_over  <- file.path(figdir, "24_promoter_overhang.pdf")
if (!force && all(file.exists(c(out_tsv, out_tgt, out_feat, out_sweep, out_def, out_over)))) {
  message("Outputs exist, skipping (use --force to rerun): ", out_tsv)
  quit(save = "no", status = 0)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

PROMOTER_COL <- "#d62828"
SCHEME_COLS  <- c(fixed = "grey45", growing = "#1d3557")

ext   <- load_cluster_extents(ctss_f, mat_f, min_total, max_dist)
cl    <- ext$clusters
genes <- load_gene_starts(genes_f)
pairs <- candidate_pairs(cl, genes, max(c(grid, chosen)))
message("Clusters: ", nrow(cl), "; same-strand gene-cluster pairs within ",
        max(c(grid, chosen)), " bp: ", nrow(pairs))

# ── sweep ─────────────────────────────────────────────────────────────────────
body_gr <- GRanges(genes$chr, IRanges(genes$start + 1L, genes$end), gene = genes$gene)

# Share of a setting's promoter bases that lie inside the body of a gene other
# than the promoter's own, on either strand. A gene's promoter is the union of
# its spans, which is one interval because every span reaches the gene start.
other_gene_share <- function(p) {
  g <- p[, .(chr = chr[1], start = min(start), end = max(end)), by = gene]
  g_gr <- GRanges(g$chr, IRanges(g$start, g$end))
  h <- findOverlaps(g_gr, body_gr, ignore.strand = TRUE)
  h <- h[g$gene[queryHits(h)] != body_gr$gene[subjectHits(h)]]
  if (length(h) == 0L) return(c(bp = sum(width(g_gr)), inside = 0))
  inter <- pintersect(g_gr[queryHits(h)], body_gr[subjectHits(h)], ignore.strand = TRUE)
  inside <- sum(sum(width(reduce(split(inter, queryHits(h))))))
  c(bp = sum(width(g_gr)), inside = inside)
}

summarise_setting <- function(p, scheme, max_w) {
  ov <- other_gene_share(p)
  span <- p[, .(len = max(end) - min(start) + 1L), by = gene]
  data.table(scheme = scheme, max_window = max_w,
             genes = uniqueN(p$gene),
             pct_genes = 100 * uniqueN(p$gene) / nrow(genes),
             clusters = uniqueN(p$cluster_id),
             pct_reads = 100 * sum(cl$reads[cl$cluster_id %in% p$cluster_id]) / sum(cl$reads),
             clusters_per_gene = nrow(p) / uniqueN(p$gene),
             median_promoter_bp = median(span$len),
             pct_bp_in_other_gene = 100 * ov[["inside"]] / ov[["bp"]])
}

sweep <- rbindlist(c(
  lapply(grid, function(w) summarise_setting(define_promoters(pairs, w), "fixed", w)),
  lapply(grid[grid >= grow_from], function(w)
    summarise_setting(define_promoters(pairs, grid[grid >= grow_from & grid <= w]), "growing", w))))
fwrite(sweep, out_tsv, sep = "\t")
message("Written: ", out_tsv)
print(sweep, digits = 3)

# ── target gene counts ────────────────────────────────────────────────────────
# The genes a setting reaches depend only on its longest window, so the n-th
# smallest nearest-cluster gap is the shortest window giving n genes a promoter.
targets <- as.integer(unlist(params$promoter_target_genes))
tfrac   <- if (is.null(params$promoter_target_fraction)) 1 else
  as.numeric(params$promoter_target_fraction)
reach   <- sort(pairs[, .(gap = min(gap)), by = gene]$gap)
tgt <- unique(data.table(target = rep(targets, each = 2L),
                         fraction = rep(c(tfrac, 1), length(targets))))
tgt[, genes_needed := as.integer(ceiling(target * fraction))]
# NA means the target is not reached within the longest grid window.
tgt[, x_bp := ifelse(genes_needed <= length(reach), reach[pmin(genes_needed, length(reach))],
                     NA_integer_)]
tgt[, genes_at_chosen := sum(reach <= max(chosen))]
fwrite(tgt, out_tgt, sep = "\t")
message("Written: ", out_tgt)
print(tgt)

# ── window features by gene group ─────────────────────────────────────────────
chosen_p <- define_promoters(pairs, chosen)
x_feat   <- max(chosen)
grp_lab  <- c(first  = sprintf("reached by %d bp", min(chosen)),
              gained = sprintf("gained at %d-%d bp", chosen[min(2L, length(chosen))], x_feat),
              none   = sprintf("no promoter at %d bp", x_feat))
dec <- chosen_p[, .(window = window[1]), by = gene]
genes[, group := fifelse(gene %in% dec[window == min(chosen)]$gene, "first",
                 fifelse(gene %in% dec$gene, "gained", "none"))]

root <- if (is.null(params$results_root)) "results" else params$results_root
enh_f <- c(get_arg("--enhancers"), file.path(root, "localization", "enhancer_regions.bed"),
           file.path(root, "tracks", "enhancers.bed"))
enh_f <- enh_f[file.exists(enh_f)][1]
enh_gr <- if (is.na(enh_f)) GRanges() else {
  e <- fread(enh_f, header = FALSE, select = 1:3)
  GRanges(e$V1, IRanges(e$V2 + 1L, e$V3))
}
message(if (is.na(enh_f)) "No enhancer BED found; the enhancer row is skipped" else
        paste0("Enhancers: ", length(enh_gr), " from ", enh_f))

upstream_gr <- function(g, x, strand = g$strand)
  GRanges(g$chr, IRanges(pmax(1L, ifelse(g$strand == "+", g$gene_start - x, g$gene_start + 1L)),
                         ifelse(g$strand == "+", g$gene_start - 1L, g$gene_start + x)),
          strand = strand)
flip  <- function(s) ifelse(s == "+", "-", "+")
cl_gr <- GRanges(cl$chr, IRanges(cl$start, cl$end), strand = cl$strand)
body_str <- GRanges(genes$chr, IRanges(genes$start + 1L, genes$end), strand = genes$strand)

# Indices of the genes whose upstream region overlaps something, other than
# their own body where the subject is a gene.
hit_genes <- function(q, subj, ignore.strand = FALSE, self = FALSE) {
  h <- findOverlaps(q, subj, ignore.strand = ignore.strand)
  if (self) h <- h[queryHits(h) != subjectHits(h)]
  unique(queryHits(h))
}
w_same <- upstream_gr(genes, x_feat)
w_anti <- upstream_gr(genes, x_feat, flip(genes$strand))
body_anti <- GRanges(seqnames(body_str), ranges(body_str), strand = flip(genes$strand))
feat_hits <- list(
  "Same-strand gene body in the window"       = hit_genes(w_same, body_str, self = TRUE),
  "Opposite-strand gene body in the window"   = hit_genes(w_same, body_anti, self = TRUE),
  "Opposite-strand TSS cluster in the window" = hit_genes(w_anti, cl_gr),
  "Opposite-strand cluster within the head-to-head range" =
    hit_genes(upstream_gr(genes, divergent, flip(genes$strand)), cl_gr))
names(feat_hits)[4] <- sprintf("Opposite-strand cluster within %d bp (head-to-head)", divergent)
if (length(enh_gr)) feat_hits[["Enhancer in the window"]] <-
  hit_genes(w_same, enh_gr, ignore.strand = TRUE)

feat <- rbindlist(lapply(names(feat_hits), function(f) {
  flag <- seq_len(nrow(genes)) %in% feat_hits[[f]]
  genes[, .(feature = f, n = sum(flag[.I]), genes = .N), by = group]
}))
feat[, pct := 100 * n / genes]
feat[, group_label := factor(sprintf("%s (%d genes)", grp_lab[group], genes),
                             levels = unique(sprintf("%s (%d genes)", grp_lab[group], genes)[
                               order(match(group, names(grp_lab)))]))]
fwrite(feat[, .(group = grp_lab[group], genes, feature, n, pct)], out_feat, sep = "\t")
message("Written: ", out_feat)
print(dcast(feat, feature ~ group, value.var = "pct"), digits = 3)

# ── sweep figure ──────────────────────────────────────────────────────────────
# Panel A: for every gene with a cluster within the longest grid window, the gap
# to its nearest cluster, which is the window length that gene first needs.
nearest <- pairs[gap <= max(grid), .(gap = min(gap)), by = gene]
p_gap <- ggplot(nearest, aes(gap)) +
  stat_ecdf(aes(y = after_stat(ecdf) * nrow(nearest) / nrow(genes) * 100), pad = FALSE) +
  geom_vline(xintercept = chosen, linetype = "dashed", colour = PROMOTER_COL, linewidth = 0.3) +
  geom_hline(data = tgt, aes(yintercept = 100 * genes_needed / nrow(genes)),
             linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_text(data = tgt, aes(x = 1, y = 100 * genes_needed / nrow(genes),
                            label = sprintf("%d genes (%g %% of %d): %s", genes_needed,
                                            100 * fraction, target,
                                            ifelse(is.na(x_bp), "not reached",
                                                   paste(x_bp, "bp")))),
            hjust = 0, vjust = -0.3, size = 2.7, colour = "grey30") +
  scale_x_log10(breaks = c(1, 5, 10, 15, 20, 30, 40, 50, 100, 300)) +
  labs(x = "Gap from gene start to the nearest same-strand cluster (bp, log scale)",
       y = "Genes reached (% of all genes)",
       title = "A. Window length each gene needs",
       subtitle = sprintf("Dashed: windows %s bp", paste(chosen, collapse = ", "))) +
  theme_minimal(base_size = 11)

metric_labels <- c(pct_genes = "Genes with a promoter (%)",
                   pct_reads = "Clustered reads in promoter clusters (%)",
                   clusters_per_gene = "Promoter clusters per gene",
                   median_promoter_bp = "Median promoter length (bp)",
                   pct_bp_in_other_gene = "Promoter bases inside another gene (%)")
long <- melt(sweep, id.vars = c("scheme", "max_window"), measure.vars = names(metric_labels))
long[, variable := factor(metric_labels[as.character(variable)], levels = metric_labels)]
p_sweep <- ggplot(long, aes(max_window, value, colour = scheme)) +
  geom_vline(xintercept = max(chosen), linetype = "dashed", colour = PROMOTER_COL, linewidth = 0.3) +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_x_log10(breaks = grid) +
  scale_colour_manual(values = SCHEME_COLS,
                      labels = c(fixed = "one fixed window",
                                 growing = sprintf("growing windows from %d bp", grow_from))) +
  labs(x = "Longest window (bp, log scale)", y = NULL, colour = NULL,
       title = "B. Window settings compared") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

decided  <- chosen_p[, .(window = window[1]), by = gene][, .N, by = window]
p_dec <- ggplot(decided, aes(factor(window), N)) +
  geom_col(fill = PROMOTER_COL, width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 3) +
  labs(x = "Deciding window (bp)", y = "Genes",
       title = "C. Window that decided each gene",
       subtitle = sprintf("Windows %s bp; %d genes", paste(chosen, collapse = ", "),
                          sum(decided$N))) +
  theme_minimal(base_size = 11)

p_feat <- ggplot(feat, aes(pct, feature, fill = group_label)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge2(width = 0.75, reverse = TRUE),
            hjust = -0.15, size = 2.7) +
  scale_fill_manual(values = setNames(c(PROMOTER_COL, "#f4a261", "grey60"), levels(feat$group_label))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Genes (%)", y = NULL, fill = NULL,
       title = sprintf("D. What lies in the %d bp upstream of each gene, by gene group", x_feat)) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

pdf(out_sweep, width = 13, height = 15)
print((p_gap | p_dec) / p_sweep / p_feat + plot_layout(heights = c(1, 1.6, 1)))
invisible(dev.off())
message("Written: ", out_sweep)

# ── definition figure ─────────────────────────────────────────────────────────
# Every locus is drawn in gene orientation, with x the offset from the gene
# start (negative upstream), so + and - genes read alike.
draw_locus <- function(ctss, clusters, windows, decided_w, promoter, title,
                       xlim, gene_label = "gene") {
  n_w <- length(windows)
  win <- data.table(w = windows, y = -0.25 - 0.13 * seq_len(n_w))
  win[, state := factor(fifelse(w < decided_w, "reaches no cluster",
                        fifelse(w == decided_w, "first to reach a cluster", "not tried")),
                        levels = c("reaches no cluster", "first to reach a cluster", "not tried"))]
  y_prom <- min(win$y) - 0.25
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_segment(data = ctss, aes(x = x, xend = x, y = 0, yend = h), linewidth = 0.6) +
    geom_segment(data = clusters, aes(x = from, xend = to, y = -0.12, yend = -0.12),
                 linewidth = 1.6, colour = "grey30") +
    # show.legend = TRUE draws a key for a state the panel does not use.
    geom_segment(data = win, aes(x = -w, xend = -1, y = y, yend = y, colour = state),
                 linewidth = 0.9, show.legend = TRUE) +
    geom_text(data = win, aes(x = 1.5, y = y, label = paste(w, "bp")),
              hjust = 0, size = 2.8, colour = "grey30") +
    geom_segment(data = promoter, aes(x = from, xend = to, y = y_prom, yend = y_prom),
                 colour = PROMOTER_COL, linewidth = 1.4) +
    annotate("text", x = mean(unlist(promoter[1, .(from, to)])), y = y_prom - 0.1,
             label = "promoter", colour = PROMOTER_COL, size = 3.3) +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 0.55, linewidth = 0.6) +
    annotate("segment", x = 0, xend = xlim[2] * 0.6, y = 0.55, yend = 0.55, linewidth = 0.6,
             arrow = arrow(length = unit(0.12, "cm"))) +
    annotate("text", x = 1, y = 0.65, label = gene_label, hjust = 0, size = 3) +
    scale_colour_manual(values = c("reaches no cluster" = "grey70",
                                   "first to reach a cluster" = "black",
                                   "not tried" = "grey88"),
                        limits = levels(win$state), drop = FALSE) +
    coord_cartesian(xlim = xlim, ylim = c(y_prom - 0.18, 1.05)) +
    labs(x = "Position relative to the gene start (bp)", y = NULL, colour = NULL,
         title = title) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "bottom")
}

# Panel A: the schematic, with invented CTSS heights.
sch_ctss <- data.table(x = c(-48, -44, -41, -38, -35, -31, -27, -24),
                       h = c(0.3, 0.13, 0.2, 0.42, 0.16, 0.9, 0.13, 0.15))
sch_cl   <- data.table(from = -48, to = -24)
sch_w    <- sort(unique(c(chosen, if (max(chosen) < 30) 30L)))
sch_dec  <- sch_w[sch_w >= 24][1]
p_sch <- draw_locus(sch_ctss, sch_cl, sch_w, sch_dec,
                    data.table(from = -48, to = -1),
                    sprintf("A. Definition: windows end at the gene start and are tried shortest first; the first that\n    reaches a same-strand TSS cluster decides, and the promoter runs from that cluster to the gene start"),
                    xlim = c(-55, 25))

# Panels B-D: real loci, chosen as the deepest examples of three situations
# among genes whose promoter clusters end within 20 bp past the gene start and
# span at most 150 bp, so that the panel shows the promoter rather than reads
# spread along a gene body.
# strand is a single value, so this is an if rather than an ifelse, which would
# return only the first element.
to_rel <- function(pos, strand, gs) if (strand == "+") pos - gs else gs - pos
real_locus <- function(g, title) {
  gp <- pairs[gene == g]
  gs <- gp$gene_start[1]; st <- gp$strand[1]
  pr <- chosen_p[gene == g]
  dw <- pr$window[1]
  # Every cluster within the longest window is drawn, so the panel shows which
  # ones the growing rule dropped.
  cls <- gp[gap <= max(chosen)]
  m <- ext$members[cluster_id %in% cls$cluster_id]
  m[, x := to_rel(pos, st, gs)]
  # Square-root heights, so a minor CTSS stays visible beside the dominant one.
  m[, h := 0.9 * sqrt(n_total / max(n_total))]
  clr <- cls[, .(a = to_rel(cluster_start, st, gs), b = to_rel(cluster_end, st, gs))][
    , .(from = pmin(a, b), to = pmax(a, b))]
  prr <- pr[, .(a = to_rel(start, st, gs), b = to_rel(end, st, gs))][
    , .(from = min(pmin(a, b)), to = max(pmax(a, b)))]
  lo <- min(clr$from, -max(chosen)) - 8
  draw_locus(m[, .(x, h)], clr, chosen, dw, prr,
             sprintf("%s  %s (%s strand), %s reads", title, g, st,
                     format(sum(pr$reads), big.mark = ",")),
             xlim = c(lo, max(clr$to, 0) + 0.3 * (0 - lo)), gene_label = g)
}

by_gene <- chosen_p[, .(window = window[1], reads = sum(reads), n = .N,
                        past  = max(fifelse(strand == "+", cluster_end - gene_start,
                                            gene_start - cluster_start)),
                        width = max(cluster_end - cluster_start + 1L)), by = gene]
by_gene <- by_gene[past <= 20L & width <= 150L]
dropped <- pairs[gap <= max(chosen)][, .(n_all = .N), by = gene]
by_gene <- merge(by_gene, dropped, by = "gene")
pick <- function(sel) by_gene[sel][order(-reads)]$gene[1]
ex <- c(B = pick(by_gene$window == min(chosen) & by_gene$n_all == 1L),
        C = pick(by_gene$window == max(chosen) & by_gene$n_all == 1L),
        D = pick(by_gene$n_all > by_gene$n))
titles <- c(B = "B. Decided by the shortest window:",
            C = "C. Decided by the longest window:",
            D = "D. A further cluster dropped by the growing rule:")
real <- lapply(names(ex)[!is.na(ex)], function(k) real_locus(ex[[k]], titles[[k]]))
message("Example loci: ", paste(names(ex), ex, sep = "=", collapse = ", "))

pdf(out_def, width = 9, height = 4.2 * (1 + length(real)))
# One legend, under the last panel: the panels share their window colours.
panels <- c(list(p_sch), real)
panels[-length(panels)] <- lapply(panels[-length(panels)],
                                  function(p) p + theme(legend.position = "none"))
print(wrap_plots(panels, ncol = 1))
invisible(dev.off())
message("Written: ", out_def)

# ── overhang figure ───────────────────────────────────────────────────────────
# A cluster that runs past the gene start keeps its full extent, so the bases it
# covers from the start onward are labelled Promoter although, with an
# annotation whose starts are start codons, they are coding sequence. The
# overhang of a gene-cluster pair is the number of those bases, and its read
# share is the fraction of the cluster's pooled reads that lie on them.
mem <- ext$members[cluster_id %in% chosen_p$cluster_id, .(cluster_id, pos, n_total)]
ov <- chosen_p[, .(cluster_id, gene, strand, gene_start, reads,
                   overhang = pmax(0L, fifelse(strand == "+", cluster_end - gene_start + 1L,
                                               gene_start - cluster_start + 1L)))]
past <- mem[ov, on = "cluster_id", allow.cartesian = TRUE][
  , .(past_reads = sum(n_total[fifelse(strand == "+", pos >= gene_start, pos <= gene_start)])),
  by = .(cluster_id, gene)]
ov <- merge(ov, past, by = c("cluster_id", "gene"))
ov[, pct_past := 100 * past_reads / reads]
crossing <- ov[overhang > 0L]
med_over <- median(crossing$overhang)
message("Pairs running past the gene start: ", nrow(crossing), " of ", nrow(ov),
        "; median overhang ", med_over, " bp; over 50 bp: ", sum(crossing$overhang > 50L),
        "; reads past the start: ",
        round(100 * sum(ov$past_reads) / sum(ov$reads), 1), " % of promoter-cluster reads")

p_hist <- ggplot(crossing, aes(overhang)) +
  geom_histogram(bins = 40, fill = "grey40") +
  geom_vline(xintercept = c(med_over, 50), linetype = "dashed",
             colour = PROMOTER_COL, linewidth = 0.3) +
  annotate("text", x = c(med_over, 50), y = Inf, vjust = c(1.5, 3.2), hjust = -0.1, size = 3,
           colour = PROMOTER_COL,
           label = c(sprintf("median %g bp", med_over),
                     sprintf("> 50 bp: %d", sum(crossing$overhang > 50L)))) +
  scale_x_log10() +
  labs(x = "Bases past the gene start (log scale)", y = "Gene-cluster pairs",
       title = "A. How far promoter clusters run past the gene start",
       subtitle = sprintf("%d of %d pairs run past the start; the other %d end before it",
                          nrow(crossing), nrow(ov), nrow(ov) - nrow(crossing))) +
  theme_minimal(base_size = 11)

p_share <- ggplot(crossing, aes(overhang, pct_past)) +
  geom_point(aes(size = reads), alpha = 0.25, stroke = 0) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = PROMOTER_COL, linewidth = 0.3) +
  scale_x_log10() + scale_size_area(max_size = 4, guide = "none") +
  labs(x = "Bases past the gene start (log scale)",
       y = "Cluster reads past the gene start (%)",
       title = "B. Share of a cluster's reads that lie past the start",
       subtitle = "Point area is the cluster's pooled depth") +
  theme_minimal(base_size = 11)

# Real loci: the deepest pair near the median overhang, the deepest beyond 50 bp
# that still fits a readable panel, and the longest overhang of all.
pick_over <- function(sel) crossing[sel][order(-reads)]$gene[1]
ex_over <- c(C = pick_over(abs(crossing$overhang - med_over) <= 2),
             D = pick_over(crossing$overhang > 50L & crossing$overhang <= 150L),
             E = crossing[which.max(overhang)]$gene)
titles_over <- c(C = "C. Near the median overhang:",
                 D = "D. More than 50 bp past the start:",
                 E = "E. The longest overhang:")
real_over <- lapply(names(ex_over)[!is.na(ex_over)], function(k) {
  g <- ex_over[[k]]
  o <- crossing[gene == g][which.max(overhang)]
  real_locus(g, titles_over[[k]]) +
    annotate("rect", xmin = 0, xmax = o$overhang - 1, ymin = -Inf, ymax = Inf,
             fill = PROMOTER_COL, alpha = 0.08) +
    labs(subtitle = sprintf("Shaded: %d bp past the gene start, holding %.1f %% of the cluster's reads",
                            o$overhang, o$pct_past))
})
message("Overhang loci: ", paste(names(ex_over), ex_over, sep = "=", collapse = ", "))

over_panels <- real_over
over_panels[-length(over_panels)] <- lapply(over_panels[-length(over_panels)],
                                            function(p) p + theme(legend.position = "none"))
pdf(out_over, width = 10, height = 4.5 + 4.2 * length(real_over))
print(((p_hist | p_share) / wrap_plots(over_panels, ncol = 1)) +
        plot_layout(heights = c(1, 1.15 * length(real_over))))
invisible(dev.off())
message("Written: ", out_over)

sessionInfo()
