#!/usr/bin/env Rscript
# Stage 08 step 3: test whether de novo sequence architectures partition by
# 5' end type, and characterise each architecture.
#
# The central test: seqArchR clustered the promoter windows without ever seeing
# the end-type labels. If M-type TSSs were simply 1Sg TSSs that failed to gain
# a non-templated G, architecture and end type would be independent. A
# significant association means M-type ends sit in a distinguishable sequence
# context.
#
# Input : results/seqarchr/seq_cluster_assignment.tsv
#         results/seqarchr/window_metadata.tsv
#         results/seqarchr/promoter_windows.fa
# Output: results/seqarchr/arch_by_endtype.tsv, arch_chisq.tsv,
#         arch_summary.tsv, results/figures/16_*.pdf

suppressMessages({
    library(optparse)
    library(Biostrings)
    library(ggplot2)
    library(ggseqlogo)
})

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.end_types_path <- (if (length(.f))
  file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
else "scripts/utils/end_types.R")
if (!file.exists(.end_types_path)) .end_types_path <- "scripts/utils/end_types.R"
source(.end_types_path)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--outdir",     type = "character"),
    make_option("--figdir",     type = "character"),
    make_option("--upstream",   type = "integer", default = 150L),
    make_option("--downstream", type = "integer", default = 50L),
    make_option("--force",      action = "store_true", default = FALSE)
)))

dir.create(opt$figdir, recursive = TRUE, showWarnings = FALSE)
CANON <- "#e63946"   # canonical-initiator colour, per project convention

assign_df <- read.delim(file.path(opt$outdir, "seq_cluster_assignment.tsv"),
                        stringsAsFactors = FALSE)
meta <- read.delim(file.path(opt$outdir, "window_metadata.tsv"),
                   stringsAsFactors = FALSE)
seqs <- readDNAStringSet(file.path(opt$outdir, "promoter_windows.fa"))

df <- merge(meta, assign_df, by = "cluster_id")
stopifnot(nrow(df) > 0)
message(sprintf("Merged %d sequences across %d architectures",
                nrow(df), length(unique(df$arch_cluster))))

# ---------------------------------------------------------------------------
# 1. Architecture x end type
# ---------------------------------------------------------------------------
# Restrict the formal test to 1Sg and M: 2Sgg/3Sggg-dominant TSSs are too few
# to contribute meaningful expected counts, and would only inflate the number
# of sparse cells.
tab_full <- table(df$arch_cluster, df$end_type)
write.table(as.data.frame.matrix(tab_full),
            file.path(opt$outdir, "arch_by_endtype.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

tab <- tab_full[, colnames(tab_full) %in% c("1Sg", "M"), drop = FALSE]
tab <- tab[rowSums(tab) >= 20, , drop = FALSE]

cs <- chisq.test(tab)
message(sprintf("Chi-square: X2 = %.1f, df = %d, p = %.3g",
                cs$statistic, cs$parameter, cs$p.value))

# Cramer's V as an effect size, so a large N does not masquerade as a strong effect
cramers_v <- sqrt(cs$statistic / (sum(tab) * (min(dim(tab)) - 1)))
message(sprintf("Cramer's V = %.3f", cramers_v))

# Per-architecture residuals: which architectures are M-enriched or depleted
resid_df <- as.data.frame.table(cs$stdres, responseName = "std_residual")
names(resid_df)[1:2] <- c("arch_cluster", "end_type")

# Per-architecture binomial test of M fraction vs the global M fraction
global_m <- sum(tab[, "M"]) / sum(tab)
per_arch <- do.call(rbind, lapply(rownames(tab), function(a) {
    n_m <- tab[a, "M"]; n_tot <- sum(tab[a, ])
    bt  <- binom.test(n_m, n_tot, p = global_m)
    data.frame(arch_cluster = a, n = n_tot, n_M = n_m,
               pct_M = 100 * n_m / n_tot,
               p_value = bt$p.value, stringsAsFactors = FALSE)
}))
per_arch$p_adj <- p.adjust(per_arch$p_value, method = "BH")
per_arch$enrichment <- ifelse(per_arch$pct_M > 100 * global_m, "M-enriched", "M-depleted")
per_arch$enrichment[per_arch$p_adj >= 0.05] <- "ns"

write.table(data.frame(
    statistic = c("chisq_X2", "chisq_df", "chisq_p", "cramers_V", "global_pct_M"),
    value     = c(cs$statistic, cs$parameter, cs$p.value, cramers_v, 100 * global_m)),
    file.path(opt$outdir, "arch_chisq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(per_arch, file.path(opt$outdir, "arch_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Figures
# ---------------------------------------------------------------------------
prop_df <- as.data.frame.table(prop.table(tab_full, 1) * 100,
                               responseName = "pct")
names(prop_df)[1:2] <- c("arch_cluster", "end_type")

pdf(file.path(opt$figdir, "16_arch_endtype_composition.pdf"), width = 9, height = 5)
print(
    ggplot(prop_df, aes(arch_cluster, pct, fill = end_type)) +
        geom_col() +
        geom_hline(yintercept = 100 * global_m, linetype = "dashed", colour = "grey20") +
        scale_fill_manual(values = END_TYPE_COLORS_BRIGHT) +
        labs(x = "seqArchR architecture", y = "% of TSSs",
             fill = "Dominant end type",
             title = "End-type composition of de novo sequence architectures",
             subtitle = sprintf("Dashed line = global M fraction (%.1f%%); chi-square p = %.3g, Cramer's V = %.3f",
                                100 * global_m, cs$p.value, cramers_v)) +
        theme_bw()
)
dev.off()

pdf(file.path(opt$figdir, "16_arch_endtype_residuals.pdf"), width = 8, height = 5)
print(
    ggplot(subset(resid_df, end_type == "M"),
           aes(reorder(arch_cluster, std_residual), std_residual,
               fill = std_residual > 0)) +
        geom_col() +
        geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#457b9d", "FALSE" = CANON), guide = "none") +
        labs(x = "Architecture", y = "Standardised residual (M-type)",
             title = "Which architectures are enriched or depleted for M-type ends",
             subtitle = "Beyond +/-2 indicates a departure from independence") +
        theme_bw()
)
dev.off()

# Sequence logo per architecture, so each cluster can be read as a motif
arch_levels <- sort(unique(df$arch_cluster))
seq_by_arch <- lapply(arch_levels, function(a) {
    ids <- df$cluster_id[df$arch_cluster == a]
    as.character(seqs[ids])
})
names(seq_by_arch) <- sprintf("%s (n=%d)", arch_levels, lengths(seq_by_arch))

# ggseqlogo indexes positions 1..width, so relabel them to TSS-relative
# coordinates and thin the breaks, otherwise 200 ticks overplot into a black bar.
win_w    <- opt$upstream + opt$downstream
rel_pos  <- seq(-opt$upstream, opt$downstream - 1L)
brk_rel  <- seq(-opt$upstream, opt$downstream - 1L, by = 25)
brk_idx  <- match(brk_rel, rel_pos)

pdf(file.path(opt$figdir, "16_arch_seqlogos.pdf"),
    width = 11, height = 2.2 * length(arch_levels))
print(
    ggseqlogo(seq_by_arch, method = "bits", ncol = 1) +
        scale_x_continuous(breaks = brk_idx, labels = brk_rel,
                           expand = expansion(mult = 0.01)) +
        geom_vline(xintercept = opt$upstream + 0.5,
                   linetype = "dashed", colour = "grey40") +
        labs(title = "Sequence architecture per seqArchR cluster",
             subtitle = "Dashed line = dominant TSS",
             x = "Position relative to dominant TSS (bp)") +
        theme(axis.text.x = element_text(size = 7))
)
dev.off()

# Zoomed core-promoter view: the initiator region is where the architectures
# actually differ, and it is invisible at full window width.
zoom_lo  <- max(1L, opt$upstream - 35L + 1L)
zoom_hi  <- min(win_w, opt$upstream + 21L)
seq_zoom <- lapply(seq_by_arch, function(s) substr(s, zoom_lo, zoom_hi))
zoom_rel <- rel_pos[zoom_lo:zoom_hi]

pdf(file.path(opt$figdir, "16_arch_seqlogos_zoom.pdf"),
    width = 13, height = 2.0 * length(arch_levels))
print(
    ggseqlogo(seq_zoom, method = "bits", ncol = 1) +
        scale_x_continuous(breaks = seq_along(zoom_rel), labels = zoom_rel) +
        geom_vline(xintercept = which(zoom_rel == 0) - 0.5,
                   linetype = "dashed", colour = "grey40") +
        labs(title = "Core promoter region per architecture (TSS -35 to +20 bp)",
             subtitle = "Dashed line = dominant TSS",
             x = "Position relative to dominant TSS (bp)") +
        theme(axis.text.x = element_text(size = 7))
)
dev.off()

# Initiator dinucleotide usage per architecture
init_tab <- prop.table(table(df$arch_cluster, df$initiator), 1) * 100
init_df  <- as.data.frame.table(init_tab, responseName = "pct")
names(init_df)[1:2] <- c("arch_cluster", "initiator")
init_df$is_yr <- init_df$initiator %in% c("CA", "CG", "TA", "TG")

pdf(file.path(opt$figdir, "16_arch_initiator_usage.pdf"), width = 10, height = 5)
print(
    ggplot(init_df, aes(initiator, pct, fill = is_yr)) +
        geom_col() +
        facet_wrap(~ arch_cluster) +
        scale_fill_manual(values = c("TRUE" = CANON, "FALSE" = "grey60"), guide = "none") +
        labs(x = "Initiator dinucleotide", y = "% of TSSs in architecture",
             title = "Initiator usage per architecture (YR initiators in red)") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
)
dev.off()

# GC content per architecture, split by end type
df$gc <- 100 * letterFrequency(seqs[df$cluster_id], "GC", as.prob = TRUE)[, 1]
pdf(file.path(opt$figdir, "16_arch_gc_content.pdf"), width = 9, height = 5)
print(
    ggplot(subset(df, end_type %in% c("1Sg", "M")),
           aes(arch_cluster, gc, fill = end_type)) +
        geom_boxplot(outlier.size = 0.3) +
        scale_fill_manual(values = c("1Sg" = CANON, "M" = "#457b9d")) +
        labs(x = "Architecture", y = "GC content of window (%)",
             fill = "Dominant end type",
             title = "GC content per architecture and end type") +
        theme_bw()
)
dev.off()

message("Wrote tables to ", opt$outdir, " and figures 16_*.pdf to ", opt$figdir)
print(per_arch)

sessionInfo()
