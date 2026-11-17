#!/usr/bin/env Rscript
# Stage 22 step 11: which other 5' end types are produced at the sites of the
# long clips.
#
# A site is an exact aligned 5' end (chr, pos, strand) in the pooled fine-typed
# table of Stage 17. Every site carrying a long clip and at least one other read
# is matched to three sites without a long clip that carry the same number of
# other reads (exact up to 50 reads, quarter-log2 bins above), drawn with
# replacement under a fixed seed. Presence (a site carries the type at all) is
# compared rather than read counts, because single hotspots dominate the read
# counts of several types (chrXIV:230672(-) alone holds 64 % of the 1St reads).
#
#   A  log2 ratio of the share of sites carrying each fine type (>= 5,000 reads
#      pooled) at long-clip sites against the matched control, by clip length;
#      a clip at least 2/3 A is marked A-rich;
#   B  share of 2-4 nt clip reads at a long-clip site whose clip is a suffix of
#      that site's most abundant long clip, by origin class group, against the
#      same with the long clips shuffled across sites of the same group.
#
# Outputs: <outdir>/lc_coprod_presence.tsv, lc_coprod_nested.tsv;
#          <figdir>/31_longclip_coproduction.{pdf,png}
#
# Usage (from project root):
#   Rscript scripts/14_longclip/11_coproduction.R \
#       [--fine results/dominance/fine_typed_ctss.tsv.gz] \
#       [--reads results/longclip/longclip_reads.tsv.gz] \
#       [--outdir results/longclip] [--figdir results/figures] [--force]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i)) args[i + 1] else default
}
fine_f  <- get_arg("--fine",   "results/dominance/fine_typed_ctss.tsv.gz")
reads_f <- get_arg("--reads",  "results/longclip/longclip_reads.tsv.gz")
outdir  <- get_arg("--outdir", "results/longclip")
figdir  <- get_arg("--figdir", "results/figures")
force   <- "--force" %in% args

fig_f <- file.path(figdir, "31_longclip_coproduction.pdf")
if (file.exists(fig_f) && !force) {
  message("Output exists, skipping (use --force to rerun): ", fig_f)
  quit(save = "no", status = 0)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "lc_helpers.R"))

set.seed(1L)
MIN_TYPE_READS <- 5000L   # fine types tested in panel A
N_CTRL         <- 3L      # control sites per long-clip site
A_RICH         <- 2 / 3   # clip fraction of A that marks a clip A-rich
LC_TYPE        <- "longclip"
UNCLIPPED      <- END_TYPES_REAL[!grepl("S", END_TYPES_REAL)]

# ── sites ─────────────────────────────────────────────────────────────────────
f <- fread(fine_f)[, .(n = sum(n)), by = .(chr, pos, strand, type)]
f[, clen := fifelse(type %in% UNCLIPPED, 0L,
                    suppressWarnings(as.integer(sub("^([0-9]+)S[acgt]+$", "\\1", type))))]

site <- f[, .(lc = sum(n[type == LC_TYPE]), oth = sum(n[type != LC_TYPE])),
          by = .(chr, pos, strand)]
lr <- fread(reads_f)
cls <- lr[, .(n = sum(n_reads)), by = .(chr, strand, pos = pos5, class)][
  order(-n)][, .SD[1], by = .(chr, strand, pos)]
site <- merge(site, cls[, .(chr, strand, pos, lc_class = class)],
              by = c("chr", "strand", "pos"), all.x = TRUE)
site[, db := fifelse(oth <= 50, as.numeric(oth), 50 + floor(log2(pmax(oth, 1) / 50) * 4))]

lcs  <- site[lc > 0 & oth > 0]
pool <- site[lc == 0]
ctl <- rbindlist(lapply(split(lcs, lcs$db), function(d) {
  p <- pool[db == d$db[1]]
  if (!nrow(p)) return(NULL)
  p[sample(.N, nrow(d) * N_CTRL, replace = TRUE)]
}))
message("Long-clip sites with other reads: ", nrow(lcs),
        " (median other reads ", median(lcs$oth), "); control sites: ", nrow(ctl),
        " (median ", median(ctl$oth), ")")

# ── A: presence enrichment ────────────────────────────────────────────────────
types <- f[type != LC_TYPE & !is.na(clen), .(n = sum(n)), by = type][n >= MIN_TYPE_READS, type]
presence <- function(s) {
  x <- f[s[, .(chr, pos, strand, id = .I)], on = .(chr, pos, strand)][type %in% types]
  x[, .(sites = uniqueN(id)), by = type][, pct := 100 * sites / nrow(s)][]
}
pr <- merge(presence(lcs)[, .(type, sites_lc = sites, pct_lc = pct)],
            presence(ctl)[, .(type, sites_ctrl = sites, pct_ctrl = pct)],
            by = "type", all = TRUE)
for (v in names(pr)[-1]) set(pr, which(is.na(pr[[v]])), v, 0)
pr[, log2_ratio := log2((pct_lc + 0.1) / (pct_ctrl + 0.1))]
pr[, clen := f[match(pr$type, f$type), clen]]
pr[, seq := toupper(sub("^[0-9]+S", "", type))]
pr[, a_frac := fifelse(clen == 0L, 0, nchar(gsub("[^A]", "", seq)) / pmax(nchar(seq), 1L))]
pr[, a_rich := clen > 0L & a_frac >= A_RICH]
pr[, `:=`(n_sites_lc = nrow(lcs), n_sites_ctrl = nrow(ctl))]
setorder(pr, -log2_ratio)
fwrite(pr[, .(type, clen, a_rich, sites_lc, n_sites_lc, pct_lc, sites_ctrl,
              n_sites_ctrl, pct_ctrl, log2_ratio)],
       file.path(outdir, "lc_coprod_presence.tsv"), sep = "\t")

# ── B: nesting of short clips in the site's long clip ─────────────────────────
top_clip <- lr[, .(n = sum(n_reads)), by = .(chr, strand, pos = pos5, clip)][
  order(-n)][, .SD[1], by = .(chr, strand, pos)][, n := NULL]
top_clip <- merge(top_clip, cls[, .(chr, strand, pos, class)], by = c("chr", "strand", "pos"))
top_clip[, grp := fifelse(grepl("^splice", class), "spliced", as.character(class))]
sh <- f[clen %in% 2:4][top_clip, on = .(chr, strand, pos), nomatch = 0]
sh[, sseq := toupper(sub("^[0-9]+S", "", type))]
sh[, nested := endsWith(clip, sseq)]
sh[, nested_shuf := endsWith(sample(clip), sseq), by = grp]   # shuffled within the class group
nest <- rbind(
  sh[, .(n_reads = sum(n), pct = 100 * sum(n[nested]) / sum(n), set = "observed"), by = .(grp, clen)],
  sh[, .(n_reads = sum(n), pct = 100 * sum(n[nested_shuf]) / sum(n), set = "long clips shuffled"),
     by = .(grp, clen)])
setorder(nest, grp, clen, set)
fwrite(nest, file.path(outdir, "lc_coprod_nested.tsv"), sep = "\t")
print(nest)

# ── figure ────────────────────────────────────────────────────────────────────
th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "top", legend.key.size = unit(0.35, "cm"))

pa <- copy(pr)[, x := factor(fifelse(clen == 0L, "0 (M)", as.character(clen)),
                             c("0 (M)", "1", "2", "3", "4"))]
pa[, colour := fifelse(a_rich, "A-rich clip", "other")]
lab <- pa[abs(log2_ratio) >= 1.4 | type == END_TYPE_CANONICAL | clen == 0L]
pA <- ggplot(pa, aes(x, log2_ratio, colour = colour)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(size = pct_lc), position = position_jitter(width = 0.15, seed = 1L), alpha = 0.85) +
  geom_text_repel(data = lab, aes(label = type), size = 2.4, show.legend = FALSE,
                  position = position_jitter(width = 0.15, seed = 1L),
                  min.segment.length = 0, max.overlaps = Inf, seed = 1L) +
  scale_colour_manual(values = c("A-rich clip" = "#d97706", other = "#1d4ed8"), name = NULL) +
  scale_size_area(max_size = 4, name = "% of long-clip\nsites carrying it") +
  labs(x = "clip length of the co-produced type (nt)",
       y = "log2(% long-clip sites / % matched control sites)",
       title = "A  Long-clip sites also carry the other long clips",
       subtitle = sprintf(paste("%s long-clip sites against %s control sites with the same number of other reads;",
                                "types with >= %s reads"),
                          format(nrow(lcs), big.mark = ","), format(nrow(ctl), big.mark = ","),
                          format(MIN_TYPE_READS, big.mark = ","))) + th

pb <- copy(nest)[, `:=`(clen = factor(paste0(clen, " nt")),
                        set = factor(set, c("observed", "long clips shuffled")),
                        grp = factor(grp, c("spliced", setdiff(LC_CLASSES, c("splice_annotated", "splice_unannotated")))))]
pb <- pb[n_reads >= 100]
pB <- ggplot(pb, aes(clen, pct, fill = set)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "white", linewidth = 0.3) +
  geom_text(data = pb[set == "observed"], aes(y = pct + 4, label = format(n_reads, big.mark = ",")),
            position = position_nudge(x = -0.2), size = 2.2, colour = "grey30") +
  facet_wrap(~grp, nrow = 1) +
  scale_fill_manual(values = c(observed = "#1d4ed8", "long clips shuffled" = "#9ca3af"), name = NULL) +
  scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
  labs(x = "length of the short clip", y = "% of its reads that are a suffix\nof the site's long clip",
       title = "B  Short clips at a long-clip site are often the same tail, shorter",
       subtitle = "number above a bar: short-clip reads; facets with < 100 reads omitted") + th

fig <- pA / pB + plot_layout(heights = c(1.2, 1))
ggsave(fig_f, fig, width = 8, height = 8.5)
ggsave(sub("\\.pdf$", ".png", fig_f), fig, width = 8, height = 8.5, dpi = 150)
message("Written: ", fig_f)
sessionInfo()
