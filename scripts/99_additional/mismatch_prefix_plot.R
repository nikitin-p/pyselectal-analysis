#!/usr/bin/env Rscript
# Draw which three-base read prefixes carry a mismatch at position 3, ranked by
# enrichment, for the focus type and the baseline type side by side.
#
# The tables come from scripts/99_additional/mismatch_profile.py. Under STAR
# --alignEndsType Local a mismatch at position 1 or 2 is soft-clipped, so a read
# that mismatches at position 3 has the pattern match-match-mismatch, and its
# first three bases are a candidate untemplated addition whose first two bases
# happened to match the genome. Panel B of mismatch_plot.R ranks the prefixes by
# the number of such reads, which favours prefixes that are merely common. This
# figure ranks them by enrichment instead and shows what each one replaced:
#   A  the focus type's prefixes, by enrichment of the position-3 mismatch;
#   B  the same for the baseline type;
#   C  the genomic base that the third read base replaced, for the prefixes of A and B.
#
# Enrichment is the share of a prefix's reads that mismatch at position 3 over
# the type's mismatch rate at positions 4-10, the same background as the dashed
# line of mismatch_plot.R panel B. Prefixes with fewer than --min_n such reads,
# pooled over libraries, are left out so that a handful of reads cannot top the
# ranking. A prefix that spells a three-base clip of the vocabulary is coloured
# as that end type.
#
# Usage (from project root):
#   Rscript scripts/99_additional/mismatch_prefix_plot.R \
#       [--indir    results/mismatch]                     \
#       [--out      results/figures/mismatch_prefix3.pdf] \
#       [--focus    M]                                    \
#       [--baseline <the canonical type of end_types.R>]  \
#       [--min_n    100]                                  \
#       [--top      20]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "..", "utils", "end_types.R"))
source(file.path(script_dir, "..", "utils", "plot_helpers.R"))

in_dir   <- get_arg("--indir",    "results/mismatch")
out_f    <- get_arg("--out",      "results/figures/mismatch_prefix3.pdf")
focus    <- get_arg("--focus",    "M")
baseline <- get_arg("--baseline", END_TYPE_CANONICAL)
min_n    <- as.integer(get_arg("--min_n", "100"))
top_n    <- as.integer(get_arg("--top", "20"))

subs <- fread(file.path(in_dir, "position_substitutions.tsv"))
cov  <- fread(file.path(in_dir, "position_coverage.tsv"))
pre  <- fread(file.path(in_dir, "prefix3.tsv"))
pref <- fread(file.path(in_dir, "prefix3_ref.tsv"))

types <- c(focus, baseline)
stopifnot(all(types %in% cov$type))

theme_set(theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95", colour = NA),
                plot.title = element_text(face = "bold", size = 10)))

# ── enrichment of the position-3 mismatch, by prefix ──────────────────────────

bg <- merge(subs[pos >= 4 & pos <= 10, .(mm = sum(n)), by = type],
            cov[pos >= 4 & pos <= 10, .(n = sum(n)), by = type], by = "type")
bg[, bg_rate := 100 * mm / n]

enr <- pre[type %in% types,
           .(n_reads = sum(n_reads), n_mm = sum(n_pos3_mismatch)),
           by = .(type, prefix)]
enr <- merge(enr, bg[, .(type, bg_rate)], by = "type")
enr[, rate := 100 * n_mm / n_reads]
enr[, enrichment := rate / bg_rate]
enr[, clip_type := paste0("3S", tolower(prefix))]
enr[!clip_type %in% END_TYPES_REAL, clip_type := "not a vocabulary clip"]
fwrite(enr[order(type, -enrichment)], file.path(in_dir, "prefix3_enrichment.tsv"),
       sep = "\t")

top <- enr[n_mm >= min_n][order(-enrichment), head(.SD, top_n), by = type]
for (t in types) {
  message(sprintf("%s: background %.3f %%; top prefixes %s", t,
                  bg[type == t, bg_rate],
                  paste(sprintf("%s (%.1fx)", top[type == t, prefix],
                                top[type == t, enrichment])[1:min(5, top[type == t, .N])],
                        collapse = ", ")))
}

clip_levels <- c(intersect(END_TYPES_REAL, top$clip_type), "not a vocabulary clip")
clip_fill   <- c(END_TYPE_COLORS_BRIGHT, "not a vocabulary clip" = "grey75")

enrich_panel <- function(t, tag) {
  d <- top[type == t]
  d[, prefix := factor(prefix, levels = rev(prefix))]
  d[, clip_type := factor(clip_type, levels = clip_levels)]
  ggplot(d, aes(enrichment, prefix, fill = clip_type)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_text(aes(label = fmt_n(n_mm)), hjust = -0.15, size = 2.6,
              colour = "grey25") +
    scale_fill_manual(values = clip_fill, limits = clip_levels, drop = FALSE,
                      name = "Prefix spells the clip of") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = sprintf("%s  %s reads: prefixes enriched for a position-3 mismatch",
                         tag, t),
         subtitle = sprintf(paste0("Top %d with at least %d such reads (count at bar end); ",
                                   "background %.2f %% at positions 4-10"),
                            nrow(d), min_n, bg[type == t, bg_rate]),
         x = "Position-3 mismatch rate over background (fold)",
         y = "First three read bases (5'->3')") +
    theme(axis.text.y = element_text(family = "mono"),
          panel.grid.major.y = element_blank())
}

p_a <- enrich_panel(focus, "A")
p_b <- enrich_panel(baseline, "B") + guides(fill = "none")  # A carries the legend

# ── genomic base replaced at position 3, for the prefixes above ───────────────

bases <- c("A", "C", "G", "T")
ref_p <- pref[type %in% types, .(n = sum(n)), by = .(type, prefix, ref)]
ref_p <- merge(ref_p, top[, .(type, prefix)], by = c("type", "prefix"))
ref_p[, pct := 100 * n / sum(n), by = .(type, prefix)]
ref_p[, type := factor(type, levels = types)]
ref_p[, row := factor(paste(type, prefix),
                      levels = rev(top[, paste(type, prefix)]))]
ref_p[, ref := factor(ref, levels = bases)]

p_c <- ggplot(ref_p, aes(pct, row, fill = ref)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(pct >= 12, as.character(ref), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.5, colour = "white", fontface = "bold") +
  facet_wrap(~type, nrow = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^\\S+ ", "", x)) +
  scale_fill_manual(values = c(A = "#109648", C = "#255c99", G = "#f7b32b",
                               T = "#d62839"),
                    name = "Genomic base") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "C  Genomic base replaced by the third read base",
       subtitle = "Same prefixes and order as A and B; the third read base itself can never be the genomic one",
       x = "Position-3 mismatches of the prefix (%)",
       y = "First three read bases (5'->3')") +
  theme(axis.text.y = element_text(family = "mono"),
        panel.grid.major.y = element_blank())

# ── assemble ──────────────────────────────────────────────────────────────────

dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
fig <- (p_a | p_b) / p_c + plot_layout(guides = "collect")
ggsave(out_f, fig, width = 14, height = 11)
message("Wrote ", out_f)

sessionInfo()
