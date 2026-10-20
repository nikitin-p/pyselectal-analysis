#!/usr/bin/env Rscript
# Draw the mismatch profile of the M class against a baseline end type.
#
# The tables come from scripts/99_additional/mismatch_profile.py, which counts
# every mismatch between an aligned read base and the genome at its position
# from the 5' end of the alignment. The figure has four panels:
#   A  the mismatch rate at each read position, the focus type over the baseline;
#   B  the read prefixes that carry the position-3 mismatches of the focus type;
#   C  the mismatching read base by position, as an excess over the plateau;
#   D  the number of mismatches per read, by sample.
#
# Positions 1 and 2 cannot mismatch under STAR --alignEndsType Local, because a
# mismatch there is soft-clipped instead of aligned, so position 3 is the first
# one at which a mismatch survives. Panel B asks what survives there: a read
# prefix that spells a three-base clip of the vocabulary is an untemplated
# addition whose first two bases happened to match the genome, and it is
# coloured as that end type.
#
# The plateau is the stretch of read positions given by --plateau, which lies
# clear of both the 5' edge effect and the low-quality cycles at the 3' end.
#
# Usage (from project root):
#   Rscript scripts/99_additional/mismatch_plot.R \
#       [--indir    results/mismatch]                     \
#       [--out      results/figures/mismatch_profile.pdf] \
#       [--focus    M]                                    \
#       [--baseline <the canonical type of end_types.R>]  \
#       [--plateau  11,28]                                \
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
out_f    <- get_arg("--out",      "results/figures/mismatch_profile.pdf")
focus    <- get_arg("--focus",    "M")
baseline <- get_arg("--baseline", END_TYPE_CANONICAL)
plateau  <- as.integer(strsplit(get_arg("--plateau", "11,28"), ",")[[1]])
top_n    <- as.integer(get_arg("--top", "20"))

subs  <- fread(file.path(in_dir, "position_substitutions.tsv"))
cov   <- fread(file.path(in_dir, "position_coverage.tsv"))
reads <- fread(file.path(in_dir, "reads_by_mismatch_count.tsv"))
pre   <- fread(file.path(in_dir, "prefix3.tsv"))

types <- c(focus, baseline)
stopifnot(all(types %in% cov$type))
for (d in list(subs, cov, reads, pre)) d[, type := factor(type, levels = types)]
subs <- subs[!is.na(type)]; cov <- cov[!is.na(type)]
reads <- reads[!is.na(type)]; pre <- pre[!is.na(type)]
type_cols <- END_TYPE_COLORS_BRIGHT[types]

theme_set(theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95", colour = NA),
                plot.title = element_text(face = "bold", size = 10)))

# ── A. mismatch rate by read position ─────────────────────────────────────────

rate_s <- merge(cov, subs[, .(mm = sum(n)), by = .(type, sample_id, pos)],
                by = c("type", "sample_id", "pos"), all.x = TRUE)
rate_s[is.na(mm), mm := 0]
rate_s[, rate := 100 * mm / n]
rate_p <- rate_s[, .(mm = sum(mm), n = sum(n)), by = .(type, pos)]
rate_p[, rate := 100 * mm / n]
fwrite(rate_p, file.path(in_dir, "position_rates.tsv"), sep = "\t")

in_plateau <- function(pos) pos >= plateau[1] & pos <= plateau[2]
plateau_rate <- rate_p[in_plateau(pos), .(rate = 100 * sum(mm) / sum(n)), by = type]
for (t in types) {
  message(sprintf("%s: position 3 %.3f %%, plateau (%d-%d) %.3f %%", t,
                  rate_p[type == t & pos == 3, rate], plateau[1], plateau[2],
                  plateau_rate[type == t, rate]))
}

y_top <- max(rate_s$rate) * 1.05
p_a <- ggplot(rate_s, aes(pos, rate, colour = type)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = "grey88") +
  annotate("text", x = 5, y = y_top * 0.98, hjust = 0, vjust = 1, size = 2.8,
           colour = "grey35", lineheight = 0.9,
           label = "positions 1-2:\na mismatch is clipped,\nnever aligned") +
  geom_line(aes(group = interaction(type, sample_id)), linewidth = 0.25,
            alpha = 0.45) +
  geom_line(data = rate_p, linewidth = 0.9) +
  geom_point(data = rate_p[pos == 3], size = 2) +
  scale_colour_manual(values = type_cols, name = "5' end type") +
  scale_x_continuous(breaks = c(1, 3, seq(5, 50, 5)), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, y_top), expand = c(0, 0)) +
  labs(title = "A  Mismatch rate by read position",
       subtitle = "Bold line: pooled libraries; thin lines: one library each",
       x = "Position from the 5' end of the alignment",
       y = "Mismatching reads (%)") +
  theme(legend.position = "inside", legend.position.inside = c(0.42, 0.85))

# ── B. read prefixes behind the position-3 mismatches ─────────────────────────

pre_f <- pre[type == focus, .(n_reads = sum(n_reads), n_mm = sum(n_pos3_mismatch)),
             by = prefix]
pre_f[, rate := 100 * n_mm / n_reads]
# A prefix that spells a three-base clip of the vocabulary takes that type's colour.
pre_f[, clip_type := paste0("3S", tolower(prefix))]
pre_f[!clip_type %in% END_TYPES_REAL, clip_type := "not a vocabulary clip"]
top <- pre_f[order(-n_mm)][seq_len(min(top_n, .N))]
top[, prefix := factor(prefix, levels = rev(prefix))]
clip_levels <- c(intersect(END_TYPES_REAL, top$clip_type), "not a vocabulary clip")
top[, clip_type := factor(clip_type, levels = clip_levels)]
bg_rate <- rate_p[type == focus & pos >= 4 & pos <= 10, 100 * sum(mm) / sum(n)]

p_b <- ggplot(top, aes(rate, prefix, fill = clip_type)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = bg_rate, linetype = "dashed", colour = "grey30") +
  geom_text(aes(label = fmt_n(n_mm)), hjust = -0.15, size = 2.6,
            colour = "grey25") +
  scale_fill_manual(values = c(END_TYPE_COLORS_BRIGHT,
                               "not a vocabulary clip" = "grey75"),
                    name = "Prefix spells the clip of") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = sprintf("B  %s reads that mismatch at position 3, by read prefix",
                       focus),
       subtitle = sprintf(paste0("The %d prefixes with the most such reads ",
                                 "(count at bar end); dashed: rate at positions 4-10"),
                          nrow(top)),
       x = "Reads with this prefix whose third base mismatches (%)",
       y = "First three read bases (5'->3')") +
  theme(axis.text.y = element_text(family = "mono"),
        panel.grid.major.y = element_blank())

# ── C. mismatching read base by position, excess over the plateau ─────────────

bases <- c("A", "C", "G", "T")
base_r <- subs[read %in% bases, .(mm = sum(n)), by = .(type, pos, read)]
cov_p  <- cov[, .(n = sum(n)), by = .(type, pos)]
grid   <- CJ(type = factor(types, levels = types), pos = 3:10, read = bases)
heat   <- merge(grid, base_r, by = c("type", "pos", "read"), all.x = TRUE)
heat[is.na(mm), mm := 0]
heat <- merge(heat, cov_p, by = c("type", "pos"))
plat <- merge(base_r[in_plateau(pos), .(mm = sum(mm)), by = .(type, read)],
              cov_p[in_plateau(pos), .(n = sum(n)), by = type], by = "type")
plat[, plat_rate := 100 * mm / n]
heat <- merge(heat, plat[, .(type, read, plat_rate)], by = c("type", "read"))
heat[, excess := 100 * mm / n - plat_rate]
lim <- max(abs(heat$excess))

p_c <- ggplot(heat, aes(factor(pos), factor(read, levels = rev(bases)),
                        fill = excess)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", excess),
                colour = abs(excess) > 0.6 * lim), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                      guide = "none") +
  facet_wrap(~type, ncol = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim),
                       name = "Excess\n(% of reads)") +
  labs(title = "C  Mismatching read base, excess over the plateau",
       subtitle = sprintf("Rate of mismatches to each base minus its mean rate at positions %d-%d",
                          plateau[1], plateau[2]),
       x = "Position from the 5' end of the alignment", y = "Read base") +
  theme(panel.grid = element_blank())

# ── D. mismatches per read, by sample ─────────────────────────────────────────

mm_levels <- c("0", "1", "2", "3+")
reads[, mm_class := factor(ifelse(n_mismatch >= 3, "3+", as.character(n_mismatch)),
                           levels = mm_levels)]
per_read <- reads[, .(n = sum(n)), by = .(type, sample_id, mm_class)]
per_read[, pct := 100 * n / sum(n), by = .(type, sample_id)]

p_d <- ggplot(per_read, aes(pct, sample_id, fill = mm_class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = per_read[mm_class == "0"],
            aes(x = pct / 2, label = sprintf("%.1f %%  (%s reads)", pct, fmt_n(n))),
            size = 2.6, colour = "grey15") +
  facet_wrap(~type, ncol = 1) +
  scale_fill_manual(values = c("0" = "#deebf7", "1" = "#9ecae1",
                               "2" = "#4292c6", "3+" = "#08519c"),
                    name = "Mismatches\nper read") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  labs(title = "D  Mismatches per read",
       subtitle = "Label: share and number of reads with no mismatch",
       x = "Reads (%)", y = NULL) +
  theme(panel.grid.major.y = element_blank())

# ── assemble ──────────────────────────────────────────────────────────────────

dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
fig <- (p_a | p_b) / (p_c | p_d) + plot_layout(heights = c(1, 1))
ggsave(out_f, fig, width = 14, height = 11)
message("Wrote ", out_f)

sessionInfo()
