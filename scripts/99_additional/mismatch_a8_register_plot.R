#!/usr/bin/env Rscript
# Draw the register of the upstream A for reads grouped by their first mismatch.
#
# The tables come from scripts/99_additional/mismatch_a8_register.py. The yeast
# A at -8 marks the register in which a group of reads starts, so its position
# tells whether a read that mismatches at aligned position k starts at its first
# base or at base k + 1, with the mismatching base sitting at -1 of the true
# initiator. The figure has five panels:
#   A  the A frequency by coordinate from the alignment start;
#   B  the same profiles after the mismatching group is re-registered on k + 1;
#   C  the re-registered profile split by the mismatching read base;
#   D  the genomic dinucleotide under the mismatch and the base after it;
#   E  the profile of every first-mismatch position, as a control.
#
# Usage (from project root):
#   Rscript scripts/99_additional/mismatch_a8_register_plot.R \
#       [--indir    results/mismatch]                         \
#       [--out      results/figures/mismatch_a8_register.pdf] \
#       [--focus    M]                                        \
#       [--baseline <the canonical type of end_types.R>]      \
#       [--k        3]                                        \
#       [--obs      -8]

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
out_f    <- get_arg("--out",      "results/figures/mismatch_a8_register.pdf")
focus    <- get_arg("--focus",    "M")
baseline <- get_arg("--baseline", END_TYPE_CANONICAL)
k        <- as.integer(get_arg("--k",   "3"))
obs      <- as.integer(get_arg("--obs", "-8"))

prof  <- fread(file.path(in_dir, "a8_register_profile.tsv"),
               colClasses = list(character = "group"))
dinuc <- fread(file.path(in_dir, "a8_register_dinuc.tsv"),
               colClasses = list(character = "group"))
nn    <- fread(file.path(in_dir, "a8_register_reads.tsv"),
               colClasses = list(character = "group"))

# The coordinate system has no zero, so shifting goes through a displacement.
to_disp  <- function(coord) ifelse(coord < 0, coord, coord - 1L)
to_coord <- function(disp)  ifelse(disp < 0, disp, disp + 1L)
shifted_obs <- to_coord(to_disp(obs) + k)   # where A lands if position k + 1 is +1

pct_a <- prof[, .(pct = 100 * sum(n[base == "A"]) / sum(n)),
              by = .(type, group, weight, coord)]
n_lab <- function(t, g) nn[type == t & group == g,
                           sprintf("%s reads, %s sites", fmt_n(n_reads), fmt_n(n_sites))]

mm_col  <- "#1d3557"
kg      <- as.character(k)
x_range <- c(-14L, 6L)
coord_axis <- function(lim) {
  br <- setdiff(seq(lim[1], lim[2]), 0L)
  scale_x_continuous(breaks = to_disp(br), labels = br, expand = c(0.01, 0))
}

theme_set(theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95", colour = NA),
                plot.title = element_text(face = "bold", size = 10)))

# ── A. A frequency by coordinate from the alignment start ─────────────────────

series <- c(sprintf("%s, no mismatch", baseline), sprintf("%s, no mismatch", focus),
            sprintf("%s, first mismatch at %d", focus, k))
series_cols <- setNames(c(END_TYPE_COLORS_BRIGHT[[baseline]],
                          END_TYPE_COLORS_BRIGHT[[focus]], mm_col), series)
three <- rbind(pct_a[type == baseline & group == "none"][, series := series[1]],
               pct_a[type == focus    & group == "none"][, series := series[2]],
               pct_a[type == focus    & group == kg][,     series := series[3]])
three[, series := factor(series, levels = series_cols |> names())]
three[, weight := factor(weight, levels = c("reads", "sites"),
                         labels = c("Every read counts", "Every 5' position counts once"))]
in_x <- function(d, lim = x_range) d[coord >= lim[1] & coord <= lim[2]]

marks <- data.table(x = to_disp(c(obs, shifted_obs)),
                    lab = c(sprintf("%d from\nposition 1", obs),
                            sprintf("%d from\nposition %d", obs, k + 1L)))
p_a <- ggplot(in_x(three), aes(to_disp(coord), pct, colour = series)) +
  geom_vline(data = marks, aes(xintercept = x), linetype = "dashed",
             colour = "grey45") +
  geom_text(data = marks, aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.2, hjust = 1.05, size = 2.5, colour = "grey30",
            lineheight = 0.9) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_wrap(~weight, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = series_cols, name = NULL) +
  coord_axis(x_range) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(title = "A  Upstream A, by coordinate from the alignment start",
       subtitle = sprintf("%s reads with a first mismatch at %d: %s", focus, k,
                          n_lab(focus, kg)),
       x = "Coordinate (+1 is the first aligned base)", y = "A at this coordinate (%)") +
  theme(legend.position = "bottom")

# ── B. the same profiles, re-registered on position k + 1 ─────────────────────

rereg <- copy(three)
mm_series <- series[3]
rereg[series == mm_series, coord := to_coord(to_disp(coord) - k)]
x_b <- c(-14L, 3L)
p_b <- ggplot(in_x(rereg, x_b), aes(to_disp(coord), pct, colour = series)) +
  geom_vline(xintercept = to_disp(obs), linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_wrap(~weight, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = series_cols, name = NULL) +
  coord_axis(x_b) +
  labs(title = sprintf("B  The same, with the mismatching reads re-registered on position %d",
                       k + 1L),
       subtitle = sprintf("If the mismatching base sits at -1, the peaks coincide at %d", obs),
       x = sprintf("Coordinate (+1 is aligned position %d for the mismatching reads)", k + 1L),
       y = "A at this coordinate (%)") +
  theme(legend.position = "bottom")

# ── C. split by the mismatching read base ─────────────────────────────────────

split_groups <- c(G = paste0(kg, "G"), A = paste0(kg, "A"), `C or T` = paste0(kg, "Y"))
split_groups <- split_groups[split_groups %in% pct_a[type == focus, group]]
if (length(split_groups) > 0) {
  sp <- pct_a[type == focus & group %in% split_groups & weight == "reads"]
  sp[, coord := to_coord(to_disp(coord) - k)]
  sp[, read_base := names(split_groups)[match(group, split_groups)]]
  sp[, lab := sprintf("%s  (%s reads)", read_base,
                      fmt_n(nn[type == focus][match(sp$group, group), n_reads]))]
  sp[, lab := factor(lab, levels = unique(lab[order(match(read_base, names(split_groups)))]))]
  p_c <- ggplot(in_x(sp, x_b), aes(to_disp(coord), pct, colour = lab)) +
    geom_vline(xintercept = to_disp(obs), linetype = "dashed", colour = "grey45") +
    geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
    scale_colour_manual(values = c("#e09f3e", "#9e2a2b", "grey55"),
                        name = sprintf("Read base at position %d", k)) +
    coord_axis(x_b) +
    labs(title = "C  Re-registered profile, by the mismatching read base",
         subtitle = "Every read counts",
         x = sprintf("Coordinate (+1 is aligned position %d)", k + 1L),
         y = "A at this coordinate (%)") +
    theme(legend.position = "inside", legend.position.inside = c(0.2, 0.82),
          legend.background = element_blank())
} else {
  p_c <- plot_spacer()
}

# ── D. the genomic dinucleotide under the mismatch ────────────────────────────

yr <- c("CA", "CG", "TA", "TG")
dn <- dinuc[type == focus & group == kg, .(n = sum(n)), by = .(dn = paste0(ref_k, ref_next))]
dn[, pct := 100 * n / sum(n)]
dn[, class := ifelse(dn %in% yr, "YR (canonical initiator)", "other")]
dn[, dn := factor(dn, levels = dn[order(pct)])]
p_d <- ggplot(dn, aes(pct, dn, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.2, size = 2.6,
            colour = "grey25") +
  scale_fill_manual(values = c("YR (canonical initiator)" = "#e63946",
                               "other" = "grey70"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = sprintf("D  Genomic dinucleotide at aligned positions %d and %d", k, k + 1L),
       subtitle = sprintf("%s reads with a first mismatch at %d; YR share %.1f %%",
                          focus, k, dn[class != "other", sum(pct)]),
       x = "Reads (%)", y = NULL) +
  theme(axis.text.y = element_text(family = "mono"),
        panel.grid.major.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.7, 0.15))

# ── E. every first-mismatch position ──────────────────────────────────────────

ks <- sort(as.integer(grep("^[0-9]+$", unique(pct_a[type == focus, group]), value = TRUE)))
ks <- ks[ks >= 3]
ctl <- pct_a[type == focus & group %in% as.character(ks)]
ctl[, kk := as.integer(group)]
strip <- function(kk) sprintf("first mismatch at %d\n%s reads", kk,
                              fmt_n(nn[type == focus][match(as.character(kk), group), n_reads]))
ctl[, facet := factor(strip(kk), levels = strip(ks))]
ctl_marks <- data.table(kk = ks)[, .(x = to_disp(c(obs, to_coord(to_disp(obs) + kk))),
                                     reg = c("position 1 is +1",
                                             "the base after the mismatch is +1")), by = kk]
ctl_marks[, facet := factor(strip(kk), levels = strip(ks))]
p_e <- ggplot(in_x(ctl, c(-14L, -1L)), aes(to_disp(coord), pct)) +
  geom_vline(data = ctl_marks, aes(xintercept = x, linetype = reg), colour = "grey35") +
  geom_line(aes(colour = weight), linewidth = 0.6) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = c(reads = mm_col, sites = "#a8b5c8"),
                      labels = c(reads = "every read", sites = "every 5' position once"),
                      name = NULL) +
  scale_linetype_manual(values = c("position 1 is +1" = "dotted",
                                   "the base after the mismatch is +1" = "dashed"),
                        name = sprintf("Expected place of A%d if", obs)) +
  scale_x_continuous(breaks = seq(-14, -2, 3), expand = c(0.02, 0)) +
  labs(title = sprintf("E  %s reads by the position of their first mismatch", focus),
       x = "Coordinate from the alignment start", y = "A at this coordinate (%)") +
  theme(legend.position = "bottom")

# ── assemble ──────────────────────────────────────────────────────────────────

dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
fig <- (p_a | p_b) / (p_c | p_d) / p_e + plot_layout(heights = c(1.5, 1.1, 0.9))
ggsave(out_f, fig, width = 14, height = 14)
message("Wrote ", out_f)

sessionInfo()
