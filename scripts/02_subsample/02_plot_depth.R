#!/usr/bin/env Rscript
# Stage 2, step 2: plot sequencing depth per library against the exclusion
# threshold and the subsampling target.
#
# The upper panel shows the uniquely mapped read count of every library, sorted
# from shallowest to deepest, with the Tukey lower fence (Q1 - 1.5*IQR) drawn as
# the exclusion threshold and the subsampling target drawn as the depth every
# retained library is reduced to. The lower panel shows the same depths as a
# boxplot on the same axis, which is what makes the fence readable as a property
# of the distribution rather than as an arbitrary number.
#
# Usage (from project root):
#   Rscript scripts/02_subsample/02_plot_depth.R \
#       --depth  results/qc/library_depth.tsv \
#       --stats  results/qc/depth_stats.tsv \
#       --outdir results/figures

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

depth_f <- get_arg("--depth",  "results/qc/library_depth.tsv")
stats_f <- get_arg("--stats",  "results/qc/depth_stats.tsv")
outdir  <- get_arg("--outdir", "results/figures")

stopifnot(file.exists(depth_f), file.exists(stats_f))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

.f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source(if (length(.f))
         file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
       else "scripts/utils/plot_helpers.R")

# ── input ─────────────────────────────────────────────────────────────────────
depth <- read_tsv(depth_f, show_col_types = FALSE) %>%
  mutate(
    retained = tolower(as.character(retained)) %in% c("yes", "true", "1"),
    status   = factor(ifelse(retained, "retained", "excluded"),
                      levels = c("retained", "excluded")),
    sample_id = factor(sample_id, levels = sample_id[order(depth)])
  )

stats <- read_tsv(stats_f, show_col_types = FALSE)
sv <- function(name) {
  v <- stats$value[stats$metric == name]
  if (length(v) == 0) NA_real_ else as.numeric(v[1])
}

fence  <- sv("tukey_lower_fence")
target <- sv("subsample_target")
q1     <- sv("q1")
q3     <- sv("q3")

status_cols <- c(retained = "#457b9d", excluded = "#adb5bd")
line_cols   <- c("Tukey lower fence (Q1 - 1.5xIQR)" = "#e63946",
                 "subsampling target"              = "#2a9d8f")

# The two reference lines carry the whole point of the figure, so the axis is
# extended to hold whichever of them falls outside the range of the data.
x_max <- max(c(depth$depth, fence, target), na.rm = TRUE) * 1.18

# A wide IQR pushes Q1 - 1.5*IQR below zero, and the stage clamps it there. The
# line then coincides with the axis and is easy to miss, so the clamping is
# stated in words rather than left to the reader to infer from the fence value.
fence_note <- if (!is.na(fence) && fence <= 0)
  " (clamped to zero, so no library could be excluded)" else ""

guide_df <- tibble::tibble(
  value = c(fence, target),
  label = names(line_cols)
) %>% filter(!is.na(value))

# ── panel A: per-library depth ────────────────────────────────────────────────
p_bars <- ggplot(depth, aes(x = sample_id, y = depth, fill = status)) +
  geom_col(width = 0.7) +
  geom_hline(data = guide_df,
             aes(yintercept = value, colour = label),
             linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = fmt_n(depth)), hjust = -0.15, size = 2.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, x_max),
                     expand = expansion(mult = c(0, 0))) +
  # drop = TRUE so that a run with no exclusion does not carry an empty
  # "excluded" key; the factor levels are fixed, so the colours stay stable.
  scale_fill_manual(values = status_cols, drop = TRUE, name = NULL) +
  scale_colour_manual(values = line_cols, name = NULL) +
  labs(x = NULL, y = "uniquely mapped reads",
       title = "Sequencing depth per library",
       subtitle = sprintf("%d of %d libraries retained; each retained library is subsampled to %s reads",
                          sum(depth$retained), nrow(depth), fmt_n(target))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())

# ── panel B: the distribution the threshold comes from ────────────────────────
p_box <- ggplot(depth, aes(x = "", y = depth)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = "#f1faee") +
  geom_jitter(aes(colour = status), width = 0.08, height = 0, size = 2, alpha = 0.85) +
  geom_hline(data = guide_df,
             aes(yintercept = value), colour = line_cols[guide_df$label],
             linetype = "dashed", linewidth = 0.6) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma, limits = c(0, x_max),
                     expand = expansion(mult = c(0, 0))) +
  scale_colour_manual(values = status_cols, guide = "none") +
  labs(x = NULL, y = "uniquely mapped reads",
       subtitle = sprintf("distribution across libraries: Q1 = %s, Q3 = %s, fence = %s%s",
                          fmt_n(q1), fmt_n(q3), fmt_n(fence), fence_note)) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Panel height follows the library count so that the bars keep a constant
# thickness whether the run holds 8 yeast libraries or 18 human ones.
pdf_h <- 2.6 + 0.28 * nrow(depth)

out_pdf <- file.path(outdir, "02_library_depth.pdf")
ggsave(out_pdf,
       p_bars / p_box + plot_layout(heights = c(nrow(depth), 3.5)),
       width = 9, height = pdf_h, limitsize = FALSE)

message("Wrote ", out_pdf)

sessionInfo()
