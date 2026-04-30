#!/usr/bin/env Rscript
# Stage 3, step 3: nucleotide composition of non-1Sg soft-clipped bases
# vs literature RT non-templated addition preferences.
#
# Usage (from project root):
#   Rscript scripts/03_five_prime/03_softclip_composition.R \
#       --counts_dir results/five_prime \
#       --samples    config/samples.tsv \
#       --outdir     results/figures \
#       --params     config/params.yaml
#
# Inputs:  results/five_prime/*_counts.tsv  (produced by 01_count.sh)
# Outputs: results/figures/03_softclip_nuc_composition.pdf
#          results/figures/03_softclip_nuc_by_sample.pdf
#          results/figures/03_softclip_nuc_composition.tsv

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(yaml)
})

# ── argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

counts_dir <- get_arg("--counts_dir", "results/five_prime")
samples_f  <- get_arg("--samples",    "config/samples.tsv")
outdir     <- get_arg("--outdir",     "results/figures")
params_f   <- get_arg("--params",     "config/params.yaml")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

params <- yaml.load_file(params_f)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Literature RT non-templated nucleotide addition preferences.
# Source: Zhu et al. 2001 (Biotechniques). Citation needs verification before publication.
# Override in params.yaml under rt_nontemplated_prefs: {A: x, C: x, G: x, T: x}
default_rt_prefs <- list(A = 0.10, C = 0.50, G = 0.25, T = 0.15)
rt_prefs_raw     <- params$rt_nontemplated_prefs %||% default_rt_prefs
rt_prefs_df <- tibble(
  nucleotide = toupper(names(rt_prefs_raw)),
  freq       = unlist(rt_prefs_raw)
)

# ── helpers ───────────────────────────────────────────────────────────────────
NUC_COLORS <- c(A = "#4DAF4A", C = "#377EB8", G = "#FF7F00", T = "#E41A1C")

# Extract the lowercase-base suffix from a soft-clip type label.
# "2Sgc" -> "gc",  "1Sg" -> "g",  "36Mctcac" -> NA
parse_softclip_bases <- function(type) {
  m <- regmatches(type, regexpr("^[0-9]+S([a-z]+)$", type, perl = TRUE))
  if (length(m) == 0) return(NA_character_)
  sub("^[0-9]+S", "", m)
}

# Expand a data frame (type, count, [sample_id]) into per-position rows.
expand_positions <- function(df) {
  lapply(seq_len(nrow(df)), function(i) {
    bases_str <- parse_softclip_bases(df$type[i])
    if (is.na(bases_str)) return(NULL)
    bases <- strsplit(bases_str, "")[[1]]
    out <- tibble(
      position   = seq_along(bases),
      nucleotide = toupper(bases),
      count      = df$count[i]
    )
    if ("sample_id" %in% colnames(df)) out$sample_id <- df$sample_id[i]
    out
  }) |> bind_rows()
}

# Summarise per-position nucleotide frequencies from an expanded data frame.
position_freqs <- function(expanded) {
  expanded |>
    group_by(across(any_of(c("sample_id", "position", "nucleotide")))) |>
    summarise(count = sum(count), .groups = "drop") |>
    group_by(across(any_of(c("sample_id", "position")))) |>
    mutate(freq = count / sum(count)) |>
    ungroup()
}

# ── load counts ───────────────────────────────────────────────────────────────
meta <- read_tsv(samples_f, comment = "#", col_types = cols(.default = "c"))

tsv_files <- list.files(counts_dir, pattern = "_counts\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0) stop("No *_counts.tsv files found in: ", counts_dir)

counts_raw <- lapply(tsv_files, function(f) {
  sid <- sub("_counts\\.tsv$", "", basename(f))
  df  <- read_tsv(f, col_types = cols(type = "c", count = "d"))
  df$sample_id <- sid
  df
}) |> bind_rows()

# ── filter: soft-clip types only, exclude 1Sg ────────────────────────────────
non1sg <- counts_raw |>
  filter(grepl("^[0-9]+S[a-z]+$", type), type != "1Sg")

if (nrow(non1sg) == 0) {
  message("No non-1Sg soft-clip types found — nothing to plot.")
  quit(save = "no", status = 0)
}

# ── pooled composition (summed across all samples) ───────────────────────────
pooled_expanded <- non1sg |>
  group_by(type) |>
  summarise(count = sum(count), .groups = "drop") |>
  expand_positions()

pooled_freq <- position_freqs(pooled_expanded)

# Reference row: treat literature as an extra pseudo-position "Lit."
rt_ref_panel <- rt_prefs_df |>
  mutate(position_label = "Literature\n(RT ref.)", count = NA_real_)

# Build combined data for plot 1
pooled_plot_df <- pooled_freq |>
  mutate(position_label = paste0("Position ", position)) |>
  bind_rows(
    rt_ref_panel |> rename(freq_col = freq) |>
      transmute(
        position   = NA_integer_,
        nucleotide,
        count      = 0,
        freq       = freq_col,
        position_label
      )
  ) |>
  mutate(
    position_label = factor(
      position_label,
      levels = c(paste0("Position ", sort(unique(pooled_freq$position))),
                 "Literature\n(RT ref.)")
    ),
    nucleotide = factor(nucleotide, levels = c("A", "C", "G", "T"))
  )

p_pooled <- ggplot(pooled_plot_df,
                   aes(x = position_label, y = freq, fill = nucleotide)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = length(unique(pooled_freq$position)) + 0.5,
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = NUC_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title    = "Nucleotide composition of non-1Sg soft-clipped bases",
    subtitle = "Pooled across all samples | dashed line separates data from literature reference",
    x        = NULL,
    y        = "Frequency",
    fill     = "Nucleotide"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank())

ggsave(file.path(outdir, "03_softclip_nuc_composition.pdf"),
       p_pooled, width = max(5, length(unique(pooled_plot_df$position_label)) * 1.2 + 2),
       height = 4)
message("Saved: 03_softclip_nuc_composition.pdf")

# ── per-sample composition ────────────────────────────────────────────────────
sample_expanded <- non1sg |> expand_positions()
sample_freq     <- position_freqs(sample_expanded) |>
  mutate(
    position_label = paste0("Position ", position),
    nucleotide     = factor(nucleotide, levels = c("A", "C", "G", "T"))
  )

n_samples <- length(unique(sample_freq$sample_id))

p_by_sample <- ggplot(sample_freq,
                      aes(x = position_label, y = freq, fill = nucleotide)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = NUC_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap(~ sample_id, ncol = 3) +
  labs(
    title = "Non-1Sg soft-clip nucleotide composition per sample",
    x     = NULL,
    y     = "Frequency",
    fill  = "Nucleotide"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "03_softclip_nuc_by_sample.pdf"),
       p_by_sample,
       width  = min(10, n_samples * 2.5 + 2),
       height = ceiling(n_samples / 3) * 3 + 2)
message("Saved: 03_softclip_nuc_by_sample.pdf")

# ── write data TSV ────────────────────────────────────────────────────────────
pooled_out <- pooled_freq |>
  select(position, nucleotide, count, freq) |>
  mutate(source = "data")

rt_out <- rt_prefs_df |>
  mutate(position = NA_integer_, count = NA_real_, source = "literature_rt") |>
  select(position, nucleotide, count, freq, source)

write_tsv(bind_rows(pooled_out, rt_out),
          file.path(outdir, "03_softclip_nuc_composition.tsv"))
message("Saved: 03_softclip_nuc_composition.tsv")

# ── session info ──────────────────────────────────────────────────────────────
message("\n--- sessionInfo ---")
sessionInfo()
