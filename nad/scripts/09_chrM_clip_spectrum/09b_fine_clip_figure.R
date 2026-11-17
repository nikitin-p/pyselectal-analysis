#!/usr/bin/env Rscript
# NAD branch step 9b: is the NAD-cap signature on CAGE reads 1Sa, or a
# two-base clip ending in A (2Saa, 2Sca, 2Sga, 2Sta)?
#
# Uses the uncollapsed fine clip types of step 9 (fine_clip_types.py), so that
# 2Sta, 2Sca and the three-base clips stay visible. Three-base clips ending in
# TA (or TG) are pooled as 3S·ta (3S·tg).
#   A  share in chrM / share in nuclear for each A-ending clip type;
#   B  clip types at the chrM promoter starts 58008 and 73704 (NAD) and 6548
#      (not an NAD start, an A6 tract);
#   C  the extra-T ladder 1S → 2St → 3S·t for nuclear 1Sg (m7G) and for 1Sa at
#      58008 and 73704;
#   D  A-ending clip shares at the nuclear NAD starts of step 6 run A (triple
#      mutant) against all nuclear reads.
#
# Outputs: <figdir>/09_nad_fine_clip_types.{pdf,png};
#          <outdir>/fine_clip_{enrichment,ladder,nad_starts}.tsv (panels A, C, D).
#
# Usage (from project root):
#   Rscript nad/scripts/09_chrM_clip_spectrum/09b_fine_clip_figure.R [--force]
#     [--fine nad/results/09_chrM_clip_spectrum/fine_types.tsv.gz]
#     [--step6 nad/results/06_cage_end_types/nad_positions_cage.tsv.gz]
#     [--outdir nad/results/09_chrM_clip_spectrum] [--figdir nad/figures]
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
fine_f  <- get_arg("--fine",   "nad/results/09_chrM_clip_spectrum/fine_types.tsv.gz")
step6_f <- get_arg("--step6",  "nad/results/06_cage_end_types/nad_positions_cage.tsv.gz")
outdir  <- get_arg("--outdir", "nad/results/09_chrM_clip_spectrum")
figdir  <- get_arg("--figdir", "nad/figures")
force   <- "--force" %in% args
out <- file.path(figdir, "09_nad_fine_clip_types")
if (file.exists(paste0(out, ".pdf")) && !force) { message("[skip] ", out, ".pdf exists (use --force)"); quit(save = "no") }
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

f <- read.delim(gzfile(fine_f), colClasses = c("character","integer","character","character","integer"))
f$fam <- ifelse(grepl("^3S[acgt]ta$", f$type), "3S·ta",
         ifelse(grepl("^3S[acgt]tg$", f$type), "3S·tg", f$type))
f$comp <- ifelse(f$chr == "chrM", "chrM", "nuclear")
TYPES <- c("1Sa", "2Sta", "3S·ta", "2Saa", "2Sca", "2Sga")
red <- "#e63946"; grey <- "grey55"
th <- theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", size = 10),
                                           plot.subtitle = element_text(size = 8, colour = "grey30"))

# A: chrM / nuclear enrichment
tot <- f %>% group_by(comp) %>% summarise(T = sum(n))
sh <- function(d) d %>% group_by(comp, fam) %>% summarise(n = sum(n), .groups = "drop") %>%
  left_join(tot, by = "comp") %>% mutate(p = n / T)
a <- sh(f) %>% filter(fam %in% TYPES)
a2 <- sh(filter(f, !(chr == "chrM" & pos == 6548 & strand == "+"))) %>% filter(fam == "2Sta") %>%
  mutate(fam = "2Sta\n(no 6548)")
A <- bind_rows(a, a2) %>% select(comp, fam, p) %>% tidyr::pivot_wider(names_from = comp, values_from = p) %>%
  mutate(ratio = chrM / nuclear, fam = factor(fam, c(TYPES[1:2], "2Sta\n(no 6548)", TYPES[3:6])),
         sig = fam %in% c("1Sa", "2Sta", "2Sta\n(no 6548)", "3S·ta"))
pA <- ggplot(A, aes(fam, ratio, fill = sig)) + geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_text(aes(label = ifelse(ratio < 1, sprintf("×%.2f", ratio), sprintf("×%.1f", ratio)),
                vjust = ifelse(ratio < 1, 1.3, -0.3)), size = 2.8) +
  scale_y_log10(limits = c(0.02, 200), expand = expansion(mult = c(0, 0.05))) + scale_fill_manual(values = c(`TRUE` = red, `FALSE` = grey), guide = "none") +
  labs(x = NULL, y = "share in chrM / share in nuclear", title = "A  Which A-ending clips are enriched on chrM?",
       subtitle = "Only 1Sa and its T-extended relatives 2Sta, 3S·ta; 2Saa, 2Sca, 2Sga are depleted") + th

# B: per-site counts
sites <- data.frame(pos = c(58008, 73704, 6548),
                    lab = c("58008 (21S)\nNAD start", "73704 (COX2)\nNAD start", "6548 (15S)\nnot an NAD start"))
B <- f %>% filter(chr == "chrM", strand == "+", pos %in% sites$pos, fam %in% TYPES) %>%
  group_by(pos, fam) %>% summarise(n = sum(n), .groups = "drop") %>%
  tidyr::complete(pos, fam = TYPES, fill = list(n = 0)) %>% left_join(sites, by = "pos") %>%
  mutate(fam = factor(fam, rev(TYPES)), lab = factor(lab, sites$lab),
         sig = fam %in% c("1Sa", "2Sta", "3S·ta") & pos != 6548)
pB <- ggplot(B, aes(n + 1, fam, fill = sig)) + geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.6) + facet_wrap(~lab, nrow = 1) +
  scale_x_log10(expand = expansion(mult = c(0, 0.3))) + scale_fill_manual(values = c(`TRUE` = red, `FALSE` = grey), guide = "none") +
  labs(x = "CAGE reads + 1 (7 libraries)", y = NULL, title = "B  Clip types at chrM promoter starts",
       subtitle = "At NAD starts 1Sa dominates; at 6548 (A6 tract, GT|AAAAAA) the clips come from an A-run shift, not a cap") + th

# C: the T ladder, m7G vs NAD
lad <- function(d, t1, t2, t3) { s <- sapply(c(t1, t2, t3), function(t) sum(d$n[d$fam == t]))
  c(s[2] / s[1], s[3] / s[2]) }
nuc <- filter(f, comp == "nuclear")
C <- rbind(
  data.frame(grp = "nuclear m7G\n1Sg → 2Stg → 3S·tg", step = c("2nd / 1st", "3rd / 2nd"), r = lad(nuc, "1Sg", "2Stg", "3S·tg")),
  data.frame(grp = "chrM 58008\n1Sa → 2Sta → 3S·ta", step = c("2nd / 1st", "3rd / 2nd"),
             r = lad(filter(f, chr == "chrM", pos == 58008, strand == "+"), "1Sa", "2Sta", "3S·ta")),
  data.frame(grp = "chrM 73704\n1Sa → 2Sta → 3S·ta", step = c("2nd / 1st", "3rd / 2nd"),
             r = lad(filter(f, chr == "chrM", pos == 73704, strand == "+"), "1Sa", "2Sta", "3S·ta")))
C$grp <- factor(C$grp, unique(C$grp))
pC <- ggplot(C, aes(grp, r, fill = step)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", r)), position = position_dodge(0.75), vjust = -0.3, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + scale_fill_manual(values = c("grey35", "grey70"), name = NULL) +
  labs(x = NULL, y = "read ratio", title = "C  1Sa carries the same extra-T ladder as 1Sg",
       subtitle = "Each step adds one 5′ base (T, then any); the ladder is present at both NAD starts, 58008 matches m7G closely") +
  th + theme(legend.position = "top")

# D: nuclear NAD starts vs nuclear background
ra <- read.delim(gzfile(step6_f), colClasses = c(chr = "character", pos = "integer", strand = "character"))
ra <- unique(ra[ra$run == "A", c("chr", "pos", "strand")])
nadset <- semi_join(nuc, ra, by = c("chr", "pos", "strand"))
D <- bind_rows(
  nadset %>% mutate(set = sprintf("NAD-enriched starts (%s)", format(nrow(ra), big.mark = ","))),
  nuc %>% mutate(set = "all nuclear reads")) %>%
  group_by(set) %>% mutate(T = sum(n)) %>% group_by(set, fam) %>% summarise(p = 100 * sum(n) / T[1], .groups = "drop") %>%
  filter(fam %in% TYPES) %>% mutate(fam = factor(fam, TYPES))
pD <- ggplot(D, aes(fam, p, colour = set)) + geom_point(position = position_dodge(0.5), size = 3) +
  scale_y_log10() + scale_colour_manual(values = c(grey, red), name = NULL) +
  labs(x = NULL, y = "% of CAGE reads (log)", title = "D  Nuclear NAD starts (triple mutant, step 6 run A)",
       subtitle = "No A-ending clip type is above nuclear background") + th + theme(legend.position = "top")

p <- (pA | pB) / (pC | pD) + plot_layout(widths = c(1, 1.3))
ggsave(paste0(out, ".pdf"), p, width = 13, height = 8.5, device = cairo_pdf)
ggsave(paste0(out, ".png"), p, width = 13, height = 8.5, dpi = 150)
write.table(A, file.path(outdir, "fine_clip_enrichment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(C, file.path(outdir, "fine_clip_ladder.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(D, file.path(outdir, "fine_clip_nad_starts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(A)); print(C); print(as.data.frame(D))
sessionInfo()
