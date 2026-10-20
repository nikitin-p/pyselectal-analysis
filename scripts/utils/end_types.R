# The 5' end-type vocabulary, and everything keyed by it.
#
# This file is the SOLE definition. Every stage sources it rather than
# hardcoding a type vector, which is what makes adding a type a one-line change
# and makes coverage checkable: a literal "2Sgg" or "3Sggg" appearing anywhere
# under scripts/ outside this file is by definition a site that was missed.
# test/run_end_types_test.sh enforces that.
#
# Source it from a script with:
#   .f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
#   source(if (length(.f))
#            file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "end_types.R")
#          else "scripts/utils/end_types.R")

# ── the vocabulary ────────────────────────────────────────────────────────────
# Canonical order, mechanistic rather than alphabetical: each type sits beside
# the sibling it shares a genome-adjacent base with. 2Stg is a 1Sg carrying one
# extra 5' T and so joins the G-ladder block; 2Sga abuts the genome with A and
# so sits with 1Sa. "other" is the residual and always sorts last.
#
# The set is every fine type above 0.1 % of reads pooled over the eight yeast
# libraries (results/five_prime/*_counts.tsv), which puts 3Sggg (0.08 %) last.
# The three-base clips ending in TG (3Sgtg, 3Sttg, 3Sctg) are 2Stg carrying one
# further 5' base and so follow it; 2Sag abuts with G and closes the G block.
# 2Saa and 3Sgaa abut with A and join 1Sa, 3Sgat abuts with T and joins 1St.
END_TYPES <- c("1Sg", "2Sgg", "3Sggg", "2Stg", "3Sgtg", "3Sttg", "3Sctg", "2Sag",
               "1Sa", "2Sga", "2Saa", "3Sgaa",
               "1St", "3Sgat",
               "1Sc",
               "M", "other")

END_TYPES_REAL <- setdiff(END_TYPES, "other")

# The canonical class: the single untemplated G that is the CAGE / oligo-capping
# signature. It is named here rather than at each use site, so that a stage
# defining "exotic" as "not the canonical type" reads it from the vocabulary
# instead of hardcoding the label.
END_TYPE_CANONICAL <- END_TYPES[1]

# ── the coarse subset drawn by Stage 3 ────────────────────────────────────────
# Stage 3 works from pyselectal's own --count output and relabels it into the
# broad categories of results/figures/03_replicate_heatmap.pdf, which resolve
# the G-appending ladder, its two two-base siblings and M by name and sweep
# every single-base non-G clip into "other_soft_clipped". This vector is that
# same set expressed in the fine vocabulary, so a figure built on the fine
# types can be restricted to exactly the types Stage 3 names.
#
# It is deliberately NOT a rule over END_TYPES (there is no clean predicate
# that keeps 2Sga while dropping 1Sa): it is the historical Stage 3 category
# list, and it is defined here rather than at a use site so that the fine-type
# literals stay confined to this file.
END_TYPES_COARSE <- c("1Sg", "2Sgg", "3Sggg", "2Stg", "2Sga", "M")

# ── the genome-adjacent clip base ─────────────────────────────────────────────
# A fine-type label is <clip length>S<clip sequence>, with the clip written in
# the read's original 5'->3' orientation (see clip5_seq() in
# scripts/10_dominant_endtype/dom_helpers.R). The base abutting the alignment is
# therefore the LAST character of the label, not the first.
#
# Under STAR --alignEndsType Local a base is soft-clipped only where it
# mismatches the genome, so that abutting base can never equal the genomic base
# at -1, which is what makes a clipped type impossible on a dinucleotide
# starting with its own clip base. M carries no clip and is unconstrained.
#
# This was verified against results/bam_subsampled/yeast_cage_ana_r1_sub.bam and
# data/sacCer3.fa over 106,874 two-base-clipped reads: the clip GA occurs on a
# genomic G at -1 in 4,865 reads and on a genomic A in zero, and the clip TG
# occurs on a genomic T in 11,631 reads and on a genomic G in zero.
.last_char <- function(x) toupper(substr(x, nchar(x), nchar(x)))
END_TYPE_CLIP_BASE <- setNames(.last_char(END_TYPES_REAL), END_TYPES_REAL)
END_TYPE_CLIP_BASE["M"] <- NA_character_

# ── palettes ──────────────────────────────────────────────────────────────────
# Two palettes are kept, deliberately. They are the two that already existed in
# the pipeline, and merging them would restyle every current figure, which is a
# change nobody asked for. Each new type is hued next to the sibling it belongs
# with.
END_TYPE_COLORS_BRIGHT <- c(
  "1Sg"   = "#4361ee", "2Sgg" = "#f72585", "3Sggg" = "#7209b7",
  "2Stg"  = "#4895ef", "3Sgtg" = "#4cc9f0", "3Sttg" = "#3a0ca3",
  "3Sctg" = "#560bad", "2Sag"  = "#b5179e",
  "1Sa"   = "#f4a261", "2Sga" = "#e76f51", "2Saa" = "#e9c46a",
  "3Sgaa" = "#bc6c25",
  "1St"   = "#90be6d", "3Sgat" = "#2a9d8f", "1Sc"  = "#b5838d",
  "M"     = "#06d6a0", "other" = "#adb5bd"
)

END_TYPE_COLORS_TAB10 <- c(
  "1Sg"   = "#1f77b4", "2Sgg" = "#ff7f0e", "3Sggg" = "#2ca02c",
  "2Stg"  = "#17becf", "3Sgtg" = "#9edae5", "3Sttg" = "#aec7e8",
  "3Sctg" = "#393b79", "2Sag"  = "#98df8a",
  "1Sa"   = "#d62728", "2Sga" = "#e377c2", "2Saa" = "#ff9896",
  "3Sgaa" = "#ffbb78",
  "1St"   = "#9467bd", "3Sgat" = "#c5b0d5", "1Sc"  = "#8c564b",
  "M"     = "#7f7f7f", "other" = "#bcbd22"
)

# ── column-name helpers ───────────────────────────────────────────────────────
end_type_cols <- function(prefix, types = END_TYPES) paste0(prefix, types)
