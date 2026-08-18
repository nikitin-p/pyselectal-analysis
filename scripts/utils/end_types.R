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
END_TYPES <- c("1Sg", "2Sgg", "3Sggg", "2Stg",
               "1Sa", "2Sga",
               "1St", "1Sc",
               "M", "other")

END_TYPES_REAL <- setdiff(END_TYPES, "other")

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
# change nobody asked for. Each simply gains its two new entries, hued next to
# the sibling the new type belongs with.
END_TYPE_COLORS_BRIGHT <- c(
  "1Sg"   = "#4361ee", "2Sgg" = "#f72585", "3Sggg" = "#7209b7",
  "2Stg"  = "#4895ef",
  "1Sa"   = "#f4a261", "2Sga" = "#e76f51",
  "1St"   = "#90be6d", "1Sc"  = "#b5838d",
  "M"     = "#06d6a0", "other" = "#adb5bd"
)

END_TYPE_COLORS_TAB10 <- c(
  "1Sg"   = "#1f77b4", "2Sgg" = "#ff7f0e", "3Sggg" = "#2ca02c",
  "2Stg"  = "#17becf",
  "1Sa"   = "#d62728", "2Sga" = "#e377c2",
  "1St"   = "#9467bd", "1Sc"  = "#8c564b",
  "M"     = "#7f7f7f", "other" = "#bcbd22"
)

# ── column-name helpers ───────────────────────────────────────────────────────
end_type_cols <- function(prefix, types = END_TYPES) paste0(prefix, types)
