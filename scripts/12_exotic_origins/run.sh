#!/usr/bin/env bash
# Stage 19: origins of the exotic 5' end types.
#
# "Exotic" means every fine 5'-end type other than 1Sg, so the unclipped class M
# is included. The stage asks what kind of start site each type's signal comes
# from, in two currencies at once (alignments and distinct CTSS positions), and
# how productive the positions each type dominates are. The categories are
# relative to the type of the bar being drawn (below the depth floor,
# 1Sg-dominated, dominated by this type itself, or dominated by neither), which
# lets 1Sg be drawn on the same axis as the reference the exotics are read
# against.
#
# Steps:
#   1. Re-run the Stage 17 dominance cascade with the CTSS position as the unit,
#      at this stage's own depth floor (exo_min_depth)
#   2. Tally every type by category and draw the composite figure
#   3. Build one genomic-context sequence motif per dominating type
#   4. Annotate the dominated positions and draw their genomic localization
#
# Step 3 reads the genome FASTA named by the sample sheet and step 4 reads the
# annotation GTF; those are the only references this stage opens, and no BAM is
# read at any step. Step 4 also reuses Stage 10's enhancer BED when it exists,
# because detecting enhancers would require the BAMs.
#
# This stage reads Stage 17's tables only. Run scripts/10_dominant_endtype/run.sh
# first; no BAM is opened here and nothing is re-typed.
#
# Usage:
#   bash scripts/12_exotic_origins/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2"; shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --force)  FORCE_FLAG="--force"; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# The r_cage environment may be absent; the system R library carries every
# package this stage needs, so fall back to it rather than failing.
CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)
    log_info "using Rscript from conda env: $CONDA_ENV"
else
    require_cmd Rscript
    RSCRIPT=$(command -v Rscript)
    log_warn "conda env '$CONDA_ENV' not found; using system Rscript"
fi

# Relative paths in the params file are resolved against the project root;
# absolute paths are used as given.
resolve_path() {
    local v="$1"
    [[ -z "$v" ]] && { printf ''; return; }
    [[ "$v" != /* ]] && v="$ROOT_DIR/$v"
    printf '%s' "$v"
}
ANNOTATION_GTF="$(resolve_path "$(yaml_get annotation_gtf "$PARAMS")")"
LOC_DIR="$ROOT_DIR/$(yaml_get results_root "$PARAMS")/localization"

DOM_DIR="$ROOT_DIR/$(yaml_get results_dominance "$PARAMS")"
EXO_DIR="$ROOT_DIR/$(yaml_get results_exotic_origins "$PARAMS")"
FIG_DIR="$ROOT_DIR/$(yaml_get results_figures "$PARAMS")"

LOG="$ROOT_DIR/logs/19_exotic_origins_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 19: origins of the exotic 5' end types"

log_info "  step 1: position-level dominance cascade"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_position_dominance.R" \
    --ctss       "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --background "$DOM_DIR/background_shares.tsv"  \
    --outdir     "$EXO_DIR" \
    --params     "$PARAMS" $FORCE_FLAG

log_info "  step 2: tallies and composite figure"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_origins_figure.R" \
    --dominance      "$EXO_DIR/position_dominance.tsv.gz" \
    --ctss           "$DOM_DIR/fine_typed_ctss.tsv.gz"    \
    --cluster_counts "$DOM_DIR/cluster_type_counts.tsv"   \
    --outdir         "$EXO_DIR" \
    --figdir         "$FIG_DIR" \
    --params         "$PARAMS" $FORCE_FLAG

log_info "  step 3: genomic context motifs per dominating type"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/03_context_motifs.R" \
    --dominance "$EXO_DIR/position_dominance.tsv.gz" \
    --samples   "$SAMPLES" \
    --outdir    "$EXO_DIR" \
    --figdir    "$FIG_DIR" \
    --params    "$PARAMS" $FORCE_FLAG

log_info "  step 4: genomic localization of the dominated positions"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/04_dominated_localization.R" \
    --dominance "$EXO_DIR/position_dominance.tsv.gz" \
    --gtf       "$ANNOTATION_GTF" \
    --enhancers "$LOC_DIR/enhancer_regions.bed" \
    --outdir    "$EXO_DIR" \
    --figdir    "$FIG_DIR" \
    --params    "$PARAMS" $FORCE_FLAG

log_info "Stage 19 complete."
log_info "  tables  → $EXO_DIR"
log_info "  figures → $FIG_DIR/20_exo_origins.pdf"
log_info "            $FIG_DIR/21_exo_context_logos.pdf"
log_info "            $FIG_DIR/21_exo_context_delta.pdf"
log_info "            $FIG_DIR/22_exo_localization.pdf"
