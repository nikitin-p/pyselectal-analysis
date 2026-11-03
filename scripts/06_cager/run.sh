#!/usr/bin/env bash
# Stage 10: CAGEr promoter clusters — normalize, distclu, IQW, annotation, plots.
#
# Steps:
#   1. Build CAGEexp from typed_ctss.tsv.gz → distclu + IQW → RDS + tagclusters TSV
#   2. Annotate consensus tag clusters against gene annotation BED
#   3. Diagnostic plots: correlation, reverse cumulatives, IQW, annotation enrichment
#   4. Browser tracks: per-library power-law TPM bigWigs (filtered CTSSs only),
#      per-sample tag clusters and consensus clusters
#
# Prerequisites:
#   - Stage 8 must have run (results/tss/typed_ctss.tsv.gz must exist)
#
# Usage:
#   bash scripts/06_cager/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
[[ "${1:-}" == "--force" ]] && FORCE_FLAG="--force"

# The r_cage environment may be absent; the system R library carries CAGEr, so
# fall back to it rather than failing.
CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    RSCRIPT="conda run -n $CONDA_ENV Rscript"
    log_info "using Rscript from conda env: $CONDA_ENV"
else
    require_cmd Rscript
    RSCRIPT="$(command -v Rscript)"
    log_warn "conda env '$CONDA_ENV' not found; using system Rscript"
fi
BEDTOBIGBED="${BEDTOBIGBED:-$(command -v bedToBigBed || true)}"
GENOME_FASTA="$(yaml_get genome_fasta "$PARAMS")"
[[ "$GENOME_FASTA" != /* ]] && GENOME_FASTA="$ROOT_DIR/$GENOME_FASTA"

CAGER_DIR="$(yaml_get results_cager "$PARAMS")"
CAGER_DIR="${CAGER_DIR:-results/cager}"
FIGURES_DIR="$(yaml_get results_figures "$PARAMS")"
FIGURES_DIR="${FIGURES_DIR:-results/figures}"

log_info "Stage 10: CAGEr promoter clusters"

log_info "  step 1: build CAGEexp + distclu + IQW"
$RSCRIPT "$SCRIPT_DIR/01_cageexp.R" \
    --ctss    "$ROOT_DIR/results/tss/typed_ctss.tsv.gz" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/$CAGER_DIR" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 2: annotate tag clusters"
$RSCRIPT "$SCRIPT_DIR/02_annotate.R" \
    --tagclusters "$ROOT_DIR/$CAGER_DIR/tagclusters.tsv.gz" \
    --samples     "$SAMPLES" \
    --outdir      "$ROOT_DIR/$CAGER_DIR" \
    --params      "$PARAMS" \
    $FORCE_FLAG

log_info "  step 3: diagnostic plots"
$RSCRIPT "$SCRIPT_DIR/03_plot.R" \
    --rds          "$ROOT_DIR/$CAGER_DIR/cageexp.rds" \
    --annot        "$ROOT_DIR/$CAGER_DIR/annotation_enrichment.tsv" \
    --tagclusters  "$ROOT_DIR/$CAGER_DIR/tagclusters_annotated.tsv" \
    --outdir       "$ROOT_DIR/$FIGURES_DIR" \
    --params       "$PARAMS" \
    $FORCE_FLAG

log_info "  step 4: browser tracks (TPM bigWigs, tag and consensus clusters)"
[[ -z "$BEDTOBIGBED" ]] && log_warn "bedToBigBed not found; only BED files are written (set BEDTOBIGBED=PATH)"
$RSCRIPT "$SCRIPT_DIR/04_export_tracks.R" \
    --rds    "$ROOT_DIR/$CAGER_DIR/cageexp.rds" \
    --fai    "$GENOME_FASTA.fai" \
    --outdir "$ROOT_DIR/$CAGER_DIR/tracks" \
    --params "$PARAMS" \
    ${BEDTOBIGBED:+--bedtobigbed "$BEDTOBIGBED"} \
    $FORCE_FLAG

log_info "Stage 10 complete."
log_info "  tracks     → $CAGER_DIR/tracks"
log_info "  cageexp    → $CAGER_DIR/cageexp.rds"
log_info "  clusters   → $CAGER_DIR/tagclusters_annotated.tsv"
log_info "  figures    → $FIGURES_DIR/10_*.pdf"
