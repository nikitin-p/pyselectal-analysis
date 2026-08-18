#!/usr/bin/env bash
# Stage 14 (human): render the HTML analysis report for the human LUHMES dataset.
#
# Renders reports/cage_5prime_report_human.qmd into a self-contained HTML
# document that regenerates every figure from the summary tables written by the
# earlier stages. No BAM file is read, so this stage is fast and can be re-run at
# any time.
#
# The report is deliberately tolerant of missing inputs: each section checks for
# its own tables and reports their absence rather than failing the render, so it
# can be rendered against a partial run and re-rendered as stages complete.
#
# Usage:
#   bash scripts/07_report/run_human.sh [--params FILE] [--samples FILE] \
#       [--compare DIR] [--compare_label NAME] [--out FILE] [--force]
#
# Defaults:
#   --params        config/params_human.yaml
#   --samples       config/samples_human.tsv
#   --compare       the yeast results tree, for the cross-dataset section
#   --compare_label yeast

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params_human.yaml"
SAMPLES="config/samples_human.tsv"
COMPARE=""
COMPARE_LABEL="yeast"
HTML=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)        PARAMS="$2";        shift 2 ;;
        --samples)       SAMPLES="$2";       shift 2 ;;
        --compare)       COMPARE="$2";       shift 2 ;;
        --compare_label) COMPARE_LABEL="$2"; shift 2 ;;
        --out)           HTML="$2";          shift 2 ;;
        --force)         FORCE=1;            shift   ;;
        *) die "Unknown argument: $1" ;;
    esac
done

QMD="$ROOT_DIR/reports/cage_5prime_report_human.qmd"
[[ -f "$QMD" ]] || die "Report source not found: $QMD"
[[ -f "$PARAMS" ]] || die "Params file not found: $PARAMS"

# The results tree comes from the params file, so that the report always reads
# the same tree the pipeline wrote.
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"
RESULTS_REL="${RESULTS_ROOT#$ROOT_DIR/}"

# The comparison tree defaults to the yeast run, which lives under the default
# results root of config/params.yaml rather than under this dataset's root.
if [[ -z "$COMPARE" ]]; then
    COMPARE="$(resolve_results_root "$ROOT_DIR/config/params.yaml" "$ROOT_DIR")"
fi
COMPARE_REL="${COMPARE#$ROOT_DIR/}"

[[ -z "$HTML" ]] && HTML="$ROOT_DIR/reports/cage_5prime_report_human.html"

log_info "Stage 14 (human): HTML analysis report"
log_info "  results  : $RESULTS_REL"
log_info "  samples  : $SAMPLES"
log_info "  compare  : $COMPARE_REL  (labelled '$COMPARE_LABEL')"

if ! command -v quarto >/dev/null 2>&1; then
    log_error "quarto not found on PATH — install it from https://quarto.org before running this stage."
    exit 1
fi

if [[ ! -d "$RESULTS_ROOT" ]]; then
    log_warn "results tree does not exist yet: $RESULTS_REL"
    log_warn "the report will render, but every section will report its inputs as missing"
fi

if [[ -f "$HTML" && $FORCE -eq 0 ]]; then
    log_info "  output exists, skipping (use --force to re-render): $HTML"
    exit 0
fi

# quarto is called without --output: renaming the output while embed-resources is
# on makes its post-processing step look for a supporting directory under the new
# name and fail. The default name is used instead, and the file is moved
# afterwards when a different destination was asked for.
DEFAULT_HTML="$ROOT_DIR/reports/$(basename "${QMD%.qmd}").html"

quarto render "$QMD" \
    --to html \
    -P root:.. \
    -P "results:$RESULTS_REL" \
    -P "samples:$SAMPLES" \
    -P "compare:$COMPARE_REL" \
    -P "compare_label:$COMPARE_LABEL"

if [[ "$DEFAULT_HTML" != "$HTML" ]]; then
    mv "$DEFAULT_HTML" "$HTML"
fi

log_info "Stage 14 (human) complete."
log_info "  report → ${HTML#$ROOT_DIR/}"
