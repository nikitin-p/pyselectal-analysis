#!/usr/bin/env bash
# Stage 07: Render the HTML analysis report
#
# Renders reports/cage_5prime_report.qmd into a self-contained HTML document that
# regenerates every figure from the summary tables written by the earlier stages.
# No BAM files are read, so this stage is fast and can be re-run at any time.
#
# Usage:
#   bash scripts/07_report/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"

QMD="$ROOT_DIR/reports/cage_5prime_report.qmd"
HTML="$ROOT_DIR/reports/cage_5prime_report.html"
FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

log_info "Stage 07: HTML analysis report"

if ! command -v quarto >/dev/null 2>&1; then
    log_error "quarto not found on PATH — install it from https://quarto.org before running this stage."
    exit 1
fi

if [[ -f "$HTML" && $FORCE -eq 0 ]]; then
    log_info "  output exists, skipping (use --force to re-render): $HTML"
    exit 0
fi

quarto render "$QMD"

log_info "Stage 07 complete."
log_info "  report → reports/cage_5prime_report.html"
