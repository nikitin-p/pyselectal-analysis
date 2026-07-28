#!/usr/bin/env bash
# Run motif analysis comparing 1Sg vs 2Sgg/3Sggg TSS
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/utils/log.sh

FORCE=""
[[ "${1:-}" == "--force" ]] && FORCE="--force"

log_info "Running motif analysis..."

Rscript scripts/05_tss_end_types/06_motif_analysis.R \
    --tss_matrix results/tss/tss_matrix.tsv \
    --params     config/params.yaml \
    --gtf        data/sacCer3.ensGene.gtf \
    --outdir     results/motifs \
    --figdir     results/figures \
    --min_reads  10 \
    --top_pct    10 \
    --upstream   50 \
    --downstream 10 \
    $FORCE

log_info "Motif analysis complete"
