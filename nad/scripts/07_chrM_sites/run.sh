#!/usr/bin/env bash
# NAD branch step 7 (cluster part): DESeq2 S vs N at the chrM ATATAAGTA starts.
#
# Runs on vavilov through sbatch, in the r_nad env of step 5:
#   sbatch --parsable -p main -J nad-chrM -c 2 --mem=8G -t 1:00:00 \
#          -o ~/logs/%x-%j.out --export=ALL \
#          --wrap "bash nad/scripts/07_chrM_sites/run.sh"
# Then copy nad/results/07_chrM_sites/*.tsv to the Mac and run
# 07_join_cage.R there (see its header).
# Spec: docs/specs/2026-11-04-nad-step7-chrM-promoter-sites.md
#
# Outputs: nad/results/07_chrM_sites/{motif_sites,controls,sanity}.tsv
# Idempotent: skips if motif_sites.tsv exists, unless --force.
#
# Usage: bash nad/scripts/07_chrM_sites/run.sh [--params config/params.yaml]
#        [--manifest ...] [--matrix ...] [--step5-dir ...] [--step5-log ...] [--outdir ...] [--dry-run] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
MATRIX="$NAD_DIR/results/04_ctss/matrix_reads.tsv.gz"
STEP5_DIR="$NAD_DIR/results/05_deseq"
RES_DIR="$NAD_DIR/results/07_chrM_sites"
FORCE=0; DRY=""; STEP5_LOG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)    PARAMS="$2"; shift 2 ;;
        --manifest)  MANIFEST="$2"; shift 2 ;;
        --matrix)    MATRIX="$2"; shift 2 ;;
        --step5-dir) STEP5_DIR="$2"; shift 2 ;;
        --step5-log) STEP5_LOG="$2"; shift 2 ;;
        --outdir)    RES_DIR="$2"; shift 2 ;;
        --dry-run)   DRY="--dry-run"; shift ;;
        --force)     FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

FASTA="$(yaml_get genome_fasta "$PARAMS")"; [[ "$FASTA" != /* ]] && FASTA="$ROOT_DIR/$FASTA"
MASK="$STEP5_DIR/ncrna_mask.bed"
[[ -s "$MATRIX" ]] || die "missing matrix: $MATRIX (run step 04_ctss first)"
[[ -s "$MASK"   ]] || die "missing mask: $MASK (run step 05_deseq first)"
[[ -s "$FASTA"  ]] || die "genome FASTA not found: $FASTA"
mkdir -p "$RES_DIR" "$NAD_DIR/logs"

# The newest step-5 log that printed size factors, for the sanity comparison.
[[ -n "$STEP5_LOG" ]] || STEP5_LOG="$(grep -l "size factors:" "$NAD_DIR"/logs/05_deseq_*.log 2>/dev/null | sort | tail -n 1 || true)"

exec > >(tee -a "$NAD_DIR/logs/07_chrM_sites_$(date +%Y%m%d_%H%M%S).log") 2>&1

if [[ -s "$RES_DIR/motif_sites.tsv" && $FORCE -eq 0 && -z "$DRY" ]]; then
    log_info "[skip] $RES_DIR/motif_sites.tsv exists (use --force)"
    column -t "$RES_DIR/sanity.tsv"; exit 0
fi

if [[ -n "${RSCRIPT:-}" ]]; then
    log_info "Rscript from RSCRIPT: $RSCRIPT"
else
    ENV="$(yaml_get conda_env_r_nad "$PARAMS")"
    [[ -n "$ENV" ]] || die "conda_env_r_nad is not set in $PARAMS"
    RSCRIPT="$(conda run -n "$ENV" which Rscript)"
    log_info "Rscript: $RSCRIPT"
fi
[[ -n "$STEP5_LOG" ]] && log_info "step-5 log for size factors: $STEP5_LOG" \
                      || log_warn "no step-5 log with size factors found; size factors will not be compared"

"$RSCRIPT" "$SCRIPT_DIR/deseq_chrM_sites.R" \
    --matrix "$MATRIX" --manifest "$MANIFEST" --mask "$MASK" --fasta "$FASTA" \
    --step5-dir "$STEP5_DIR" ${STEP5_LOG:+--step5-log "$STEP5_LOG"} \
    --outdir "$RES_DIR" $DRY
