#!/usr/bin/env bash
# NAD branch step 12: replicate concordance of NAD captureSeq libraries.
#
# Runs on vavilov through sbatch, in the r_nad env of step 5 (from project root):
#   sbatch --parsable -p main -J nad-repqc -c 2 --mem=8G -t 2:00:00 \
#          -o ~/logs/%x-%j.out --export=ALL \
#          --wrap "bash nad/scripts/12_replicate_qc/run.sh --self-test && bash nad/scripts/12_replicate_qc/run.sh"
# Spec: docs/specs/2026-11-16-nad-step12-replicate-qc.md
#
# Outputs: nad/results/12_replicate_qc/{read_classes,pairs,enrich_pairs,pca,cooks,loo,verdict}.tsv,
#          drivers_<strain>.tsv; nad/figures/12_{pairs_<strain>,distance}.pdf
#          (12_pca.pdf and 12_loo.pdf: Rscript nad/scripts/12_replicate_qc/12b_figures.R on the Mac)
#          --self-test: the same tables for a synthetic strain in nad/results/12_replicate_qc/selftest/
# Idempotent: skips if verdict.tsv exists, unless --force (the self-test always runs).
# RSCRIPT in the environment overrides the Rscript of the r_nad env.
#
# Usage: bash nad/scripts/12_replicate_qc/run.sh [--params config/params.yaml]
#        [--manifest ...] [--matrix ...] [--step5-dir ...] [--outdir ...] [--figdir ...]
#        [--self-test] [--force]

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
RES_DIR="$NAD_DIR/results/12_replicate_qc"
FIG_DIR="$NAD_DIR/figures"
FORCE=0; SELF=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)    PARAMS="$2"; shift 2 ;;
        --manifest)  MANIFEST="$2"; shift 2 ;;
        --matrix)    MATRIX="$2"; shift 2 ;;
        --step5-dir) STEP5_DIR="$2"; shift 2 ;;
        --outdir)    RES_DIR="$2"; shift 2 ;;
        --figdir)    FIG_DIR="$2"; shift 2 ;;
        --self-test) SELF="--self-test"; shift ;;
        --force)     FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -s "$MATRIX" ]] || die "missing matrix: $MATRIX (run step 04_ctss first)"
for st in $(awk -F'\t' 'NR>1{print $2}' "$MANIFEST" | sort -u); do
    [[ -s "$STEP5_DIR/${st}_S_vs_N.tsv.gz" ]] || die "missing $STEP5_DIR/${st}_S_vs_N.tsv.gz (run step 05_deseq first)"
done
mkdir -p "$RES_DIR" "$FIG_DIR" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/12_replicate_qc${SELF:+_selftest}_$(date +%Y%m%d_%H%M%S).log") 2>&1

if [[ -z "$SELF" && -s "$RES_DIR/verdict.tsv" && $FORCE -eq 0 ]]; then
    log_info "[skip] $RES_DIR/verdict.tsv exists (use --force)"
    column -t "$RES_DIR/verdict.tsv"; exit 0
fi

if [[ -z "${RSCRIPT:-}" ]]; then
    ENV="$(yaml_get conda_env_r_nad "$PARAMS")"
    [[ -n "$ENV" ]] || die "conda_env_r_nad is not set in $PARAMS"
    RSCRIPT="$(conda run -n "$ENV" which Rscript)"
fi
log_info "Rscript: $RSCRIPT"

"$RSCRIPT" "$SCRIPT_DIR/12_replicate_qc.R" \
    --matrix "$MATRIX" --manifest "$MANIFEST" --step5-dir "$STEP5_DIR" \
    --outdir "$RES_DIR" --figdir "$FIG_DIR" $SELF

[[ -n "$SELF" ]] && exit 0
log_info "verdict:"
column -t "$RES_DIR/verdict.tsv"
OUT="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="outlier") c=i; next} $c=="TRUE"{print $2}' "$RES_DIR/verdict.tsv")"
if [[ -n "$OUT" ]]; then
    log_warn "outlier run(s): $OUT. Step 5 is unchanged. To refit without them:"
    EXCL="$NAD_DIR/scripts/00_download/priority1_runs_no_outliers.tsv"
    echo "  grep -v -F -e $(echo $OUT | sed 's/ / -e /g') $MANIFEST > $EXCL"
    echo "  bash nad/scripts/05_deseq/run.sh --manifest $EXCL --outdir nad/results/05_deseq_no_outliers --figdir nad/figures/05_no_outliers"
else
    log_info "no outlier run under the step-12 rule"
fi
