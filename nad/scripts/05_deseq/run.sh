#!/usr/bin/env bash
# NAD branch step 5: S-versus-N enrichment of NAD captureSeq 5' ends (DESeq2).
#
# Downloads the Ensembl R64-1-1 GTF once (data/sacCer3.ensGene.gtf has no
# tRNA/rRNA/snoRNA genes, which the mask needs), then runs deseq_s_vs_n.R on
# the raw-read matrix of step 4. Spec: docs/specs/2026-11-03-nad-step5-deseq2.md
#
# Outputs: nad/results/05_deseq/{ncrna_mask.bed,<strain>_S_vs_N.tsv.gz,
#          <strain>_enriched.tsv,summary.tsv,context_enriched.tsv}
#          nad/figures/05_<strain>_ma.pdf, nad/figures/05_first_base.pdf
# Idempotent: skips the download if the GTF exists and the fit if summary.tsv
# exists, unless --force.
#
# Usage: bash nad/scripts/05_deseq/run.sh [--params config/params.yaml]
#        [--manifest nad/scripts/00_download/priority1_runs.tsv]
#        [--matrix nad/results/04_ctss/matrix_reads.tsv.gz] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
MATRIX="$NAD_DIR/results/04_ctss/matrix_reads.tsv.gz"
RES_DIR="$NAD_DIR/results/05_deseq"
FIG_DIR="$NAD_DIR/figures"
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)   PARAMS="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --matrix)   MATRIX="$2"; shift 2 ;;
        --outdir)   RES_DIR="$2"; shift 2 ;;
        --figdir)   FIG_DIR="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

FASTA="$(yaml_get genome_fasta "$PARAMS")"; [[ "$FASTA" != /* ]] && FASTA="$ROOT_DIR/$FASTA"
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.113.gtf.gz"
GTF="$NAD_DIR/data/annotation/$(basename "$GTF_URL")"
[[ -s "$MATRIX" ]] || die "missing matrix: $MATRIX (run step 04_ctss first)"
[[ -s "$FASTA"  ]] || die "genome FASTA not found: $FASTA"
mkdir -p "$RES_DIR" "$FIG_DIR" "$NAD_DIR/logs" "$(dirname "$GTF")"
exec > >(tee -a "$NAD_DIR/logs/05_deseq_$(date +%Y%m%d_%H%M%S).log") 2>&1

if [[ -s "$RES_DIR/summary.tsv" && $FORCE -eq 0 ]]; then
    log_info "[skip] $RES_DIR/summary.tsv exists (use --force)"
    column -t "$RES_DIR/summary.tsv"; exit 0
fi

if [[ ! -s "$GTF" ]]; then
    log_info "downloading $GTF_URL"
    curl -fsSL --retry 3 -o "$GTF.tmp" "$GTF_URL"
    gzip -t "$GTF.tmp" && mv "$GTF.tmp" "$GTF"
fi

ENV="$(yaml_get conda_env_r_nad "$PARAMS")"
[[ -n "$ENV" ]] || die "conda_env_r_nad is not set in $PARAMS"
RSCRIPT="$(conda run -n "$ENV" which Rscript)"
log_info "Rscript: $RSCRIPT"

"$RSCRIPT" "$SCRIPT_DIR/deseq_s_vs_n.R" \
    --matrix "$MATRIX" --manifest "$MANIFEST" --gtf "$GTF" --fasta "$FASTA" \
    --outdir "$RES_DIR" --figdir "$FIG_DIR"

log_info "summary:"
column -t "$RES_DIR/summary.tsv"
