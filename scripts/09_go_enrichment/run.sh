#!/usr/bin/env bash
# Stage 16: GO Biological Process enrichment by 5' end type and by initiator.
#
# Asks whether promoters that favour a given 5' end type, or that carry a given
# initiator dinucleotide, serve distinct biological processes. Both families are
# tested against the same universe of genes with an assigned, expressed TSS, so
# a hit means enrichment relative to the detected transcriptome rather than
# relative to the genome.
#
# Usage:
#   bash scripts/09_go_enrichment/run.sh [--force]

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
        --params)  PARAMS="$2";  shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --force)   FORCE_FLAG="--force"; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"

OUTDIR="$ROOT_DIR/$(yaml_get results_go "$PARAMS")"
FIGDIR="$ROOT_DIR/$(yaml_get results_figures "$PARAMS")"
# The generic annotation_bed / go_orgdb keys are read first; the older
# *_yeast keys remain as a fallback for params files that predate them.
GENES_BED="$(yaml_get annotation_bed "$PARAMS")"
[[ -z "$GENES_BED" ]] && GENES_BED="$(yaml_get annotation_bed_yeast "$PARAMS")"
[[ -z "$GENES_BED" ]] && die "No gene BED configured: set annotation_bed in $PARAMS"
[[ "$GENES_BED" != /* ]] && GENES_BED="$ROOT_DIR/$GENES_BED"

QUANTILE="$(yaml_get go_quantile          "$PARAMS")"
MIN_TYPE="$(yaml_get go_min_type_reads    "$PARAMS")"
MIN_TOTAL="$(yaml_get go_min_total        "$PARAMS")"
UPSTREAM="$(yaml_get go_upstream_window   "$PARAMS")"
DOWNSTREAM="$(yaml_get go_downstream_window "$PARAMS")"
PVAL="$(yaml_get go_pvalue_cutoff         "$PARAMS")"
QVAL="$(yaml_get go_qvalue_cutoff         "$PARAMS")"
SIMPLIFY="$(yaml_get go_simplify_cutoff   "$PARAMS")"
MIN_SET="$(yaml_get go_min_set_size       "$PARAMS")"
ORGDB="$(yaml_get go_orgdb "$PARAMS")"
[[ -z "$ORGDB" ]] && ORGDB="$(yaml_get go_orgdb_yeast "$PARAMS")"
KEYTYPE="$(yaml_get go_keytype            "$PARAMS")"

# clusterProfiler is installed in the system R library, not in conda_env_r.
RSCRIPT="$(which Rscript)"

[[ -f "$GENES_BED" ]] || {
    log_info "Gene BED missing; deriving it from the GTF"
    bash "$ROOT_DIR/scripts/00_setup/make_yeast_gene_bed.sh"
}

log_info "Stage 16 step 0: checking dependencies"
$RSCRIPT "$SCRIPT_DIR/00_install.R"

log_info "Stage 16 step 1: assigning TSS clusters to genes"
$RSCRIPT "$SCRIPT_DIR/01_assign_genes.R" \
    --tss_matrix    "$RESULTS_ROOT/tss/tss_matrix.tsv"    \
    --tss_clustered "$RESULTS_ROOT/tss/tss_clustered.tsv" \
    --genes_bed     "$GENES_BED"    \
    --outdir        "$OUTDIR"       \
    --upstream      "$UPSTREAM"     \
    --downstream    "$DOWNSTREAM"   \
    --min_total     "$MIN_TOTAL"    \
    $FORCE_FLAG

log_info "Stage 16 step 2: enrichment by 5' end type"
$RSCRIPT "$SCRIPT_DIR/02_enrich_endtype.R" \
    --outdir         "$OUTDIR"   \
    --quantile       "$QUANTILE" \
    --min_type_reads "$MIN_TYPE" \
    --min_set        "$MIN_SET"  \
    --pvalue         "$PVAL"     \
    --qvalue         "$QVAL"     \
    --simplify       "$SIMPLIFY" \
    --orgdb          "$ORGDB"    \
    --keytype        "$KEYTYPE"  \
    $FORCE_FLAG

log_info "Stage 16 step 3: enrichment by initiator dinucleotide"
$RSCRIPT "$SCRIPT_DIR/03_enrich_initiator.R" \
    --outdir   "$OUTDIR"   \
    --min_set  "$MIN_SET"  \
    --pvalue   "$PVAL"     \
    --qvalue   "$QVAL"     \
    --simplify "$SIMPLIFY" \
    --orgdb    "$ORGDB"    \
    --keytype  "$KEYTYPE"  \
    $FORCE_FLAG

log_info "Stage 16 step 4: figures"
$RSCRIPT "$SCRIPT_DIR/04_plot.R" \
    --outdir "$OUTDIR" --figdir "$FIGDIR" $FORCE_FLAG

log_info "Stage 16 complete."
log_info "  tables  → $OUTDIR/"
log_info "  figures → $FIGDIR/17_go_*.pdf"
