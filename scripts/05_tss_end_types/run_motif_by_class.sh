#!/usr/bin/env bash
# Stage 8 step 10: motifs per 5' end type, split by initiator dinucleotide class.
#
# 3Sggg is left out of the default type list because its 16,369 reads do not survive a
# four-way split into readable logos. Pass a different --types list to override.
#
# The window is asymmetric, running from -35 to +20, for the same reason the seqArchR
# stage uses an asymmetric window: the yeast core promoter extends considerably further
# upstream of the start site than downstream of it, so a symmetric window spends half its
# width on transcript body sequence that carries no promoter signal.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/utils/log.sh

FORCE=""
[[ "${1:-}" == "--force" ]] && FORCE="--force"

log_info "Running motif analysis by initiator dinucleotide class..."

Rscript scripts/05_tss_end_types/10_motif_by_dinuc_class.R \
    --ctss       results/tss/typed_ctss.tsv.gz \
    --samples    config/samples.tsv \
    --genome     data/sacCer3.fa \
    --outdir     results/motifs \
    --figdir     results/figures \
    --types      1Sg,M,2Sgg,1Sa,1Sc,1St \
    --upstream   35 \
    --downstream 20 \
    --min_reads  200 \
    --min_ess    100 \
    $FORCE

log_info "Motif analysis by dinucleotide class complete"
