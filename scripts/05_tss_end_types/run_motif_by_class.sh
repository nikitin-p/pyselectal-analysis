#!/usr/bin/env bash
# Stage 8 step 10: motifs per 5' end type, split by initiator dinucleotide class.
#
# 3Sggg is left out of the default type list (eight of the nine real end types) because
# its 16,369 reads do not survive a four-way split into readable logos. Pass a different
# --types list to override.
#
# The window is symmetric, running from -15 to +15, because the signal this stage is
# about sits immediately around the start site: the initiator dinucleotide, the bases
# that the reverse transcriptase reads through, and the first few templated positions.
# The wider upstream reach that the seqArchR stage needs in order to see the TATA box is
# not useful here, and it dilutes the logos it is added to.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

source scripts/utils/log.sh

FORCE=""
[[ "${1:-}" == "--force" ]] && FORCE="--force"

log_info "Running motif analysis by initiator dinucleotide class..."

Rscript scripts/05_tss_end_types/10_motif_by_dinuc_class.R \
    --level      ctss \
    --tss        results/tss/tss_matrix.tsv \
    --ctss       results/tss/typed_ctss.tsv.gz \
    --samples    config/samples.tsv \
    --genome     data/sacCer3.fa \
    --outdir     results/motifs \
    --figdir     results/figures \
    --types      1Sg,M,2Sgg,2Stg,1Sa,2Sga,1Sc,1St \
    --upstream   15 \
    --downstream 15 \
    --min_reads  200 \
    --min_ess    100 \
    $FORCE

log_info "Motif analysis by dinucleotide class complete"
