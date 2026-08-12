#!/usr/bin/env bash
# Derive a gene-level BED6 from the sacCer3 Ensembl GTF.
#
# The GTF is transcript-level, so a gene's span is taken as the minimum start
# and the maximum end over all of its transcripts. The name column carries the
# systematic ORF name, which is the key that org.Sc.sgd.db indexes under
# keyType = "ORF", so downstream GO analysis needs no identifier translation.
#
# Usage:
#   bash scripts/00_setup/make_yeast_gene_bed.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"

GTF="$ROOT_DIR/data/sacCer3.ensGene.gtf"
BED="$ROOT_DIR/data/sacCer3.genes.bed"

[[ -f "$GTF" ]] || die "GTF not found: $GTF"

if [[ -f "$BED" && "${1:-}" != "--force" ]]; then
    log_info "[skip] $BED exists; use --force to regenerate"
    exit 0
fi

log_info "Deriving gene spans from $GTF"

awk -F'\t' '
    $1 ~ /^#/       { next }
    $3 != "transcript" { next }
    {
        if (match($9, /gene_id "[^"]+"/) == 0) next
        gid = substr($9, RSTART + 9, RLENGTH - 10)
        key = gid SUBSEP $1 SUBSEP $7
        s = $4 - 1                       # GTF is 1-based inclusive; BED is 0-based
        if (!(key in lo) || s  < lo[key]) lo[key] = s
        if (!(key in hi) || $5 > hi[key]) hi[key] = $5
    }
    END {
        for (k in lo) {
            split(k, a, SUBSEP)
            printf "%s\t%d\t%d\t%s\t0\t%s\n", a[2], lo[k], hi[k], a[1], a[3]
        }
    }
' "$GTF" | sort -k1,1 -k2,2n > "$BED"

log_info "Wrote $(wc -l < "$BED" | tr -d ' ') genes to $BED"
