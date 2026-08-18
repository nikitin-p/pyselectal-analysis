#!/usr/bin/env bash
# Derive a gene-level BED6 from a GENCODE GRCh38 GTF.
#
# The counterpart of make_yeast_gene_bed.sh, with three differences that follow
# from GENCODE rather than from Ensembl yeast annotation.
#
# First, GENCODE emits explicit "gene" feature lines, so a gene's span is read
# directly rather than reduced over its transcripts.
#
# Second, the name column carries gene_name, which is the SYMBOL that
# org.Hs.eg.db indexes, so Stage 16 needs no identifier translation. This is the
# reason the BED is not built on gene_id: GENCODE gene_id values carry a version
# suffix and are Ensembl accessions, neither of which SYMBOL lookups accept.
# Genes with no gene_name fall back to their gene_id, which keeps them available
# for the nearest-gene assignment of Stage 9 even though Stage 16 will not match
# them.
#
# Third, pseudoautosomal genes are annotated twice, once on chrX and once on
# chrY with a _PAR_Y suffix on the gene_id. The chrY copies are dropped, because
# they are duplicates of the chrX entries under an identical symbol, and keeping
# both would let one promoter be assigned to two genes that are in fact one.
#
# Usage:
#   bash scripts/00_setup/make_human_gene_bed.sh [--params FILE] [--gtf FILE] \
#       [--out FILE] [--gene_types "protein_coding,lncRNA"] [--force]
#
# Defaults:
#   --params  config/params_human.yaml   (read for annotation_gtf)
#   --out     data/GRCh38.genes.bed
#
# --gene_types restricts the output to the listed GENCODE gene_type values. The
# default is every type, which matches the yeast BED and avoids silently
# discarding real promoters. Restricting to "protein_coding,lncRNA" makes the
# nearest-gene assignment of Stage 9 cleaner, at the cost of reassigning any TSS
# whose true target is a pseudogene or a small RNA to a more distant neighbour.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params_human.yaml"
GTF=""
BED=""
GENE_TYPES=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)     PARAMS="$2";     shift 2 ;;
        --gtf)        GTF="$2";        shift 2 ;;
        --out)        BED="$2";        shift 2 ;;
        --gene_types) GENE_TYPES="$2"; shift 2 ;;
        --force)      FORCE=1;         shift   ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -z "$GTF" ]] && GTF="$(yaml_get annotation_gtf "$PARAMS")"
[[ -z "$GTF" ]] && die "No GTF given and annotation_gtf is empty in $PARAMS"
[[ "$GTF" != /* ]] && GTF="$ROOT_DIR/$GTF"

[[ -z "$BED" ]] && BED="$ROOT_DIR/data/GRCh38.genes.bed"
[[ "$BED" != /* ]] && BED="$ROOT_DIR/$BED"

[[ -f "$GTF" ]] || die "GTF not found: $GTF"

if [[ -f "$BED" && $FORCE -eq 0 ]]; then
    log_info "[skip] $BED exists; use --force to regenerate"
    exit 0
fi

mkdir -p "$(dirname "$BED")"

log_info "Deriving gene spans from $GTF"
[[ -n "$GENE_TYPES" ]] && log_info "  restricting to gene types: $GENE_TYPES"

TMP_BED="$(mktemp "${TMPDIR:-/tmp}/human_gene_bed.XXXXXX")"
trap 'rm -f "$TMP_BED"' EXIT

awk -F'\t' -v types="$GENE_TYPES" '
    BEGIN {
        n_types = split(types, t, ",")
        for (i = 1; i <= n_types; i++) { gsub(/^[ \t]+|[ \t]+$/, "", t[i]); want[t[i]] = 1 }
    }
    $1 ~ /^#/    { next }
    $3 != "gene" { next }
    {
        if (match($9, /gene_id "[^"]+"/) == 0) next
        gid = substr($9, RSTART + 9, RLENGTH - 10)

        # Pseudoautosomal duplicates on chrY carry a _PAR_Y suffix; the chrX
        # copy is retained and this one is dropped.
        if (gid ~ /_PAR_Y$/) { n_par++; next }

        gtype = ""
        if (match($9, /gene_type "[^"]+"/) > 0)
            gtype = substr($9, RSTART + 11, RLENGTH - 12)
        if (n_types > 0 && !(gtype in want)) { n_filtered++; next }

        # gene_name is the SYMBOL that org.Hs.eg.db indexes; gene_id is only a
        # fallback so that unnamed genes still take part in Stage 9.
        name = gid
        if (match($9, /gene_name "[^"]+"/) > 0) {
            name = substr($9, RSTART + 11, RLENGTH - 12)
        } else {
            n_unnamed++
        }

        printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, $4 - 1, $5, name, $7
        n_out++
    }
    END {
        printf "genes written: %d\n", n_out                > "/dev/stderr"
        printf "PAR_Y dropped: %d\n", n_par + 0            > "/dev/stderr"
        printf "type-filtered: %d\n", n_filtered + 0       > "/dev/stderr"
        printf "no gene_name : %d\n", n_unnamed + 0        > "/dev/stderr"
    }
' "$GTF" | sort -k1,1 -k2,2n > "$TMP_BED"

[[ -s "$TMP_BED" ]] || die "No genes extracted from $GTF — check that it has 'gene' feature lines"

mv "$TMP_BED" "$BED"
trap - EXIT

N_GENES=$(wc -l < "$BED" | tr -d ' ')
N_UNIQ=$(cut -f4 "$BED" | sort -u | wc -l | tr -d ' ')
log_info "Wrote $N_GENES genes to $BED"

# Duplicate symbols survive on purpose: a handful of GENCODE symbols are used on
# more than one contig, and collapsing them would merge unrelated loci into one
# span. Stage 16 deduplicates by symbol anyway, so the count is reported rather
# than acted upon.
if [[ "$N_GENES" -ne "$N_UNIQ" ]]; then
    log_warn "$(( N_GENES - N_UNIQ )) rows share a name with another row ($N_UNIQ distinct names)"
fi
