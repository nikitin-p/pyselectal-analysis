#!/usr/bin/env bash
# NAD branch step 4: strand-aware 5'-end counts per position from step-3 BAMs.
#
# Outputs:
#   nad/results/04_ctss/<run>.ctss.tsv.gz   chr pos strand base_p1 reads umis mm5 g0
#   nad/results/04_ctss/matrix_umis.tsv.gz  position x run, UMI-collapsed
#   nad/results/04_ctss/matrix_reads.tsv.gz position x run, raw reads
#   nad/results/04_ctss/first_base.tsv      per run: share of reads by base_p1
# The matrices are small and are what step 5 (S vs N enrichment, R) reads.
# Idempotent: per-run tables are skipped if present unless --force.
#
# Usage: bash nad/scripts/04_ctss/run.sh [--params config/params.yaml]
#        [--manifest nad/scripts/00_download/priority1_runs.tsv] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)   PARAMS="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

FASTA="$(yaml_get genome_fasta "$PARAMS")"; [[ "$FASTA" != /* ]] && FASTA="$ROOT_DIR/$FASTA"
BAM_DIR="$NAD_DIR/data/bam"
RES_DIR="$NAD_DIR/results/04_ctss"
mkdir -p "$RES_DIR" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/04_ctss_$(date +%Y%m%d_%H%M%S).log") 2>&1

PYENV="$(yaml_get conda_env_pyselectal "$PARAMS")"
PY=(conda run --no-capture-output -n "$PYENV" python -I)

# pysam.FastaFile needs a .fai beside the FASTA; build it if missing.
[[ -s "$FASTA.fai" ]] || "${PY[@]}" -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$FASTA"

# Runs are independent; with one CPU per run they finish in parallel.
pids=()
for run in $(tail -n +2 "$MANIFEST" | cut -f1); do
    out="$RES_DIR/$run.ctss.tsv.gz"
    [[ -s "$BAM_DIR/$run.bam" ]] || die "missing BAM: $BAM_DIR/$run.bam (run step 03_map first)"
    if [[ -s "$out" && $FORCE -eq 0 ]]; then log_info "[skip] $run"; continue; fi
    log_info "[count] $run"
    ( "${PY[@]}" "$SCRIPT_DIR/count_5p.py" --bam "$BAM_DIR/$run.bam" --fasta "$FASTA" --out "$out.tmp" \
        && mv "$out.tmp" "$out" ) &
    pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail == 0 )) || die "at least one run failed"

"${PY[@]}" "$SCRIPT_DIR/make_matrix.py" --manifest "$MANIFEST" --ctss-dir "$RES_DIR" --outdir "$RES_DIR"

# First-base composition per run: the NAD captures should be A-rich in S over N.
{
    printf "run\tbase_p1\treads\tfrac\n"
    for run in $(tail -n +2 "$MANIFEST" | cut -f1); do
        zcat "$RES_DIR/$run.ctss.tsv.gz" | awk -v r="$run" -F'\t' 'NR>1{n[$4]+=$5; t+=$5}
            END{for(b in n) printf "%s\t%s\t%d\t%.4f\n", r, b, n[b], n[b]/t}' | sort -k2,2
    done
} > "$RES_DIR/first_base.tsv"
log_info "first-base composition (A share per run):"
awk -F'\t' '$2=="A"' "$RES_DIR/first_base.tsv" | column -t
