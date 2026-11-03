#!/usr/bin/env bash
# NAD branch step 2: trim NAD captureSeq reads down to the RNA insert.
#
# For every run in the manifest:
#   cutadapt  - NextSeq poly-G quality trim, then remove the 3' adapter
#               AGATCGGAAGAGC... (reads without it are kept: long inserts
#               are not read through to the adapter);
#   umi_grun.awk - cut C + 6-nt UMI before the adapter, move the UMI into the
#               read name, strip the leading G run from the TdT C-tail, record
#               its length in the read name, drop inserts < 18 nt.
# Outputs:
#   nad/data/trimmed/<run>.fastq.gz
#   nad/results/02_trim/<run>.cutadapt.json   cutadapt report
#   nad/results/02_trim/<run>.grun.tsv        leading-G run length histogram
#   nad/results/02_trim/<run>.stats.tsv       reads in / out / with UMI / too short
#   nad/results/02_trim/summary.tsv           one row per run
# Idempotent: a run whose trimmed FASTQ exists is skipped unless --force.
#
# Usage (from anywhere; submit through Slurm on vavilov):
#   bash nad/scripts/02_trim/run.sh [--params config/params.yaml]
#        [--manifest nad/scripts/00_download/priority1_runs.tsv] [--threads 8] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
FORCE=0
MINLEN=18
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)   PARAMS="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

RAW_DIR="$NAD_DIR/data/nad_captureseq"
OUT_DIR="$NAD_DIR/data/trimmed"
RES_DIR="$NAD_DIR/results/02_trim"
mkdir -p "$OUT_DIR" "$RES_DIR" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/02_trim_$(date +%Y%m%d_%H%M%S).log") 2>&1

CONDA_TOOLS="$(yaml_get conda_env_tools "$PARAMS")"
CUTADAPT="$(conda run -n "$CONDA_TOOLS" which cutadapt)"
[[ -x "$CUTADAPT" ]] || die "cutadapt not found in conda env $CONDA_TOOLS"
log_info "cutadapt: $CUTADAPT ($("$CUTADAPT" --version))"
log_info "threads: $THREADS, min insert: $MINLEN nt"

tail -n +2 "$MANIFEST" | cut -f1 | while read -r run; do
    in="$RAW_DIR/$run.fastq.gz"
    out="$OUT_DIR/$run.fastq.gz"
    [[ -s "$in" ]] || die "missing raw FASTQ: $in (run step 00_download first)"
    if [[ -s "$out" && $FORCE -eq 0 ]]; then
        log_info "[skip] $run: $out exists"
        continue
    fi
    log_info "[trim] $run"
    tmp="$out.tmp"
    # cutadapt applies --nextseq-trim before adapter removal; -O 5 keeps short
    # partial adapters at the read end from being missed, -e 0.1 is the default.
    "$CUTADAPT" -j "$THREADS" --nextseq-trim=20 \
        -a "ad=$ADAPTER" -O 5 -e 0.1 \
        --rename '{id} ad={adapter_name}' \
        --json "$RES_DIR/$run.cutadapt.json" \
        -o - "$in" 2> "$RES_DIR/$run.cutadapt.log" \
      | awk -v MINLEN="$MINLEN" \
            -v GRUN_TSV="$RES_DIR/$run.grun.tsv" \
            -v STATS_TSV="$RES_DIR/$run.stats.tsv" \
            -f "$SCRIPT_DIR/umi_grun.awk" \
      | gzip -c > "$tmp"
    mv "$tmp" "$out"
    log_info "[ok]   $run: $(tail -n 1 "$RES_DIR/$run.stats.tsv")"
done

# Summary over all runs in the manifest.
{
    printf "run\tstrain\tlibrary\treplicate\treads_in\treads_out\treads_out_with_umi\treads_too_short\n"
    tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run strain lib rep _rest; do
        s="$RES_DIR/$run.stats.tsv"
        [[ -s "$s" ]] && printf "%s\t%s\t%s\t%s\t%s\n" "$run" "$strain" "$lib" "$rep" "$(tail -n 1 "$s")"
    done
} > "$RES_DIR/summary.tsv"
log_info "summary: $RES_DIR/summary.tsv"
column -t "$RES_DIR/summary.tsv"
