#!/usr/bin/env bash
# Download NAD captureSeq FASTQ files (GSE146368, Zhang et al. 2020) from ENA
# over HTTPS and verify each against the ENA md5. Idempotent: a file whose md5
# already matches is skipped. Runs sequentially; ENA throttles parallel pulls.
#
# Usage (from the repo root):
#   bash nad/scripts/00_download/download.sh [MANIFEST] [OUTDIR]
# Defaults: nad/scripts/00_download/priority1_runs.tsv, nad/data/nad_captureseq
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
manifest="${1:-$here/priority1_runs.tsv}"
outdir="${2:-$repo/nad/data/nad_captureseq}"
logdir="$repo/nad/logs"
mkdir -p "$outdir" "$logdir"
log="$logdir/00_download_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1

echo "manifest: $manifest"
echo "outdir:   $outdir"

failed=0
# columns: run strain library replicate read_count fastq_bytes md5 url
while IFS=$'\t' read -r run strain lib rep reads bytes md5 url; do
  [[ "$run" == "run" ]] && continue
  dest="$outdir/$run.fastq.gz"
  if [[ -s "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
    echo "[skip] $run ($strain $lib $rep) already present, md5 ok"
    continue
  fi
  echo "[get]  $run ($strain $lib $rep), $bytes bytes"
  ok=0
  for attempt in 1 2 3; do
    if wget -q -c --tries=5 --timeout=60 -O "$dest" "$url" \
       && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
      ok=1; break
    fi
    echo "       attempt $attempt failed for $run, retrying from scratch"
    rm -f "$dest"
  done
  if (( ok )); then echo "[ok]   $run md5 verified"; else echo "[FAIL] $run"; failed=$((failed+1)); fi
done < "$manifest"

# Human-readable sample table beside the data, for later stages.
cut -f1-5 "$manifest" > "$outdir/samples.tsv"

echo "done, failures: $failed"
exit $(( failed > 0 ))
