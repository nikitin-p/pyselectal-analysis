#!/usr/bin/env bash
# Logging helpers. Source this file: source scripts/utils/log.sh
# Usage: log_info "message"  |  log_warn "message"  |  log_error "message"

LOG_TS() { date '+%Y-%m-%d %H:%M:%S'; }

log_info()  { echo "[$(LOG_TS)] INFO  $*" >&2; }
log_warn()  { echo "[$(LOG_TS)] WARN  $*" >&2; }
log_error() { echo "[$(LOG_TS)] ERROR $*" >&2; }

# die MESSAGE — print error and exit 1
die() { log_error "$*"; exit 1; }

# require_cmd CMD — exit if command not found
require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

# log_to_file LOGFILE COMMAND [ARGS...] — run command, tee stderr to logfile
log_to_file() {
    local logfile="$1"; shift
    mkdir -p "$(dirname "$logfile")"
    "$@" 2> >(tee -a "$logfile" >&2)
}
