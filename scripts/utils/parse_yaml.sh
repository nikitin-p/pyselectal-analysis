#!/usr/bin/env bash
# Minimal YAML key reader (flat keys only, no nesting).
# Usage: yaml_get KEY FILE
# Returns the value for "KEY: value" lines; strips quotes and inline comments.

yaml_get() {
    local key="$1" file="$2"
    grep -E "^${key}:" "$file" \
        | head -1 \
        | sed -E "s/^${key}:[[:space:]]*//" \
        | sed 's/#.*//' \
        | sed -E 's/^[[:space:]]*"?([^"]*)"?[[:space:]]*/\1/' \
        | tr -d '\r'
}

# Resolve the root of the results tree for a given params file.
# Usage: RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"
#
# Every stage writes beneath this root, so that one dataset's run cannot
# overwrite another's: the per-sample outputs are keyed by sample_id, but the
# aggregate tables and the figures are not, and they would otherwise collide.
# Params files that predate the key fall back to "results", which is where the
# yeast run already lives.
resolve_results_root() {
    local params="$1" root_dir="$2" value
    value="$(yaml_get results_root "$params")"
    [[ -z "$value" ]] && value="results"
    [[ "$value" != /* ]] && value="$root_dir/$value"
    printf '%s' "$value"
}
