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
