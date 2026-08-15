#!/usr/bin/env bash
set -Eeuo pipefail
mode="${FLUID_ST_DDROM_MODE:-online}"
case "$mode" in offline|online) ;; *) echo "ERROR: FLUID_ST_DDROM_MODE must be offline or online" >&2; exit 2;; esac
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/st_pipeline.sh" "$mode"
