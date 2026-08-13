#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "$repo_root"

export FLUID_ST_DDROM_MAX_MODES="${FLUID_ST_DDROM_MAX_MODES:-10}"
export FLUID_ST_DDROM_MAX_VELOCITY_MODES="${FLUID_ST_DDROM_MAX_VELOCITY_MODES:-$FLUID_ST_DDROM_MAX_MODES}"
export FLUID_ST_DDROM_MAX_PRESSURE_MODES="${FLUID_ST_DDROM_MAX_PRESSURE_MODES:-$FLUID_ST_DDROM_MAX_MODES}"

exec bash shell/fluid_sups_hrom_cavity/run_ST_DDROM.sh "$PWD"
