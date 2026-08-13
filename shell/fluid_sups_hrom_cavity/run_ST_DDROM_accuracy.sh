#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "$repo_root"

export FLUID_ST_DDROM_MODE="${FLUID_ST_DDROM_MODE:-online}"
export FLUID_ST_DDROM_VALIDATE="${FLUID_ST_DDROM_VALIDATE:-1}"
export FLUID_ST_DDROM_ACCURACY_CSV="${FLUID_ST_DDROM_ACCURACY_CSV:-$PWD/result_fluid_sups_cavity/10-${FLUID_ST_MPI_RANKS:-4}-1/fluid_st_ddrom_accuracy.csv}"

exec bash shell/fluid_sups_hrom_cavity/run_ST_DDROM.sh "$PWD"
