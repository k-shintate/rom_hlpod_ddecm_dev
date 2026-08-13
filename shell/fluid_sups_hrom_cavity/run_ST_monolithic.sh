#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}"

export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"
export FLUID_ST_SLABS_PER_WINDOW="${FLUID_ST_SLABS_PER_WINDOW:-4}"
export FLUID_ST_WRITE_WINDOW_BINARY="${FLUID_ST_WRITE_WINDOW_BINARY:-1}"
export FLUID_ST_SOLVE_MODE=monolithic
export FLUID_ST_PREVIOUS_JVP_EPS="${FLUID_ST_PREVIOUS_JVP_EPS:-1.0e-6}"

bash shell/fluid_sups_hrom_cavity/run_ST.sh
