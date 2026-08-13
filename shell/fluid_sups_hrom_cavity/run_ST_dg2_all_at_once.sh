#!/bin/bash
set -Eeuo pipefail
repo_root="${1:-$(pwd)}"
cd "${repo_root}"

export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"
export FLUID_ST_TIME_DEGREE=2
export FLUID_ST_SLABS_PER_WINDOW="${FLUID_ST_SLABS_PER_WINDOW:-4}"
export FLUID_ST_WRITE_WINDOW_BINARY="${FLUID_ST_WRITE_WINDOW_BINARY:-1}"
export FLUID_ST_SOLVE_MODE=all_at_once
export FLUID_ST_DG_JACOBIAN_FD_EPS="${FLUID_ST_DG_JACOBIAN_FD_EPS:-1.0e-7}"

bash shell/fluid_sups_hrom_cavity/run_ST.sh
