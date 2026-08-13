#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}"

export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"
export FLUID_ST_SOLVE_MODE=all_at_once
export FLUID_ST_TIME_DEGREE="${FLUID_ST_TIME_DEGREE:-1}"
export FLUID_ST_SLABS_PER_WINDOW="${FLUID_ST_SLABS_PER_WINDOW:-4}"

# Analytic kernel is used, but the first sampled element kernels are compared
# directly against finite differences of the linked verified residual provider.
export FLUID_ST_JOLD_MODE=verify
export FLUID_ST_JOLD_VERIFY_SAMPLES="${FLUID_ST_JOLD_VERIFY_SAMPLES:-64}"
export FLUID_ST_JOLD_VERIFY_TOL="${FLUID_ST_JOLD_VERIFY_TOL:-5.0e-6}"
export FLUID_ST_DG_JACOBIAN_FD_EPS="${FLUID_ST_DG_JACOBIAN_FD_EPS:-1.0e-7}"

bash shell/fluid_sups_hrom_cavity/run_ST.sh
