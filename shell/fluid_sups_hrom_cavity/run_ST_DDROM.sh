#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "$repo_root"

np="${FLUID_ST_MPI_RANKS:-4}"
nm="${FLUID_ST_RESULT_MODES_TAG:-10}"
nd="${FLUID_ST_RESULT_DD_TAG:-1}"
result_dir="${repo_root}/result_fluid_sups_cavity/${nm}-${np}-${nd}"

exe_src="${repo_root}/solvers/fluid_sups/fluid_sups_st_ddrom"
exe_dst="${result_dir}/fluid_sups_st_ddrom"

[[ -x "$exe_src" ]] || {
    echo "ERROR: executable not found: $exe_src" >&2
    echo "Build with: make -f Makefile_ST_cavity dev_build" >&2
    exit 1
}

[[ -d "$result_dir" ]] || {
    echo "ERROR: result directory not found: $result_dir" >&2
    echo "Run the GEDATSU FOM preparation/collection first." >&2
    exit 1
}

cp -f "$exe_src" "$exe_dst"
cd "$result_dir"

export FLUID_ST_SOLVE_MODE="${FLUID_ST_SOLVE_MODE:-all_at_once}"
export FLUID_ST_JOLD_MODE="${FLUID_ST_JOLD_MODE:-fd}"

# Keep the DDROM initial state identical to the cavity FOM collection run.
# fluid_sups_st_model defaults to the Karman-vortex initial state when this
# variable is absent, which corrupts only the first-window causal trace and
# then contaminates the entire ROM march.
export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"

export FLUID_ST_DDROM_MODE="${FLUID_ST_DDROM_MODE:-online}"
export FLUID_ST_DDROM_MAX_MODES="${FLUID_ST_DDROM_MAX_MODES:-10}"
export FLUID_ST_DDROM_POD_ENERGY="${FLUID_ST_DDROM_POD_ENERGY:-1.0}"

export FLUID_ST_DDROM_SNAPSHOT_DIR="${FLUID_ST_DDROM_SNAPSHOT_DIR:-${result_dir}/fluid_st_windows}"
export FLUID_ST_DDROM_BASIS_DIR="${FLUID_ST_DDROM_BASIS_DIR:-${result_dir}/fluid_st_ddrom_basis}"
export FLUID_ST_DDROM_WINDOW_DIR="${FLUID_ST_DDROM_WINDOW_DIR:-${result_dir}/fluid_st_ddrom_windows}"
export FLUID_ST_DDROM_COEFFICIENT_CSV="${FLUID_ST_DDROM_COEFFICIENT_CSV:-${result_dir}/fluid_st_ddrom_coefficients.csv}"

export FLUID_ST_DDROM_WRITE_OUTPUT="${FLUID_ST_DDROM_WRITE_OUTPUT:-1}"
export FLUID_ST_DDROM_WRITE_WINDOW_BINARY="${FLUID_ST_DDROM_WRITE_WINDOW_BINARY:-1}"

log="fluid_st_ddrom_${FLUID_ST_DDROM_MODE}_np${np}.log"

echo "ST-DDROM runner initial condition: ${FLUID_ST_INITIAL_CONDITION}"

mpirun -np "$np" ./fluid_sups_st_ddrom ./ \
    2>&1 | tee "$log"
