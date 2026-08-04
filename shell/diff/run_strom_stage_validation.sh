#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 7 ]; then
    cat >&2 <<'USAGE'
Usage:
  bash run_strom_stage_validation.sh \
      <execution_ST_DDROM_MPI.sh arguments...>

Example:
  bash run_strom_stage_validation.sh 20 1.0 10 4 4 8 0

Environment:
  ST_DDROM_EXEC_SCRIPT       execution wrapper path
                              default: shell/diff/execution_ST_DDROM_MPI.sh
  ST_DDROM_STAGE_ROOT        output root
                              default: ./stddrom_strom_stage_validation
  ST_DDROM_MAX_MODES         modes/rank, default: 10
  ST_DDROM_NUM_SNAPSHOTS     snapshot sets, default: 1
  ST_DDROM_STAGE_SLABS       multiple-slab count, default: 4
  ST_DDROM_STROM_COMPARE_TOL comparison tolerance, default: 1e-8
  ST_DDROM_REBUILD_STAGE     1: delete old stage data, default: 1
USAGE
    exit 2
fi

EXEC_SCRIPT="${ST_DDROM_EXEC_SCRIPT:-shell/diff/execution_ST_DDROM_MPI.sh}"
ROOT="${ST_DDROM_STAGE_ROOT:-./stddrom_strom_stage_validation}"
MAX_MODES="${ST_DDROM_MAX_MODES:-10}"
NUM_SNAPSHOTS="${ST_DDROM_NUM_SNAPSHOTS:-1}"
MULTI_SLABS="${ST_DDROM_STAGE_SLABS:-4}"
COMPARE_TOL="${ST_DDROM_STROM_COMPARE_TOL:-1.0e-8}"
REBUILD="${ST_DDROM_REBUILD_STAGE:-1}"

mkdir -p "$ROOT"
# The execution wrapper may change the working directory to the case directory.
# Pass absolute paths so CSV/snapshot/basis locations remain valid after that cd.
ROOT="$(cd "$ROOT" && pwd -P)"

run_stage() {
    local label="$1"
    local degree="$2"
    local slabs="$3"
    shift 3
    local stage_dir="$ROOT/$label"
    local snapshots="$stage_dir/snapshots"
    local basis="$stage_dir/basis_r${MAX_MODES}"
    local compare_csv="$stage_dir/stddrom_vs_sequential_strom.csv"
    local log="$stage_dir/run.log"

    if [ "$REBUILD" = "1" ]; then
        rm -rf "$stage_dir"
    fi
    mkdir -p "$stage_dir" "$snapshots" "$basis"

    echo
    echo "================================================================"
    echo "Stage: $label  dG($degree), slabs/window=$slabs"
    echo "================================================================"

    ST_DDROM_ACTION=collect \
    ST_DDROM_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_HLPOD_KERNELS=0 \
    ST_DDROM_COMPARE_STROM=0 \
    ST_DDROM_TIME_DEGREE="$degree" \
    ST_DDROM_SLABS_PER_WINDOW="$slabs" \
    ST_DDROM_FOM_FORMULATION=current_state \
    ST_DDROM_NUM_SNAPSHOTS="$NUM_SNAPSHOTS" \
    ST_DDROM_MAX_MODES="$MAX_MODES" \
    ST_DDROM_SNAPSHOT_DIR="$snapshots" \
    ST_DDROM_BASIS_DIR="$basis" \
    ST_DDROM_WRITE_OUTPUT=0 \
    ST_DDROM_WRITE_MODE_VTK=0 \
    bash "$EXEC_SCRIPT" "$@" 2>&1 | tee "$stage_dir/collect.log"

    ST_DDROM_ACTION=offline \
    ST_DDROM_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_HLPOD_KERNELS=0 \
    ST_DDROM_COMPARE_STROM=0 \
    ST_DDROM_TIME_DEGREE="$degree" \
    ST_DDROM_SLABS_PER_WINDOW="$slabs" \
    ST_DDROM_FOM_FORMULATION=current_state \
    ST_DDROM_NUM_SNAPSHOTS="$NUM_SNAPSHOTS" \
    ST_DDROM_MAX_MODES="$MAX_MODES" \
    ST_DDROM_SNAPSHOT_DIR="$snapshots" \
    ST_DDROM_BASIS_DIR="$basis" \
    ST_DDROM_WRITE_OUTPUT=0 \
    ST_DDROM_WRITE_MODE_VTK=0 \
    bash "$EXEC_SCRIPT" "$@" 2>&1 | tee "$stage_dir/offline.log"

    ST_DDROM_ACTION=online \
    ST_DDROM_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_STANDARD_DDROM=0 \
    ST_DDROM_COMPARE_HLPOD_KERNELS=0 \
    ST_DDROM_COMPARE_STROM=1 \
    ST_DDROM_TIME_DEGREE="$degree" \
    ST_DDROM_SLABS_PER_WINDOW="$slabs" \
    ST_DDROM_FOM_FORMULATION=current_state \
    ST_DDROM_MAX_MODES="$MAX_MODES" \
    ST_DDROM_SNAPSHOT_DIR="$snapshots" \
    ST_DDROM_BASIS_DIR="$basis" \
    ST_DDROM_STROM_COMPARE_CSV="$compare_csv" \
    ST_DDROM_STROM_COMPARE_TOL="$COMPARE_TOL" \
    ST_DDROM_REDUCED_EPSILON="${ST_DDROM_REDUCED_EPSILON:-1.0e-10}" \
    ST_DDROM_REDUCED_MAX_ITER="${ST_DDROM_REDUCED_MAX_ITER:-10000}" \
    ST_DDROM_VALIDATE="${ST_DDROM_VALIDATE:-1}" \
    ST_DDROM_WRITE_COMPARISON_VTK=0 \
    ST_DDROM_WRITE_MODE_VTK=0 \
    ST_DDROM_WRITE_OUTPUT=0 \
    bash "$EXEC_SCRIPT" "$@" 2>&1 | tee "$log"
}

run_stage dg0_s1 0 1 "$@"
run_stage dg0_s${MULTI_SLABS} 0 "$MULTI_SLABS" "$@"
run_stage dg1_s1 1 1 "$@"
run_stage dg1_s${MULTI_SLABS} 1 "$MULTI_SLABS" "$@"

python3 "$(dirname "$0")/summarize_strom_stages.py" "$ROOT"

echo
echo "Stage summary: $ROOT/stage_summary.csv"
