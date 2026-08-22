#!/usr/bin/env bash
# Fair timing comparison for diffusion ST-DDROM.
# 1) Build/train once (offline).
# 2) Run online all-at-once using the stored basis/operators.
# 3) Run online window-marching using exactly the same basis/operators/RHS setup.
# VTK output is disabled by default to reduce I/O contamination.

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
    echo "Usage: $0 <e> <ep> <nm> <nd> <np> <pa> <st>" >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
launcher="${ST_DDROM_LAUNCHER:-$script_dir/execution_ST_spaceblock_MPI.sh}"
if [[ ! -x "$launcher" ]]; then
    echo "ERROR: launcher not found/executable: $launcher" >&2
    exit 1
fi

write_output="${DIFF_ST_DDROM_WRITE_OUTPUT:-0}"
initial_rebuild="${ST_FOM_REBUILD:-1}"

printf '\n=== 1/3 offline POD/operator generation ===\n'
env \
    ST_FOM_REBUILD="$initial_rebuild" \
    ST_CLEAN_RESULT=1 \
    ST_DDROM_RUN_LABEL="stddrom_offline_compare_np${5}" \
    DIFF_ST_DDROM_MODE=offline \
    DIFF_ST_DDROM_SOLVER=all-at-once \
    DIFF_ST_DDROM_WRITE_OUTPUT="$write_output" \
    "$launcher" "$@"

printf '\n=== 2/3 online all-at-once ===\n'
env \
    ST_FOM_REBUILD=0 \
    ST_CLEAN_RESULT=0 \
    ST_DDROM_RUN_LABEL="stddrom_online_all_at_once_np${5}" \
    DIFF_ST_DDROM_MODE=online \
    DIFF_ST_DDROM_SOLVER=all-at-once \
    DIFF_ST_DDROM_WRITE_OUTPUT="$write_output" \
    "$launcher" "$@"

printf '\n=== 3/3 online window-marching ===\n'
env \
    ST_FOM_REBUILD=0 \
    ST_CLEAN_RESULT=0 \
    ST_DDROM_RUN_LABEL="stddrom_online_window_marching_np${5}" \
    DIFF_ST_DDROM_MODE=online \
    DIFF_ST_DDROM_SOLVER=window-marching \
    DIFF_ST_DDROM_WRITE_OUTPUT="$write_output" \
    "$launcher" "$@"

cat <<'MSG'

Comparison complete.
The two online runs use the same offline basis and reduced operators.
Compare the two online rows in stddrom_timing.csv, especially:
  all_at_once_system_assembly + all_at_once_solve
versus:
  window_marching_solve
and compare Total ST-DDROM time for end-to-end online cost.
MSG
