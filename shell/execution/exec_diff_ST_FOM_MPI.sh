#!/usr/bin/env bash

# =============================================================================
# Normal spatially-partitioned Space-Time FOM launcher for diffusion.
#
# This script does NOT use ROM/POD/ST-DDROM.  It:
#   1. builds diff_ST_MPI from Makefile_ST,
#   2. creates the ordinary spatial diffusion case,
#   3. partitions only the spatial mesh with GEDATSU,
#   4. runs the normal all-at-once ST-FOM by MPI,
#   5. checks st_accuracy.csv against the manufactured solution.
#
# Usage:
#   exec_diff_ST_FOM_MPI.sh <e> <ep> <np> <degree> <slabs/window> <dt> <finish_time>
#
# Example:
#   bash shell/execution/exec_diff_ST_FOM_MPI.sh 20 1.0 4 1 4 0.01 1.0
#
# Environment overrides:
#   ST_PROJECT_ROOT=/path/to/rom_hlpod_ddecm_dev
#   TEST_THERMAL_ROOT=/path/to/test_thermal
#   ST_FOM_REBUILD=1|0              default: 1
#   ST_CLEAN_RESULT=1|0             default: 1
#   ST_PREPARE_ONLY=1|0             default: 0
#   ST_FOM_OUTPUT_INTERVAL=N         default: 1
#   ST_FOM_NUM_IP_EACH_AXIS=N        default: 3
#   ST_FOM_MAT_EPSILON=value         default: 1.0e-10
#   ST_FOM_MAT_MAX_ITER=N            default: 10000
#   ST_FOM_RESULT_TAG=name           optional result-directory tag
#   ST_FOM_ACCURACY_CSV=name.csv     default: st_accuracy.csv
#   ST_FOM_RUN_LABEL=name            default: fom_st_np${np}_dg${degree}
#   OMP_NUM_THREADS=N                default: 1
# =============================================================================

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  exec_diff_ST_FOM_MPI.sh <e> <ep> <np> <degree> <slabs/window> <dt> <finish_time>

Arguments:
  e              number of hexahedral elements in each spatial direction
  ep             domain length in each spatial direction
  np             MPI process count / spatial partition count
  degree         dG time degree: 0, 1, or 2
  slabs/window   number of time slabs solved together in one ST window
  dt             time-step size
  finish_time    final physical time

Example:
  bash shell/execution/exec_diff_ST_FOM_MPI.sh 20 1.0 4 1 4 0.01 1.0

This solves the normal full-order space-time problem and writes
st_accuracy.csv containing relative L2, absolute L2, and Linf errors against
the manufactured solution.
USAGE
}

if [[ "$#" -ne 7 ]]; then
    usage
    exit 2
fi

e="$1"
ep="$2"
np="$3"
degree="$4"
slabs="$5"
dt="$6"
finish_time="$7"

is_positive_int() {
    [[ "$1" =~ ^[0-9]+$ ]] && (( 10#$1 > 0 ))
}

is_nonnegative_int() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

is_positive_number() {
    [[ "$1" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]] &&
    awk -v x="$1" 'BEGIN { exit !(x > 0.0) }'
}

if ! is_positive_int "$e"; then
    echo "ERROR: e must be a positive integer: $e" >&2
    exit 2
fi
if ! is_positive_int "$np"; then
    echo "ERROR: np must be a positive integer: $np" >&2
    exit 2
fi
if ! is_nonnegative_int "$degree" || (( degree < 0 || degree > 2 )); then
    echo "ERROR: degree must be 0, 1, or 2: $degree" >&2
    exit 2
fi
if ! is_positive_int "$slabs"; then
    echo "ERROR: slabs/window must be a positive integer: $slabs" >&2
    exit 2
fi
if ! is_positive_number "$ep"; then
    echo "ERROR: ep must be a positive number: $ep" >&2
    exit 2
fi
if ! is_positive_number "$dt"; then
    echo "ERROR: dt must be a positive number: $dt" >&2
    exit 2
fi
if ! is_positive_number "$finish_time"; then
    echo "ERROR: finish_time must be a positive number: $finish_time" >&2
    exit 2
fi

# Verify that finish_time/dt is integral and that complete ST windows fit exactly.
if ! total_steps="$(
    awk -v tf="$finish_time" -v h="$dt" -v sw="$slabs" '
        BEGIN {
            steps = tf / h;
            n = int(steps + 0.5);
            scale = (steps < 0 ? -steps : steps);
            if (scale < 1.0) scale = 1.0;
            diff = steps - n;
            if (diff < 0) diff = -diff;
            if (diff > 1.0e-9 * scale) exit 10;
            if (n <= 0) exit 11;
            if ((n % sw) != 0) exit 12;
            printf "%d", n;
        }'
)"; then
    echo "ERROR: finish_time/dt must be an integer and total steps must be divisible by slabs/window." >&2
    echo "       finish_time=$finish_time dt=$dt slabs/window=$slabs" >&2
    exit 2
fi
num_windows=$(( total_steps / slabs ))

script_path="${BASH_SOURCE[0]:-$0}"
script_dir="$(cd "$(dirname "$script_path")" && pwd)"

find_project_root() {
    local candidate
    for candidate in \
        "${ST_PROJECT_ROOT:-}" \
        "$(pwd)" \
        "${script_dir}/../.." \
        "${script_dir}/../../.."; do
        [[ -n "$candidate" ]] || continue
        candidate="$(cd "$candidate" 2>/dev/null && pwd || true)"
        [[ -n "$candidate" ]] || continue
        if [[ -f "$candidate/shell/install.sh" && \
              -f "$candidate/solvers/diff/Makefile_ST" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

project_root="$(find_project_root || true)"
if [[ -z "$project_root" ]]; then
    echo "ERROR: project root was not found." >&2
    echo "Set ST_PROJECT_ROOT to rom_hlpod_ddecm_dev." >&2
    exit 1
fi

# Keep the project's normal environment initialization.
# shellcheck source=/dev/null
source "$project_root/shell/install.sh"

solver_dir="$project_root/solvers/diff"

rebuild_solver="${ST_FOM_REBUILD:-1}"
clean_result="${ST_CLEAN_RESULT:-1}"
prepare_only="${ST_PREPARE_ONLY:-0}"
omp_threads="${OMP_NUM_THREADS:-1}"
output_interval="${ST_FOM_OUTPUT_INTERVAL:-1}"
num_ip_each_axis="${ST_FOM_NUM_IP_EACH_AXIS:-3}"
mat_epsilon="${ST_FOM_MAT_EPSILON:-1.0e-10}"
mat_max_iter="${ST_FOM_MAT_MAX_ITER:-10000}"
accuracy_csv="${ST_FOM_ACCURACY_CSV:-st_accuracy.csv}"
solve_mode="${ST_FOM_SOLVE_MODE:-causal}"
graph_source="${ST_FOM_GRAPH_SOURCE:-partition}"
run_label="${ST_FOM_RUN_LABEL:-fom_st_np${np}_dg${degree}}"
result_tag="${ST_FOM_RESULT_TAG:-stfom-e${e}-np${np}-dg${degree}-sw${slabs}-dt${dt}-tf${finish_time}}"
result_dir="$project_root/result_diff/$result_tag"

if ! is_positive_int "$output_interval"; then
    echo "ERROR: ST_FOM_OUTPUT_INTERVAL must be a positive integer." >&2
    exit 2
fi
if ! is_positive_int "$num_ip_each_axis"; then
    echo "ERROR: ST_FOM_NUM_IP_EACH_AXIS must be a positive integer." >&2
    exit 2
fi
if ! is_positive_int "$mat_max_iter"; then
    echo "ERROR: ST_FOM_MAT_MAX_ITER must be a positive integer." >&2
    exit 2
fi
if ! is_positive_number "$mat_epsilon"; then
    echo "ERROR: ST_FOM_MAT_EPSILON must be positive." >&2
    exit 2
fi
if ! is_positive_int "$omp_threads"; then
    echo "ERROR: OMP_NUM_THREADS must be a positive integer." >&2
    exit 2
fi

case "$solve_mode" in
    causal|window|all-at-once) ;;
    *)
        echo "ERROR: ST_FOM_SOLVE_MODE must be causal or window." >&2
        exit 2
        ;;
esac
case "$graph_source" in
    partition|rebuild) ;;
    *)
        echo "ERROR: ST_FOM_GRAPH_SOURCE must be partition or rebuild." >&2
        exit 2
        ;;
esac

thermal_root="${TEST_THERMAL_ROOT:-$project_root/../test_thermal}"
cmd2cond="${CMD2COND:-$thermal_root/bin/cmd2cond}"
meshgen_hex="${MESHGEN_HEX:-$thermal_root/bin/meshgen_hex}"
surf_dbc_all="${SURF_DBC_ALL:-$thermal_root/bin/surf_dbc_all}"
gedatsu_bin="${GEDATSU_BIN:-$thermal_root/submodule/monolis/submodule/gedatsu/bin}"
mesh2graph="${GEDATSU_MESH2GRAPH:-$gedatsu_bin/gedatsu_simple_mesh2graph_convertor}"
mesh_partitioner="${GEDATSU_MESH_PARTITIONER:-$gedatsu_bin/gedatsu_simple_mesh_partitioner}"
graph_partitioner="${GEDATSU_GRAPH_PARTITIONER:-$gedatsu_bin/gedatsu_nodal_graph_partitioner}"
bc_partitioner="${GEDATSU_BC_PARTITIONER:-$gedatsu_bin/gedatsu_bc_partitioner_R}"

for tool in \
    "$cmd2cond" \
    "$meshgen_hex" \
    "$surf_dbc_all" \
    "$mesh2graph" \
    "$mesh_partitioner" \
    "$graph_partitioner" \
    "$bc_partitioner"; do
    if [[ ! -x "$tool" ]]; then
        echo "ERROR: required executable was not found: $tool" >&2
        exit 1
    fi
done

if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: mpirun was not found in PATH." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Build normal ST-FOM.  No ROM/POD objects are linked by this Makefile.
# -----------------------------------------------------------------------------
if [[ "$rebuild_solver" -eq 1 || ! -x "$solver_dir/diff_ST_MPI" ]]; then
    echo "Building normal spatially partitioned ST-FOM..."
    (
        cd "$solver_dir"
        make -f Makefile_ST clean
        make -f Makefile_ST \
            ST_DEGREE="$degree" \
            ST_SLABS_PER_WINDOW="$slabs"
    )
fi

if [[ ! -x "$solver_dir/diff_ST_MPI" ]]; then
    echo "ERROR: build did not create $solver_dir/diff_ST_MPI" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Prepare ordinary spatial diffusion mesh and partition it only in space.
# -----------------------------------------------------------------------------
if [[ "$clean_result" -eq 1 ]]; then
    rm -rf "$result_dir"
fi

mkdir -p \
    "$result_dir" \
    "$result_dir/fem_solver_prm" \
    "$result_dir/calctime"

cp -f "$solver_dir/diff_ST_MPI" "$result_dir/diff_ST_MPI"
cd "$result_dir"

# Generate the normal ST-FOM calculation conditions.
rm -f cond.dat
"$cmd2cond" \
    "#num_ip_each_axis" int 1 "$num_ip_each_axis" \
    "#mat_epsilon" double 1 "$mat_epsilon" \
    "#mat_max_iter" int 1 "$mat_max_iter" \
    "#time_spacing" double 1 "$dt" \
    "#output_interval" int 1 "$output_interval" \
    "#finish_time" double 1 "$finish_time" \
    "#st_time_degree" int 1 "$degree" \
    "#st_slabs_per_window" int 1 "$slabs"

if [[ "$clean_result" -eq 1 || ! -f node.dat || ! -f elem.dat || ! -f D_bc.dat ]]; then
    "$meshgen_hex" "$e" "$e" "$e" "$ep" "$ep" "$ep"
    "$surf_dbc_all" 1 1.0
fi

if [[ "$clean_result" -eq 1 || ! -f graph.dat ]]; then
    "$mesh2graph" -i elem.dat -o graph.dat
fi

partition_ready=1
for ((rank=0; rank<np; rank++)); do
    for required in \
        "parted.0/node.dat.${rank}" \
        "parted.0/elem.dat.${rank}" \
        "parted.0/graph.dat.${rank}" \
        "parted.0/D_bc.dat.${rank}"; do
        [[ -f "$required" ]] || partition_ready=0
    done
done

if [[ "$clean_result" -eq 1 || "$partition_ready" -eq 0 ]]; then
    rm -rf parted.0
    echo "Partitioning ordinary spatial mesh into ${np} domains..."
    "$mesh_partitioner" -n "$np"
    "$graph_partitioner" -n "$np" -i graph.dat
    "$bc_partitioner" -n "$np" -i D_bc.dat -ig node.dat
fi

for ((rank=0; rank<np; rank++)); do
    for required in \
        "parted.0/node.dat.${rank}" \
        "parted.0/elem.dat.${rank}" \
        "parted.0/graph.dat.${rank}" \
        "parted.0/D_bc.dat.${rank}"; do
        if [[ ! -f "$required" ]]; then
            echo "ERROR: missing spatial partition file: $result_dir/$required" >&2
            exit 1
        fi
    done
done

echo
echo "Normal ST-FOM preparation completed."
echo "  case directory : $result_dir"
echo "  mesh            : ${e} x ${e} x ${e}, length=${ep}"
echo "  MPI domains     : $np"
echo "  time scheme     : dG(${degree})"
echo "  dt              : $dt"
echo "  finish time     : $finish_time"
echo "  total steps     : $total_steps"
echo "  slabs/window    : $slabs"
echo "  ST windows      : $num_windows"
echo "  solve mode      : $solve_mode"
echo "  graph source    : $graph_source"
echo "  ROM/POD         : disabled"
echo "  accuracy CSV    : $accuracy_csv"

if [[ "$prepare_only" -eq 1 ]]; then
    echo "ST_PREPARE_ONLY=1: stopping before mpirun."
    exit 0
fi

# -----------------------------------------------------------------------------
# Run normal ST-FOM.
# -----------------------------------------------------------------------------
log_file="${run_label}.log"
rm -f "$log_file" "$accuracy_csv"

command=(
    mpirun
    -np "$np"
    ./diff_ST_MPI
    ./
)

printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n\n'

set +e
env \
    OMP_NUM_THREADS="$omp_threads" \
    ST_FOM_TIME_DEGREE="$degree" \
    ST_FOM_SLABS_PER_WINDOW="$slabs" \
    ST_FOM_ACCURACY_CSV="$accuracy_csv" \
    ST_FOM_SOLVE_MODE="$solve_mode" \
    ST_FOM_GRAPH_SOURCE="$graph_source" \
    "${command[@]}" 2>&1 | tee "$log_file"
run_status=${PIPESTATUS[0]}
set -e

if [[ "$run_status" -ne 0 ]]; then
    echo "ERROR: ST-FOM solver exited with status $run_status." >&2
    echo "Check: $result_dir/$log_file" >&2
    exit "$run_status"
fi

if ! grep -q '\*\* Total time:' "$log_file"; then
    echo "ERROR: solver ended without the normal total-time message." >&2
    echo "Check: $result_dir/$log_file" >&2
    exit 1
fi

if [[ ! -s "$accuracy_csv" ]]; then
    echo "ERROR: accuracy CSV was not created: $result_dir/$accuracy_csv" >&2
    exit 1
fi

# The first line is the header; the last line is the final slab.
final_accuracy="$(tail -n 1 "$accuracy_csv")"

echo
echo "Normal spatially partitioned ST-FOM run completed."
echo "  run log      : $result_dir/$log_file"
echo "  accuracy CSV : $result_dir/$accuracy_csv"
echo "  final row    : $final_accuracy"
grep '\*\* Total time:' "$log_file" | tail -n 1
