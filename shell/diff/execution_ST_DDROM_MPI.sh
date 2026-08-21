#!/usr/bin/env bash

# Space-time domain-decomposed ROM launcher.
#
# Spatial mesh: ordinary GEDATSU/MONOLIS partition.
# ROM metanode: (time window, MPI spatial rank).
# No st_meta.dat, owner.dat, or ST mesh postprocessor is used.

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  execution_ST_DDROM_MPI.sh <e> <ep> <nm> <nd> <np> <pa> <st>

Arguments:
  e   : elements in each spatial direction
  ep  : domain length in each spatial direction
  nm  : default maximum local POD modes and result-directory tag
  nd  : legacy option / result-directory tag
  np  : MPI spatial partitions and process count
  pa  : legacy option
  st  : legacy solver option

Environment:
  ST_DDROM_ACTION=full|fom|collect|offline|online   default: full
  ST_DDROM_MAX_MODES=<int>                         default: nm
  ST_DDROM_EPSILON=<real>                          default: 1.0
  ST_DDROM_NUM_SNAPSHOTS=<int>                     default: 1
  ST_DDROM_SNAPSHOT_ID=<int>                       default: 0
  ST_DDROM_VALIDATE=0|1                            default: 1
  ST_DDROM_WRITE_OUTPUT=0|1                        default: 0
  ST_DDROM_REDUCED_SOLVER=window|all_at_once       default: window
  ST_DDROM_REDUCED_EPSILON=<real>                  default: 1.0e-10
  ST_DDROM_CHECK_FOM_RESIDUAL=0|1                  primitive residual diagnostic, default: 1
  ST_DDROM_STRICT_FOM_RESIDUAL=0|1                 reject on primitive free residual, default: 0
  ST_DDROM_FOM_RESIDUAL_TOL=<real>                 comparison threshold, default: 1.0e-6
  ST_DDROM_REBUILD_SNAPSHOTS=0|1                   default: 1 for full
  ST_DDROM_REBUILD_BASIS=0|1                       default: 1 for full/offline
  ST_CLEAN_RESULT=0|1
  ST_FOM_REBUILD=0|1                               default: 1

Example:
  bash shell/diff/execution_ST_DDROM_MPI.sh 20 1.0 10 4 4 8 0
USAGE
}

if [[ "$#" -ne 7 ]]; then
    usage
    exit 2
fi

e="$1"
ep="$2"
nm="$3"
nd="$4"
np="$5"
pa="$6"
st="$7"

for pair in "e:$e" "nm:$nm" "nd:$nd" "np:$np" "pa:$pa"; do
    name="${pair%%:*}"
    value="${pair#*:}"
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: ${name} must be a nonnegative integer: ${value}" >&2
        exit 2
    fi
done

if (( e <= 0 || np <= 0 || nm <= 0 )); then
    echo "ERROR: e, np, and nm must be positive." >&2
    exit 2
fi

if ! [[ "$ep" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
    echo "ERROR: ep must be a positive number: ${ep}" >&2
    exit 2
fi

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

# shellcheck source=/dev/null
source "$project_root/shell/install.sh"

solver_dir="$project_root/solvers/diff"
result_dir="$project_root/result_diff/${nm}-${np}-${nd}"
action="${ST_DDROM_ACTION:-full}"
max_modes="${ST_DDROM_MAX_MODES:-$nm}"
epsilon="${ST_DDROM_EPSILON:-1.0}"
num_snapshots="${ST_DDROM_NUM_SNAPSHOTS:-1}"
snapshot_id="${ST_DDROM_SNAPSHOT_ID:-0}"
validate="${ST_DDROM_VALIDATE:-1}"
write_output="${ST_DDROM_WRITE_OUTPUT:-0}"
reduced_solver="${ST_DDROM_REDUCED_SOLVER:-window}"
reduced_epsilon="${ST_DDROM_REDUCED_EPSILON:-1.0e-10}"
reduced_max_iter="${ST_DDROM_REDUCED_MAX_ITER:-10000}"
fom_residual_check="${ST_DDROM_CHECK_FOM_RESIDUAL:-1}"
fom_residual_tolerance="${ST_DDROM_FOM_RESIDUAL_TOL:-1.0e-6}"
rebuild_solver="${ST_FOM_REBUILD:-1}"
omp_threads="${OMP_NUM_THREADS:-1}"

case "$action" in
    full|fom|collect)
        clean_result="${ST_CLEAN_RESULT:-1}"
        ;;
    offline|online)
        clean_result="${ST_CLEAN_RESULT:-0}"
        ;;
    *)
        echo "ERROR: unknown ST_DDROM_ACTION=$action" >&2
        exit 2
        ;;
esac

case "$action" in
    full)
        strict_fom_residual="${ST_DDROM_STRICT_FOM_RESIDUAL:-0}"
        rebuild_snapshots="${ST_DDROM_REBUILD_SNAPSHOTS:-1}"
        rebuild_basis="${ST_DDROM_REBUILD_BASIS:-1}"
        ;;
    collect)
        strict_fom_residual="${ST_DDROM_STRICT_FOM_RESIDUAL:-0}"
        rebuild_snapshots="${ST_DDROM_REBUILD_SNAPSHOTS:-0}"
        rebuild_basis="${ST_DDROM_REBUILD_BASIS:-0}"
        ;;
    offline)
        strict_fom_residual="${ST_DDROM_STRICT_FOM_RESIDUAL:-0}"
        rebuild_snapshots="${ST_DDROM_REBUILD_SNAPSHOTS:-0}"
        rebuild_basis="${ST_DDROM_REBUILD_BASIS:-1}"
        ;;
    *)
        strict_fom_residual="${ST_DDROM_STRICT_FOM_RESIDUAL:-0}"
        rebuild_snapshots="${ST_DDROM_REBUILD_SNAPSHOTS:-0}"
        rebuild_basis="${ST_DDROM_REBUILD_BASIS:-0}"
        ;;
esac

for boolean_name in \
    clean_result \
    rebuild_solver \
    fom_residual_check \
    strict_fom_residual \
    rebuild_snapshots \
    rebuild_basis; do
    boolean_value="${!boolean_name}"
    if [[ "$boolean_value" != 0 && "$boolean_value" != 1 ]]; then
        echo "ERROR: ${boolean_name} must be 0 or 1: ${boolean_value}" >&2
        exit 2
    fi
done

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

if [[ "$rebuild_solver" -eq 1 || ! -x "$result_dir/hlpod_diff_st_ddrom" ]]; then
    echo "Building ST-DDROM executable..."
    (
        cd "$solver_dir"
        make -f Makefile_ST clean
        make -f Makefile_ST
    )

    if [[ ! -x "$solver_dir/hlpod_diff_st_ddrom" ]]; then
        echo "ERROR: build did not create $solver_dir/hlpod_diff_st_ddrom" >&2
        exit 1
    fi
fi

if [[ "$clean_result" -eq 1 ]]; then
    rm -rf "$result_dir"
fi

mkdir -p \
    "$result_dir" \
    "$result_dir/fem_solver_prm" \
    "$result_dir/calctime"

cp -f "$solver_dir/hlpod_diff_st_ddrom" "$result_dir/hlpod_diff_st_ddrom"
cd "$result_dir"

if [[ ! -f cond.dat ]]; then
    rm -f cond.dat rom_cond.dat

    "$cmd2cond" \
        "#snapshot_interval" int 1 1 \
        "#rom_finish_time" double 1 4.0 \
        "#rom_output_interval" int 1 1
    mv cond.dat rom_cond.dat

    "$cmd2cond" \
        "#time_spacing" double 1 0.01 \
        "#output_interval" int 1 1 \
        "#finish_time" double 1 1.0
fi

if [[ ! -f node.dat || ! -f elem.dat ]]; then
    "$meshgen_hex" "$e" "$e" "$e" "$ep" "$ep" "$ep"
    "$surf_dbc_all" 1 1.0
fi

if [[ ! -f graph.dat ]]; then
    "$mesh2graph" -i elem.dat -o graph.dat
fi

partition_ready=1
for required in \
    "parted.0/node.dat.0" \
    "parted.0/elem.dat.0" \
    "parted.0/graph.dat.0" \
    "parted.0/D_bc.dat.0"; do
    [[ -f "$required" ]] || partition_ready=0
done

if [[ "$partition_ready" -eq 0 ]]; then
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

snapshot_dir="$result_dir/stddrom_snapshots"
basis_dir="$result_dir/stddrom_basis_r${max_modes}"

if [[ "$rebuild_snapshots" -eq 1 ]]; then
    if [[ "$action" == full ]]; then
        echo "Removing old ST-DDROM snapshots: $snapshot_dir"
        rm -rf "$snapshot_dir"
    elif [[ "$action" == collect ]]; then
        snapshot_tag="$(printf '%04d' "$snapshot_id")"
        echo "Removing old snapshot ID ${snapshot_id} from: $snapshot_dir"
        rm -f "$snapshot_dir/stddrom_snapshot_${snapshot_tag}_rank_"*.bin
    fi
fi

if [[ "$rebuild_basis" -eq 1 && \
      ( "$action" == full || "$action" == offline ) ]]; then
    echo "Removing old ST-DDROM basis/operators: $basis_dir"
    rm -rf "$basis_dir"
fi

run_mode() {
    local mode="$1"
    local label="$2"
    local log_file="${label}.log"

    echo
    echo "Running ST-DDROM mode: ${mode}"

    env \
        OMP_NUM_THREADS="$omp_threads" \
        OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" \
        ST_DDROM_MODE="$mode" \
        ST_DDROM_SNAPSHOT_ID="$snapshot_id" \
        ST_DDROM_NUM_SNAPSHOTS="$num_snapshots" \
        ST_DDROM_MAX_MODES="$max_modes" \
        ST_DDROM_EPSILON="$epsilon" \
        ST_DDROM_VALIDATE="$validate" \
        ST_DDROM_WRITE_OUTPUT="$write_output" \
        ST_DDROM_REDUCED_SOLVER="$reduced_solver" \
        ST_DDROM_REDUCED_EPSILON="$reduced_epsilon" \
        ST_DDROM_REDUCED_MAX_ITER="$reduced_max_iter" \
        ST_DDROM_CHECK_FOM_RESIDUAL="$fom_residual_check" \
        ST_DDROM_STRICT_FOM_RESIDUAL="$strict_fom_residual" \
        ST_DDROM_FOM_RESIDUAL_TOL="$fom_residual_tolerance" \
        ST_DDROM_SNAPSHOT_DIR="$snapshot_dir" \
        ST_DDROM_BASIS_DIR="$basis_dir" \
        mpirun -np "$np" \
            ./hlpod_diff_st_ddrom ./ \
            -nd "$nd" \
            -nm "$nm" \
            -pa "$pa" \
            -st "$st" \
        2>&1 | tee "$log_file"

    if ! grep -q '\*\* Total time:' "$log_file"; then
        echo "ERROR: ${mode} ended without the normal total-time message." >&2
        exit 1
    fi
}

case "$action" in
    full)
        run_mode collect "stddrom_collect_${snapshot_id}_np${np}"
        run_mode offline "stddrom_offline_np${np}_r${max_modes}"
        run_mode online "stddrom_online_np${np}_r${max_modes}"
        ;;
    fom)
        run_mode fom "stddrom_fom_np${np}"
        ;;
    collect)
        run_mode collect "stddrom_collect_${snapshot_id}_np${np}"
        ;;
    offline)
        run_mode offline "stddrom_offline_np${np}_r${max_modes}"
        ;;
    online)
        run_mode online "stddrom_online_np${np}_r${max_modes}"
        ;;
esac

echo
echo "ST-DDROM action completed: $action"
echo "Case directory : $result_dir"
echo "Snapshots      : $snapshot_dir"
echo "Basis/operators: $basis_dir"
