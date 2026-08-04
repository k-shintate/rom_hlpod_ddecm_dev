#!/usr/bin/env bash

# Spatially partitioned space-time block solver launcher.
# This path does not use meshgen_ST, st_meta.dat, owner.dat, or
# st_partition_postprocessor.  Only the ordinary spatial mesh is partitioned.

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
    cat >&2 <<'USAGE'
Usage:
  execution_ST_spaceblock_MPI.sh <e> <ep> <nm> <nd> <np> <pa> <st>

Arguments:
  e   : number of hexahedral elements in each spatial direction
  ep  : domain length in each spatial direction
  nm  : result-directory tag / legacy ROM option
  nd  : legacy domain option (not used for spatial partitioning)
  np  : MPI process count and spatial partition count
  pa  : legacy POD-threshold option
  st  : legacy solver option

Example:
  bash shell/diff/execution_ST_spaceblock_MPI.sh 20 1.0 10 4 4 8 0
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

if (( e <= 0 || np <= 0 )); then
    echo "ERROR: e and np must be positive." >&2
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

rebuild_solver="${ST_FOM_REBUILD:-1}"
clean_result="${ST_CLEAN_RESULT:-1}"
prepare_only="${ST_PREPARE_ONLY:-0}"
omp_threads="${OMP_NUM_THREADS:-1}"
run_label="${ST_FOM_RUN_LABEL:-fom_st_spaceblock_np${np}}"

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

# ---------------------------------------------------------------------------
# Build the solver.  Makefile_ST must select the space-block main and include
# core_FOM_ST_spaceblock.c in its object/source list.
# ---------------------------------------------------------------------------

if [[ "$rebuild_solver" -eq 1 || ! -x "$result_dir/diff_ST_MPI" ]]; then
    echo "Building spatially partitioned ST block solver..."
    (
        cd "$solver_dir"
        make -f Makefile_ST clean
        make -f Makefile_ST
    )

    if [[ ! -x "$solver_dir/diff_ST_MPI" ]]; then
        echo "ERROR: build did not create $solver_dir/diff_ST_MPI" >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Create the ordinary spatial case and partition only the spatial mesh.
# ---------------------------------------------------------------------------

if [[ "$clean_result" -eq 1 ]]; then
    rm -rf "$result_dir"
fi

mkdir -p \
    "$result_dir" \
    "$result_dir/fem_solver_prm" \
    "$result_dir/calctime"

cp -f "$solver_dir/diff_ST_MPI" "$result_dir/diff_ST_MPI"

cd "$result_dir"

if [[ "$clean_result" -eq 1 || ! -f cond.dat ]]; then
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

if [[ "$clean_result" -eq 1 || ! -f node.dat || ! -f elem.dat ]]; then
    "$meshgen_hex" "$e" "$e" "$e" "$ep" "$ep" "$ep"
    "$surf_dbc_all" 1 1.0
fi

if [[ "$clean_result" -eq 1 || ! -f graph.dat ]]; then
    "$mesh2graph" -i elem.dat -o graph.dat
fi

# Recreate partition data whenever the case is cleaned or rank-0 files are absent.
partition_ready=1
for required in \
    "parted.0/node.dat.0" \
    "parted.0/elem.dat.0" \
    "parted.0/graph.dat.0" \
    "parted.0/D_bc.dat.0"; do
    [[ -f "$required" ]] || partition_ready=0
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
 echo "Spatial-only partition preparation completed."
echo "  case directory : $result_dir"
echo "  partition data : $result_dir/parted.0"
echo "  MPI domains    : $np"
echo "  ST degree      : compile-time STSB_DEGREE (default 1)"
echo "  window slabs   : compile-time STSB_NUM_WINDOW_SLABS (default 4)"
echo "  no st_meta.dat / owner.dat / ST postprocessor is used"

if [[ "$prepare_only" -eq 1 ]]; then
    echo "ST_PREPARE_ONLY=1: stopping before mpirun."
    exit 0
fi

# ---------------------------------------------------------------------------
# Run the spatially partitioned ST block solver.
# ---------------------------------------------------------------------------

cat > gdb_cmd <<GDBEOF
run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}
backtrace
quit
GDBEOF

log_file="${run_label}.log"
rm -f "$log_file"

command=(
    mpirun
    -np "$np"
    ./diff_ST_MPI
    ./
    -nd "$nd"
    -nm "$nm"
    -pa "$pa"
    -st "$st"
)

printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'

env OMP_NUM_THREADS="$omp_threads" \
    "${command[@]}" 2>&1 | tee "$log_file"

if ! grep -q '\*\* Total time:' "$log_file"; then
    echo "ERROR: solver ended without the normal total-time message." >&2
    echo "Check: $result_dir/$log_file" >&2
    exit 1
fi

echo
echo "Spatially partitioned ST block run completed."
echo "Run log: $result_dir/$log_file"
grep '\*\* Total time:' "$log_file" | tail -n 1
