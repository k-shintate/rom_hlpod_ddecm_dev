#!/usr/bin/env bash

# Fluid all-at-once space-time FOM launcher with the same spatial
# GEDATSU/MONOLIS decomposition used by the convection-diffusion ST-DDROM.
#
# Spatial mesh: ordinary GEDATSU/MONOLIS partition.
# ST unknowns : all slab/time/physical DOFs are a block on each spatial node.
# No space-time mesh partition, st_meta.dat, owner.dat, or ST mesh
# postprocessor is used.
#
# Future fluid ST-DDROM metanode:
#   (time window, MPI spatial rank)

if [[ -z "${BASH_VERSION:-}" ]]; then
    exec /usr/bin/env bash "$0" "$@"
fi

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2; exit "$rc"' ERR

usage() {
    cat >&2 <<'USAGE'
Usage:
  execution_ST_GEDATSU_MPI.sh <e> <ep> <nm> <nd> <np> <pa> <st>

Arguments:
  e   : elements in each spatial direction
  ep  : domain length
  nm  : result-directory tag
  nd  : legacy/result-directory tag
  np  : ordinary spatial GEDATSU partitions and MPI process count
  pa  : legacy option
  st  : legacy solver option

Environment:
  FLUID_ST_TIME_DEGREE=0|1|2               default: 1
  FLUID_ST_SLABS_PER_WINDOW=<int>          default: 4
  FLUID_ST_SOLVE_MODE=all_at_once          default: all_at_once
  FLUID_ST_JOLD_MODE=fd|verify|analytic    default: fd
  FLUID_ST_REPARTITION=0|1                 default: 1
  FLUID_ST_MAKEFILE=<name>                 default: Makefile_ST_cavity
  FLUID_ST_TARGET=<target>                 default: fluid_sups_st_fom
  FLUID_ST_EXECUTABLE=<name>               default: fluid_sups_st_fom
  GEDATSU_BIN=<directory>                  optional override
  TEST_THERMAL_ROOT=<directory>            optional override

Example:
  bash shell/fluid_sups_hrom_cavity/execution_ST_GEDATSU_MPI.sh \
      20 1 10 1 4 0 3
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
              -d "$candidate/solvers/fluid_sups" ]]; then
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

result_dir="$project_root/result_fluid_sups_cavity/${nm}-${np}-${nd}"
solver_dir="$project_root/solvers/fluid_sups"

thermal_root="${TEST_THERMAL_ROOT:-$project_root/../test_thermal}"
gedatsu_bin="${GEDATSU_BIN:-$thermal_root/submodule/monolis/submodule/gedatsu/bin}"

mesh2graph="${GEDATSU_MESH2GRAPH:-$gedatsu_bin/gedatsu_simple_mesh2graph_convertor}"
mesh_partitioner="${GEDATSU_MESH_PARTITIONER:-$gedatsu_bin/gedatsu_simple_mesh_partitioner}"
graph_partitioner="${GEDATSU_GRAPH_PARTITIONER:-$gedatsu_bin/gedatsu_nodal_graph_partitioner}"
bc_partitioner="${GEDATSU_BC_PARTITIONER:-$gedatsu_bin/gedatsu_bc_partitioner_R}"

for tool in \
    "$mesh2graph" \
    "$mesh_partitioner" \
    "$graph_partitioner" \
    "$bc_partitioner"; do
    if [[ ! -x "$tool" ]]; then
        echo "ERROR: required GEDATSU executable was not found: $tool" >&2
        exit 1
    fi
done

export FLUID_ST_INITIAL_CONDITION="${FLUID_ST_INITIAL_CONDITION:-cavity}"
export FLUID_ST_TIME_DEGREE="${FLUID_ST_TIME_DEGREE:-1}"
export FLUID_ST_SLABS_PER_WINDOW="${FLUID_ST_SLABS_PER_WINDOW:-4}"
export FLUID_ST_SOLVE_MODE="${FLUID_ST_SOLVE_MODE:-all_at_once}"
export FLUID_ST_WRITE_WINDOW_BINARY="${FLUID_ST_WRITE_WINDOW_BINARY:-1}"

# Keep MPI validation independent from the analytic-J_old development path.
export FLUID_ST_JOLD_MODE="${FLUID_ST_JOLD_MODE:-fd}"

repartition="${FLUID_ST_REPARTITION:-1}"
if [[ "$repartition" != 0 && "$repartition" != 1 ]]; then
    echo "ERROR: FLUID_ST_REPARTITION must be 0 or 1: $repartition" >&2
    exit 2
fi

echo
echo "============================================================"
echo "Fluid ST-FOM ordinary spatial GEDATSU partition"
echo "  mesh divisions        : $e"
echo "  domain length         : $ep"
echo "  spatial MPI ranks     : $np"
echo "  dG degree             : $FLUID_ST_TIME_DEGREE"
echo "  slabs/window          : $FLUID_ST_SLABS_PER_WINDOW"
echo "  solve mode            : $FLUID_ST_SOLVE_MODE"
echo "  J_old mode            : $FLUID_ST_JOLD_MODE"
echo "  result directory      : $result_dir"
echo "============================================================"

# ----------------------------------------------------------------------
# 1. Generate the ordinary, unpartitioned fluid cavity mesh/BC using the
#    existing project mesh workflow.
# ----------------------------------------------------------------------
bash "$project_root/shell/fluid_sups_hrom_cavity/meshgen.sh" \
    "$e" "$ep" "$nm" "$nd" "$np" "$pa"

if [[ ! -d "$result_dir" ]]; then
    echo "ERROR: meshgen did not create result directory: $result_dir" >&2
    exit 1
fi

cd "$result_dir"

for required in node.dat elem.dat D_bc_v.dat; do
    if [[ ! -f "$required" ]]; then
        echo "ERROR: ordinary global fluid mesh/BC file is missing: $result_dir/$required" >&2
        echo "The GEDATSU decomposition intentionally starts from node.dat, elem.dat, and D_bc_v.dat." >&2
        exit 1
    fi
done

# ----------------------------------------------------------------------
# 2. Build the ordinary spatial nodal graph.
# ----------------------------------------------------------------------
echo
echo "Building ordinary spatial graph.dat from elem.dat..."
rm -f graph.dat
"$mesh2graph" -i elem.dat -o graph.dat

if [[ ! -f graph.dat ]]; then
    echo "ERROR: graph.dat was not generated." >&2
    exit 1
fi

# ----------------------------------------------------------------------
# 3. Partition ONLY the ordinary spatial mesh.
#
#    This is the same design as the supplied convection-diffusion launcher:
#      node/elem  -> simple mesh partitioner
#      graph      -> nodal graph partitioner
#      D_bc_v     -> real-valued BC partitioner
#
#    Time windows/slabs are not GEDATSU partitions.
# ----------------------------------------------------------------------
partition_marker="parted.0/.fluid_st_spatial_np"
partition_ready=1

if [[ "$repartition" -eq 1 ]]; then
    partition_ready=0
elif [[ ! -f "$partition_marker" ]] || \
     [[ "$(cat "$partition_marker" 2>/dev/null || true)" != "$np" ]]; then
    partition_ready=0
else
    for ((rank=0; rank<np; rank++)); do
        for required in \
            "parted.0/node.dat.${rank}" \
            "parted.0/elem.dat.${rank}" \
            "parted.0/graph.dat.${rank}" \
            "parted.0/D_bc_v.dat.${rank}"; do
            [[ -f "$required" ]] || partition_ready=0
        done
    done
fi

if [[ "$partition_ready" -eq 0 ]]; then
    rm -rf parted.0

    echo
    echo "Partitioning ordinary fluid spatial mesh into ${np} GEDATSU domains..."

    "$mesh_partitioner" -n "$np"
    "$graph_partitioner" -n "$np" -i graph.dat
    "$bc_partitioner" -n "$np" -i D_bc_v.dat -ig node.dat

    mkdir -p parted.0
    printf '%s\n' "$np" > "$partition_marker"
fi

# ----------------------------------------------------------------------
# 4. Strict partition-file validation.
# ----------------------------------------------------------------------
for ((rank=0; rank<np; rank++)); do
    for required in \
        "parted.0/node.dat.${rank}" \
        "parted.0/elem.dat.${rank}" \
        "parted.0/graph.dat.${rank}" \
        "parted.0/D_bc_v.dat.${rank}"; do
        if [[ ! -f "$required" ]]; then
            echo "ERROR: missing spatial partition file: $result_dir/$required" >&2
            exit 1
        fi
    done
done

echo
echo "Spatial GEDATSU partition files are ready:"
echo "  $result_dir/parted.0/node.dat.<rank>"
echo "  $result_dir/parted.0/elem.dat.<rank>"
echo "  $result_dir/parted.0/graph.dat.<rank>"
echo "  $result_dir/parted.0/D_bc_v.dat.<rank>"
echo
echo "No time-window partition files are generated."

# ----------------------------------------------------------------------
# 5. Build and run through the existing ST execution path.
#    The executable reads the rank-local files through
#    monolis_get_global_input_file_name(), so its physical MONOLIS
#    communicator is exactly this ordinary spatial decomposition.
# ----------------------------------------------------------------------
cd "$project_root"

bash "$project_root/shell/fluid_sups_hrom_cavity/execution_ST.sh" \
    "$e" "$ep" "$nm" "$nd" "$np" "$pa" "$st"
