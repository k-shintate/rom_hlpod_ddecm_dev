#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Partitioned windowed space-time FOM execution script
#
# Usage:
#   bash shell/diff/execution_ST_MPI.sh <e> <ep> <nm> <nd> <np> <pa> <st>
#
# Example:
#   bash shell/diff/execution_ST_MPI.sh 20 1.0 10 4 4 8 0
#
# Directory layout produced by this script:
#
#   result_diff/<nm>-<np>-<nd>/
#     cond.dat
#     node.dat
#     elem.dat
#     graph.dat
#     D_bc.dat
#     parted.0/               ordinary spatial partition
#     st/
#       st_meta.dat
#       graph.dat             symmetric ST partition graph
#       matrix_graph.dat      exact directed ST matrix graph
#       st_elem.dat or elem.dat
#       st_elem_map.dat
#       st_node_map.dat
#       parted.0/             partitioned ST data read by the solver
#
# The compiled C code must use:
#   ST_MPI_TOP_DIR  = "st"
#   ST_MPI_PART_DIR = "parted"
# so MONOLIS resolves the latter as st/parted.0.
# -----------------------------------------------------------------------------

if [[ -z "${BASH_VERSION:-}" ]]; then
    echo "ERROR: this script must be executed with Bash." >&2
    echo "Run: bash $0 ..." >&2
    exit 2
fi

set -eo pipefail

if [[ "$#" -ne 7 ]]; then
    echo "Usage: $0 <e> <ep> <nm> <nd> <np> <pa> <st>" >&2
    exit 2
fi

e="$1"
ep="$2"
nm="$3"
nd="$4"
np="$5"
pa="$6"
st="$7"

for name in e nm nd np pa; do
    value="${!name}"
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $name must be a non-negative integer: $value" >&2
        exit 2
    fi
done

if (( e <= 0 || nd <= 0 || np <= 0 )); then
    echo "ERROR: e, nd and np must be positive." >&2
    exit 2
fi

if ! [[ "$ep" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; then
    echo "ERROR: ep must be a positive real value: $ep" >&2
    exit 2
fi

# Enable nounset only after install.sh has been sourced.  Existing environment
# setup scripts often inspect optional variables without ${name:-} guards.
script_path="${BASH_SOURCE[0]}"
script_dir="$(cd "$(dirname "$script_path")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

if [[ ! -f "$project_root/shell/install.sh" ]]; then
    echo "ERROR: install.sh was not found:" >&2
    echo "  $project_root/shell/install.sh" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$project_root/shell/install.sh"
set -u

solver_dir="$project_root/solvers/diff"
result_dir="$project_root/result_diff/${nm}-${np}-${nd}"
thermal_root="${TEST_THERMAL_ROOT:-$project_root/../test_thermal}"

makefile="$solver_dir/Makefile_ST"
cmd2cond="$thermal_root/bin/cmd2cond"
meshgen_hex="$thermal_root/bin/meshgen_hex"
surf_dbc_all="$thermal_root/bin/surf_dbc_all"
gedatsu_bin="$thermal_root/submodule/monolis/submodule/gedatsu/bin"
mesh_to_graph="$gedatsu_bin/gedatsu_simple_mesh2graph_convertor"
mesh_partitioner="$gedatsu_bin/gedatsu_simple_mesh_partitioner"
graph_partitioner="$gedatsu_bin/gedatsu_nodal_graph_partitioner"
bc_partitioner="$gedatsu_bin/gedatsu_bc_partitioner_R"

st_degree="${ST_MPI_DEGREE:-1}"
st_window_slabs="${ST_MPI_WINDOW_SLABS:-4}"
omp_threads="${OMP_NUM_THREADS:-1}"
clean_result="${ST_CLEAN_RESULT:-1}"
rebuild_solver="${ST_FOM_REBUILD:-1}"
prepare_only="${ST_PREPARE_ONLY:-0}"
run_label="${ST_FOM_RUN_LABEL:-fom_st_partitioned_np${np}}"

# These names must match core_FOM_ST_MPI.h.
st_dir="$result_dir/st"
st_part_base="$st_dir/parted"
st_part_dir="$st_part_base.0"
spatial_part_dir="$result_dir/parted.0"

require_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: required file was not found:" >&2
        echo "  $path" >&2
        exit 1
    fi
}

require_executable() {
    local path="$1"
    if [[ ! -x "$path" ]]; then
        echo "ERROR: required executable was not found:" >&2
        echo "  $path" >&2
        exit 1
    fi
}

for executable in \
    "$cmd2cond" \
    "$meshgen_hex" \
    "$surf_dbc_all" \
    "$mesh_to_graph" \
    "$mesh_partitioner" \
    "$graph_partitioner" \
    "$bc_partitioner"; do
    require_executable "$executable"
done

require_file "$makefile"

# -----------------------------------------------------------------------------
# Build diff_ST_MPI, meshgen_ST and st_partition_postprocessor.
# -----------------------------------------------------------------------------

if [[ "$rebuild_solver" -eq 1 || \
      ! -x "$solver_dir/diff_ST_MPI" || \
      ! -x "$solver_dir/meshgen_ST" || \
      ! -x "$solver_dir/st_partition_postprocessor" ]]; then
    echo "Building partitioned ST tools and solver..."
    (
        cd "$solver_dir"
        make -f Makefile_ST clean
        make -f Makefile_ST
    )
fi

require_executable "$solver_dir/diff_ST_MPI"
require_executable "$solver_dir/meshgen_ST"
require_executable "$solver_dir/st_partition_postprocessor"

# -----------------------------------------------------------------------------
# Create the ordinary spatial problem.
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

if [[ "$clean_result" -eq 1 || \
      ! -f cond.dat || \
      ! -f node.dat || \
      ! -f elem.dat || \
      ! -f D_bc.dat || \
      ! -f graph.dat ]]; then

    echo "Generating ordinary spatial input..."

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

    "$meshgen_hex" "$e" "$e" "$e" "$ep" "$ep" "$ep"
    "$surf_dbc_all" 1 1.0
    "$mesh_to_graph" -i elem.dat -o graph.dat
fi

for input_file in cond.dat node.dat elem.dat D_bc.dat graph.dat; do
    require_file "$result_dir/$input_file"
done

# -----------------------------------------------------------------------------
# Generate the global one-window ST mesh under result_dir/st.
# meshgen_ST in the current project writes st_meta.dat, matrix_graph.dat,
# st_elem_map.dat and st_node_map.dat under the requested output directory.
# Depending on its local revision, connectivity may be named st_elem.dat or
# elem.dat; both are supported below.
# -----------------------------------------------------------------------------

echo "Generating global one-window ST mesh..."
rm -rf "$st_dir"
mkdir -p "$st_dir"

"$solver_dir/meshgen_ST" \
    "$st_degree" \
    "$st_window_slabs" \
    "$result_dir/node.dat" \
    "$result_dir/elem.dat" \
    "$st_dir"

st_meta="$st_dir/st_meta.dat"
st_matrix_graph="$st_dir/matrix_graph.dat"
st_partition_graph="$st_dir/graph.dat"
st_elem_map="$st_dir/st_elem_map.dat"
st_node_map="$st_dir/st_node_map.dat"

if [[ -f "$st_dir/st_elem.dat" ]]; then
    st_elem="$st_dir/st_elem.dat"
elif [[ -f "$st_dir/elem.dat" ]]; then
    st_elem="$st_dir/elem.dat"
else
    echo "ERROR: meshgen_ST created neither st/st_elem.dat nor st/elem.dat." >&2
    exit 1
fi

for st_file in \
    "$st_meta" \
    "$st_matrix_graph" \
    "$st_partition_graph" \
    "$st_elem" \
    "$st_elem_map" \
    "$st_node_map"; do
    require_file "$st_file"
done

# -----------------------------------------------------------------------------
# Partition the ordinary spatial mesh.  The postprocessor uses the spatial
# ownership information to assign all temporal replicas of a space node to the
# same rank.
# -----------------------------------------------------------------------------

echo "Partitioning the ordinary spatial mesh into $np domains..."
rm -rf "$spatial_part_dir"

"$mesh_partitioner" -n "$np"
"$graph_partitioner" -n "$np" -i graph.dat
"$bc_partitioner" -n "$np" -i D_bc.dat -ig node.dat

# st_partition_postprocessor expects a spatial owner map.  Preserve the
# project's existing owner-file convention, but locate it robustly.
if [[ -n "${ST_OWNER_FILE:-}" ]]; then
    owner_file="$ST_OWNER_FILE"
    if [[ "$owner_file" != /* ]]; then
        owner_file="$result_dir/$owner_file"
    fi
else
    owner_file=""
    for candidate in \
        "$result_dir/owner.dat" \
        "$spatial_part_dir/owner.dat" \
        "$result_dir/node_owner.dat" \
        "$spatial_part_dir/node_owner.dat"; do
        if [[ -f "$candidate" ]]; then
            owner_file="$candidate"
            break
        fi
    done
fi

if [[ -z "$owner_file" || ! -f "$owner_file" ]]; then
    echo "ERROR: the spatial owner map required by st_partition_postprocessor" >&2
    echo "was not found." >&2
    echo "Set it explicitly, for example:" >&2
    echo "  ST_OWNER_FILE=/absolute/path/to/owner.dat ..." >&2
    echo "Searched:" >&2
    echo "  $result_dir/owner.dat" >&2
    echo "  $spatial_part_dir/owner.dat" >&2
    echo "  $result_dir/node_owner.dat" >&2
    echo "  $spatial_part_dir/node_owner.dat" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Build partitioned ST graph/connectivity/maps.
# The output base is st/parted, and the program creates st/parted.0.
# -----------------------------------------------------------------------------

echo "Postprocessing the ST partition for $np domains..."
rm -rf "$st_part_dir"

# Files read by ST_partition_read() through MONOLIS parted-file naming.
for rank in $(seq 0 $((np - 1))); do
    require_file "$st_part_dir/graph.dat.$rank"
    require_file "$st_part_dir/st_node_map.dat.$rank"
    require_file "$st_part_dir/elem.dat.$rank"
    require_file "$st_part_dir/st_elem_map.dat.$rank"
done

cat <<SUMMARY
ST preprocessing completed.
  result directory : $result_dir
  spatial partition: $spatial_part_dir
  ST metadata      : $st_meta
  ST matrix graph  : $st_matrix_graph
  ST owner map     : $owner_file
  ST partition     : $st_part_dir
  degree/slabs     : $st_degree/$st_window_slabs
  MPI domains      : $np
SUMMARY

if [[ "$prepare_only" -eq 1 ]]; then
    echo "ST_PREPARE_ONLY=1: stopping before mpirun."
    exit 0
fi

# -----------------------------------------------------------------------------
# Execute the partitioned windowed ST FOM.
# -----------------------------------------------------------------------------

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
echo "Partitioned FOM-ST run completed."
echo "Run log: $result_dir/$log_file"
grep '\*\* Total time:' "$log_file" | tail -n 1
