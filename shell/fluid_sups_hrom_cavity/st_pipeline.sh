#!/usr/bin/env bash
set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]}:${LINENO}: ${BASH_COMMAND} (exit=${rc})" >&2; exit "$rc"' ERR

usage() {
    cat >&2 <<'EOF'
usage: st_pipeline.sh <prepare|fom|offline|online|verify|validated|validated_full|status>

  prepare        generate and partition the cavity case
  fom            run ST-FOM and collect window binaries
  offline        build ST-DDROM POD basis/operators from existing FOM windows
  online         run ST-DDROM from an existing basis
  verify         verify the validated online result
  validated      reuse existing mesh + FOM windows: offline -> online -> verify
  validated_full end-to-end: prepare -> FOM -> offline -> online -> verify
  status         print which prerequisite files are present
EOF
    exit 2
}
[[ "$#" -eq 1 ]] || usage
action="$1"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
# shellcheck source=/dev/null
source "${script_dir}/st_case.env"
# shellcheck source=/dev/null
source "${project_root}/shell/install.sh"

E="$FLUID_ST_MESH_DIVISIONS"; EP="$FLUID_ST_DOMAIN_SIZE"
NM="$FLUID_ST_RESULT_MODES_TAG"; ND="$FLUID_ST_RESULT_DD_TAG"
NP="$FLUID_ST_MPI_RANKS"; PA="$FLUID_ST_PA_TAG"; ST="$FLUID_ST_STAGE_TAG"

default_result_dir="${project_root}/result_fluid_sups_cavity/${NM}-${NP}-${ND}"
result_dir="${FLUID_ST_RESULT_DIR:-$default_result_dir}"
solver_dir="${project_root}/solvers/fluid_sups"
test_thermal="${TEST_THERMAL_ROOT:-${project_root}/../test_thermal}"
gedatsu_bin="${GEDATSU_BIN:-${test_thermal}/submodule/monolis/submodule/gedatsu/bin}"
makefile="${FLUID_ST_MAKEFILE:-Makefile_ST_cavity}"

window_dir() {
    printf '%s\n' "${FLUID_ST_WINDOW_DIR:-${result_dir}/fluid_st_windows_s${FLUID_ST_SLABS_PER_WINDOW}}"
}

ddrom_snapshot_dir() {
    printf '%s\n' "${FLUID_ST_DDROM_SNAPSHOT_DIR:-$(window_dir)}"
}

ddrom_basis_dir() {
    printf '%s\n' "${FLUID_ST_DDROM_BASIS_DIR:-${result_dir}/fluid_st_ddrom_basis}"
}

case_inputs_ready() {
    local rank f
    [[ -f "${result_dir}/cond.dat" ]] || return 1
    [[ -d "${result_dir}/parted.0" ]] || return 1
    for ((rank=0; rank<NP; rank++)); do
        for f in node.dat elem.dat graph.dat D_bc_v.dat; do
            [[ -f "${result_dir}/parted.0/${f}.${rank}" ]] || return 1
        done
    done
    return 0
}

snapshot_count_for_rank() {
    local rank="$1"
    local dir count name
    dir="$(ddrom_snapshot_dir)"
    count=0
    while :; do
        printf -v name 'fluid_st_window_%04d_rank_%06d.bin' "$count" "$rank"
        [[ -f "${dir}/${name}" ]] || break
        count=$((count + 1))
    done
    printf '%d\n' "$count"
}

fom_snapshots_ready() {
    local rank count
    local required="$FLUID_ST_DDROM_TRAINING_WINDOWS"
    for ((rank=0; rank<NP; rank++)); do
        count="$(snapshot_count_for_rank "$rank")"
        (( count >= required )) || return 1
    done
    return 0
}

basis_ready() {
    local rank dir name
    dir="$(ddrom_basis_dir)"
    for ((rank=0; rank<NP; rank++)); do
        printf -v name 'stddrom_basis_rank_%06d.bin' "$rank"
        [[ -f "${dir}/${name}" ]] || return 1
    done
    return 0
}

print_status() {
    local rank count
    echo "Fluid ST validated-case status"
    echo "  project root      : $project_root"
    echo "  result directory  : $result_dir"
    echo "  expected MPI ranks: $NP"
    echo "  dG/slabs          : ${FLUID_ST_TIME_DEGREE}/${FLUID_ST_SLABS_PER_WINDOW}"
    echo "  training/online   : ${FLUID_ST_DDROM_TRAINING_WINDOWS}/${FLUID_ST_DDROM_ONLINE_WINDOWS} windows"
    echo
    if case_inputs_ready; then
        echo "  case mesh/BC      : READY"
    else
        echo "  case mesh/BC      : MISSING"
    fi
    echo "  snapshot directory: $(ddrom_snapshot_dir)"
    for ((rank=0; rank<NP; rank++)); do
        count="$(snapshot_count_for_rank "$rank")"
        printf '    rank %d snapshots: %d\n' "$rank" "$count"
    done
    if fom_snapshots_ready; then
        echo "  FOM snapshots     : READY"
    else
        echo "  FOM snapshots     : MISSING/INSUFFICIENT"
    fi
    if basis_ready; then
        echo "  DDROM basis       : READY"
    else
        echo "  DDROM basis       : MISSING"
    fi
}

require_case_inputs() {
    if case_inputs_ready; then
        return 0
    fi
    cat >&2 <<EOF
ERROR: the ST cavity case is not prepared.

Expected:
  ${result_dir}/cond.dat
  ${result_dir}/parted.0/node.dat.0 ... node.dat.$((NP-1))
  ${result_dir}/parted.0/elem.dat.*
  ${result_dir}/parted.0/graph.dat.*
  ${result_dir}/parted.0/D_bc_v.dat.*

This is an input-data problem, not an ST-DDROM numerical failure.

Create the case with:
  bash "${script_dir}/st_pipeline.sh" prepare

For a complete FOM -> ROM regression from scratch:
  bash "${script_dir}/st_pipeline.sh" validated_full

To point at an existing prepared result directory:
  export FLUID_ST_RESULT_DIR=/absolute/path/to/result_fluid_sups_cavity/10-4-1
EOF
    exit 3
}

require_fom_snapshots() {
    if fom_snapshots_ready; then
        return 0
    fi
    local rank count
    echo "ERROR: insufficient FOM ST window snapshots." >&2
    echo "Expected at least ${FLUID_ST_DDROM_TRAINING_WINDOWS} consecutive windows on every rank in:" >&2
    echo "  $(ddrom_snapshot_dir)" >&2
    for ((rank=0; rank<NP; rank++)); do
        count="$(snapshot_count_for_rank "$rank")"
        echo "  rank ${rank}: ${count}" >&2
    done
    cat >&2 <<EOF

Collect them with:
  bash "${script_dir}/st_pipeline.sh" fom

Or run the complete validated workflow:
  bash "${script_dir}/st_pipeline.sh" validated_full
EOF
    exit 4
}

require_basis() {
    if basis_ready; then
        return 0
    fi
    cat >&2 <<EOF
ERROR: ST-DDROM basis files are missing from:
  $(ddrom_basis_dir)

Build them with:
  bash "${script_dir}/st_pipeline.sh" offline
EOF
    exit 5
}

build_target() {
    local target="$1"
    cd "$solver_dir"
    if [[ "${FLUID_ST_CLEAN_BUILD:-0}" == 1 ]]; then
        make -f "$makefile" clean_st
    fi
    make -f "$makefile" "$target"
}

prepare_case() {
    # meshgen intentionally recreates result_dir. Do not call it implicitly
    # from offline/online because that could erase valid FOM snapshots.
    FLUID_ST_RESULT_DIR="$result_dir" \
      bash "${script_dir}/meshgen.sh" "$E" "$EP" "$NM" "$ND" "$NP" "$PA"

    cd "$result_dir"
    local mesh2graph="${gedatsu_bin}/gedatsu_simple_mesh2graph_convertor"
    local mesh_part="${gedatsu_bin}/gedatsu_simple_mesh_partitioner"
    local graph_part="${gedatsu_bin}/gedatsu_nodal_graph_partitioner"
    local bc_part="${gedatsu_bin}/gedatsu_bc_partitioner_R"
    for tool in "$mesh2graph" "$mesh_part" "$graph_part" "$bc_part"; do
        [[ -x "$tool" ]] || { echo "ERROR: missing GEDATSU tool: $tool" >&2; exit 1; }
    done
    "$mesh2graph" -i elem.dat -o graph.dat
    rm -rf parted.0
    "$mesh_part" -n "$NP"
    "$graph_part" -n "$NP" -i graph.dat
    "$bc_part" -n "$NP" -i D_bc_v.dat -ig node.dat
    mkdir -p parted.0
    printf '%s\n' "$NP" > parted.0/.fluid_st_spatial_np
    local rank f
    for ((rank=0; rank<NP; rank++)); do
        for f in node.dat elem.dat graph.dat D_bc_v.dat; do
            [[ -f "parted.0/${f}.${rank}" ]] || { echo "ERROR: missing parted.0/${f}.${rank}" >&2; exit 1; }
        done
    done
    echo "Prepared ST cavity case: $result_dir"
}

run_fom() {
    require_case_inputs
    build_target st_fom
    cp -f "${solver_dir}/fluid_sups_st_fom" "$result_dir/"
    local windows
    windows="$(window_dir)"
    mkdir -p "$windows"
    cd "$result_dir"
    env \
      OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}" OMP_NUM_THREADS=1 \
      FLUID_ST_INITIAL_CONDITION="$FLUID_ST_INITIAL_CONDITION" \
      FLUID_ST_TIME_DEGREE="$FLUID_ST_TIME_DEGREE" \
      FLUID_ST_SLABS_PER_WINDOW="$FLUID_ST_SLABS_PER_WINDOW" \
      FLUID_ST_SOLVE_MODE="$FLUID_ST_SOLVE_MODE" \
      FLUID_ST_JOLD_MODE="$FLUID_ST_JOLD_MODE" \
      FLUID_ST_WRITE_WINDOW_BINARY="$FLUID_ST_WRITE_WINDOW_BINARY" \
      FLUID_ST_WINDOW_DIR="$windows" \
      mpirun -np "$NP" ./fluid_sups_st_fom ./ \
      2>&1 | tee "fluid_st_fom_np${NP}_s${FLUID_ST_SLABS_PER_WINDOW}.log"
}

run_ddrom() {
    local mode="$1"
    require_case_inputs
    if [[ "$mode" == "offline" ]]; then
        require_fom_snapshots
    elif [[ "$mode" == "online" ]]; then
        require_basis
        if [[ "$FLUID_ST_DDROM_VALIDATE" == 1 ]]; then
            require_fom_snapshots
        fi
    fi

    build_target st_ddrom
    cp -f "${solver_dir}/fluid_sups_st_ddrom" "$result_dir/"
    local snapshots basis
    snapshots="$(ddrom_snapshot_dir)"
    basis="$(ddrom_basis_dir)"
    cd "$result_dir"
    env \
      OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}" OMP_NUM_THREADS=1 \
      FLUID_ST_INITIAL_CONDITION="$FLUID_ST_INITIAL_CONDITION" \
      FLUID_ST_TIME_DEGREE="$FLUID_ST_TIME_DEGREE" \
      FLUID_ST_SLABS_PER_WINDOW="$FLUID_ST_SLABS_PER_WINDOW" \
      FLUID_ST_SOLVE_MODE="$FLUID_ST_SOLVE_MODE" \
      FLUID_ST_JOLD_MODE="$FLUID_ST_JOLD_MODE" \
      FLUID_ST_DDROM_MODE="$mode" \
      FLUID_ST_DDROM_TRAINING_WINDOWS="$FLUID_ST_DDROM_TRAINING_WINDOWS" \
      FLUID_ST_DDROM_ONLINE_WINDOWS="$FLUID_ST_DDROM_ONLINE_WINDOWS" \
      FLUID_ST_DDROM_MAX_VELOCITY_MODES="$FLUID_ST_DDROM_MAX_VELOCITY_MODES" \
      FLUID_ST_DDROM_MAX_PRESSURE_MODES="$FLUID_ST_DDROM_MAX_PRESSURE_MODES" \
      FLUID_ST_DDROM_VELOCITY_POD_ENERGY="$FLUID_ST_DDROM_VELOCITY_POD_ENERGY" \
      FLUID_ST_DDROM_PRESSURE_POD_ENERGY="$FLUID_ST_DDROM_PRESSURE_POD_ENERGY" \
      FLUID_ST_DDROM_VALIDATE="$FLUID_ST_DDROM_VALIDATE" \
      FLUID_ST_DDROM_SNAPSHOT_DIR="$snapshots" \
      FLUID_ST_DDROM_BASIS_DIR="$basis" \
      FLUID_ST_DDROM_WINDOW_DIR="${result_dir}/fluid_st_ddrom_windows" \
      FLUID_ST_DDROM_COEFFICIENT_CSV="${result_dir}/fluid_st_ddrom_coefficients.csv" \
      FLUID_ST_DDROM_ACCURACY_CSV="${result_dir}/fluid_st_ddrom_accuracy.csv" \
      FLUID_ST_DDROM_WRITE_OUTPUT="$FLUID_ST_DDROM_WRITE_OUTPUT" \
      FLUID_ST_DDROM_WRITE_WINDOW_BINARY="$FLUID_ST_DDROM_WRITE_WINDOW_BINARY" \
      mpirun -np "$NP" ./fluid_sups_st_ddrom ./ \
      2>&1 | tee "fluid_st_ddrom_${mode}_np${NP}.log"
}

verify_case() {
    require_case_inputs
    require_fom_snapshots
    require_basis
    python3 "${script_dir}/verify_validated_case.py" "$result_dir"
}

case "$action" in
  prepare) prepare_case ;;
  fom) run_fom ;;
  offline) run_ddrom offline ;;
  online) run_ddrom online ;;
  verify) verify_case ;;
  validated)
    require_case_inputs
    require_fom_snapshots
    run_ddrom offline
    run_ddrom online
    verify_case
    ;;
  validated_full)
    prepare_case
    run_fom
    run_ddrom offline
    run_ddrom online
    verify_case
    ;;
  status) print_status ;;
  *) usage ;;
esac
