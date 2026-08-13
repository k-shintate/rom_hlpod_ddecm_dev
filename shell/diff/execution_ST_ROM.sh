#!/bin/bash

set -euo pipefail

if [ "$#" -lt 7 ]; then
    echo "Usage: $0 <arg1> <arg2> <nm> <nd> <np> <pa> <st>" >&2
    echo "  nm: number of shared POD modes" >&2
    echo "  nd: number of POD subdomains" >&2
    echo "  np: number used in the result-directory name" >&2
    echo "  pa: POD threshold exponent" >&2
    echo "  st: legacy solver type" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Online reduced solver parameters
# ---------------------------------------------------------------------------

rom_solver="${ST_ROM_SOLVER:-verify_monolis}"

monolis_epsilon="${ST_ROM_MONOLIS_EPSILON:-1.0e-12}"
monolis_max_iter="${ST_ROM_MONOLIS_MAX_ITER:-10000}"
monolis_verify_tolerance="${ST_ROM_MONOLIS_VERIFY_TOL:-1.0e-10}"

blas_threads="${OPENBLAS_NUM_THREADS:-4}"
rom_validate="${ST_ROM_VALIDATE:-1}"
rom_write_output="${ST_ROM_WRITE_OUTPUT:-0}"

# PODモード数
nm="$3"

# POD計算領域数
nd="$4"

# 結果ディレクトリ名に使用する並列領域数
np="$5"

# 基底本数可変の閾値 1.0E-{pa}
pa="$6"

# solver type
st="$7"

project_root="$(pwd)"
directory="${project_root}/result_diff/${nm}-${np}-${nd}"

# 外側の環境変数で上書き可能
blas_threads="${OPENBLAS_NUM_THREADS:-4}"
rom_validate="${ST_ROM_VALIDATE:-1}"
rom_write_output="${ST_ROM_WRITE_OUTPUT:-0}"

# 1ならofflineを再構築する
rebuild_offline="${ST_ROM_REBUILD_OFFLINE:-1}"

# 現在のserial ST-ROM実装では1プロセス
rom_mpi_np=1

source "${project_root}/shell/install.sh"

mkdir -p "${directory}"

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

cd "${project_root}/solvers/diff"

make -f Makefile_ST clean
make -f Makefile_ST

cp -f diff_ST_MPI "${directory}/"

cd "${directory}"

mkdir -p \
    pod_modes_vtk \
    pod_modes \
    fem_solver_prm \
    pod_solver_prm \
    calctime \
    DDECM

for ((i=0; i<nd; i++)); do
    mkdir -p "pod_modes/subdomain${i}"
done

# ---------------------------------------------------------------------------
# GDB command
# ---------------------------------------------------------------------------

fname="gdb_cmd"

cat > "${fname}" <<GDBEOF
run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}
backtrace
exit
GDBEOF

# ---------------------------------------------------------------------------
# Paths
#
# collect/offline/onlineの全段階で同じ絶対パスを使用する。
# ---------------------------------------------------------------------------

snapshot_dir="${directory}/strom_snapshots_mesh8000"
basis_dir="${directory}/strom_basis_shared_mesh8000_r${nm}"
validation_csv="${directory}/strom_validation_shared_r${nm}.csv"

snapshot0="${snapshot_dir}/snapshot_0000_window_0000.bin"

basis_file="${basis_dir}/all_windows_shared_pod.bin"

operator_file="${basis_dir}/all_windows_shared_reduced_operators.bin"

mkdir -p "${snapshot_dir}"

# ---------------------------------------------------------------------------
# 1. FOM snapshot collection
# ---------------------------------------------------------------------------

if [ ! -f "${snapshot0}" ]; then
    echo "Collecting FOM snapshots..."

    env \
        ST_ROM_MODE=collect \
        ST_ROM_SNAPSHOT_ID=0 \
        ST_ROM_SNAPSHOT_DIR="${snapshot_dir}" \
        mpirun -np "${rom_mpi_np}" ./diff_ST_MPI \
        2>&1 | tee collect.log
else
    echo "Using existing snapshots: ${snapshot_dir}"
fi

snapshot_count=$(
    find "${snapshot_dir}" \
        -maxdepth 1 \
        -type f \
        -name 'snapshot_0000_window_*.bin' |
    wc -l |
    tr -d ' '
)

if [ "${snapshot_count}" -ne 25 ]; then
    echo \
        "ERROR: expected 25 window snapshots, found ${snapshot_count}" \
        >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. All-window shared POD + reduced D/L construction
#
# basisとD/Lを同一のoffline実行で生成し、
# 基底と縮約演算子の不整合を防ぐ。
#
# direct couplingの検証が完了するまではlegacy経路を使用する。
# ---------------------------------------------------------------------------

if [ "${rebuild_offline}" -eq 1 ] ||
   [ ! -f "${basis_file}" ] ||
   [ ! -f "${operator_file}" ]; then

    echo \
        "Building all-window shared POD basis and reduced operators..."

    # 古い基底と縮約演算子の組合せを完全に削除
    rm -rf "${basis_dir}"
    mkdir -p "${basis_dir}"

    env \
        OPENBLAS_NUM_THREADS="${blas_threads}" \
        OMP_NUM_THREADS=1 \
        ST_ROM_MODE=offline \
        ST_ROM_NUM_SNAPSHOTS=1 \
        ST_ROM_EPSILON=1.0 \
        ST_ROM_MAX_MODES="${nm}" \
        ST_ROM_OPERATOR_REUSE=1 \
        ST_ROM_DIRECT_COUPLING=0 \
        ST_ROM_VERIFY_DIRECT_B=0 \
        ST_ROM_SNAPSHOT_DIR="${snapshot_dir}" \
        ST_ROM_BASIS_DIR="${basis_dir}" \
        mpirun -np "${rom_mpi_np}" ./diff_ST_MPI \
        2>&1 | tee "offline_shared_r${nm}.log"
else
    echo \
        "Using existing shared basis and reduced operators: ${basis_dir}"
fi

# ---------------------------------------------------------------------------
# Offline output checks
# ---------------------------------------------------------------------------

if [ ! -f "${basis_file}" ]; then
    echo \
        "ERROR: shared basis was not created: ${basis_file}" \
        >&2
    exit 1
fi

if [ ! -f "${operator_file}" ]; then
    echo \
        "ERROR: reduced operator file was not created: ${operator_file}" \
        >&2
    exit 1
fi

ls -lh \
    "${basis_file}" \
    "${operator_file}"

# ---------------------------------------------------------------------------
# 3. Online all-at-once solve
#
# offlineで生成したD/Lを読み込む。
# onlineでは演算子を再構築しない。
# ---------------------------------------------------------------------------

echo \
    "Running and validating online all-window shared ST-ROM..."

env \
    OPENBLAS_NUM_THREADS="${blas_threads}" \
    OMP_NUM_THREADS=1 \
    ST_ROM_MODE=online \
    ST_ROM_SOLVER="${rom_solver}" \
    ST_ROM_OPERATOR_REUSE=1 \
    ST_ROM_REBUILD_OPERATORS=0 \
    ST_ROM_DIRECT_COUPLING=0 \
    ST_ROM_VERIFY_DIRECT_B=0 \
    ST_ROM_BASIS_DIR="${basis_dir}" \
    ST_ROM_SNAPSHOT_DIR="${snapshot_dir}" \
    ST_ROM_VALIDATE="${rom_validate}" \
    ST_ROM_WRITE_OUTPUT="${rom_write_output}" \
    ST_ROM_VALIDATION_SNAPSHOT_ID=0 \
    ST_ROM_VALIDATION_CSV="${validation_csv}" \
    ST_ROM_MONOLIS_EPSILON="${monolis_epsilon}" \
    ST_ROM_MONOLIS_MAX_ITER="${monolis_max_iter}" \
    ST_ROM_MONOLIS_VERIFY_TOL="${monolis_verify_tolerance}" \
    mpirun -np "${rom_mpi_np}" ./diff_ST_MPI \
    2>&1 | tee "online_shared_r${nm}.log"

echo
echo "Shared basis:       ${basis_file}"
echo "Reduced operators:  ${operator_file}"
echo "Validation output:  ${validation_csv}"
echo "Online log:         ${directory}/online_shared_r${nm}.log"

cd "${project_root}"