#!/usr/bin/env bash
set -euo pipefail

# ST-FOM v10 temporal accuracy sweep
# Fixed spatial mesh; compare dG(0), dG(1), dG(2)
# at dt = 0.1, 0.05, 0.01 using final-time accuracy.
#
# Usage:
#   bash shell/execution/exec_compare_ST_FOM_accuracy_v10_dt3.sh \
#        <e> <ep> <np> [finish_time] [slabs_per_window]
#
# Recommended for tf=1.0 and dt={0.1,0.05,0.01}:
#   slabs_per_window=2
#
# Example:
#   bash shell/execution/exec_compare_ST_FOM_accuracy_v10_dt3.sh \
#        40 1.0 8 1.0 2

if (( $# < 3 || $# > 5 )); then
    echo "Usage: $0 <e> <ep> <np> [finish_time] [slabs_per_window]" >&2
    exit 2
fi

E="$1"
EP="$2"
NP="$3"
TF="${4:-1.0}"
SW="${5:-2}"

DEGREES=(0 1 2)
DTS=(0.1 0.05 0.01)

script_path="${BASH_SOURCE[0]:-$0}"
script_dir="$(cd "$(dirname "$script_path")" && pwd)"
ROOT_DIR="$(cd "$script_dir/../.." && pwd)"
RUNNER="$ROOT_DIR/shell/execution/exec_diff_ST_FOM_MPI.sh"
PLOTTER="$ROOT_DIR/shell/plot_st_accuracy_compare_v10_dt3.py"
SOLVER_SRC="$ROOT_DIR/solvers/diff/diff_core/core_FOM_ST_spaceblock.c"

if [[ ! -f "$RUNNER" ]]; then
    echo "ERROR: ST-FOM runner not found: $RUNNER" >&2
    exit 1
fi
if [[ ! -f "$PLOTTER" ]]; then
    echo "ERROR: plotter not found: $PLOTTER" >&2
    exit 1
fi
if [[ ! -f "$SOLVER_SRC" ]]; then
    echo "ERROR: solver source not found: $SOLVER_SRC" >&2
    exit 1
fi

# v10: matrix values must be cleared before every reassembly.
if ! grep -q "monolis_clear_mat_value_R" "$SOLVER_SRC"; then
    echo "ERROR: v10 matrix-clear fix does not appear to be installed." >&2
    echo "       Expected monolis_clear_mat_value_R in:" >&2
    echo "       $SOLVER_SRC" >&2
    exit 1
fi

calc_steps() {
    python3 - "$1" "$2" <<'PY'
import sys
tf = float(sys.argv[1])
dt = float(sys.argv[2])
n = tf / dt
nr = int(round(n))
if abs(n - nr) > 1.0e-10:
    raise SystemExit(2)
print(nr)
PY
}

# Check all dt values before doing expensive calculations.
for dt in "${DTS[@]}"; do
    steps="$(calc_steps "$TF" "$dt")" || {
        echo "ERROR: finish_time/dt is not an integer: tf=$TF dt=$dt" >&2
        exit 1
    }
    if (( steps % SW != 0 )); then
        echo "ERROR: total steps $steps is not divisible by slabs/window $SW for dt=$dt." >&2
        echo "       For tf=1.0 and dt={0.1,0.05,0.01}, use slabs/window=2." >&2
        exit 1
    fi
done

compare_tag="stfom-compare-e${E}-np${NP}-sw${SW}-tf${TF}-dt3"
COMPARE_DIR="$ROOT_DIR/result_diff/$compare_tag"
mkdir -p "$COMPARE_DIR"

SUMMARY_CSV="$COMPARE_DIR/final_accuracy_summary.csv"
cat > "$SUMMARY_CSV" <<'CSV'
mesh_e,length,num_mpi,finish_time,slabs_per_window,degree,dt,final_step,final_time,relative_l2,absolute_l2,linf,case_directory,accuracy_csv
CSV

echo "============================================================"
echo " ST-FOM v10 final-time temporal accuracy comparison"
echo "------------------------------------------------------------"
echo " spatial mesh       : ${E} x ${E} x ${E}"
echo " domain length      : ${EP}"
echo " MPI ranks          : ${NP}"
echo " finish time        : ${TF}"
echo " slabs/window       : ${SW}"
echo " time degrees       : dG(0), dG(1), dG(2)"
echo " dt                 : 0.1, 0.05, 0.01"
echo " solve mode         : causal"
echo " graph source       : partition"
echo " VTK during sweep   : disabled"
echo " result directory   : ${COMPARE_DIR}"
echo "============================================================"

for dg in "${DEGREES[@]}"; do
    for dt in "${DTS[@]}"; do
        steps="$(calc_steps "$TF" "$dt")"

        echo
        echo "============================================================"
        echo " Running dG(${dg}), dt=${dt}, steps=${steps}"
        echo "============================================================"

        ST_FOM_SOLVE_MODE=causal \
        ST_FOM_GRAPH_SOURCE=partition \
        ST_FOM_WRITE_VTK=0 \
        bash "$RUNNER" "$E" "$EP" "$NP" "$dg" "$SW" "$dt" "$TF"

        CASE_DIR="$ROOT_DIR/result_diff/stfom-e${E}-np${NP}-dg${dg}-sw${SW}-dt${dt}-tf${TF}"
        ACC_CSV="$CASE_DIR/st_accuracy.csv"

        if [[ ! -s "$ACC_CSV" ]]; then
            echo "ERROR: accuracy CSV missing or empty: $ACC_CSV" >&2
            exit 1
        fi

        FINAL_ROW="$(tail -n 1 "$ACC_CSV")"
        IFS=',' read -r final_step final_time csv_dt csv_degree csv_sw rel_l2 abs_l2 linf <<< "$FINAL_ROW"

        python3 - "$final_time" "$TF" <<'PY'
import sys
actual = float(sys.argv[1])
expected = float(sys.argv[2])
if abs(actual - expected) > 1.0e-10:
    raise SystemExit(f"final time mismatch: actual={actual}, expected={expected}")
PY

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$E" "$EP" "$NP" "$TF" "$SW" \
            "$csv_degree" "$csv_dt" "$final_step" "$final_time" \
            "$rel_l2" "$abs_l2" "$linf" "$CASE_DIR" "$ACC_CSV" \
            >> "$SUMMARY_CSV"

        echo "Final accuracy: dG(${csv_degree}) dt=${csv_dt} rel_L2=${rel_l2} abs_L2=${abs_l2} Linf=${linf}"
    done
done

python3 "$PLOTTER" "$SUMMARY_CSV" "$COMPARE_DIR"

echo
echo "============================================================"
echo " Comparison completed"
echo " summary CSV       : $SUMMARY_CSV"
echo " pairwise orders   : $COMPARE_DIR/observed_time_order_pairwise.csv"
echo " relative L2 plot  : $COMPARE_DIR/final_relative_l2_by_degree.png"
echo " dt convergence    : $COMPARE_DIR/final_relative_l2_vs_dt.png"
echo " Linf convergence  : $COMPARE_DIR/final_linf_vs_dt.png"
echo " table             : $COMPARE_DIR/final_accuracy_table.md"
echo "============================================================"
