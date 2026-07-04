#!/usr/bin/env bash

N_PARTS=8
CASE_DIR="result_mag/5-8-128"

GEDATSU_BIN="./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin"
UTIL_BIN="./../../utils/bin"
MAG_SHELL="../../shell/mag"

# elem_phi.dat.* 作成
for i in $(seq 0 $((N_PARTS - 1))); do
  python3 "${MAG_SHELL}/find_part_elem.py" \
    "parted.0/graph_nedelec_elem.dat.${i}" \
    "parted.0/graph_elem_conn.dat.${i}" \
    -o "parted.0/elem_phi.dat.${i}"
done

echo "Done."

