#!/usr/bin/env bash

N_PARTS=128
CASE_DIR="result_mag/5-64-128"

# 結果ディレクトリへ移動
#cd "${CASE_DIR}"

GEDATSU_BIN="./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin"
UTIL_BIN="./../../utils/bin"
MAG_SHELL="../../shell/mag"

# Nedelec 要素グラフを変換
python3 "${MAG_SHELL}/elemtograph.py" \
  graph_nedelec_elem.dat \
  graph_nedelec_elem_test.dat

# nodal graph partition
"${GEDATSU_BIN}/gedatsu_nodal_graph_partitioner" \
  -n "${N_PARTS}" \
  -i graph_nedelec_elem_test.dat

# element connectivity graph 作成
"${UTIL_BIN}/elem2graph" \
  -ie part_elem_conn.dat \
  -og graph_elem_conn.dat

# connectivity graph partition
"${GEDATSU_BIN}/gedatsu_connectivity_graph_partitioner" \
  -n "${N_PARTS}" \
  -i graph_nedelec_elem.dat \
  -ig graph_nedelec_elem_test.dat

# connectivity graph partition
"${GEDATSU_BIN}/gedatsu_connectivity_graph_partitioner" \
  -n "${N_PARTS}" \
  -i graph_elem_conn.dat \
  -ig graph_nedelec_elem_test.dat

"${GEDATSU_BIN}/gedatsu_dist_val_partitioner_I" \
  -n "${N_PARTS}" \
  -i elem_bool.dat \
  -ig graph_nedelec_elem.dat

# elem_phi.dat.* 作成
for i in $(seq 0 $((N_PARTS - 1))); do
  python3 "${MAG_SHELL}/find_part_elem.py" \
    "parted.0/graph_nedelec_elem.dat.${i}" \
    "parted.0/graph_elem_conn.dat.${i}" \
    -o "parted.0/elem_phi.dat.${i}"
done

# node coordinate partition
"${GEDATSU_BIN}/gedatsu_dist_val_partitioner_R" \
  -n "${N_PARTS}" \
  -i node_coordinate_elem.dat \
  -ig graph_nedelec_elem.dat

# nedelec element partition
"${GEDATSU_BIN}/gedatsu_dist_val_partitioner_I" \
  -n "${N_PARTS}" \
  -i nedelec_elem.dat \
  -ig graph_nedelec_elem.dat

# nedelec edge sign partition
"${GEDATSU_BIN}/gedatsu_dist_val_partitioner_I" \
  -n "${N_PARTS}" \
  -i nedelec_edge_sign.dat \
  -ig graph_nedelec_elem.dat

# D boundary condition partition
"${GEDATSU_BIN}/gedatsu_bc_partitioner_R" \
  -n "${N_PARTS}" \
  -i D_bc_ned.dat \
  -ig graph_nedelec_elem_test.dat

echo "Done."


