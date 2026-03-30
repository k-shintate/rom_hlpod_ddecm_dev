#!/bin/bash

#mesh
#一方向分割数
e=$1
#解析領域の大きさ
ep=$2

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6

# 実行ディレクトリ
directory="result_mag/${nm}-${np}-${nd}"

#rm -r $directory
mkdir -p $directory
cd $directory

mv parted.0/D_bc.dat.0 parted.0/D_bc.dat.14_tmp
mv parted.0/D_bc.dat.1 parted.0/D_bc.dat.15_tmp
mv parted.0/D_bc.dat.2 parted.0/D_bc.dat.1_tmp
mv parted.0/D_bc.dat.3 parted.0/D_bc.dat.2_tmp
mv parted.0/D_bc.dat.4 parted.0/D_bc.dat.3_tmp
mv parted.0/D_bc.dat.5 parted.0/D_bc.dat.7_tmp
mv parted.0/D_bc.dat.6 parted.0/D_bc.dat.0_tmp
mv parted.0/D_bc.dat.7 parted.0/D_bc.dat.6_tmp
mv parted.0/D_bc.dat.8 parted.0/D_bc.dat.12_tmp
mv parted.0/D_bc.dat.9 parted.0/D_bc.dat.13_tmp
mv parted.0/D_bc.dat.10 parted.0/D_bc.dat.4_tmp
mv parted.0/D_bc.dat.11 parted.0/D_bc.dat.10_tmp
mv parted.0/D_bc.dat.12 parted.0/D_bc.dat.5_tmp
mv parted.0/D_bc.dat.13 parted.0/D_bc.dat.11_tmp
mv parted.0/D_bc.dat.14 parted.0/D_bc.dat.8_tmp
mv parted.0/D_bc.dat.15 parted.0/D_bc.dat.9_tmp

mv parted.0/D_bc.dat.0_tmp parted.0/D_bc.dat.0
mv parted.0/D_bc.dat.1_tmp parted.0/D_bc.dat.1
mv parted.0/D_bc.dat.2_tmp parted.0/D_bc.dat.2
mv parted.0/D_bc.dat.3_tmp parted.0/D_bc.dat.3
mv parted.0/D_bc.dat.4_tmp parted.0/D_bc.dat.4
mv parted.0/D_bc.dat.5_tmp parted.0/D_bc.dat.5
mv parted.0/D_bc.dat.6_tmp parted.0/D_bc.dat.6
mv parted.0/D_bc.dat.7_tmp parted.0/D_bc.dat.7
mv parted.0/D_bc.dat.8_tmp parted.0/D_bc.dat.8
mv parted.0/D_bc.dat.9_tmp parted.0/D_bc.dat.9
mv parted.0/D_bc.dat.10_tmp parted.0/D_bc.dat.10
mv parted.0/D_bc.dat.11_tmp parted.0/D_bc.dat.11
mv parted.0/D_bc.dat.12_tmp parted.0/D_bc.dat.12
mv parted.0/D_bc.dat.13_tmp parted.0/D_bc.dat.13
mv parted.0/D_bc.dat.14_tmp parted.0/D_bc.dat.14
mv parted.0/D_bc.dat.15_tmp parted.0/D_bc.dat.15


mv parted.0/graph.dat.n_internal.0 parted.0/graph.dat.n_internal.14_tmp
mv parted.0/graph.dat.n_internal.1 parted.0/graph.dat.n_internal.15_tmp
mv parted.0/graph.dat.n_internal.2 parted.0/graph.dat.n_internal.1_tmp
mv parted.0/graph.dat.n_internal.3 parted.0/graph.dat.n_internal.2_tmp
mv parted.0/graph.dat.n_internal.4 parted.0/graph.dat.n_internal.3_tmp
mv parted.0/graph.dat.n_internal.5 parted.0/graph.dat.n_internal.7_tmp
mv parted.0/graph.dat.n_internal.6 parted.0/graph.dat.n_internal.0_tmp
mv parted.0/graph.dat.n_internal.7 parted.0/graph.dat.n_internal.6_tmp
mv parted.0/graph.dat.n_internal.8 parted.0/graph.dat.n_internal.12_tmp
mv parted.0/graph.dat.n_internal.9 parted.0/graph.dat.n_internal.13_tmp
mv parted.0/graph.dat.n_internal.10 parted.0/graph.dat.n_internal.4_tmp
mv parted.0/graph.dat.n_internal.11 parted.0/graph.dat.n_internal.10_tmp
mv parted.0/graph.dat.n_internal.12 parted.0/graph.dat.n_internal.5_tmp
mv parted.0/graph.dat.n_internal.13 parted.0/graph.dat.n_internal.11_tmp
mv parted.0/graph.dat.n_internal.14 parted.0/graph.dat.n_internal.8_tmp
mv parted.0/graph.dat.n_internal.15 parted.0/graph.dat.n_internal.9_tmp

mv parted.0/graph.dat.n_internal.0_tmp parted.0/graph.dat.n_internal.0
mv parted.0/graph.dat.n_internal.1_tmp parted.0/graph.dat.n_internal.1
mv parted.0/graph.dat.n_internal.2_tmp parted.0/graph.dat.n_internal.2
mv parted.0/graph.dat.n_internal.3_tmp parted.0/graph.dat.n_internal.3
mv parted.0/graph.dat.n_internal.4_tmp parted.0/graph.dat.n_internal.4
mv parted.0/graph.dat.n_internal.5_tmp parted.0/graph.dat.n_internal.5
mv parted.0/graph.dat.n_internal.6_tmp parted.0/graph.dat.n_internal.6
mv parted.0/graph.dat.n_internal.7_tmp parted.0/graph.dat.n_internal.7
mv parted.0/graph.dat.n_internal.8_tmp parted.0/graph.dat.n_internal.8
mv parted.0/graph.dat.n_internal.9_tmp parted.0/graph.dat.n_internal.9
mv parted.0/graph.dat.n_internal.10_tmp parted.0/graph.dat.n_internal.10
mv parted.0/graph.dat.n_internal.11_tmp parted.0/graph.dat.n_internal.11
mv parted.0/graph.dat.n_internal.12_tmp parted.0/graph.dat.n_internal.12
mv parted.0/graph.dat.n_internal.13_tmp parted.0/graph.dat.n_internal.13
mv parted.0/graph.dat.n_internal.14_tmp parted.0/graph.dat.n_internal.14
mv parted.0/graph.dat.n_internal.15_tmp parted.0/graph.dat.n_internal.15


mv parted.0/elem_bool.dat.0 parted.0/elem_bool.dat.14_tmp
mv parted.0/elem_bool.dat.1 parted.0/elem_bool.dat.15_tmp
mv parted.0/elem_bool.dat.2 parted.0/elem_bool.dat.1_tmp
mv parted.0/elem_bool.dat.3 parted.0/elem_bool.dat.2_tmp
mv parted.0/elem_bool.dat.4 parted.0/elem_bool.dat.3_tmp
mv parted.0/elem_bool.dat.5 parted.0/elem_bool.dat.7_tmp
mv parted.0/elem_bool.dat.6 parted.0/elem_bool.dat.0_tmp
mv parted.0/elem_bool.dat.7 parted.0/elem_bool.dat.6_tmp
mv parted.0/elem_bool.dat.8 parted.0/elem_bool.dat.12_tmp
mv parted.0/elem_bool.dat.9 parted.0/elem_bool.dat.13_tmp
mv parted.0/elem_bool.dat.10 parted.0/elem_bool.dat.4_tmp
mv parted.0/elem_bool.dat.11 parted.0/elem_bool.dat.10_tmp
mv parted.0/elem_bool.dat.12 parted.0/elem_bool.dat.5_tmp
mv parted.0/elem_bool.dat.13 parted.0/elem_bool.dat.11_tmp
mv parted.0/elem_bool.dat.14 parted.0/elem_bool.dat.8_tmp
mv parted.0/elem_bool.dat.15 parted.0/elem_bool.dat.9_tmp

mv parted.0/elem_bool.dat.0_tmp parted.0/elem_bool.dat.0
mv parted.0/elem_bool.dat.1_tmp parted.0/elem_bool.dat.1
mv parted.0/elem_bool.dat.2_tmp parted.0/elem_bool.dat.2
mv parted.0/elem_bool.dat.3_tmp parted.0/elem_bool.dat.3
mv parted.0/elem_bool.dat.4_tmp parted.0/elem_bool.dat.4
mv parted.0/elem_bool.dat.5_tmp parted.0/elem_bool.dat.5
mv parted.0/elem_bool.dat.6_tmp parted.0/elem_bool.dat.6
mv parted.0/elem_bool.dat.7_tmp parted.0/elem_bool.dat.7
mv parted.0/elem_bool.dat.8_tmp parted.0/elem_bool.dat.8
mv parted.0/elem_bool.dat.9_tmp parted.0/elem_bool.dat.9
mv parted.0/elem_bool.dat.10_tmp parted.0/elem_bool.dat.10
mv parted.0/elem_bool.dat.11_tmp parted.0/elem_bool.dat.11
mv parted.0/elem_bool.dat.12_tmp parted.0/elem_bool.dat.12
mv parted.0/elem_bool.dat.13_tmp parted.0/elem_bool.dat.13
mv parted.0/elem_bool.dat.14_tmp parted.0/elem_bool.dat.14
mv parted.0/elem_bool.dat.15_tmp parted.0/elem_bool.dat.15


mv parted.0/graph.dat.0 parted.0/graph.dat.14_tmp
mv parted.0/graph.dat.1 parted.0/graph.dat.15_tmp
mv parted.0/graph.dat.2 parted.0/graph.dat.1_tmp
mv parted.0/graph.dat.3 parted.0/graph.dat.2_tmp
mv parted.0/graph.dat.4 parted.0/graph.dat.3_tmp
mv parted.0/graph.dat.5 parted.0/graph.dat.7_tmp
mv parted.0/graph.dat.6 parted.0/graph.dat.0_tmp
mv parted.0/graph.dat.7 parted.0/graph.dat.6_tmp
mv parted.0/graph.dat.8 parted.0/graph.dat.12_tmp
mv parted.0/graph.dat.9 parted.0/graph.dat.13_tmp
mv parted.0/graph.dat.10 parted.0/graph.dat.4_tmp
mv parted.0/graph.dat.11 parted.0/graph.dat.10_tmp
mv parted.0/graph.dat.12 parted.0/graph.dat.5_tmp
mv parted.0/graph.dat.13 parted.0/graph.dat.11_tmp
mv parted.0/graph.dat.14 parted.0/graph.dat.8_tmp
mv parted.0/graph.dat.15 parted.0/graph.dat.9_tmp

mv parted.0/graph.dat.0_tmp parted.0/graph.dat.0
mv parted.0/graph.dat.1_tmp parted.0/graph.dat.1
mv parted.0/graph.dat.2_tmp parted.0/graph.dat.2
mv parted.0/graph.dat.3_tmp parted.0/graph.dat.3
mv parted.0/graph.dat.4_tmp parted.0/graph.dat.4
mv parted.0/graph.dat.5_tmp parted.0/graph.dat.5
mv parted.0/graph.dat.6_tmp parted.0/graph.dat.6
mv parted.0/graph.dat.7_tmp parted.0/graph.dat.7
mv parted.0/graph.dat.8_tmp parted.0/graph.dat.8
mv parted.0/graph.dat.9_tmp parted.0/graph.dat.9
mv parted.0/graph.dat.10_tmp parted.0/graph.dat.10
mv parted.0/graph.dat.11_tmp parted.0/graph.dat.11
mv parted.0/graph.dat.12_tmp parted.0/graph.dat.12
mv parted.0/graph.dat.13_tmp parted.0/graph.dat.13
mv parted.0/graph.dat.14_tmp parted.0/graph.dat.14
mv parted.0/graph.dat.15_tmp parted.0/graph.dat.15


mv parted.0/graph.dat.recv.0 parted.0/graph.dat.recv.14_tmp
mv parted.0/graph.dat.recv.1 parted.0/graph.dat.recv.15_tmp
mv parted.0/graph.dat.recv.2 parted.0/graph.dat.recv.1_tmp
mv parted.0/graph.dat.recv.3 parted.0/graph.dat.recv.2_tmp
mv parted.0/graph.dat.recv.4 parted.0/graph.dat.recv.3_tmp
mv parted.0/graph.dat.recv.5 parted.0/graph.dat.recv.7_tmp
mv parted.0/graph.dat.recv.6 parted.0/graph.dat.recv.0_tmp
mv parted.0/graph.dat.recv.7 parted.0/graph.dat.recv.6_tmp
mv parted.0/graph.dat.recv.8 parted.0/graph.dat.recv.12_tmp
mv parted.0/graph.dat.recv.9 parted.0/graph.dat.recv.13_tmp
mv parted.0/graph.dat.recv.10 parted.0/graph.dat.recv.4_tmp
mv parted.0/graph.dat.recv.11 parted.0/graph.dat.recv.10_tmp
mv parted.0/graph.dat.recv.12 parted.0/graph.dat.recv.5_tmp
mv parted.0/graph.dat.recv.13 parted.0/graph.dat.recv.11_tmp
mv parted.0/graph.dat.recv.14 parted.0/graph.dat.recv.8_tmp
mv parted.0/graph.dat.recv.15 parted.0/graph.dat.recv.9_tmp

mv parted.0/graph.dat.recv.0_tmp parted.0/graph.dat.recv.0
mv parted.0/graph.dat.recv.1_tmp parted.0/graph.dat.recv.1
mv parted.0/graph.dat.recv.2_tmp parted.0/graph.dat.recv.2
mv parted.0/graph.dat.recv.3_tmp parted.0/graph.dat.recv.3
mv parted.0/graph.dat.recv.4_tmp parted.0/graph.dat.recv.4
mv parted.0/graph.dat.recv.5_tmp parted.0/graph.dat.recv.5
mv parted.0/graph.dat.recv.6_tmp parted.0/graph.dat.recv.6
mv parted.0/graph.dat.recv.7_tmp parted.0/graph.dat.recv.7
mv parted.0/graph.dat.recv.8_tmp parted.0/graph.dat.recv.8
mv parted.0/graph.dat.recv.9_tmp parted.0/graph.dat.recv.9
mv parted.0/graph.dat.recv.10_tmp parted.0/graph.dat.recv.10
mv parted.0/graph.dat.recv.11_tmp parted.0/graph.dat.recv.11
mv parted.0/graph.dat.recv.12_tmp parted.0/graph.dat.recv.12
mv parted.0/graph.dat.recv.13_tmp parted.0/graph.dat.recv.13
mv parted.0/graph.dat.recv.14_tmp parted.0/graph.dat.recv.14
mv parted.0/graph.dat.recv.15_tmp parted.0/graph.dat.recv.15

mv parted.0/graph.dat.send.0 parted.0/graph.dat.send.14_tmp
mv parted.0/graph.dat.send.1 parted.0/graph.dat.send.15_tmp
mv parted.0/graph.dat.send.2 parted.0/graph.dat.send.1_tmp
mv parted.0/graph.dat.send.3 parted.0/graph.dat.send.2_tmp
mv parted.0/graph.dat.send.4 parted.0/graph.dat.send.3_tmp
mv parted.0/graph.dat.send.5 parted.0/graph.dat.send.7_tmp
mv parted.0/graph.dat.send.6 parted.0/graph.dat.send.0_tmp
mv parted.0/graph.dat.send.7 parted.0/graph.dat.send.6_tmp
mv parted.0/graph.dat.send.8 parted.0/graph.dat.send.12_tmp
mv parted.0/graph.dat.send.9 parted.0/graph.dat.send.13_tmp
mv parted.0/graph.dat.send.10 parted.0/graph.dat.send.4_tmp
mv parted.0/graph.dat.send.11 parted.0/graph.dat.send.10_tmp
mv parted.0/graph.dat.send.12 parted.0/graph.dat.send.5_tmp
mv parted.0/graph.dat.send.13 parted.0/graph.dat.send.11_tmp
mv parted.0/graph.dat.send.14 parted.0/graph.dat.send.8_tmp
mv parted.0/graph.dat.send.15 parted.0/graph.dat.send.9_tmp

mv parted.0/graph.dat.send.0_tmp parted.0/graph.dat.send.0
mv parted.0/graph.dat.send.1_tmp parted.0/graph.dat.send.1
mv parted.0/graph.dat.send.2_tmp parted.0/graph.dat.send.2
mv parted.0/graph.dat.send.3_tmp parted.0/graph.dat.send.3
mv parted.0/graph.dat.send.4_tmp parted.0/graph.dat.send.4
mv parted.0/graph.dat.send.5_tmp parted.0/graph.dat.send.5
mv parted.0/graph.dat.send.6_tmp parted.0/graph.dat.send.6
mv parted.0/graph.dat.send.7_tmp parted.0/graph.dat.send.7
mv parted.0/graph.dat.send.8_tmp parted.0/graph.dat.send.8
mv parted.0/graph.dat.send.9_tmp parted.0/graph.dat.send.9
mv parted.0/graph.dat.send.10_tmp parted.0/graph.dat.send.10
mv parted.0/graph.dat.send.11_tmp parted.0/graph.dat.send.11
mv parted.0/graph.dat.send.12_tmp parted.0/graph.dat.send.12
mv parted.0/graph.dat.send.13_tmp parted.0/graph.dat.send.13
mv parted.0/graph.dat.send.14_tmp parted.0/graph.dat.send.14
mv parted.0/graph.dat.send.15_tmp parted.0/graph.dat.send.15

mv parted.0/graph_nedelec_elem.dat.0 parted.0/graph_nedelec_elem.dat.14_tmp
mv parted.0/graph_nedelec_elem.dat.1 parted.0/graph_nedelec_elem.dat.15_tmp
mv parted.0/graph_nedelec_elem.dat.2 parted.0/graph_nedelec_elem.dat.1_tmp
mv parted.0/graph_nedelec_elem.dat.3 parted.0/graph_nedelec_elem.dat.2_tmp
mv parted.0/graph_nedelec_elem.dat.4 parted.0/graph_nedelec_elem.dat.3_tmp
mv parted.0/graph_nedelec_elem.dat.5 parted.0/graph_nedelec_elem.dat.7_tmp
mv parted.0/graph_nedelec_elem.dat.6 parted.0/graph_nedelec_elem.dat.0_tmp
mv parted.0/graph_nedelec_elem.dat.7 parted.0/graph_nedelec_elem.dat.6_tmp
mv parted.0/graph_nedelec_elem.dat.8 parted.0/graph_nedelec_elem.dat.12_tmp
mv parted.0/graph_nedelec_elem.dat.9 parted.0/graph_nedelec_elem.dat.13_tmp
mv parted.0/graph_nedelec_elem.dat.10 parted.0/graph_nedelec_elem.dat.4_tmp
mv parted.0/graph_nedelec_elem.dat.11 parted.0/graph_nedelec_elem.dat.10_tmp
mv parted.0/graph_nedelec_elem.dat.12 parted.0/graph_nedelec_elem.dat.5_tmp
mv parted.0/graph_nedelec_elem.dat.13 parted.0/graph_nedelec_elem.dat.11_tmp
mv parted.0/graph_nedelec_elem.dat.14 parted.0/graph_nedelec_elem.dat.8_tmp
mv parted.0/graph_nedelec_elem.dat.15 parted.0/graph_nedelec_elem.dat.9_tmp

mv parted.0/graph_nedelec_elem.dat.0_tmp parted.0/graph_nedelec_elem.dat.0
mv parted.0/graph_nedelec_elem.dat.1_tmp parted.0/graph_nedelec_elem.dat.1
mv parted.0/graph_nedelec_elem.dat.2_tmp parted.0/graph_nedelec_elem.dat.2
mv parted.0/graph_nedelec_elem.dat.3_tmp parted.0/graph_nedelec_elem.dat.3
mv parted.0/graph_nedelec_elem.dat.4_tmp parted.0/graph_nedelec_elem.dat.4
mv parted.0/graph_nedelec_elem.dat.5_tmp parted.0/graph_nedelec_elem.dat.5
mv parted.0/graph_nedelec_elem.dat.6_tmp parted.0/graph_nedelec_elem.dat.6
mv parted.0/graph_nedelec_elem.dat.7_tmp parted.0/graph_nedelec_elem.dat.7
mv parted.0/graph_nedelec_elem.dat.8_tmp parted.0/graph_nedelec_elem.dat.8
mv parted.0/graph_nedelec_elem.dat.9_tmp parted.0/graph_nedelec_elem.dat.9
mv parted.0/graph_nedelec_elem.dat.10_tmp parted.0/graph_nedelec_elem.dat.10
mv parted.0/graph_nedelec_elem.dat.11_tmp parted.0/graph_nedelec_elem.dat.11
mv parted.0/graph_nedelec_elem.dat.12_tmp parted.0/graph_nedelec_elem.dat.12
mv parted.0/graph_nedelec_elem.dat.13_tmp parted.0/graph_nedelec_elem.dat.13
mv parted.0/graph_nedelec_elem.dat.14_tmp parted.0/graph_nedelec_elem.dat.14
mv parted.0/graph_nedelec_elem.dat.15_tmp parted.0/graph_nedelec_elem.dat.15

mv parted.0/nedelec_edge_sign.dat.0 parted.0/nedelec_edge_sign.dat.14_tmp
mv parted.0/nedelec_edge_sign.dat.1 parted.0/nedelec_edge_sign.dat.15_tmp
mv parted.0/nedelec_edge_sign.dat.2 parted.0/nedelec_edge_sign.dat.1_tmp
mv parted.0/nedelec_edge_sign.dat.3 parted.0/nedelec_edge_sign.dat.2_tmp
mv parted.0/nedelec_edge_sign.dat.4 parted.0/nedelec_edge_sign.dat.3_tmp
mv parted.0/nedelec_edge_sign.dat.5 parted.0/nedelec_edge_sign.dat.7_tmp
mv parted.0/nedelec_edge_sign.dat.6 parted.0/nedelec_edge_sign.dat.0_tmp
mv parted.0/nedelec_edge_sign.dat.7 parted.0/nedelec_edge_sign.dat.6_tmp
mv parted.0/nedelec_edge_sign.dat.8 parted.0/nedelec_edge_sign.dat.12_tmp
mv parted.0/nedelec_edge_sign.dat.9 parted.0/nedelec_edge_sign.dat.13_tmp
mv parted.0/nedelec_edge_sign.dat.10 parted.0/nedelec_edge_sign.dat.4_tmp
mv parted.0/nedelec_edge_sign.dat.11 parted.0/nedelec_edge_sign.dat.10_tmp
mv parted.0/nedelec_edge_sign.dat.12 parted.0/nedelec_edge_sign.dat.5_tmp
mv parted.0/nedelec_edge_sign.dat.13 parted.0/nedelec_edge_sign.dat.11_tmp
mv parted.0/nedelec_edge_sign.dat.14 parted.0/nedelec_edge_sign.dat.8_tmp
mv parted.0/nedelec_edge_sign.dat.15 parted.0/nedelec_edge_sign.dat.9_tmp

mv parted.0/nedelec_edge_sign.dat.0_tmp parted.0/nedelec_edge_sign.dat.0
mv parted.0/nedelec_edge_sign.dat.1_tmp parted.0/nedelec_edge_sign.dat.1
mv parted.0/nedelec_edge_sign.dat.2_tmp parted.0/nedelec_edge_sign.dat.2
mv parted.0/nedelec_edge_sign.dat.3_tmp parted.0/nedelec_edge_sign.dat.3
mv parted.0/nedelec_edge_sign.dat.4_tmp parted.0/nedelec_edge_sign.dat.4
mv parted.0/nedelec_edge_sign.dat.5_tmp parted.0/nedelec_edge_sign.dat.5
mv parted.0/nedelec_edge_sign.dat.6_tmp parted.0/nedelec_edge_sign.dat.6
mv parted.0/nedelec_edge_sign.dat.7_tmp parted.0/nedelec_edge_sign.dat.7
mv parted.0/nedelec_edge_sign.dat.8_tmp parted.0/nedelec_edge_sign.dat.8
mv parted.0/nedelec_edge_sign.dat.9_tmp parted.0/nedelec_edge_sign.dat.9
mv parted.0/nedelec_edge_sign.dat.10_tmp parted.0/nedelec_edge_sign.dat.10
mv parted.0/nedelec_edge_sign.dat.11_tmp parted.0/nedelec_edge_sign.dat.11
mv parted.0/nedelec_edge_sign.dat.12_tmp parted.0/nedelec_edge_sign.dat.12
mv parted.0/nedelec_edge_sign.dat.13_tmp parted.0/nedelec_edge_sign.dat.13
mv parted.0/nedelec_edge_sign.dat.14_tmp parted.0/nedelec_edge_sign.dat.14
mv parted.0/nedelec_edge_sign.dat.15_tmp parted.0/nedelec_edge_sign.dat.15

mv parted.0/node.dat.0 parted.0/node.dat.14_tmp
mv parted.0/node.dat.1 parted.0/node.dat.15_tmp
mv parted.0/node.dat.2 parted.0/node.dat.1_tmp
mv parted.0/node.dat.3 parted.0/node.dat.2_tmp
mv parted.0/node.dat.4 parted.0/node.dat.3_tmp
mv parted.0/node.dat.5 parted.0/node.dat.7_tmp
mv parted.0/node.dat.6 parted.0/node.dat.0_tmp
mv parted.0/node.dat.7 parted.0/node.dat.6_tmp
mv parted.0/node.dat.8 parted.0/node.dat.12_tmp
mv parted.0/node.dat.9 parted.0/node.dat.13_tmp
mv parted.0/node.dat.10 parted.0/node.dat.4_tmp
mv parted.0/node.dat.11 parted.0/node.dat.10_tmp
mv parted.0/node.dat.12 parted.0/node.dat.5_tmp
mv parted.0/node.dat.13 parted.0/node.dat.11_tmp
mv parted.0/node.dat.14 parted.0/node.dat.8_tmp
mv parted.0/node.dat.15 parted.0/node.dat.9_tmp

mv parted.0/node.dat.0_tmp parted.0/node.dat.0
mv parted.0/node.dat.1_tmp parted.0/node.dat.1
mv parted.0/node.dat.2_tmp parted.0/node.dat.2
mv parted.0/node.dat.3_tmp parted.0/node.dat.3
mv parted.0/node.dat.4_tmp parted.0/node.dat.4
mv parted.0/node.dat.5_tmp parted.0/node.dat.5
mv parted.0/node.dat.6_tmp parted.0/node.dat.6
mv parted.0/node.dat.7_tmp parted.0/node.dat.7
mv parted.0/node.dat.8_tmp parted.0/node.dat.8
mv parted.0/node.dat.9_tmp parted.0/node.dat.9
mv parted.0/node.dat.10_tmp parted.0/node.dat.10
mv parted.0/node.dat.11_tmp parted.0/node.dat.11
mv parted.0/node.dat.12_tmp parted.0/node.dat.12
mv parted.0/node.dat.13_tmp parted.0/node.dat.13
mv parted.0/node.dat.14_tmp parted.0/node.dat.14
mv parted.0/node.dat.15_tmp parted.0/node.dat.15

mv parted.0/graph.dat.id.0 parted.0/graph.dat.id.14_tmp
mv parted.0/graph.dat.id.1 parted.0/graph.dat.id.15_tmp
mv parted.0/graph.dat.id.2 parted.0/graph.dat.id.1_tmp
mv parted.0/graph.dat.id.3 parted.0/graph.dat.id.2_tmp
mv parted.0/graph.dat.id.4 parted.0/graph.dat.id.3_tmp
mv parted.0/graph.dat.id.5 parted.0/graph.dat.id.7_tmp
mv parted.0/graph.dat.id.6 parted.0/graph.dat.id.0_tmp
mv parted.0/graph.dat.id.7 parted.0/graph.dat.id.6_tmp
mv parted.0/graph.dat.id.8 parted.0/graph.dat.id.12_tmp
mv parted.0/graph.dat.id.9 parted.0/graph.dat.id.13_tmp
mv parted.0/graph.dat.id.10 parted.0/graph.dat.id.4_tmp
mv parted.0/graph.dat.id.11 parted.0/graph.dat.id.10_tmp
mv parted.0/graph.dat.id.12 parted.0/graph.dat.id.5_tmp
mv parted.0/graph.dat.id.13 parted.0/graph.dat.id.11_tmp
mv parted.0/graph.dat.id.14 parted.0/graph.dat.id.8_tmp
mv parted.0/graph.dat.id.15 parted.0/graph.dat.id.9_tmp

mv parted.0/graph.dat.id.0_tmp parted.0/graph.dat.id.0
mv parted.0/graph.dat.id.1_tmp parted.0/graph.dat.id.1
mv parted.0/graph.dat.id.2_tmp parted.0/graph.dat.id.2
mv parted.0/graph.dat.id.3_tmp parted.0/graph.dat.id.3
mv parted.0/graph.dat.id.4_tmp parted.0/graph.dat.id.4
mv parted.0/graph.dat.id.5_tmp parted.0/graph.dat.id.5
mv parted.0/graph.dat.id.6_tmp parted.0/graph.dat.id.6
mv parted.0/graph.dat.id.7_tmp parted.0/graph.dat.id.7
mv parted.0/graph.dat.id.8_tmp parted.0/graph.dat.id.8
mv parted.0/graph.dat.id.9_tmp parted.0/graph.dat.id.9
mv parted.0/graph.dat.id.10_tmp parted.0/graph.dat.id.10
mv parted.0/graph.dat.id.11_tmp parted.0/graph.dat.id.11
mv parted.0/graph.dat.id.12_tmp parted.0/graph.dat.id.12
mv parted.0/graph.dat.id.13_tmp parted.0/graph.dat.id.13
mv parted.0/graph.dat.id.14_tmp parted.0/graph.dat.id.14
mv parted.0/graph.dat.id.15_tmp parted.0/graph.dat.id.15



cd ../..
