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

rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 10
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.000001 "#output_interval" int 1 10  "#finish_time" double 1 1.0

python3 ./../../shell/mag/node2bc.py ./node.dat > node_bc.dat
python3 ./../../shell/mag/node2dist_val.py ./node.dat node_distval.dat graph.dat

#python3 ./../../shell/mag/node2bc.py ./nedelec_node.dat > nedelec_node_bc.dat
#python3 ./../../shell/mag/node2dist_val.py ./nedelec_node.dat nedelec_node_distval.dat 773300

mkdir -p parted.0/parted.1

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $nd -i graph.dat -d ./parted.0/parted.1
cp -r parted.0/parted.1/metagraph.dat ./parted.0
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -ie elem.dat -o nodal_graph.dat

./../../utils/bin/elem2graph -ie nedelec_elem.dat -og graph_nedelec_elem.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $nd -i graph_nedelec_elem.dat -ig graph.dat -d ./parted.0/parted.1

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_I -n $nd -i nedelec_edge_sign.dat -ig graph_nedelec_elem.dat -d ./parted.0/parted.1

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_I -n $nd -i elem_bool.dat -ig graph_nedelec_elem.dat -d ./parted.0/parted.1
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n $nd -i node_distval.dat -ig graph.dat -d ./parted.0/parted.1
#./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n $nd -i nedelec_node_distval.dat -ig graph.dat 

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc.dat -ig graph.dat -d ./parted.0/parted.1
#./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i InnerSphere_quad_bc.dat -ig graph.dat


cd parted.0/parted.1
./../../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $np -i ../parted.0/metagraph.dat -d ./
cd ../..

cd ../..
