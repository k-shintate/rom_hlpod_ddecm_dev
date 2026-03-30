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
./../../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 10
mv cond.dat rom_cond.dat

./../../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.000001 "#output_interval" int 1 1000  "#finish_time" double 1 2.0


python3 ./../../shell/mag/node2bc.py ./node.dat > node_bc.dat
python3 ./../../shell/mag/node2dist_val.py ./node.dat node_distval.dat 773833
#python3 ./../../shell/mag/node2dist_val.py ./node.dat node_distval.dat 1602391

#python3 ./../../shell/mag/node2bc.py ./nedelec_node.dat > nedelec_node_bc.dat
#python3 ./../../shell/mag/node2dist_val.py ./nedelec_node.dat nedelec_node_distval.dat 194562

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $np -i graph.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -ie elem.dat -o nodal_graph.dat

./../../utils/bin/elem2graph -ie elem.dat -og graph_elem.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i graph_elem.dat -ig graph.dat
./../../utils/bin/elem2graph -ie nedelec_elem_only.dat -og graph_nedelec_elem_only.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i graph_nedelec_elem_only.dat -ig graph.dat
./../../utils/bin/elem2graph -ie nedelec_elem.dat -og graph_nedelec_elem.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i graph_nedelec_elem.dat -ig graph.dat

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_I -n $np -i nedelec_edge_sign.dat -ig graph_nedelec_elem.dat 

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_I -n $np -i elem_bool.dat -ig graph_nedelec_elem.dat 
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n $np -i node_distval.dat -ig graph.dat 
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n $np -i nedelec_node_distval.dat -ig graph.dat 

for i in $(seq 0 $((np-1)))
do
    python3 ../../shell/mag/distval2node.py ./parted.0/node_distval.dat.$i  ./parted.0/node.dat.$i
    #python3 ../../shell/mag/distval2node.py ./parted.0/nedelec_node_distval.dat.$i  ./parted.0/nedelec_node.dat.$i
done

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $np -i D_bc.dat -ig graph.dat

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $np -i InnerSphere_quad_bc.dat -ig graph.dat


cd ../..
