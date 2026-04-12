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

./../../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.0001 "#output_interval" int 1 1  "#finish_time" double 1 2.0

###
python3 ../../solvers/mag/threephase_transformer/mesh_io/add_circuit_dof.py graph_nedelec_elem.dat graph_nedelec_elem_circuit.dat graph.dat 
python3 ../../shell/mag/merge_graph.py ./ --elem graph_elem.dat --nedelec graph_nedelec_elem_circuit.dat --out graph.dat

python3 ../../solvers/mag/threephase_transformer/mesh_io/add_circuit_dof_graph.py graph.dat 


cd ../..
