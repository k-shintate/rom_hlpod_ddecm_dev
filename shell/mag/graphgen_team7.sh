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
./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 1
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.0001 "#output_interval" int 1 100  "#finish_time" double 1 1.0

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i elem.dat -o graph_elem.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i nedelec_elem.dat -o graph_nedelec_elem.dat

#python3 ./../../shell/mag/bool_elem.py ./elem_widing1.dat ./elem_widing2.dat ./elem_widing3.dat ./elem_iron.dat ./elem_air.dat  ./elem.dat ./elem_bool.dat
python3 ./../../shell/mag/bool_elem.py ./elem_widing1.dat ./elem_widing2.dat ./elem_widing3.dat ./elem_iron.dat ./elem_widing4.dat ./elem_air.dat  ./elem.dat ./elem_bool.dat

cd ../..
