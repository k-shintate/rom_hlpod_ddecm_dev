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
directory="result_diff/${nm}-${np}-${nd}"

rm -r $directory
mkdir -p $directory
cd $directory

rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 1
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.01 "#output_interval" int 1 1  "#finish_time" double 1 1.0

./../../../test_thermal/bin/meshgen_hex $e $e $e $ep $ep $ep
./../../../test_thermal/bin/surf_dbc_all 1 1.0

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n $nd
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc.dat -ig node.dat

cd ../..