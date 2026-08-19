#!/bin/bash

#mesh
#一方向分割数
e=$1
#解析領域の大きさ
ep=$2
np=$3

# 実行ディレクトリ
directory="result_diff/${nm}-${np}-${nd}"

cd solvers/fluid_sups/karman_vortex
. shell/meshgen_thermal_fin.sh
cd ../../..

rm -r $directory
mkdir -p $directory
cd $directory

mv ../../solvers/fluid_sups/karman_vortex/merged_mesh/parted.0 ./
mv ../../solvers/fluid_sups/karman_vortex/merged_mesh/D_bc.dat ./
mv ../../solvers/fluid_sups/karman_vortex/merged_mesh/elem.dat ./


# for parameteric study
rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 1
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.01 "#output_interval" int 1 1  "#finish_time" double 1 0.2

cd ../..