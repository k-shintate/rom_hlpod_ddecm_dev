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
directory="result_fluid_sups_karman_vortex/${nm}-${np}-${nd}"

cd solvers/fluid_sups/karman_vortex
. shell/meshgen_karman_vortex_2d.sh
cd ../../..

rm -r $directory
mkdir -p $directory
cd $directory

mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/node.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/elem.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/D_bc_v.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/D_bc_p.dat ./

# for parameteric study
rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#density_array" double 1 1.0
mv cond.dat density.dat
./../../../test_thermal/bin/cmd2cond "#viscosity_array" double 1 0.01
mv cond.dat viscosity.dat

./../../../test_thermal/bin/cmd2cond "#target_density" double 1 1.0
mv cond.dat target_density.dat
./../../../test_thermal/bin/cmd2cond "#target_viscosity" double 1 0.01
mv cond.dat target_viscosity.dat

./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 200 "#rom_output_interval" int 1 100
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.005 "#output_interval" int 1 100  "#finish_time" double 1 200

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n $nd
#./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc_p.dat -ig node.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc_v.dat -ig node.dat

cd ../..
