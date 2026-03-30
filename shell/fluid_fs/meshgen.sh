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
directory="result_fluid_fs/${nm}-${np}-${nd}"

rm -r $directory
mkdir -p $directory
cd $directory

# for parameteric study
rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#density_array" double 1 100.0
mv cond.dat density.dat
./../../../test_thermal/bin/cmd2cond "#viscosity_array" double 1 1
mv cond.dat viscosity.dat

./../../../test_thermal/bin/cmd2cond "#target_density" double 1 100.0
mv cond.dat target_density.dat
./../../../test_thermal/bin/cmd2cond "#target_viscosity" double 1 1
mv cond.dat target_viscosity.dat

./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 1.0 "#rom_output_interval" int 1 1
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.01 "#output_interval" int 1 1  "#finish_time" double 1 1.0

./../../../test_thermal/bin/meshgen_hex $e $e $e 1.0 1.0 1.0

./../../../test_thermal/submodule/monolis/bin/monolis_extract_all_surf_hex
# -> surf.datを出力

# 一様流速境界部およびnonslip境界部の切り出し
./../../../test_thermal/bin/mesh_surf_extract 0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_z.dat -ov surf_z.vtk
./../../../test_thermal/bin/mesh_surf_remove  0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_wall.dat -ov surf_wall.vtk

# 速度のDirechlet B.C.ファイルの作成
./../../../test_thermal/bin/surf_dbc 3 1.0 0.0 0.0 -ie surf_z.dat -o D_bc_z.dat
./../../../test_thermal/bin/surf_dbc 3 0.0 0.0 0.0 -ie surf_wall.dat -o D_bc_wall.dat
#./../../../../test_thermal/bin/surf_dbc 3 0.0 0.0 0.0 -ie surf_wall.dat -o D_bc_wall.dat
./../../../test_thermal/bin/surf_bc_merge D_bc_z.dat D_bc_wall.dat -o D_bc_v.dat

# 圧力に関する境界部の切り出し
./../../../test_thermal/bin/mesh_surf_extract 0.48 0.48 0.0 0.52 0.52 0.0 -oe surf_p.dat -ov surf_p.vtk

# 圧力のDirichlet B.C.ファイルの作成
./../../../test_thermal/bin/surf_dbc 1 0.0 -ie surf_p.dat -o D_bc_p.dat

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n $nd
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc_p.dat -ig node.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $nd -i D_bc_v.dat -ig node.dat

cd ../..