#!/bin/bash

#mesh
#一方向分割数
e=$1
#解析領域の大きさ
ep=$2
np=$3

# 実行ディレクトリ
directory="result_fluid_sups_karman_vortex/FOM_${e}-${ep}-${np}"

cd solvers/fluid_sups/karman_vortex
. shell/meshgen_karman_vortex_2d.sh
cd ../../..

rm -r $directory
mkdir -p $directory
cd $directory

mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/node.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/elem.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/D_bc_v.dat ./

mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/surf.dat ./

#cp -r ../../node.dat ./
#cp -r ../../elem.dat ./
#cp -r ../../D_bc_v.dat ./
#cp -r ../../D_bc_p.dat ./

# for parameteric study
rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#density_array" double 1 200
mv cond.dat density.dat
./../../../test_thermal/bin/cmd2cond "#viscosity_array" double 1 1
mv cond.dat viscosity.dat

./../../../test_thermal/bin/cmd2cond "#target_density" double 1 200
mv cond.dat target_density.dat
./../../../test_thermal/bin/cmd2cond "#target_viscosity" double 1 1
mv cond.dat target_viscosity.dat

./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 200 "#rom_output_interval" int 1 1000
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.001 "#output_interval" int 1 1000  "#finish_time" double 1 200 "#mat_max_iter" int 1 10000 "#mat_epsilon" double 1 1.0e-8

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n $np
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $np
./../../utils/bin/elem2graph -ie surf.dat -og surf_graph.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i surf_graph.dat -ig graph.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $np -i D_bc_v.dat -ig node.dat

cd ../..