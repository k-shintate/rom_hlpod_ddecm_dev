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
. run_karman_mixed.sh   
#. shell/meshgen_karman_vortex_2d.sh
cd ../../..

rm -r $directory
mkdir -p $directory
cd $directory

mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/graph.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/node.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/elem_prism.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/elem_tet.dat ./
#mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/D_bc_v.dat ./
mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/surf.dat ./



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

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.001 "#output_interval" int 1 1  "#finish_time" double 1 200 "#mat_max_iter" int 1 10000 "#mat_epsilon" double 1 1.0e-8

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $np
./../../utils/bin/elem2graph -ie surf.dat -og surf_graph.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i surf_graph.dat -ig graph.dat
./../../utils/bin/elem2graph -ie elem_prism.dat -og graph_elem_prism.dat
./../../utils/bin/elem2graph -ie elem_tet.dat -og graph_elem_tet.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i graph_elem_tet.dat -ig graph.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n $np -i graph_elem_prism.dat -ig graph.dat

python3 ./../../shell/mag/node2dist_val.py ./node.dat node_distval.dat ./graph.dat

python3 ./../../shell/mag/make_karman_mixed_bc.py ./node.dat D_bc_v.dat

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n $np -i node_distval.dat -ig graph.dat

./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n $np -i D_bc_v.dat -ig graph.dat

for i in $(seq 0 $((np-1)))
do
    python3 ../../shell/mag/distval2node.py ./parted.0/node_distval.dat.$i  ./parted.0/node.dat.$i
done

cd ../..