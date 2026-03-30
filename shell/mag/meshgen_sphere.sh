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

rm -r $directory
mkdir -p $directory
cd $directory

rm -r cond.dat
./../../../test_thermal/bin/cmd2cond "#snapshot_interval" int 1 1 "#rom_finish_time" double 1 4.0 "#rom_output_interval" int 1 10
mv cond.dat rom_cond.dat

./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 0.001 "#output_interval" int 1 1  "#finish_time" double 1 1.0

gmsh -3 ../../plot/gmsh/sphere_air.geo -o sphere.msh
mv sphere.msh sphere.unv
gmsh -format msh2 -0 sphere.unv -o sphere2.msh
python3 ../../plot/mesh_io/save_physical_groups.py sphere2.msh
python3 ../../plot/mesh_io/karman_vortex_mesh_io.py sphere2.msh

mkdir mesh_tmp

python3 ../../plot/mesh_io/elem_to_node.py ./mesh_tmp/OuterBoundary_quad_connectivity.dat ./mesh_tmp/OuterSphere_quad_node.dat
python3 ../../plot/mesh_io/elem_to_node.py ./mesh_tmp/OuterSphere_quad_connectivity.dat ./mesh_tmp/InnerSphere_quad_node.dat

mv ./mesh_tmp/All_node_coordinates.dat ./node.dat
mv ./mesh_tmp/All_hexahedron_connectivity.dat ./elem.dat
mv ./mesh_tmp/Shell_hexahedron_connectivity.dat ./elem_inner.dat
mv ./mesh_tmp/OuterAir_hexahedron_connectivity.dat ./elem_air.dat

python3 ../../plot/mesh_io/node_to_bc.py ./mesh_tmp/OuterSphere_quad_node.dat ./D_bc.dat 1 0
python3 ../../plot/mesh_io/node_to_bc.py ./mesh_tmp/InnerSphere_quad_node.dat ./InnerSphere_quad_bc.dat 1 0

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i elem.dat -o graph_elem.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i nedelec_elem.dat -o graph_nedelec_elem.dat

cd ../..
