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

# 1. 直接 MSH2 形式でメッシュ生成 (-format msh2 を追加)
gmsh -3 -format msh2 ../../solvers/mag/threephase_transformer/gmsh/team21c.geo -o sphere2.msh
#gmsh sphere2.msh

# 2. Python スクリプトを実行
python3 ../../solvers/mag/threephase_transformer/mesh_io/save_physical_groups.py sphere2.msh
python3 ../../solvers/mag/threephase_transformer/mesh_io/karman_vortex_mesh_io.py sphere2.msh

mkdir mesh_tmp

python3 ../../solvers/mag/threephase_transformer/mesh_io/elem_to_node.py ./mesh_tmp/AIR_OUTER_WALLS_triangle_connectivity.dat ./mesh_tmp/OuterSphere_triangle_node.dat
python3 ../../solvers/mag/threephase_transformer/mesh_io/elem_to_node.py ./mesh_tmp/INTERFACE_AIR_SOLID_triangle_connectivity.dat ./mesh_tmp/InnerSphere_triangle_node.dat

mv ./mesh_tmp/All_node_coordinates.dat ./node.dat
mv ./mesh_tmp/All_tetra_connectivity.dat ./elem.dat
mv ./mesh_tmp/DOMAIN_tetra_connectivity.dat ./elem_air.dat
mv ./mesh_tmp/MAGNETIC_STEEL_tetra_connectivity.dat ./elem_iron.dat
mv ./mesh_tmp/SHIELD_tetra_connectivity.dat ./elem_shield.dat

mv ./mesh_tmp/COIL_1_tetra_connectivity.dat ./elem_widing1.dat
mv ./mesh_tmp/COIL_2_tetra_connectivity.dat ./elem_widing2.dat


python3 ../../solvers/mag/threephase_transformer/mesh_io/node_to_bc.py ./mesh_tmp/OuterSphere_triangle_node.dat ./D_bc.dat 1 0
python3 ../../solvers/mag/threephase_transformer/mesh_io/node_to_bc.py ./mesh_tmp/InnerSphere_triangle_node.dat ./InnerSphere_quad_bc.dat 1 0

cd ../..
