 
#!/bin/bash

# inlet : constant velocity (1, 0, 0)
# outlet : traction free (fractional step method, sups method) + constant pressure (fractional step method)
# cylinder wall : non-slip wall
# right wall, left wall, top wall, bottom wall : slip wall

mkdir -p mesh_tmp
mkdir -p mesh_karman_vortex

gmsh -3 -format msh2 -0 gmsh/backstep.geo -o ./mesh_karman_vortex/backstep.msh
#gmsh ./mesh_karman_vortex/backstep.msh

python3 mesh_io/save_physical_groups.py ./mesh_karman_vortex/backstep.msh

### for node.dat
mv ./mesh_tmp/Fluid_node_coordinates.dat ./mesh_karman_vortex/node.dat                      #node.dat
### for elem.dat
mv ./mesh_tmp/Fluid_hexahedron_connectivity.dat ./mesh_karman_vortex/elem.dat               #elem.dat
cp  ./mesh_tmp/Nitsche_walls_quad_connectivity.dat ./mesh_karman_vortex/
mv  ./mesh_karman_vortex/Nitsche_walls_quad_connectivity.dat ./mesh_karman_vortex/surf.dat  #surf.dat

### for bc  要素コネクティビティファイルを境界条件ファイルに変換
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Inlet_quad_connectivity.dat ./mesh_tmp/Inlet_quad_bc.dat 3 1.0 0 0                 #流入境界条件
#python3 mesh_io/elem_to_bc.py ./mesh_tmp/Cylinder_wall_quad_connectivity.dat  ./mesh_tmp/Cylinder_wall_quad_bc.dat 3 0 0 0  #nonslip wall
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Top_wall_quad_connectivity.dat  ./mesh_tmp/Top_wall_quad_bc.dat 3 0 0 0        #slip wall #一時的にnonslip wall の処置
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Upstream_lower_wall_quad_connectivity.dat  ./mesh_tmp/Upstream_lower_wall_quad_bc.dat 3 0 0 0          #slip wall
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Bottom_wall_quad_connectivity.dat  ./mesh_tmp/Bottom_wall_quad_bc.dat 3 0 0 0      #slip wall
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Step_wall_quad_connectivity.dat  ./mesh_tmp/Step_wall_quad_bc.dat 3 0 0 0            #slip wall
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Front_quad_connectivity.dat  ./mesh_tmp/Front_quad_bc.dat 3 0 0 0      #slip wall
python3 mesh_io/elem_to_bc.py ./mesh_tmp/Back_quad_connectivity.dat  ./mesh_tmp/Back_quad_bc.dat 3 0 0 0            #slip wall

# 境界条件の削除
python3 mesh_io/fileter_index.py ./mesh_tmp/Top_wall_quad_bc.dat ./mesh_tmp/filtered_Top_wall_quad_bc.dat 0 2           #slip wall の場合、接線方向の境界条件を削除
python3 mesh_io/fileter_index.py ./mesh_tmp/Upstream_lower_wall_quad_bc.dat ./mesh_tmp/filtered_Upstream_lower_wall_quad_bc.dat 0 2
python3 mesh_io/fileter_index.py ./mesh_tmp/Bottom_wall_quad_bc.dat ./mesh_tmp/filtered_Bottom_wall_quad_bc.dat 0 2
python3 mesh_io/fileter_index.py ./mesh_tmp/Step_wall_quad_bc.dat ./mesh_tmp/filtered_Step_wall_quad_bc.dat 1 2
python3 mesh_io/fileter_index.py ./mesh_tmp/Front_quad_bc.dat ./mesh_tmp/filtered_Front_quad_bc.dat 0 1
python3 mesh_io/fileter_index.py ./mesh_tmp/Back_quad_bc.dat ./mesh_tmp/filtered_Back_quad_bc.dat 0 1

# b.c. のマージ (境界条件ファイル間で共有する節点がある場合、後から入力した節点の情報を優先する)
#python3 mesh_io/merge_bc.py ./mesh_tmp/filtered_Top_wall_quad_bc.dat ./mesh_tmp/filtered_Upstream_lower_wall_quad_bc.dat ./mesh_tmp/filtered_Bottom_wall_quad_bc.dat \
#    ./mesh_tmp/filtered_Step_wall_quad_bc.dat ./mesh_tmp/Inlet_quad_bc.dat -d ./ -o mesh_tmp/merged_bc.dat

#python3 mesh_io/merge_bc.py ./mesh_tmp/filtered_Top_wall_quad_bc.dat ./mesh_tmp/filtered_Upstream_lower_wall_quad_bc.dat ./mesh_tmp/filtered_Bottom_wall_quad_bc.dat \
#    ./mesh_tmp/filtered_Step_wall_quad_bc.dat ./mesh_tmp/filtered_Front_quad_bc.dat ./mesh_tmp/filtered_Back_quad_bc.dat ./mesh_tmp/Inlet_quad_bc.dat -d ./ -o mesh_tmp/merged_bc.dat


mv ./mesh_tmp/Inlet_quad_bc.dat mesh_tmp/merged_bc.dat

#python3 mesh_io/merge_bc.py ./mesh_tmp/filtered_Front_quad_bc.dat ./mesh_tmp/filtered_Back_quad_bc.dat ./mesh_tmp/Inlet_quad_bc.dat -d ./ -o mesh_tmp/merged_bc.dat
#python3 mesh_io/merge_bc.py ./mesh_tmp/filtered_Front_quad_bc.dat ./mesh_tmp/filtered_Back_quad_bc.dat ./mesh_tmp/Inlet_quad_bc.dat -d ./ -o mesh_tmp/merged_bc.dat


mv ./mesh_tmp/merged_bc.dat ./mesh_karman_vortex/D_bc_v.dat

