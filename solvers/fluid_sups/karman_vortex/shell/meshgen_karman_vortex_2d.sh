 
#!/bin/bash

# inlet : constant velocity (1, 0, 0)
# outlet : traction free (fractional step method, sups method) + constant pressure (fractional step method)
# cylinder wall : non-slip wall
# right wall, left wall, top wall, bottom wall : slip wall

mkdir -p mesh_tmp
mkdir -p mesh_karman_vortex

gmsh -3 -format msh2 -0 gmsh/karman_vortex_2d.geo -o ./mesh_karman_vortex/karman_vortex_2d.msh
gmsh ./mesh_karman_vortex/karman_vortex_2d.msh

#python3 mesh_io/save_physical_groups.py ./mesh_tmp/karman_vortex_2d_v2.msh
python3 mesh_io/karman_vortex_mesh_io.py ./mesh_karman_vortex/karman_vortex_2d.msh

### for node.dat
mv ./mesh_tmp/Fluid_node_coordinates.dat ./mesh_karman_vortex/node.dat
### for elem.dat
mv ./mesh_tmp/Fluid_hexahedron_connectivity.dat ./mesh_karman_vortex/elem.dat
cp  ./mesh_tmp/Cylinder_wall_quad_connectivity.dat ./mesh_karman_vortex/
mv  ./mesh_karman_vortex/Cylinder_wall_quad_connectivity.dat ./mesh_karman_vortex/surf.dat

### for D_bc
python3 mesh_io/elem_to_node.py ./mesh_tmp/Inlet_quad_connectivity.dat ./mesh_tmp/Inlet_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Right_wall_quad_connectivity.dat ./mesh_tmp/Right_wall_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Bottom_wall_quad_connectivity.dat ./mesh_tmp/Bottom_wall_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Left_wall_quad_connectivity.dat ./mesh_tmp/Left_wall_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Top_wall_quad_connectivity.dat ./mesh_tmp/Top_wall_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Outlet_quad_connectivity.dat ./mesh_tmp/Outlet_quad_node.dat
python3 mesh_io/elem_to_node.py ./mesh_tmp/Cylinder_wall_quad_connectivity.dat ./mesh_tmp/Cylinder_wall_quad_node.dat

### for D_bc_v.dat
python3 mesh_io/node_to_bc.py ./mesh_tmp/Inlet_quad_node.dat ./mesh_tmp/Inlet_quad_bc.dat 3 1.0 0 0
python3 mesh_io/node_to_bc.py ./mesh_tmp/Cylinder_wall_quad_node.dat  ./mesh_tmp/Cylinder_wall_quad_bc.dat 3 0 0 0
python3 mesh_io/node_to_bc.py ./mesh_tmp/Right_wall_quad_node.dat  ./mesh_tmp/Right_wall_quad_bc.dat 3 0 0 0
python3 mesh_io/node_to_bc.py ./mesh_tmp/Left_wall_quad_node.dat  ./mesh_tmp/Left_wall_quad_bc.dat 3 0 0 0
python3 mesh_io/node_to_bc.py ./mesh_tmp/Bottom_wall_quad_node.dat  ./mesh_tmp/Bottom_wall_quad_bc.dat 3 0 0 0
python3 mesh_io/node_to_bc.py ./mesh_tmp/Top_wall_quad_node.dat  ./mesh_tmp/Top_wall_quad_bc.dat 3 0 0 0

#slip wall
python3 mesh_io/fileter_index.py ./mesh_tmp/Right_wall_quad_bc.dat ./mesh_tmp/filtered_Right_wall_quad_bc.dat 0 2
python3 mesh_io/fileter_index.py ./mesh_tmp/Left_wall_quad_bc.dat ./mesh_tmp/filtered_Left_wall_quad_bc.dat 0 2
python3 mesh_io/fileter_index.py ./mesh_tmp/Bottom_wall_quad_bc.dat ./mesh_tmp/filtered_Bottom_wall_quad_bc.dat 0 1
python3 mesh_io/fileter_index.py ./mesh_tmp/Top_wall_quad_bc.dat ./mesh_tmp/filtered_Top_wall_quad_bc.dat 0 1

# b.c. のマージ (境界条件ファイル間で共有する節点がある場合、後から入力した節点の情報を優先する)
python3 mesh_io/merge_bc.py ./mesh_tmp/filtered_Right_wall_quad_bc.dat \
    ./mesh_tmp/Cylinder_wall_quad_bc.dat ./mesh_tmp/filtered_Left_wall_quad_bc.dat ./mesh_tmp/filtered_Bottom_wall_quad_bc.dat \
    ./mesh_tmp/filtered_Top_wall_quad_bc.dat ./mesh_tmp/Inlet_quad_bc.dat -d ./ -o mesh_tmp/merged_bc.dat

mv ./mesh_tmp/merged_bc.dat ./mesh_karman_vortex/D_bc_v.dat

### for D_bc_p.dat
python3 mesh_io/node_to_bc.py ./mesh_tmp/Outlet_quad_node.dat  ./mesh_tmp/Outlet_quad_bc.dat 1 1
mv ./mesh_tmp/Outlet_quad_bc.dat ./mesh_karman_vortex/D_bc_p.dat

#rm -r mesh_tmp
