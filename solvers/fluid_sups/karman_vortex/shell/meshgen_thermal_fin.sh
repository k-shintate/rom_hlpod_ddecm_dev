 
#!/bin/bash

# inlet : constant velocity (1, 0, 0)
# outlet : traction free (fractional step method, sups method) + constant pressure (fractional step method)
# cylinder wall : non-slip wall
# right wall, left wall, top wall, bottom wall : slip wall

mkdir -p mesh_tmp
mkdir -p mesh_karman_vortex

#gmsh -3 -format msh2 -0 gmsh/thermal_fin_fig3.geo -o ./mesh_karman_vortex/thermal_fin.msh
#gmsh ./mesh_karman_vortex/thermal_fin.msh

#python3 mesh_io/save_physical_groups.py ./mesh_karman_vortex/thermal_fin.msh

### for node.dat
#mv ./mesh_tmp/Fluid_node_coordinates.dat ./mesh_karman_vortex/node.dat                      #node.dat
### for elem.dat
#mv ./mesh_tmp/Fluid_hexahedron_connectivity.dat ./mesh_karman_vortex/elem.dat               #elem.dat
#cp  ./mesh_tmp/Cylinder_wall_quad_connectivity.dat ./mesh_karman_vortex/
#mv  ./mesh_karman_vortex/Cylinder_wall_quad_connectivity.dat ./mesh_karman_vortex/surf.dat  #surf.dat

# --------------------------------------------------
# Repeat mesh + given_subdomain + velocity BC
# --------------------------------------------------

rm -rf merged_mesh

python3 mesh_io/repeat_dat_mesh_with_bc.py \
  --geo gmsh/thermal_fin_fig3.geo \
  --physical-group-script mesh_io/save_physical_groups.py \
  --copies 8 \
  --depth-copies 4 \
  --axis x \
  --gap 0.0 \
  --symmetric \
  --block-length 1 \
  --bc xmin:1 \
  --bc-value 0.0 \
  --out-dir merged_mesh

cd merged_mesh

cp D_bc.dat D_bc.dat.bak

awk '
NR==1 {print; next}
{print $1, 0, $3}
' D_bc.dat.bak > D_bc.dat

head D_bc.dat

cd ..

# --------------------------------------------------
# Resolve GEDATSU paths BEFORE cd
# --------------------------------------------------

GEDATSU_DIR=$(
    realpath \
    ../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin
)

MESH_PARTITIONER="$GEDATSU_DIR/gedatsu_simple_mesh_partitioner"

# 実際のファイル名を確認すること
BC_PARTITIONER="$GEDATSU_DIR/gedatsu_bc_partitioner_R"


# --------------------------------------------------
# Partition mesh
# --------------------------------------------------

cd merged_mesh

mkdir -p parted.0

"$MESH_PARTITIONER" \
    -n "32" \
    -in node.dat \
    -ie elem.dat \
    --given_subdomain given_subdomain.dat \
    -d parted.0


# --------------------------------------------------
# Partition velocity BC
# --------------------------------------------------

"$BC_PARTITIONER" \
    -n "32" \
    -i D_bc.dat \
    -ig node.dat

cd ..


