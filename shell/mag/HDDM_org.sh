#!/bin/bash
# set -eu

# variable
domain1=16
domain2=4
echo "ndomain_basis : ${domain1}   ndomain_parallel : ${domain2}"
# command
if [ -e ./parted.0 ]; then
  if [ -n "$(ls ./parted.0)" ]; then
    rm parted.0/*
  fi
  if [ -n "$(ls ./parted.0/parted.1)" ]; then
    rm parted.0/parted.1/*
  fi
fi

./../../../../test_thermal/bin/meshgen_hex 30 30 30 1 1 1
./../../../../test_thermal/bin/surf_dbc_all 3 1.0 1.0 1.0
mv D_bc.dat bc.dat
./../../../../test_thermal/bin/surf_dbc_all 3 0.0 0.0 0.0
mv D_bc.dat load.dat

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n ${domain1} -d ./parted.0/parted.1
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor
mv node.dat node_orig.dat
python3 ./../../shell/mag/node2dist_val.py ./node_orig.dat node.dat 29791
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
./../../utils/bin/elem2graph -ie elem.dat -og connectivity.dat
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n ${domain1} -i connectivity.dat -ig graph.dat -d ./parted.0/parted.1
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n ${domain1} -i node.dat -d ./parted.0/parted.1
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${domain1} -i bc.dat -d ./parted.0/parted.1
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${domain1} -i load.dat -d ./parted.0/parted.1

cd parted.0/parted.1
./../../../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${domain2} -i metagraph.dat -d ../
cd ../..

#if [ -e a.out ]; then
#  make clean
#fi
#make

mpirun -np ${domain2} ./../../utils/load_balancing/hddm ${domain1}

#mpirun -np ${domain2} ~/monolis_solid_r_shimizu/bin/monolis_solid -HDDM

