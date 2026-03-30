#!/bin/bash
# set -eu

# variable
domain1=${1}   # 基底取得領域数
domain2=${2}   # 並列計算領域数
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
gedatsu_simple_mesh2graph_convertor
gedatsu_nodal_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
gedatsu_connectivity_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
gedatsu_dist_val_partitioner_R -n ${domain1} -i node.dat -d ./parted.0/parted.1
gedatsu_bc_partitioner_R -n ${domain1} -i bc.dat -d ./parted.0/parted.1
gedatsu_bc_partitioner_R -n ${domain1} -i load.dat -d ./parted.0/parted.1

cd parted.0/parted.1
gedatsu_nodal_graph_partitioner -n ${domain2} -i metagraph.dat -d ../
cd ../..

if [ -e a.out ]; then
  make clean
fi
make
mpirun -np ${domain2} ./a.out ${domain1}

#mpirun -np ${domain2} ~/monolis_solid_r_shimizu/bin/monolis_solid -HDDM

