#!/bin/bash
# set -eu

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

#rm -r $directory
mkdir -p $directory
cd $directory
cp -r ../../../../utils/load_balancing/hddm ./

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

#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n ${domain1} -d ./parted.0/parted.1
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_dist_val_partitioner_R -n ${domain1} -i node.dat -d ./parted.0/parted.1
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${domain1} -i bc.dat -d ./parted.0/parted.1
#./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${domain1} -i load.dat -d ./parted.0/parted.1

cd parted.0/parted.1
./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${domain2} -i metagraph.dat -d ../
cd ../..

if [ -e a.out ]; then
  make clean
fi
make

fname="gdb_cmd"
rm $fname
touch $fname
#echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> $fname
echo "run ./ 4" >> $fname
echo "backtrace" >> $fname
echo "exit" >> $fname

#mpirun -np $np  ./hlpod_fluid_sups_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st

#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_offline_ROM

mpirun -np ${domain2} ./../../hddm ${domain1}

cd ../..

#mpirun -np ${domain2} ~/monolis_solid_r_shimizu/bin/monolis_solid -HDDM

