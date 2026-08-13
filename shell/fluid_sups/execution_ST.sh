#!/bin/bash

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6
#solver type
st=$7

# 実行ディレクトリ
directory="result_fluid_sups/${nm}-${np}-${nd}"

. shell/install.sh

cd solvers/fluid_sups

make
make

cp -r hlpod_fluid_sups_offline ./../../$directory
cp -r hlpod_fluid_sups_online ./../../$directory

cd ./../../$directory

mkdir -p {pod_modes_vtk,pod_modes_v,pod_modes_p,fem_solver_prm,pod_solver_prm,calctime}
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes_v/subdomain${i}"
done
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes_p/subdomain${i}"
done

mpirun -np $np  ./hlpod_fluid_sups_offline ./ -nd $nd -nm $nm -pa 0 -st $st
mpirun -np $np  ./hlpod_fluid_sups_online ./ -nd $nd -nm $nm -pa $pa -st $st

cd ../..
