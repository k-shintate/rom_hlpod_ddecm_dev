#!/bin/bash

#一方向分割数
e=$1
#解析領域の大きさ
ep=$2
#並列計算領域数 (=並列数)
np=$3
#計算ノード数
N_node=$4
#計算ノード当たりのCPU数
N_cpu=$5

# 実行ディレクトリ
directory="result_fluid_sups_karman_vortex/FOM_${e}-${ep}-${np}"

cd solvers/fluid_sups
make -f Makefile_FOM_weakbc clean
make -f Makefile_FOM_weakbc
cd ../..

cd solvers/fluid_sups
cp -r hlpod_fluid_sups_karman_vortex_FOM_weakbc ./../../$directory
cd ./../../$directory
mkdir -p {fem_solver_prm,calctime,hot_start}

mpirun  -np ${np} ./hlpod_fluid_sups_karman_vortex_FOM_weakbc ./ --wall-mode nitsche-noslip --wall-gamma 100 --wall-cdt 0.1
#mpirun  -np ${np} ./hlpod_fluid_sups_karman_vortex_FOM_weakbc ./ --wall-mode nitsche-spalding --wall-exchange opposite-node --wall-gamma 100 --wall-cdt 0.1

cd ../..
