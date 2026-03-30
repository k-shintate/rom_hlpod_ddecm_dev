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
directory="result_diff/${nm}-${np}-${nd}"

. shell/install.sh

cd solvers/diff

make
make

cp -r hlpod_diff ./../../$directory

cd ./../../$directory

mkdir -p {pod_modes_vtk,pod_modes,fem_solver_prm,pod_solver_prm,calctime}
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes/subdomain${i}"
done

mpirun -np $np  ./hlpod_diff ./ -nd $nd -nm $nm -pa $pa -st $st
cd ../..