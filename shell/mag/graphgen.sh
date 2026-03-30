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

#rm -r $directory
mkdir -p $directory
cd $directory

./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i elem.dat -o graph_elem.dat
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh2graph_convertor -i nedelec_elem.dat -o graph_nedelec_elem.dat

cd ../..