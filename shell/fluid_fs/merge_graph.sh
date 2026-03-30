#!/bin/bash

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6

# 実行ディレクトリ
directory="result_fluid_fs/${nm}-${np}-${nd}"

cd $directory

#重みなし
cp -r ./parted.0/metagraph.dat ./
./../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n $np -i metagraph.dat -d metagraph_parted.0

mkdir merged_graph
mpirun -np $np ./../../utils/bin/merge_graph ./ -np1 $nd -np2 $np
mpirun -np $np ./../../utils/bin/merge_graph_bc ./ -np1 $nd -np2 $np -nbc 1 -fbc D_bc_p
mpirun -np $np ./../../utils/bin/merge_graph_bc ./ -np1 $nd -np2 $np -nbc 3 -fbc D_bc_v

# POD計算領域の分割ファイルをparted.1 として扱う
mv parted.0 parted.1
# マージ後の分割ファイル (並列計算領域に相当) をparted.0 として扱う (test_thermalやmonolisのデフォルトの分割ファイル名であるため)
mv merged_graph parted.0

cd ../..
