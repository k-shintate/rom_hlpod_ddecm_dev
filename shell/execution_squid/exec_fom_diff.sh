#!/bin/bash

#mesh
#一方向分割数
e=20
#解析領域の大きさ
ep=5

#並列計算領域数 (=並列数)
num_parallel=(2)

#計算ノード数
N_node1=1
#計算ノード当たりのCPU数
N_cpu1=2

for np in "${num_parallel[@]}"
do
#    . shell/diff_squid/meshgen_FOM.sh $e $ep $np $N_node1 $N_cpu1
    . shell/diff_squid/execution_FOM.sh $e $ep $np $N_node1 $N_cpu1
done
