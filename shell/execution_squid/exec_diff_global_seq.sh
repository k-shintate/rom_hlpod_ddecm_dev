#!/bin/bash

#mesh
#一方向分割数
e=30
#解析領域の大きさ
ep=5

#podモード数
num_modes=(10 20)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=1

#計算ノード数
N_node1=1
#計算ノード当たりのCPU数
N_cpu1=12
#計算ノード数
N_node2=1
#計算ノード当たりのCPU数
N_cpu2=12


for nm in "${num_modes[@]}"
do
    . shell/diff_squid/meshgen.sh $e $ep $nm 1 1 $pa $N_node1 $N_cpu1
    . shell/diff_squid/execution.sh $e $ep $nm 1 1 $pa $st $N_node1 $N_cpu1
done
