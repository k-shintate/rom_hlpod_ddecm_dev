#!/bin/bash

#mesh
#一方向分割数
e=15
#解析領域の大きさ
ep=5

#podモード数
num_modes=(5)
#POD計算領域数
num_1stdd=(3)
#並列計算領域数 (=並列数)
num_parallel=(1)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=2

#計算ノード数
N_node1=1
#計算ノード当たりのCPU数
N_cpu1=12

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
	
        . shell/diff_squid/meshgen.sh $e $ep $nm $nd $np $pa
	    . shell/diff_hrom_squid/execution_offline.sh $e $ep $nm $nd $np $pa $st
        . shell/diff_hrom_squid/execution_online.sh $e $ep $nm $nd $np $pa $st

        done
	done
done

