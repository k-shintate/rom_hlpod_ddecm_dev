#!/bin/bash

#mesh
#一方向分割数
e=20

#解析領域の大きさ
ep=5

#podモード数
num_modes=(20)
#POD計算領域数
num_1stdd=(6)
#並列計算領域数 (=並列数)
num_parallel=(6)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=3

#負荷分散後並列計算領域数
np2=(4)

#計算ノード数
N_node1=1
#計算ノード当たりのCPU数
N_cpu1=6
#計算ノード数
N_node2=1
#計算ノード当たりのCPU数
N_cpu2=6

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
	
#       . shell/diff_squid/meshgen.sh $e $ep $nm $nd $np $pa $N_node1 $N_cpu1
#	. shell/diff_squid/merge_graph.sh $e $ep $nm $nd $np $pa $N_node1 $N_cpu1
#	. shell/diff_hrom_squid/execution_offline_FOM.sh $e $ep $nm $nd $np $pa $st $N_node1 $N_cpu1
	. shell/diff_hrom_squid/execution_offline_ROM.sh $e $ep $nm $nd $np $pa $st $N_node1 $N_cpu1
#	. shell/diff_squid/merge_graph_lb.sh $e $ep $nm $nd $np $np2 $pa $N_node1 $N_cpu1
#	. shell/diff_hrom_squid/execution_online.sh $e $ep $nm $nd $np2 $pa $st $N_node1 $N_cpu1

	done
	done
done
