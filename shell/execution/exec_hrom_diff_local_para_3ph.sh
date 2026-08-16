#!/bin/bash

#mesh
#一方向分割数
e=20

#解析領域の大きさ
ep=5

#podモード数
num_modes=(10)
#POD計算領域数
num_1stdd=(72)
#並列計算領域数 (=並列数)
num_parallel=(16)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=3

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
	
        . shell/diff_hrom/meshgen.sh $e $ep $nm $nd $np $pa
	    . shell/diff_hrom/merge_graph.sh $e $ep $nm $nd $np $pa
	    . shell/diff_hrom/execution_offline_3ph.sh $e $ep $nm $nd $np $pa $st

        done
	done
done
