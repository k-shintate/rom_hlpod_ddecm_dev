#!/bin/bash

#mesh
#一方向分割数
e=20
#解析領域の大きさ
ep=1

#podモード数
num_modes=(20)
#POD計算領域数
num_1stdd=(5)
#並列計算領域数 (=並列数)
num_parallel=(5)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type (POD計算領域数 = 並列計算領域数)
st=4 

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
            . shell/fluid_sups/meshgen.sh $e $ep $nm $nd $np $pa
            #. shell/fluid_sups/merge_graph.sh $e $ep $nm $nd $np $pa
            . shell/fluid_sups/execution.sh $e $ep $nm $nd $np $pa $st
        done
	done
done

