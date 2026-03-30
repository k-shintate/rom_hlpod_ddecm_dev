#!/bin/bash

#mesh
#一方向分割数
e=40
#解析領域の大きさ
ep=1

#podモード数
num_modes=(5)
#POD計算領域数
num_1stdd=(32)
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
    
        #. shell/mag/meshgen_threephase_tet.sh $e $ep $nm $nd $np $pa
        #. shell/mag/execution.sh $e $ep $nm $nd $np $pa $st
        #. shell/mag/graphgen_threephase.sh $e $ep $nm $nd $np $pa
        #python3 ./shell/mag/merge_graph.py result_mag/$nm-$np-$nd --elem graph_elem.dat --nedelec graph_nedelec_elem.dat --out graph.dat
        ##. shell/mag/partitioner.sh $e $ep $nm $nd $np $pa $st
        #. shell/mag/partitioner_hddm.sh $e $ep $nm $nd $np $pa $st
        #. shell/mag/partitioner_hddm2.sh $e $ep $nm $nd $np $pa $st
        ##. shell/mag/mesh_convert.sh $e $ep $nm $nd $np $pa $st

        . shell/mag/execution_threephase.sh $e $ep $nm $nd $np $pa $st
        #. shell/mag/execution_threephase_online.sh $e $ep $nm $nd $np $pa $st
        
        done
	done
done
