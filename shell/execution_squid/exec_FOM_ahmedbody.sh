#!/bin/bash

#並列計算領域数 (=並列数)
num_parallel=(72)

#計算ノード数
N_node1=1
#計算ノード当たりのCPU数
N_cpu1=72

for np in "${num_parallel[@]}"
do
    . shell/fluid_sups_hrom_karman_vortex_squid/meshgen_FOM.sh $np $N_node1 $N_cpu1
    . shell/fluid_sups_hrom_karman_vortex_squid/execution_FOM_ahmedbody.sh $np $N_node1 $N_cpu1
done
