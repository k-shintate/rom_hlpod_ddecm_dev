#!/bin/bash

#mesh
#一方向分割数
e=5
#解析領域の大きさ
ep=5

#並列計算領域数 (=並列数)
num_parallel=(8)

for np in "${num_parallel[@]}"
do
	. shell/fluid_sups_hrom_karman_vortex/meshgen_FOM_ahmedbody.sh $e $ep $np
	. shell/fluid_sups_hrom_karman_vortex/execution_FOM.sh $e $ep $np
done
