#!/bin/bash

#mesh
#一方向分割数
e=30
#解析領域の大きさ
ep=5

#podモード数
num_modes=(10)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=1

for nm in "${num_modes[@]}"
do
    . shell/diff/meshgen_ST.sh $e $ep $nm 4 4 $pa
    . shell/diff/execution_ST.sh $e $ep $nm 4 4 $pa $st
done
