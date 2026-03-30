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

for nm in "${num_modes[@]}"
do
    . shell/diff/meshgen.sh $e $ep $nm 1 1 $pa
    . shell/diff/execution.sh $e $ep $nm 1 1 $pa $st
done

