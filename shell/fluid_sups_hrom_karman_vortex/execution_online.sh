#!/bin/bash

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6
#solver type
st=$7

# 実行ディレクトリ
directory_offline="result_fluid_sups_karman_vortex/${nm}-${np}-${nd}"
directory_online="result_fluid_sups_karman_vortex/online_${nm}-${np}-${nd}"

rm -r $directory_online
#mkdir -p $directory_offline
mkdir -p result_fluid_sups_karman_vortex/tmp
cp -r result_fluid_sups_karman_vortex/${nm}-${np}-${nd} result_fluid_sups_karman_vortex/tmp
mv result_fluid_sups_karman_vortex/tmp/${nm}-${np}-${nd} $directory_online

cd solvers/fluid_sups

make -f Makefile_HROM_karman_vortex clean
make -f Makefile_HROM_karman_vortex

cp -r hlpod_fluid_sups_online_HROM ./../../$directory_online

cd ../../$directory_online

fname="gdb_cmd"
rm $fname
touch $fname
echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> $fname
echo "backtrace" >> $fname
echo "exit" >> $fname

#mpirun -np $np  ./hlpod_fluid_sups_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_offline_ROM
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_online_HROM
#mpirun -np $np  ./hlpod_fluid_sups_offline_ROM ./ -nd $nd -nm $nm -pa $pa -st $st
mpirun -np $np  ./hlpod_fluid_sups_online_HROM ./ -nd $nd -nm $nm -pa $pa -st $st

cd ../..