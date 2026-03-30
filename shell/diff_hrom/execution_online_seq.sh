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
directory_offline="result_diff/${nm}-${np}-${nd}"
directory_online="result_diff/online_${nm}-${np}-${nd}"

rm -r $directory_online
#mkdir -p $directory_offline
mkdir -p result_diff/tmp
cp -r result_diff/20-3-${nd} result_diff/tmp
mv result_diff/tmp/20-3-${nd} $directory_online

cd $directory_online
rm -r parted.0
mv parted.1 parted.0
cd ../..

cd solvers/diff

make -f Makefile_HROM clean
make -f Makefile_HROM

cp -r hlpod_diff_online_HROM ./../../$directory_online

cd ../../$directory_online

fname="gdb_cmd"
rm $fname
touch $fname
echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> $fname
echo "backtrace" >> $fname
echo "exit" >> $fname

#mpirun -np $np  ./hlpod_diff_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_diff_offline_ROM
mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_diff_online_HROM
#mpirun -np $np  ./hlpod_diff_offline_ROM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np $np  ./hlpod_diff_online_HROM ./ -nd $nd -nm $nm -pa $pa -st $st

cd ../..