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
directory="result_fluid_sups_karman_vortex/${nm}-${np}-${nd}"

. shell/install.sh

cd solvers/fluid_sups

make -f Makefile_HROM_karman_vortex clean
make -f Makefile_HROM_karman_vortex

cp -r hlpod_fluid_sups_offline_FOM ./../../$directory
cp -r hlpod_fluid_sups_offline_ROM ./../../$directory
cp -r hlpod_fluid_sups_online_HROM ./../../$directory

cd ./../../$directory

mkdir -p {pod_modes_vtk,pod_modes_v,pod_modes_p,fem_solver_prm,pod_solver_prm,hr_solver_prm,calctime,DDECM,hr_prm,hot_start}
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes_v/subdomain${i}"
done
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes_p/subdomain${i}"
done

fname="gdb_cmd"
rm $fname
touch $fname
echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> $fname
echo "backtrace" >> $fname
echo "exit" >> $fname

#./hlpod_fluid_sups_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st
mpirun -np $np  ./hlpod_fluid_sups_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_offline_FOM
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_offline_ROM
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_fluid_sups_online_HROM
#mpirun -np $np  ./hlpod_fluid_sups_offline_ROM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np $np  ./hlpod_fluid_sups_online_HROM ./ -nd $nd -nm $nm -pa $pa -st $st

cd ../..
