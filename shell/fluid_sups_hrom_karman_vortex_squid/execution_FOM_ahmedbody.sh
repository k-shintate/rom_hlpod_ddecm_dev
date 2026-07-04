#!/bin/bash

#並列計算領域数 (=並列数)
np=$1
#計算ノード数
N_node=$2
#計算ノード当たりのCPU数
N_cpu=$3

# 実行ディレクトリ
directory="result_fluid_sups_ahmedbody/FOM-${np}"

cd solvers/fluid_sups
make -f Makefile_FOM_squid clean
make -f Makefile_FOM_squid
cd ../..

fname="execution_hrom_offline_FOM_ahmedbody.sh"
rm $fname
touch $fname
echo "#!/bin/bash" >> $fname
echo "#------- qsub option -----------" >>$fname
echo "#PBS -q SQUID" >> $fname                          # debug == DBG, regular == SQUID
echo "#PBS --group=jh240017" >> $fname
echo "#PBS -b ${N_node}" >>$fname                       # node DBG == 2, SQUID == max 512
echo "#PBS -l elapstim_req=1:00:00" >>$fname            # debug == max 0:10:00, regular == max 24:00:00
echo "#PBS -l cpunum_job=${N_cpu}" >> $fname            # CPU/node == max 72
echo "#PBS -l memsz_job=248GB" >> $fname
echo "#PBS -T intmpi" >> $fname
echo "#------- Program execution -----------" >> $fname
echo "module load BaseCPU/2025" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname

echo "cd solvers/fluid_sups" >> $fname
echo "cp -r hlpod_fluid_sups_ahmedbody_FOM ./../../$directory" >> $fname
echo "cd ./../../$directory" >> $fname
echo "mkdir -p {pod_modes_vtk,pod_modes,fem_solver_prm,pod_solver_prm,hr_solver_prm,calctime,DDECM,hr_prm,hot_start}" >> $fname

echo "mpirun \${NQSV_MPIOPTS} -np ${np} ./hlpod_fluid_sups_ahmedbody_FOM ./  > ../../result/result_FOM-${np}.txt" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
