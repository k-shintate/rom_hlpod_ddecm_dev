#!/bin/bash

#podモード数
ne=$1
#POD計算領域数
e=$2
#並列計算領域数 (=並列数)
np=$5
#計算ノード数
N_node=$8
#計算ノード当たりのCPU数
N_cpu=$9

# 実行ディレクトリ
directory="result_diff/FOM_${e}-${ep}-${np}"

cd solvers/diff
make -f Makefile_FOM_squid clean
make -f Makefile_FOM_squid
cd ../..

fname="execution_online_seq_hrom.sh"
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
echo "module load BaseCPU/2023" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname

echo "cd solvers/diff" >> $fname
echo "cp -r hlpod_diff_FOM ./../../$directory" >> $fname
echo "cd ./../../$directory" >> $fname

echo "mpirun \${NQSV_MPIOPTS} -np ${np} ./hlpod_diff_FOM ./ > ../../result/result_FOM_${nd}-${nm}.txt" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
