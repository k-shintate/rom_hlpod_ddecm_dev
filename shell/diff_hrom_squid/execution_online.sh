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
#計算ノード数
N_node=$8
#計算ノード当たりのCPU数
N_cpu=$9

# 実行ディレクトリ
directory_offline="result_diff/${nm}-${np}-${nd}"
directory_online="result_diff/online_${nm}-${np}-${nd}"

cd solvers/diff
make -f Makefile_HROM_squid clean
make -f Makefile_HROM_squid
cd ../..

fname="execution_online_hrom.sh"
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

#echo "rm -r $directory_online" >> $fname
#echo "mkdir -p result_diff/tmp" >> $fname
#echo "cp -r result_diff/${nm}-${np}-${nd} result_diff/tmp" >> $fname
#echo "mv result_diff/tmp/${nm}-${np}-${nd} $directory_online" >> $fname

echo "cd solvers/diff" >> $fname
echo "cp -r hlpod_diff_online_HROM ./../../$directory_online" >> $fname
echo "cd ../../$directory_online" >> $fname
echo "mkdir hr_solver_prm" >> $fname
echo "mpirun \${NQSV_MPIOPTS} -np ${np}  ./hlpod_diff_online_HROM ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}  > ../../result/result_online_${nd}-${nm}.txt" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
