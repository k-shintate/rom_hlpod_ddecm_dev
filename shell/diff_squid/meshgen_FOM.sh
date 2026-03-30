#!/bin/bash

#mesh
#一方向分割数
e=$1
#解析領域の大きさ
ep=$2
#並列計算領域数 (=並列数)
np=$3
#計算ノード数
N_node=$4
#計算ノード当たりのCPU数
N_cpu=$5

# 実行ディレクトリ
directory="result_diff/FOM_${e}-${ep}-${np}"

fname="meshgen.sh"
rm $fname
touch $fname
echo "#!/bin/bash" >> $fname
echo "#------- qsub option -----------" >>$fname
echo "#PBS -q SQUID" >> $fname                          # debug == DBG, regular == SQUID
echo "#PBS --group=jh240017" >> $fname
echo "#PBS -b 1" >>$fname                       # node DBG == 2, SQUID == max 512
echo "#PBS -l elapstim_req=1:00:00" >>$fname            # debug == max 0:10:00, regular == max 24:00:00
echo "#PBS -l cpunum_job=1" >> $fname            # CPU/node == max 72
echo "#PBS -l memsz_job=248GB" >> $fname
echo "#PBS -T intmpi" >> $fname
echo "#------- Program execution -----------" >> $fname
echo "module load BaseCPU/2023" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname


echo "rm -r $directory" >> $fname
echo "mkdir -p $directory" >> $fname
echo "cd $directory" >> $fname

echo "rm -r cond.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#snapshot_interval"' int 1 1 '"#rom_finish_time"' double 1 4.0 '"#rom_output_interval"' int 1 1" >> $fname
echo "mv cond.dat rom_cond.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#time_spacing"' double 1 0.01 '"#output_interval"' int 1 1 '"#finish_time"' double 1 1.0" >> $fname

echo "./../../../test_thermal/bin/meshgen_hex ${e} ${e} ${e} ${ep} ${ep} ${ep}" >> $fname
echo "./../../../test_thermal/bin/surf_dbc_all 1 1.0" >> $fname

echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n ${np}" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${np} -i D_bc.dat -ig node.dat" >> $fname

echo "cp -r parted.0 ./../../" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
