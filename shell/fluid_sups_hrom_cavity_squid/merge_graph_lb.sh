#!/bin/bash

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$7

np2=$6
#計算ノード数
N_node=$8
#計算ノード当たりのCPU数
N_cpu=$9

# 実行ディレクトリ
directory_offline="result_fluid_sups_cavity/${nm}-${np}-${nd}"
directory_online="result_fluid_sups_cavity/online_${nm}-${np2}-${nd}"

fname="merge_graph.sh"
rm $fname
touch $fname
echo "#!/bin/bash" >> $fname
echo "#------- qsub option -----------" >>$fname
echo "#PBS -q DBG" >> $fname                          # debug == DBG, regular == SQUID
echo "#PBS --group=jh240017" >> $fname
echo "#PBS -b ${N_node}" >>$fname                       # node DBG == 2, SQUID == max 512
echo "#PBS -l elapstim_req=0:10:00" >>$fname            # debug == max 0:10:00, regular == max 24:00:00
echo "#PBS -l cpunum_job=${N_cpu}" >> $fname            # CPU/node == max 72
echo "#PBS -l memsz_job=248GB" >> $fname
echo "#PBS -T intmpi" >> $fname
echo "#------- Program execution -----------" >> $fname
echo "module load BaseCPU/2023" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname

echo "rm -r $directory_online" >> $fname
echo "mkdir -p result_fluid_sups_cavity/tmp" >> $fname
echo "cp -r result_fluid_sups_cavity/${nm}-${np}-${nd} result_fluid_sups_cavity/tmp" >> $fname
echo "mv result_fluid_sups_cavity/tmp/${nm}-${np}-${nd} $directory_online" >> $fname

echo "cd solvers/fluid_sups" >> $fname
echo "cp -r hlpod_fluid_sups_online_HROM ./../../$directory_online" >> $fname
echo "cd ../../$directory_online" >> $fname

# ノード重みあり
echo "mkdir merged_graph" >> $fname
echo "./../../utils/bin/merge_graph ./ -np1 ${nd} -np2 ${np2}" >> $fname
echo "cp -r ./merged_graph/node_weight.dat ./" >> $fname
echo "cp -r ./parted.0/metagraph.dat ./" >> $fname
echo "rm -r metagraph_parted.0" >> $fname

#echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${np2} -i metagraph.dat -d metagraph_parted.0 -inw node_weight.dat" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${np2} -i metagraph.dat -d metagraph_parted.0" >> $fname

echo "rm -r parted.0" >> $fname
echo "mv parted.1 parted.0" >> $fname

echo "mpirun \${NQSV_MPIOPTS} -np ${np2} ./../../utils/bin/merge_graph ./ -np1 ${nd} -np2 ${np2}" >> $fname

# POD計算領域の分割ファイルをparted.1 として扱う
echo "mv parted.0 parted.1" >> $fname
# マージ後の分割ファイル (並列計算領域に相当) をparted.0 として扱う (test_thermalやmonolisのデフォルトの分割ファイル名であるため)
echo "mv merged_graph parted.0" >> $fname

echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${np2} -i D_bc_v.dat -ig node.dat" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${np2} -i D_bc_p.dat -ig node.dat" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
