#!/bin/bash

nd=$1
#計算ノード数
N_node=$2
#計算ノード当たりのCPU数
N_cpu=$3

# 実行ディレクトリ
directory="result_fluid_sups_ahmedbody/FOM-${np}"

fname="meshgen.sh"
rm $fname
touch $fname
echo "#!/bin/bash" >> $fname
echo "#------- qsub option -----------" >>$fname
echo "#PBS -q SQUID" >> $fname                          # debug == DBG, regular == SQUID
echo "#PBS --group=jh240017" >> $fname
echo "#PBS -b 1" >>$fname                       # node DBG == 2, SQUID == max 512
echo "#PBS -l elapstim_req=0:30:00" >>$fname            # debug == max 0:10:00, regular == max 24:00:00
echo "#PBS -l cpunum_job=1" >> $fname            # CPU/node == max 72
echo "#PBS -l memsz_job=248GB" >> $fname
echo "#PBS -T intmpi" >> $fname
echo "#------- Program execution -----------" >> $fname
echo "module load BaseCPU/2025" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname

#echo "cd solvers/fluid_sups/karman_vortex" >> $fname
#echo ". shell/meshgen_karman_vortex_2d.sh" >> $fname
#echo "cd ../../.." >> $fname

echo "rm -r $directory" >> $fname
echo "mkdir -p $directory" >> $fname
echo "cd $directory" >> $fname

#echo "mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/node.dat ./" >> $fname
#echo "mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/elem.dat ./" >> $fname
#echo "mv ../../solvers/fluid_sups/karman_vortex/mesh_karman_vortex/D_bc_v.dat ./" >> $fname

echo "cp -r ../../node.dat ./" >> $fname
echo "cp -r ../../elem.dat ./" >> $fname
echo "cp -r ../../D_bc_v.dat ./" >> $fname
echo "cp -r ../../surf.dat ./" >> $fname
echo "cp -r ../../surf_graph.dat ./" >> $fname
echo "cp -r ../../graph.dat ./" >> $fname

# for parameteric study
echo "rm -r cond.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#density_array"' double 1 0.68435754189e5" >> $fname
echo "mv cond.dat density.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#viscosity_array"' double 1 1" >> $fname
echo "mv cond.dat viscosity.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#target_density"' double 1 0.68435754189e5" >> $fname
echo "mv cond.dat target_density.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#target_viscosity"' double 1 1" >> $fname
echo "mv cond.dat target_viscosity.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#snapshot_interval"' int 1 1000 '"#rom_finish_time"' double 1 0.5 '"#rom_output_interval"' int 1 1000" >> $fname
echo "mv cond.dat rom_cond.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond "#time_spacing" double 1 8.0e-5 "#output_interval" int 1 1  "#finish_time" double 1 0.5 "#mat_max_iter" int 1 100000 "#mat_epsilon" double 1 1.0e-8 "#density" double 1 0.68435754189e5 "#viscosity" double 1 1
" >> $fname

echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n ${nd}" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_nodal_graph_partitioner -n ${nd}" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_connectivity_graph_partitioner -n ${nd} -i surf_graph.dat -ig graph.dat" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n ${nd}" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${nd} -i D_bc_v.dat -ig node.dat" >> $fname


echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname