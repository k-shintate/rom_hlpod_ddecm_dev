#!/bin/bash

#mesh
#一方向分割数
e=$1
#解析領域の大きさ
ep=$2

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6

#計算ノード数
N_node=$7
#計算ノード当たりのCPU数
N_cpu=$8

# 実行ディレクトリ
directory="result_fluid_sups_cavity/${nm}-${np}-${nd}"

fname="meshgen.sh"
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

echo "rm -r $directory" >> $fname
echo "mkdir -p $directory" >> $fname
echo "cd $directory" >> $fname

# for parameteric study
echo "rm -r cond.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#density_array"' double 1 100.0" >> $fname
echo "mv cond.dat density.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#viscosity_array"' double 1 1" >> $fname
echo "mv cond.dat viscosity.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#target_density"' double 1 100.0" >> $fname
echo "mv cond.dat target_density.dat" >> $fname
echo "./../../../test_thermal/bin/cmd2cond '"#target_viscosity"' double 1 1" >> $fname
echo "mv cond.dat target_viscosity.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#snapshot_interval"' int 1 1 '"#rom_finish_time"' double 1 1.0 '"#rom_output_interval"' int 1 10" >> $fname
echo "mv cond.dat rom_cond.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#inc_svd_interval"' int 1 50" >> $fname
echo "mv cond.dat hrom_cond.dat" >> $fname

echo "./../../../test_thermal/bin/cmd2cond '"#time_spacing"' double 1 0.005 '"#output_interval"' int 1 10  '"#finish_time"' double 1 1.0" >> $fname

echo "./../../../test_thermal/bin/meshgen_hex ${e} ${e} ${e} 1.0 1.0 1.0" >> $fname

echo "./../../../test_thermal/submodule/monolis/bin/monolis_extract_all_surf_hex" >> $fname
# -> surf.datを出力

echo "./../../../test_thermal/bin/mesh_surf_extract 0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_z.dat -ov surf_z.vtk" >> $fname
echo "./../../../test_thermal/bin/mesh_surf_remove  0.0 0.0 1.0 1.0 1.0 1.0 -oe surf_wall.dat -ov surf_wall.vtk" >> $fname

# 速度のDirechlet B.C.ファイルの作成
echo "./../../../test_thermal/bin/surf_dbc 3 1.0 0.0 0.0 -ie surf_z.dat -o D_bc_z.dat" >> $fname
echo "./../../../test_thermal/bin/surf_dbc 3 0.0 0.0 0.0 -ie surf_wall.dat -o D_bc_wall.dat" >> $fname
echo "./../../../test_thermal/bin/surf_bc_merge D_bc_z.dat D_bc_wall.dat -o D_bc_v.dat" >> $fname

# 圧力に関する境界部の切り出し
echo "./../../../test_thermal/bin/mesh_surf_extract 0.48 0.48 0.0 0.52 0.52 0.0 -oe surf_p.dat -ov surf_p.vtk" >> $fname

# 圧力のDirichlet B.C.ファイルの作成
echo "./../../../test_thermal/bin/surf_dbc 1 0.0 -ie surf_p.dat -o D_bc_p.dat" >> $fname

echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_simple_mesh_partitioner -n ${nd}" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${nd} -i D_bc_p.dat -ig node.dat" >> $fname
echo "./../../../../test_thermal/submodule/monolis/submodule/gedatsu/bin/gedatsu_bc_partitioner_R -n ${nd} -i D_bc_v.dat -ig node.dat" >> $fname

echo "cd ../.." >> $fname
echo "" >> $fname

qsub $fname
