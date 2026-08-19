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
directory="result_diff/${nm}-${np}-${nd}"

. shell/install.sh

cd solvers/diff

make -f Makefile_HROM clean
make -f Makefile_HROM

cp -r hlpod_diff_offline_FOM ./../../$directory
cp -r hlpod_diff_offline_ROM ./../../$directory
cp -r hlpod_diff_online_HROM ./../../$directory

cp -r hrom_stage3_serial ./../../$directory
cp -r hrom_stage3_to_parallel ./../../$directory

cd ./../../$directory

mkdir -p {pod_modes_vtk,pod_modes,fem_solver_prm,hr_solver_prm,pod_solver_prm,calctime,DDECM,hr_prm}
for ((i=0; i<nd; i++))
do
    mkdir -p "pod_modes/subdomain${i}"
done

fname="gdb_cmd"
rm $fname
touch $fname
echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> $fname
echo "backtrace" >> $fname
echo "exit" >> $fname

mpirun -np $np  ./hlpod_diff_offline_FOM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_diff_offline_FOM
mpirun -np $np  ./hlpod_diff_offline_ROM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_diff_offline_ROM

#export MONOLIS_IADMM_NU0=0.10
#export MONOLIS_IADMM_POWER=1.0

#./hrom_stage3_serial ./ $np  10000 1.0e-8 1.0e-8 1.0e-8  200 100.0 1.0e-10 1.0e-6  1 3

operator_signature_rtol=1.0e-8
mode5_seed_tol=1.0e-4
mode5_defect_tol=1.0e-1
mode5_max_round=32
mode5_add_per_round=1
mode5_hard_fallback=0

for prm1 in 1.0e-3 1.0e-4 1.0e-5 1.0e-6 1.0e-7 1.0e-8
do
    for prm2 in 0 4 5
    do
        if [ "${prm2}" -eq 0 ]; then
            # Classical sparse-NNLS baseline: no output threshold/refit.
            weight_tol=0.0
            coverage_tau=1.0e-8   # unused by mode 0, but must remain positive
            threshold_refit=0
        else
            # Adaptive mode 5: the output support is judged using weight_tol.
            weight_tol=1.0e-8
            coverage_tau=1.0e-8   # used only if hard fallback is enabled
            threshold_refit=1
        fi

        ./hrom_stage3_serial ./ "${np}" \
            20000 \
            "${prm1}" \
            "${weight_tol}" \
            "${coverage_tau}" \
            200 \
            100.0 \
            1.0e-10 \
            1.0e-6 \
            1 \
            "${prm2}" \
            2 \
            "${operator_signature_rtol}" \
            "${threshold_refit}" \
            32 \
            "${mode5_seed_tol}" \
            "${mode5_defect_tol}" \
            "${mode5_max_round}" \
            "${mode5_add_per_round}" \
            "${mode5_hard_fallback}"

        save_dir="hrom_stage3_results/prm1_${prm1}_prm2_${prm2}"
        mkdir -p "${save_dir}/DDECM"

        cp DDECM/stage_reduction_summary.txt "${save_dir}/DDECM/"
        cp DDECM/stage3_selected_elem.txt "${save_dir}/DDECM/"
        cp DDECM/stage3_NNLS_residual.dat "${save_dir}/DDECM/"

        rm -f l2_error_rom.txt l2_error_fem_hrom.txt

        mpirun -np "${np}" ./hlpod_diff_online_HROM ./ \
            -nd "${nd}" -nm "${nm}" -pa "${pa}" -st "${st}"

        cp l2_error_rom.txt "${save_dir}/"
        cp l2_error_fem_hrom.txt "${save_dir}/"
    done
done




#./hrom_stage3_to_parallel ./ 3 parted.1 metagraph_parted.0

#mpirun -np ${np}  gdb --command=gdb_cmd ./hlpod_diff_online_HROM
#mpirun -np $np  ./hlpod_diff_offline_ROM ./ -nd $nd -nm $nm -pa $pa -st $st
#mpirun -np $np  ./hlpod_diff_online_HROM ./ -nd $nd -nm $nm -pa $pa -st $st

cd ../..