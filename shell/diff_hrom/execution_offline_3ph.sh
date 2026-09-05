#!/bin/bash

# Stage 3 mode 0 uses MPI-parallel MONOLIS NNLS, but its output responsibility
# is kept identical to the conventional serial Stage 3:
#   DDECM/lb_selected_elem*.txt
# The existing online-HROM path performs the subsequent rank-wise conversion.

nm=$3
nd=$4
np=$5
pa=$6
st=$7

directory="result_diff/${nm}-${np}-${nd}"

. shell/install.sh

cd solvers/diff
make -f Makefile_HROM clean
make -f Makefile_HROM

cp -r hlpod_diff_offline_FOM ./../../"$directory"
cp -r hlpod_diff_offline_ROM ./../../"$directory"
cp -r hlpod_diff_online_HROM ./../../"$directory"
cp -r hrom_stage3_serial ./../../"$directory"
cp -r hrom_stage3_parallel_mode0 ./../../"$directory"

cd ./../../"$directory"

mkdir -p {pod_modes_vtk,pod_modes,fem_solver_prm,hr_solver_prm,pod_solver_prm,calctime,DDECM,hr_prm}
for ((i=0; i<nd; i++)); do
    mkdir -p "pod_modes/subdomain${i}"
done

fname="gdb_cmd"
rm -f "$fname"
touch "$fname"
echo "run ./ -nd ${nd} -nm ${nm} -pa ${pa} -st ${st}" >> "$fname"
echo "backtrace" >> "$fname"
echo "exit" >> "$fname"

mpirun -np "$np" ./hlpod_diff_offline_FOM ./ -nd "$nd" -nm "$nm" -pa "$pa" -st "$st"
mpirun -np "$np" ./hlpod_diff_offline_ROM ./ -nd "$nd" -nm "$nm" -pa "$pa" -st "$st"

operator_signature_rtol=1.0e-8
mode5_seed_tol=1.0e-4
mode5_defect_tol=1.0e-1
mode5_max_round=32
mode5_add_per_round=1
mode5_hard_fallback=0

for prm1 in 1.0e-12
do
    for prm2 in 0
    do
        stage3_log="stage3_prm1_${prm1}_prm2_${prm2}.log"

        if [ "${prm2}" -eq 0 ]; then
            weight_tol=0.0

            #mpirun -np "${np}" ./hrom_stage3_parallel_mode0 \
            #    ./ \
            #    "${np}" \
            #    20000 \
            #    "${prm1}" \
            #    0 \
            #    "${weight_tol}" \
            #    2>&1 | tee "${stage3_log}"

            stage3_rc=${PIPESTATUS[0]}
            if [ "${stage3_rc}" -ne 0 ]; then
                echo "WARNING: Stage 3 parallel mpirun returned rc=${stage3_rc}; script continues." >&2
            fi

            if grep -q '\[STAGE3-MODE0-LEGACY-DONE\]' "${stage3_log}" 2>/dev/null; then
                echo "Stage 3 legacy lb_selected_elem output: OK"
            else
                echo "WARNING: Stage 3 legacy output completion log was not found." >&2
            fi
        else
            weight_tol=1.0e-8
            coverage_tau=1.0e-8
            threshold_refit=1

            ./hrom_stage3_serial ./ "${np}" \
                20000 "${prm1}" "${weight_tol}" "${coverage_tau}" \
                200 100.0 1.0e-10 1.0e-6 1 "${prm2}" 2 \
                "${operator_signature_rtol}" "${threshold_refit}" 32 \
                "${mode5_seed_tol}" "${mode5_defect_tol}" \
                "${mode5_max_round}" "${mode5_add_per_round}" \
                "${mode5_hard_fallback}" \
                2>&1 | tee "${stage3_log}"
        fi

        save_dir="hrom_stage3_results/prm1_${prm1}_prm2_${prm2}"
        mkdir -p "${save_dir}/DDECM"

        for f in DDECM/stage_reduction_summary.txt DDECM/stage3_NNLS_residual.dat
        do
            [ -f "$f" ] && cp "$f" "${save_dir}/DDECM/"
        done

        # Save exactly the Stage-3 legacy outputs used by the conventional flow.
        for f in DDECM/lb_selected_elem.*.txt DDECM/lb_selected_elem_D_bc.*.txt
        do
            [ -f "$f" ] && cp "$f" "${save_dir}/DDECM/"
        done

        [ -f "${stage3_log}" ] && cp "${stage3_log}" "${save_dir}/"

        rm -f l2_error_rom.txt l2_error_fem_hrom.txt

        # Important: remove rank-wise files from earlier runs so the online HROM
        # cannot accidentally use stale data.  The conventional online path must
        # regenerate them from the new lb_selected_elem*.txt files.
        rm -f DDECM/selected_elem.[0-9]*.txt
        rm -f DDECM/selected_elem_D_bc.[0-9]*.txt

        echo "Starting conventional online HROM from Stage-3 lb_selected_elem files"
        mpirun -np "${np}" ./hlpod_diff_online_HROM ./ \
            -nd "${nd}" -nm "${nm}" -pa "${pa}" -st "${st}"

        online_rc=$?
        if [ "${online_rc}" -ne 0 ]; then
            echo "WARNING: online HROM returned rc=${online_rc}; script continues." >&2
        fi

        # The conventional online path should have generated these rank-wise files.
        rank_files_ok=1
        for ((r=0; r<np; r++)); do
            if [ ! -f "DDECM/selected_elem.${r}.txt" ]; then
                echo "WARNING: online path did not create DDECM/selected_elem.${r}.txt" >&2
                rank_files_ok=0
            fi
            if [ ! -f "DDECM/selected_elem_D_bc.${r}.txt" ]; then
                echo "WARNING: online path did not create DDECM/selected_elem_D_bc.${r}.txt" >&2
                rank_files_ok=0
            fi
        done
        if [ "${rank_files_ok}" -eq 1 ]; then
            echo "Conventional selected_elem rank-file generation: OK"
        fi

        # Save files after online conversion too.
        for f in DDECM/selected_elem.[0-9]*.txt DDECM/selected_elem_D_bc.[0-9]*.txt
        do
            [ -f "$f" ] && cp "$f" "${save_dir}/DDECM/"
        done

        [ -f l2_error_rom.txt ] && cp l2_error_rom.txt "${save_dir}/"
        [ -f l2_error_fem_hrom.txt ] && cp l2_error_fem_hrom.txt "${save_dir}/"
    done
done

cd ../..
