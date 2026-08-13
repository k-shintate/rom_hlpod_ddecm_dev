#!/bin/bash
set -euo pipefail

project_root="${1:-$(pwd)}"

if [ ! -d "${project_root}/solvers/fluid_sups" ]; then
    echo "ERROR: project root does not contain solvers/fluid_sups: ${project_root}" >&2
    exit 2
fi

solver_dir="${project_root}/solvers/fluid_sups"
test_thermal="${project_root}/../test_thermal"

symbols=(
    BBFE_elemmat_fluid_sups_coef_metric_tensor
    BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC
    BBFE_elemmat_fluid_sups_mat_NR
    BBFE_elemmat_fluid_sups_vec_NR
    BBFE_elemmat_fluid_sups_vec_NR_mass_con
    BBFE_elemmat_fluid_sups_vec_NR_momentum
)

echo "============================================================"
echo "Fluid NR provider locator"
echo "project root : ${project_root}"
echo "solver dir   : ${solver_dir}"
echo "test_thermal : ${test_thermal}"
echo "============================================================"
echo

echo "[1] Source references under the current fluid solver"
for s in "${symbols[@]}"; do
    echo
    echo "---- ${s}"
    grep -R -n \
        --include='*.c' \
        --include='*.cc' \
        --include='*.cpp' \
        --include='*.cxx' \
        --include='*.h' \
        --include='*.hpp' \
        "${s}" "${solver_dir}" 2>/dev/null || true
done

echo
echo "[2] Source references under test_thermal"
if [ -d "${test_thermal}" ]; then
    for s in "${symbols[@]}"; do
        echo
        echo "---- ${s}"
        grep -R -n \
            --include='*.c' \
            --include='*.cc' \
            --include='*.cpp' \
            --include='*.cxx' \
            --include='*.h' \
            --include='*.hpp' \
            "${s}" "${test_thermal}" 2>/dev/null || true
    done
else
    echo "test_thermal directory not found: ${test_thermal}"
fi

echo
echo "[3] Existing object files under solvers/fluid_sups"
found_object_provider=0
while IFS= read -r obj; do
    matched=0
    for s in "${symbols[@]}"; do
        if nm -g "${obj}" 2>/dev/null | awk '{print $NF}' | grep -Eq "^_?${s}$|${s}"; then
            if [ "${matched}" -eq 0 ]; then
                echo
                echo "OBJECT PROVIDER CANDIDATE: ${obj}"
                matched=1
                found_object_provider=1
            fi
            nm -g "${obj}" 2>/dev/null | grep "${s}" || true
        fi
    done
done < <(find "${solver_dir}" -type f -name '*.o' -print 2>/dev/null)

if [ "${found_object_provider}" -eq 0 ]; then
    echo "No current .o file exports the requested names."
fi

echo
echo "[4] Libraries/archives under project tree and test_thermal"
search_roots=("${project_root}")
if [ -d "${test_thermal}" ]; then
    search_roots+=("${test_thermal}")
fi

tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT

found_lib_provider=0
while IFS= read -r lib; do
    : > "${tmp}"
    nm -g "${lib}" > "${tmp}" 2>/dev/null || nm "${lib}" > "${tmp}" 2>/dev/null || true
    matched=0
    for s in "${symbols[@]}"; do
        if grep -q "${s}" "${tmp}"; then
            if [ "${matched}" -eq 0 ]; then
                echo
                echo "LIBRARY PROVIDER CANDIDATE: ${lib}"
                matched=1
                found_lib_provider=1
            fi
            grep "${s}" "${tmp}" | head -n 5 || true
        fi
    done
done < <(
    find "${search_roots[@]}" -type f \
        \( -name '*.a' -o -name '*.dylib' -o -name '*.so' \) \
        -print 2>/dev/null
)

if [ "${found_lib_provider}" -eq 0 ]; then
    echo "No library/archive under the searched roots exports the requested names."
fi

echo
echo "[5] Check the local core sources that are present in the legacy HROM build"
for f in \
    "${solver_dir}/fluid_sups_core/core_NR.c" \
    "${solver_dir}/fluid_sups_core/core_FOM.c" \
    "${solver_dir}/fluid_sups_core/fluid_core.c"
do
    if [ -f "${f}" ]; then
        echo
        echo "FILE: ${f}"
        grep -n -E \
            'BBFE_elemmat_fluid_sups_(coef_metric_tensor|coef_metric_tensor_LSIC|mat_NR|vec_NR|vec_NR_mass_con|vec_NR_momentum)' \
            "${f}" || true
    fi
done

echo
echo "============================================================"
echo "Interpretation"
echo "============================================================"
echo "A) If core_NR.o/core_NR.c defines the six BBFE_elemmat_* functions:"
echo "   add that provider object/source to ST_FOM_OBJ."
echo
echo "B) If another libBBFE_elemmat.* contains them:"
echo "   point Makefile_ST_cavity to that validated library revision."
echo
echo "C) If only call sites are found and no definition/provider is found:"
echo "   the validated NR code depends on a BBFE source revision that is not"
echo "   present in the current build tree. Recover/rebuild that exact source."
echo
echo "Do NOT replace *_NR with the older non-NR SUPS routines."
