#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
solver="${repo_root}/solvers/fluid_sups"
header="${solver}/fluid_sups_core/fluid_sups_nr_api_compat.h"
caller="${solver}/fluid_sups_core/fluid_sups_st_model.o"
provider="${solver}/fluid_sups_core/core_FOM_st_provider.o"
makefile="${solver}/Makefile_ST_cavity"

symbols=(
  BBFE_elemmat_fluid_sups_coef_metric_tensor
  BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC
  BBFE_elemmat_fluid_sups_mat_NR
  BBFE_elemmat_fluid_sups_vec_NR
)

for f in "${header}" "${provider}" "${makefile}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: missing: ${f}" >&2
        exit 2
    fi
done

echo "============================================================"
echo "[1] Header actually present on disk"
echo "============================================================"
grep -n -A18 -E \
'^(double|void)[[:space:]]+BBFE_elemmat_fluid_sups_(coef_metric_tensor|coef_metric_tensor_LSIC|mat_NR|vec_NR)' \
"${header}" || true

echo
echo "const lines:"
grep -n 'const double' "${header}" || true

if ! grep -q 'const double J_inv' "${header}"; then
    echo
    echo "ERROR: the on-disk compatibility header is still non-const." >&2
    exit 3
fi

# IMPORTANT:
# Do not grep for the literal text extern "C" anywhere in the file because
# comments such as 'Do NOT add extern "C"' are harmless and caused a false
# positive in the previous probe.
if grep -n -E '^[[:space:]]*extern[[:space:]]+"C"([[:space:]]*\{|[[:space:]]+)' "${header}"; then
    echo
    echo 'ERROR: an active extern "C" declaration/block is present,' >&2
    echo 'but core_FOM_st_provider.o uses C++ linkage.' >&2
    exit 4
fi

echo
echo 'Active extern "C" declarations: none'

echo
echo "============================================================"
echo "[2] Remove the stale caller object"
echo "============================================================"
rm -f "${caller}"
ls -l "${caller}" 2>/dev/null || echo "caller object removed"

echo
echo "============================================================"
echo "[3] Force rebuild ONLY fluid_sups_st_model.o"
echo "============================================================"
cd "${solver}"

if ! make -B -f Makefile_ST_cavity ./fluid_sups_core/fluid_sups_st_model.o; then
    make -B -f Makefile_ST_cavity fluid_sups_core/fluid_sups_st_model.o
fi

if [ ! -f "${caller}" ]; then
    echo "ERROR: caller object was not produced." >&2
    exit 5
fi

echo
echo "============================================================"
echo "[4] Raw provider/caller ABI comparison after rebuild"
echo "============================================================"

bad=0
for s in "${symbols[@]}"; do
    p="$(nm -g "${provider}" 2>/dev/null | grep "${s}" | grep ' T ' | head -n 1 | awk '{print $NF}')"
    c="$(nm -g "${caller}" 2>/dev/null | grep "${s}" | grep ' U ' | head -n 1 | awk '{print $NF}')"

    echo
    echo "---- ${s}"
    echo "provider: ${p}"
    echo "caller  : ${c}"

    if [ -n "${p}" ] && [ "${p}" = "${c}" ]; then
        echo "STATUS  : EXACT MATCH"
    else
        echo "STATUS  : MISMATCH"
        bad=1
        echo "provider demangled:"
        printf '%s\n' "${p}" | c++filt || true
        echo "caller demangled:"
        printf '%s\n' "${c}" | c++filt || true
    fi
done

if [ "${bad}" -eq 0 ]; then
    echo
    echo "============================================================"
    echo "All four caller symbols now exactly match the provider."
    echo "Proceed with:"
    echo "  cd ${solver}"
    echo "  make -f Makefile_ST_cavity dev_build"
    echo "============================================================"
    exit 0
fi

echo
echo "============================================================"
echo "[5] Mismatch persists: locate competing declarations"
echo "============================================================"

grep -R -n \
    --include='*.h' \
    --include='*.hpp' \
    --include='*.c' \
    --include='*.cpp' \
    -E 'BBFE_elemmat_fluid_sups_(coef_metric_tensor|coef_metric_tensor_LSIC|mat_NR|vec_NR)[[:space:]]*\(' \
    "${solver}/fluid_sups_core" \
    "${solver}/fluid_sups_matvec" || true

echo
echo "============================================================"
echo "[6] Includes used by fluid_sups_st_model.c"
echo "============================================================"
grep -n '^[[:space:]]*#[[:space:]]*include' \
    "${solver}/fluid_sups_core/fluid_sups_st_model.c" || true

exit 10
