#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}/solvers/fluid_sups"

provider="./fluid_sups_core/core_FOM_st_provider.o"
caller="./fluid_sups_core/fluid_sups_st_model.o"

symbols=(
  BBFE_elemmat_fluid_sups_coef_metric_tensor
  BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC
  BBFE_elemmat_fluid_sups_mat_NR
  BBFE_elemmat_fluid_sups_vec_NR
)

bad=0

for s in "${symbols[@]}"; do
    p="$(nm -g "${provider}" 2>/dev/null | grep "${s}" | grep ' T ' | head -n 1 | awk '{print $NF}')"
    c="$(nm -g "${caller}" 2>/dev/null | grep "${s}" | grep ' U ' | head -n 1 | awk '{print $NF}')"

    echo "---- ${s}"
    echo "provider: ${p}"
    echo "caller  : ${c}"

    if [ -z "${p}" ] || [ -z "${c}" ]; then
        echo "STATUS  : missing"
        bad=1
    elif [ "${p}" = "${c}" ]; then
        echo "STATUS  : EXACT MATCH"
    else
        echo "STATUS  : MISMATCH"
        echo "provider demangled:"
        printf '%s\n' "${p}" | c++filt
        echo "caller demangled:"
        printf '%s\n' "${c}" | c++filt
        bad=1
    fi
    echo
done

exit "${bad}"
