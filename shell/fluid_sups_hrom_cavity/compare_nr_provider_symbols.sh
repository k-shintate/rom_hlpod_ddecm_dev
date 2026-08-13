#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
cd "${repo_root}/solvers/fluid_sups"

provider="./fluid_sups_core/core_FOM_st_provider.o"
caller="./fluid_sups_core/fluid_sups_st_model.o"

if [ ! -f "${provider}" ] || [ ! -f "${caller}" ]; then
    echo "ERROR: build both provider and caller objects first." >&2
    exit 2
fi

symbols=(
  BBFE_elemmat_fluid_sups_coef_metric_tensor
  BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC
  BBFE_elemmat_fluid_sups_mat_NR
  BBFE_elemmat_fluid_sups_vec_NR
)

echo "============================================================"
echo "NR provider/caller symbol comparison"
echo "============================================================"

for s in "${symbols[@]}"; do
    echo
    echo "---- ${s}"
    echo "[provider raw]"
    nm -g "${provider}" | grep "${s}" || true
    echo "[provider demangled]"
    nm -g "${provider}" | grep "${s}" | c++filt || true
    echo "[caller raw]"
    nm -g "${caller}" | grep "${s}" || true
    echo "[caller demangled]"
    nm -g "${caller}" | grep "${s}" | c++filt || true
done
