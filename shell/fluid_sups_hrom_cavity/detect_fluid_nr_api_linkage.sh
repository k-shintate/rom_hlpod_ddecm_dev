#!/bin/bash
set -euo pipefail

mode_only=0
if [ "${1:-}" = "--mode" ]; then
    mode_only=1
    shift
fi

root="${1:-../../../test_thermal}"
libdir="${root}/lib"

symbols=(
  BBFE_elemmat_fluid_sups_coef_metric_tensor
  BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC
  BBFE_elemmat_fluid_sups_mat_NR
  BBFE_elemmat_fluid_sups_vec_NR
  BBFE_elemmat_fluid_sups_vec_NR_mass_con
  BBFE_elemmat_fluid_sups_vec_NR_momentum
)

libs=()
while IFS= read -r f; do
    libs+=("$f")
done < <(
    find "${libdir}" -maxdepth 1 -type f \
      \( -name 'libBBFE_elemmat.a' \
      -o -name 'libBBFE_elemmat.dylib' \
      -o -name 'libBBFE_elemmat.so' \) -print 2>/dev/null
)

if [ "${#libs[@]}" -eq 0 ]; then
    [ "${mode_only}" -eq 1 ] && echo "missing"
    echo "ERROR: libBBFE_elemmat was not found in ${libdir}" >&2
    exit 3
fi

tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT

for lib in "${libs[@]}"; do
    nm -g "${lib}" >> "${tmp}" 2>/dev/null || nm "${lib}" >> "${tmp}" 2>/dev/null || true
done

num_c=0
num_cpp=0
num_missing=0

for s in "${symbols[@]}"; do
    # Exact C symbol. Mach-O prepends one underscore; ELF generally does not.
    if awk '{print $NF}' "${tmp}" | grep -Eq "^_?${s}$"; then
        num_c=$((num_c + 1))
        continue
    fi

    # A C++ mangled symbol still contains the original function identifier.
    if grep -Fq "${s}" "${tmp}"; then
        num_cpp=$((num_cpp + 1))
        continue
    fi

    num_missing=$((num_missing + 1))
done

if [ "${num_missing}" -eq 0 ] && [ "${num_c}" -eq "${#symbols[@]}" ]; then
    mode="c"
elif [ "${num_missing}" -eq 0 ] && [ "${num_cpp}" -eq "${#symbols[@]}" ]; then
    mode="cpp"
elif [ "${num_missing}" -gt 0 ]; then
    mode="missing"
else
    mode="mixed"
fi

if [ "${mode_only}" -eq 1 ]; then
    echo "${mode}"
    exit 0
fi

echo "BBFE SUPS/NR API linkage diagnosis"
echo "  library dir : ${libdir}"
echo "  C symbols   : ${num_c}/${#symbols[@]}"
echo "  C++ symbols : ${num_cpp}/${#symbols[@]}"
echo "  missing     : ${num_missing}/${#symbols[@]}"
echo "  mode        : ${mode}"
echo

for s in "${symbols[@]}"; do
    echo "---- ${s}"
    grep -F "${s}" "${tmp}" | head -n 5 || echo "NOT FOUND"
done

if [ "${mode}" = "missing" ]; then
    echo >&2
    echo "The linked libBBFE_elemmat does not contain all verified NR functions." >&2
    echo "Search for the exact definitions with:" >&2
    echo "  grep -R -n -E '^[[:space:]]*(void|double)[[:space:]]+BBFE_elemmat_fluid_sups_(mat_NR|vec_NR|vec_NR_mass_con|vec_NR_momentum|coef_metric_tensor)' ${root}" >&2
    exit 4
fi

if [ "${mode}" = "mixed" ]; then
    echo "ERROR: mixed C/C++ linkage was detected; inspect the symbols above." >&2
    exit 5
fi
