#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
dst="${repo_root}/solvers/fluid_sups/fluid_sups_core/fluid_sups_nr_api_compat.h"

if [ ! -d "$(dirname "${dst}")" ]; then
    echo "ERROR: solver directory not found: $(dirname "${dst}")" >&2
    exit 2
fi

if [ -f "${dst}" ] && [ ! -f "${dst}.before_const_exact_fix" ]; then
    cp "${dst}" "${dst}.before_const_exact_fix"
fi

cat > "${dst}" <<'EOF'
#pragma once

/*
 * Exact C++ declarations matching core_FOM_st_provider.o.
 * Do NOT add extern "C".
 */

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    double Jacobian,
    double density,
    double viscosity,
    const double* velocity,
    double dt);

double BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
    const double J_inv[3][3],
    double Jacobian,
    double density,
    double viscosity,
    const double* velocity,
    double dt);

void BBFE_elemmat_fluid_sups_mat_NR(
    double element_matrix[4][4],
    const double J_inv[3][3],
    double N_i,
    double N_j,
    const double* grad_N_i,
    const double* grad_N_j,
    const double* velocity,
    double** grad_velocity,
    const double* grad_pressure,
    double density,
    double viscosity,
    double tau,
    double tau_c,
    double dt,
    const double* delta_velocity_time);

void BBFE_elemmat_fluid_sups_vec_NR(
    double* element_vector,
    double N_i,
    const double* grad_N_i,
    const double* velocity,
    const double* velocity_old,
    double** grad_velocity,
    double pressure,
    const double* grad_pressure,
    double density,
    double viscosity,
    double tau,
    double tau_c,
    double dt,
    const double* delta_velocity_time);
EOF

echo "patched: ${dst}"
echo
echo "No extern C:"
grep -n 'extern "C"' "${dst}" && {
    echo "ERROR: extern C unexpectedly present" >&2
    exit 3
} || true

echo
echo "const-qualified declarations:"
grep -n 'const double' "${dst}"
