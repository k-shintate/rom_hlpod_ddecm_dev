#pragma once

/*
 * AUTO-GENERATED from fluid_sups_core/core_FOM.c.
 * Do not hand-edit this file.
 *
 * The exact parameter spelling is copied from the verified provider
 * definitions so C++ name mangling cannot drift from the provider.
 */

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

double BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

void BBFE_elemmat_fluid_sups_mat_NR(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3]);

void BBFE_elemmat_fluid_sups_vec_NR(
    double         vec[4],
    const double   N_i,
    const double   grad_N_i[3],
    const double   v[3],
    const double   u_old[3],
    double**       grad_u,
    const double   p_cur,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3]);
