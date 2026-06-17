#pragma once

#include "convdiff_core.h"
#include "nedelec_core.h"

void BBFE_std_shapefunc_hex1st_nedelec_get_val(
    const double xi[3],
    double** N_edge,
    const double J_inv[3][3]);

void BBFE_std_shapefunc_hex1st_nedelec_get_curl(
    const double xi[3],
    double **curl_N_edge,
    const double J[3][3],
    const double detJ);

/* --- Nedelec(1st kind, lowest order) の値 --- */
void BBFE_std_shapefunc_tet1st_nedelec_get_val(
    const double xi[3],
    double **N_edge,              /* [6][3] */
    const double J_inv[3][3]);

/* --- Nedelec(1st kind, lowest order) の curl --- */
void BBFE_std_shapefunc_tet1st_nedelec_get_curl(
    const double xi[3],           /* 形だけ合わせる（実際は一定） */
    double **curl_N_edge,         /* [6][3] */
    const double J[3][3],         /* [a(ref)][i(phys)] = ∂x_i/∂ξ_a */
    const double detJ);

static void tet2_W_ij(
    int i,
    int j,
    const double lam[4],
    const double glam[4][3],
    double W[3],
    double curlW[3]);

static void curl_scalar_times_vec_2nd(
    const double gradf[3],
    double f,
    const double V[3],
    const double curlV[3],
    double out[3]);

void BBFE_std_shapefunc_tet2nd_nedelec_get_ref(
    const double xi[3],
    const int elem_conn[4],
    double N[20][3],
    double curlN[20][3]);

void BBFE_std_shapefunc_tet2nd_nedelec_get_val_curl(
    const double xi[3],
    const int elem_conn[4],
    const double J[3][3],
    const double J_inv[3][3],
    const double detJ,
    double **N_phys,
    double **curlN_phys);
