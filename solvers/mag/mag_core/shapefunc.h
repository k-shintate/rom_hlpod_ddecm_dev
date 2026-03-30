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