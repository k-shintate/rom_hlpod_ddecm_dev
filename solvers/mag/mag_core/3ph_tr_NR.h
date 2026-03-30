#pragma once

#include "mag_dataset.h"
#include "convdiff_core.h"
#include "nedelec_core.h"
#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include <limits.h>


/* ============================================================
 * Constants & Material Properties (SI Units)
 * ============================================================ */
static const double SIGMA_COPPER = 5.77e7; /* [S/m] */
static const double SIGMA_CORE   = 3.72e3; /* [S/m] */
static const double FREQ_HZ      = 50.0;   /* [Hz] */
static const double CURRENT_AMP  = 1;    /* [A] */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double MU0    = 4.0 * M_PI * 1.0e-7;
static const double NU_LIN = 1.0 / (4.0 * M_PI * 1.0e-7); /* 1/mu0 */


/* ============================================================
 * Helper Functions & Structures
 * ============================================================ */

/* dot */
static inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* norm */
static inline double norm3(const double a[3]) {
    return sqrt(dot3(a,a));
}

/* --- Coil Configuration Structure --- */
typedef struct {
    double axis[3];
    double center[3];
    double turns;
    double area;
    double phase_shift;
} COIL_INFO;

static const double MM_TO_M = 0.001;

/* ============================================================
 * Coil info (prop 1..3 are coils)
 * ============================================================ */
static inline int get_coil_info(int elem_prop, COIL_INFO* info) {
    info->axis[0] = 0.0; info->axis[1] = 1.0; info->axis[2] = 0.0; /* Y-axis */
    info->turns = 100.0;

    /* PI*(9^2 - 4.5^2) = 190.85 mm^2 -> m^2 */
    info->area  = 190.85 * (MM_TO_M * MM_TO_M);

    double cy = 30.0 * MM_TO_M;
    double cz = 20.0 * MM_TO_M;

    if (elem_prop == 1) {
        info->center[0] = 2.5 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = 0.0;
        return 1;
    } else if (elem_prop == 2) {
        info->center[0] = 35.0 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = -2.0 * M_PI / 3.0;
        return 1;
    } else if (elem_prop == 3) {
        info->center[0] = 67.5 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = +2.0 * M_PI / 3.0;
        return 1;
    }
    return 0;
}

void eval_nu_and_dnudB(
    double Bmag,
    double* nu,
    double* dnudB);

double get_reluctivity_nu(
    double Bmag);

void debug_max_B_and_nu_core(
    FE_SYSTEM* sys,
    const double* x_curr,
    int n_dof_total,
    int it, int step, double t,
    const char* directory);

void update_Aphi_NR(
    double* x_curr,
    const double* delta,
    int n_dof_total,
    double relaxation);

void get_sigmas_for_prop(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi);

void compute_B_ip(
    const NEDELEC* ned, int e, int p,
    const double* x_curr, int nEdge,
    double B[3]);

void get_interp_coords(
    int e, int p,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    double x_ip[3]);

static inline double dot3(
    const double a[3],
    const double b[3]);

static inline double norm3(
    const double a[3]);

void set_element_mat_NR_Aphi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt);

void set_element_vec_NR_Aphi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time);

void apply_dirichlet_bc_for_A_and_phi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);

void log_accuracy_metrics(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    int step,
    double t,
    double dt,
    double P_input);
