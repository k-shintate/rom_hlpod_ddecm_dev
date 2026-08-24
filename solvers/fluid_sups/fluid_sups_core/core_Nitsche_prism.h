#pragma once

#include "core_FOM_mixed.h"
#include "core_NR.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * PRI6/TRI3 wall treatment for a boundary-layer prism mesh.
 *
 * PRI6 reference node convention used by core_FOM_mixed.c:
 *
 *   zeta=-1 : 0 -- 1 -- 2   (triangular face 0)
 *   zeta=+1 : 3 -- 4 -- 5   (triangular face 1)
 *
 * The three quadrilateral side faces are intentionally not candidates for
 * TRI3 ownership.  This matches a wall-normal prism layer where the wall is
 * one of the two triangular end faces.
 */

enum {
    PRI6_WALL_STRONG_NOSLIP          = 0,
    PRI6_WALL_NITSCHE_NOSLIP         = 1,
    PRI6_WALL_NITSCHE_FREESLIP       = 2,
    PRI6_WALL_NITSCHE_SPALDING       = 3,
    PRI6_WALL_STRONG_NORMAL_FREESLIP = 4,
    PRI6_WALL_STRONG_NORMAL_SPALDING = 5
};

enum {
    PRI6_WALL_EXCHANGE_WALL_Q        = 0,
    PRI6_WALL_EXCHANGE_OPPOSITE_FACE = 1
};

typedef struct {
    int wall_mode;
    int exchange_mode;

    double gamma_n;
    double dt_penalty_coeff;

    double kappa;
    double B;
    double ut_eps;
} PRI6WallOptions;

typedef struct {
    long long surface_faces;
    long long mapped_faces;
    long long skipped_no_local_owner;
    long long assembled_faces;
    long long quadrature_points;

    double area;
    double h_n_min;
    double h_n_max;
    double beta_min;
    double beta_max;

    double rhs_abs_sum;
    double matrix_abs_sum;
    double max_abs_normal_residual;
} PRI6NitscheDiagnostics;

/* Defaults mirror the existing TET4 wall implementation. */
void BBFE_pri6_wall_options_default(PRI6WallOptions* opt);

/*
 * Build TRI3 -> owner PRI6 map.
 * owner_elem[es] is local prism element id, or -1 if this MPI rank has only a
 * surface copy and does not own the corresponding prism.
 * owner_lface[es] is 0 (nodes 0,1,2), 1 (nodes 3,4,5), or -1.
 * Caller frees owner_elem/owner_lface.
 */
int BBFE_pri6_tri3_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    int** owner_elem,
    int** owner_lface,
    int* num_unmapped_local);

/*
 * Assemble symmetric Nitsche wall terms on TRI3 faces owned by PRI6 cells.
 * The volume SUPS/VMS terms must already have been assembled into monolis.
 */
int set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_pri,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    PRI6NitscheDiagnostics* diag);

/*
 * Mixed TET4 + PRI6 Newton/VMS solver for the common case where the selected
 * wall surface is entirely on the prism boundary layer.
 */
void solver_fom_NR_mixed_prism_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_wall,
    const BBFE_BASIS* basis_wall,
    const PRI6WallOptions* wall_opt,
    double t,
    int step);


/* -------------------------------------------------------------------------- */
/* Karman drag/lift diagnostics on TRI3 wall faces owned by PRI6 cells.       */
/* -------------------------------------------------------------------------- */

typedef struct {
    double area;
    double face_count;
    double face_qp_count;

    /* Body-force contributions in global x/y/z coordinates. */
    double force_pressure[3];

    /*
     * Legacy FOM viscous-force convention, consistent with the existing
     * SUPS volume viscous form.
     */
    double force_viscous[3];

    /* Symmetric physical Cauchy viscous-stress comparison. */
    double force_cauchy_viscous[3];

    /* Numerical Nitsche penalty reaction; intentionally reported separately. */
    double force_nitsche[3];

    /* Tangential Spalding wall-law reaction, when that mode is selected. */
    double force_walllaw[3];
} PRI6ForceDiagnostics;

int BBFE_pri6_force_diagnostics_global(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    MONOLIS_COM* mono_com,
    PRI6ForceDiagnostics* diag);

void BBFE_pri6_reset_karman_force_files(
    const char* directory);

void BBFE_pri6_output_karman_force_diagnostics(
    const PRI6ForceDiagnostics* diag,
    double t,
    double rho,
    double U_ref,
    double A_ref,
    const double eU[3],
    const double eL[3],
    const char* directory);

#ifdef __cplusplus
}
#endif
