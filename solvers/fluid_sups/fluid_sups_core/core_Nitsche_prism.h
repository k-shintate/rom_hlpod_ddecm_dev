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


/*
 * MPI surface-copy policy
 * -----------------------
 * Nitsche assembly routines below accept overlap/duplicate surface copies.
 * SUM-type physical postprocessing (wall area, force, Cd, Cl) must instead use
 * a unique-face surface or an equivalent ownership mask so each physical face
 * contributes exactly once.
 */

typedef struct {
    long long surface_faces;
    long long mapped_faces;
    long long skipped_no_local_owner;
    long long assembled_faces;
    long long quadrature_points;

    /* Assembly-copy area when the input surface contains overlap copies. */
    double area;
    double h_n_min;
    double h_n_max;
    double beta_min;
    double beta_mean;
    double beta_max;

    /* Split penalty diagnostics.  These are assembly-copy weighted stats. */
    double beta_mu_min;
    double beta_mu_mean;
    double beta_mu_max;

    double beta_dt_min;
    double beta_dt_mean;
    double beta_dt_max;

    double beta_dt_over_beta_mu_mean;
    double beta_dt_over_beta_mu_max;

    /* Internal weighted sums used by the MPI reducer. */
    double beta_area_sum;
    double beta_mu_area_sum;
    double beta_dt_area_sum;
    double beta_ratio_area_sum;

    double rhs_abs_sum;
    double matrix_abs_sum;
    double max_abs_normal_residual;
} PRI6NitscheDiagnostics;

/*
 * Physical/internal wall-surface diagnostics.
 *
 * IMPORTANT:
 *   The input surfaces to BBFE_mixed_wall_physical_surface_diagnostics_global()
 *   must be globally non-duplicated: each physical face exactly once.
 *   This is the surface class intended for SUM-type physical quantities.
 */
typedef struct {
    long long pri_faces;
    long long tet_faces;
    long long total_faces;

    double pri_area;
    double tet_area;
    double total_area;
} MixedWallPhysicalSurfaceDiagnostics;

int BBFE_mixed_wall_physical_surface_diagnostics_global(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_DATA* fe_pri,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_DATA* fe_tet,
    MONOLIS_COM* mono_com,
    MixedWallPhysicalSurfaceDiagnostics* diag);


/* -------------------------------------------------------------------------- */
/* Mixed Ahmed-wall physical diagnostics                                      */
/* -------------------------------------------------------------------------- */

/*
 * These diagnostics MUST use the non-duplicated/internal wall surfaces.
 * They are physical SUM/statistics and must not be evaluated on overlap copies.
 *
 * y_plus definition used here:
 *   tau_w  = | P_t [ mu (grad(u) + grad(u)^T) n ] |
 *   u_tau  = sqrt(tau_w / rho)
 *   y_plus = rho * u_tau * h_n / mu
 *
 * h_n is the wall-normal distance to the opposite PRI6 triangular face or the
 * opposite TET4 vertex.  Thus y_plus is the first-off-wall/exchange-layer y+.
 */
typedef struct {
    long long faces;
    long long quadrature_points;

    double area;

    double h_n_min;
    double h_n_mean;
    double h_n_max;

    double tau_w_min;
    double tau_w_mean;
    double tau_w_max;

    double u_tau_min;
    double u_tau_mean;
    double u_tau_max;

    double y_plus_min;
    double y_plus_mean;
    double y_plus_max;

    double abs_u_normal_mean;
    double abs_u_normal_max;

    double u_tangent_mean;
    double u_tangent_max;

    double pressure_min;
    double pressure_mean;
    double pressure_max;
} MixedWallFamilyDiagnostics;

typedef struct {
    MixedWallFamilyDiagnostics pri6;
    MixedWallFamilyDiagnostics tet4;
    MixedWallFamilyDiagnostics total;
} MixedWallDiagnostics;

int BBFE_mixed_wall_diagnostics_global(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_DATA* fe_tet,
    const BBFE_BASIS* basis_surf_tet,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    MONOLIS_COM* mono_com,
    MixedWallDiagnostics* diag);

void BBFE_mixed_reset_wall_diagnostics_file(
    const char* directory);

void BBFE_mixed_output_wall_diagnostics(
    const MixedWallDiagnostics* diag,
    double t,
    const char* directory);

/* Defaults mirror the existing TET4 wall implementation. */
void BBFE_pri6_wall_options_default(PRI6WallOptions* opt);

/* Runtime wall-mode helpers used by the Ahmed mixed executable. */
const char* BBFE_pri6_wall_mode_name(int mode);
int BBFE_pri6_wall_mode_from_string(const char* text, int* mode_out);
const char* BBFE_pri6_wall_exchange_name(int mode);
int BBFE_pri6_wall_exchange_from_string(const char* text, int* mode_out);
int BBFE_pri6_wall_mode_uses_nitsche(int mode);
int BBFE_pri6_wall_mode_uses_spalding(int mode);
int BBFE_pri6_wall_mode_uses_projected_strong_normal(int mode);
int BBFE_pri6_wall_mode_is_strong_noslip(int mode);

/* Feature angle used to cluster nodal normals for projected-strong modes. */
int BBFE_mixed_wall_set_feature_angle_deg(double feature_angle_deg);
double BBFE_mixed_wall_get_feature_angle_deg(void);

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
/* TET4/TRI3 wall support                                                     */
/* -------------------------------------------------------------------------- */

/*
 * Build TRI3 -> owner TET4 map.
 * owner_elem[es] is the local tetrahedron element id, or -1 when this rank
 * does not contain the owner tetrahedron.  owner_lface[es] is 0..3, or -1.
 */
int BBFE_tet4_tri3_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_tet,
    int** owner_elem,
    int** owner_lface,
    int* num_unmapped_local);

/*
 * Assemble the same symmetric no-slip/free-slip Nitsche operator used by the
 * PRI6 wall code, but on TRI3 faces owned by TET4 cells.
 */
int set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe_tet,
    const BBFE_BASIS* basis_tet,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    PRI6NitscheDiagnostics* diag);

/*
 * Mixed TET4+PRI6 Newton/VMS solver with the Ahmed wall split into
 *   - PRI6-owned TRI3 faces, and
 *   - TET4-owned TRI3 faces.
 * Both wall families are assembled into the same Monolis matrix/RHS before
 * the strong (non-Ahmed-wall) Dirichlet rows are applied.
 */

/* VMS/Newton counterpart of the dual-surface solver. Returns 0 only when the
 * nonlinear residual converges; on failure v_old is intentionally not advanced. */
int solver_fom_VMS_mixed_dual_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet,
    const BBFE_BASIS* basis_surf_tet,
    int enable_pri6_nitsche,
    const PRI6WallOptions* wall_opt_pri,
    int enable_tet4_nitsche,
    const PRI6WallOptions* wall_opt_tet,
    double t,
    int step);

void solver_fom_NR_mixed_dual_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet,
    const BBFE_BASIS* basis_surf_tet,
    int enable_pri6_nitsche,
    const PRI6WallOptions* wall_opt_pri,
    int enable_tet4_nitsche,
    const PRI6WallOptions* wall_opt_tet,
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

/*
 * Physical force SUM: surf MUST be a globally non-duplicated/internal PRI6
 * surface.  Do not pass an overlap-capable Nitsche assembly surface here.
 */
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
