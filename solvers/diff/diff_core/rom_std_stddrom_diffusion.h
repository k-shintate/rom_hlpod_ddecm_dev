#ifndef ROM_STD_STDDROM_DIFFUSION_H
#define ROM_STD_STDDROM_DIFFUSION_H

/*
 * dG(0) all-at-once ST-DDROM adapter for the unsteady diffusion equation
 *
 *     M dT/dt + K T = F.
 *
 * One reduced unknown is attached to each time window.  Inside a window the
 * dG(0) slabs are retained as a single space-time vector in node-major order
 *
 *     [space node][slab].
 *
 * The full-order window operators are
 *
 *     A_w = block lower bidiagonal(D, -L),
 *     B_w : previous-window last slab -> current-window first slab,
 *
 *     D = M/dt + K,  L = M/dt.
 *
 * The reduced all-window equation solved by rom_std_stddrom is
 *
 *     D_r a_w - L_r a_{w-1} = g_r,w.
 *
 * Time-dependent Dirichlet data are handled by an affine lift T = ell + q;
 * the POD basis is homogeneous on Dirichlet nodes.
 */

#include "core_FOM.h"
#include "rom_std_stddrom.h"
#include "st_fom_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    FE_SYSTEM* fom;
    const ST_FOM_WINDOW_LAYOUT* layout;
    ROM_STDD_SYSTEM* reduced;

    MONOLIS window_operator;
    MONOLIS temporal_operator;
    MONOLIS_COM self_com;
    int operators_initialized;

    double* operator_input;
    double* lift;
    double* previous_input;
    double* source;
    double* work_A;
    double* work_B;
    double* effective_rhs;
} ROM_STDD_DIFFUSION_CONTEXT;

int ROM_std_stdd_diffusion_context_initialize(
    ROM_STDD_DIFFUSION_CONTEXT* context,
    FE_SYSTEM* fom,
    const ST_FOM_WINDOW_LAYOUT* layout,
    ROM_STDD_SYSTEM* reduced);

void ROM_std_stdd_diffusion_context_finalize(
    ROM_STDD_DIFFUSION_CONTEXT* context);

/* Build time-invariant dG(0) full-order A and B operators. */
int ROM_std_stdd_diffusion_build_operators(
    ROM_STDD_DIFFUSION_CONTEXT* context);

/* Matrix-free callbacks used by ROM_std_stdd_project_operators_rank_sweep. */
int ROM_std_stdd_diffusion_apply_A(
    const double* x,
    double* y,
    void* user_context);

int ROM_std_stdd_diffusion_apply_B(
    const double* x,
    double* y,
    void* user_context);

/* Make the POD basis exactly homogeneous on Dirichlet rows. */
int ROM_std_stdd_diffusion_enforce_homogeneous_basis(
    ROM_STDD_DIFFUSION_CONTEXT* context);

/*
 * Fill system->reduced_rhs for every online window using
 *
 *   g_0 = f_0 - A ell_0 + B u_initial,
 *   g_w = f_w - A ell_w + B ell_{w-1},  w > 0.
 */
int ROM_std_stdd_diffusion_project_all_rhs(
    ROM_STDD_DIFFUSION_CONTEXT* context);

/* Build the nodal Dirichlet lift for one window. */
int ROM_std_stdd_diffusion_build_lift(
    const ROM_STDD_DIFFUSION_CONTEXT* context,
    int window_id,
    double* lift_out);

/*
 * Reconstruct one complete local window T = ell + Phi a.  Halo values are
 * exchanged before the lift is added.
 */
int ROM_std_stdd_diffusion_reconstruct_window(
    ROM_STDD_DIFFUSION_CONTEXT* context,
    int window_id,
    double* full_window_out);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STDDROM_DIFFUSION_H */
