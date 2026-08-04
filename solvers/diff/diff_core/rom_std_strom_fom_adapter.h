#ifndef ROM_STD_STROM_FOM_ADAPTER_H
#define ROM_STD_STROM_FOM_ADAPTER_H

#include "core_ROM.h"
#include "core_FOM_ST.h"
#include "rom_std_strom.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Serial adapter between the existing windowed time-dG FOM and the
 * window-local space-time ROM.
 *
 * The previous-window quantity expected by set_element_vec_ST_dG_window()
 * is not a single nodal right-trace vector.  It is the complete temporal
 * state of the previous window's last slab:
 *
 *     T_prev_st[alpha * num_space_nodes + i],
 *
 * whose length is num_time_dofs * num_space_nodes.
 */
typedef struct {
    FE_SYSTEM* sys;
    const ST_TIME_ELEMENT* te;

    int num_windows;
    int num_window_slabs;
    int num_space_nodes;
    int num_time_dofs;
    int previous_state_dof;
    int st_dof;

    double* zero_previous_state;
    double* previous_state;
    double* rhs_with_previous;
    double* rhs_without_previous;

    /* Preintegrated element mass matrices used by the direct B action. */
    int local_num_nodes;
    double* element_mass;
} ROM_STROM_FOM_CONTEXT;

int ROM_std_strom_fom_context_initialize(
    ROM_STROM_FOM_CONTEXT* context,
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows);

void ROM_std_strom_fom_context_finalize(
    ROM_STROM_FOM_CONTEXT* context);

/* Construct the boundary-modified matrix A_w used by the serial FOM. */
int ROM_std_strom_fom_prepare_operator(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id);

/* Callback wrappers matching ROM_STROM_PREPARE/APPLY_PREPARED signatures. */
int ROM_std_strom_fom_prepare_A_callback(
    int window_id,
    void* context);

int ROM_std_strom_fom_apply_prepared_A_callback(
    const double* x,
    double* y,
    void* context);

/*
 * Build the nonhomogeneous Dirichlet lift U_D,w.
 *
 * Boundary entries contain the prescribed value at each slab/time DOF;
 * all unconstrained entries are zero.  The lifted ROM variable is
 *
 *     q_w = U_w - U_D,w,
 *
 * and therefore has homogeneous Dirichlet values.
 */
int ROM_std_strom_fom_build_dirichlet_lift(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    double* lift);

/*
 * Apply the matrix currently stored in context->sys->monolis without
 * rebuilding it.  Use this immediately after build_full_rhs() when the
 * matrix and RHS must come from exactly the same boundary-condition pass.
 */
int ROM_std_strom_fom_apply_prepared_A(
    ROM_STROM_FOM_CONTEXT* context,
    const double* x,
    double* y);

/* Callback y = A_w x for ROM_std_strom_calc_reduced_diag(). */
int ROM_std_strom_fom_apply_A(
    int window_id,
    const double* x,
    double* y,
    void* context);

/*
 * Callback y = B_w x_prev.
 *
 * B_w includes extraction of the previous window's last-slab temporal state
 * followed by the time-dG inflow assembly into the current window.
 */
int ROM_std_strom_fom_apply_B(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* context);

/*
 * Fast direct dG inflow action.  It assembles only the previous-window mass
 * contribution and zeroes Dirichlet rows, avoiding two complete RHS builds.
 */
int ROM_std_strom_fom_apply_B_direct(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* context);

/*
 * Assemble the full boundary-modified FOM right-hand side for one window.
 * previous_state has length num_time_dofs * num_space_nodes.
 */
int ROM_std_strom_fom_build_full_rhs(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_state,
    double* full_rhs);


/*
 * Squared norms used to aggregate residuals over multiple windows without
 * averaging already-normalized per-window values.
 */
typedef struct {
    double residual_norm2;
    double lhs_norm2;
    double rhs_norm2;
} ROM_STROM_RESIDUAL_NORMS;

/*
 * Residual diagnostics for the lifted representation
 *
 *     U_w = U_D,w + q_w.
 *
 * The boundary-modified full system is
 *
 *     A_tilde_w U_w = b_tilde_w(U_{w-1}),
 *
 * and the equivalent homogeneous system is
 *
 *     A_tilde_w q_w
 *       = b_tilde_w(U_{w-1}) - A_tilde_w U_D,w.
 */
typedef struct {
    ROM_STROM_RESIDUAL_NORMS all;
    ROM_STROM_RESIDUAL_NORMS free_dofs;
    ROM_STROM_RESIDUAL_NORMS dirichlet_dofs;

    /* ||r_full - r_hom||_2^2 */
    double full_hom_difference_norm2;

    /* Scale for the identity check, accumulated from both formulations. */
    double full_hom_scale_norm2;
} ROM_STROM_LIFTED_RESIDUAL;

/* Return ||r|| / (||lhs|| + ||rhs||). */
double ROM_std_strom_relative_residual(
    const ROM_STROM_RESIDUAL_NORMS* norms);

/* Add squared norms, for global accumulation over windows. */
void ROM_std_strom_accumulate_residual_norms(
    ROM_STROM_RESIDUAL_NORMS* destination,
    const ROM_STROM_RESIDUAL_NORMS* source);

/*
 * Evaluate both algebraically equivalent residuals using the SAME
 * Dirichlet-modified matrix produced by build_full_rhs():
 *
 *   r_hom  = A_tilde q - (b_tilde - A_tilde U_D),
 *   r_full = A_tilde (U_D + q) - b_tilde.
 *
 * previous_full_state is the full last-slab state of the previous window,
 * including its lift.  homogeneous_solution is q_w, not U_w.
 */
int ROM_std_strom_fom_calculate_lifted_residual(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_full_state,
    const double* homogeneous_solution,
    const double* current_lift,
    double* full_rhs,
    double* matrix_times_lift,
    double* homogeneous_rhs,
    double* matrix_times_homogeneous,
    double* homogeneous_residual,
    double* full_residual,
    ROM_STROM_LIFTED_RESIDUAL* result);

/*
 * Evaluate the full FOM residual
 *
 *     residual = A_w * window_solution - f_w(previous_state).
 *
 * The three optional squared norms are useful for a scale-safe relative
 * residual ||r|| / (||A_w U|| + ||f_w||).
 */
int ROM_std_strom_fom_calculate_residual(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_state,
    const double* window_solution,
    double* full_rhs,
    double* matrix_times_solution,
    double* residual,
    double* residual_norm2,
    double* matrix_times_solution_norm2,
    double* full_rhs_norm2);

int ROM_std_strom_fom_write_snapshot(
    const char* directory,
    int snapshot_id,
    int window_id,
    const double* values,
    int n);

int ROM_std_strom_fom_read_snapshot(
    const char* directory,
    int snapshot_id,
    int window_id,
    double* values,
    int n);

int ROM_std_strom_fom_write_basis(
    const char* directory,
    const ROM_STROM_WINDOW* window);

int ROM_std_strom_fom_read_basis(
    const char* directory,
    ROM_STROM_WINDOW* window);

/*
 * A basis built from the pooled snapshot matrix
 * [q_0^(0), ..., q_{Nw-1}^(0), q_0^(1), ...].
 */
int ROM_std_strom_fom_write_shared_basis(
    const char* directory,
    const ROM_STROM_WINDOW* basis_holder);

int ROM_std_strom_fom_read_shared_basis(
    const char* directory,
    ROM_STROM_WINDOW* window);

/* Persist/load projected D and L blocks for complete offline/online split. */
int ROM_std_strom_fom_write_reduced_operators(
    const char* directory,
    const ROM_STROM_SYSTEM* system,
    int operator_reuse,
    double dt,
    int st_dof);

int ROM_std_strom_fom_read_reduced_operators(
    const char* directory,
    ROM_STROM_SYSTEM* system,
    int operator_reuse,
    double dt,
    int st_dof);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STROM_FOM_ADAPTER_H */
