#ifndef ROM_STD_STDDROM_STSB_H
#define ROM_STD_STDDROM_STSB_H

#include "core_FOM_ST_spaceblock.h"
#include "rom_std_stddrom.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double residual_all2;
    double effective_rhs_all2;

    double residual_free2;
    double effective_rhs_free2;

    double residual_dirichlet2;
    double effective_rhs_dirichlet2;

    double homogeneous_dirichlet2;
    double homogeneous_all2;
} ROM_STDD_STSB_RESIDUAL_SUM;

typedef struct {
    FE_SYSTEM* sys;
    const ST_TIME_ELEMENT* time_element;

    int num_window_slabs;
    int num_windows;
    int st_dof_per_space_node;
    int local_st_dof;

    double* zero_previous_trace;
    double* previous_trace;
    double* previous_state_work;
    double* rhs_with_previous;
    double* rhs_without_previous;
    double* current_lift;
    double* previous_lift;
    double* operator_product;
} ROM_STDD_STSB_CONTEXT;

int ROM_std_stdd_stsb_context_initialize(
    ROM_STDD_STSB_CONTEXT* context,
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int num_windows);

void ROM_std_stdd_stsb_context_finalize(
    ROM_STDD_STSB_CONTEXT* context);

/*
 * Prepare the primitive, non-Dirichlet-eliminated space-time operator A.
 * Dirichlet data are represented only through lifting.
 */
int ROM_std_stdd_stsb_prepare_operator(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id);

/* Callback compatible with ROM_STDD_APPLY_OPERATOR. */
int ROM_std_stdd_stsb_apply_prepared_A(
    const double* x,
    double* y,
    void* user_context);

/* Window-invariant B action used during offline rank-sweep projection. */
int ROM_std_stdd_stsb_apply_B(
    const double* x_previous_window,
    double* y_current_window,
    void* user_context);

int ROM_std_stdd_stsb_build_dirichlet_lift(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    double* lift);

/*
 * Enforce exact homogeneous Dirichlet rows in the owned local basis and
 * re-orthonormalize the retained modes.  Call this after POD and before
 * writing/projecting the basis.
 */
int ROM_std_stdd_stsb_enforce_homogeneous_basis(
    ROM_STDD_STSB_CONTEXT* context,
    ROM_STDD_SYSTEM* system);

/* Build the primitive RHS b_w(previous_trace), without BC elimination. */
int ROM_std_stdd_stsb_build_raw_rhs(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* previous_trace,
    double* raw_rhs);

/*
 * Build
 *
 *   b_h,w = b_w(previous full/lift trace) - A U_D,w.
 *
 * For w=0, initial_trace is used.  For w>0, the previous-window lift trace
 * is used; the previous homogeneous state is represented by the reduced
 * temporal block in the all-at-once matrix.
 */
int ROM_std_stdd_stsb_build_effective_rhs(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* initial_trace,
    double* effective_rhs);

/*
 * Evaluate the primitive full-state residual
 *
 *   r_w = A U_w - b_w(U_{w-1}^+),
 *
 * split into owned free and Dirichlet spatial rows.
 */
int ROM_std_stdd_stsb_evaluate_raw_full_state_residual(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* previous_trace,
    const double* full_window_state,
    ROM_STDD_STSB_RESIDUAL_SUM* residual_sum);

/*
 * Evaluate the lifted FOM residual for a stored homogeneous reference state:
 *
 *   r_w = A q_w
 *       - [b_w(U_{w-1}^{full}) - A U_{D,w}],
 *
 * where A and b are primitive and have not undergone Dirichlet elimination.
 *
 * For window 0, initial_trace is used.  For later windows,
 * previous_reference_internal is communicated, the previous lift is added,
 * and its complete right trace is used.  direct_projected_defect receives
 * Phi^T r_w on the local reduced row block.
 */
int ROM_std_stdd_stsb_evaluate_lifted_reference_residual(
    ROM_STDD_STSB_CONTEXT* context,
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* initial_trace,
    const double* current_reference_internal,
    const double* previous_reference_internal,
    double* direct_projected_defect,
    ROM_STDD_STSB_RESIDUAL_SUM* residual_sum);

/* Copy owned rows of q = U-U_D into snapshot storage. */
int ROM_std_stdd_stsb_store_internal_homogeneous(
    const ROM_STDD_SYSTEM* system,
    const double* full_window_state,
    const double* lift,
    double* internal_snapshot);

/*
 * Reconstruct a full local window state.  The ROM homogeneous vector is
 * communicated over the ordinary spatial partition before the lift is added.
 */
int ROM_std_stdd_stsb_reconstruct_full_window(
    ROM_STDD_STSB_CONTEXT* context,
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* homogeneous_work,
    double* full_window_state);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STDDROM_STSB_H */
