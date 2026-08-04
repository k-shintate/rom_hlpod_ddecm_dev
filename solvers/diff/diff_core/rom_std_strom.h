#ifndef ROM_STD_STROM_H
#define ROM_STD_STROM_H

#include <stddef.h>

/*
 * Replace this include with the header that defines HLPOD_VALUES and HLPOD_MAT
 * in your project.
 */
#include "core_ROM.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    ROM_STROM_SUCCESS = 0,
    ROM_STROM_ERR_ARGUMENT = 1,
    ROM_STROM_ERR_ALLOCATION = 2,
    ROM_STROM_ERR_OPERATOR = 3,
    ROM_STROM_ERR_SOLVER = 4,
    ROM_STROM_ERR_IO = 5,
    ROM_STROM_ERR_DIMENSION = 6
} ROM_STROM_STATUS;

/* y = A_w x */
typedef int (*ROM_STROM_APPLY_WINDOW_OPERATOR)(
    int window_id,
    const double* x,
    double* y,
    void* context);

/* Prepare A_w once, then reuse it for every POD mode in that window. */
typedef int (*ROM_STROM_PREPARE_WINDOW_OPERATOR)(
    int window_id,
    void* context);

typedef int (*ROM_STROM_APPLY_PREPARED_OPERATOR)(
    const double* x,
    double* y,
    void* context);

/* y_cur = B_w x_prev = G_w C_{+,w-1} x_prev */
typedef int (*ROM_STROM_APPLY_WINDOW_COUPLING)(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* context);

/* y = M_x x, or another spatial operator used in the dG trace term. */
typedef int (*ROM_STROM_APPLY_SPATIAL_OPERATOR)(
    const double* x,
    double* y,
    void* context);

/* Solve A x = b. The implementation must not assume that A is symmetric. */
typedef int (*ROM_STROM_DENSE_SOLVER)(
    int n,
    double** A,
    const double* b,
    double* x,
    void* context);

/* Callbacks for assembling an all-at-once sparse reduced matrix. */
typedef int (*ROM_STROM_ADD_MATRIX_ENTRY)(
    int row,
    int column,
    double value,
    void* context);

typedef int (*ROM_STROM_SET_VECTOR_ENTRY)(
    int row,
    double value,
    void* context);

typedef struct {
    int window_id;

    int total_num_nodes;
    int dof;
    int num_slabs;
    int num_time_basis;

    int spatial_dof;
    int st_dof;

    int num_snapshots;
    int num_modes_max;
    int num_modes;

    HLPOD_VALUES* hlpod_vals;
    HLPOD_MAT* hlpod_mat;

    /* Singular values retained from the offline SVD. */
    double* singular_values;
    int num_singular_values;

    /*
     * Shared all-window bases may be attached without a deep copy.  These
     * flags prevent double-free while preserving the legacy owning path.
     */
    int owns_pod_modes;
    int owns_singular_values;

    /* D_w = Phi_w^T A_w Phi_w. */
    double** reduced_diag;

    /* L_w = Phi_w^T B_w Phi_{w-1}. Null for w = 0. */
    double** reduced_coupling_prev;

    /* g_w and a_w. */
    double* reduced_rhs;
    double* reduced_coef;

    /* Optional reference state U_ref,w. */
    double* reference_state;
} ROM_STROM_WINDOW;

typedef struct {
    int num_windows;
    int total_reduced_dof;

    ROM_STROM_WINDOW* windows;
    int* mode_offset;
} ROM_STROM_SYSTEM;

/*
 * Data for the reference implementation of
 * B_w = G_w C_{+,w-1}.
 *
 * The coefficient arrays contain the evaluation of the temporal trial basis
 * at the previous-window right endpoint and the current-window left-boundary
 * dG test/jump coefficients, respectively.
 */
typedef struct {
    int spatial_dof;
    int prev_num_slabs;
    int prev_num_time_basis;
    int cur_num_slabs;
    int cur_num_time_basis;

    const double* theta_plus_prev;
    const double* jump_left_cur;

    ROM_STROM_APPLY_SPATIAL_OPERATOR apply_spatial_trace_operator;
    void* spatial_operator_context;
} ROM_STROM_DG_TRANSITION;

typedef struct {
    int num_windows;
    ROM_STROM_DG_TRANSITION* transition; /* transition[w], w >= 1 */
} ROM_STROM_DG_COUPLING_CONTEXT;

/* ---------- Basic indexing and memory lifecycle ---------- */

int ROM_std_strom_index(
    int node,
    int component,
    int slab,
    int time_basis,
    int total_num_nodes,
    int dof,
    int num_time_basis);

int ROM_std_strom_window_initialize(
    ROM_STROM_WINDOW* window,
    int window_id,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT* hlpod_mat,
    int total_num_nodes,
    int dof,
    int num_slabs,
    int num_time_basis,
    int num_snapshots,
    int num_modes_max,
    int allocate_reference_state);

void ROM_std_strom_window_finalize(ROM_STROM_WINDOW* window);

int ROM_std_strom_system_initialize(
    ROM_STROM_SYSTEM* system,
    int num_windows);

void ROM_std_strom_system_finalize(ROM_STROM_SYSTEM* system);

int ROM_std_strom_update_mode_offsets(ROM_STROM_SYSTEM* system);

/* ---------- Offline snapshot and POD operations ---------- */

int ROM_std_strom_store_snapshot(
    ROM_STROM_WINDOW* window,
    const double* window_solution,
    int snapshot_id,
    int subtract_reference_state);

int ROM_std_strom_set_pod_modes(
    ROM_STROM_WINDOW* window,
    double rom_epsilon);

int ROM_std_strom_write_pod_modes_binary(
    const ROM_STROM_WINDOW* window,
    const char* filename);

int ROM_std_strom_read_pod_modes_binary(
    ROM_STROM_WINDOW* window,
    const char* filename);

/*
 * Deep-copy an already constructed POD basis into another window object.
 * The destination keeps its own snapshot matrix and owns an independent
 * copy of the basis and singular values.
 */
int ROM_std_strom_copy_pod_basis(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source);

/*
 * Attach the source basis by reference.  Source must outlive destination.
 * This avoids num_windows deep copies of the all-window shared basis.
 */
int ROM_std_strom_attach_pod_basis(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source);

/* ---------- Projection of window and coupling operators ---------- */

int ROM_std_strom_calc_reduced_diag(
    ROM_STROM_WINDOW* window,
    ROM_STROM_APPLY_WINDOW_OPERATOR apply_A,
    void* operator_context);

/*
 * Fast path: prepare A_w once, form Y=A_w Phi by sparse matvecs, then use
 * BLAS DGEMM for D=Phi^T Y when ROM_STROM_USE_CBLAS is enabled.
 */
int ROM_std_strom_calc_reduced_diag_prepared(
    ROM_STROM_WINDOW* window,
    ROM_STROM_PREPARE_WINDOW_OPERATOR prepare_A,
    ROM_STROM_APPLY_PREPARED_OPERATOR apply_prepared_A,
    void* operator_context);

int ROM_std_strom_calc_reduced_coupling(
    ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context);

/* Form Y=B_w Phi_prev, then use BLAS DGEMM for L=Phi_cur^T Y. */
int ROM_std_strom_calc_reduced_coupling_fast(
    ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context);

int ROM_std_strom_copy_reduced_diag(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source);

int ROM_std_strom_copy_reduced_coupling(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source);

int ROM_std_strom_set_reduced_diag_data(
    ROM_STROM_WINDOW* window,
    const double* row_major_data);

int ROM_std_strom_set_reduced_coupling_data(
    ROM_STROM_WINDOW* window,
    const double* row_major_data,
    int previous_num_modes);

int ROM_std_strom_project_rhs(
    ROM_STROM_WINDOW* window,
    const double* effective_full_rhs);

int ROM_std_strom_build_effective_rhs(
    const ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    const double* full_rhs,
    double* effective_rhs,
    ROM_STROM_APPLY_WINDOW_OPERATOR apply_A,
    void* operator_context,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context);

/* ---------- Online solves ---------- */

int ROM_std_strom_solve_dense_gauss(
    int n,
    double** A,
    const double* b,
    double* x,
    void* context);

/* LAPACK DGESV path, enabled with ROM_STROM_USE_LAPACK. */
int ROM_std_strom_solve_dense_lapack(
    int n,
    double** A,
    const double* b,
    double* x,
    void* context);

int ROM_std_strom_solve_sequential(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_DENSE_SOLVER dense_solver,
    void* solver_context);

int ROM_std_strom_assemble_all_at_once(
    const ROM_STROM_SYSTEM* system,
    ROM_STROM_ADD_MATRIX_ENTRY add_matrix_entry,
    ROM_STROM_SET_VECTOR_ENTRY set_rhs_entry,
    void* assembly_context);

/* Dense all-at-once solver intended for verification and small ROMs. */
int ROM_std_strom_solve_all_at_once_dense(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_DENSE_SOLVER dense_solver,
    void* solver_context);

/* ---------- Reconstruction and trace operations ---------- */

int ROM_std_strom_reconstruct_window(
    const ROM_STROM_WINDOW* window,
    double* window_solution);


/* Euclidean orthogonal projection onto U_ref,w + span(Phi_w). */
int ROM_std_strom_project_window_state(
    const ROM_STROM_WINDOW* window,
    const double* window_state,
    double* projected_state);

int ROM_std_strom_extract_right_trace(
    const double* window_vector,
    double* trace,
    int spatial_dof,
    int num_slabs,
    int num_time_basis,
    const double* theta_plus);

int ROM_std_strom_inject_left_trace(
    const double* spatial_trace,
    double* current_window_vector,
    int spatial_dof,
    int num_time_basis,
    const double* jump_left_coefficients);

/* Reference implementation of y_cur = G_w C_{+,w-1} x_prev. */
int ROM_std_strom_apply_dg_coupling(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* context);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STROM_H */
