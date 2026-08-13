#ifndef ROM_STD_STDDROM_H
#define ROM_STD_STDDROM_H

/*
 * Space-time domain-decomposed reduced-order model.
 *
 * Version 1 deliberately uses one spatial ROM subdomain per MPI rank and one
 * shared local POD basis per rank for every time window.  The global reduced
 * metanode is
 *
 *     (window, spatial-rank).
 *
 * The reduced graph contains
 *
 *   - same-window spatial-neighbour blocks, and
 *   - previous-window temporal blocks for self and spatial neighbours.
 *
 * This is the distributed counterpart of the serial windowed ST-ROM
 * all-at-once system.  A fixed common number of modes is used on all ranks so
 * that MONOLIS BCSR can be used safely.  Local POD ranks estimated from an
 * energy tolerance are promoted to the global maximum rank.
 */

#include <stddef.h>
#include <stdint.h>

#include "core_ROM.h"
#include "core_FOM_ST.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    ROM_STDD_SUCCESS = 0,
    ROM_STDD_ERR_ARGUMENT = 1,
    ROM_STDD_ERR_ALLOCATION = 2,
    ROM_STDD_ERR_IO = 3,
    ROM_STDD_ERR_DIMENSION = 4,
    ROM_STDD_ERR_OPERATOR = 5,
    ROM_STDD_ERR_SOLVER = 6,
    ROM_STDD_ERR_INCOMPATIBLE = 7
} ROM_STDD_STATUS;

/* y = A x.  The callback is collective over the spatial MPI communicator. */
typedef int (*ROM_STDD_APPLY_OPERATOR)(
    const double* x,
    double* y,
    void* context);

/* y = B x_prev.  This callback is also collective. */
typedef int (*ROM_STDD_APPLY_COUPLING)(
    const double* x_prev,
    double* y,
    void* context);

typedef struct {
    int num_windows;
    int num_modes;
    int num_modes_capacity;

    int num_space_nodes;
    int num_internal_space_nodes;
    int st_dof_per_space_node;
    int local_st_dof;
    int internal_st_dof;

    int num_neighbor_ranks;
    int* neighbor_ranks; /* receive-neighbour order from the spatial COM. */

    /*
     * Shared local basis on owned spatial nodes only.
     * Row-major layout:
     *
     *   basis[internal_st_row * num_modes_capacity + mode]
     */
    double* basis;
    double* singular_values;
    int num_singular_values;

    /*
     * Reduced operator blocks.  slot 0 is self; slot n+1 is neighbour n.
     * Each block is row-major num_modes x num_modes.
     * The operators are assumed window invariant in version 1.
     */
    double* diagonal_blocks;
    double* temporal_blocks;

    /* Internal reduced right-hand side and solution, window-major. */
    double* reduced_rhs;
    double* reduced_coef;
} ROM_STDD_SYSTEM;

typedef struct {
    MONOLIS monolis;
    MONOLIS_COM comm;

    int num_windows;
    int block_size;
    int num_internal_meta_nodes;
    int num_total_meta_nodes;
    int num_neighbor_ranks;

    int* graph_index;
    int* graph_item;
    int num_graph_items;

    double* solution;
    int max_iter;
    double tolerance;
} ROM_STDD_REDUCED_SOLVER;

/* ---------------------------- lifecycle -------------------------------- */

int ROM_std_stdd_system_initialize(
    ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int num_windows,
    int num_space_nodes,
    int num_internal_space_nodes,
    int st_dof_per_space_node,
    int num_modes_capacity);

void ROM_std_stdd_system_finalize(
    ROM_STDD_SYSTEM* system);

/* -------------------------- snapshot I/O ------------------------------- */

int ROM_std_stdd_write_snapshot_rank(
    const ROM_STDD_SYSTEM* system,
    const double* window_snapshots,
    int snapshot_id,
    const char* directory);

int ROM_std_stdd_read_snapshot_rank(
    const ROM_STDD_SYSTEM* system,
    double* window_snapshots,
    int snapshot_id,
    const char* directory);

/* ------------------------------ POD ------------------------------------ */

int ROM_std_stdd_compute_shared_local_pod(
    ROM_STDD_SYSTEM* system,
    const double* snapshot_matrix,
    int num_snapshot_columns,
    double energy_tolerance,
    int communicator);

int ROM_std_stdd_write_basis_rank(
    const ROM_STDD_SYSTEM* system,
    int num_snapshot_columns,
    const char* directory);

int ROM_std_stdd_read_basis_rank(
    ROM_STDD_SYSTEM* system,
    const char* directory);

/* ----------------------- distributed projection ------------------------ */

int ROM_std_stdd_project_operators_rank_sweep(
    ROM_STDD_SYSTEM* system,
    ROM_STDD_APPLY_OPERATOR apply_A,
    ROM_STDD_APPLY_COUPLING apply_B,
    void* operator_context);

int ROM_std_stdd_write_operators_rank(
    const ROM_STDD_SYSTEM* system,
    const char* directory);

int ROM_std_stdd_read_operators_rank(
    ROM_STDD_SYSTEM* system,
    const char* directory);

/* g_local = Phi_local^T b_local over owned ST rows. */
int ROM_std_stdd_project_local_rhs(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* effective_full_rhs,
    double* reduced_rhs_out);

/* q_local = Phi_local a on owned rows; halo rows are left zero. */
int ROM_std_stdd_reconstruct_local_homogeneous(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* local_homogeneous_state);

/* ---------------------- reduced graph and solve ------------------------ */

int ROM_std_stdd_reduced_solver_initialize(
    ROM_STDD_REDUCED_SOLVER* solver,
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int max_iter,
    double tolerance);

int ROM_std_stdd_reduced_solver_assemble(
    ROM_STDD_REDUCED_SOLVER* solver,
    const ROM_STDD_SYSTEM* system);

int ROM_std_stdd_reduced_solver_solve(
    ROM_STDD_REDUCED_SOLVER* solver,
    ROM_STDD_SYSTEM* system);

/*
 * Solve the block-lower-triangular all-window system by exact window
 * marching.  Each window solve is still a distributed spatial reduced
 * solve over all MPI ranks.  The local dense self block is used as a
 * left block-Jacobi scaling before MONOLIS BiCGSTAB.
 */
int ROM_std_stdd_solve_window_marching(
    ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int max_iter,
    double tolerance);

/*
 * Relative residual of the distributed reduced equations
 *
 *   D a_w - L a_{w-1} = g_w
 *
 * evaluated directly from the stored reduced blocks.
 */
double ROM_std_stdd_system_reduced_residual(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com);

/*
 * Evaluate the same distributed reduced equations for an arbitrary
 * window-major coefficient trajectory.  The coefficient layout is
 *
 *   coefficients[window * num_modes_capacity + mode].
 *
 * When window_relative_residuals is non-NULL, one relative residual is
 * returned for each window.  The return value is the global all-window
 * relative residual.
 */
double ROM_std_stdd_system_reduced_residual_from_coefficients(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    const double* coefficients,
    double* window_relative_residuals);

/*
 * Evaluate one local reduced row block for an arbitrary coefficient
 * trajectory.  Every output vector has system->num_modes entries.
 *
 *   defect = self_action + neighbor_action
 *          - temporal_action - reduced_rhs.
 *
 * temporal_action is zero for window 0.
 */
int ROM_std_stdd_evaluate_local_reduced_equation(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    const double* coefficients,
    int window_id,
    double* defect,
    double* self_action,
    double* neighbor_action,
    double* temporal_action);

int ROM_std_stdd_write_reduced_graph_rank(
    const ROM_STDD_REDUCED_SOLVER* solver,
    const char* directory);

void ROM_std_stdd_reduced_solver_finalize(
    ROM_STDD_REDUCED_SOLVER* solver);

/* Relative ||A_r a - g||_2 / ||g||_2 over all ranks. */
double ROM_std_stdd_reduced_algebraic_residual(
    ROM_STDD_REDUCED_SOLVER* solver);

/* ---------------------------- validation ------------------------------- */

double ROM_std_stdd_global_relative_error(
    const double* value,
    const double* reference,
    int local_length,
    int communicator);

int ROM_std_stdd_project_reference_coefficients(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* reference_homogeneous,
    double* coefficients);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STDDROM_H */
