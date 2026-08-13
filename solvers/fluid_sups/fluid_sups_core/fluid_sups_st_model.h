#pragma once

#include "fluid_core.h"
#include "st_fom_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FLUID_SUPS_PHYSICAL_DOF 4

typedef struct {
    int num_ip_each_axis;
    double mat_epsilon;
    int mat_max_iter;

    double dt;
    double finish_time;
    int output_interval;

    double density;
    double viscosity;

    double nr_epsilon;
    int nr_max_iter;

    int st_time_degree;
    int st_slabs_per_window;
    int write_window_binary;

    /*
     * Finite-difference scale used only to assemble the previous-slab
     * Jacobian block dR_s/dU_{s-1}.  The nonlinear residual and diagonal
     * Newton Jacobian remain the verified SUPS/NR kernels.
     */
    double st_previous_jacobian_fd_epsilon;

    /*
     * Finite-difference scale used for the previous-state directional
     * derivative needed by the dG(1)/dG(2) SUPS temporal Jacobian.
     */
    double st_dg_jacobian_fd_epsilon;

    double** v;
    double* p;
    double** delta_v;
    double* delta_p;
    double** v_old;
    double* p_old;
} FLUID_SUPS_VALUES;

typedef struct {
    const char* directory;
} FLUID_SUPS_CONDITIONS;

typedef struct {
    BBFE_BASIS basis;
    BBFE_DATA fe;
    FLUID_SUPS_CONDITIONS cond;
    FLUID_SUPS_VALUES vals;
    BBFE_BC bc;
    BBFE_BC bc_NR;

    /* Verified one-step spatial NR system (4 dof/node). */
    MONOLIS monolis;
    MONOLIS_COM mono_com;

    /*
     * True all-at-once dG window Newton system.
     *
     * block size = 4 * num_slabs * (time_degree + 1) dof/node.
     *
     * Examples for four slabs:
     *   dG(0): 16 dof/node
     *   dG(1): 32 dof/node
     *   dG(2): 48 dof/node
     */
    BBFE_BC bc_NR_window;
    MONOLIS monolis_window;

    /*
     * Same design as the verified convection-diffusion STSB FOM:
     *
     *   - one physical spatial MPI decomposition,
     *   - one enlarged space-time block DOF on each spatial node,
     *   - the original spatial MONOLIS communicator is reused for the
     *     window solve and all halo exchanges.
     */
    ST_FOM_SPATIAL_GRAPH st_spatial_graph;
    int num_owned_space_nodes;
    int window_dof;
} FLUID_SUPS_SYSTEM;

typedef struct {
    int converged;
    int nonlinear_iterations;
    double relative_total_residual;
    double relative_momentum_residual;
    double relative_mass_residual;
} FLUID_SUPS_NR_REPORT;

int fluid_sups_st_initialize(
    FLUID_SUPS_SYSTEM* system,
    int argc,
    char* argv[]);

void fluid_sups_st_finalize(FLUID_SUPS_SYSTEM* system);

int fluid_sups_st_solve_step(
    FLUID_SUPS_SYSTEM* system,
    int step,
    double time_old,
    double time_new,
    FLUID_SUPS_NR_REPORT* report);

/*
 * Solve one complete dG(0) window as ONE nonlinear Newton system.
 *
 * Each Newton iteration assembles the full lower-bidiagonal space-time
 * Jacobian with 4*num_slabs dof/node and calls MONOLIS exactly once:
 *
 *   [J0             ] [dU0]   [-R0]
 *   [C1 J1          ] [dU1] = [-R1]
 *   [   C2 J2       ] [dU2]   [-R2]
 *   [      ...  ... ] [...]   [...]
 *
 * J_s is the verified current-state NR Jacobian.
 * C_s = dR_s/dU_{s-1} is assembled elementwise by finite differences of the
 * verified residual with respect to previous-slab velocity.  Previous pressure
 * has no coupling because the verified residual kernel does not take p_old.
 */
int fluid_sups_st_solve_window_all_at_once_dg0(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report);

/*
 * Generic all-at-once dG window solve.
 *
 * degree = 0 dispatches to the already validated dG(0) implementation.
 * degree = 1 or 2 uses the temporal weak residual:
 *
 *   int_0^1 psi_alpha R_SUPS(U, dU/dt) d tau
 *   + e_minus_alpha * J_jump(U^- - U_prev^+) = 0.
 *
 * The point residual is the verified SUPS/NR residual kernel.  At a temporal
 * quadrature point,
 *
 *   U(tau)       = sum_beta psi_beta(tau) U_beta,
 *   Delta U_time = sum_beta dpsi_beta/dtau U_beta,
 *   U_old,eff    = U(tau) - Delta U_time,
 *
 * so the existing kernel receives exactly the physical dU/dt represented by
 * the higher-order temporal polynomial.
 *
 * The jump operator is obtained from the same verified kernel through its
 * previous-state directional derivative, rather than assuming a simplified
 * -M/dt coupling.  This keeps the SUPS/PSPG transient contribution.
 */
int fluid_sups_st_solve_window_all_at_once_dg(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report);


/* ------------------------------------------------------------------------- */
/* Non-hyper-reduced ST-DDROM adapter API                                    */
/* ------------------------------------------------------------------------- */

/*
 * Build the affine window lift used by the ROM:
 *
 *   U = U_D + Phi a.
 *
 * U_D contains the physical velocity Dirichlet values and the pressure-gauge
 * values, and is synchronized over the spatial MONOLIS communicator.
 */
int fluid_sups_st_rom_build_window_lift(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* lift);

/*
 * Build the constant-in-time window predictor from the previous right trace.
 * This uses the exact same node-major ST ordering as the all-at-once FOM.
 */
int fluid_sups_st_rom_build_window_predictor(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* predictor);

/*
 * Enforce physical Dirichlet values, the pressure gauge, and halo consistency
 * on a reconstructed full window state.
 */
int fluid_sups_st_rom_enforce_window_constraints(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* window_state);

/*
 * Assemble the exact full-order Newton pair at the supplied physical state.
 *
 * On return:
 *   monolis_window.mat.R.A = J(U)
 *   monolis_window.mat.R.B = -R(U)
 *
 * including the same homogeneous correction BC and pressure gauge as the
 * validated all-at-once FOM.
 */
int fluid_sups_st_rom_assemble_window_newton(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state);

/* Assemble only -R(U), for post-update convergence evaluation. */
int fluid_sups_st_rom_assemble_window_residual(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state);

/* Global norms of the currently assembled monolis_window.mat.R.B. */
int fluid_sups_st_rom_get_window_rhs_norms(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* total_norm,
    double* momentum_norm,
    double* mass_norm);

int fluid_sups_st_export_state(
    const FLUID_SUPS_SYSTEM* system,
    double* state);

int fluid_sups_st_import_state(
    FLUID_SUPS_SYSTEM* system,
    const double* state,
    int update_previous_state);

int fluid_sups_st_synchronize_state(FLUID_SUPS_SYSTEM* system);

int fluid_sups_st_write_output(
    FLUID_SUPS_SYSTEM* system,
    int file_number,
    double time);

/*
 * Same output payload as fluid_sups_st_write_output(), with a caller-selected
 * filename prefix.  The actual rank suffix remains ".vtk.<rank>".
 */
int fluid_sups_st_write_output_named(
    FLUID_SUPS_SYSTEM* system,
    const char* prefix,
    int file_number,
    double time);

void fluid_sups_st_print_conditions(const FLUID_SUPS_SYSTEM* system);

#ifdef __cplusplus
}
#endif
