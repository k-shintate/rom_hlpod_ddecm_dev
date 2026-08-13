#pragma once

#include "fluid_sups_st_model.h"
#include "rom_std_stddrom.h"
#include "rom_std_stddrom_hlpod_bridge.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int converged;
    int nonlinear_iterations;

    double initial_full_residual;
    double final_full_residual;
    double relative_full_residual;

    double final_momentum_residual;
    double final_mass_residual;

    double initial_projected_residual;
    double final_projected_residual;
    double relative_projected_residual;
    double final_coefficient_correction;
} ROM_STDD_FLUID_WINDOW_REPORT;

typedef struct {
    FLUID_SUPS_SYSTEM* fom;
    const ST_FOM_WINDOW_LAYOUT* layout;
    ROM_STDD_SYSTEM* reduced;

    /* Communication-free matvec communicator after explicit halo exchange. */
    MONOLIS_COM mono_com0;

    int max_newton_iterations;
    double nonlinear_tolerance;
    double nonlinear_absolute_tolerance;

    int reduced_linear_max_iter;
    double reduced_linear_tolerance;

    double* initial_trace;
    double* previous_trace;

    double* current_window;
    double* lift;
    double* predictor;
    double* homogeneous;
    double* operator_input;
    double* reduced_correction;
} ROM_STDD_FLUID_CONTEXT;

int ROM_std_stdd_fluid_context_initialize(
    ROM_STDD_FLUID_CONTEXT* context,
    FLUID_SUPS_SYSTEM* fom,
    const ST_FOM_WINDOW_LAYOUT* layout,
    ROM_STDD_SYSTEM* reduced,
    int max_newton_iterations,
    double nonlinear_tolerance,
    double nonlinear_absolute_tolerance,
    int reduced_linear_max_iter,
    double reduced_linear_tolerance);

void ROM_std_stdd_fluid_context_finalize(
    ROM_STDD_FLUID_CONTEXT* context);

/*
 * Zero every physical velocity-Dirichlet and pressure-gauge row in the local
 * basis, then perform a two-pass local modified Gram--Schmidt.
 */
int ROM_std_stdd_fluid_enforce_homogeneous_basis(
    ROM_STDD_FLUID_CONTEXT* context);

/* HLPOD import followed by the exact homogeneous-basis enforcement above. */
int ROM_std_stdd_fluid_import_hlpod_basis(
    ROM_STDD_FLUID_CONTEXT* context,
    const HLPOD_MAT* hlpod_mat,
    int num_modes);

/*
 * Store one homogeneous owned-row ST window into an already allocated HLPOD
 * snapshot matrix:
 *
 *   snapmat[row][snapshot_id] = U_window[row] - U_D[row].
 */
int ROM_std_stdd_fluid_store_hlpod_snapshot(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* full_window_state,
    HLPOD_MAT* hlpod_mat,
    int snapshot_id,
    int snapshot_capacity);

/* Solve one nonlinear reduced window and update context->previous_trace. */
int ROM_std_stdd_fluid_solve_window(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    ROM_STDD_FLUID_WINDOW_REPORT* report);

/*
 * Write right-trace VTK files for every output-interval slab of the current
 * reconstructed window.
 */
int ROM_std_stdd_fluid_write_current_window_vtk(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    int* file_number);

/* Write the owned part of the current ROM window in the FOM-compatible format. */
int ROM_std_stdd_fluid_write_current_window_binary(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    const char* directory);

/* Relative owned-state error against a FOM window binary/reference array. */
int ROM_std_stdd_fluid_owned_window_error(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* reference_owned_window,
    double* relative_error);

#ifdef __cplusplus
}
#endif
