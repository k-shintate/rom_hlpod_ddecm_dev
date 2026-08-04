/*
 * Minimal integration skeleton for the ST-ROM implementation.
 * Replace project-specific include names and assembly routines as necessary.
 */
#include "rom_std_strom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "monolis.h"
#include "monolis_utils.h"

typedef struct {
    MONOLIS** window_matrix;
    MONOLIS_COM** window_com;
} EXAMPLE_WINDOW_OPERATOR_CONTEXT;

typedef struct {
    MONOLIS* spatial_mass_matrix;
    MONOLIS_COM* spatial_mass_com;
} EXAMPLE_SPATIAL_MASS_CONTEXT;

static int example_apply_window_operator(
    int window_id,
    const double* x,
    double* y,
    void* context)
{
    EXAMPLE_WINDOW_OPERATOR_CONTEXT* ctx =
        (EXAMPLE_WINDOW_OPERATOR_CONTEXT*)context;

    if(ctx == NULL || ctx->window_matrix == NULL ||
       ctx->window_com == NULL){
        return 1;
    }

    monolis_matvec_product_R(
        ctx->window_matrix[window_id],
        ctx->window_com[window_id],
        (double*)x,
        y);

    return 0;
}

static int example_apply_spatial_mass(
    const double* x,
    double* y,
    void* context)
{
    EXAMPLE_SPATIAL_MASS_CONTEXT* ctx =
        (EXAMPLE_SPATIAL_MASS_CONTEXT*)context;

    if(ctx == NULL || ctx->spatial_mass_matrix == NULL ||
       ctx->spatial_mass_com == NULL){
        return 1;
    }

    monolis_matvec_product_R(
        ctx->spatial_mass_matrix,
        ctx->spatial_mass_com,
        (double*)x,
        y);

    return 0;
}

/*
 * Offline construction.
 *
 * training_window_solution[w][m] points to an ST vector of length st_dof(w).
 * full_window_rhs[w] points to the FOM window RHS f_w.
 */
int example_build_strom(
    ROM_STROM_SYSTEM* system,
    HLPOD_VALUES* pod_values,
    HLPOD_MAT* pod_matrices,
    int num_windows,
    int total_num_nodes,
    int dof,
    const int* num_slabs,
    const int* num_time_basis,
    int num_snapshots,
    int num_modes_max,
    double rom_epsilon,
    double*** training_window_solution,
    const double** full_window_rhs,
    EXAMPLE_WINDOW_OPERATOR_CONTEXT* A_context,
    ROM_STROM_DG_COUPLING_CONTEXT* B_context)
{
    int status;

    status = ROM_std_strom_system_initialize(system, num_windows);
    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    for(int w = 0; w < num_windows; w++){
        status = ROM_std_strom_window_initialize(
            &system->windows[w],
            w,
            &pod_values[w],
            &pod_matrices[w],
            total_num_nodes,
            dof,
            num_slabs[w],
            num_time_basis[w],
            num_snapshots,
            num_modes_max,
            0); /* set to 1 when U_ref,w is used */

        if(status != ROM_STROM_SUCCESS){
            return status;
        }

        for(int m = 0; m < num_snapshots; m++){
            status = ROM_std_strom_store_snapshot(
                &system->windows[w],
                training_window_solution[w][m],
                m,
                0);

            if(status != ROM_STROM_SUCCESS){
                return status;
            }
        }

        status = ROM_std_strom_set_pod_modes(
            &system->windows[w],
            rom_epsilon);

        if(status != ROM_STROM_SUCCESS){
            return status;
        }

        status = ROM_std_strom_calc_reduced_diag(
            &system->windows[w],
            example_apply_window_operator,
            A_context);

        if(status != ROM_STROM_SUCCESS){
            return status;
        }
    }

    for(int w = 1; w < num_windows; w++){
        status = ROM_std_strom_calc_reduced_coupling(
            &system->windows[w],
            &system->windows[w - 1],
            ROM_std_strom_apply_dg_coupling,
            B_context);

        if(status != ROM_STROM_SUCCESS){
            return status;
        }
    }

    for(int w = 0; w < num_windows; w++){
        double* effective_rhs = (double*)calloc(
            (size_t)system->windows[w].st_dof,
            sizeof(double));

        if(effective_rhs == NULL){
            return ROM_STROM_ERR_ALLOCATION;
        }

        status = ROM_std_strom_build_effective_rhs(
            &system->windows[w],
            w > 0 ? &system->windows[w - 1] : NULL,
            full_window_rhs[w],
            effective_rhs,
            example_apply_window_operator,
            A_context,
            ROM_std_strom_apply_dg_coupling,
            B_context);

        if(status == ROM_STROM_SUCCESS){
            status = ROM_std_strom_project_rhs(
                &system->windows[w],
                effective_rhs);
        }

        free(effective_rhs);

        if(status != ROM_STROM_SUCCESS){
            return status;
        }
    }

    return ROM_std_strom_update_mode_offsets(system);
}

/* Online sequential solve and reconstruction. */
int example_run_strom(
    ROM_STROM_SYSTEM* system,
    double** reconstructed_window_solution)
{
    int status = ROM_std_strom_solve_sequential(
        system,
        NULL, /* use the included pivoted Gaussian solver */
        NULL);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    for(int w = 0; w < system->num_windows; w++){
        status = ROM_std_strom_reconstruct_window(
            &system->windows[w],
            reconstructed_window_solution[w]);

        if(status != ROM_STROM_SUCCESS){
            return status;
        }
    }

    return ROM_STROM_SUCCESS;
}

/*
 * Verification procedure:
 * 1. Save coefficients from example_run_strom().
 * 2. Call ROM_std_strom_solve_all_at_once_dense().
 * 3. Compare both coefficient vectors. They should agree to roundoff.
 */
