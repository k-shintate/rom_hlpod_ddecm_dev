#include "rom_std_stddrom_fluid.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t fluid_rom_basis_offset(
    const ROM_STDD_SYSTEM* system,
    int row,
    int mode)
{
    return
        (size_t)row * (size_t)system->num_modes_capacity
        + (size_t)mode;
}

static int fluid_rom_validate(
    const ROM_STDD_FLUID_CONTEXT* context)
{
    if(context == NULL || context->fom == NULL ||
       context->layout == NULL || context->reduced == NULL ||
       context->reduced->basis == NULL ||
       context->reduced->num_modes <= 0 ||
       context->reduced->num_modes >
           context->reduced->num_modes_capacity ||
       context->reduced->num_space_nodes !=
           context->fom->fe.total_num_nodes ||
       context->reduced->num_internal_space_nodes !=
           context->fom->num_owned_space_nodes ||
       context->reduced->st_dof_per_space_node !=
           context->fom->window_dof ||
       context->reduced->local_st_dof !=
           (int)ST_fom_window_vector_length(context->layout) ||
       context->reduced->internal_st_dof !=
           context->fom->num_owned_space_nodes
               * context->fom->window_dof)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    return ROM_STDD_SUCCESS;
}

static void fluid_rom_zero_report(
    ROM_STDD_FLUID_WINDOW_REPORT* report)
{
    if(report == NULL){
        return;
    }

    memset(report, 0, sizeof(*report));
    report->initial_full_residual = INFINITY;
    report->final_full_residual = INFINITY;
    report->relative_full_residual = INFINITY;
    report->final_momentum_residual = INFINITY;
    report->final_mass_residual = INFINITY;
    report->initial_projected_residual = INFINITY;
    report->final_projected_residual = INFINITY;
    report->relative_projected_residual = INFINITY;
    report->final_coefficient_correction = INFINITY;
}

int ROM_std_stdd_fluid_context_initialize(
    ROM_STDD_FLUID_CONTEXT* context,
    FLUID_SUPS_SYSTEM* fom,
    const ST_FOM_WINDOW_LAYOUT* layout,
    ROM_STDD_SYSTEM* reduced,
    int max_newton_iterations,
    double nonlinear_tolerance,
    double nonlinear_absolute_tolerance,
    int reduced_linear_max_iter,
    double reduced_linear_tolerance)
{
    const size_t local_st_dof =
        reduced != NULL ? (size_t)reduced->local_st_dof : 0u;
    const size_t trace_dof =
        fom != NULL
        ? (size_t)fom->fe.total_num_nodes
            * (size_t)FLUID_SUPS_PHYSICAL_DOF
        : 0u;

    if(context == NULL || fom == NULL || layout == NULL ||
       reduced == NULL || max_newton_iterations <= 0 ||
       !isfinite(nonlinear_tolerance) ||
       nonlinear_tolerance <= 0.0 ||
       !isfinite(nonlinear_absolute_tolerance) ||
       nonlinear_absolute_tolerance <= 0.0 ||
       reduced_linear_max_iter <= 0 ||
       !isfinite(reduced_linear_tolerance) ||
       reduced_linear_tolerance <= 0.0 ||
       reduced->num_space_nodes != fom->fe.total_num_nodes ||
       reduced->num_internal_space_nodes !=
           fom->num_owned_space_nodes ||
       reduced->st_dof_per_space_node != fom->window_dof ||
       local_st_dof != ST_fom_window_vector_length(layout))
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(context, 0, sizeof(*context));

    context->fom = fom;
    context->layout = layout;
    context->reduced = reduced;
    context->max_newton_iterations = max_newton_iterations;
    context->nonlinear_tolerance = nonlinear_tolerance;
    context->nonlinear_absolute_tolerance =
        nonlinear_absolute_tolerance;
    context->reduced_linear_max_iter = reduced_linear_max_iter;
    context->reduced_linear_tolerance =
        reduced_linear_tolerance;

    monolis_com_initialize_by_self(&context->mono_com0);

    context->initial_trace =
        (double*)calloc(trace_dof, sizeof(double));
    context->previous_trace =
        (double*)calloc(trace_dof, sizeof(double));

    context->current_window =
        (double*)calloc(local_st_dof, sizeof(double));
    context->lift =
        (double*)calloc(local_st_dof, sizeof(double));
    context->predictor =
        (double*)calloc(local_st_dof, sizeof(double));
    context->homogeneous =
        (double*)calloc(local_st_dof, sizeof(double));
    context->operator_input =
        (double*)calloc(local_st_dof, sizeof(double));

    context->reduced_correction =
        (double*)calloc(
            (size_t)reduced->num_modes_capacity,
            sizeof(double));

    if(context->initial_trace == NULL ||
       context->previous_trace == NULL ||
       context->current_window == NULL ||
       context->lift == NULL ||
       context->predictor == NULL ||
       context->homogeneous == NULL ||
       context->operator_input == NULL ||
       context->reduced_correction == NULL)
    {
        ROM_std_stdd_fluid_context_finalize(context);
        return ROM_STDD_ERR_ALLOCATION;
    }

    if(fluid_sups_st_export_state(
        fom,
        context->initial_trace) != 0)
    {
        ROM_std_stdd_fluid_context_finalize(context);
        return ROM_STDD_ERR_OPERATOR;
    }

    monolis_mpi_update_R(
        &fom->mono_com,
        fom->fe.total_num_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        context->initial_trace);

    memcpy(
        context->previous_trace,
        context->initial_trace,
        trace_dof * sizeof(double));

    return ROM_STDD_SUCCESS;
}

void ROM_std_stdd_fluid_context_finalize(
    ROM_STDD_FLUID_CONTEXT* context)
{
    if(context == NULL){
        return;
    }

    free(context->initial_trace);
    free(context->previous_trace);
    free(context->current_window);
    free(context->lift);
    free(context->predictor);
    free(context->homogeneous);
    free(context->operator_input);
    free(context->reduced_correction);

    memset(context, 0, sizeof(*context));
}

int ROM_std_stdd_fluid_enforce_homogeneous_basis(
    ROM_STDD_FLUID_CONTEXT* context)
{
    const double minimum_norm =
        128.0 * DBL_EPSILON;
    ROM_STDD_SYSTEM* reduced;

    if(context == NULL || context->fom == NULL ||
       context->reduced == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    reduced = context->reduced;

    if(reduced->basis == NULL || reduced->num_modes <= 0 ||
       reduced->num_modes > reduced->num_modes_capacity)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /*
     * The FOM correction BC contains both physical velocity Dirichlet rows
     * and the unique pressure-gauge rows.  The ROM homogeneous basis must
     * vanish on exactly the same rows.
     */
    for(int row = 0; row < reduced->internal_st_dof; row++){
        if(!context->fom->bc_NR_window.D_bc_exists[row]){
            continue;
        }

        for(int mode = 0; mode < reduced->num_modes; mode++){
            reduced->basis[
                fluid_rom_basis_offset(
                    reduced,
                    row,
                    mode)] = 0.0;
        }
    }

    /* Two-pass local modified Gram--Schmidt. */
    for(int mode = 0; mode < reduced->num_modes; mode++){
        for(int pass = 0; pass < 2; pass++){
            for(int previous = 0;
                previous < mode;
                previous++)
            {
                double dot = 0.0;

                for(int row = 0;
                    row < reduced->internal_st_dof;
                    row++)
                {
                    dot +=
                        reduced->basis[
                            fluid_rom_basis_offset(
                                reduced,
                                row,
                                previous)]
                        *
                        reduced->basis[
                            fluid_rom_basis_offset(
                                reduced,
                                row,
                                mode)];
                }

                for(int row = 0;
                    row < reduced->internal_st_dof;
                    row++)
                {
                    reduced->basis[
                        fluid_rom_basis_offset(
                            reduced,
                            row,
                            mode)] -=
                        dot
                        * reduced->basis[
                            fluid_rom_basis_offset(
                                reduced,
                                row,
                                previous)];
                }
            }
        }

        {
            double norm2 = 0.0;

            for(int row = 0;
                row < reduced->internal_st_dof;
                row++)
            {
                const double value =
                    reduced->basis[
                        fluid_rom_basis_offset(
                            reduced,
                            row,
                            mode)];
                norm2 += value * value;
            }

            if(!isfinite(norm2) ||
               norm2 <= minimum_norm * minimum_norm)
            {
                return ROM_STDD_ERR_DIMENSION;
            }

            {
                const double inverse_norm =
                    1.0 / sqrt(norm2);

                for(int row = 0;
                    row < reduced->internal_st_dof;
                    row++)
                {
                    reduced->basis[
                        fluid_rom_basis_offset(
                            reduced,
                            row,
                            mode)] *= inverse_norm;
                }
            }
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_fluid_import_hlpod_basis(
    ROM_STDD_FLUID_CONTEXT* context,
    const HLPOD_MAT* hlpod_mat,
    int num_modes)
{
    int status;

    if(context == NULL || context->reduced == NULL){
        return ROM_STDD_ERR_ARGUMENT;
    }

    status = ROM_std_stdd_import_hlpod_local_basis(
        context->reduced,
        hlpod_mat,
        num_modes);

    if(status != ROM_STDD_SUCCESS){
        return status;
    }

    return ROM_std_stdd_fluid_enforce_homogeneous_basis(
        context);
}

int ROM_std_stdd_fluid_store_hlpod_snapshot(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* full_window_state,
    HLPOD_MAT* hlpod_mat,
    int snapshot_id,
    int snapshot_capacity)
{
    if(context == NULL || full_window_state == NULL ||
       hlpod_mat == NULL || hlpod_mat->snapmat == NULL ||
       snapshot_id < 0 ||
       snapshot_id >= snapshot_capacity)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(fluid_sups_st_rom_build_window_lift(
        context->fom,
        context->layout,
        context->lift) != 0)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    for(int row = 0;
        row < context->reduced->internal_st_dof;
        row++)
    {
        hlpod_mat->snapmat[row][snapshot_id] =
            full_window_state[row]
            - context->lift[row];
    }

    return ROM_STDD_SUCCESS;
}

static int fluid_rom_apply_current_jacobian(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_FLUID_CONTEXT* context =
        (ROM_STDD_FLUID_CONTEXT*)user_context;

    if(context == NULL || x == NULL || y == NULL ||
       fluid_rom_validate(context) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memcpy(
        context->operator_input,
        x,
        (size_t)context->reduced->local_st_dof
            * sizeof(double));

    /*
     * Same verified convention as the convection--diffusion ST-DDROM:
     *
     *   1. exchange source-rank basis values explicitly,
     *   2. evaluate J*x with a communication-free communicator.
     */
    monolis_mpi_update_R(
        &context->fom->mono_com,
        context->fom->fe.total_num_nodes,
        context->fom->window_dof,
        context->operator_input);

    memset(
        y,
        0,
        (size_t)context->reduced->local_st_dof
            * sizeof(double));

    monolis_matvec_product_R(
        &context->fom->monolis_window,
        &context->mono_com0,
        context->operator_input,
        y);

    return ROM_STDD_SUCCESS;
}

static int fluid_rom_apply_zero_temporal_coupling(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_FLUID_CONTEXT* context =
        (ROM_STDD_FLUID_CONTEXT*)user_context;

    (void)x;

    if(context == NULL || y == NULL ||
       context->reduced == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(
        y,
        0,
        (size_t)context->reduced->local_st_dof
            * sizeof(double));

    return ROM_STDD_SUCCESS;
}

static int fluid_rom_project_predictor_coefficients(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id)
{
    double* coefficient;

    if(fluid_sups_st_rom_build_window_lift(
        context->fom,
        context->layout,
        context->lift) != 0 ||
       fluid_sups_st_rom_build_window_predictor(
        context->fom,
        context->layout,
        context->previous_trace,
        context->predictor) != 0)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    coefficient =
        context->reduced->reduced_coef
        + (size_t)window_id
            * (size_t)context->reduced->num_modes_capacity;

    memset(
        coefficient,
        0,
        (size_t)context->reduced->num_modes_capacity
            * sizeof(double));

    /*
     * Each spatial rank owns an independent local HLPOD space, so the
     * block-diagonal L2 projection is local.
     */
    for(int mode = 0;
        mode < context->reduced->num_modes;
        mode++)
    {
        double value = 0.0;

        for(int row = 0;
            row < context->reduced->internal_st_dof;
            row++)
        {
            value +=
                context->reduced->basis[
                    fluid_rom_basis_offset(
                        context->reduced,
                        row,
                        mode)]
                *
                (context->predictor[row]
                 - context->lift[row]);
        }

        coefficient[mode] = value;
    }

    return ROM_STDD_SUCCESS;
}

static int fluid_rom_reconstruct_window(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id)
{
    const double* coefficient =
        context->reduced->reduced_coef
        + (size_t)window_id
            * (size_t)context->reduced->num_modes_capacity;

    memset(
        context->homogeneous,
        0,
        (size_t)context->reduced->local_st_dof
            * sizeof(double));

    for(int row = 0;
        row < context->reduced->internal_st_dof;
        row++)
    {
        double value = 0.0;

        for(int mode = 0;
            mode < context->reduced->num_modes;
            mode++)
        {
            value +=
                context->reduced->basis[
                    fluid_rom_basis_offset(
                        context->reduced,
                        row,
                        mode)]
                * coefficient[mode];
        }

        context->homogeneous[row] = value;
    }

    monolis_mpi_update_R(
        &context->fom->mono_com,
        context->fom->fe.total_num_nodes,
        context->fom->window_dof,
        context->homogeneous);

    for(int row = 0;
        row < context->reduced->local_st_dof;
        row++)
    {
        context->current_window[row] =
            context->lift[row]
            + context->homogeneous[row];
    }

    if(fluid_sups_st_rom_enforce_window_constraints(
        context->fom,
        context->layout,
        context->current_window) != 0)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    return ROM_STDD_SUCCESS;
}

static double fluid_rom_relative_norm(
    double value,
    double initial)
{
    return initial > 1.0e-300
        ? value / initial
        : value;
}

static double fluid_rom_projected_rhs_norm(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* reduced_rhs)
{
    double norm2 = 0.0;

    for(int mode = 0;
        mode < context->reduced->num_modes;
        mode++)
    {
        norm2 += reduced_rhs[mode] * reduced_rhs[mode];
    }

    monolis_allreduce_R(
        1,
        &norm2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);

    return sqrt(fmax(norm2, 0.0));
}

static double fluid_rom_coefficient_correction_norm(
    ROM_STDD_FLUID_CONTEXT* context)
{
    double norm2 = 0.0;

    for(int mode = 0;
        mode < context->reduced->num_modes;
        mode++)
    {
        const double value =
            context->reduced_correction[mode];
        norm2 += value * value;
    }

    monolis_allreduce_R(
        1,
        &norm2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);

    return sqrt(fmax(norm2, 0.0));
}

static int fluid_rom_solve_reduced_correction(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id)
{
    ROM_STDD_SYSTEM one_window =
        *context->reduced;

    one_window.num_windows = 1;
    one_window.reduced_rhs =
        context->reduced->reduced_rhs
        + (size_t)window_id
            * (size_t)context->reduced->num_modes_capacity;
    one_window.reduced_coef =
        context->reduced_correction;

    memset(
        context->reduced_correction,
        0,
        (size_t)context->reduced->num_modes_capacity
            * sizeof(double));

    return ROM_std_stdd_solve_window_marching(
        &one_window,
        &context->fom->mono_com,
        context->reduced_linear_max_iter,
        context->reduced_linear_tolerance);
}

int ROM_std_stdd_fluid_solve_window(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    ROM_STDD_FLUID_WINDOW_REPORT* report)
{
    double initial_total = 0.0;
    double initial_projected = 0.0;
    double total = INFINITY;
    double momentum = INFINITY;
    double mass = INFINITY;
    double projected = INFINITY;
    double correction_norm = INFINITY;
    int converged = 0;

    fluid_rom_zero_report(report);

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(fluid_rom_project_predictor_coefficients(
        context,
        window_id) != ROM_STDD_SUCCESS ||
       fluid_rom_reconstruct_window(
        context,
        window_id) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    for(int iteration = 0;
        iteration < context->max_newton_iterations;
        iteration++)
    {
        double* reduced_rhs =
            context->reduced->reduced_rhs
            + (size_t)window_id
                * (size_t)context->reduced->num_modes_capacity;

        if(fluid_sups_st_rom_assemble_window_newton(
            context->fom,
            context->layout,
            context->previous_trace,
            context->current_window) != 0 ||
           fluid_sups_st_rom_get_window_rhs_norms(
            context->fom,
            context->layout,
            &total,
            &momentum,
            &mass) != 0 ||
           ROM_std_stdd_project_local_rhs(
            context->reduced,
            window_id,
            context->fom->monolis_window.mat.R.B,
            reduced_rhs) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        projected =
            fluid_rom_projected_rhs_norm(
                context,
                reduced_rhs);

        if(iteration == 0){
            initial_total = total;
            initial_projected = projected;
        }

        /*
         * Galerkin ROM convergence is defined by Phi^T R = 0.
         *
         * The full residual is retained as a model-error diagnostic and will
         * generally not vanish when the reduced basis is truncated.
         */
        if(projected <= context->nonlinear_absolute_tolerance ||
           fluid_rom_relative_norm(
               projected,
               initial_projected)
               <= context->nonlinear_tolerance)
        {
            converged = 1;
            correction_norm = 0.0;

            if(report != NULL){
                report->nonlinear_iterations = iteration;
            }
            break;
        }

        /*
         * Nonlinear NS differs from the linear convection--diffusion bridge:
         * J_r(a) must be rebuilt at every reduced Newton iteration.
         */
        if(ROM_std_stdd_project_operators_rank_sweep(
            context->reduced,
            fluid_rom_apply_current_jacobian,
            fluid_rom_apply_zero_temporal_coupling,
            context) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        if(fluid_rom_solve_reduced_correction(
            context,
            window_id) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        correction_norm =
            fluid_rom_coefficient_correction_norm(
                context);

        {
            double* coefficient =
                context->reduced->reduced_coef
                + (size_t)window_id
                    * (size_t)context->reduced->num_modes_capacity;

            for(int mode = 0;
                mode < context->reduced->num_modes;
                mode++)
            {
                coefficient[mode] +=
                    context->reduced_correction[mode];
            }
        }

        if(fluid_rom_reconstruct_window(
            context,
            window_id) != ROM_STDD_SUCCESS ||
           fluid_sups_st_rom_assemble_window_residual(
            context->fom,
            context->layout,
            context->previous_trace,
            context->current_window) != 0 ||
           fluid_sups_st_rom_get_window_rhs_norms(
            context->fom,
            context->layout,
            &total,
            &momentum,
            &mass) != 0 ||
           ROM_std_stdd_project_local_rhs(
            context->reduced,
            window_id,
            context->fom->monolis_window.mat.R.B,
            reduced_rhs) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        projected =
            fluid_rom_projected_rhs_norm(
                context,
                reduced_rhs);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "fluid ST-DDROM window %d Newton %d: "
                "projected rel=%e abs=%e ||da||=%e; "
                "full rel=%e abs=%e momentum=%e mass=%e\n",
                window_id + 1,
                iteration + 1,
                fluid_rom_relative_norm(
                    projected,
                    initial_projected),
                projected,
                correction_norm,
                fluid_rom_relative_norm(
                    total,
                    initial_total),
                total,
                momentum,
                mass);
        }

        if(projected <= context->nonlinear_absolute_tolerance ||
           fluid_rom_relative_norm(
               projected,
               initial_projected)
               <= context->nonlinear_tolerance)
        {
            converged = 1;

            if(report != NULL){
                report->nonlinear_iterations =
                    iteration + 1;
            }
            break;
        }
    }

    if(!converged){
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(
                stderr,
                "ERROR: fluid ST-DDROM projected Newton residual "
                "did not converge in window %d: "
                "projected=%e relative=%e full=%e\n",
                window_id + 1,
                projected,
                fluid_rom_relative_norm(
                    projected,
                    initial_projected),
                total);
        }
        return ROM_STDD_ERR_OPERATOR;
    }

    if(ST_fom_extract_last_right_trace(
        context->layout,
        context->current_window,
        context->previous_trace) != ST_FOM_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    monolis_mpi_update_R(
        &context->fom->mono_com,
        context->fom->fe.total_num_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        context->previous_trace);

    if(report != NULL){
        report->converged = 1;
        report->initial_full_residual =
            initial_total;
        report->final_full_residual =
            total;
        report->relative_full_residual =
            fluid_rom_relative_norm(
                total,
                initial_total);
        report->final_momentum_residual =
            momentum;
        report->final_mass_residual =
            mass;
        report->initial_projected_residual =
            initial_projected;
        report->final_projected_residual =
            projected;
        report->relative_projected_residual =
            fluid_rom_relative_norm(
                projected,
                initial_projected);
        report->final_coefficient_correction =
            correction_norm;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_fluid_write_current_window_vtk(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    int* file_number)
{
    double* trace = NULL;
    int status = ROM_STDD_SUCCESS;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       file_number == NULL ||
       window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    trace = (double*)calloc(
        (size_t)context->fom->fe.total_num_nodes
            * FLUID_SUPS_PHYSICAL_DOF,
        sizeof(double));
    if(trace == NULL){
        return ROM_STDD_ERR_ALLOCATION;
    }

    for(int slab = 0;
        slab < context->layout->num_slabs;
        slab++)
    {
        const int global_step =
            window_id * context->layout->num_slabs
            + slab + 1;
        const double time =
            (double)global_step
            * context->fom->vals.dt;

        if(global_step
               % context->fom->vals.output_interval
           != 0)
        {
            continue;
        }

        if(ST_fom_extract_slab_right_trace(
            context->layout,
            slab,
            context->current_window,
            trace) != ST_FOM_SUCCESS ||
           fluid_sups_st_import_state(
            context->fom,
            trace,
            0) != 0 ||
           fluid_sups_st_synchronize_state(
            context->fom) != 0 ||
           fluid_sups_st_write_output_named(
            context->fom,
            "st_ddrom_result",
            *file_number,
            time) != 0)
        {
            status = ROM_STDD_ERR_IO;
            goto cleanup;
        }

        (*file_number)++;
    }

    /*
     * Restore the right trace as the persistent physical state for the next
     * window and for post-processing after this output loop.
     */
    if(fluid_sups_st_import_state(
        context->fom,
        context->previous_trace,
        1) != 0 ||
       fluid_sups_st_synchronize_state(
        context->fom) != 0)
    {
        status = ROM_STDD_ERR_OPERATOR;
    }

cleanup:
    free(trace);
    return status;
}

int ROM_std_stdd_fluid_write_current_window_binary(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    const char* directory)
{
    char filename[4096];
    const int rank =
        monolis_mpi_get_global_my_rank();
    const double time_start =
        (double)(window_id
            * context->layout->num_slabs)
        * context->fom->vals.dt;
    const double time_end =
        time_start
        + (double)context->layout->num_slabs
            * context->fom->vals.dt;
    int written;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       directory == NULL || directory[0] == '\0' ||
       window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ST_fom_make_directory_recursive(directory)
       != ST_FOM_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    written = snprintf(
        filename,
        sizeof(filename),
        "%s/fluid_st_ddrom_window_%04d_rank_%06d.bin",
        directory,
        window_id,
        rank);

    if(written < 0 ||
       (size_t)written >= sizeof(filename))
    {
        return ROM_STDD_ERR_IO;
    }

    return ST_fom_write_owned_window_binary(
        filename,
        window_id,
        time_start,
        time_end,
        context->layout,
        context->current_window) == ST_FOM_SUCCESS
        ? ROM_STDD_SUCCESS
        : ROM_STDD_ERR_IO;
}

int ROM_std_stdd_fluid_owned_window_error(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* reference_owned_window,
    double* relative_error)
{
    double difference2 = 0.0;
    double reference2 = 0.0;
    double approximation2 = 0.0;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       reference_owned_window == NULL ||
       relative_error == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ST_fom_owned_difference_squares(
        context->layout,
        reference_owned_window,
        context->current_window,
        &difference2,
        &reference2,
        &approximation2) != ST_FOM_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    monolis_allreduce_R(
        1,
        &difference2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);
    monolis_allreduce_R(
        1,
        &reference2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);

    *relative_error =
        sqrt(fmax(difference2, 0.0))
        /
        fmax(
            sqrt(fmax(reference2, 0.0)),
            1.0e-300);

    (void)approximation2;
    return ROM_STDD_SUCCESS;
}
