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

