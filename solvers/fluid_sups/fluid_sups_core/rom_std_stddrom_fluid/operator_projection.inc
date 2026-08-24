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

