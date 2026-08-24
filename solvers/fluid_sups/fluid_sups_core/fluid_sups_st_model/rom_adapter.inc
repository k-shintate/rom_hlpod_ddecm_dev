/* ========================================================================== */
/* Public nonlinear FOM adapter used by the NS ST-DDROM                       */
/* ========================================================================== */

static int fluid_validate_rom_layout(
    const FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    if(system == NULL || layout == NULL ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->num_owned_space_points != system->num_owned_space_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       ST_fom_window_dof_per_space_point(layout) != system->window_dof)
    {
        return 1;
    }

    return 0;
}

int fluid_sups_st_rom_build_window_lift(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* lift)
{
    const size_t length =
        layout != NULL ? ST_fom_window_vector_length(layout) : 0u;

    if(fluid_validate_rom_layout(system, layout) != 0 ||
       lift == NULL || length == 0u)
    {
        return 1;
    }

    memset(lift, 0, length * sizeof(double));

    fluid_apply_physical_bc_to_window(
        system,
        layout,
        lift);
    fluid_apply_pressure_gauge_to_window_state(
        system,
        layout,
        lift);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        lift);

    return 0;
}

int fluid_sups_st_rom_build_window_predictor(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* predictor)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || predictor == NULL)
    {
        return 1;
    }

    fluid_dg_constant_predictor(
        system,
        layout,
        previous_trace,
        predictor);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        predictor);

    return 0;
}

int fluid_sups_st_rom_enforce_window_constraints(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       window_state == NULL)
    {
        return 1;
    }

    fluid_apply_physical_bc_to_window(
        system,
        layout,
        window_state);
    fluid_apply_pressure_gauge_to_window_state(
        system,
        layout,
        window_state);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        window_state);

    return 0;
}

int fluid_sups_st_rom_assemble_window_newton(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_assemble_window_newton_system(
            system,
            layout,
            previous_trace,
            window_state);
    }

    return fluid_dg_assemble_window(
        system,
        layout,
        previous_trace,
        window_state,
        1);
}

int fluid_sups_st_rom_assemble_window_residual(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_assemble_window_residual_only(
            system,
            layout,
            previous_trace,
            window_state);
    }

    return fluid_dg_assemble_window(
        system,
        layout,
        previous_trace,
        window_state,
        0);
}

int fluid_sups_st_rom_get_window_rhs_norms(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* total_norm,
    double* momentum_norm,
    double* mass_norm)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       total_norm == NULL || momentum_norm == NULL || mass_norm == NULL)
    {
        return 1;
    }

    fluid_window_rhs_norms(
        system,
        layout,
        total_norm,
        momentum_norm,
        mass_norm,
        NULL,
        NULL,
        NULL);

    return 0;
}


