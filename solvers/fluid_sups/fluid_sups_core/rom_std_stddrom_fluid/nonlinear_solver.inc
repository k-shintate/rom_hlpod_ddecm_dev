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

