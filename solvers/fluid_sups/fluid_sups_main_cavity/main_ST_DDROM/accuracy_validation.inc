static double ddrom_safe_ratio_from_squares(
    double numerator2,
    double denominator2)
{
    return sqrt(fmax(numerator2, 0.0))
        / fmax(sqrt(fmax(denominator2, 0.0)), 1.0e-300);
}

static double ddrom_reduced_vector_global_norm(
    const ROM_STDD_FLUID_CONTEXT* adapter,
    const double* vector)
{
    double norm2 = 0.0;

    for(int mode = 0;
        mode < adapter->reduced->num_modes;
        mode++)
    {
        norm2 += vector[mode] * vector[mode];
    }

    monolis_allreduce_R(
        1,
        &norm2,
        MONOLIS_MPI_SUM,
        adapter->fom->mono_com.comm);

    return sqrt(fmax(norm2, 0.0));
}

static int ddrom_reconstruct_reference_full(
    ROM_STDD_FLUID_CONTEXT* adapter,
    const double* reference_owned,
    double* reference_full)
{
    if(adapter == NULL || reference_owned == NULL ||
       reference_full == NULL)
    {
        return 1;
    }

    memset(
        reference_full,
        0,
        (size_t)adapter->reduced->local_st_dof
            * sizeof(double));

    /*
     * MONOLIS ownership contract for this distributed fluid path:
     * owned nodes are the local prefix, followed by halo nodes.
     */
    memcpy(
        reference_full,
        reference_owned,
        (size_t)adapter->reduced->internal_st_dof
            * sizeof(double));

    monolis_mpi_update_R(
        &adapter->fom->mono_com,
        adapter->fom->fe.total_num_nodes,
        adapter->fom->window_dof,
        reference_full);

    return 0;
}

static int ddrom_validate_current_window(
    ROM_STDD_FLUID_CONTEXT* adapter,
    int window_id,
    const double* reference_owned,
    double* reference_full,
    double* reference_homogeneous,
    double* projected_homogeneous,
    double* reference_coefficients,
    double* reference_previous_trace,
    double* reference_trace,
    double* rom_trace,
    double* reference_projected_rhs,
    DDROM_WINDOW_ACCURACY* accuracy,
    DDROM_ACCURACY_SUMMARY* summary)
{
    ROM_STDD_SYSTEM* reduced;
    double local[20] = {0.0};
    double reference_total_residual = 0.0;
    double reference_momentum_residual = 0.0;
    double reference_mass_residual = 0.0;
    const double* rom_coefficients;

    if(adapter == NULL || adapter->reduced == NULL ||
       adapter->layout == NULL || adapter->fom == NULL ||
       reference_owned == NULL || reference_full == NULL ||
       reference_homogeneous == NULL ||
       projected_homogeneous == NULL ||
       reference_coefficients == NULL ||
       reference_previous_trace == NULL ||
       reference_trace == NULL || rom_trace == NULL ||
       reference_projected_rhs == NULL ||
       accuracy == NULL || summary == NULL ||
       window_id < 0 ||
       window_id >= adapter->reduced->num_windows)
    {
        return 1;
    }

    reduced = adapter->reduced;
    memset(accuracy, 0, sizeof(*accuracy));

    if(fluid_sups_st_rom_build_window_lift(
           adapter->fom,
           adapter->layout,
           adapter->lift) != 0 ||
       ddrom_reconstruct_reference_full(
           adapter,
           reference_owned,
           reference_full) != 0)
    {
        return 1;
    }

    memset(
        reference_homogeneous,
        0,
        (size_t)reduced->internal_st_dof
            * sizeof(double));
    memset(
        projected_homogeneous,
        0,
        (size_t)reduced->internal_st_dof
            * sizeof(double));
    memset(
        reference_coefficients,
        0,
        (size_t)reduced->num_modes_capacity
            * sizeof(double));

    for(int row = 0; row < reduced->internal_st_dof; row++){
        reference_homogeneous[row] =
            reference_owned[row] - adapter->lift[row];
    }

    if(ROM_std_stdd_project_reference_coefficients(
           reduced,
           window_id,
           reference_homogeneous,
           reference_coefficients) != ROM_STDD_SUCCESS)
    {
        return 1;
    }

    for(int row = 0; row < reduced->internal_st_dof; row++){
        double value = 0.0;

        for(int mode = 0; mode < reduced->num_modes; mode++){
            value += reduced->basis[
                (size_t)row
                    * (size_t)reduced->num_modes_capacity
                + (size_t)mode]
                * reference_coefficients[mode];
        }

        projected_homogeneous[row] = value;
    }

    rom_coefficients =
        reduced->reduced_coef
        + (size_t)window_id
            * (size_t)reduced->num_modes_capacity;

    for(int row = 0; row < reduced->internal_st_dof; row++){
        const int physical_component =
            (row % reduced->st_dof_per_space_node)
                % FLUID_SUPS_PHYSICAL_DOF;
        const double q_ref =
            reference_homogeneous[row];
        const double q_projection =
            projected_homogeneous[row];
        const double q_rom =
            adapter->current_window[row]
            - adapter->lift[row];

        const double projection_difference =
            q_projection - q_ref;
        const double state_difference =
            q_rom - q_ref;
        const double dynamic_difference =
            q_rom - q_projection;

        local[0] += q_ref * q_ref;
        local[1] +=
            projection_difference
            * projection_difference;
        local[2] +=
            state_difference
            * state_difference;
        local[3] +=
            dynamic_difference
            * dynamic_difference;

        if(physical_component < 3){
            local[4] += q_ref * q_ref;
            local[5] +=
                projection_difference
                * projection_difference;
            local[6] +=
                state_difference
                * state_difference;
            local[7] +=
                dynamic_difference
                * dynamic_difference;
        }
        else{
            local[8] += q_ref * q_ref;
            local[9] +=
                projection_difference
                * projection_difference;
            local[10] +=
                state_difference
                * state_difference;
            local[11] +=
                dynamic_difference
                * dynamic_difference;
        }
    }

    for(int mode = 0; mode < reduced->num_modes; mode++){
        const double reference_value =
            reference_coefficients[mode];
        const double difference =
            rom_coefficients[mode]
            - reference_value;

        local[12] +=
            reference_value * reference_value;
        local[13] += difference * difference;
    }

    if(ST_fom_extract_last_right_trace(
           adapter->layout,
           reference_full,
           reference_trace) != ST_FOM_SUCCESS ||
       ST_fom_extract_last_right_trace(
           adapter->layout,
           adapter->current_window,
           rom_trace) != ST_FOM_SUCCESS)
    {
        return 1;
    }

    for(int node = 0;
        node < adapter->fom->num_owned_space_nodes;
        node++)
    {
        for(int component = 0;
            component < FLUID_SUPS_PHYSICAL_DOF;
            component++)
        {
            const int index =
                FLUID_SUPS_PHYSICAL_DOF * node
                + component;
            const double reference_value =
                reference_trace[index];
            const double difference =
                rom_trace[index]
                - reference_value;

            local[14] +=
                reference_value * reference_value;
            local[15] +=
                difference * difference;

            if(component < 3){
                local[16] +=
                    reference_value * reference_value;
                local[17] +=
                    difference * difference;
            }
            else{
                local[18] +=
                    reference_value * reference_value;
                local[19] +=
                    difference * difference;
            }
        }
    }

    monolis_allreduce_R(
        20,
        local,
        MONOLIS_MPI_SUM,
        adapter->fom->mono_com.comm);

    accuracy->projection_error =
        ddrom_safe_ratio_from_squares(
            local[1], local[0]);
    accuracy->state_error =
        ddrom_safe_ratio_from_squares(
            local[2], local[0]);
    accuracy->dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[3], local[0]);

    accuracy->velocity_projection_error =
        ddrom_safe_ratio_from_squares(
            local[5], local[4]);
    accuracy->velocity_state_error =
        ddrom_safe_ratio_from_squares(
            local[6], local[4]);
    accuracy->velocity_dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[7], local[4]);

    accuracy->pressure_projection_error =
        ddrom_safe_ratio_from_squares(
            local[9], local[8]);
    accuracy->pressure_state_error =
        ddrom_safe_ratio_from_squares(
            local[10], local[8]);
    accuracy->pressure_dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[11], local[8]);

    accuracy->coefficient_error =
        ddrom_safe_ratio_from_squares(
            local[13], local[12]);

    accuracy->trace_error =
        ddrom_safe_ratio_from_squares(
            local[15], local[14]);
    accuracy->trace_velocity_error =
        ddrom_safe_ratio_from_squares(
            local[17], local[16]);
    accuracy->trace_pressure_error =
        ddrom_safe_ratio_from_squares(
            local[19], local[18]);

    /*
     * Evaluate the stored FOM reference with its own causal previous trace,
     * not with the ROM previous trace.  This verifies that the binary
     * reference itself remains a solution of the all-at-once residual.
     */
    if(fluid_sups_st_rom_assemble_window_residual(
           adapter->fom,
           adapter->layout,
           reference_previous_trace,
           reference_full) != 0 ||
       fluid_sups_st_rom_get_window_rhs_norms(
           adapter->fom,
           adapter->layout,
           &reference_total_residual,
           &reference_momentum_residual,
           &reference_mass_residual) != 0 ||
       ROM_std_stdd_project_local_rhs(
           reduced,
           window_id,
           adapter->fom->monolis_window.mat.R.B,
           reference_projected_rhs) != ROM_STDD_SUCCESS)
    {
        return 1;
    }

    accuracy->reference_full_residual =
        reference_total_residual;
    accuracy->reference_momentum_residual =
        reference_momentum_residual;
    accuracy->reference_mass_residual =
        reference_mass_residual;
    accuracy->reference_projected_defect =
        ddrom_reduced_vector_global_norm(
            adapter,
            reference_projected_rhs);

    memcpy(
        reference_previous_trace,
        reference_trace,
        (size_t)adapter->fom->fe.total_num_nodes
            * (size_t)FLUID_SUPS_PHYSICAL_DOF
            * sizeof(double));

    summary->num_validated_windows++;

    summary->ref2 += local[0];
    summary->projection_diff2 += local[1];
    summary->state_diff2 += local[2];
    summary->dynamic_diff2 += local[3];

    summary->velocity_ref2 += local[4];
    summary->velocity_projection_diff2 += local[5];
    summary->velocity_state_diff2 += local[6];
    summary->velocity_dynamic_diff2 += local[7];

    summary->pressure_ref2 += local[8];
    summary->pressure_projection_diff2 += local[9];
    summary->pressure_state_diff2 += local[10];
    summary->pressure_dynamic_diff2 += local[11];

    summary->coefficient_ref2 += local[12];
    summary->coefficient_diff2 += local[13];

    summary->trace_ref2 += local[14];
    summary->trace_diff2 += local[15];
    summary->trace_velocity_ref2 += local[16];
    summary->trace_velocity_diff2 += local[17];
    summary->trace_pressure_ref2 += local[18];
    summary->trace_pressure_diff2 += local[19];

    summary->max_projection_error =
        fmax(summary->max_projection_error,
             accuracy->projection_error);
    summary->max_state_error =
        fmax(summary->max_state_error,
             accuracy->state_error);
    summary->max_dynamic_error =
        fmax(summary->max_dynamic_error,
             accuracy->dynamic_error);
    summary->max_velocity_state_error =
        fmax(summary->max_velocity_state_error,
             accuracy->velocity_state_error);
    summary->max_pressure_state_error =
        fmax(summary->max_pressure_state_error,
             accuracy->pressure_state_error);
    summary->max_coefficient_error =
        fmax(summary->max_coefficient_error,
             accuracy->coefficient_error);
    summary->max_trace_error =
        fmax(summary->max_trace_error,
             accuracy->trace_error);
    summary->max_reference_full_residual =
        fmax(summary->max_reference_full_residual,
             accuracy->reference_full_residual);
    summary->max_reference_projected_defect =
        fmax(summary->max_reference_projected_defect,
             accuracy->reference_projected_defect);

    return 0;
}

static void ddrom_print_accuracy_summary(
    const DDROM_ACCURACY_SUMMARY* summary)
{
    if(summary == NULL ||
       monolis_mpi_get_global_my_rank() != 0)
    {
        return;
    }

    printf(
        "\n"
        "============================================================\n"
        "Fluid ST-DDROM accuracy summary\n"
        "  validated windows                    : %d\n"
        "  projection error ||qF-PhiPhiTqF||/||qF|| : %.15e\n"
        "  ROM state error  ||qF-qROM||/||qF||      : %.15e\n"
        "  dynamic error    ||PhiPhiTqF-qROM||/||qF||: %.15e\n"
        "  velocity projection error            : %.15e\n"
        "  velocity ROM state error             : %.15e\n"
        "  velocity dynamic error               : %.15e\n"
        "  pressure projection error            : %.15e\n"
        "  pressure ROM state error             : %.15e\n"
        "  pressure dynamic error               : %.15e\n"
        "  reduced coefficient error            : %.15e\n"
        "  right-trace error                    : %.15e\n"
        "  right-trace velocity error           : %.15e\n"
        "  right-trace pressure error           : %.15e\n"
        "  max window projection error          : %.15e\n"
        "  max window ROM state error           : %.15e\n"
        "  max window dynamic error             : %.15e\n"
        "  max window velocity state error      : %.15e\n"
        "  max window pressure state error      : %.15e\n"
        "  max window coefficient error         : %.15e\n"
        "  max window right-trace error         : %.15e\n"
        "  max FOM-reference full residual      : %.15e\n"
        "  max projected defect at FOM reference: %.15e\n"
        "============================================================\n",
        summary->num_validated_windows,
        ddrom_safe_ratio_from_squares(
            summary->projection_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->state_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->dynamic_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_projection_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_state_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_dynamic_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_projection_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_state_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_dynamic_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->coefficient_diff2,
            summary->coefficient_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_diff2,
            summary->trace_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_velocity_diff2,
            summary->trace_velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_pressure_diff2,
            summary->trace_pressure_ref2),
        summary->max_projection_error,
        summary->max_state_error,
        summary->max_dynamic_error,
        summary->max_velocity_state_error,
        summary->max_pressure_state_error,
        summary->max_coefficient_error,
        summary->max_trace_error,
        summary->max_reference_full_residual,
        summary->max_reference_projected_defect);
}

static int ddrom_write_accuracy_header(FILE* fp)
{
    if(fp == NULL){
        return 1;
    }

    fprintf(
        fp,
        "window,"
        "projection_error,state_error,dynamic_error,"
        "velocity_projection_error,velocity_state_error,"
        "velocity_dynamic_error,"
        "pressure_projection_error,pressure_state_error,"
        "pressure_dynamic_error,"
        "coefficient_error,"
        "trace_error,trace_velocity_error,trace_pressure_error,"
        "rom_projected_residual,rom_relative_projected_residual,"
        "rom_full_residual,rom_relative_full_residual,"
        "reference_full_residual,reference_momentum_residual,"
        "reference_mass_residual,reference_projected_defect\n");

    return 0;
}

static int ddrom_write_coefficient_header(
    FILE* fp,
    int num_modes)
{
    if(fp == NULL || num_modes <= 0){
        return 1;
    }

    fprintf(fp, "window");
    for(int mode = 0; mode < num_modes; mode++){
        fprintf(fp, ",a_%04d", mode);
    }
    fprintf(
        fp,
        ",relative_full_residual,projected_residual,"
        "coefficient_correction\n");

    return 0;
}

