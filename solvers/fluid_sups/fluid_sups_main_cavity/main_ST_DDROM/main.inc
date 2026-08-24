int main(int argc, char* argv[])
{
    FLUID_SUPS_SYSTEM fom;
    ST_FOM_WINDOW_LAYOUT layout;
    ROM_STDD_SYSTEM reduced;
    ROM_STDD_FLUID_CONTEXT adapter;
    FLUID_ST_DDROM_MODE mode;

    const char* snapshot_directory;
    const char* basis_directory;
    const char* rom_window_directory;
    const char* coefficient_csv_name;
    const char* accuracy_csv_name;

    int total_steps;
    int num_windows;
    int online_windows = 0;
    int available_snapshot_windows = 0;
    int num_training_windows = 0;
    int requested_training_windows = 0;
    int max_modes;
    int max_velocity_modes;
    int max_pressure_modes;
    int total_mode_capacity;
    int velocity_modes = 0;
    int pressure_modes = 0;
    int file_number = 0;
    int write_output;
    int write_window_binary;
    int validate_accuracy;
    double pod_energy;
    double velocity_pod_energy;
    double pressure_pod_energy;

    FILE* coefficient_csv = NULL;
    FILE* accuracy_csv = NULL;
    double* snapshot_matrix = NULL;
    double* reference_owned = NULL;
    double* reference_full = NULL;
    double* reference_homogeneous = NULL;
    double* projected_homogeneous = NULL;
    double* reference_coefficients = NULL;
    double* reference_previous_trace = NULL;
    double* reference_trace = NULL;
    double* rom_trace = NULL;
    double* reference_projected_rhs = NULL;
    DDROM_ACCURACY_SUMMARY accuracy_summary;

    int status = EXIT_FAILURE;

    memset(&fom, 0, sizeof(fom));
    memset(&layout, 0, sizeof(layout));
    memset(&reduced, 0, sizeof(reduced));
    memset(&adapter, 0, sizeof(adapter));
    memset(&accuracy_summary, 0, sizeof(accuracy_summary));

    monolis_global_initialize();

    if(fluid_sups_st_initialize(
        &fom,
        argc,
        argv) != 0)
    {
        ddrom_fail("fluid ST-DDROM FOM initialization failed");
    }

    total_steps = ddrom_total_steps(
        fom.vals.finish_time,
        fom.vals.dt);

    if(total_steps % fom.vals.st_slabs_per_window != 0){
        ddrom_fail(
            "total steps must be divisible by slabs/window");
    }

    num_windows =
        total_steps / fom.vals.st_slabs_per_window;

    if(ST_fom_layout_initialize(
        &layout,
        fom.fe.total_num_nodes,
        fom.num_owned_space_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        fom.vals.st_slabs_per_window,
        fom.vals.st_time_degree) != ST_FOM_SUCCESS)
    {
        ddrom_fail("ST-DDROM layout initialization failed");
    }

    /*
     * Backward-compatible fallback:
     *
     * FLUID_ST_DDROM_MAX_MODES=N means N velocity modes AND N pressure
     * modes unless the field-specific limits below are provided.
     */
    max_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_MODES",
        10);

    max_velocity_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_VELOCITY_MODES",
        max_modes);

    max_pressure_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_PRESSURE_MODES",
        max_modes);

    if(max_modes <= 0 ||
       max_velocity_modes <= 0 ||
       max_pressure_modes <= 0)
    {
        ddrom_fail(
            "ST-DDROM velocity/pressure mode limits must be positive");
    }

    if(max_velocity_modes > INT32_MAX - max_pressure_modes){
        ddrom_fail("ST-DDROM total mode capacity overflow");
    }

    total_mode_capacity =
        max_velocity_modes + max_pressure_modes;

    if(ROM_std_stdd_system_initialize(
        &reduced,
        &fom.mono_com,
        num_windows,
        fom.fe.total_num_nodes,
        fom.num_owned_space_nodes,
        fom.window_dof,
        total_mode_capacity) != ROM_STDD_SUCCESS)
    {
        ddrom_fail("ROM_STDD_SYSTEM initialization failed");
    }

    if(ROM_std_stdd_fluid_context_initialize(
        &adapter,
        &fom,
        &layout,
        &reduced,
        ddrom_env_int(
            "FLUID_ST_DDROM_NR_MAX_ITER",
            fom.vals.nr_max_iter),
        ddrom_env_double(
            "FLUID_ST_DDROM_NR_EPSILON",
            fom.vals.nr_epsilon),
        ddrom_env_double(
            "FLUID_ST_DDROM_NR_ABS_EPSILON",
            1.0e-12),
        ddrom_env_int(
            "FLUID_ST_DDROM_LINEAR_MAX_ITER",
            fom.vals.mat_max_iter),
        ddrom_env_double(
            "FLUID_ST_DDROM_LINEAR_EPSILON",
            fom.vals.mat_epsilon))
       != ROM_STDD_SUCCESS)
    {
        ddrom_fail("fluid ST-DDROM adapter initialization failed");
    }

    mode = ddrom_read_mode();

    snapshot_directory = ddrom_env_string(
        "FLUID_ST_DDROM_SNAPSHOT_DIR",
        "./fluid_st_windows");

    available_snapshot_windows =
        ddrom_count_common_training_windows(
            snapshot_directory,
            num_windows,
            fom.mono_com.comm);

    requested_training_windows =
        ddrom_env_int(
            "FLUID_ST_DDROM_TRAINING_WINDOWS",
            0);

    if(requested_training_windows < 0){
        ddrom_fail(
            "FLUID_ST_DDROM_TRAINING_WINDOWS must be >= 0");
    }

    if(requested_training_windows == 0){
        num_training_windows =
            available_snapshot_windows;
    }
    else{
        num_training_windows =
            requested_training_windows;
    }

    if((mode == FLUID_ST_DDROM_OFFLINE ||
        mode == FLUID_ST_DDROM_FULL) &&
       (num_training_windows <= 0 ||
        num_training_windows > available_snapshot_windows))
    {
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(
                stderr,
                "ERROR: insufficient common FOM training windows.\n"
                "  requested training windows : %d\n"
                "  available on every rank    : %d\n"
                "  online/runtime windows     : %d\n"
                "Set FLUID_ST_DDROM_TRAINING_WINDOWS to a value "
                "between 1 and the available count, or regenerate "
                "the missing FOM window binaries.\n",
                num_training_windows,
                available_snapshot_windows,
                num_windows);
        }
        ddrom_fail(
            "invalid ST-DDROM training-window count");
    }

    online_windows = ddrom_env_int(
        "FLUID_ST_DDROM_ONLINE_WINDOWS",
        num_windows);

    if(online_windows <= 0 ||
       online_windows > num_windows)
    {
        ddrom_fail(
            "FLUID_ST_DDROM_ONLINE_WINDOWS must be in [1,num_windows]");
    }

    basis_directory = ddrom_env_string(
        "FLUID_ST_DDROM_BASIS_DIR",
        "./fluid_st_ddrom_basis");

    rom_window_directory = ddrom_env_string(
        "FLUID_ST_DDROM_WINDOW_DIR",
        "./fluid_st_ddrom_windows");

    coefficient_csv_name = ddrom_env_string(
        "FLUID_ST_DDROM_COEFFICIENT_CSV",
        "./fluid_st_ddrom_coefficients.csv");

    accuracy_csv_name = ddrom_env_string(
        "FLUID_ST_DDROM_ACCURACY_CSV",
        "./fluid_st_ddrom_accuracy.csv");

    validate_accuracy = ddrom_env_int(
        "FLUID_ST_DDROM_VALIDATE",
        1);

    if(validate_accuracy != 0 &&
       validate_accuracy != 1)
    {
        ddrom_fail(
            "FLUID_ST_DDROM_VALIDATE must be 0 or 1");
    }

    pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_POD_ENERGY",
        1.0);

    velocity_pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_VELOCITY_POD_ENERGY",
        pod_energy);

    pressure_pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_PRESSURE_POD_ENERGY",
        pod_energy);

    if(velocity_pod_energy <= 0.0 ||
       velocity_pod_energy > 1.0 ||
       pressure_pod_energy <= 0.0 ||
       pressure_pod_energy > 1.0)
    {
        ddrom_fail(
            "velocity/pressure POD energies must be in (0,1]");
    }

    write_output = ddrom_env_int(
        "FLUID_ST_DDROM_WRITE_OUTPUT",
        1);

    write_window_binary = ddrom_env_int(
        "FLUID_ST_DDROM_WRITE_WINDOW_BINARY",
        1);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "\nfluid ST-DDROM configuration\n"
            "  mode                    : %s\n"
            "  dG degree               : %d\n"
            "  slabs/window            : %d\n"
            "  runtime window capacity : %d\n"
            "  online solve windows    : %d\n"
            "  FOM snapshots available : %d\n"
            "  POD training windows    : %d\n"
            "  physical MPI ranks      : %d\n"
            "  local/owned nodes rank0 : %d/%d\n"
            "  ST dof/node             : %d\n"
            "  basis structure         : separated velocity/pressure\n"
            "  max velocity modes/rank : %d\n"
            "  max pressure modes/rank : %d\n"
            "  total mode capacity/rank: %d\n"
            "  velocity POD energy     : %.15e\n"
            "  pressure POD energy     : %.15e\n"
            "  snapshot directory      : %s\n"
            "  basis directory         : %s\n"
            "  accuracy validation     : %s\n",
            mode == FLUID_ST_DDROM_OFFLINE
                ? "offline"
                : (mode == FLUID_ST_DDROM_ONLINE
                    ? "online"
                    : "full"),
            layout.time.degree,
            layout.num_slabs,
            num_windows,
            online_windows,
            available_snapshot_windows,
            num_training_windows,
            monolis_mpi_get_global_comm_size(),
            fom.fe.total_num_nodes,
            fom.num_owned_space_nodes,
            fom.window_dof,
            max_velocity_modes,
            max_pressure_modes,
            total_mode_capacity,
            velocity_pod_energy,
            pressure_pod_energy,
            snapshot_directory,
            basis_directory,
            validate_accuracy ? "enabled" : "disabled");
    }

    if(mode == FLUID_ST_DDROM_OFFLINE ||
       mode == FLUID_ST_DDROM_FULL)
    {
        if(ddrom_build_snapshot_matrix(
            &adapter,
            snapshot_directory,
            num_training_windows,
            &snapshot_matrix) != 0)
        {
            ddrom_fail(
                "failed to build the homogeneous ST snapshot matrix");
        }

        if(ddrom_compute_separated_velocity_pressure_pod(
            &adapter,
            snapshot_matrix,
            num_training_windows,
            max_velocity_modes,
            max_pressure_modes,
            velocity_pod_energy,
            pressure_pod_energy,
            fom.mono_com.comm,
            &velocity_modes,
            &pressure_modes) != ROM_STDD_SUCCESS)
        {
            ddrom_fail(
                "separated velocity/pressure space-time POD failed");
        }

        if(ROM_std_stdd_fluid_enforce_homogeneous_basis(
            &adapter) != ROM_STDD_SUCCESS)
        {
            ddrom_fail(
                "homogeneous fluid ST basis enforcement failed");
        }

        if(ROM_std_stdd_write_basis_rank(
            &reduced,
            num_training_windows,
            basis_directory) != ROM_STDD_SUCCESS)
        {
            ddrom_fail("ST-DDROM basis write failed");
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "fluid ST-DDROM offline completed: "
                "velocity modes/rank=%d, pressure modes/rank=%d, "
                "total modes/rank=%d, training windows=%d "
                "(runtime windows=%d)\n",
                velocity_modes,
                pressure_modes,
                reduced.num_modes,
                num_training_windows,
                num_windows);
        }
    }

    if(mode == FLUID_ST_DDROM_ONLINE ||
       mode == FLUID_ST_DDROM_FULL)
    {
        if(mode == FLUID_ST_DDROM_ONLINE){
            if(ROM_std_stdd_read_basis_rank(
                &reduced,
                basis_directory) != ROM_STDD_SUCCESS)
            {
                ddrom_fail("ST-DDROM basis read failed");
            }

            if(ROM_std_stdd_fluid_enforce_homogeneous_basis(
                &adapter) != ROM_STDD_SUCCESS)
            {
                ddrom_fail(
                    "online homogeneous basis enforcement failed");
            }
        }

        if(ddrom_classify_split_basis_modes(
               &reduced,
               &velocity_modes,
               &pressure_modes) != 0)
        {
            ddrom_fail(
                "loaded basis is not a separated velocity/pressure basis; "
                "rerun ST-DDROM offline with the split-basis implementation");
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "fluid ST-DDROM loaded split basis: "
                "velocity modes/rank=%d pressure modes/rank=%d "
                "total=%d\n",
                velocity_modes,
                pressure_modes,
                reduced.num_modes);
        }

        {
            char rank_csv[4096];
            const int rank =
                monolis_mpi_get_global_my_rank();
            int written = snprintf(
                rank_csv,
                sizeof(rank_csv),
                "%s.rank_%06d",
                coefficient_csv_name,
                rank);

            if(written < 0 ||
               (size_t)written >= sizeof(rank_csv))
            {
                ddrom_fail("coefficient CSV path is too long");
            }

            coefficient_csv = fopen(rank_csv, "w");
            if(coefficient_csv == NULL){
                ddrom_fail("cannot open coefficient CSV");
            }

            ddrom_write_coefficient_header(
                coefficient_csv,
                reduced.num_modes);
        }

        reference_owned = (double*)calloc(
            (size_t)reduced.internal_st_dof,
            sizeof(double));

        if(reference_owned == NULL){
            ddrom_fail("reference-window allocation failed");
        }

        if(validate_accuracy && available_snapshot_windows > 0){
            const size_t local_st_dof =
                (size_t)reduced.local_st_dof;
            const size_t internal_st_dof =
                (size_t)reduced.internal_st_dof;
            const size_t trace_dof =
                (size_t)fom.fe.total_num_nodes
                * (size_t)FLUID_SUPS_PHYSICAL_DOF;

            reference_full =
                (double*)calloc(
                    local_st_dof,
                    sizeof(double));
            reference_homogeneous =
                (double*)calloc(
                    internal_st_dof,
                    sizeof(double));
            projected_homogeneous =
                (double*)calloc(
                    internal_st_dof,
                    sizeof(double));
            reference_coefficients =
                (double*)calloc(
                    (size_t)reduced.num_modes_capacity,
                    sizeof(double));
            reference_previous_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            reference_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            rom_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            reference_projected_rhs =
                (double*)calloc(
                    (size_t)reduced.num_modes_capacity,
                    sizeof(double));

            if(reference_full == NULL ||
               reference_homogeneous == NULL ||
               projected_homogeneous == NULL ||
               reference_coefficients == NULL ||
               reference_previous_trace == NULL ||
               reference_trace == NULL ||
               rom_trace == NULL ||
               reference_projected_rhs == NULL)
            {
                ddrom_fail(
                    "ST-DDROM accuracy-validation allocation failed");
            }

            memcpy(
                reference_previous_trace,
                adapter.initial_trace,
                trace_dof * sizeof(double));

            if(monolis_mpi_get_global_my_rank() == 0){
                accuracy_csv = fopen(
                    accuracy_csv_name,
                    "w");
                if(accuracy_csv == NULL){
                    ddrom_fail(
                        "cannot open ST-DDROM accuracy CSV");
                }
                ddrom_write_accuracy_header(
                    accuracy_csv);
            }
        }

        for(int window = 0;
            window < online_windows;
            window++)
        {
            ROM_STDD_FLUID_WINDOW_REPORT report;
            double relative_reference_error = NAN;
            const double* coefficient;

            if(ROM_std_stdd_fluid_solve_window(
                &adapter,
                window,
                &report) != ROM_STDD_SUCCESS)
            {
                fprintf(
                    stderr,
                    "ERROR: fluid ST-DDROM failed in window %d\n",
                    window + 1);
                goto cleanup;
            }

            coefficient =
                reduced.reduced_coef
                + (size_t)window
                    * (size_t)reduced.num_modes_capacity;

            fprintf(coefficient_csv, "%d", window);
            for(int mode_id = 0;
                mode_id < reduced.num_modes;
                mode_id++)
            {
                fprintf(
                    coefficient_csv,
                    ",%.17e",
                    coefficient[mode_id]);
            }
            fprintf(
                coefficient_csv,
                ",%.17e,%.17e,%.17e\n",
                report.relative_full_residual,
                report.final_projected_residual,
                report.final_coefficient_correction);
            fflush(coefficient_csv);

            /*
             * Detailed FOM/ROM accuracy diagnostics are evaluated whenever a
             * matching reference window exists.  The reference uses its own
             * previous FOM trace for the residual check.
             */
            if(window < available_snapshot_windows &&
               ddrom_read_owned_window(
                   snapshot_directory,
                   window,
                   &layout,
                   reference_owned) == 0)
            {
                ROM_std_stdd_fluid_owned_window_error(
                    &adapter,
                    reference_owned,
                    &relative_reference_error);

                if(validate_accuracy){
                    DDROM_WINDOW_ACCURACY window_accuracy;

                    if(ddrom_validate_current_window(
                           &adapter,
                           window,
                           reference_owned,
                           reference_full,
                           reference_homogeneous,
                           projected_homogeneous,
                           reference_coefficients,
                           reference_previous_trace,
                           reference_trace,
                           rom_trace,
                           reference_projected_rhs,
                           &window_accuracy,
                           &accuracy_summary) != 0)
                    {
                        ddrom_fail(
                            "detailed ST-DDROM accuracy validation failed");
                    }

                    if(monolis_mpi_get_global_my_rank() == 0){
                        printf(
                            "fluid ST-DDROM accuracy window %d: "
                            "projection=%e state=%e dynamic=%e "
                            "velocity(state)=%e pressure(state)=%e "
                            "coef=%e trace=%e "
                            "FOM residual=%e FOM projected defect=%e\n",
                            window + 1,
                            window_accuracy.projection_error,
                            window_accuracy.state_error,
                            window_accuracy.dynamic_error,
                            window_accuracy.velocity_state_error,
                            window_accuracy.pressure_state_error,
                            window_accuracy.coefficient_error,
                            window_accuracy.trace_error,
                            window_accuracy.reference_full_residual,
                            window_accuracy.reference_projected_defect);

                        if(accuracy_csv != NULL){
                            fprintf(
                                accuracy_csv,
                                "%d,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,"
                                "%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,%.17e\n",
                                window,
                                window_accuracy.projection_error,
                                window_accuracy.state_error,
                                window_accuracy.dynamic_error,
                                window_accuracy.velocity_projection_error,
                                window_accuracy.velocity_state_error,
                                window_accuracy.velocity_dynamic_error,
                                window_accuracy.pressure_projection_error,
                                window_accuracy.pressure_state_error,
                                window_accuracy.pressure_dynamic_error,
                                window_accuracy.coefficient_error,
                                window_accuracy.trace_error,
                                window_accuracy.trace_velocity_error,
                                window_accuracy.trace_pressure_error,
                                report.final_projected_residual,
                                report.relative_projected_residual,
                                report.final_full_residual,
                                report.relative_full_residual,
                                window_accuracy.reference_full_residual,
                                window_accuracy.reference_momentum_residual,
                                window_accuracy.reference_mass_residual,
                                window_accuracy.reference_projected_defect);
                            fflush(accuracy_csv);
                        }
                    }
                }
            }

            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "fluid ST-DDROM window %d completed: "
                    "Newton=%d projected rel residual=%e "
                    "full rel residual=%e "
                    "FOM window error=%e\n",
                    window + 1,
                    report.nonlinear_iterations,
                    report.relative_projected_residual,
                    report.relative_full_residual,
                    relative_reference_error);
            }

            if(write_output){
                if(ROM_std_stdd_fluid_write_current_window_vtk(
                    &adapter,
                    window,
                    &file_number) != ROM_STDD_SUCCESS)
                {
                    ddrom_fail("ROM VTK output failed");
                }
            }

            if(write_window_binary){
                if(ROM_std_stdd_fluid_write_current_window_binary(
                    &adapter,
                    window,
                    rom_window_directory) != ROM_STDD_SUCCESS)
                {
                    ddrom_fail("ROM window binary output failed");
                }
            }
        }

        if(validate_accuracy &&
           accuracy_summary.num_validated_windows > 0)
        {
            ddrom_print_accuracy_summary(
                &accuracy_summary);

            if(monolis_mpi_get_global_my_rank() == 0 &&
               accuracy_summary.num_validated_windows <
                   online_windows)
            {
                printf(
                    "NOTE: accuracy was evaluated for %d/%d windows because "
                    "only %d common FOM reference windows are available.\n",
                    accuracy_summary.num_validated_windows,
                    online_windows,
                    available_snapshot_windows);
            }
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\nfluid ST-DDROM online completed successfully.\n");
        }
    }

    status = EXIT_SUCCESS;

cleanup:
    if(coefficient_csv != NULL){
        fclose(coefficient_csv);
    }
    if(accuracy_csv != NULL){
        fclose(accuracy_csv);
    }

    free(snapshot_matrix);
    free(reference_owned);
    free(reference_full);
    free(reference_homogeneous);
    free(projected_homogeneous);
    free(reference_coefficients);
    free(reference_previous_trace);
    free(reference_trace);
    free(rom_trace);
    free(reference_projected_rhs);

    ROM_std_stdd_fluid_context_finalize(&adapter);
    ROM_std_stdd_system_finalize(&reduced);
    fluid_sups_st_finalize(&fom);
    monolis_global_finalize();

    return status;
}
