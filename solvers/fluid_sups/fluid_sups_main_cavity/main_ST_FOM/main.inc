static int checked_total_steps(double finish_time, double dt)
{
    const int total_steps = (int)llround(finish_time / dt);
    const double reconstructed = (double)total_steps * dt;
    const double tolerance = 1.0e-12 * fmax(1.0, fabs(finish_time));

    if(total_steps <= 0 || fabs(reconstructed - finish_time) > tolerance){
        fprintf(stderr,
            "ERROR: finish_time/dt must be a positive integer; "
            "finish_time=%.17e dt=%.17e\n",
            finish_time, dt);
        exit(EXIT_FAILURE);
    }
    return total_steps;
}

int main(int argc, char* argv[])
{
    FLUID_SUPS_SYSTEM system;
    FLUID_ST_DRIVER_CONTEXT driver;
    ST_FOM_CAUSAL_OPTIONS options;
    ST_FOM_CAUSAL_CALLBACKS callbacks;
    int total_steps;
    int num_windows;
    int status;
    double start_time;
    double end_time;
    const char* directory_override;
    FLUID_ST_SOLVE_MODE solve_mode;

    memset(&system, 0, sizeof(system));
    memset(&driver, 0, sizeof(driver));
    memset(&options, 0, sizeof(options));
    memset(&callbacks, 0, sizeof(callbacks));

    monolis_global_initialize();
    start_time = monolis_get_time();

    if(fluid_sups_st_initialize(&system, argc, argv) != 0){
        fprintf(stderr, "ERROR: fluid ST-FOM initialization failed.\n");
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    solve_mode =
        fluid_st_read_solve_mode(system.vals.st_time_degree);

    if(solve_mode == FLUID_ST_SOLVE_CAUSAL &&
       system.vals.st_time_degree != 0)
    {
        fprintf(
            stderr,
            "ERROR: causal callback runner is only the validated dG(0) "
            "reference.  Use FLUID_ST_SOLVE_MODE=all_at_once for dG(%d).\n",
            system.vals.st_time_degree);
        fluid_sups_st_finalize(&system);
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    total_steps = checked_total_steps(
        system.vals.finish_time, system.vals.dt);
    if(total_steps % system.vals.st_slabs_per_window != 0){
        fprintf(stderr,
            "ERROR: total time steps (%d) must be divisible by "
            "st_slabs_per_window (%d).\n",
            total_steps, system.vals.st_slabs_per_window);
        fluid_sups_st_finalize(&system);
        monolis_global_finalize();
        return EXIT_FAILURE;
    }
    num_windows = total_steps / system.vals.st_slabs_per_window;

    directory_override = getenv("FLUID_ST_WINDOW_DIR");
    if(directory_override != NULL && directory_override[0] != '\0'){
        snprintf(driver.window_directory,
            sizeof(driver.window_directory), "%s", directory_override);
    }
    else{
        snprintf(driver.window_directory,
            sizeof(driver.window_directory), "%s/fluid_st_windows",
            system.cond.directory);
    }

    driver.system = &system;
    driver.write_window_binary = system.vals.write_window_binary;

    if(driver.write_window_binary &&
       ST_fom_make_directory_recursive(driver.window_directory)
           != ST_FOM_SUCCESS)
    {
        fprintf(stderr, "ERROR: cannot create %s\n", driver.window_directory);
        fluid_sups_st_finalize(&system);
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        char csv_filename[4096];
        snprintf(csv_filename, sizeof(csv_filename),
            "%s/fluid_st_window_checksums.csv", system.cond.directory);
        driver.checksum_csv = fopen(csv_filename, "w");
        if(driver.checksum_csv == NULL){
            fprintf(stderr, "ERROR: cannot open %s\n", csv_filename);
            fluid_sups_st_finalize(&system);
            monolis_global_finalize();
            return EXIT_FAILURE;
        }
        fprintf(driver.checksum_csv,
            "window,time_start,time_end,window_norm,velocity_norm,"
            "pressure_norm,right_trace_norm\n");

        snprintf(csv_filename, sizeof(csv_filename),
            "%s/fluid_st_step_checksums.csv", system.cond.directory);
        driver.step_checksum_csv = fopen(csv_filename, "w");
        if(driver.step_checksum_csv == NULL){
            fprintf(stderr, "ERROR: cannot open %s\n", csv_filename);
            fclose(driver.checksum_csv);
            fluid_sups_st_finalize(&system);
            monolis_global_finalize();
            return EXIT_FAILURE;
        }
        fprintf(driver.step_checksum_csv,
            "step,time,state_norm,velocity_norm,pressure_norm,"
            "relative_total_residual,relative_momentum_residual,"
            "relative_mass_residual\n");
    }

    status = ST_fom_layout_initialize(
        &options.layout,
        system.fe.total_num_nodes,
        system.num_owned_space_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        system.vals.st_slabs_per_window,
        system.vals.st_time_degree);
    if(status != ST_FOM_SUCCESS){
        fprintf(stderr, "ERROR: invalid fluid ST layout.\n");
        if(driver.checksum_csv != NULL){
            fclose(driver.checksum_csv);
        }
        if(driver.step_checksum_csv != NULL){
            fclose(driver.step_checksum_csv);
        }
        fluid_sups_st_finalize(&system);
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    if(options.layout.num_owned_space_points !=
           system.monolis_window.mat.N ||
       options.layout.num_space_points !=
           system.monolis_window.mat.NP)
    {
        fprintf(
            stderr,
            "ERROR [rank %d]: ST layout/window MONOLIS ownership mismatch: "
            "layout owned/local=%d/%d, matrix N/NP=%d/%d\n",
            monolis_mpi_get_global_my_rank(),
            options.layout.num_owned_space_points,
            options.layout.num_space_points,
            system.monolis_window.mat.N,
            system.monolis_window.mat.NP);

        if(driver.checksum_csv != NULL){
            fclose(driver.checksum_csv);
        }
        if(driver.step_checksum_csv != NULL){
            fclose(driver.step_checksum_csv);
        }
        fluid_sups_st_finalize(&system);
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    options.num_windows = num_windows;
    options.dt = system.vals.dt;
    options.output_interval = system.vals.output_interval;
    options.write_output = 1;

    callbacks.user_context = &driver;
    callbacks.solve_step = fluid_step_callback;
    callbacks.export_state = fluid_export_callback;
    callbacks.write_output = fluid_output_callback;
    callbacks.window_complete = fluid_window_callback;

    if(monolis_mpi_get_global_my_rank() == 0){
        if(solve_mode == FLUID_ST_SOLVE_ALL_AT_ONCE){
            printf(
                "%s Starting TRUE ALL-AT-ONCE window Newton ST-FOM: "
                "dG(%d), windows=%d, slabs/window=%d, "
                "time dof/slab=%d, window dof/node=%d\n",
                CODENAME,
                options.layout.time.degree,
                num_windows,
                system.vals.st_slabs_per_window,
                options.layout.time.n_dof,
                system.window_dof);
            printf(
                "%s Spatial MPI decomposition is inherited from the physical "
                "fluid mesh; the full window block is distributed by spatial "
                "node exactly as in the convection-diffusion STSB FOM.\n",
                CODENAME);
            printf(
                "%s Each Newton iteration assembles the complete "
                "space-time window Jacobian and performs ONE MONOLIS solve. "
                "dG(1)/dG(2) use temporal Gauss integration of the verified "
                "SUPS residual plus the upwind jump operator.\n",
                CODENAME);
        }
        else{
            printf(
                "%s Starting causal space-time FOM: dG(0), windows=%d, "
                "slabs/window=%d, physical dof/node=%d\n",
                CODENAME,
                num_windows,
                system.vals.st_slabs_per_window,
                FLUID_SUPS_PHYSICAL_DOF);
            printf(
                "%s The causal solve is exact block forward substitution of "
                "the dG(0) nonlinear ST window system and reuses the verified "
                "NR one-step kernel without changing its algebra.\n",
                CODENAME);
        }
    }

    status = solve_mode == FLUID_ST_SOLVE_ALL_AT_ONCE
        ? fluid_run_all_at_once_windows(&driver, &options)
        : ST_fom_run_causal_dg0_windows(&options, &callbacks);

    if(driver.checksum_csv != NULL){
        fclose(driver.checksum_csv);
    }
    if(driver.step_checksum_csv != NULL){
        fclose(driver.step_checksum_csv);
    }

    fluid_sups_st_finalize(&system);
    end_time = monolis_get_time();

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("** Total time: %f\n", end_time - start_time);
    }

    monolis_global_finalize();
    return status == ST_FOM_SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
