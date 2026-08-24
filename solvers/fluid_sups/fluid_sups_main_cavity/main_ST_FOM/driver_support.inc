typedef struct {
    FLUID_SUPS_SYSTEM* system;
    char window_directory[4096];
    FILE* checksum_csv;
    FILE* step_checksum_csv;
    int write_window_binary;
} FLUID_ST_DRIVER_CONTEXT;

typedef enum {
    FLUID_ST_SOLVE_CAUSAL = 0,
    FLUID_ST_SOLVE_ALL_AT_ONCE = 1
} FLUID_ST_SOLVE_MODE;

static FLUID_ST_SOLVE_MODE fluid_st_read_solve_mode(int time_degree)
{
    const char* value = getenv("FLUID_ST_SOLVE_MODE");

    if(value == NULL || value[0] == '\0'){
        return time_degree == 0
            ? FLUID_ST_SOLVE_CAUSAL
            : FLUID_ST_SOLVE_ALL_AT_ONCE;
    }

    if(strcmp(value, "causal") == 0 ||
       strcmp(value, "sequential") == 0)
    {
        return FLUID_ST_SOLVE_CAUSAL;
    }

    if(strcmp(value, "all_at_once") == 0 ||
       strcmp(value, "monolithic") == 0 ||
       strcmp(value, "window") == 0)
    {
        return FLUID_ST_SOLVE_ALL_AT_ONCE;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_SOLVE_MODE must be causal or all_at_once; got %s\n",
        value);
    exit(EXIT_FAILURE);
}


static void fluid_compute_state_norms(
    FLUID_SUPS_SYSTEM* system,
    double values[3])
{
    int node;
    const int owned = system->num_owned_space_nodes;

    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
    for(node = 0; node < owned; node++){
        const double vx = system->vals.v[node][0];
        const double vy = system->vals.v[node][1];
        const double vz = system->vals.v[node][2];
        const double pressure = system->vals.p[node];
        values[1] += vx * vx + vy * vy + vz * vz;
        values[2] += pressure * pressure;
    }
    values[0] = values[1] + values[2];
    monolis_allreduce_R(3, values, MONOLIS_MPI_SUM, system->mono_com.comm);
}
static int fluid_step_callback(
    void* user_context,
    int global_step,
    double time_old,
    double time_new)
{
    FLUID_ST_DRIVER_CONTEXT* context =
        (FLUID_ST_DRIVER_CONTEXT*)user_context;
    FLUID_SUPS_NR_REPORT report;
    int status;

    status = fluid_sups_st_solve_step(
        context->system,
        global_step,
        time_old,
        time_new,
        &report);

    {
        double norms[3];
        fluid_compute_state_norms(context->system, norms);
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s ST step %d completed: NR iterations=%d, "
                "total=%e, momentum=%e, mass=%e, ||state||=%e\n",
                CODENAME,
                global_step,
                report.nonlinear_iterations,
                report.relative_total_residual,
                report.relative_momentum_residual,
                report.relative_mass_residual,
                sqrt(fmax(norms[0], 0.0)));
            if(context->step_checksum_csv != NULL){
                fprintf(
                    context->step_checksum_csv,
                    "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                    global_step,
                    time_new,
                    sqrt(fmax(norms[0], 0.0)),
                    sqrt(fmax(norms[1], 0.0)),
                    sqrt(fmax(norms[2], 0.0)),
                    report.relative_total_residual,
                    report.relative_momentum_residual,
                    report.relative_mass_residual);
                fflush(context->step_checksum_csv);
            }
            fflush(stdout);
        }
    }

    return status;
}

static int fluid_export_callback(
    void* user_context,
    double* state)
{
    FLUID_ST_DRIVER_CONTEXT* context =
        (FLUID_ST_DRIVER_CONTEXT*)user_context;
    return fluid_sups_st_export_state(context->system, state);
}

static int fluid_output_callback(
    void* user_context,
    int file_number,
    int global_step,
    double time)
{
    FLUID_ST_DRIVER_CONTEXT* context =
        (FLUID_ST_DRIVER_CONTEXT*)user_context;
    (void)global_step;
    return fluid_sups_st_write_output(
        context->system, file_number, time);
}

static int fluid_window_callback(
    void* user_context,
    int window_id,
    double time_start,
    double time_end,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    const double* right_trace)
{
    FLUID_ST_DRIVER_CONTEXT* context =
        (FLUID_ST_DRIVER_CONTEXT*)user_context;
    double local_values[4] = {0.0, 0.0, 0.0, 0.0};
    int node;
    int slab;
    int alpha;
    int component;
    const int rank = monolis_mpi_get_global_my_rank();

    for(node = 0; node < layout->num_owned_space_points; node++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                for(component = 0;
                    component < layout->physical_dof;
                    component++)
                {
                    const size_t index = ST_fom_vector_index(
                        layout, node, slab, alpha, component);
                    const double value = window_state[index];
                    local_values[0] += value * value;
                    if(component < 3){
                        local_values[1] += value * value;
                    }
                    else{
                        local_values[2] += value * value;
                    }
                }
            }
        }
        for(component = 0; component < layout->physical_dof; component++){
            const double value = right_trace[
                (size_t)node * (size_t)layout->physical_dof
                + (size_t)component];
            local_values[3] += value * value;
        }
    }

    monolis_allreduce_R(
        4,
        local_values,
        MONOLIS_MPI_SUM,
        context->system->mono_com.comm);

    if(rank == 0){
        printf(
            "%s ST window %d [%e, %e]: "
            "||U_window||=%e ||v_window||=%e ||p_window||=%e "
            "||right trace||=%e\n",
            CODENAME,
            window_id + 1,
            time_start,
            time_end,
            sqrt(fmax(local_values[0], 0.0)),
            sqrt(fmax(local_values[1], 0.0)),
            sqrt(fmax(local_values[2], 0.0)),
            sqrt(fmax(local_values[3], 0.0)));

        if(context->checksum_csv != NULL){
            fprintf(
                context->checksum_csv,
                "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                window_id,
                time_start,
                time_end,
                sqrt(fmax(local_values[0], 0.0)),
                sqrt(fmax(local_values[1], 0.0)),
                sqrt(fmax(local_values[2], 0.0)),
                sqrt(fmax(local_values[3], 0.0)));
            fflush(context->checksum_csv);
        }
    }

    if(context->write_window_binary){
        char filename[4096];
        const int written = snprintf(
            filename,
            sizeof(filename),
            "%s/fluid_st_window_%04d_rank_%06d.bin",
            context->window_directory,
            window_id,
            rank);
        int status;

        if(written < 0 || (size_t)written >= sizeof(filename)){
            return 1;
        }
        status = ST_fom_write_owned_window_binary(
            filename,
            window_id,
            time_start,
            time_end,
            layout,
            window_state);
        if(status != ST_FOM_SUCCESS){
            fprintf(stderr,
                "ERROR [rank %d]: cannot write ST window file %s\n",
                rank, filename);
            return 1;
        }
    }

    return 0;
}


static int fluid_record_all_at_once_step(
    FLUID_ST_DRIVER_CONTEXT* context,
    const double* state,
    const FLUID_SUPS_NR_REPORT* report,
    int global_step,
    double time,
    int* file_number,
    int output_interval)
{
    double norms[3];

    if(fluid_sups_st_import_state(
        context->system, state, 0) != 0 ||
       fluid_sups_st_synchronize_state(
        context->system) != 0)
    {
        return 1;
    }

    fluid_compute_state_norms(context->system, norms);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "%s all-at-once reconstructed step %d: "
            "window Newton=%d total=%e momentum=%e mass=%e "
            "||state||=%e\n",
            CODENAME,
            global_step,
            report->nonlinear_iterations,
            report->relative_total_residual,
            report->relative_momentum_residual,
            report->relative_mass_residual,
            sqrt(fmax(norms[0], 0.0)));

        if(context->step_checksum_csv != NULL){
            fprintf(
                context->step_checksum_csv,
                "%d,%.17e,%.17e,%.17e,%.17e,"
                "%.17e,%.17e,%.17e\n",
                global_step,
                time,
                sqrt(fmax(norms[0], 0.0)),
                sqrt(fmax(norms[1], 0.0)),
                sqrt(fmax(norms[2], 0.0)),
                report->relative_total_residual,
                report->relative_momentum_residual,
                report->relative_mass_residual);
            fflush(context->step_checksum_csv);
        }
    }

    if(global_step % output_interval == 0){
        if(fluid_sups_st_write_output(
            context->system,
            *file_number,
            time) != 0)
        {
            return 1;
        }
        (*file_number)++;
    }

    return 0;
}

static int fluid_run_all_at_once_windows(
    FLUID_ST_DRIVER_CONTEXT* driver,
    const ST_FOM_CAUSAL_OPTIONS* options)
{
    FLUID_SUPS_SYSTEM* system = driver->system;
    const ST_FOM_WINDOW_LAYOUT* layout =
        &options->layout;
    const int slabs = layout->num_slabs;
    const int nstate =
        system->fe.total_num_nodes
        * FLUID_SUPS_PHYSICAL_DOF;
    const size_t window_length =
        ST_fom_window_vector_length(layout);
    double* previous_trace = NULL;
    double* window_state = NULL;
    double* slab_state = NULL;
    FLUID_SUPS_NR_REPORT* slab_reports = NULL;
    int file_number = 0;
    int window;
    int status = ST_FOM_ERR_CALLBACK;

    previous_trace =
        (double*)calloc((size_t)nstate, sizeof(double));
    window_state =
        (double*)calloc(window_length, sizeof(double));
    slab_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    slab_reports =
        (FLUID_SUPS_NR_REPORT*)calloc(
            (size_t)slabs,
            sizeof(*slab_reports));

    if(previous_trace == NULL ||
       window_state == NULL ||
       slab_state == NULL ||
       slab_reports == NULL)
    {
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    if(fluid_sups_st_export_state(
        system, previous_trace) != 0)
    {
        goto cleanup;
    }

    for(window = 0;
        window < options->num_windows;
        window++)
    {
        const double time_start =
            (double)(window * slabs) * options->dt;
        const double time_end =
            time_start + (double)slabs * options->dt;
        FLUID_SUPS_NR_REPORT window_report;
        int slab;

        memset(
            window_state,
            0,
            window_length * sizeof(double));
        memset(
            slab_reports,
            0,
            (size_t)slabs * sizeof(*slab_reports));
        memset(
            &window_report,
            0,
            sizeof(window_report));

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\n%s ################ ALL-AT-ONCE window %d/%d "
                "[%e, %e] ################\n",
                CODENAME,
                window + 1,
                options->num_windows,
                time_start,
                time_end);
        }

        if(fluid_sups_st_solve_window_all_at_once_dg(
            system,
            window + 1,
            layout,
            previous_trace,
            window_state,
            slab_reports,
            &window_report) != 0)
        {
            goto cleanup;
        }

        /*
         * Only after the complete window Newton has converged do we unpack
         * individual slab states for VTK/checksum output.
         */
        for(slab = 0; slab < slabs; slab++){
            const int global_step =
                window * slabs + slab + 1;
            const double time =
                (double)global_step * options->dt;

            if(ST_fom_extract_slab_right_trace(
                layout,
                slab,
                window_state,
                slab_state) != ST_FOM_SUCCESS)
            {
                goto cleanup;
            }

            if(fluid_record_all_at_once_step(
                driver,
                slab_state,
                &slab_reports[slab],
                global_step,
                time,
                &file_number,
                options->output_interval) != 0)
            {
                goto cleanup;
            }
        }

        if(ST_fom_extract_last_right_trace(
            layout,
            window_state,
            previous_trace) != ST_FOM_SUCCESS)
        {
            goto cleanup;
        }

        if(fluid_sups_st_import_state(
            system,
            previous_trace,
            1) != 0 ||
           fluid_sups_st_synchronize_state(
            system) != 0)
        {
            goto cleanup;
        }

        if(fluid_window_callback(
            driver,
            window,
            time_start,
            time_end,
            layout,
            window_state,
            previous_trace) != 0)
        {
            goto cleanup;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s ALL-AT-ONCE window %d completed: "
                "Newton=%d total=%e momentum=%e mass=%e\n",
                CODENAME,
                window + 1,
                window_report.nonlinear_iterations,
                window_report.relative_total_residual,
                window_report.relative_momentum_residual,
                window_report.relative_mass_residual);
        }
    }

    status = ST_FOM_SUCCESS;

cleanup:
    free(previous_trace);
    free(window_state);
    free(slab_state);
    free(slab_reports);
    return status;
}

