#include "fluid_sups_st_model.h"
#include "st_fom_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    FLUID_SUPS_SYSTEM* system;
    char window_directory[4096];
    FILE* checksum_csv;
    FILE* step_checksum_csv;
    int write_window_binary;
} FLUID_ST_DRIVER_CONTEXT;

static int fluid_owned_space_points(const FLUID_SUPS_SYSTEM* system)
{
    int owned = 0;
    if(system == NULL ||
       ST_fom_resolve_owned_space_points(
           system->fe.total_num_nodes,
           monolis_mpi_get_global_comm_size(),
           system->mono_com.n_internal_vertex,
           &owned) != ST_FOM_SUCCESS)
    {
        return 0;
    }
    return owned;
}


static void fluid_compute_state_norms(
    FLUID_SUPS_SYSTEM* system,
    double values[3])
{
    int node;
    const int owned = fluid_owned_space_points(system);

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
    int component;
    const int rank = monolis_mpi_get_global_my_rank();

    for(node = 0; node < layout->num_owned_space_points; node++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            for(component = 0; component < layout->physical_dof; component++){
                const size_t index = ST_fom_vector_index(
                    layout, node, slab, 0, component);
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
        fluid_owned_space_points(&system),
        FLUID_SUPS_PHYSICAL_DOF,
        system.vals.st_slabs_per_window,
        0);
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
        printf(
            "%s Starting causal space-time FOM: dG(0), windows=%d, "
            "slabs/window=%d, physical dof/node=%d\n",
            CODENAME,
            num_windows,
            system.vals.st_slabs_per_window,
            FLUID_SUPS_PHYSICAL_DOF);
        printf(
            "%s The causal solve is exact block forward substitution of the "
            "dG(0) nonlinear ST window system and reuses the verified NR "
            "one-step kernel without changing its algebra.\n",
            CODENAME);
    }

    status = ST_fom_run_causal_dg0_windows(&options, &callbacks);

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
