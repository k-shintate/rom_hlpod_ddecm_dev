#include "fluid_sups_st_model.h"
#include "st_fom_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static void compute_state_norms(
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

int main(int argc, char* argv[])
{
    FLUID_SUPS_SYSTEM system;
    int file_number = 0;
    int total_steps;
    int step;
    double start_time;
    double end_time;
    FILE* checksum_csv = NULL;
    int exit_status = EXIT_SUCCESS;

    memset(&system, 0, sizeof(system));
    monolis_global_initialize();
    start_time = monolis_get_time();

    if(fluid_sups_st_initialize(&system, argc, argv) != 0){
        fprintf(stderr, "ERROR: fluid baseline initialization failed.\n");
        monolis_global_finalize();
        return EXIT_FAILURE;
    }

    total_steps = checked_total_steps(
        system.vals.finish_time, system.vals.dt);

    if(monolis_mpi_get_global_my_rank() == 0){
        char filename[4096];
        snprintf(filename, sizeof(filename),
            "%s/fluid_nr_step_checksums.csv", system.cond.directory);
        checksum_csv = fopen(filename, "w");
        if(checksum_csv == NULL){
            fprintf(stderr, "ERROR: cannot open %s\n", filename);
            exit_status = EXIT_FAILURE;
            goto cleanup;
        }
        fprintf(checksum_csv,
            "step,time,state_norm,velocity_norm,pressure_norm,"
            "relative_total_residual,relative_momentum_residual,"
            "relative_mass_residual\n");
    }

    for(step = 1; step <= total_steps; step++){
        const double time_old = (double)(step - 1) * system.vals.dt;
        const double time_new = (double)step * system.vals.dt;
        FLUID_SUPS_NR_REPORT report;
        double norms[3];

        if(fluid_sups_st_solve_step(
            &system, step, time_old, time_new, &report) != 0)
        {
            exit_status = EXIT_FAILURE;
            goto cleanup;
        }

        compute_state_norms(&system, norms);
        if(checksum_csv != NULL){
            fprintf(checksum_csv,
                "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                step, time_new,
                sqrt(fmax(norms[0], 0.0)),
                sqrt(fmax(norms[1], 0.0)),
                sqrt(fmax(norms[2], 0.0)),
                report.relative_total_residual,
                report.relative_momentum_residual,
                report.relative_mass_residual);
            fflush(checksum_csv);
        }

        if(step % system.vals.output_interval == 0){
            if(fluid_sups_st_write_output(
                &system, file_number, time_new) != 0)
            {
                exit_status = EXIT_FAILURE;
                goto cleanup;
            }
            file_number++;
        }
    }

cleanup:
    if(checksum_csv != NULL){
        fclose(checksum_csv);
    }
    fluid_sups_st_finalize(&system);
    end_time = monolis_get_time();
    if(monolis_mpi_get_global_my_rank() == 0){
        printf("** Total time: %f\n", end_time - start_time);
    }
    monolis_global_finalize();
    return exit_status;
}
