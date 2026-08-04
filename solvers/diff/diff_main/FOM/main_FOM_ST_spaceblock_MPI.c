#include "core_ROM.h"
#include "core_FOM_ST_spaceblock.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef STSB_DEGREE
#define STSB_DEGREE 1
#endif

#ifndef STSB_NUM_WINDOW_SLABS
#define STSB_NUM_WINDOW_SLABS 4
#endif

static int env_flag(const char* name, int default_value)
{
    const char* text = getenv(name);

    if(text == NULL || text[0] == '\0') {
        return default_value;
    }
    if(strcmp(text, "0") == 0 || strcmp(text, "false") == 0 ||
       strcmp(text, "FALSE") == 0 || strcmp(text, "off") == 0) {
        return 0;
    }
    return 1;
}

static double relative_norm(double numerator2, double denominator2)
{
    return sqrt(fmax(numerator2, 0.0)) /
        fmax(sqrt(fmax(denominator2, 0.0)), 1.0e-300);
}

int main(int argc, char* argv[])
{
    printf("\n");

    FE_SYSTEM sys;
    ST_TIME_ELEMENT te;
    STSB_SPATIAL_GRAPH graph;
    int num_windows = 0;
    int num_internal_space_nodes = 0;
    const int degree = STSB_DEGREE;
    const int num_window_slabs = STSB_NUM_WINDOW_SLABS;

    monolis_global_initialize();

    const int my_rank = monolis_mpi_get_global_my_rank();
    const double time_start_total = monolis_get_time();

    /*
     * This is the same initialization routine used by ST-DDROM collect.
     * It reconstructs the spatial graph from FE connectivity and assembles
     * the same space-time block matrix in both executables.
     */
    STSB_initialize_reference_fom(
        &sys,
        &te,
        &graph,
        argc,
        argv,
        degree,
        num_window_slabs,
        &num_windows,
        &num_internal_space_nodes);

    FILE* fp = NULL;
    fp = BBFE_sys_write_fopen(
        fp,
        "l2_error.txt",
        sys.cond.directory);
    fclose(fp);

    if(my_rank == 0) {
        fp = ROM_BB_write_fopen(
            fp,
            "calctime/time_fem.txt",
            sys.cond.directory);
        fclose(fp);

        ROM_std_hlpod_write_solver_prm_fopen(
            "fem_solver_prm",
            sys.cond.directory);
    }

    const int ndof = STSB_get_num_dofs_per_node(
        &te,
        num_window_slabs);

    const size_t num_window_values =
        (size_t)sys.fe.total_num_nodes * (size_t)ndof;

    if(num_window_values > (size_t)INT_MAX) {
        fprintf(
            stderr,
            "ERROR: BB_std allocation count exceeds INT_MAX.\n"
            "  requested values = %zu\n",
            num_window_values);
        exit(EXIT_FAILURE);
    }

    double* window_solution = BB_std_calloc_1d_double(
        window_solution,
        (int)num_window_values);

    double* previous_trace = BB_std_calloc_1d_double(
        previous_trace,
        sys.fe.total_num_nodes);

    if(window_solution == NULL || previous_trace == NULL) {
        fprintf(stderr, "ERROR: solution-vector allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys.fe,
        sys.vals.T,
        previous_trace);

    const int check_external_residual =
        env_flag("STSB_CHECK_EXTERNAL_RESIDUAL", 1);
    STSB_CONSTRAINED_RESIDUAL_SUM cumulative_residual;
    memset(&cumulative_residual, 0, sizeof(cumulative_residual));

    int file_num = 0;

    for(int window = 0; window < num_windows; window++) {
        const double window_start_time =
            (double)(window * num_window_slabs) * sys.vals.dt;

        memset(
            window_solution,
            0,
            num_window_values * sizeof(double));

        STSB_solve_window(
            &sys,
            &te,
            num_window_slabs,
            window_start_time,
            window + 1,
            previous_trace,
            window_solution);

        /* Read-only diagnostic, shared with ST-DDROM collect. */
        if(check_external_residual) {
            STSB_CONSTRAINED_RESIDUAL_SUM window_residual;
            memset(&window_residual, 0, sizeof(window_residual));

            if(STSB_evaluate_constrained_residual(
                &sys,
                ndof,
                window_solution,
                &window_residual) != 0) {
                fprintf(stderr, "ERROR: external residual evaluation failed.\n");
                exit(EXIT_FAILURE);
            }

            cumulative_residual.residual_all2 += window_residual.residual_all2;
            cumulative_residual.rhs_all2 += window_residual.rhs_all2;
            cumulative_residual.residual_free2 += window_residual.residual_free2;
            cumulative_residual.rhs_free2 += window_residual.rhs_free2;
            cumulative_residual.residual_dirichlet2 +=
                window_residual.residual_dirichlet2;
            cumulative_residual.rhs_dirichlet2 +=
                window_residual.rhs_dirichlet2;

            if(my_rank == 0) {
                const double all = relative_norm(
                    window_residual.residual_all2,
                    window_residual.rhs_all2);
                const double free_residual = relative_norm(
                    window_residual.residual_free2,
                    window_residual.rhs_free2);
                const double bc_all_rhs =
                    sqrt(fmax(window_residual.residual_dirichlet2, 0.0)) /
                    fmax(
                        sqrt(fmax(window_residual.rhs_all2, 0.0)),
                        1.0e-300);

                printf(
                    "STSB external residual window %d: "
                    "all=%.15e free=%.15e BC/all-rhs=%.15e\n",
                    window,
                    all,
                    free_residual,
                    bc_all_rhs);
            }
        }

        for(int local_slab = 0;
            local_slab < num_window_slabs;
            local_slab++) {

            const int global_step =
                window * num_window_slabs + local_slab + 1;

            const double time =
                (double)global_step * sys.vals.dt;

            if(global_step % sys.vals.output_interval != 0) {
                continue;
            }

            STSB_extract_slab_right_trace(
                &te,
                num_window_slabs,
                local_slab,
                sys.fe.total_num_nodes,
                window_solution,
                sys.vals.T);

            output_files(
                &sys,
                file_num,
                time);

            file_num++;
        }

        STSB_extract_last_right_trace(
            &te,
            num_window_slabs,
            sys.fe.total_num_nodes,
            window_solution,
            previous_trace);

        {
            double solution_l2 = 0.0;
            double trace_l2 = 0.0;

            if(STSB_compute_window_checksum(
                &sys,
                ndof,
                window_solution,
                previous_trace,
                &solution_l2,
                &trace_l2) != 0) {
                fprintf(stderr, "ERROR: FOM checksum evaluation failed.\n");
                exit(EXIT_FAILURE);
            }

            if(my_rank == 0) {
                printf(
                    "STSB reference checksum window %d: "
                    "solution=%.15e trace=%.15e\n",
                    window,
                    solution_l2,
                    trace_l2);
            }
        }
    }

    STSB_extract_last_right_trace(
        &te,
        num_window_slabs,
        sys.fe.total_num_nodes,
        window_solution,
        sys.vals.T);

    if(check_external_residual && my_rank == 0) {
        printf(
            "\nSTSB external residual summary\n"
            "  constrained residual all        : %.15e\n"
            "  constrained residual free       : %.15e\n"
            "  constrained residual BC/all-rhs : %.15e\n",
            relative_norm(
                cumulative_residual.residual_all2,
                cumulative_residual.rhs_all2),
            relative_norm(
                cumulative_residual.residual_free2,
                cumulative_residual.rhs_free2),
            sqrt(fmax(cumulative_residual.residual_dirichlet2, 0.0)) /
                fmax(
                    sqrt(fmax(cumulative_residual.rhs_all2, 0.0)),
                    1.0e-300));
    }

    BB_std_free_1d_double(
        previous_trace,
        sys.fe.total_num_nodes);

    BB_std_free_1d_double(
        window_solution,
        (int)num_window_values);

    STSB_finalize_reference_fom(&sys, &graph);

    const double time_end_total = monolis_get_time();

    if(my_rank == 0) {
        printf(
            "** Total time: %f\n",
            time_end_total - time_start_total);
        printf("\n");
    }

    monolis_global_finalize();
    return 0;
}
