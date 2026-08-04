#include "core_ROM.h"
#include "core_FOM_ST_spaceblock.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef STSB_DEGREE
#define STSB_DEGREE 1
#endif

#ifndef STSB_NUM_WINDOW_SLABS
#define STSB_NUM_WINDOW_SLABS 4
#endif

#ifndef STSB_GRAPH_LABEL
#define STSB_GRAPH_LABEL "graph.dat"
#endif

int main(int argc, char* argv[])
{
    printf("\n");

    FE_SYSTEM sys;

    monolis_global_initialize();

    const double time_start_total = monolis_get_time();

    sys.cond.directory = BBFE_convdiff_get_directory_name(
        argc,
        argv,
        CODENAME);

    read_calc_conditions(
        &sys.vals,
        sys.cond.directory);

    /*
     * node.dat, elem.dat, graph.dat, D_bc.dat and MONOLIS communication
     * files must all correspond to the same ordinary spatial partition.
     */
    BBFE_convdiff_pre(
        &sys.fe,
        &sys.basis,
        &sys.bc,
        &sys.monolis,
        &sys.monolis_com,
        argc,
        argv,
        sys.cond.directory,
        sys.vals.num_ip_each_axis,
        true);

    memory_allocation_nodal_values(
        &sys.vals,
        sys.fe.total_num_nodes);

    manusol_set_init_value(
        &sys.fe,
        sys.vals.T);

    FILE* fp = NULL;
    fp = BBFE_sys_write_fopen(
        fp,
        "l2_error.txt",
        sys.cond.directory);
    fclose(fp);

    if(monolis_mpi_get_global_my_rank() == 0) {
        fp = ROM_BB_write_fopen(
            fp,
            "calctime/time_fem.txt",
            sys.cond.directory);
        fclose(fp);

        ROM_std_hlpod_write_solver_prm_fopen(
            "fem_solver_prm",
            sys.cond.directory);
    }

    BBFE_elemmat_set_Jacobi_mat(
        &sys.fe,
        &sys.basis);

    BBFE_elemmat_set_shapefunc_derivative(
        &sys.fe,
        &sys.basis);

    const int degree = STSB_DEGREE;
    const int num_window_slabs = STSB_NUM_WINDOW_SLABS;

    ST_TIME_ELEMENT te;
    ST_TimeElement_init_dG(
        &te,
        degree);

    STSB_validate_configuration(
        &sys.fe,
        &sys.vals,
        &te,
        num_window_slabs);

    const int num_windows = STSB_get_num_windows(
        sys.vals.finish_time,
        sys.vals.dt,
        num_window_slabs);

    const int ndof = STSB_get_num_dofs_per_node(
        &te,
        num_window_slabs);

    /*
     * BBFE_convdiff_pre has already initialized the ordinary spatial
     * MONOLIS matrix.  Preserve its owned-node count before replacing its
     * scalar matrix by the ST block matrix.
     */
    const int num_internal_space_nodes = sys.monolis.mat.N;

    if(num_internal_space_nodes <= 0 ||
       num_internal_space_nodes > sys.fe.total_num_nodes) {
        fprintf(
            stderr,
            "ERROR: invalid internal-node count from spatial MONOLIS data.\n"
            "  internal nodes = %d\n"
            "  total nodes    = %d\n",
            num_internal_space_nodes,
            sys.fe.total_num_nodes);
        exit(EXIT_FAILURE);
    }

    STSB_SPATIAL_GRAPH graph;
    STSB_spatial_graph_initialize(&graph);

    STSB_spatial_graph_build_from_mesh(
		&graph,
		&sys.fe);

    STSB_validate_spatial_graph_against_mesh(
		&graph,
		&sys.fe);

    monolis_initialize(&sys.monolis0);

    STSB_set_nonzero_pattern_from_spatial_graph(
        &sys.monolis0,
        &graph,
        &te,
        num_window_slabs,
        num_internal_space_nodes);

    monolis_clear_mat_value_rhs_R(&sys.monolis0);

    STSB_assemble_window_matrix(
        &sys.monolis0,
        &sys.fe,
        &sys.basis,
        &sys.vals,
        &te,
        num_window_slabs);

    /* Copy complete block structure and values into the working matrix. */
    monolis_copy_mat_R(
        &sys.monolis0,
        &sys.monolis);

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

    STSB_set_previous_trace_from_nodal_value(
        &sys.fe,
        sys.vals.T,
        previous_trace);

    int file_num = 0;

    for(int window=0; window<num_windows; window++) {
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

        for(int local_slab=0;
            local_slab<num_window_slabs;
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

        /* No ST-to-spatial redistribution is necessary. */
        STSB_extract_last_right_trace(
            &te,
            num_window_slabs,
            sys.fe.total_num_nodes,
            window_solution,
            previous_trace);
    }

    STSB_extract_last_right_trace(
        &te,
        num_window_slabs,
        sys.fe.total_num_nodes,
        window_solution,
        sys.vals.T);

    BB_std_free_1d_double(
        previous_trace,
        sys.fe.total_num_nodes);

    BB_std_free_1d_double(
        window_solution,
        (int)num_window_values);

    STSB_spatial_graph_finalize(&graph);

    BBFE_convdiff_finalize(
        &sys.fe,
        &sys.basis,
        &sys.bc);

    monolis_finalize(&sys.monolis);
    monolis_finalize(&sys.monolis0);

    monolis_global_finalize();

    const double time_end_total = monolis_get_time();

    printf(
        "** Total time: %f\n",
        time_end_total - time_start_total);

    printf("\n");

    return 0;
}
