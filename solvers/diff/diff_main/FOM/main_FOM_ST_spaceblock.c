#include "core_FOM_ST_spaceblock.h"
#include "st_fom_common.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int env_int(const char* name, int default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    long value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }
    errno = 0;
    value = strtol(text, &end, 10);
    if(errno != 0 || end == text || *end != '\0' || value < 0){
        fprintf(stderr, "ERROR: invalid integer %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static double env_double(const char* name, double default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    double value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }
    errno = 0;
    value = strtod(text, &end);
    if(errno != 0 || end == text || *end != '\0' || !isfinite(value)){
        fprintf(stderr, "ERROR: invalid real %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }
    return value;
}

static void write_compare_vtk(
    FE_SYSTEM* system,
    int window,
    const double* monolithic,
    const double* causal)
{
    char filename[4096];
    double* absolute_error;
    FILE* fp = NULL;
    int node;
    const int rank = monolis_mpi_get_global_my_rank();

    absolute_error = (double*)calloc(
        (size_t)system->fe.total_num_nodes, sizeof(double));
    if(absolute_error == NULL){
        fprintf(stderr, "ERROR: VTK comparison allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    for(node = 0; node < system->fe.total_num_nodes; node++){
        absolute_error[node] = fabs(monolithic[node] - causal[node]);
    }

    snprintf(filename, sizeof(filename),
        "convdiff_stfom_dg0_compare_%06d.vtk.%d", window + 1, rank);
    fp = BBFE_sys_write_fopen(fp, filename, system->cond.directory);

    if(system->fe.local_num_nodes == 4){
        BBFE_sys_write_vtk_shape(fp, &system->fe, TYPE_VTK_TETRA);
    }
    else if(system->fe.local_num_nodes == 8){
        BBFE_sys_write_vtk_shape(fp, &system->fe, TYPE_VTK_HEXAHEDRON);
    }
    else{
        fprintf(stderr, "ERROR: unsupported element for comparison VTK.\n");
        fclose(fp);
        free(absolute_error);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "POINT_DATA %d\n", system->fe.total_num_nodes);
    BB_vtk_write_point_vals_scalar(
        fp, (double*)monolithic, system->fe.total_num_nodes,
        "monolithic_dG0_right_trace");
    BB_vtk_write_point_vals_scalar(
        fp, (double*)causal, system->fe.total_num_nodes,
        "causal_dG0_right_trace");
    BB_vtk_write_point_vals_scalar(
        fp, absolute_error, system->fe.total_num_nodes,
        "absolute_difference");
    fclose(fp);
    free(absolute_error);
}

int main(int argc, char* argv[])
{
    FE_SYSTEM monolithic_system;
    FE_SYSTEM causal_system;
    ST_TIME_ELEMENT monolithic_time;
    ST_TIME_ELEMENT causal_time;
    STSB_SPATIAL_GRAPH monolithic_graph;
    STSB_SPATIAL_GRAPH causal_graph;
    ST_FOM_WINDOW_LAYOUT common_layout;
    int monolithic_windows = 0;
    int causal_windows = 0;
    int monolithic_owned = 0;
    int causal_owned = 0;
    int slabs = env_int("CONVDIFF_STFOM_SLABS_PER_WINDOW", 4);
    int write_vtk = env_int("CONVDIFF_STFOM_WRITE_COMPARE_VTK", 1);
    double tolerance = env_double("CONVDIFF_STFOM_COMPARE_TOL", 1.0e-10);
    double* monolithic_previous = NULL;
    double* monolithic_trace = NULL;
    double* causal_previous = NULL;
    double* causal_one_slab = NULL;
    double* monolithic_window = NULL;
    double* causal_window = NULL;
    FILE* csv = NULL;
    double max_window_error = 0.0;
    double max_trace_error = 0.0;
    int window;
    int status = EXIT_SUCCESS;

    if(slabs <= 0 || tolerance <= 0.0){
        fprintf(stderr, "ERROR: invalid validation options.\n");
        return EXIT_FAILURE;
    }

    memset(&monolithic_system, 0, sizeof(monolithic_system));
    memset(&causal_system, 0, sizeof(causal_system));
    memset(&monolithic_time, 0, sizeof(monolithic_time));
    memset(&causal_time, 0, sizeof(causal_time));
    STSB_spatial_graph_initialize(&monolithic_graph);
    STSB_spatial_graph_initialize(&causal_graph);

    monolis_global_initialize();

    STSB_initialize_reference_fom(
        &monolithic_system,
        &monolithic_time,
        &monolithic_graph,
        argc,
        argv,
        0,
        slabs,
        &monolithic_windows,
        &monolithic_owned);

    STSB_initialize_reference_fom(
        &causal_system,
        &causal_time,
        &causal_graph,
        argc,
        argv,
        0,
        1,
        &causal_windows,
        &causal_owned);

    if(monolithic_system.fe.total_num_nodes !=
           causal_system.fe.total_num_nodes ||
       monolithic_owned != causal_owned ||
       causal_windows != monolithic_windows * slabs)
    {
        fprintf(stderr,
            "ERROR: monolithic/causal dG(0) FOM dimensions are inconsistent.\n");
        status = EXIT_FAILURE;
        goto cleanup;
    }

    if(ST_fom_layout_initialize(
        &common_layout,
        monolithic_system.fe.total_num_nodes,
        monolithic_owned,
        1,
        slabs,
        0) != ST_FOM_SUCCESS)
    {
        fprintf(stderr, "ERROR: common dG(0) layout initialization failed.\n");
        status = EXIT_FAILURE;
        goto cleanup;
    }

    monolithic_previous = (double*)calloc(
        (size_t)monolithic_system.fe.total_num_nodes, sizeof(double));
    monolithic_trace = (double*)calloc(
        (size_t)monolithic_system.fe.total_num_nodes, sizeof(double));
    causal_previous = (double*)calloc(
        (size_t)causal_system.fe.total_num_nodes, sizeof(double));
    causal_one_slab = (double*)calloc(
        (size_t)causal_system.fe.total_num_nodes, sizeof(double));
    monolithic_window = (double*)calloc(
        ST_fom_window_vector_length(&common_layout), sizeof(double));
    causal_window = (double*)calloc(
        ST_fom_window_vector_length(&common_layout), sizeof(double));

    if(monolithic_previous == NULL || monolithic_trace == NULL ||
       causal_previous == NULL || causal_one_slab == NULL ||
       monolithic_window == NULL || causal_window == NULL)
    {
        fprintf(stderr, "ERROR: validation work allocation failed.\n");
        status = EXIT_FAILURE;
        goto cleanup;
    }

    STSB_set_previous_trace_from_nodal_value(
        &monolithic_system.fe,
        monolithic_system.vals.T,
        monolithic_previous);
    STSB_set_previous_trace_from_nodal_value(
        &causal_system.fe,
        causal_system.vals.T,
        causal_previous);

    if(monolis_mpi_get_global_my_rank() == 0){
        char filename[4096];
        snprintf(filename, sizeof(filename),
            "%s/convdiff_stfom_dg0_validation.csv",
            monolithic_system.cond.directory);
        csv = fopen(filename, "w");
        if(csv == NULL){
            fprintf(stderr, "ERROR: cannot open %s\n", filename);
            status = EXIT_FAILURE;
            goto cleanup;
        }
        fprintf(csv,
            "window,time_start,time_end,window_relative_error,"
            "right_trace_relative_error\n");
    }

    for(window = 0; window < monolithic_windows; window++){
        const double time_start =
            (double)(window * slabs) * monolithic_system.vals.dt;
        const double time_end = time_start
            + (double)slabs * monolithic_system.vals.dt;
        double values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double window_relative;
        double trace_relative;
        int slab;
        int node;

        memset(monolithic_window, 0,
            ST_fom_window_vector_length(&common_layout) * sizeof(double));
        memset(causal_window, 0,
            ST_fom_window_vector_length(&common_layout) * sizeof(double));

        STSB_solve_window(
            &monolithic_system,
            &monolithic_time,
            slabs,
            time_start,
            window + 1,
            monolithic_previous,
            monolithic_window);
        STSB_extract_last_right_trace(
            &monolithic_time,
            slabs,
            monolithic_system.fe.total_num_nodes,
            monolithic_window,
            monolithic_trace);

        for(slab = 0; slab < slabs; slab++){
            const double slab_start = time_start
                + (double)slab * causal_system.vals.dt;
            const int global_step = window * slabs + slab + 1;

            memset(causal_one_slab, 0,
                (size_t)causal_system.fe.total_num_nodes * sizeof(double));
            STSB_solve_window(
                &causal_system,
                &causal_time,
                1,
                slab_start,
                global_step,
                causal_previous,
                causal_one_slab);

            if(ST_fom_pack_dg0_slab(
                &common_layout,
                slab,
                causal_one_slab,
                causal_window) != ST_FOM_SUCCESS)
            {
                fprintf(stderr, "ERROR: common dG(0) slab packing failed.\n");
                status = EXIT_FAILURE;
                goto cleanup;
            }

            STSB_extract_last_right_trace(
                &causal_time,
                1,
                causal_system.fe.total_num_nodes,
                causal_one_slab,
                causal_previous);
        }

        if(ST_fom_owned_difference_squares(
            &common_layout,
            monolithic_window,
            causal_window,
            &values[0],
            &values[1],
            &values[2]) != ST_FOM_SUCCESS)
        {
            fprintf(stderr, "ERROR: common dG(0) window comparison failed.\n");
            status = EXIT_FAILURE;
            goto cleanup;
        }

        for(node = 0; node < monolithic_owned; node++){
            const double difference =
                monolithic_trace[node] - causal_previous[node];
            values[3] += difference * difference;
            values[4] += monolithic_trace[node] * monolithic_trace[node];
        }

        monolis_allreduce_R(
            5,
            values,
            MONOLIS_MPI_SUM,
            monolithic_system.monolis_com.comm);

        window_relative = sqrt(fmax(values[0], 0.0)) /
            fmax(sqrt(fmax(values[1], 0.0)), 1.0e-300);
        trace_relative = sqrt(fmax(values[3], 0.0)) /
            fmax(sqrt(fmax(values[4], 0.0)), 1.0e-300);
        if(window_relative > max_window_error){
            max_window_error = window_relative;
        }
        if(trace_relative > max_trace_error){
            max_trace_error = trace_relative;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "convdiff dG(0) ST validation window %d: "
                "window=%e right_trace=%e\n",
                window, window_relative, trace_relative);
            fprintf(csv, "%d,%.17e,%.17e,%.17e,%.17e\n",
                window, time_start, time_end,
                window_relative, trace_relative);
            fflush(csv);
        }

        if(write_vtk){
            write_compare_vtk(
                &monolithic_system,
                window,
                monolithic_trace,
                causal_previous);
        }

        memcpy(monolithic_previous, monolithic_trace,
            (size_t)monolithic_system.fe.total_num_nodes * sizeof(double));
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("\nConvection-diffusion common ST-FOM validation\n");
        printf("  formulation                 : dG(0)\n");
        printf("  monolithic slabs/window     : %d\n", slabs);
        printf("  causal reference            : repeated one-slab dG(0) solves\n");
        printf("  max window relative error   : %.15e\n", max_window_error);
        printf("  max trace relative error    : %.15e\n", max_trace_error);
        printf("  tolerance                   : %.15e\n", tolerance);
        printf("  status                      : %s\n",
            (max_window_error <= tolerance && max_trace_error <= tolerance)
                ? "PASS" : "FAIL");
    }

    if(max_window_error > tolerance || max_trace_error > tolerance){
        status = EXIT_FAILURE;
    }

cleanup:
    if(csv != NULL){
        fclose(csv);
    }
    free(monolithic_previous);
    free(monolithic_trace);
    free(causal_previous);
    free(causal_one_slab);
    free(monolithic_window);
    free(causal_window);
    STSB_finalize_reference_fom(&causal_system, &causal_graph);
    STSB_finalize_reference_fom(&monolithic_system, &monolithic_graph);
    monolis_global_finalize();
    return status;
}
