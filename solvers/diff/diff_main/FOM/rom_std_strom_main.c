/*
 * Complete reference main program for the window-local ST-ROM package.
 *
 * This program exercises the complete workflow
 *
 *   training snapshots
 *      -> window-local ST-POD
 *      -> D_w = Phi_w^T A_w Phi_w
 *      -> L_w = Phi_w^T B_w Phi_{w-1}
 *      -> sequential reduced solve
 *      -> dense all-at-once verification
 *      -> reconstruction and error output.
 *
 * The operator and training data used below are deterministic demonstration
 * data.  In the production application, replace
 *
 *   demo_apply_window_operator(),
 *   demo_generate_training_snapshot(), and
 *   demo_build_online_exact_solution()
 *
 * by adapters to the time-dG FOM, while retaining the main control flow.
 */

#include "rom_std_strom.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "monolis.h"
#include "monolis_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEMO_NUM_WINDOWS       3
#define DEMO_TOTAL_NUM_NODES   4
#define DEMO_DOF               1
#define DEMO_NUM_SLABS         2
#define DEMO_NUM_TIME_BASIS    2
#define DEMO_NUM_SNAPSHOTS     6
#define DEMO_NUM_MODES_MAX     4
#define DEMO_ROM_EPSILON       1.0

typedef struct {
    ROM_STROM_SYSTEM* system;
    double diagonal_base;
    double time_subdiagonal;
} DEMO_WINDOW_OPERATOR_CONTEXT;

static const char* strom_status_string(int status)
{
    switch(status){
    case ROM_STROM_SUCCESS:
        return "success";
    case ROM_STROM_ERR_ARGUMENT:
        return "invalid argument";
    case ROM_STROM_ERR_ALLOCATION:
        return "allocation failure";
    case ROM_STROM_ERR_OPERATOR:
        return "operator application failure";
    case ROM_STROM_ERR_SOLVER:
        return "reduced solver failure";
    case ROM_STROM_ERR_IO:
        return "I/O failure";
    case ROM_STROM_ERR_DIMENSION:
        return "invalid dimension";
    default:
        return "unknown error";
    }
}

static int make_output_directory(const char* directory)
{
    if(directory == NULL || directory[0] == '\0'){
        return 1;
    }

    if(mkdir(directory, 0775) == 0){
        return 0;
    }

    if(errno == EEXIST){
        return 0;
    }

    fprintf(
        stderr,
        "ERROR: cannot create output directory '%s': %s\n",
        directory,
        strerror(errno));
    return 1;
}

static double** allocate_window_vectors(
    const ROM_STROM_SYSTEM* system)
{
    double** vectors;

    if(system == NULL || system->num_windows <= 0){
        return NULL;
    }

    vectors = (double**)calloc(
        (size_t)system->num_windows,
        sizeof(double*));

    if(vectors == NULL){
        return NULL;
    }

    for(int w = 0; w < system->num_windows; w++){
        vectors[w] = (double*)calloc(
            (size_t)system->windows[w].st_dof,
            sizeof(double));

        if(vectors[w] == NULL){
            for(int q = 0; q < w; q++){
                free(vectors[q]);
            }
            free(vectors);
            return NULL;
        }
    }

    return vectors;
}

static void free_window_vectors(
    double*** vectors,
    int num_windows)
{
    if(vectors == NULL || *vectors == NULL){
        return;
    }

    for(int w = 0; w < num_windows; w++){
        free((*vectors)[w]);
    }

    free(*vectors);
    *vectors = NULL;
}

/*
 * Matrix-free demonstration window operator.
 *
 * The time block index is q = slab * num_time_basis + alpha.  For each
 * spatial degree of freedom p,
 *
 *   (A_w x)_{q,p}
 *       = d_{w,q,p} x_{q,p} - beta x_{q-1,p}, q > 0.
 *
 * This lower-triangular form mimics a causal time-dG window operator while
 * remaining independent of a particular finite-element application.
 */
static int demo_apply_window_operator(
    int window_id,
    const double* x,
    double* y,
    void* context)
{
    DEMO_WINDOW_OPERATOR_CONTEXT* operator_context;
    const ROM_STROM_WINDOW* window;
    int num_time_blocks;

    if(x == NULL || y == NULL || context == NULL){
        return 1;
    }

    operator_context = (DEMO_WINDOW_OPERATOR_CONTEXT*)context;

    if(operator_context->system == NULL ||
       window_id < 0 ||
       window_id >= operator_context->system->num_windows){
        return 1;
    }

    window = &operator_context->system->windows[window_id];
    num_time_blocks = window->num_slabs * window->num_time_basis;

    for(int q = 0; q < num_time_blocks; q++){
        for(int p = 0; p < window->spatial_dof; p++){
            const int row = q * window->spatial_dof + p;
            const double diagonal =
                operator_context->diagonal_base
                + 0.20 * (double)window_id
                + 0.05 * (double)q
                + 0.01 * (double)p;

            double value = diagonal * x[row];

            if(q > 0){
                const int previous_row =
                    (q - 1) * window->spatial_dof + p;

                value -= operator_context->time_subdiagonal
                       * x[previous_row];
            }

            y[row] = value;
        }
    }

    return 0;
}

static double demo_pattern_value(
    int pattern_id,
    int window_id,
    int row,
    int st_dof)
{
    const double xi =
        ((double)row + 1.0) / ((double)st_dof + 1.0);

    switch(pattern_id){
    case 0:
        return 1.0 + 0.05 * (double)window_id;
    case 1:
        return xi;
    case 2:
        return sin(M_PI * xi + 0.17 * (double)window_id);
    default:
        return cos(2.0 * M_PI * xi - 0.11 * (double)window_id);
    }
}

/*
 * Generates snapshots whose rank is at most DEMO_NUM_MODES_MAX.
 * Snapshot zero is later reused as the exact online trajectory, so it lies
 * in the retained trial space when all demonstration modes are retained.
 */
static void demo_generate_training_snapshot(
    const ROM_STROM_WINDOW* window,
    int snapshot_id,
    double* snapshot)
{
    for(int row = 0; row < window->st_dof; row++){
        double value = 0.0;

        for(int k = 0; k < DEMO_NUM_MODES_MAX; k++){
            const double coefficient =
                cos(0.37 * (double)(snapshot_id + 1) * (double)(k + 1))
                + 0.20 * sin(0.23 * (double)(snapshot_id + k + 2));

            value += coefficient * demo_pattern_value(
                k,
                window->window_id,
                row,
                window->st_dof);
        }

        snapshot[row] = value;
    }
}

static void demo_build_online_exact_solution(
    const ROM_STROM_WINDOW* window,
    double* exact_solution)
{
    demo_generate_training_snapshot(window, 0, exact_solution);
}

static double relative_l2_error(
    const double* reference,
    const double* approximation,
    int n)
{
    double numerator = 0.0;
    double denominator = 0.0;

    for(int i = 0; i < n; i++){
        const double difference = reference[i] - approximation[i];
        numerator += difference * difference;
        denominator += reference[i] * reference[i];
    }

    if(denominator <= 0.0){
        return sqrt(numerator);
    }

    return sqrt(numerator / denominator);
}

static int write_window_solution_csv(
    const char* directory,
    int window_id,
    const double* exact_solution,
    const double* rom_solution,
    int n)
{
    char filename[1024];
    FILE* fp;

    if(snprintf(
        filename,
        sizeof(filename),
        "%s/window_%04d_solution.csv",
        directory,
        window_id) >= (int)sizeof(filename)){
        return 1;
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        return 1;
    }

    fprintf(fp, "row,exact,rom,error\n");
    for(int i = 0; i < n; i++){
        fprintf(
            fp,
            "%d,%.17e,%.17e,%.17e\n",
            i,
            exact_solution[i],
            rom_solution[i],
            rom_solution[i] - exact_solution[i]);
    }

    fclose(fp);
    return 0;
}

static int write_pod_data(
    const char* directory,
    const ROM_STROM_WINDOW* window)
{
    char filename[1024];

    if(snprintf(
        filename,
        sizeof(filename),
        "%s/window_%04d_pod.bin",
        directory,
        window->window_id) >= (int)sizeof(filename)){
        return ROM_STROM_ERR_IO;
    }

    return ROM_std_strom_write_pod_modes_binary(window, filename);
}

static void print_usage(const char* program_name)
{
    printf(
        "Usage: %s [output_directory]\n"
        "\n"
        "Runs the complete window-local ST-ROM reference workflow.\n"
        "The default output directory is 'strom_output'.\n",
        program_name);
}

int main(int argc, char* argv[])
{
    const char* output_directory = "strom_output";
    const double theta_plus[DEMO_NUM_TIME_BASIS] = {0.0, 1.0};
    const double jump_left[DEMO_NUM_TIME_BASIS] = {1.0, 0.0};

    ROM_STROM_SYSTEM system;
    HLPOD_VALUES* pod_values = NULL;
    HLPOD_MAT* pod_matrices = NULL;
    ROM_STROM_DG_COUPLING_CONTEXT coupling_context;
    DEMO_WINDOW_OPERATOR_CONTEXT operator_context;

    double** exact_solution = NULL;
    double** full_rhs = NULL;
    double** reconstructed_solution = NULL;
    double* snapshot_work = NULL;
    double* operator_work = NULL;
    double* coupling_work = NULL;
    double* sequential_coefficients = NULL;

    int status = ROM_STROM_SUCCESS;
    int initialized_system = 0;
    int exit_code = EXIT_FAILURE;
    int my_rank = 0;

    memset(&system, 0, sizeof(system));
    memset(&coupling_context, 0, sizeof(coupling_context));
    memset(&operator_context, 0, sizeof(operator_context));

    if(argc > 2){
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    if(argc == 2){
        if(strcmp(argv[1], "--help") == 0 ||
           strcmp(argv[1], "-h") == 0){
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        }
        output_directory = argv[1];
    }

    monolis_global_initialize();
    my_rank = monolis_mpi_get_global_my_rank();

    if(my_rank == 0 && make_output_directory(output_directory) != 0){
        goto cleanup;
    }

    pod_values = (HLPOD_VALUES*)calloc(
        DEMO_NUM_WINDOWS,
        sizeof(HLPOD_VALUES));

    pod_matrices = (HLPOD_MAT*)calloc(
        DEMO_NUM_WINDOWS,
        sizeof(HLPOD_MAT));

    if(pod_values == NULL || pod_matrices == NULL){
        fprintf(stderr, "ERROR: POD structure allocation failed.\n");
        goto cleanup;
    }

    status = ROM_std_strom_system_initialize(
        &system,
        DEMO_NUM_WINDOWS);

    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: system initialization failed: %s\n",
            strom_status_string(status));
        goto cleanup;
    }
    initialized_system = 1;

    for(int w = 0; w < DEMO_NUM_WINDOWS; w++){
        status = ROM_std_strom_window_initialize(
            &system.windows[w],
            w,
            &pod_values[w],
            &pod_matrices[w],
            DEMO_TOTAL_NUM_NODES,
            DEMO_DOF,
            DEMO_NUM_SLABS,
            DEMO_NUM_TIME_BASIS,
            DEMO_NUM_SNAPSHOTS,
            DEMO_NUM_MODES_MAX,
            0);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: window %d initialization failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }
    }

    operator_context.system = &system;
    operator_context.diagonal_base = 2.0;
    operator_context.time_subdiagonal = 0.25;

    snapshot_work = (double*)calloc(
        (size_t)system.windows[0].st_dof,
        sizeof(double));

    if(snapshot_work == NULL){
        fprintf(stderr, "ERROR: snapshot work allocation failed.\n");
        goto cleanup;
    }

    /* -------------------------------------------------------------- */
    /* Offline stage: snapshot storage, POD and diagonal projections. */
    /* -------------------------------------------------------------- */

    for(int w = 0; w < system.num_windows; w++){
        ROM_STROM_WINDOW* window = &system.windows[w];

        for(int m = 0; m < window->num_snapshots; m++){
            demo_generate_training_snapshot(window, m, snapshot_work);

            status = ROM_std_strom_store_snapshot(
                window,
                snapshot_work,
                m,
                0);

            if(status != ROM_STROM_SUCCESS){
                fprintf(
                    stderr,
                    "ERROR: storing snapshot (%d,%d) failed: %s\n",
                    w,
                    m,
                    strom_status_string(status));
                goto cleanup;
            }
        }

        status = ROM_std_strom_set_pod_modes(
            window,
            DEMO_ROM_EPSILON);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: POD for window %d failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }

        status = ROM_std_strom_calc_reduced_diag(
            window,
            demo_apply_window_operator,
            &operator_context);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: D_%d construction failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }

        if(my_rank == 0){
            status = write_pod_data(output_directory, window);
            if(status != ROM_STROM_SUCCESS){
                fprintf(
                    stderr,
                    "ERROR: POD output for window %d failed.\n",
                    w);
                goto cleanup;
            }
        }
    }

    /* -------------------------------------------------------------- */
    /* Time-dG transition data and L_w projections.                    */
    /* -------------------------------------------------------------- */

    coupling_context.num_windows = system.num_windows;
    coupling_context.transition =
        (ROM_STROM_DG_TRANSITION*)calloc(
            (size_t)system.num_windows,
            sizeof(ROM_STROM_DG_TRANSITION));

    if(coupling_context.transition == NULL){
        fprintf(stderr, "ERROR: transition allocation failed.\n");
        goto cleanup;
    }

    for(int w = 1; w < system.num_windows; w++){
        ROM_STROM_DG_TRANSITION* transition =
            &coupling_context.transition[w];

        transition->spatial_dof = system.windows[w].spatial_dof;
        transition->prev_num_slabs = system.windows[w - 1].num_slabs;
        transition->prev_num_time_basis =
            system.windows[w - 1].num_time_basis;
        transition->cur_num_slabs = system.windows[w].num_slabs;
        transition->cur_num_time_basis =
            system.windows[w].num_time_basis;
        transition->theta_plus_prev = theta_plus;
        transition->jump_left_cur = jump_left;

        /* Identity spatial trace operator in this demonstration. */
        transition->apply_spatial_trace_operator = NULL;
        transition->spatial_operator_context = NULL;

        status = ROM_std_strom_calc_reduced_coupling(
            &system.windows[w],
            &system.windows[w - 1],
            ROM_std_strom_apply_dg_coupling,
            &coupling_context);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: L_%d construction failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }
    }

    status = ROM_std_strom_update_mode_offsets(&system);
    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: reduced offsets failed: %s\n",
            strom_status_string(status));
        goto cleanup;
    }

    exact_solution = allocate_window_vectors(&system);
    full_rhs = allocate_window_vectors(&system);
    reconstructed_solution = allocate_window_vectors(&system);

    if(exact_solution == NULL || full_rhs == NULL ||
       reconstructed_solution == NULL){
        fprintf(stderr, "ERROR: window-vector allocation failed.\n");
        goto cleanup;
    }

    operator_work = (double*)calloc(
        (size_t)system.windows[0].st_dof,
        sizeof(double));
    coupling_work = (double*)calloc(
        (size_t)system.windows[0].st_dof,
        sizeof(double));

    if(operator_work == NULL || coupling_work == NULL){
        fprintf(stderr, "ERROR: operator work allocation failed.\n");
        goto cleanup;
    }

    /* -------------------------------------------------------------- */
    /* Build a consistent online FOM right-hand side                   */
    /* A_w U_w - B_w U_{w-1} = f_w.                                   */
    /* -------------------------------------------------------------- */

    for(int w = 0; w < system.num_windows; w++){
        ROM_STROM_WINDOW* window = &system.windows[w];
        double* effective_rhs;

        demo_build_online_exact_solution(
            window,
            exact_solution[w]);

        memset(
            operator_work,
            0,
            (size_t)window->st_dof * sizeof(double));

        if(demo_apply_window_operator(
            w,
            exact_solution[w],
            operator_work,
            &operator_context) != 0){
            fprintf(stderr, "ERROR: online A_%d application failed.\n", w);
            goto cleanup;
        }

        memcpy(
            full_rhs[w],
            operator_work,
            (size_t)window->st_dof * sizeof(double));

        if(w > 0){
            memset(
                coupling_work,
                0,
                (size_t)window->st_dof * sizeof(double));

            if(ROM_std_strom_apply_dg_coupling(
                w,
                exact_solution[w - 1],
                coupling_work,
                &coupling_context) != ROM_STROM_SUCCESS){
                fprintf(stderr, "ERROR: online B_%d application failed.\n", w);
                goto cleanup;
            }

            for(int i = 0; i < window->st_dof; i++){
                full_rhs[w][i] -= coupling_work[i];
            }
        }

        effective_rhs = (double*)calloc(
            (size_t)window->st_dof,
            sizeof(double));

        if(effective_rhs == NULL){
            fprintf(stderr, "ERROR: effective RHS allocation failed.\n");
            goto cleanup;
        }

        status = ROM_std_strom_build_effective_rhs(
            window,
            w > 0 ? &system.windows[w - 1] : NULL,
            full_rhs[w],
            effective_rhs,
            demo_apply_window_operator,
            &operator_context,
            ROM_std_strom_apply_dg_coupling,
            &coupling_context);

        if(status == ROM_STROM_SUCCESS){
            status = ROM_std_strom_project_rhs(
                window,
                effective_rhs);
        }

        free(effective_rhs);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: reduced RHS for window %d failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }
    }

    /* -------------------------------------------------------------- */
    /* Sequential reduced solve and reconstruction.                    */
    /* -------------------------------------------------------------- */

    status = ROM_std_strom_solve_sequential(
        &system,
        NULL,
        NULL);

    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: sequential ROM solve failed: %s\n",
            strom_status_string(status));
        goto cleanup;
    }

    sequential_coefficients = (double*)calloc(
        (size_t)system.total_reduced_dof,
        sizeof(double));

    if(sequential_coefficients == NULL){
        fprintf(stderr, "ERROR: coefficient allocation failed.\n");
        goto cleanup;
    }

    for(int w = 0; w < system.num_windows; w++){
        ROM_STROM_WINDOW* window = &system.windows[w];
        const int offset = system.mode_offset[w];

        for(int i = 0; i < window->num_modes; i++){
            sequential_coefficients[offset + i] =
                window->reduced_coef[i];
        }

        status = ROM_std_strom_reconstruct_window(
            window,
            reconstructed_solution[w]);

        if(status != ROM_STROM_SUCCESS){
            fprintf(
                stderr,
                "ERROR: reconstruction of window %d failed: %s\n",
                w,
                strom_status_string(status));
            goto cleanup;
        }
    }

    /* -------------------------------------------------------------- */
    /* Dense all-at-once verification.                                 */
    /* -------------------------------------------------------------- */

    status = ROM_std_strom_solve_all_at_once_dense(
        &system,
        NULL,
        NULL);

    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: all-at-once ROM solve failed: %s\n",
            strom_status_string(status));
        goto cleanup;
    }

    if(my_rank == 0){
        double coefficient_difference_squared = 0.0;
        double coefficient_reference_squared = 0.0;

        printf("\nWindow-local ST-ROM reference run\n");
        printf("---------------------------------\n");
        printf("num_windows       = %d\n", system.num_windows);
        printf("st_dof/window      = %d\n", system.windows[0].st_dof);
        printf("total reduced dof  = %d\n", system.total_reduced_dof);

        for(int w = 0; w < system.num_windows; w++){
            const ROM_STROM_WINDOW* window = &system.windows[w];
            const double error = relative_l2_error(
                exact_solution[w],
                reconstructed_solution[w],
                window->st_dof);

            printf(
                "window %d: modes = %d, relative reconstruction error = %.6e\n",
                w,
                window->num_modes,
                error);

            if(write_window_solution_csv(
                output_directory,
                w,
                exact_solution[w],
                reconstructed_solution[w],
                window->st_dof) != 0){
                fprintf(
                    stderr,
                    "ERROR: failed to write window %d solution CSV.\n",
                    w);
                goto cleanup;
            }

            for(int i = 0; i < window->num_modes; i++){
                const int global_index = system.mode_offset[w] + i;
                const double difference =
                    window->reduced_coef[i]
                    - sequential_coefficients[global_index];

                coefficient_difference_squared +=
                    difference * difference;
                coefficient_reference_squared +=
                    sequential_coefficients[global_index]
                    * sequential_coefficients[global_index];
            }
        }

        printf(
            "sequential/all-at-once coefficient difference = %.6e\n",
            coefficient_reference_squared > 0.0
                ? sqrt(
                    coefficient_difference_squared
                    / coefficient_reference_squared)
                : sqrt(coefficient_difference_squared));
        printf("output directory   = %s\n\n", output_directory);
    }

    exit_code = EXIT_SUCCESS;

cleanup:
    free(snapshot_work);
    free(operator_work);
    free(coupling_work);
    free(sequential_coefficients);

    free_window_vectors(&exact_solution, DEMO_NUM_WINDOWS);
    free_window_vectors(&full_rhs, DEMO_NUM_WINDOWS);
    free_window_vectors(&reconstructed_solution, DEMO_NUM_WINDOWS);

    free(coupling_context.transition);
    coupling_context.transition = NULL;

    if(initialized_system){
        ROM_std_strom_system_finalize(&system);
    }

    free(pod_values);
    free(pod_matrices);

    monolis_global_finalize();
    return exit_code;
}
