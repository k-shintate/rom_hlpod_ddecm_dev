#include "core_ROM.h"
#include "core_FOM_ST.h"
#include "rom_std_strom.h"
#include "rom_std_strom_fom_adapter.h"
#include "rom_std_strom_monolis.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef ST_ROM_DEGREE
#define ST_ROM_DEGREE 1
#endif

#ifndef ST_ROM_WINDOW_SLABS
#define ST_ROM_WINDOW_SLABS 4
#endif

/*
 * Runtime mode is selected by environment variables so that the existing
 * BBFE command-line directory parser remains unchanged.
 *
 * ST_ROM_MODE=fom      : ordinary serial windowed ST-FOM
 * ST_ROM_MODE=collect  : ST-FOM + write homogeneous lifted snapshots
 * ST_ROM_MODE=offline  : pool every window snapshot and construct one
 *                        all-window shared POD basis
 * ST_ROM_MODE=online   : read bases and solve the serial window ST-ROM
 *
 * Additional variables:
 * ST_ROM_SNAPSHOT_ID       integer used in collect mode
 * ST_ROM_NUM_SNAPSHOTS     number of training cases in offline mode
 * ST_ROM_SNAPSHOT_DIR      default: <case>/strom_snapshots
 * ST_ROM_BASIS_DIR         default: <case>/strom_basis
 * ST_ROM_EPSILON           default: 0.999999
 * ST_ROM_MAX_MODES         offline default: num_windows*num_snapshots
 *                          online default capacity: 256
 * ST_ROM_VALIDATE          online validation, default: 1
 * ST_ROM_VALIDATION_SNAPSHOT_ID  validation snapshot, default: 0
 * ST_ROM_VALIDATION_CSV    default: <case>/strom_validation.csv
 * ST_ROM_SOLVER            all_at_once (default), monolis, verify_monolis,
 *                          sequential, or verify
 * ST_ROM_OPERATOR_REUSE     1 (default) when A_w and B_w are window-invariant
 * ST_ROM_DIRECT_COUPLING    1 (default): preintegrated direct B action
 * ST_ROM_VERIFY_DIRECT_B    0 (default): compare direct and legacy B once
 * ST_ROM_REBUILD_OPERATORS  0 (default): rebuild instead of loading D/L
 * ST_ROM_WRITE_OUTPUT       1 (default): write VTK/temperature/L2 output
 * ST_ROM_MONOLIS_EPSILON    reduced BiCGSTAB tolerance, default: 1.0e-12
 * ST_ROM_MONOLIS_MAX_ITER   reduced BiCGSTAB maximum iterations, default: 10000
 * ST_ROM_MONOLIS_VERIFY_TOL dense/MONOLIS coefficient tolerance, default: 1.0e-10
 * ST_ROM_WRITE_REDUCED_GRAPH 1 (default): write the retained window metagraph
 * ST_ROM_REDUCED_GRAPH_FILE optional graph output path
 *
 * Snapshot format version 2 stores q_w = U_w - U_D,w.  Old version-1
 * full-state snapshots are deliberately rejected and must be recollected.
 */

typedef enum {
    ST_RUN_FOM = 0,
    ST_RUN_COLLECT,
    ST_RUN_OFFLINE,
    ST_RUN_ONLINE
} ST_RUN_MODE;

typedef enum {
    ST_ROM_SOLVER_ALL_AT_ONCE = 0,
    ST_ROM_SOLVER_MONOLIS,
    ST_ROM_SOLVER_VERIFY_MONOLIS,
    ST_ROM_SOLVER_SEQUENTIAL,
    ST_ROM_SOLVER_VERIFY
} ST_ROM_ONLINE_SOLVER;

typedef struct {
    ST_RUN_MODE mode;
    ST_ROM_ONLINE_SOLVER online_solver;
    int snapshot_id;
    int num_snapshots;
    int max_modes;
    int max_modes_explicit;
    int validate_online;
    int validation_snapshot_id;
    int operator_reuse;
    int direct_coupling;
    int verify_direct_b;
    int rebuild_operators;
    int write_output;
    int monolis_max_iter;
    int write_reduced_graph;
    double epsilon;
    double monolis_tolerance;
    double monolis_verify_tolerance;
    char snapshot_dir[4096];
    char basis_dir[4096];
    char validation_csv[4096];
    char reduced_graph_file[4096];
} ST_ROM_OPTIONS;

static void fail_message(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static int checked_product3_int(int a, int b, int c)
{
    const long long value = (long long)a * b * c;

    if(a <= 0 || b <= 0 || c <= 0 || value > INT_MAX){
        fail_message("invalid or overflowing ST dimension");
    }

    return (int)value;
}

static void require_serial_runtime(void)
{
    int comm_size = monolis_mpi_get_global_comm_size();

    if(comm_size != 1){
        fail_message(
            "this implementation is serial; run with mpirun -np 1");
    }
}

static int make_directory(const char* directory)
{
    if(directory == NULL || directory[0] == '\0'){
        return 1;
    }

    if(mkdir(directory, 0775) == 0 || errno == EEXIST){
        return 0;
    }

    fprintf(
        stderr,
        "ERROR: mkdir(%s): %s\n",
        directory,
        strerror(errno));

    return 1;
}

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

    if(errno != 0 ||
       end == text ||
       *end != '\0' ||
       value < 0 ||
       value > INT_MAX)
    {
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

    if(errno != 0 ||
       end == text ||
       *end != '\0' ||
       !isfinite(value))
    {
        fprintf(stderr, "ERROR: invalid real %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return value;
}

static void join_default_directory(
    char* output,
    size_t output_size,
    const char* case_directory,
    const char* suffix,
    const char* environment_name)
{
    const char* override = getenv(environment_name);

    if(override != NULL && override[0] != '\0'){
        if(snprintf(output, output_size, "%s", override)
           >= (int)output_size)
        {
            fail_message("ROM directory path is too long");
        }
        return;
    }

    if(snprintf(
        output,
        output_size,
        "%s/%s",
        case_directory,
        suffix) >= (int)output_size)
    {
        fail_message("ROM directory path is too long");
    }
}

static ST_RUN_MODE read_mode(void)
{
    const char* mode = getenv("ST_ROM_MODE");

    if(mode == NULL || strcmp(mode, "fom") == 0){
        return ST_RUN_FOM;
    }
    if(strcmp(mode, "collect") == 0){
        return ST_RUN_COLLECT;
    }
    if(strcmp(mode, "offline") == 0){
        return ST_RUN_OFFLINE;
    }
    if(strcmp(mode, "online") == 0){
        return ST_RUN_ONLINE;
    }

    fprintf(stderr, "ERROR: unknown ST_ROM_MODE=%s\n", mode);
    exit(EXIT_FAILURE);
}

static const char* mode_name(ST_RUN_MODE mode)
{
    switch(mode){
    case ST_RUN_FOM:
        return "fom";
    case ST_RUN_COLLECT:
        return "collect";
    case ST_RUN_OFFLINE:
        return "offline";
    case ST_RUN_ONLINE:
        return "online";
    default:
        return "unknown";
    }
}

static ST_ROM_ONLINE_SOLVER read_online_solver(void)
{
    const char* solver = getenv("ST_ROM_SOLVER");

    if(solver == NULL || solver[0] == '\0' ||
       strcmp(solver, "all_at_once") == 0 ||
       strcmp(solver, "all-at-once") == 0)
    {
        return ST_ROM_SOLVER_ALL_AT_ONCE;
    }

    if(strcmp(solver, "monolis") == 0 ||
       strcmp(solver, "all_at_once_monolis") == 0 ||
       strcmp(solver, "all-at-once-monolis") == 0)
    {
        return ST_ROM_SOLVER_MONOLIS;
    }

    if(strcmp(solver, "verify_monolis") == 0 ||
       strcmp(solver, "verify-monolis") == 0)
    {
        return ST_ROM_SOLVER_VERIFY_MONOLIS;
    }

    if(strcmp(solver, "sequential") == 0){
        return ST_ROM_SOLVER_SEQUENTIAL;
    }

    if(strcmp(solver, "verify") == 0){
        return ST_ROM_SOLVER_VERIFY;
    }

    fprintf(stderr, "ERROR: unknown ST_ROM_SOLVER=%s\n", solver);
    exit(EXIT_FAILURE);
}

static const char* online_solver_name(ST_ROM_ONLINE_SOLVER solver)
{
    switch(solver){
    case ST_ROM_SOLVER_ALL_AT_ONCE:
        return "all_at_once_dense";
    case ST_ROM_SOLVER_MONOLIS:
        return "all_at_once_monolis_bcsr";
    case ST_ROM_SOLVER_VERIFY_MONOLIS:
        return "verify_dense_vs_monolis_bcsr";
    case ST_ROM_SOLVER_SEQUENTIAL:
        return "sequential";
    case ST_ROM_SOLVER_VERIFY:
        return "verify_sequential_vs_all_at_once";
    default:
        return "unknown";
    }
}

static void read_rom_options(
    ST_ROM_OPTIONS* options,
    const char* case_directory)
{
    memset(options, 0, sizeof(*options));

    options->mode = read_mode();
    options->online_solver = read_online_solver();
    options->snapshot_id = env_int("ST_ROM_SNAPSHOT_ID", 0);
    options->num_snapshots = env_int("ST_ROM_NUM_SNAPSHOTS", 1);
    options->epsilon = env_double("ST_ROM_EPSILON", 0.999999);
    options->validate_online = env_int(
        "ST_ROM_VALIDATE",
        options->mode == ST_RUN_ONLINE ? 1 : 0);
    options->validation_snapshot_id = env_int(
        "ST_ROM_VALIDATION_SNAPSHOT_ID",
        0);
    options->operator_reuse = env_int("ST_ROM_OPERATOR_REUSE", 1);
    options->direct_coupling = env_int("ST_ROM_DIRECT_COUPLING", 1);
    options->verify_direct_b = env_int("ST_ROM_VERIFY_DIRECT_B", 0);
    options->rebuild_operators = env_int("ST_ROM_REBUILD_OPERATORS", 0);
    options->write_output = env_int("ST_ROM_WRITE_OUTPUT", 1);
    options->monolis_max_iter = env_int(
        "ST_ROM_MONOLIS_MAX_ITER",
        10000);
    options->monolis_tolerance = env_double(
        "ST_ROM_MONOLIS_EPSILON",
        1.0e-12);
    options->monolis_verify_tolerance = env_double(
        "ST_ROM_MONOLIS_VERIFY_TOL",
        1.0e-10);
    options->write_reduced_graph = env_int(
        "ST_ROM_WRITE_REDUCED_GRAPH",
        1);

    if(options->validate_online != 0 &&
       options->validate_online != 1)
    {
        fail_message("ST_ROM_VALIDATE must be 0 or 1");
    }

    if((options->operator_reuse != 0 && options->operator_reuse != 1) ||
       (options->direct_coupling != 0 && options->direct_coupling != 1) ||
       (options->verify_direct_b != 0 && options->verify_direct_b != 1) ||
       (options->rebuild_operators != 0 && options->rebuild_operators != 1) ||
       (options->write_output != 0 && options->write_output != 1) ||
       (options->write_reduced_graph != 0 && options->write_reduced_graph != 1))
    {
        fail_message("ST-ROM boolean options must be 0 or 1");
    }

    options->max_modes_explicit =
        getenv("ST_ROM_MAX_MODES") != NULL;

    if(options->max_modes_explicit){
        options->max_modes = env_int("ST_ROM_MAX_MODES", 1);
    }
    else if(options->mode == ST_RUN_ONLINE){
        /* Capacity used while loading the shared basis file. */
        options->max_modes = 256;
    }
    else{
        /* Offline chooses num_windows*num_snapshots after partitioning. */
        options->max_modes = 0;
    }

    if(options->num_snapshots <= 0 ||
       (options->max_modes_explicit && options->max_modes <= 0) ||
       options->epsilon < 0.0 ||
       options->monolis_max_iter <= 0 ||
       options->monolis_tolerance <= 0.0 ||
       options->monolis_verify_tolerance <= 0.0)
    {
        fail_message("invalid ST-ROM options");
    }

    join_default_directory(
        options->snapshot_dir,
        sizeof(options->snapshot_dir),
        case_directory,
        "strom_snapshots",
        "ST_ROM_SNAPSHOT_DIR");

    join_default_directory(
        options->basis_dir,
        sizeof(options->basis_dir),
        case_directory,
        "strom_basis",
        "ST_ROM_BASIS_DIR");

    join_default_directory(
        options->validation_csv,
        sizeof(options->validation_csv),
        case_directory,
        "strom_validation.csv",
        "ST_ROM_VALIDATION_CSV");

    {
        const char* graph_override = getenv("ST_ROM_REDUCED_GRAPH_FILE");

        if(graph_override != NULL && graph_override[0] != '\0'){
            if(snprintf(
                options->reduced_graph_file,
                sizeof(options->reduced_graph_file),
                "%s",
                graph_override) >= (int)sizeof(options->reduced_graph_file))
            {
                fail_message("reduced graph path is too long");
            }
        }
        else if(snprintf(
            options->reduced_graph_file,
            sizeof(options->reduced_graph_file),
            "%s/all_windows_shared_reduced_metagraph.dat",
            options->basis_dir) >= (int)sizeof(options->reduced_graph_file))
        {
            fail_message("reduced graph path is too long");
        }
    }
}

static void initialize_output(FE_SYSTEM* sys, int my_rank)
{
    FILE* fp = NULL;

    fp = BBFE_sys_write_fopen(
        fp,
        "l2_error.txt",
        sys->cond.directory);
    fclose(fp);

    if(my_rank == 0){
        fp = ROM_BB_write_fopen(
            fp,
            "calctime/time_fem.txt",
            sys->cond.directory);
        fclose(fp);

        ROM_std_hlpod_write_solver_prm_fopen(
            "fem_solver_prm",
            sys->cond.directory);
    }
}

static void validate_time_partition(
    const FE_SYSTEM* sys,
    int num_window_slabs,
    int* num_windows)
{
    int total_num_steps;
    double reconstructed_finish_time;
    double time_tolerance;

    if(sys->vals.dt <= 0.0 ||
       sys->vals.finish_time <= 0.0 ||
       sys->vals.output_interval <= 0)
    {
        fail_message("invalid dt, finish_time, or output_interval");
    }

    total_num_steps =
        (int)llround(sys->vals.finish_time / sys->vals.dt);

    reconstructed_finish_time =
        (double)total_num_steps * sys->vals.dt;

    time_tolerance =
        1.0e-12 * fmax(1.0, fabs(sys->vals.finish_time));

    if(fabs(reconstructed_finish_time - sys->vals.finish_time)
       > time_tolerance)
    {
        fail_message("finish_time / dt must be an integer");
    }

    if(total_num_steps % num_window_slabs != 0){
        fail_message(
            "total_num_steps must be divisible by ST_ROM_WINDOW_SLABS");
    }

    *num_windows = total_num_steps / num_window_slabs;
}

static void initialize_fom(
    FE_SYSTEM* sys,
    ST_TIME_ELEMENT* te,
    int argc,
    char** argv,
    int* num_windows)
{
    const int degree = ST_ROM_DEGREE;
    const int num_window_slabs = ST_ROM_WINDOW_SLABS;

    sys->cond.directory = BBFE_convdiff_get_directory_name(
        argc,
        argv,
        CODENAME);

    read_calc_conditions(&sys->vals, sys->cond.directory);

    BBFE_convdiff_pre(
        &sys->fe,
        &sys->basis,
        &sys->bc,
        &sys->monolis,
        &sys->monolis_com,
        argc,
        argv,
        sys->cond.directory,
        sys->vals.num_ip_each_axis,
        true);

    memory_allocation_nodal_values(
        &sys->vals,
        sys->fe.total_num_nodes);

    manusol_set_init_value(
        &sys->fe,
        sys->vals.T);

    BBFE_elemmat_set_Jacobi_mat(
        &sys->fe,
        &sys->basis);

    BBFE_elemmat_set_shapefunc_derivative(
        &sys->fe,
        &sys->basis);

    ST_TimeElement_init_dG(te, degree);

    validate_time_partition(
        sys,
        num_window_slabs,
        num_windows);

    /*
     * BBFE_convdiff_pre() constructed a spatial MONOLIS matrix.  Replace
     * only the matrix by the serial space-time window graph.  The existing
     * MONOLIS_COM is retained; with one MPI rank it has no halo exchange.
     */
    monolis_finalize(&sys->monolis);
    monolis_initialize(&sys->monolis);
    monolis_initialize(&sys->monolis0);

    set_nonzero_pattern_space_time(
        &sys->monolis0,
        "graph.dat",
        sys->cond.directory,
        te,
        num_window_slabs);

    monolis_clear_mat_value_rhs_R(&sys->monolis0);

    set_element_mat_ST_dG_window(
        &sys->monolis0,
        &sys->fe,
        &sys->basis,
        &sys->vals,
        te,
        num_window_slabs);

    monolis_copy_mat_R(
        &sys->monolis0,
        &sys->monolis);
}

static void finalize_fom(FE_SYSTEM* sys)
{
    BBFE_convdiff_finalize(
        &sys->fe,
        &sys->basis,
        &sys->bc);

    monolis_finalize(&sys->monolis);
    monolis_finalize(&sys->monolis0);
}

static void initialize_strom_system(
    ROM_STROM_SYSTEM* system,
    HLPOD_VALUES** values_out,
    HLPOD_MAT** matrices_out,
    const FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    int num_snapshots,
    int max_modes)
{
    HLPOD_VALUES* values;
    HLPOD_MAT* matrices;
    const int expected_st_dof = checked_product3_int(
        sys->fe.total_num_nodes,
        num_window_slabs,
        te->n_dof);
    int status;

    status = ROM_std_strom_system_initialize(system, num_windows);
    if(status != ROM_STROM_SUCCESS){
        fail_message("ROM_std_strom_system_initialize failed");
    }

    values = (HLPOD_VALUES*)calloc(
        (size_t)num_windows,
        sizeof(HLPOD_VALUES));

    matrices = (HLPOD_MAT*)calloc(
        (size_t)num_windows,
        sizeof(HLPOD_MAT));

    if(values == NULL || matrices == NULL){
        free(values);
        free(matrices);
        ROM_std_strom_system_finalize(system);
        fail_message("ST-ROM HLPOD allocation failed");
    }

    for(int w = 0; w < num_windows; w++){
        status = ROM_std_strom_window_initialize(
            &system->windows[w],
            w,
            &values[w],
            &matrices[w],
            sys->fe.total_num_nodes,
            1,
            num_window_slabs,
            te->n_dof,
            num_snapshots,
            max_modes,
            0);

        if(status != ROM_STROM_SUCCESS){
            fail_message("ROM_std_strom_window_initialize failed");
        }

        if(system->windows[w].st_dof != expected_st_dof){
            fail_message("ST-ROM/FOM vector dimensions differ");
        }
    }

    *values_out = values;
    *matrices_out = matrices;
}

static void finalize_strom_system(
    ROM_STROM_SYSTEM* system,
    HLPOD_VALUES** values,
    HLPOD_MAT** matrices)
{
    ROM_std_strom_system_finalize(system);

    free(*values);
    free(*matrices);

    *values = NULL;
    *matrices = NULL;
}

static void run_fom_or_collect(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    const ST_ROM_OPTIONS* options,
    int write_standard_output)
{
    const int nx = sys->fe.total_num_nodes;
    const int nt = te->n_dof;
    const int previous_state_dof = checked_product3_int(nx, nt, 1);
    const int st_dof = checked_product3_int(nx, nt, num_window_slabs);

    ROM_STROM_FOM_CONTEXT collect_context;
    int collect_context_initialized = 0;

    double* T_window = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    double* T_prev_st = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    double* T_last_st = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    double* dirichlet_lift = NULL;
    double* homogeneous_snapshot = NULL;
    int file_num = 0;

    memset(&collect_context, 0, sizeof(collect_context));

    if(T_window == NULL || T_prev_st == NULL || T_last_st == NULL){
        free(T_window);
        free(T_prev_st);
        free(T_last_st);
        fail_message("FOM vector allocation failed");
    }

    ST_TimeElement_set_prev_trace_from_nodal_value(
        te,
        nx,
        sys->vals.T,
        T_prev_st);

    if(options->mode == ST_RUN_COLLECT){
        int status;

        if(make_directory(options->snapshot_dir) != 0){
            fail_message("cannot create snapshot directory");
        }

        status = ROM_std_strom_fom_context_initialize(
            &collect_context,
            sys,
            te,
            num_window_slabs,
            num_windows);

        if(status != ROM_STROM_SUCCESS){
            fail_message("collect FOM adapter initialization failed");
        }

        collect_context_initialized = 1;

        dirichlet_lift = (double*)calloc(
            (size_t)st_dof,
            sizeof(double));

        homogeneous_snapshot = (double*)calloc(
            (size_t)st_dof,
            sizeof(double));

        if(dirichlet_lift == NULL || homogeneous_snapshot == NULL){
            fail_message("homogeneous snapshot allocation failed");
        }
    }

    for(int w = 0; w < num_windows; w++){
        const double t_window_start =
            (double)(w * num_window_slabs) * sys->vals.dt;

        memset(T_window, 0, (size_t)st_dof * sizeof(double));

        solver_fom_space_time(
            sys,
            te,
            num_window_slabs,
            t_window_start,
            w + 1,
            T_prev_st,
            T_window,
            T_last_st);

        if(options->mode == ST_RUN_COLLECT){
            int status = ROM_std_strom_fom_build_dirichlet_lift(
                &collect_context,
                w,
                dirichlet_lift);

            if(status != ROM_STROM_SUCCESS){
                fail_message("Dirichlet lift construction failed in collect mode");
            }

            for(int i = 0; i < st_dof; i++){
                homogeneous_snapshot[i] =
                    T_window[i] - dirichlet_lift[i];
            }

            status = ROM_std_strom_fom_write_snapshot(
                options->snapshot_dir,
                options->snapshot_id,
                w,
                homogeneous_snapshot,
                st_dof);

            if(status != ROM_STROM_SUCCESS){
                fail_message("homogeneous snapshot write failed");
            }
        }

        if(write_standard_output){
            for(int slab = 0; slab < num_window_slabs; slab++){
                const int global_step =
                    w * num_window_slabs + slab + 1;

                const double time =
                    (double)global_step * sys->vals.dt;

                if(global_step % sys->vals.output_interval == 0){
                    ST_Window_extract_slab_right_trace_to_nodal_value(
                        te,
                        nx,
                        num_window_slabs,
                        slab,
                        T_window,
                        sys->vals.T);

                    output_files(sys, file_num, time);
                    file_num++;
                }
            }
        }

        memcpy(
            T_prev_st,
            T_last_st,
            (size_t)previous_state_dof * sizeof(double));
    }

    ST_Window_extract_right_trace_to_nodal_value(
        te,
        nx,
        num_window_slabs,
        T_window,
        sys->vals.T);

    if(collect_context_initialized){
        ROM_std_strom_fom_context_finalize(&collect_context);
    }

    free(dirichlet_lift);
    free(homogeneous_snapshot);
    free(T_window);
    free(T_prev_st);
    free(T_last_st);
}

static void write_shared_basis_manifest(
    const ST_ROM_OPTIONS* options,
    const FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    int shared_snapshot_count,
    int st_dof,
    int num_modes,
    const double* singular_values,
    int num_singular_values)
{
    char filename[4096];
    FILE* fp;

    if(snprintf(
        filename,
        sizeof(filename),
        "%s/all_windows_shared_pod.info",
        options->basis_dir) >= (int)sizeof(filename))
    {
        fail_message("shared basis manifest path is too long");
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(
            stderr,
            "ERROR: cannot open shared basis manifest %s: %s\n",
            filename,
            strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "basis_scope=all_windows_shared\n");
    fprintf(fp, "snapshot_representation=homogeneous_lifted\n");
    fprintf(fp, "num_space_nodes=%d\n", sys->fe.total_num_nodes);
    fprintf(fp, "num_time_dofs=%d\n", te->n_dof);
    fprintf(fp, "num_slabs_per_window=%d\n", num_window_slabs);
    fprintf(fp, "dt=%.17g\n", sys->vals.dt);
    fprintf(fp, "finish_time=%.17g\n", sys->vals.finish_time);
    fprintf(fp, "num_windows=%d\n", num_windows);
    fprintf(fp, "num_training_cases=%d\n", options->num_snapshots);
    fprintf(fp, "num_snapshot_columns=%d\n", shared_snapshot_count);
    fprintf(fp, "st_dof_per_window=%d\n", st_dof);
    fprintf(fp, "num_modes=%d\n", num_modes);
    fprintf(fp, "pod_epsilon=%.17g\n", options->epsilon);
    fprintf(fp, "snapshot_column_order=training_case_major_then_window\n");

    for(int i = 0; i < num_singular_values; i++){
        fprintf(fp, "sigma_%04d=%.17e\n", i, singular_values[i]);
    }

    fclose(fp);
}

static void validate_shared_basis_manifest(
    const ST_ROM_OPTIONS* options,
    const FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    int loaded_num_modes)
{
    char filename[4096];
    char line[8192];
    char basis_scope[128] = "";
    char snapshot_representation[128] = "";
    int num_space_nodes = -1;
    int num_time_dofs = -1;
    int slabs_per_window = -1;
    int manifest_num_windows = -1;
    int st_dof_per_window = -1;
    int manifest_num_modes = -1;
    double dt = NAN;
    double finish_time = NAN;
    FILE* fp;

    if(snprintf(
        filename,
        sizeof(filename),
        "%s/all_windows_shared_pod.info",
        options->basis_dir) >= (int)sizeof(filename))
    {
        fail_message("shared basis manifest path is too long");
    }

    fp = fopen(filename, "r");
    if(fp == NULL){
        fprintf(
            stderr,
            "ERROR: cannot open shared basis manifest %s: %s\n",
            filename,
            strerror(errno));
        exit(EXIT_FAILURE);
    }

    while(fgets(line, sizeof(line), fp) != NULL){
        char* equals = strchr(line, '=');
        char* value;
        char* newline;

        if(equals == NULL){
            continue;
        }

        *equals = '\0';
        value = equals + 1;
        newline = strpbrk(value, "\r\n");
        if(newline != NULL){
            *newline = '\0';
        }

        if(strcmp(line, "basis_scope") == 0){
            snprintf(basis_scope, sizeof(basis_scope), "%s", value);
        }
        else if(strcmp(line, "snapshot_representation") == 0){
            snprintf(
                snapshot_representation,
                sizeof(snapshot_representation),
                "%s",
                value);
        }
        else if(strcmp(line, "num_space_nodes") == 0){
            num_space_nodes = atoi(value);
        }
        else if(strcmp(line, "num_time_dofs") == 0){
            num_time_dofs = atoi(value);
        }
        else if(strcmp(line, "num_slabs_per_window") == 0){
            slabs_per_window = atoi(value);
        }
        else if(strcmp(line, "dt") == 0){
            dt = strtod(value, NULL);
        }
        else if(strcmp(line, "finish_time") == 0){
            finish_time = strtod(value, NULL);
        }
        else if(strcmp(line, "num_windows") == 0){
            manifest_num_windows = atoi(value);
        }
        else if(strcmp(line, "st_dof_per_window") == 0){
            st_dof_per_window = atoi(value);
        }
        else if(strcmp(line, "num_modes") == 0){
            manifest_num_modes = atoi(value);
        }
    }

    fclose(fp);

    const int expected_st_dof = checked_product3_int(
        sys->fe.total_num_nodes,
        te->n_dof,
        num_window_slabs);

    const double dt_tolerance =
        1.0e-12 * fmax(1.0, fabs(sys->vals.dt));

    const double finish_tolerance =
        1.0e-12 * fmax(1.0, fabs(sys->vals.finish_time));

    if(strcmp(basis_scope, "all_windows_shared") != 0 ||
       strcmp(snapshot_representation, "homogeneous_lifted") != 0 ||
       num_space_nodes != sys->fe.total_num_nodes ||
       num_time_dofs != te->n_dof ||
       slabs_per_window != num_window_slabs ||
       manifest_num_windows != num_windows ||
       st_dof_per_window != expected_st_dof ||
       manifest_num_modes != loaded_num_modes ||
       !isfinite(dt) ||
       fabs(dt - sys->vals.dt) > dt_tolerance ||
       !isfinite(finish_time) ||
       fabs(finish_time - sys->vals.finish_time) > finish_tolerance)
    {
        fprintf(
            stderr,
            "ERROR: shared POD basis metadata does not match this run.\n"
            "       basis scope            : %s\n"
            "       snapshot representation : %s\n"
            "       nodes                   : file=%d run=%d\n"
            "       time dofs               : file=%d run=%d\n"
            "       slabs/window            : file=%d run=%d\n"
            "       windows                 : file=%d run=%d\n"
            "       st dof/window           : file=%d run=%d\n"
            "       modes                   : file=%d loaded=%d\n"
            "       dt                      : file=%.17g run=%.17g\n"
            "       finish time             : file=%.17g run=%.17g\n",
            basis_scope,
            snapshot_representation,
            num_space_nodes,
            sys->fe.total_num_nodes,
            num_time_dofs,
            te->n_dof,
            slabs_per_window,
            num_window_slabs,
            manifest_num_windows,
            num_windows,
            st_dof_per_window,
            expected_st_dof,
            manifest_num_modes,
            loaded_num_modes,
            dt,
            sys->vals.dt,
            finish_time,
            sys->vals.finish_time);
        exit(EXIT_FAILURE);
    }
}


static ROM_STROM_APPLY_WINDOW_COUPLING select_coupling_action(
    const ST_ROM_OPTIONS* options)
{
    return options->direct_coupling
        ? ROM_std_strom_fom_apply_B_direct
        : ROM_std_strom_fom_apply_B;
}

static void verify_direct_coupling_once(
    const ST_ROM_OPTIONS* options,
    ROM_STROM_SYSTEM* system,
    ROM_STROM_FOM_CONTEXT* context)
{
    double* x = NULL;
    double* direct_value = NULL;
    double* legacy_value = NULL;
    double difference_norm2 = 0.0;
    double scale_norm2 = 0.0;
    const ROM_STROM_WINDOW* basis_window;
    int status;

    if(!options->verify_direct_b ||
       system->num_windows < 2 ||
       !options->direct_coupling)
    {
        return;
    }

    basis_window = &system->windows[0];
    x = (double*)calloc((size_t)basis_window->st_dof, sizeof(double));
    direct_value = (double*)calloc((size_t)basis_window->st_dof, sizeof(double));
    legacy_value = (double*)calloc((size_t)basis_window->st_dof, sizeof(double));

    if(x == NULL || direct_value == NULL || legacy_value == NULL){
        free(x);
        free(direct_value);
        free(legacy_value);
        fail_message("direct B verification allocation failed");
    }

    for(int p = 0; p < basis_window->st_dof; p++){
        x[p] = basis_window->hlpod_mat->pod_modes[p][0];
    }

    status = ROM_std_strom_fom_apply_B_direct(
        1,
        x,
        direct_value,
        context);
    if(status != ROM_STROM_SUCCESS){
        fail_message("direct B verification action failed");
    }

    status = ROM_std_strom_fom_apply_B(
        1,
        x,
        legacy_value,
        context);
    if(status != ROM_STROM_SUCCESS){
        fail_message("legacy B verification action failed");
    }

    for(int i = 0; i < basis_window->st_dof; i++){
        const double difference = direct_value[i] - legacy_value[i];
        difference_norm2 += difference * difference;
        scale_norm2 += direct_value[i] * direct_value[i]
            + legacy_value[i] * legacy_value[i];
    }

    const double relative_difference =
        sqrt(difference_norm2) / fmax(sqrt(scale_norm2), 1.0e-30);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("direct/legacy B relative difference: %.15e\n",
            relative_difference);
        fflush(stdout);
    }

    free(x);
    free(direct_value);
    free(legacy_value);

    if(relative_difference > 1.0e-10){
        fail_message(
            "direct B action does not match the legacy RHS-difference action; "
            "set ST_ROM_DIRECT_COUPLING=0");
    }
}

static void build_reduced_operators_fast(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_FOM_CONTEXT* context,
    const ST_ROM_OPTIONS* options)
{
    const double start_time = monolis_get_time();
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B =
        select_coupling_action(options);
    int status;

    verify_direct_coupling_once(options, system, context);

    if(options->operator_reuse){
        system->windows[0].window_id = 0;
        status = ROM_std_strom_calc_reduced_diag_prepared(
            &system->windows[0],
            ROM_std_strom_fom_prepare_A_callback,
            ROM_std_strom_fom_apply_prepared_A_callback,
            context);
        if(status != ROM_STROM_SUCCESS){
            fail_message("shared D projection failed");
        }

        for(int w = 1; w < system->num_windows; w++){
            status = ROM_std_strom_copy_reduced_diag(
                &system->windows[w],
                &system->windows[0]);
            if(status != ROM_STROM_SUCCESS){
                fail_message("shared D copy failed");
            }
        }

        if(system->num_windows > 1){
            system->windows[1].window_id = 1;
            status = ROM_std_strom_calc_reduced_coupling_fast(
                &system->windows[1],
                &system->windows[0],
                apply_B,
                context);
            if(status != ROM_STROM_SUCCESS){
                fail_message("shared L projection failed");
            }

            for(int w = 2; w < system->num_windows; w++){
                status = ROM_std_strom_copy_reduced_coupling(
                    &system->windows[w],
                    &system->windows[1]);
                if(status != ROM_STROM_SUCCESS){
                    fail_message("shared L copy failed");
                }
            }
        }
    }
    else{
        for(int w = 0; w < system->num_windows; w++){
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("projecting reduced operators: window %d/%d\r",
                    w + 1,
                    system->num_windows);
                fflush(stdout);
            }

            status = ROM_std_strom_calc_reduced_diag_prepared(
                &system->windows[w],
                ROM_std_strom_fom_prepare_A_callback,
                ROM_std_strom_fom_apply_prepared_A_callback,
                context);
            if(status != ROM_STROM_SUCCESS){
                fail_message("D_w projection failed");
            }

            if(w > 0){
                status = ROM_std_strom_calc_reduced_coupling_fast(
                    &system->windows[w],
                    &system->windows[w - 1],
                    apply_B,
                    context);
                if(status != ROM_STROM_SUCCESS){
                    fail_message("L_w projection failed");
                }
            }
        }
        if(monolis_mpi_get_global_my_rank() == 0){
            printf("\n");
        }
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "reduced operator projection completed: reuse=%d direct_B=%d time=%.6f s\n",
            options->operator_reuse,
            options->direct_coupling,
            monolis_get_time() - start_time);
        fflush(stdout);
    }
}

static void run_offline(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    const ST_ROM_OPTIONS* options)
{
    ROM_STROM_WINDOW shared_basis_window = {0};
    ROM_STROM_SYSTEM operator_system = {0};
    ROM_STROM_FOM_CONTEXT fom_context = {0};
    HLPOD_VALUES shared_values = {0};
    HLPOD_MAT shared_matrix = {0};
    HLPOD_VALUES* operator_values = NULL;
    HLPOD_MAT* operator_matrices = NULL;

    const int st_dof = checked_product3_int(
        sys->fe.total_num_nodes,
        te->n_dof,
        num_window_slabs);

    const int shared_snapshot_count = checked_product3_int(
        options->num_snapshots,
        num_windows,
        1);

    const int requested_max_modes =
        options->max_modes_explicit
        ? options->max_modes
        : shared_snapshot_count;

    const int shared_max_modes =
        requested_max_modes < shared_snapshot_count
        ? requested_max_modes
        : shared_snapshot_count;

    double* snapshot = NULL;
    int status;

    if(shared_max_modes <= 0){
        fail_message("invalid shared POD mode capacity");
    }

    if(make_directory(options->basis_dir) != 0){
        fail_message("cannot create basis directory");
    }

    /*
     * All q_w snapshots use the same local window coordinate system and have
     * identical size.  Pool them as columns of one matrix:
     *
     *   S_shared = [q_0^(0), ..., q_{Nw-1}^(0),
     *               q_0^(1), ..., q_{Nw-1}^(1), ...].
     *
     * This is not vertical concatenation of the full trajectory.  Each
     * window remains one column, so one training trajectory already supplies
     * num_windows POD samples.
     */
    status = ROM_std_strom_window_initialize(
        &shared_basis_window,
        0,
        &shared_values,
        &shared_matrix,
        sys->fe.total_num_nodes,
        1,
        num_window_slabs,
        te->n_dof,
        shared_snapshot_count,
        shared_max_modes,
        0);

    if(status != ROM_STROM_SUCCESS){
        fail_message("shared POD holder initialization failed");
    }

    snapshot = (double*)calloc((size_t)st_dof, sizeof(double));
    if(snapshot == NULL){
        fail_message("shared snapshot work allocation failed");
    }

    for(int m = 0; m < options->num_snapshots; m++){
        for(int w = 0; w < num_windows; w++){
            const int shared_column = m * num_windows + w;

            status = ROM_std_strom_fom_read_snapshot(
                options->snapshot_dir,
                m,
                w,
                snapshot,
                st_dof);

            if(status != ROM_STROM_SUCCESS){
                fprintf(
                    stderr,
                    "ERROR: cannot read snapshot m=%d w=%d from %s\n",
                    m,
                    w,
                    options->snapshot_dir);
                exit(EXIT_FAILURE);
            }

            status = ROM_std_strom_store_snapshot(
                &shared_basis_window,
                snapshot,
                shared_column,
                0);

            if(status != ROM_STROM_SUCCESS){
                fail_message("shared snapshot store failed");
            }
        }
    }

    status = ROM_std_strom_set_pod_modes(
        &shared_basis_window,
        options->epsilon);

    if(status != ROM_STROM_SUCCESS){
        fail_message("all-window shared ST-POD failed");
    }

    status = ROM_std_strom_fom_write_shared_basis(
        options->basis_dir,
        &shared_basis_window);

    if(status != ROM_STROM_SUCCESS){
        fail_message("shared POD basis write failed");
    }

    write_shared_basis_manifest(
        options,
        sys,
        te,
        num_window_slabs,
        num_windows,
        shared_snapshot_count,
        st_dof,
        shared_basis_window.num_modes,
        shared_basis_window.singular_values,
        shared_basis_window.num_singular_values);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("offline POD basis scope: all_windows_shared\n");
        printf("shared snapshot columns: %d (%d training cases x %d windows)\n",
            shared_snapshot_count,
            options->num_snapshots,
            num_windows);
        printf("shared POD modes: %d\n", shared_basis_window.num_modes);
        printf("shared sigma_1: %e\n",
            shared_basis_window.num_singular_values > 0
                ? shared_basis_window.singular_values[0]
                : 0.0);
        printf("shared basis file: %s/all_windows_shared_pod.bin\n",
            options->basis_dir);
    }

    /*
     * Complete offline/online split:
     *   - attach one shared basis without deep copies,
     *   - prepare A once per required block,
     *   - form A*Phi and B*Phi,
     *   - project with BLAS,
     *   - persist D/L for online loading.
     */
    initialize_strom_system(
        &operator_system,
        &operator_values,
        &operator_matrices,
        sys,
        te,
        num_window_slabs,
        num_windows,
        1,
        shared_basis_window.num_modes);

    for(int w = 0; w < num_windows; w++){
        status = ROM_std_strom_attach_pod_basis(
            &operator_system.windows[w],
            &shared_basis_window);
        if(status != ROM_STROM_SUCCESS){
            fail_message("cannot attach shared basis for offline projection");
        }
    }

    status = ROM_std_strom_update_mode_offsets(&operator_system);
    if(status != ROM_STROM_SUCCESS){
        fail_message("offline mode offset construction failed");
    }

    status = ROM_std_strom_fom_context_initialize(
        &fom_context,
        sys,
        te,
        num_window_slabs,
        num_windows);
    if(status != ROM_STROM_SUCCESS){
        fail_message("offline FOM adapter initialization failed");
    }

    build_reduced_operators_fast(
        &operator_system,
        &fom_context,
        options);

    status = ROM_std_strom_fom_write_reduced_operators(
        options->basis_dir,
        &operator_system,
        options->operator_reuse,
        sys->vals.dt,
        st_dof);
    if(status != ROM_STROM_SUCCESS){
        fail_message("reduced operator write failed");
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "reduced operator file: %s/all_windows_shared_reduced_operators.bin\n",
            options->basis_dir);
    }

    ROM_std_strom_fom_context_finalize(&fom_context);
    finalize_strom_system(
        &operator_system,
        &operator_values,
        &operator_matrices);
    free(snapshot);
    ROM_std_strom_window_finalize(&shared_basis_window);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("All-window shared POD/operator offline completed.\n");
        fflush(stdout);
    }

}


/* -------------------------------------------------------------------------
 * Online accuracy validation and reduced-equation diagnostics
 * ------------------------------------------------------------------------- */

typedef struct {
    double projection_num2;
    double projection_den2;

    double state_num2;
    double state_den2;

    ROM_STROM_RESIDUAL_NORMS rom_residual_all;
    ROM_STROM_RESIDUAL_NORMS rom_residual_free;
    ROM_STROM_RESIDUAL_NORMS rom_residual_dirichlet;

    ROM_STROM_RESIDUAL_NORMS reference_residual_all;
    ROM_STROM_RESIDUAL_NORMS reference_residual_free;
    ROM_STROM_RESIDUAL_NORMS reference_residual_dirichlet;

    double rom_full_hom_difference_num2;
    double rom_full_hom_scale_num2;
    double reference_full_hom_difference_num2;
    double reference_full_hom_scale_num2;

    double reduced_defect_num2;
    double reduced_defect_lhs2;
    double reduced_defect_rhs2;

    double coefficient_num2;
    double coefficient_den2;

    double trace_num2;
    double trace_den2;

    double max_projection_error;
    double max_state_error;

    double max_rom_residual_all;
    double max_rom_residual_free;
    double max_rom_residual_dirichlet;

    double max_reference_residual_all;
    double max_reference_residual_free;
    double max_reference_residual_dirichlet;

    double max_rom_full_hom_identity_error;
    double max_reference_full_hom_identity_error;

    double max_reduced_defect;
    double max_coefficient_error;
    double max_trace_error;
} ST_ROM_ACCURACY;

typedef struct {
    int enabled;
    int snapshot_id;
    int st_dof;
    int num_space_nodes;
    int previous_state_dof;
    int max_modes;
    int previous_reference_modes;

    FILE* csv;

    /* Version-2 snapshot q_ref = U_ref - U_D. */
    double* reference_homogeneous;
    double* projected_homogeneous;
    double* reference_full;

    /* Shared lifted-residual workspace. */
    double* full_rhs;
    double* matrix_times_lift;
    double* homogeneous_rhs;
    double* matrix_times_homogeneous;
    double* homogeneous_residual;
    double* full_residual;

    double* fom_trace;
    double* rom_trace;

    /* Full last-slab states, including the Dirichlet lift. */
    double* previous_rom_state;
    double* previous_reference_state;

    double* reference_coef;
    double* previous_reference_coef;
    double* reduced_lhs;
    double* reduced_rhs_effective;

    ST_ROM_ACCURACY accuracy;
} ST_ROM_VALIDATION_WORK;

static double safe_relative_from_squares(
    double numerator2,
    double denominator2)
{
    const double denominator = sqrt(fmax(denominator2, 0.0));

    return sqrt(fmax(numerator2, 0.0))
        / fmax(denominator, 1.0e-30);
}

static void calculate_error_squares(
    const double* reference,
    const double* approximation,
    int n,
    double* numerator2,
    double* denominator2)
{
    double num2 = 0.0;
    double den2 = 0.0;

    for(int i = 0; i < n; i++){
        const double difference = reference[i] - approximation[i];
        num2 += difference * difference;
        den2 += reference[i] * reference[i];
    }

    *numerator2 = num2;
    *denominator2 = den2;
}

static double residual_relative_error(
    double residual2,
    double lhs2,
    double rhs2)
{
    const double scale =
        sqrt(fmax(lhs2, 0.0)) + sqrt(fmax(rhs2, 0.0));

    return sqrt(fmax(residual2, 0.0))
        / fmax(scale, 1.0e-30);
}

static void update_maximum(double* maximum, double value)
{
    if(value > *maximum){
        *maximum = value;
    }
}

static void calculate_window_coefficients(
    const ROM_STROM_WINDOW* window,
    const double* state,
    double* coefficients)
{
    if(window == NULL ||
       state == NULL ||
       coefficients == NULL ||
       window->num_modes <= 0 ||
       window->hlpod_mat == NULL ||
       window->hlpod_mat->pod_modes == NULL)
    {
        fail_message("invalid coefficient projection arguments");
    }

    for(int j = 0; j < window->num_modes; j++){
        double value = 0.0;

        for(int p = 0; p < window->st_dof; p++){
            const double reference_value =
                window->reference_state != NULL
                ? window->reference_state[p]
                : 0.0;

            value +=
                window->hlpod_mat->pod_modes[p][j]
                * (state[p] - reference_value);
        }

        coefficients[j] = value;
    }
}

static void calculate_reduced_reference_defect(
    const ROM_STROM_WINDOW* window,
    const double* reference_coef,
    const double* previous_reference_coef,
    int previous_reference_modes,
    double* reduced_lhs,
    double* reduced_rhs_effective,
    double* defect_num2,
    double* lhs_num2,
    double* rhs_num2)
{
    double defect2 = 0.0;
    double left2 = 0.0;
    double right2 = 0.0;

    if(window == NULL ||
       reference_coef == NULL ||
       reduced_lhs == NULL ||
       reduced_rhs_effective == NULL ||
       defect_num2 == NULL ||
       lhs_num2 == NULL ||
       rhs_num2 == NULL ||
       window->reduced_diag == NULL ||
       window->reduced_rhs == NULL)
    {
        fail_message("invalid reduced reference defect arguments");
    }

    if(window->window_id > 0 &&
       (previous_reference_coef == NULL ||
        previous_reference_modes <= 0 ||
        window->reduced_coupling_prev == NULL))
    {
        fail_message("missing previous-window reduced diagnostic data");
    }

    for(int i = 0; i < window->num_modes; i++){
        double lhs = 0.0;
        double rhs = window->reduced_rhs[i];

        for(int j = 0; j < window->num_modes; j++){
            lhs += window->reduced_diag[i][j] * reference_coef[j];
        }

        if(window->window_id > 0){
            for(int j = 0; j < previous_reference_modes; j++){
                rhs +=
                    window->reduced_coupling_prev[i][j]
                    * previous_reference_coef[j];
            }
        }

        reduced_lhs[i] = lhs;
        reduced_rhs_effective[i] = rhs;

        const double difference = lhs - rhs;
        defect2 += difference * difference;
        left2 += lhs * lhs;
        right2 += rhs * rhs;
    }

    *defect_num2 = defect2;
    *lhs_num2 = left2;
    *rhs_num2 = right2;
}

static void initialize_validation_work(
    ST_ROM_VALIDATION_WORK* work,
    const ST_ROM_OPTIONS* options,
    const FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int max_modes)
{
    const int nx = sys->fe.total_num_nodes;
    const int previous_state_dof = checked_product3_int(nx, te->n_dof, 1);
    const int st_dof = checked_product3_int(
        nx,
        te->n_dof,
        num_window_slabs);

    memset(work, 0, sizeof(*work));

    work->enabled = options->validate_online;
    work->snapshot_id = options->validation_snapshot_id;
    work->st_dof = st_dof;
    work->num_space_nodes = nx;
    work->previous_state_dof = previous_state_dof;
    work->max_modes = max_modes;

    if(!work->enabled){
        return;
    }

    if(max_modes <= 0){
        fail_message("invalid validation max_modes");
    }

    work->reference_homogeneous = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->projected_homogeneous = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->reference_full = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->full_rhs = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->matrix_times_lift = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->homogeneous_rhs = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->matrix_times_homogeneous = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->homogeneous_residual = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->full_residual = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    work->fom_trace = (double*)calloc(
        (size_t)nx,
        sizeof(double));

    work->rom_trace = (double*)calloc(
        (size_t)nx,
        sizeof(double));

    work->previous_rom_state = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    work->previous_reference_state = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    work->reference_coef = (double*)calloc(
        (size_t)max_modes,
        sizeof(double));

    work->previous_reference_coef = (double*)calloc(
        (size_t)max_modes,
        sizeof(double));

    work->reduced_lhs = (double*)calloc(
        (size_t)max_modes,
        sizeof(double));

    work->reduced_rhs_effective = (double*)calloc(
        (size_t)max_modes,
        sizeof(double));

    if(work->reference_homogeneous == NULL ||
       work->projected_homogeneous == NULL ||
       work->reference_full == NULL ||
       work->full_rhs == NULL ||
       work->matrix_times_lift == NULL ||
       work->homogeneous_rhs == NULL ||
       work->matrix_times_homogeneous == NULL ||
       work->homogeneous_residual == NULL ||
       work->full_residual == NULL ||
       work->fom_trace == NULL ||
       work->rom_trace == NULL ||
       work->previous_rom_state == NULL ||
       work->previous_reference_state == NULL ||
       work->reference_coef == NULL ||
       work->previous_reference_coef == NULL ||
       work->reduced_lhs == NULL ||
       work->reduced_rhs_effective == NULL)
    {
        fail_message("validation work allocation failed");
    }

    work->csv = fopen(options->validation_csv, "w");
    if(work->csv == NULL){
        fprintf(
            stderr,
            "ERROR: cannot open validation CSV %s: %s\n",
            options->validation_csv,
            strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(
        work->csv,
        "scope,window,time_start,time_end,num_modes,"
        "projection_error,state_error,"
        "rom_residual_all,rom_residual_free,rom_residual_dirichlet,"
        "reference_residual_all,reference_residual_free,"
        "reference_residual_dirichlet,"
        "rom_full_hom_identity_error,"
        "reference_full_hom_identity_error,"
        "reduced_reference_defect,coefficient_error,trace_error\n");

    /* Both are full physical initial states for window zero. */
    ST_TimeElement_set_prev_trace_from_nodal_value(
        te,
        nx,
        sys->vals.T,
        work->previous_rom_state);

    ST_TimeElement_set_prev_trace_from_nodal_value(
        te,
        nx,
        sys->vals.T,
        work->previous_reference_state);
}

static void finalize_validation_work(
    ST_ROM_VALIDATION_WORK* work)
{
    if(work == NULL){
        return;
    }

    if(work->csv != NULL){
        fclose(work->csv);
    }

    free(work->reference_homogeneous);
    free(work->projected_homogeneous);
    free(work->reference_full);

    free(work->full_rhs);
    free(work->matrix_times_lift);
    free(work->homogeneous_rhs);
    free(work->matrix_times_homogeneous);
    free(work->homogeneous_residual);
    free(work->full_residual);

    free(work->fom_trace);
    free(work->rom_trace);

    free(work->previous_rom_state);
    free(work->previous_reference_state);

    free(work->reference_coef);
    free(work->previous_reference_coef);
    free(work->reduced_lhs);
    free(work->reduced_rhs_effective);

    memset(work, 0, sizeof(*work));
}

static void print_scalar_reduced_diagnostics(
    const ROM_STROM_WINDOW* window,
    const ST_ROM_VALIDATION_WORK* work,
    double reference_residual_all,
    double reference_residual_free,
    double reference_residual_dirichlet,
    double reduced_reference_defect,
    double coefficient_error)
{
    if(window->window_id != 0 || window->num_modes != 1){
        return;
    }

    const double D00 = window->reduced_diag[0][0];
    const double g0 = window->reduced_rhs[0];
    const double a_reference = work->reference_coef[0];
    const double a_rom = window->reduced_coef[0];
    const double direct_scalar =
        fabs(D00) > 1.0e-30
        ? g0 / D00
        : NAN;

    printf("\nwindow-0 scalar reduced diagnostics\n");
    printf("  reference residual all    : %.15e\n", reference_residual_all);
    printf("  reference residual free   : %.15e\n", reference_residual_free);
    printf("  reference residual BC     : %.15e\n", reference_residual_dirichlet);
    printf("  D00                       : %.15e\n", D00);
    printf("  g0                        : %.15e\n", g0);
    printf("  a_reference               : %.15e\n", a_reference);
    printf("  a_rom                     : %.15e\n", a_rom);
    printf("  D00 * a_reference         : %.15e\n", D00 * a_reference);
    printf("  direct scalar g0 / D00    : %.15e\n", direct_scalar);
    printf("  reduced reference defect  : %.15e\n", reduced_reference_defect);
    printf("  coefficient error         : %.15e\n\n", coefficient_error);
}

static void validate_online_window(
    ST_ROM_VALIDATION_WORK* work,
    const ST_ROM_OPTIONS* options,
    ROM_STROM_FOM_CONTEXT* fom_context,
    const ROM_STROM_WINDOW* window,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double dt,
    const double* current_lift,
    const double* reconstructed_homogeneous,
    const double* reconstructed_full)
{
    const int w = window->window_id;
    const double t_window_start =
        (double)(w * num_window_slabs) * dt;
    const double t_window_end =
        t_window_start + (double)num_window_slabs * dt;

    double projection_num2;
    double projection_den2;
    double state_num2;
    double state_den2;

    double reduced_defect_num2;
    double reduced_defect_lhs2;
    double reduced_defect_rhs2;

    double coefficient_num2;
    double coefficient_den2;

    double trace_num2;
    double trace_den2;

    double projection_error;
    double state_error;

    double rom_residual_all;
    double rom_residual_free;
    double rom_residual_dirichlet;

    double reference_residual_all;
    double reference_residual_free;
    double reference_residual_dirichlet;

    double rom_full_hom_identity_error;
    double reference_full_hom_identity_error;

    double reduced_reference_defect;
    double coefficient_error;
    double trace_error;

    ROM_STROM_LIFTED_RESIDUAL rom_lifted_residual;
    ROM_STROM_LIFTED_RESIDUAL reference_lifted_residual;
    int status;

    if(!work->enabled){
        return;
    }

    if(current_lift == NULL ||
       reconstructed_homogeneous == NULL ||
       reconstructed_full == NULL)
    {
        fail_message("missing lifted validation state");
    }

    if(window->num_modes > work->max_modes){
        fail_message("validation mode workspace is too small");
    }

    memset(
        work->reference_coef,
        0,
        (size_t)work->max_modes * sizeof(double));

    memset(
        work->reduced_lhs,
        0,
        (size_t)work->max_modes * sizeof(double));

    memset(
        work->reduced_rhs_effective,
        0,
        (size_t)work->max_modes * sizeof(double));

    status = ROM_std_strom_fom_read_snapshot(
        options->snapshot_dir,
        work->snapshot_id,
        w,
        work->reference_homogeneous,
        work->st_dof);

    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: cannot read homogeneous validation snapshot "
            "id=%d window=%d from %s\n",
            work->snapshot_id,
            w,
            options->snapshot_dir);
        exit(EXIT_FAILURE);
    }

    status = ROM_std_strom_project_window_state(
        window,
        work->reference_homogeneous,
        work->projected_homogeneous);

    if(status != ROM_STROM_SUCCESS){
        fail_message("homogeneous POD projection validation failed");
    }

    for(int i = 0; i < work->st_dof; i++){
        work->reference_full[i] =
            current_lift[i] + work->reference_homogeneous[i];
    }

    /* POD representation is assessed in q. */
    calculate_error_squares(
        work->reference_homogeneous,
        work->projected_homogeneous,
        work->st_dof,
        &projection_num2,
        &projection_den2);

    /* Physical state error is assessed in U = U_D + q. */
    calculate_error_squares(
        work->reference_full,
        reconstructed_full,
        work->st_dof,
        &state_num2,
        &state_den2);

    /*
     * Both residual formulations use the same matched pair
     * (A_tilde, b_tilde) inside each call.
     */
    status = ROM_std_strom_fom_calculate_lifted_residual(
        fom_context,
        w,
        work->previous_rom_state,
        reconstructed_homogeneous,
        current_lift,
        work->full_rhs,
        work->matrix_times_lift,
        work->homogeneous_rhs,
        work->matrix_times_homogeneous,
        work->homogeneous_residual,
        work->full_residual,
        &rom_lifted_residual);

    if(status != ROM_STROM_SUCCESS){
        fail_message("consistent lifted ROM residual validation failed");
    }

    status = ROM_std_strom_fom_calculate_lifted_residual(
        fom_context,
        w,
        work->previous_reference_state,
        work->reference_homogeneous,
        current_lift,
        work->full_rhs,
        work->matrix_times_lift,
        work->homogeneous_rhs,
        work->matrix_times_homogeneous,
        work->homogeneous_residual,
        work->full_residual,
        &reference_lifted_residual);

    if(status != ROM_STROM_SUCCESS){
        fail_message("consistent lifted reference residual validation failed");
    }

    calculate_window_coefficients(
        window,
        work->reference_homogeneous,
        work->reference_coef);

    calculate_reduced_reference_defect(
        window,
        work->reference_coef,
        work->previous_reference_coef,
        work->previous_reference_modes,
        work->reduced_lhs,
        work->reduced_rhs_effective,
        &reduced_defect_num2,
        &reduced_defect_lhs2,
        &reduced_defect_rhs2);

    calculate_error_squares(
        work->reference_coef,
        window->reduced_coef,
        window->num_modes,
        &coefficient_num2,
        &coefficient_den2);

    ST_Window_extract_right_trace_to_nodal_value(
        te,
        work->num_space_nodes,
        num_window_slabs,
        work->reference_full,
        work->fom_trace);

    ST_Window_extract_right_trace_to_nodal_value(
        te,
        work->num_space_nodes,
        num_window_slabs,
        reconstructed_full,
        work->rom_trace);

    calculate_error_squares(
        work->fom_trace,
        work->rom_trace,
        work->num_space_nodes,
        &trace_num2,
        &trace_den2);

    projection_error = safe_relative_from_squares(
        projection_num2,
        projection_den2);

    state_error = safe_relative_from_squares(
        state_num2,
        state_den2);

    rom_residual_all = ROM_std_strom_relative_residual(
        &rom_lifted_residual.all);
    rom_residual_free = ROM_std_strom_relative_residual(
        &rom_lifted_residual.free_dofs);
    rom_residual_dirichlet = ROM_std_strom_relative_residual(
        &rom_lifted_residual.dirichlet_dofs);

    reference_residual_all = ROM_std_strom_relative_residual(
        &reference_lifted_residual.all);
    reference_residual_free = ROM_std_strom_relative_residual(
        &reference_lifted_residual.free_dofs);
    reference_residual_dirichlet = ROM_std_strom_relative_residual(
        &reference_lifted_residual.dirichlet_dofs);

    rom_full_hom_identity_error =
        sqrt(fmax(
            rom_lifted_residual.full_hom_difference_norm2,
            0.0))
        / fmax(
            sqrt(fmax(
                rom_lifted_residual.full_hom_scale_norm2,
                0.0)),
            1.0e-30);

    reference_full_hom_identity_error =
        sqrt(fmax(
            reference_lifted_residual.full_hom_difference_norm2,
            0.0))
        / fmax(
            sqrt(fmax(
                reference_lifted_residual.full_hom_scale_norm2,
                0.0)),
            1.0e-30);

    reduced_reference_defect = residual_relative_error(
        reduced_defect_num2,
        reduced_defect_lhs2,
        reduced_defect_rhs2);

    coefficient_error = safe_relative_from_squares(
        coefficient_num2,
        coefficient_den2);

    trace_error = safe_relative_from_squares(
        trace_num2,
        trace_den2);

    work->accuracy.projection_num2 += projection_num2;
    work->accuracy.projection_den2 += projection_den2;

    work->accuracy.state_num2 += state_num2;
    work->accuracy.state_den2 += state_den2;

    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.rom_residual_all,
        &rom_lifted_residual.all);
    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.rom_residual_free,
        &rom_lifted_residual.free_dofs);
    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.rom_residual_dirichlet,
        &rom_lifted_residual.dirichlet_dofs);

    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.reference_residual_all,
        &reference_lifted_residual.all);
    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.reference_residual_free,
        &reference_lifted_residual.free_dofs);
    ROM_std_strom_accumulate_residual_norms(
        &work->accuracy.reference_residual_dirichlet,
        &reference_lifted_residual.dirichlet_dofs);

    work->accuracy.rom_full_hom_difference_num2 +=
        rom_lifted_residual.full_hom_difference_norm2;
    work->accuracy.rom_full_hom_scale_num2 +=
        rom_lifted_residual.full_hom_scale_norm2;

    work->accuracy.reference_full_hom_difference_num2 +=
        reference_lifted_residual.full_hom_difference_norm2;
    work->accuracy.reference_full_hom_scale_num2 +=
        reference_lifted_residual.full_hom_scale_norm2;

    work->accuracy.reduced_defect_num2 += reduced_defect_num2;
    work->accuracy.reduced_defect_lhs2 += reduced_defect_lhs2;
    work->accuracy.reduced_defect_rhs2 += reduced_defect_rhs2;

    work->accuracy.coefficient_num2 += coefficient_num2;
    work->accuracy.coefficient_den2 += coefficient_den2;

    work->accuracy.trace_num2 += trace_num2;
    work->accuracy.trace_den2 += trace_den2;

    update_maximum(
        &work->accuracy.max_projection_error,
        projection_error);
    update_maximum(
        &work->accuracy.max_state_error,
        state_error);

    update_maximum(
        &work->accuracy.max_rom_residual_all,
        rom_residual_all);
    update_maximum(
        &work->accuracy.max_rom_residual_free,
        rom_residual_free);
    update_maximum(
        &work->accuracy.max_rom_residual_dirichlet,
        rom_residual_dirichlet);

    update_maximum(
        &work->accuracy.max_reference_residual_all,
        reference_residual_all);
    update_maximum(
        &work->accuracy.max_reference_residual_free,
        reference_residual_free);
    update_maximum(
        &work->accuracy.max_reference_residual_dirichlet,
        reference_residual_dirichlet);

    update_maximum(
        &work->accuracy.max_rom_full_hom_identity_error,
        rom_full_hom_identity_error);
    update_maximum(
        &work->accuracy.max_reference_full_hom_identity_error,
        reference_full_hom_identity_error);

    update_maximum(
        &work->accuracy.max_reduced_defect,
        reduced_reference_defect);
    update_maximum(
        &work->accuracy.max_coefficient_error,
        coefficient_error);
    update_maximum(
        &work->accuracy.max_trace_error,
        trace_error);

    printf(
        "validation window %d: modes=%d projection=%e state=%e "
        "rom_res_all=%e rom_res_free=%e rom_res_bc=%e "
        "ref_res_all=%e ref_res_free=%e ref_res_bc=%e "
        "rom_full_hom=%e ref_full_hom=%e "
        "red_def=%e coef=%e trace=%e\n",
        w,
        window->num_modes,
        projection_error,
        state_error,
        rom_residual_all,
        rom_residual_free,
        rom_residual_dirichlet,
        reference_residual_all,
        reference_residual_free,
        reference_residual_dirichlet,
        rom_full_hom_identity_error,
        reference_full_hom_identity_error,
        reduced_reference_defect,
        coefficient_error,
        trace_error);

    fprintf(
        work->csv,
        "window,%d,%.15e,%.15e,%d,"
        "%.15e,%.15e,"
        "%.15e,%.15e,%.15e,"
        "%.15e,%.15e,%.15e,"
        "%.15e,%.15e,"
        "%.15e,%.15e,%.15e\n",
        w,
        t_window_start,
        t_window_end,
        window->num_modes,
        projection_error,
        state_error,
        rom_residual_all,
        rom_residual_free,
        rom_residual_dirichlet,
        reference_residual_all,
        reference_residual_free,
        reference_residual_dirichlet,
        rom_full_hom_identity_error,
        reference_full_hom_identity_error,
        reduced_reference_defect,
        coefficient_error,
        trace_error);

    print_scalar_reduced_diagnostics(
        window,
        work,
        reference_residual_all,
        reference_residual_free,
        reference_residual_dirichlet,
        reduced_reference_defect,
        coefficient_error);

    /* Pass the full physical state U = U_D + q to the next window. */
    ST_Window_extract_last_slab_state(
        te,
        work->num_space_nodes,
        num_window_slabs,
        reconstructed_full,
        work->previous_rom_state);

    ST_Window_extract_last_slab_state(
        te,
        work->num_space_nodes,
        num_window_slabs,
        work->reference_full,
        work->previous_reference_state);

    memset(
        work->previous_reference_coef,
        0,
        (size_t)work->max_modes * sizeof(double));

    memcpy(
        work->previous_reference_coef,
        work->reference_coef,
        (size_t)window->num_modes * sizeof(double));

    work->previous_reference_modes = window->num_modes;
}

static void write_validation_summary(
    ST_ROM_VALIDATION_WORK* work,
    const FE_SYSTEM* sys,
    const ROM_STROM_SYSTEM* system)
{
    double projection_error;
    double state_error;

    double rom_residual_all;
    double rom_residual_free;
    double rom_residual_dirichlet;

    double reference_residual_all;
    double reference_residual_free;
    double reference_residual_dirichlet;

    double rom_full_hom_identity_error;
    double reference_full_hom_identity_error;

    double reduced_reference_defect;
    double coefficient_error;
    double trace_error;

    if(!work->enabled){
        return;
    }

    projection_error = safe_relative_from_squares(
        work->accuracy.projection_num2,
        work->accuracy.projection_den2);

    state_error = safe_relative_from_squares(
        work->accuracy.state_num2,
        work->accuracy.state_den2);

    rom_residual_all = ROM_std_strom_relative_residual(
        &work->accuracy.rom_residual_all);
    rom_residual_free = ROM_std_strom_relative_residual(
        &work->accuracy.rom_residual_free);
    rom_residual_dirichlet = ROM_std_strom_relative_residual(
        &work->accuracy.rom_residual_dirichlet);

    reference_residual_all = ROM_std_strom_relative_residual(
        &work->accuracy.reference_residual_all);
    reference_residual_free = ROM_std_strom_relative_residual(
        &work->accuracy.reference_residual_free);
    reference_residual_dirichlet = ROM_std_strom_relative_residual(
        &work->accuracy.reference_residual_dirichlet);

    rom_full_hom_identity_error =
        sqrt(fmax(
            work->accuracy.rom_full_hom_difference_num2,
            0.0))
        / fmax(
            sqrt(fmax(
                work->accuracy.rom_full_hom_scale_num2,
                0.0)),
            1.0e-30);

    reference_full_hom_identity_error =
        sqrt(fmax(
            work->accuracy.reference_full_hom_difference_num2,
            0.0))
        / fmax(
            sqrt(fmax(
                work->accuracy.reference_full_hom_scale_num2,
                0.0)),
            1.0e-30);

    reduced_reference_defect = residual_relative_error(
        work->accuracy.reduced_defect_num2,
        work->accuracy.reduced_defect_lhs2,
        work->accuracy.reduced_defect_rhs2);

    coefficient_error = safe_relative_from_squares(
        work->accuracy.coefficient_num2,
        work->accuracy.coefficient_den2);

    trace_error = safe_relative_from_squares(
        work->accuracy.trace_num2,
        work->accuracy.trace_den2);

    fprintf(
        work->csv,
        "global,-1,%.15e,%.15e,%d,"
        "%.15e,%.15e,"
        "%.15e,%.15e,%.15e,"
        "%.15e,%.15e,%.15e,"
        "%.15e,%.15e,"
        "%.15e,%.15e,%.15e\n",
        0.0,
        sys->vals.finish_time,
        system->total_reduced_dof,
        projection_error,
        state_error,
        rom_residual_all,
        rom_residual_free,
        rom_residual_dirichlet,
        reference_residual_all,
        reference_residual_free,
        reference_residual_dirichlet,
        rom_full_hom_identity_error,
        reference_full_hom_identity_error,
        reduced_reference_defect,
        coefficient_error,
        trace_error);

    fflush(work->csv);

    printf("\nST-ROM validation summary\n");
    printf("  validation snapshot id        : %d\n", work->snapshot_id);
    printf("  global projection error       : %e\n", projection_error);
    printf("  global state error            : %e\n", state_error);

    printf("  global ROM residual all       : %e\n", rom_residual_all);
    printf("  global ROM residual free      : %e\n", rom_residual_free);
    printf("  global ROM residual BC        : %e\n", rom_residual_dirichlet);

    printf("  global reference residual all : %e\n", reference_residual_all);
    printf("  global reference residual free: %e\n", reference_residual_free);
    printf("  global reference residual BC  : %e\n", reference_residual_dirichlet);

    printf("  global ROM full/hom identity  : %e\n", rom_full_hom_identity_error);
    printf("  global ref full/hom identity  : %e\n", reference_full_hom_identity_error);

    printf("  global reduced defect         : %e\n", reduced_reference_defect);
    printf("  global coefficient error      : %e\n", coefficient_error);
    printf("  global trace error            : %e\n", trace_error);

    printf("  max projection error          : %e\n", work->accuracy.max_projection_error);
    printf("  max state error               : %e\n", work->accuracy.max_state_error);

    printf("  max ROM residual all          : %e\n", work->accuracy.max_rom_residual_all);
    printf("  max ROM residual free         : %e\n", work->accuracy.max_rom_residual_free);
    printf("  max ROM residual BC           : %e\n", work->accuracy.max_rom_residual_dirichlet);

    printf("  max reference residual all    : %e\n", work->accuracy.max_reference_residual_all);
    printf("  max reference residual free   : %e\n", work->accuracy.max_reference_residual_free);
    printf("  max reference residual BC     : %e\n", work->accuracy.max_reference_residual_dirichlet);

    printf("  max ROM full/hom identity     : %e\n", work->accuracy.max_rom_full_hom_identity_error);
    printf("  max ref full/hom identity     : %e\n", work->accuracy.max_reference_full_hom_identity_error);

    printf("  max reduced defect            : %e\n", work->accuracy.max_reduced_defect);
    printf("  max coefficient error         : %e\n", work->accuracy.max_coefficient_error);
    printf("  max trace error               : %e\n", work->accuracy.max_trace_error);
}

static void pack_reduced_coefficients(
    const ROM_STROM_SYSTEM* system,
    double* packed)
{
    for(int w = 0; w < system->num_windows; w++){
        const ROM_STROM_WINDOW* window = &system->windows[w];
        const int offset = system->mode_offset[w];

        for(int i = 0; i < window->num_modes; i++){
            packed[offset + i] = window->reduced_coef[i];
        }
    }
}

static double relative_reduced_solution_difference(
    const double* reference,
    const double* comparison,
    int size)
{
    double numerator2 = 0.0;
    double denominator2 = 0.0;

    for(int i = 0; i < size; i++){
        const double difference = comparison[i] - reference[i];
        numerator2 += difference * difference;
        denominator2 += reference[i] * reference[i];
    }

    return sqrt(numerator2) /
        fmax(sqrt(denominator2), 1.0e-30);
}

static int solve_reduced_system(
    ROM_STROM_SYSTEM* system,
    ST_ROM_ONLINE_SOLVER solver,
    MONOLIS_COM* monolis_com,
    const ST_ROM_OPTIONS* options)
{
    int status;

    if(system == NULL || monolis_com == NULL || options == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(solver == ST_ROM_SOLVER_SEQUENTIAL){
        return ROM_std_strom_solve_sequential(system, NULL, NULL);
    }

    if(solver == ST_ROM_SOLVER_ALL_AT_ONCE){
        return ROM_std_strom_solve_all_at_once_dense(system, NULL, NULL);
    }

    if(solver == ST_ROM_SOLVER_MONOLIS ||
       solver == ST_ROM_SOLVER_VERIFY_MONOLIS)
    {
        ROM_STROM_MONOLIS_SOLVER monolis_solver;
        double* dense_solution = NULL;
        double* monolis_solution = NULL;

        memset(&monolis_solver, 0, sizeof(monolis_solver));

        if(solver == ST_ROM_SOLVER_VERIFY_MONOLIS){
            dense_solution = (double*)calloc(
                (size_t)system->total_reduced_dof,
                sizeof(double));

            monolis_solution = (double*)calloc(
                (size_t)system->total_reduced_dof,
                sizeof(double));

            if(dense_solution == NULL || monolis_solution == NULL){
                free(dense_solution);
                free(monolis_solution);
                return ROM_STROM_ERR_ALLOCATION;
            }

            status = ROM_std_strom_solve_all_at_once_dense(
                system,
                NULL,
                NULL);

            if(status != ROM_STROM_SUCCESS){
                free(dense_solution);
                free(monolis_solution);
                return status;
            }

            pack_reduced_coefficients(system, dense_solution);
        }

        status = ROM_std_strom_monolis_initialize(
            &monolis_solver,
            system,
            monolis_com,
            options->monolis_max_iter,
            options->monolis_tolerance);

        if(status != ROM_STROM_SUCCESS){
            free(dense_solution);
            free(monolis_solution);
            return status;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("reduced BCSR meta nodes: %d\n",
                monolis_solver.num_windows);
            printf("reduced BCSR block size: %d\n",
                monolis_solver.block_size);
            printf("reduced BCSR block nnz: %d\n",
                monolis_solver.num_graph_items);
            printf("reduced scalar dof: %d\n",
                monolis_solver.total_reduced_dof);
            printf("reduced scalar storage entries: %d\n",
                monolis_solver.num_graph_items
                * monolis_solver.block_size
                * monolis_solver.block_size);
            printf("MONOLIS solver: BiCGSTAB\n");
            printf("MONOLIS tolerance: %.15e\n",
                monolis_solver.tolerance);
            printf("MONOLIS max iterations: %d\n",
                monolis_solver.max_iter);
        }

        if(options->write_reduced_graph &&
           monolis_mpi_get_global_my_rank() == 0)
        {
            status = ROM_std_strom_write_bcsr_graph(
                &monolis_solver,
                options->reduced_graph_file);

            if(status != ROM_STROM_SUCCESS){
                ROM_std_strom_monolis_finalize(&monolis_solver);
                free(dense_solution);
                free(monolis_solution);
                return status;
            }

            printf("reduced BCSR graph: %s\n",
                options->reduced_graph_file);
        }

        status = ROM_std_strom_solve_all_at_once_monolis(
            system,
            &monolis_solver);

        if(status != ROM_STROM_SUCCESS){
            ROM_std_strom_monolis_finalize(&monolis_solver);
            free(dense_solution);
            free(monolis_solution);
            return status;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("MONOLIS reduced algebraic residual: %.15e\n",
                ROM_std_strom_reduced_residual_norm(system));
        }

        if(solver == ST_ROM_SOLVER_VERIFY_MONOLIS){
            double difference;

            pack_reduced_coefficients(system, monolis_solution);

            difference = relative_reduced_solution_difference(
                dense_solution,
                monolis_solution,
                system->total_reduced_dof);

            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "dense/MONOLIS reduced coefficient difference: "
                    "%.15e\n",
                    difference);
            }

            if(!isfinite(difference) ||
               difference > options->monolis_verify_tolerance)
            {
                fprintf(
                    stderr,
                    "ERROR: dense/MONOLIS coefficient difference "
                    "%.15e exceeds tolerance %.15e\n",
                    difference,
                    options->monolis_verify_tolerance);

                ROM_std_strom_monolis_finalize(&monolis_solver);
                free(dense_solution);
                free(monolis_solution);
                return ROM_STROM_ERR_SOLVER;
            }
        }

        ROM_std_strom_monolis_finalize(&monolis_solver);
        free(dense_solution);
        free(monolis_solution);

        return ROM_STROM_SUCCESS;
    }

    if(solver == ST_ROM_SOLVER_VERIFY){
        double* sequential_solution = (double*)calloc(
            (size_t)system->total_reduced_dof,
            sizeof(double));
        double* all_at_once_solution = (double*)calloc(
            (size_t)system->total_reduced_dof,
            sizeof(double));

        if(sequential_solution == NULL || all_at_once_solution == NULL){
            free(sequential_solution);
            free(all_at_once_solution);
            return ROM_STROM_ERR_ALLOCATION;
        }

        status = ROM_std_strom_solve_sequential(system, NULL, NULL);
        if(status != ROM_STROM_SUCCESS){
            free(sequential_solution);
            free(all_at_once_solution);
            return status;
        }

        pack_reduced_coefficients(system, sequential_solution);

        status = ROM_std_strom_solve_all_at_once_dense(system, NULL, NULL);
        if(status != ROM_STROM_SUCCESS){
            free(sequential_solution);
            free(all_at_once_solution);
            return status;
        }

        pack_reduced_coefficients(system, all_at_once_solution);

        printf(
            "sequential/all-at-once reduced coefficient difference: %.15e\n",
            relative_reduced_solution_difference(
                sequential_solution,
                all_at_once_solution,
                system->total_reduced_dof));

        free(sequential_solution);
        free(all_at_once_solution);
        return ROM_STROM_SUCCESS;
    }

    return ROM_STROM_ERR_ARGUMENT;
}

static void run_online(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows,
    const ST_ROM_OPTIONS* options)
{
    ROM_STROM_SYSTEM system;
    ROM_STROM_FOM_CONTEXT fom_context;
    ST_ROM_VALIDATION_WORK validation;
    HLPOD_VALUES* values = NULL;
    HLPOD_MAT* matrices = NULL;

    const int nx = sys->fe.total_num_nodes;
    const int previous_state_dof = checked_product3_int(nx, te->n_dof, 1);
    const int st_dof = checked_product3_int(
        nx,
        te->n_dof,
        num_window_slabs);

    double* initial_previous_state = NULL;
    double* previous_lift_state = NULL;
    double* current_lift = NULL;
    double* previous_lift = NULL;
    double* full_rhs = NULL;
    double* matrix_times_lift = NULL;
    double* homogeneous_rhs = NULL;
    double* reconstructed_homogeneous = NULL;
    double* reconstructed_full = NULL;
    int file_num = 0;
    int status;
    double stage_start;
    double operator_time = 0.0;
    double rhs_time = 0.0;
    double solve_time = 0.0;
    double reconstruction_time = 0.0;

    memset(&validation, 0, sizeof(validation));

    initialize_strom_system(
        &system,
        &values,
        &matrices,
        sys,
        te,
        num_window_slabs,
        num_windows,
        1,
        options->max_modes);

    status = ROM_std_strom_fom_read_shared_basis(
        options->basis_dir,
        &system.windows[0]);

    if(status != ROM_STROM_SUCCESS){
        fprintf(
            stderr,
            "ERROR: cannot read all-window shared basis from %s\n"
            "       expected file: all_windows_shared_pod.bin\n"
            "       rerun ST_ROM_MODE=offline with the shared-basis source.\n",
            options->basis_dir);
        exit(EXIT_FAILURE);
    }

    validate_shared_basis_manifest(
        options,
        sys,
        te,
        num_window_slabs,
        num_windows,
        system.windows[0].num_modes);

    for(int w = 1; w < num_windows; w++){
        status = ROM_std_strom_attach_pod_basis(
            &system.windows[w],
            &system.windows[0]);

        if(status != ROM_STROM_SUCCESS){
            fail_message("cannot attach all-window shared basis");
        }
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("online POD basis scope: all_windows_shared\n");
        printf("shared POD modes/window: %d\n",
            system.windows[0].num_modes);
    }

    status = ROM_std_strom_update_mode_offsets(&system);
    if(status != ROM_STROM_SUCCESS){
        fail_message("mode offset construction failed");
    }

    status = ROM_std_strom_fom_context_initialize(
        &fom_context,
        sys,
        te,
        num_window_slabs,
        num_windows);

    if(status != ROM_STROM_SUCCESS){
        fail_message("FOM adapter initialization failed");
    }

    /*
     * Prefer offline-projected D/L.  Rebuild only when requested or when the
     * file is absent/incompatible.  The rebuild path prepares A once per
     * block and uses BLAS for the dense projections.
     */
    stage_start = monolis_get_time();
    status = options->rebuild_operators
        ? ROM_STROM_ERR_IO
        : ROM_std_strom_fom_read_reduced_operators(
            options->basis_dir,
            &system,
            options->operator_reuse,
            sys->vals.dt,
            st_dof);

    if(status == ROM_STROM_SUCCESS){
        if(monolis_mpi_get_global_my_rank() == 0){
            printf("loaded reduced D/L operators from offline file.\n");
        }
    }
    else{
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "reduced operator file unavailable/incompatible; rebuilding.\n");
        }

        build_reduced_operators_fast(
            &system,
            &fom_context,
            options);

        status = ROM_std_strom_fom_write_reduced_operators(
            options->basis_dir,
            &system,
            options->operator_reuse,
            sys->vals.dt,
            st_dof);
        if(status != ROM_STROM_SUCCESS){
            fail_message("rebuilt reduced operator write failed");
        }
    }
    operator_time = monolis_get_time() - stage_start;

    initial_previous_state = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    previous_lift_state = (double*)calloc(
        (size_t)previous_state_dof,
        sizeof(double));

    current_lift = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    previous_lift = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    full_rhs = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    matrix_times_lift = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    homogeneous_rhs = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    reconstructed_homogeneous = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    reconstructed_full = (double*)calloc(
        (size_t)st_dof,
        sizeof(double));

    if(initial_previous_state == NULL ||
       previous_lift_state == NULL ||
       current_lift == NULL ||
       previous_lift == NULL ||
       full_rhs == NULL ||
       matrix_times_lift == NULL ||
       homogeneous_rhs == NULL ||
       reconstructed_homogeneous == NULL ||
       reconstructed_full == NULL)
    {
        fail_message("lifted online work allocation failed");
    }

    ST_TimeElement_set_prev_trace_from_nodal_value(
        te,
        nx,
        sys->vals.T,
        initial_previous_state);

    /*
     * U_w = U_D,w + q_w.
     *
     * For w=0:
     *   A_w q_w = b_w(U_initial) - A_w U_D,w.
     *
     * For w>0 the affine RHS contains the previous lift:
     *   A_w q_w = b_w(C_last U_D,w-1) - A_w U_D,w
     *             + B_w q_w-1.
     *
     * The final B_w q_w-1 term is supplied by L_w a_w-1 in the sequential
     * reduced solve.
     */
    stage_start = monolis_get_time();
    for(int w = 0; w < num_windows; w++){
        const double* previous_affine_state;

        status = ROM_std_strom_fom_build_dirichlet_lift(
            &fom_context,
            w,
            current_lift);

        if(status != ROM_STROM_SUCCESS){
            fail_message("current Dirichlet lift construction failed");
        }

        if(w == 0){
            previous_affine_state = initial_previous_state;
        }
        else{
            ST_Window_extract_last_slab_state(
                te,
                nx,
                num_window_slabs,
                previous_lift,
                previous_lift_state);

            previous_affine_state = previous_lift_state;
        }

        status = ROM_std_strom_fom_build_full_rhs(
            &fom_context,
            w,
            previous_affine_state,
            full_rhs);

        if(status != ROM_STROM_SUCCESS){
            fail_message("lifted full FOM RHS assembly failed");
        }

        /* Use the same prepared matrix that produced full_rhs. */
        status = ROM_std_strom_fom_apply_prepared_A(
            &fom_context,
            current_lift,
            matrix_times_lift);

        if(status != ROM_STROM_SUCCESS){
            fail_message("Dirichlet lift operator application failed");
        }

        for(int i = 0; i < st_dof; i++){
            homogeneous_rhs[i] =
                full_rhs[i] - matrix_times_lift[i];
        }

        status = ROM_std_strom_project_rhs(
            &system.windows[w],
            homogeneous_rhs);

        if(status != ROM_STROM_SUCCESS){
            fail_message("homogeneous reduced RHS projection failed");
        }

        memcpy(
            previous_lift,
            current_lift,
            (size_t)st_dof * sizeof(double));
    }
    rhs_time = monolis_get_time() - stage_start;

    stage_start = monolis_get_time();
    status = solve_reduced_system(
        &system,
        options->online_solver,
        &sys->monolis_com,
        options);

    if(status != ROM_STROM_SUCCESS){
        fail_message("lifted ST-ROM reduced solve failed");
    }
    solve_time = monolis_get_time() - stage_start;

    initialize_validation_work(
        &validation,
        options,
        sys,
        te,
        num_window_slabs,
        options->max_modes);

    stage_start = monolis_get_time();
    for(int w = 0; w < num_windows; w++){
        status = ROM_std_strom_reconstruct_window(
            &system.windows[w],
            reconstructed_homogeneous);

        if(status != ROM_STROM_SUCCESS){
            fail_message("homogeneous ST-ROM reconstruction failed");
        }

        status = ROM_std_strom_fom_build_dirichlet_lift(
            &fom_context,
            w,
            current_lift);

        if(status != ROM_STROM_SUCCESS){
            fail_message("online Dirichlet lift reconstruction failed");
        }

        for(int i = 0; i < st_dof; i++){
            reconstructed_full[i] =
                current_lift[i] + reconstructed_homogeneous[i];
        }

        validate_online_window(
            &validation,
            options,
            &fom_context,
            &system.windows[w],
            te,
            num_window_slabs,
            sys->vals.dt,
            current_lift,
            reconstructed_homogeneous,
            reconstructed_full);

        for(int slab = 0; slab < num_window_slabs; slab++){
            const int global_step =
                w * num_window_slabs + slab + 1;

            const double time =
                (double)global_step * sys->vals.dt;

            if(options->write_output &&
               global_step % sys->vals.output_interval == 0){
                ST_Window_extract_slab_right_trace_to_nodal_value(
                    te,
                    nx,
                    num_window_slabs,
                    slab,
                    reconstructed_full,
                    sys->vals.T);

                output_files(sys, file_num, time);
                file_num++;
            }
        }
    }

    reconstruction_time = monolis_get_time() - stage_start;

    ST_Window_extract_right_trace_to_nodal_value(
        te,
        nx,
        num_window_slabs,
        reconstructed_full,
        sys->vals.T);

    write_validation_summary(
        &validation,
        sys,
        &system);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "Lifted ST-ROM online solve completed: windows=%d, reduced dof=%d\n",
            num_windows,
            system.total_reduced_dof);
        printf(
            "ST-ROM timings [s]: operators=%.6f rhs=%.6f solve=%.6f "
            "reconstruct+validation+output=%.6f\n",
            operator_time,
            rhs_time,
            solve_time,
            reconstruction_time);

        if(validation.enabled){
            printf("validation CSV: %s\n", options->validation_csv);
        }
    }

    finalize_validation_work(&validation);

    free(initial_previous_state);
    free(previous_lift_state);
    free(current_lift);
    free(previous_lift);
    free(full_rhs);
    free(matrix_times_lift);
    free(homogeneous_rhs);
    free(reconstructed_homogeneous);
    free(reconstructed_full);

    ROM_std_strom_fom_context_finalize(&fom_context);
    finalize_strom_system(&system, &values, &matrices);
}

int main(int argc, char* argv[])
{
    FE_SYSTEM sys = {0};
    ST_TIME_ELEMENT te;
    ST_ROM_OPTIONS options;
    int num_windows = 0;
    int my_rank;
    double time_start;

    monolis_global_initialize();

    require_serial_runtime();

    my_rank = monolis_mpi_get_global_my_rank();
    time_start = monolis_get_time();

    initialize_fom(
        &sys,
        &te,
        argc,
        argv,
        &num_windows);

    read_rom_options(
        &options,
        sys.cond.directory);

    if(options.mode == ST_RUN_FOM ||
       options.mode == ST_RUN_COLLECT ||
       options.mode == ST_RUN_ONLINE)
    {
        initialize_output(&sys, my_rank);
    }

    if(my_rank == 0){
        printf("ST run mode: %s\n", mode_name(options.mode));
        printf("time dG degree: %d\n", ST_ROM_DEGREE);
        printf("slabs/window: %d\n", ST_ROM_WINDOW_SLABS);
        printf("num windows: %d\n", num_windows);
        if(options.mode == ST_RUN_OFFLINE ||
           options.mode == ST_RUN_ONLINE)
        {
            printf("POD basis scope: all_windows_shared\n");
        }
        if(options.mode == ST_RUN_OFFLINE ||
           options.mode == ST_RUN_ONLINE)
        {
            printf("operator reuse: %d\n", options.operator_reuse);
            printf("direct coupling: %d\n", options.direct_coupling);
        }
        if(options.mode == ST_RUN_ONLINE){
            printf("online reduced solver: %s\n",
                online_solver_name(options.online_solver));
            printf("write standard output: %d\n", options.write_output);
        }
    }

    switch(options.mode){
    case ST_RUN_FOM:
        run_fom_or_collect(
            &sys,
            &te,
            ST_ROM_WINDOW_SLABS,
            num_windows,
            &options,
            1);
        break;

    case ST_RUN_COLLECT:
        run_fom_or_collect(
            &sys,
            &te,
            ST_ROM_WINDOW_SLABS,
            num_windows,
            &options,
            1);
        break;

    case ST_RUN_OFFLINE:
        run_offline(
            &sys,
            &te,
            ST_ROM_WINDOW_SLABS,
            num_windows,
            &options);
        break;

    case ST_RUN_ONLINE:
        run_online(
            &sys,
            &te,
            ST_ROM_WINDOW_SLABS,
            num_windows,
            &options);
        break;

    default:
        fail_message("unreachable ST-ROM mode");
    }

    finalize_fom(&sys);

    if(my_rank == 0){
        printf(
            "** Total time: %f\n\n",
            monolis_get_time() - time_start);
    }

    monolis_global_finalize();
    return 0;
}
