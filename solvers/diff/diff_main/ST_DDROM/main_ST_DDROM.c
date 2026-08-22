#include "core_ROM.h"
#include "rom_std_stddrom.h"
#include "rom_std_stddrom_diffusion.h"
#include "st_fom_common.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    DIFF_STDD_MODE_OFFLINE = 0,
    DIFF_STDD_MODE_ONLINE = 1,
    DIFF_STDD_MODE_FULL = 2
} DIFF_STDD_MODE;

typedef enum {
    DIFF_STDD_SOLVER_ALL_AT_ONCE = 0,
    DIFF_STDD_SOLVER_WINDOW_MARCHING = 1,
    DIFF_STDD_SOLVER_COMPARE = 2
} DIFF_STDD_SOLVER_STRATEGY;

typedef struct {
    double snapshot_collection;
    double offline_pod;
    double reduced_operator_assembly;
    double offline_io;
    double online_data_load;
    double rhs_projection;
    double all_at_once_system_assembly;
    double all_at_once_solve;
    double window_marching_solve;
    double reconstruction;
    double reconstruction_output_stage;
    double total_stddrom;
} DIFF_STDD_TIMING;

static void diff_stdd_fail(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static int diff_stdd_env_int(const char* name, int default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    long value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }

    errno = 0;
    value = strtol(text, &end, 10);
    if(errno != 0 || end == text || *end != '\0' ||
       value < INT_MIN || value > INT_MAX)
    {
        fprintf(stderr, "ERROR: invalid integer %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return (int)value;
}

static double diff_stdd_env_double(const char* name, double default_value)
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

static const char* diff_stdd_env_string(
    const char* name,
    const char* default_value)
{
    const char* text = getenv(name);
    return text != NULL && text[0] != '\0' ? text : default_value;
}

static DIFF_STDD_MODE diff_stdd_read_mode(void)
{
    const char* text = diff_stdd_env_string("DIFF_ST_DDROM_MODE", "full");

    if(strcmp(text, "offline") == 0){
        return DIFF_STDD_MODE_OFFLINE;
    }
    if(strcmp(text, "online") == 0){
        return DIFF_STDD_MODE_ONLINE;
    }
    if(strcmp(text, "full") == 0){
        return DIFF_STDD_MODE_FULL;
    }

    fprintf(
        stderr,
        "ERROR: DIFF_ST_DDROM_MODE must be offline, online, or full; got %s\n",
        text);
    exit(EXIT_FAILURE);
}

static DIFF_STDD_SOLVER_STRATEGY diff_stdd_read_solver_strategy(void)
{
    const char* text = diff_stdd_env_string(
        "DIFF_ST_DDROM_SOLVER",
        "all-at-once");

    if(strcmp(text, "all-at-once") == 0 || strcmp(text, "aoa") == 0){
        return DIFF_STDD_SOLVER_ALL_AT_ONCE;
    }
    if(strcmp(text, "window-marching") == 0 ||
       strcmp(text, "marching") == 0){
        return DIFF_STDD_SOLVER_WINDOW_MARCHING;
    }
    if(strcmp(text, "compare") == 0 || strcmp(text, "both") == 0){
        return DIFF_STDD_SOLVER_COMPARE;
    }

    fprintf(
        stderr,
        "ERROR: DIFF_ST_DDROM_SOLVER must be all-at-once, "
        "window-marching, or compare; got %s\n",
        text);
    exit(EXIT_FAILURE);
}

static const char* diff_stdd_solver_strategy_name(
    DIFF_STDD_SOLVER_STRATEGY strategy)
{
    switch(strategy){
        case DIFF_STDD_SOLVER_ALL_AT_ONCE:
            return "all-at-once";
        case DIFF_STDD_SOLVER_WINDOW_MARCHING:
            return "window-marching";
        case DIFF_STDD_SOLVER_COMPARE:
            return "compare";
        default:
            return "unknown";
    }
}

static int diff_stdd_total_steps(double finish_time, double dt)
{
    const double value = finish_time / dt;
    const long long rounded = llround(value);

    if(dt <= 0.0 || finish_time <= 0.0 || rounded <= 0 ||
       rounded > INT_MAX ||
       fabs(value - (double)rounded) > 1.0e-10 * fmax(1.0, fabs(value)))
    {
        diff_stdd_fail("finish_time/dt must be a positive integer");
    }

    return (int)rounded;
}

static int diff_stdd_join_path(
    char* output,
    size_t capacity,
    const char* directory,
    const char* leaf)
{
    int written;

    if(output == NULL || capacity == 0u || directory == NULL || leaf == NULL){
        return -1;
    }

    if(strcmp(directory, ".") == 0){
        written = snprintf(output, capacity, "./%s", leaf);
    }
    else{
        const size_t length = strlen(directory);
        const int has_slash = length > 0u && directory[length - 1u] == '/';
        written = snprintf(
            output,
            capacity,
            has_slash ? "%s%s" : "%s/%s",
            directory,
            leaf);
    }

    return written < 0 || (size_t)written >= capacity ? -1 : 0;
}

static int diff_stdd_prepare_fom(FE_SYSTEM* sys, int argc, char* argv[])
{
    int manufactured_rank;

    memset(sys, 0, sizeof(*sys));

    sys->cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
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

    manufactured_rank = diff_stdd_env_int("DIFF_MANUFACTURED_RANK", 0);
    if(manufactured_rank > 0){
        const int total_steps =
            diff_stdd_total_steps(sys->vals.finish_time, sys->vals.dt);
        const int slabs_per_window =
            diff_stdd_env_int("DIFF_ST_DDROM_SLABS_PER_WINDOW", 4);
        const double diffusivity =
            diff_stdd_env_double("DIFF_MANUFACTURED_KAPPA", 1.0);
        const double total_amplitude =
            diff_stdd_env_double("DIFF_MANUFACTURED_AMPLITUDE", 1.0);
        double domain_min[3];
        double domain_max[3];
        double negative_min[3];

        if(slabs_per_window <= 0 ||
           total_steps % slabs_per_window != 0 ||
           diffusivity <= 0.0 || total_amplitude <= 0.0)
        {
            diff_stdd_fail("invalid rank-controlled manufactured configuration");
        }

        for(int d = 0; d < 3; d++){
            domain_min[d] = sys->fe.x[0][d];
            domain_max[d] = sys->fe.x[0][d];
            for(int node = 1; node < sys->fe.total_num_nodes; node++){
                if(sys->fe.x[node][d] < domain_min[d]){
                    domain_min[d] = sys->fe.x[node][d];
                }
                if(sys->fe.x[node][d] > domain_max[d]){
                    domain_max[d] = sys->fe.x[node][d];
                }
            }
            negative_min[d] = -domain_min[d];
        }

        monolis_allreduce_R(
            3, domain_max, MONOLIS_MPI_MAX, sys->monolis_com.comm);
        monolis_allreduce_R(
            3, negative_min, MONOLIS_MPI_MAX, sys->monolis_com.comm);
        for(int d = 0; d < 3; d++){
            domain_min[d] = -negative_min[d];
        }

        if(manusol_configure_rank_benchmark(
            manufactured_rank,
            slabs_per_window,
            total_steps / slabs_per_window,
            sys->vals.dt,
            diffusivity,
            total_amplitude,
            domain_min,
            domain_max) != 0)
        {
            diff_stdd_fail(
                "DIFF_MANUFACTURED_RANK must be in [1,64] and the mesh "
                "must have positive extent in x/y/z");
        }
    }
    else{
        const double dummy_min[3] = {0.0, 0.0, 0.0};
        const double dummy_max[3] = {1.0, 1.0, 1.0};
        (void)manusol_configure_rank_benchmark(
            0, 1, 1, 1.0, 1.0, 1.0, dummy_min, dummy_max);
    }

    memory_allocation_nodal_values(&sys->vals, sys->fe.total_num_nodes);
    manusol_set_init_value(&sys->fe, sys->vals.T);

    BBFE_elemmat_set_Jacobi_mat(&sys->fe, &sys->basis);
    BBFE_elemmat_set_shapefunc_derivative(&sys->fe, &sys->basis);

    monolis_initialize(&sys->monolis0);
    monolis_get_nonzero_pattern_by_simple_mesh_R(
        &sys->monolis0,
        sys->fe.total_num_nodes,
        sys->fe.local_num_nodes,
        1,
        sys->fe.total_num_elems,
        sys->fe.conn);

    set_element_mat(&sys->monolis0, &sys->fe, &sys->basis, &sys->vals);
    monolis_copy_mat_R(&sys->monolis0, &sys->monolis);

    return 0;
}

static void diff_stdd_finalize_fom(FE_SYSTEM* sys)
{
    BBFE_convdiff_finalize(&sys->fe, &sys->basis, &sys->bc);
    monolis_finalize(&sys->monolis);
    monolis_finalize(&sys->monolis0);
}

static int diff_stdd_collect_training_snapshots(
    FE_SYSTEM* sys,
    const ST_FOM_WINDOW_LAYOUT* layout,
    int num_training_windows,
    double* snapshot_matrix)
{
    const int ns = layout->num_slabs;
    const int owned = layout->num_owned_space_points;
    const int columns = num_training_windows;
    int step = 0;

    if(sys == NULL || layout == NULL || snapshot_matrix == NULL ||
       num_training_windows <= 0 || layout->time.degree != 0)
    {
        return -1;
    }

    manusol_set_init_value(&sys->fe, sys->vals.T);

    if(manusol_rank_benchmark_is_active()){
        /*
         * For the rank-controlled benchmark, collect the prescribed discrete
         * trajectory directly.  The online RHS is manufactured with the same
         * ST operator, so this is also the exact discrete FOM trajectory.
         */
        for(int window = 0; window < num_training_windows; window++){
            for(int slab = 0; slab < ns; slab++){
                const int global_step = window * ns + slab + 1;
                const double time = (double)global_step * sys->vals.dt;

                for(int node = 0; node < owned; node++){
                    const int row = node * ns + slab;
                    double lift = 0.0;
                    const double target = manusol_get_sol(
                        sys->fe.x[node][0],
                        sys->fe.x[node][1],
                        sys->fe.x[node][2],
                        time);

                    if(sys->bc.D_bc_exists[node]){
                        lift = target;
                    }

                    snapshot_matrix[
                        (size_t)row * (size_t)columns + (size_t)window] =
                        target - lift;
                }
            }
        }
        return 0;
    }

    for(int window = 0; window < num_training_windows; window++){
        for(int slab = 0; slab < ns; slab++){
            double time;

            step++;
            time = (double)step * sys->vals.dt;
            solver_fom(*sys, time, step);

            for(int node = 0; node < owned; node++){
                const int row = node * ns + slab;
                double lift = 0.0;

                if(sys->bc.D_bc_exists[node]){
                    lift = manusol_get_sol(
                        sys->fe.x[node][0],
                        sys->fe.x[node][1],
                        sys->fe.x[node][2],
                        time);
                }

                snapshot_matrix[(size_t)row * (size_t)columns + (size_t)window] =
                    sys->vals.T[node] - lift;
            }
        }
    }

    return 0;
}

static int diff_stdd_write_coefficients(
    const ROM_STDD_SYSTEM* reduced,
    const char* filename)
{
    FILE* fp;

    if(reduced == NULL || filename == NULL || reduced->reduced_coef == NULL){
        return -1;
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        return -1;
    }

    fprintf(fp, "window");
    for(int mode = 0; mode < reduced->num_modes; mode++){
        fprintf(fp, ",a_%d", mode);
    }
    fputc('\n', fp);

    for(int window = 0; window < reduced->num_windows; window++){
        fprintf(fp, "%d", window);
        for(int mode = 0; mode < reduced->num_modes; mode++){
            fprintf(
                fp,
                ",%.17e",
                reduced->reduced_coef[
                    (size_t)window * (size_t)reduced->num_modes_capacity
                    + (size_t)mode]);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    return 0;
}

static int diff_stdd_output_solution(
    ROM_STDD_DIFFUSION_CONTEXT* adapter,
    int output_interval,
    int write_output,
    const char* error_filename,
    double* reconstruction_time)
{
    FE_SYSTEM* sys = adapter->fom;
    ROM_STDD_SYSTEM* reduced = adapter->reduced;
    const int ns = adapter->layout->num_slabs;
    double* window_state = NULL;
    FILE* fp = NULL;
    int file_number = 0;
    int step = 0;
    double local_reconstruction_time = 0.0;

    window_state = (double*)calloc(
        (size_t)reduced->local_st_dof,
        sizeof(double));
    if(window_state == NULL){
        return -1;
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        fp = fopen(error_filename, "w");
        if(fp == NULL){
            free(window_state);
            return -1;
        }
        fprintf(fp, "# time relative_L2_error\n");
    }

    for(int window = 0; window < reduced->num_windows; window++){
        const double reconstruction_t0 = monolis_get_time();
        const int reconstruction_status =
            ROM_std_stdd_diffusion_reconstruct_window(
                adapter,
                window,
                window_state);
        local_reconstruction_time +=
            monolis_get_time() - reconstruction_t0;

        if(reconstruction_status != ROM_STDD_SUCCESS)
        {
            if(fp != NULL){
                fclose(fp);
            }
            free(window_state);
            return -1;
        }

        for(int slab = 0; slab < ns; slab++){
            double time;
            double l2_error;

            step++;
            time = (double)step * sys->vals.dt;

            for(int node = 0; node < sys->fe.total_num_nodes; node++){
                sys->vals.T[node] = window_state[
                    (size_t)node * (size_t)ns + (size_t)slab];
            }

            l2_error = BBFE_convdiff_equivval_relative_L2_error_scalar(
                &sys->fe,
                &sys->basis,
                &sys->monolis_com,
                time,
                sys->vals.T,
                manusol_get_sol);

            if(monolis_mpi_get_global_my_rank() == 0){
                fprintf(fp, "%.17e %.17e\n", time, l2_error);
                printf(
                    "diffusion ST-DDROM step=%d time=%.8e relative L2=%.8e\n",
                    step,
                    time,
                    l2_error);
            }

            if(write_output && output_interval > 0 && step % output_interval == 0){
                char local_name[256];
                const char* output_name;
                int written = snprintf(
                    local_name,
                    sizeof(local_name),
                    "stddrom_result_%06d.vtk",
                    file_number);

                if(written < 0 || (size_t)written >= sizeof(local_name)){
                    if(fp != NULL){
                        fclose(fp);
                    }
                    free(window_state);
                    return -1;
                }

                output_name = monolis_get_global_output_file_name(
                    MONOLIS_DEFAULT_TOP_DIR,
                    "./",
                    local_name);

                output_result_file_vtk(
                    &sys->fe,
                    &sys->vals,
                    output_name,
                    sys->cond.directory,
                    time);
                file_number++;
            }
        }
    }

    if(fp != NULL){
        fclose(fp);
    }

    monolis_allreduce_R(
        1,
        &local_reconstruction_time,
        MONOLIS_MPI_MAX,
        sys->monolis_com.comm);
    if(reconstruction_time != NULL){
        *reconstruction_time = local_reconstruction_time;
    }

    free(window_state);
    return 0;
}

static double diff_stdd_coefficient_relative_difference(
    const ROM_STDD_SYSTEM* reduced,
    const double* first,
    const double* second)
{
    if(reduced == NULL || first == NULL || second == NULL){
        return -1.0;
    }

    return ROM_std_stdd_global_relative_error(
        first,
        second,
        reduced->num_windows * reduced->num_modes_capacity,
        monolis_mpi_get_global_comm());
}

static int diff_stdd_write_timing_csv(
    const char* filename,
    DIFF_STDD_MODE mode,
    DIFF_STDD_SOLVER_STRATEGY strategy,
    int mpi_size,
    int total_steps,
    int slabs_per_window,
    int num_windows,
    int num_modes,
    const DIFF_STDD_TIMING* timing,
    double coefficient_difference,
    double all_at_once_residual,
    double window_marching_residual)
{
    FILE* fp;
    FILE* probe;
    int write_header = 0;

    if(filename == NULL || timing == NULL){
        return -1;
    }

    probe = fopen(filename, "r");
    if(probe == NULL){
        write_header = 1;
    }
    else{
        fclose(probe);
    }

    fp = fopen(filename, "a");
    if(fp == NULL){
        return -1;
    }

    if(write_header){
        fprintf(
            fp,
            "mode,solver,mpi_size,total_steps,slabs_per_window,num_windows,"
            "manufactured_rank,num_modes,snapshot_collection,offline_pod,"
            "reduced_operator_assembly,offline_io,online_data_load,"
            "rhs_projection,all_at_once_system_assembly,"
            "all_at_once_solve,window_marching_solve,reconstruction,"
            "reconstruction_output_stage,total_stddrom,"
            "coefficient_relative_difference,all_at_once_residual,"
            "window_marching_residual\n");
    }

    fprintf(
        fp,
        "%s,%s,%d,%d,%d,%d,%d,%d,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
        mode == DIFF_STDD_MODE_OFFLINE
            ? "offline"
            : (mode == DIFF_STDD_MODE_ONLINE ? "online" : "full"),
        diff_stdd_solver_strategy_name(strategy),
        mpi_size,
        total_steps,
        slabs_per_window,
        num_windows,
        manusol_rank_benchmark_rank(),
        num_modes,
        timing->snapshot_collection,
        timing->offline_pod,
        timing->reduced_operator_assembly,
        timing->offline_io,
        timing->online_data_load,
        timing->rhs_projection,
        timing->all_at_once_system_assembly,
        timing->all_at_once_solve,
        timing->window_marching_solve,
        timing->reconstruction,
        timing->reconstruction_output_stage,
        timing->total_stddrom,
        coefficient_difference,
        all_at_once_residual,
        window_marching_residual);

    fclose(fp);
    return 0;
}

static void diff_stdd_print_timing_summary(
    DIFF_STDD_SOLVER_STRATEGY strategy,
    const DIFF_STDD_TIMING* timing,
    double coefficient_difference,
    double all_at_once_residual,
    double window_marching_residual)
{
    if(monolis_mpi_get_global_my_rank() != 0 || timing == NULL){
        return;
    }

    printf("\nST-DDROM timing summary [s]\n");
    printf("  Snapshot collection time          : %.9e\n", timing->snapshot_collection);
    printf("  Offline POD time                  : %.9e\n", timing->offline_pod);
    printf("  Reduced operator assembly time    : %.9e\n", timing->reduced_operator_assembly);
    printf("  Reduced RHS projection time       : %.9e\n", timing->rhs_projection);
    printf("  All-at-once system assembly time  : %.9e\n", timing->all_at_once_system_assembly);
    printf("  All-at-once solve time            : %.9e\n", timing->all_at_once_solve);
    printf("  Window-marching solve time        : %.9e\n", timing->window_marching_solve);
    printf("  Reconstruction time               : %.9e\n", timing->reconstruction);
    printf("  Reconstruction/output stage time  : %.9e\n", timing->reconstruction_output_stage);
    printf("  Total ST-DDROM time               : %.9e\n", timing->total_stddrom);

    if(strategy == DIFF_STDD_SOLVER_COMPARE){
        const double aoa_path =
            timing->all_at_once_system_assembly + timing->all_at_once_solve;
        const double marching_path = timing->window_marching_solve;

        printf("\nST-DDROM solver comparison\n");
        printf("  Coefficient relative difference   : %.15e\n", coefficient_difference);
        printf("  All-at-once reduced residual      : %.15e\n", all_at_once_residual);
        printf("  Window-marching reduced residual  : %.15e\n", window_marching_residual);
        if(timing->all_at_once_solve > 0.0){
            printf(
                "  Pure solve speedup (marching/AOA) : %.9e\n",
                timing->window_marching_solve / timing->all_at_once_solve);
        }
        if(aoa_path > 0.0){
            printf(
                "  Solve-path speedup (marching/AOA) : %.9e\n",
                marching_path / aoa_path);
        }
    }
    fflush(stdout);
}

int main(int argc, char* argv[])
{
    FE_SYSTEM sys;
    ST_FOM_WINDOW_LAYOUT layout;
    ROM_STDD_SYSTEM reduced;
    ROM_STDD_DIFFUSION_CONTEXT adapter;
    ROM_STDD_REDUCED_SOLVER reduced_solver;
    DIFF_STDD_MODE mode;
    DIFF_STDD_SOLVER_STRATEGY solver_strategy;
    DIFF_STDD_TIMING timing;

    int total_steps;
    int slabs_per_window;
    int num_windows;
    int num_training_windows;
    int owned_nodes;
    int max_modes;
    int write_output;
    int linear_max_iter;
    int reduced_solver_initialized = 0;
    double pod_energy;
    double linear_tolerance;
    double total_t0 = 0.0;
    double coefficient_difference = -1.0;
    double all_at_once_residual = -1.0;
    double window_marching_residual = -1.0;

    double* snapshot_matrix = NULL;
    double* all_at_once_coefficients = NULL;
    char basis_directory[4096];
    char operator_directory[4096];
    char coefficient_filename[4096];
    char error_filename[4096];
    char timing_filename[4096];
    int status = EXIT_FAILURE;

    memset(&sys, 0, sizeof(sys));
    memset(&layout, 0, sizeof(layout));
    memset(&reduced, 0, sizeof(reduced));
    memset(&adapter, 0, sizeof(adapter));
    memset(&reduced_solver, 0, sizeof(reduced_solver));
    memset(&timing, 0, sizeof(timing));

    monolis_global_initialize();

    if(diff_stdd_prepare_fom(&sys, argc, argv) != 0){
        diff_stdd_fail("diffusion FOM initialization failed");
    }

    total_steps = diff_stdd_total_steps(sys.vals.finish_time, sys.vals.dt);
    slabs_per_window = diff_stdd_env_int("DIFF_ST_DDROM_SLABS_PER_WINDOW", 4);

    if(slabs_per_window <= 0 || total_steps % slabs_per_window != 0){
        diff_stdd_fail(
            "DIFF_ST_DDROM_SLABS_PER_WINDOW must be positive and divide total steps");
    }

    num_windows = total_steps / slabs_per_window;

    if(ST_fom_resolve_owned_space_points(
        sys.fe.total_num_nodes,
        monolis_mpi_get_global_comm_size(),
        sys.monolis_com.n_internal_vertex,
        &owned_nodes) != ST_FOM_SUCCESS)
    {
        diff_stdd_fail("cannot resolve owned spatial nodes");
    }

    if(ST_fom_layout_initialize(
        &layout,
        sys.fe.total_num_nodes,
        owned_nodes,
        1,
        slabs_per_window,
        0) != ST_FOM_SUCCESS)
    {
        diff_stdd_fail("dG(0) space-time layout initialization failed");
    }

    max_modes = diff_stdd_env_int("DIFF_ST_DDROM_MAX_MODES", 10);
    pod_energy = diff_stdd_env_double("DIFF_ST_DDROM_POD_ENERGY", 0.999999);
    linear_max_iter = diff_stdd_env_int(
        "DIFF_ST_DDROM_LINEAR_MAX_ITER",
        sys.vals.mat_max_iter);
    linear_tolerance = diff_stdd_env_double(
        "DIFF_ST_DDROM_LINEAR_EPSILON",
        sys.vals.mat_epsilon);
    write_output = diff_stdd_env_int("DIFF_ST_DDROM_WRITE_OUTPUT", 1);

    if(max_modes <= 0 || pod_energy <= 0.0 || pod_energy > 1.0 ||
       linear_max_iter <= 0 || linear_tolerance <= 0.0 ||
       (write_output != 0 && write_output != 1))
    {
        diff_stdd_fail("invalid ST-DDROM environment parameter");
    }
    if(manusol_rank_benchmark_is_active() &&
       max_modes < manusol_rank_benchmark_rank())
    {
        diff_stdd_fail(
            "rank-controlled benchmark requires "
            "DIFF_ST_DDROM_MAX_MODES >= DIFF_MANUFACTURED_RANK");
    }

    num_training_windows = diff_stdd_env_int(
        "DIFF_ST_DDROM_TRAINING_WINDOWS",
        num_windows);
    if(num_training_windows <= 0 || num_training_windows > num_windows){
        diff_stdd_fail("DIFF_ST_DDROM_TRAINING_WINDOWS must be in [1,num_windows]");
    }
    if(manusol_rank_benchmark_is_active() &&
       num_training_windows < manusol_rank_benchmark_rank())
    {
        diff_stdd_fail(
            "rank-controlled benchmark requires "
            "DIFF_ST_DDROM_TRAINING_WINDOWS >= DIFF_MANUFACTURED_RANK");
    }

    if(ROM_std_stdd_system_initialize(
        &reduced,
        &sys.monolis_com,
        num_windows,
        sys.fe.total_num_nodes,
        owned_nodes,
        slabs_per_window,
        max_modes) != ROM_STDD_SUCCESS)
    {
        diff_stdd_fail("ROM_STDD_SYSTEM initialization failed");
    }

    if(ROM_std_stdd_diffusion_context_initialize(
        &adapter,
        &sys,
        &layout,
        &reduced) != ROM_STDD_SUCCESS)
    {
        diff_stdd_fail("diffusion ST-DDROM adapter initialization failed");
    }

    if(diff_stdd_join_path(
        basis_directory,
        sizeof(basis_directory),
        sys.cond.directory,
        diff_stdd_env_string("DIFF_ST_DDROM_BASIS_SUBDIR", "stddrom_basis")) != 0 ||
       diff_stdd_join_path(
        operator_directory,
        sizeof(operator_directory),
        sys.cond.directory,
        diff_stdd_env_string("DIFF_ST_DDROM_OPERATOR_SUBDIR", "stddrom_operators")) != 0 ||
       diff_stdd_join_path(
        coefficient_filename,
        sizeof(coefficient_filename),
        sys.cond.directory,
        "stddrom_coefficients.csv") != 0 ||
       diff_stdd_join_path(
        error_filename,
        sizeof(error_filename),
        sys.cond.directory,
        "stddrom_l2_error.csv") != 0 ||
       diff_stdd_join_path(
        timing_filename,
        sizeof(timing_filename),
        sys.cond.directory,
        "stddrom_timing.csv") != 0)
    {
        diff_stdd_fail("ST-DDROM output path too long");
    }

    mode = diff_stdd_read_mode();
    solver_strategy = diff_stdd_read_solver_strategy();
    total_t0 = monolis_get_time_global_sync();

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "\nDiffusion ST-DDROM\n"
            "  mode             : %s\n"
            "  time scheme      : dG(0)\n"
            "  total steps      : %d\n"
            "  slabs/window     : %d\n"
            "  windows          : %d\n"
            "  training windows : %d\n"
            "  local/owned nodes: %d/%d\n"
            "  max modes/rank   : %d\n"
            "  POD retained E  : %.8e\n"
            "  basis directory  : %s\n"
            "  operator dir     : %s\n"
            "  solver strategy  : %s\n"
            "  manufactured rank: %d%s\n",
            mode == DIFF_STDD_MODE_OFFLINE
                ? "offline"
                : (mode == DIFF_STDD_MODE_ONLINE ? "online" : "full"),
            total_steps,
            slabs_per_window,
            num_windows,
            num_training_windows,
            sys.fe.total_num_nodes,
            owned_nodes,
            max_modes,
            pod_energy,
            basis_directory,
            operator_directory,
            diff_stdd_solver_strategy_name(solver_strategy),
            manusol_rank_benchmark_rank(),
            manusol_rank_benchmark_is_active()
                ? " (rank-controlled dG(0) benchmark)"
                : " (legacy manufactured case)");

        if(manusol_rank_benchmark_is_active()){
            printf(
                "  benchmark kappa : %.8e\n"
                "  benchmark amp   : %.8e (total modal amplitude)\n",
                manusol_rank_benchmark_diffusivity(),
                manusol_rank_benchmark_total_amplitude());
        }
    }

    if(mode == DIFF_STDD_MODE_OFFLINE || mode == DIFF_STDD_MODE_FULL){
        snapshot_matrix = (double*)calloc(
            (size_t)reduced.internal_st_dof * (size_t)num_training_windows,
            sizeof(double));
        if(snapshot_matrix == NULL){
            diff_stdd_fail("snapshot matrix allocation failed");
        }

        {
            const double t0 = monolis_get_time_global_sync();
            if(diff_stdd_collect_training_snapshots(
                &sys,
                &layout,
                num_training_windows,
                snapshot_matrix) != 0)
            {
                diff_stdd_fail("FOM training snapshot collection failed");
            }
            timing.snapshot_collection =
                monolis_get_time_global_sync() - t0;
        }

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_compute_shared_local_pod(
                &reduced,
                snapshot_matrix,
                num_training_windows,
                pod_energy,
                sys.monolis_com.comm) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("local space-time POD failed");
            }

            if(ROM_std_stdd_diffusion_enforce_homogeneous_basis(&adapter)
               != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("Dirichlet-homogeneous basis enforcement failed");
            }
            timing.offline_pod = monolis_get_time_global_sync() - t0;
        }

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_project_operators_rank_sweep(
                &reduced,
                ROM_std_stdd_diffusion_apply_A,
                ROM_std_stdd_diffusion_apply_B,
                &adapter) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("diffusion reduced-operator projection failed");
            }
            timing.reduced_operator_assembly =
                monolis_get_time_global_sync() - t0;
        }

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_write_basis_rank(
                &reduced,
                num_training_windows,
                basis_directory) != ROM_STDD_SUCCESS ||
               ROM_std_stdd_write_operators_rank(
                &reduced,
                operator_directory) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("ST-DDROM offline data write failed");
            }
            timing.offline_io = monolis_get_time_global_sync() - t0;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "offline projection complete: modes/rank=%d\n",
                reduced.num_modes);
        }
    }

    free(snapshot_matrix);
    snapshot_matrix = NULL;

    if(mode == DIFF_STDD_MODE_OFFLINE){
        timing.total_stddrom = monolis_get_time_global_sync() - total_t0;
        diff_stdd_print_timing_summary(
            solver_strategy,
            &timing,
            coefficient_difference,
            all_at_once_residual,
            window_marching_residual);
        if(monolis_mpi_get_global_my_rank() == 0){
            if(diff_stdd_write_timing_csv(
                timing_filename,
                mode,
                solver_strategy,
                monolis_mpi_get_global_comm_size(),
                total_steps,
                slabs_per_window,
                num_windows,
                reduced.num_modes,
                &timing,
                coefficient_difference,
                all_at_once_residual,
                window_marching_residual) != 0)
            {
                fprintf(stderr, "WARNING: cannot write %s\n", timing_filename);
            }
        }
        status = EXIT_SUCCESS;
        goto cleanup;
    }

    if(mode == DIFF_STDD_MODE_ONLINE){
        const double t0 = monolis_get_time_global_sync();
        if(ROM_std_stdd_read_basis_rank(&reduced, basis_directory)
           != ROM_STDD_SUCCESS)
        {
            diff_stdd_fail("cannot read diffusion ST-DDROM basis");
        }

        /* The stored basis was homogenized before operator projection.
         * Do not re-orthogonalize it here: the stored reduced operators must
         * remain paired with exactly the basis that generated them. */
        if(ROM_std_stdd_read_operators_rank(&reduced, operator_directory)
           != ROM_STDD_SUCCESS)
        {
            diff_stdd_fail("cannot read diffusion ST-DDROM operators");
        }
        timing.online_data_load = monolis_get_time_global_sync() - t0;
    }

    {
        const double t0 = monolis_get_time_global_sync();
        if(ROM_std_stdd_diffusion_project_all_rhs(&adapter) != ROM_STDD_SUCCESS){
            diff_stdd_fail("all-window reduced RHS projection failed");
        }
        timing.rhs_projection = monolis_get_time_global_sync() - t0;
    }

    if(solver_strategy == DIFF_STDD_SOLVER_ALL_AT_ONCE ||
       solver_strategy == DIFF_STDD_SOLVER_COMPARE)
    {
        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_reduced_solver_initialize(
                &reduced_solver,
                &reduced,
                &sys.monolis_com,
                linear_max_iter,
                linear_tolerance) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("all-at-once reduced solver initialization failed");
            }
            reduced_solver_initialized = 1;

            if(ROM_std_stdd_reduced_solver_assemble(&reduced_solver, &reduced)
               != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("all-at-once reduced system assembly failed");
            }
            timing.all_at_once_system_assembly =
                monolis_get_time_global_sync() - t0;
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "solving all windows simultaneously: windows=%d modes/rank=%d "
                "global reduced dof=%d\n",
                reduced.num_windows,
                reduced.num_modes,
                reduced.num_windows * reduced.num_modes
                    * monolis_mpi_get_global_comm_size());
        }

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_reduced_solver_solve(&reduced_solver, &reduced)
               != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("all-at-once reduced solve failed");
            }
            timing.all_at_once_solve =
                monolis_get_time_global_sync() - t0;
        }

        all_at_once_residual =
            ROM_std_stdd_reduced_algebraic_residual(&reduced_solver);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "all-at-once algebraic residual = %.15e\n",
                all_at_once_residual);
            fflush(stdout);
        }
    }

    if(solver_strategy == DIFF_STDD_SOLVER_COMPARE){
        const size_t coefficient_count =
            (size_t)reduced.num_windows * (size_t)reduced.num_modes_capacity;

        all_at_once_coefficients = (double*)calloc(
            coefficient_count,
            sizeof(double));
        if(all_at_once_coefficients == NULL){
            diff_stdd_fail("all-at-once coefficient backup allocation failed");
        }
        memcpy(
            all_at_once_coefficients,
            reduced.reduced_coef,
            coefficient_count * sizeof(double));

        memset(
            reduced.reduced_coef,
            0,
            coefficient_count * sizeof(double));

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_solve_window_marching(
                &reduced,
                &sys.monolis_com,
                linear_max_iter,
                linear_tolerance) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("window-marching reduced solve failed");
            }
            timing.window_marching_solve =
                monolis_get_time_global_sync() - t0;
        }

        window_marching_residual =
            ROM_std_stdd_system_reduced_residual(
                &reduced,
                &sys.monolis_com);
        coefficient_difference =
            diff_stdd_coefficient_relative_difference(
                &reduced,
                all_at_once_coefficients,
                reduced.reduced_coef);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "window-marching reduced residual = %.15e\n",
                window_marching_residual);
            printf(
                "all-at-once/window-marching coefficient relative difference "
                "= %.15e\n",
                coefficient_difference);
            fflush(stdout);
        }

        /* In compare mode, use the all-at-once trajectory for reconstruction
         * and VTK/error output. */
        memcpy(
            reduced.reduced_coef,
            all_at_once_coefficients,
            coefficient_count * sizeof(double));
    }
    else if(solver_strategy == DIFF_STDD_SOLVER_WINDOW_MARCHING){
        const size_t coefficient_count =
            (size_t)reduced.num_windows * (size_t)reduced.num_modes_capacity;

        memset(
            reduced.reduced_coef,
            0,
            coefficient_count * sizeof(double));

        {
            const double t0 = monolis_get_time_global_sync();
            if(ROM_std_stdd_solve_window_marching(
                &reduced,
                &sys.monolis_com,
                linear_max_iter,
                linear_tolerance) != ROM_STDD_SUCCESS)
            {
                diff_stdd_fail("window-marching reduced solve failed");
            }
            timing.window_marching_solve =
                monolis_get_time_global_sync() - t0;
        }

        window_marching_residual =
            ROM_std_stdd_system_reduced_residual(
                &reduced,
                &sys.monolis_com);
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "window-marching reduced residual = %.15e\n",
                window_marching_residual);
            fflush(stdout);
        }
    }

    {
        char rank_coefficients[4096];
        const int rank = monolis_mpi_get_global_my_rank();
        int written = snprintf(
            rank_coefficients,
            sizeof(rank_coefficients),
            "%s.rank_%06d",
            coefficient_filename,
            rank);
        if(written < 0 || (size_t)written >= sizeof(rank_coefficients) ||
           diff_stdd_write_coefficients(&reduced, rank_coefficients) != 0)
        {
            diff_stdd_fail("cannot write ST-DDROM coefficients");
        }
    }

    {
        const double t0 = monolis_get_time_global_sync();
        if(diff_stdd_output_solution(
            &adapter,
            sys.vals.output_interval,
            write_output,
            error_filename,
            &timing.reconstruction) != 0)
        {
            diff_stdd_fail("ST-DDROM reconstruction/output failed");
        }
        timing.reconstruction_output_stage =
            monolis_get_time_global_sync() - t0;
    }

    timing.total_stddrom = monolis_get_time_global_sync() - total_t0;

    diff_stdd_print_timing_summary(
        solver_strategy,
        &timing,
        coefficient_difference,
        all_at_once_residual,
        window_marching_residual);

    if(monolis_mpi_get_global_my_rank() == 0){
        if(diff_stdd_write_timing_csv(
            timing_filename,
            mode,
            solver_strategy,
            monolis_mpi_get_global_comm_size(),
            total_steps,
            slabs_per_window,
            num_windows,
            reduced.num_modes,
            &timing,
            coefficient_difference,
            all_at_once_residual,
            window_marching_residual) != 0)
        {
            fprintf(stderr, "WARNING: cannot write %s\n", timing_filename);
        }

        printf(
            "diffusion ST-DDROM completed: steps=%d final_time=%.8e "
            "solver=%s\n",
            total_steps,
            sys.vals.finish_time,
            diff_stdd_solver_strategy_name(solver_strategy));
        fflush(stdout);
    }

    status = EXIT_SUCCESS;

cleanup:
    free(snapshot_matrix);
    free(all_at_once_coefficients);
    if(reduced_solver_initialized){
        ROM_std_stdd_reduced_solver_finalize(&reduced_solver);
    }
    ROM_std_stdd_diffusion_context_finalize(&adapter);
    ROM_std_stdd_system_finalize(&reduced);
    diff_stdd_finalize_fom(&sys);
    monolis_global_finalize();

    return status;
}
