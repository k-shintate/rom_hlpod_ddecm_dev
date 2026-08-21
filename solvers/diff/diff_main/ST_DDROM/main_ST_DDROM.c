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
    const char* error_filename)
{
    FE_SYSTEM* sys = adapter->fom;
    ROM_STDD_SYSTEM* reduced = adapter->reduced;
    const int ns = adapter->layout->num_slabs;
    double* window_state = NULL;
    FILE* fp = NULL;
    int file_number = 0;
    int step = 0;

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
        if(ROM_std_stdd_diffusion_reconstruct_window(
            adapter,
            window,
            window_state) != ROM_STDD_SUCCESS)
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
    free(window_state);
    return 0;
}

int main(int argc, char* argv[])
{
    FE_SYSTEM sys;
    ST_FOM_WINDOW_LAYOUT layout;
    ROM_STDD_SYSTEM reduced;
    ROM_STDD_DIFFUSION_CONTEXT adapter;
    ROM_STDD_REDUCED_SOLVER reduced_solver;
    DIFF_STDD_MODE mode;

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

    double* snapshot_matrix = NULL;
    char basis_directory[4096];
    char operator_directory[4096];
    char coefficient_filename[4096];
    char error_filename[4096];
    int status = EXIT_FAILURE;

    memset(&sys, 0, sizeof(sys));
    memset(&layout, 0, sizeof(layout));
    memset(&reduced, 0, sizeof(reduced));
    memset(&adapter, 0, sizeof(adapter));
    memset(&reduced_solver, 0, sizeof(reduced_solver));

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

    num_training_windows = diff_stdd_env_int(
        "DIFF_ST_DDROM_TRAINING_WINDOWS",
        num_windows);
    if(num_training_windows <= 0 || num_training_windows > num_windows){
        diff_stdd_fail("DIFF_ST_DDROM_TRAINING_WINDOWS must be in [1,num_windows]");
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
        "stddrom_l2_error.csv") != 0)
    {
        diff_stdd_fail("ST-DDROM output path too long");
    }

    mode = diff_stdd_read_mode();

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "\nDiffusion all-at-once ST-DDROM\n"
            "  mode             : %s\n"
            "  time scheme      : dG(0)\n"
            "  total steps      : %d\n"
            "  slabs/window     : %d\n"
            "  windows          : %d\n"
            "  training windows : %d\n"
            "  local/owned nodes: %d/%d\n"
            "  max modes/rank   : %d\n"
            "  POD energy       : %.8e\n"
            "  basis directory  : %s\n"
            "  operator dir     : %s\n",
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
            operator_directory);
    }

    if(mode == DIFF_STDD_MODE_OFFLINE || mode == DIFF_STDD_MODE_FULL){
        snapshot_matrix = (double*)calloc(
            (size_t)reduced.internal_st_dof * (size_t)num_training_windows,
            sizeof(double));
        if(snapshot_matrix == NULL){
            diff_stdd_fail("snapshot matrix allocation failed");
        }

        if(diff_stdd_collect_training_snapshots(
            &sys,
            &layout,
            num_training_windows,
            snapshot_matrix) != 0)
        {
            diff_stdd_fail("FOM training snapshot collection failed");
        }

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

        if(ROM_std_stdd_project_operators_rank_sweep(
            &reduced,
            ROM_std_stdd_diffusion_apply_A,
            ROM_std_stdd_diffusion_apply_B,
            &adapter) != ROM_STDD_SUCCESS)
        {
            diff_stdd_fail("diffusion reduced-operator projection failed");
        }

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

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "offline projection complete: modes/rank=%d\n",
                reduced.num_modes);
        }
    }

    free(snapshot_matrix);
    snapshot_matrix = NULL;

    if(mode == DIFF_STDD_MODE_OFFLINE){
        status = EXIT_SUCCESS;
        goto cleanup;
    }

    if(mode == DIFF_STDD_MODE_ONLINE){
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
    }

    /* The all-at-once RHS is assembled for every window before any solve. */
    if(ROM_std_stdd_diffusion_project_all_rhs(&adapter) != ROM_STDD_SUCCESS){
        diff_stdd_fail("all-window reduced RHS projection failed");
    }

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

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "solving all windows simultaneously: windows=%d modes/rank=%d "
            "global reduced dof=%d\n",
            reduced.num_windows,
            reduced.num_modes,
            reduced.num_windows * reduced.num_modes
                * monolis_mpi_get_global_comm_size());
    }

    if(ROM_std_stdd_reduced_solver_solve(&reduced_solver, &reduced)
       != ROM_STDD_SUCCESS)
    {
        diff_stdd_fail("all-at-once reduced solve failed");
    }

    {
        /*
         * ROM_std_stdd_reduced_algebraic_residual() performs halo exchange
         * and MPI all-reductions internally, so every MPI rank must enter
         * this routine.  Only the printing itself is restricted to rank 0.
         */
        const double algebraic_residual =
            ROM_std_stdd_reduced_algebraic_residual(&reduced_solver);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "all-at-once algebraic residual = %.15e\n",
                algebraic_residual);
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

    if(diff_stdd_output_solution(
        &adapter,
        sys.vals.output_interval,
        write_output,
        error_filename) != 0)
    {
        diff_stdd_fail("ST-DDROM reconstruction/output failed");
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "diffusion ST-DDROM completed: steps=%d final_time=%.8e\n",
            total_steps,
            sys.vals.finish_time);
        fflush(stdout);
    }

    status = EXIT_SUCCESS;

cleanup:
    free(snapshot_matrix);
    if(reduced_solver_initialized){
        ROM_std_stdd_reduced_solver_finalize(&reduced_solver);
    }
    ROM_std_stdd_diffusion_context_finalize(&adapter);
    ROM_std_stdd_system_finalize(&reduced);
    diff_stdd_finalize_fom(&sys);
    monolis_global_finalize();

    return status;
}