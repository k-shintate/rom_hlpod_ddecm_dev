#include "core_ROM.h"
#include "core_FOM_ST_spaceblock.h"
#include "rom_std_stddrom.h"
#include "rom_std_stddrom_stsb.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef ST_DDROM_DEGREE
#define ST_DDROM_DEGREE 1
#endif

#ifndef ST_DDROM_WINDOW_SLABS
#define ST_DDROM_WINDOW_SLABS 4
#endif

typedef enum {
    ST_DDROM_RUN_FOM = 0,
    ST_DDROM_RUN_COLLECT,
    ST_DDROM_RUN_OFFLINE,
    ST_DDROM_RUN_ONLINE
} ST_DDROM_RUN_MODE;

typedef enum {
    ST_DDROM_REDUCED_WINDOW_MARCHING = 0,
    ST_DDROM_REDUCED_ALL_AT_ONCE
} ST_DDROM_REDUCED_SOLVER_MODE;

typedef struct {
    ST_DDROM_RUN_MODE mode;
    int snapshot_id;
    int num_snapshots;
    int max_modes;
    int validate;
    int validation_snapshot_id;
    int write_output;
    int check_fom_residual;
    int strict_fom_residual;
    ST_DDROM_REDUCED_SOLVER_MODE reduced_solver_mode;
    int reduced_max_iter;
    double epsilon;
    double reduced_tolerance;
    double fom_residual_tolerance;
    char snapshot_dir[4096];
    char basis_dir[4096];
} ST_DDROM_OPTIONS;

typedef struct {
    double projection_difference2;
    double reference2;
    double state_difference2;
    double coefficient_difference2;
    double coefficient_reference2;
    double trace_difference2;
    double trace_reference2;
} ST_DDROM_VALIDATION_SUM;

typedef STSB_CONSTRAINED_RESIDUAL_SUM ST_DDROM_FOM_SOLVE_RESIDUAL_SUM;

typedef struct {
    double fom_residual_all2;
    double fom_rhs_all2;
    double fom_residual_free2;
    double fom_rhs_free2;
    double fom_residual_dirichlet2;
    double fom_rhs_dirichlet2;
    double homogeneous_dirichlet2;
    double homogeneous_all2;

    double direct_projected_defect2;
    double assembled_defect2;
    double identity_difference2;
    double reduced_rhs2;

    double self_action2;
    double neighbor_action2;
    double temporal_action2;
} ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM;

static void fail_message(const char* message)
{
    const int rank = monolis_mpi_get_global_my_rank();

    fprintf(
        stderr,
        "ERROR [rank %d]: %s\n",
        rank,
        (message != NULL) ? message : "unknown failure");
    fflush(stderr);

    /*
     * Do not call a low-level communicator abort here.  A nonzero process status is
     * sufficient for mpirun to terminate the job, and keeps this source on
     * the MONOLIS communication interface only.
     */
    exit(EXIT_FAILURE);
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
    if(errno != 0 || end == text || *end != '\0' ||
       value < 0 || value > INT_MAX)
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
    if(errno != 0 || end == text || *end != '\0' || !isfinite(value)){
        fprintf(stderr, "ERROR: invalid real %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return value;
}

static ST_DDROM_RUN_MODE read_mode(void)
{
    const char* text = getenv("ST_DDROM_MODE");

    if(text == NULL || strcmp(text, "fom") == 0){
        return ST_DDROM_RUN_FOM;
    }
    if(strcmp(text, "collect") == 0){
        return ST_DDROM_RUN_COLLECT;
    }
    if(strcmp(text, "offline") == 0){
        return ST_DDROM_RUN_OFFLINE;
    }
    if(strcmp(text, "online") == 0){
        return ST_DDROM_RUN_ONLINE;
    }

    fprintf(stderr, "ERROR: unknown ST_DDROM_MODE=%s\n", text);
    exit(EXIT_FAILURE);
}

static const char* mode_name(ST_DDROM_RUN_MODE mode)
{
    switch(mode){
    case ST_DDROM_RUN_FOM:
        return "fom";
    case ST_DDROM_RUN_COLLECT:
        return "collect";
    case ST_DDROM_RUN_OFFLINE:
        return "offline";
    case ST_DDROM_RUN_ONLINE:
        return "online";
    }
    return "unknown";
}

static ST_DDROM_REDUCED_SOLVER_MODE read_reduced_solver_mode(void)
{
    const char* text = getenv("ST_DDROM_REDUCED_SOLVER");

    if(text == NULL || text[0] == '\0' ||
       strcmp(text, "window") == 0 ||
       strcmp(text, "window_marching") == 0)
    {
        return ST_DDROM_REDUCED_WINDOW_MARCHING;
    }

    if(strcmp(text, "all_at_once") == 0 ||
       strcmp(text, "monolis_all_at_once") == 0)
    {
        return ST_DDROM_REDUCED_ALL_AT_ONCE;
    }

    fprintf(
        stderr,
        "ERROR: unknown ST_DDROM_REDUCED_SOLVER=%s\n"
        "       supported values: window, all_at_once\n",
        text);
    exit(EXIT_FAILURE);
}

static const char* reduced_solver_mode_name(
    ST_DDROM_REDUCED_SOLVER_MODE mode)
{
    switch(mode){
    case ST_DDROM_REDUCED_WINDOW_MARCHING:
        return "window_marching_block_scaled";
    case ST_DDROM_REDUCED_ALL_AT_ONCE:
        return "all_at_once_bicgstab_diag";
    }
    return "unknown";
}

static void join_path(
    char* output,
    size_t output_size,
    const char* first,
    const char* second)
{
    if(snprintf(output, output_size, "%s/%s", first, second)
       >= (int)output_size)
    {
        fail_message("path is too long");
    }
}

static void read_options(
    ST_DDROM_OPTIONS* options,
    const char* case_directory)
{
    const char* snapshot_override;
    const char* basis_override;

    if(options == NULL || case_directory == NULL){
        fail_message("NULL option argument");
    }

    memset(options, 0, sizeof(*options));
    options->mode = read_mode();
    options->snapshot_id = env_int("ST_DDROM_SNAPSHOT_ID", 0);
    options->num_snapshots = env_int("ST_DDROM_NUM_SNAPSHOTS", 1);
    options->max_modes = env_int("ST_DDROM_MAX_MODES", 10);
    options->validate = env_int("ST_DDROM_VALIDATE", 1);
    options->validation_snapshot_id =
        env_int("ST_DDROM_VALIDATION_SNAPSHOT_ID", 0);
    options->write_output = env_int("ST_DDROM_WRITE_OUTPUT", 1);
    options->check_fom_residual =
        env_int("ST_DDROM_CHECK_FOM_RESIDUAL", 1);
    /*
     * The adapter now assembles the same strong-Dirichlet system as the
     * reference FOM.  Reject inconsistent snapshots by default so an
     * operator/RHS drift cannot silently contaminate the POD basis again.
     */
    options->strict_fom_residual =
        env_int("ST_DDROM_STRICT_FOM_RESIDUAL", 1);
    options->reduced_solver_mode = read_reduced_solver_mode();
    options->reduced_max_iter =
        env_int("ST_DDROM_REDUCED_MAX_ITER", 10000);
    options->epsilon = env_double("ST_DDROM_EPSILON", 1.0);
    options->reduced_tolerance =
        env_double("ST_DDROM_REDUCED_EPSILON", 1.0e-10);
    options->fom_residual_tolerance =
        env_double("ST_DDROM_FOM_RESIDUAL_TOL", 1.0e-6);

    if(options->num_snapshots <= 0 || options->max_modes <= 0 ||
       options->reduced_max_iter <= 0 || options->epsilon <= 0.0 ||
       options->epsilon > 1.0 || options->reduced_tolerance <= 0.0 ||
       options->fom_residual_tolerance <= 0.0 ||
       options->check_fom_residual > 1 ||
       options->strict_fom_residual > 1)
    {
        fail_message("invalid ST-DDROM runtime option");
    }

    snapshot_override = getenv("ST_DDROM_SNAPSHOT_DIR");
    basis_override = getenv("ST_DDROM_BASIS_DIR");

    if(snapshot_override != NULL && snapshot_override[0] != '\0'){
        snprintf(
            options->snapshot_dir,
            sizeof(options->snapshot_dir),
            "%s",
            snapshot_override);
    }
    else{
        join_path(
            options->snapshot_dir,
            sizeof(options->snapshot_dir),
            case_directory,
            "stddrom_snapshots");
    }

    if(basis_override != NULL && basis_override[0] != '\0'){
        snprintf(
            options->basis_dir,
            sizeof(options->basis_dir),
            "%s",
            basis_override);
    }
    else{
        join_path(
            options->basis_dir,
            sizeof(options->basis_dir),
            case_directory,
            "stddrom_basis");
    }
}

static void make_directory_collective(const char* directory)
{
    int local_status = 0;
    int global_status;

    if(directory == NULL || directory[0] == '\0'){
        local_status = 1;
    }
    else if(mkdir(directory, 0775) != 0 && errno != EEXIST){
        fprintf(
            stderr,
            "ERROR [rank %d]: mkdir(%s): %s\n",
            monolis_mpi_get_global_my_rank(),
            directory,
            strerror(errno));
        fflush(stderr);
        local_status = 1;
    }

    /*
     * Every rank executes mkdir.  Concurrent mkdir is safe here because
     * EEXIST is accepted.  The allreduce both propagates an error from any
     * rank and synchronizes all ranks before the directory is used.
     */
    global_status = local_status;
    monolis_allreduce_I(
        1,
        &global_status,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    if(global_status != 0){
        fail_message("directory creation failed on at least one rank");
    }
}

static void initialize_output(FE_SYSTEM* sys, int rank)
{
    FILE* fp = NULL;

    fp = BBFE_sys_write_fopen(
        fp,
        "l2_error.txt",
        sys->cond.directory);
    fclose(fp);

    if(rank == 0){
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

static double safe_relative_norm(double numerator2, double denominator2)
{
    if(numerator2 < 0.0 || denominator2 < 0.0 ||
       !isfinite(numerator2) || !isfinite(denominator2))
    {
        return INFINITY;
    }

    return sqrt(numerator2) /
        fmax(sqrt(denominator2), 1.0e-300);
}

static void output_window(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int window_id,
    const double* full_window_state,
    int* file_number)
{
    for(int local_slab = 0;
        local_slab < num_window_slabs;
        local_slab++)
    {
        const int global_step =
            window_id * num_window_slabs + local_slab + 1;
        const double time = (double)global_step * sys->vals.dt;

        if(global_step % sys->vals.output_interval != 0){
            continue;
        }

        STSB_extract_slab_right_trace(
            time_element,
            num_window_slabs,
            local_slab,
            sys->fe.total_num_nodes,
            full_window_state,
            sys->vals.T);

        output_files(sys, *file_number, time);
        (*file_number)++;
    }
}

static void run_fom_or_collect(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int num_windows,
    int num_internal_space_nodes,
    const ST_DDROM_OPTIONS* options,
    ROM_STDD_STSB_CONTEXT* adapter,
    ROM_STDD_SYSTEM* rom_system)
{
    const int ndof = adapter->st_dof_per_space_node;
    const int local_st_dof = adapter->local_st_dof;
    const int internal_st_dof = num_internal_space_nodes * ndof;
    double* window_solution = NULL;
    double* previous_trace = NULL;
    double* lift = NULL;
    double* snapshots = NULL;
    ST_DDROM_FOM_SOLVE_RESIDUAL_SUM cumulative_residual;
    int file_number = 0;
    int residual_violation = 0;

    memset(&cumulative_residual, 0, sizeof(cumulative_residual));

    window_solution = (double*)calloc((size_t)local_st_dof, sizeof(double));
    previous_trace = (double*)calloc(
        (size_t)sys->fe.total_num_nodes,
        sizeof(double));
    lift = (double*)calloc((size_t)local_st_dof, sizeof(double));

    if(options->mode == ST_DDROM_RUN_COLLECT){
        snapshots = (double*)calloc(
            (size_t)num_windows * (size_t)internal_st_dof,
            sizeof(double));
    }

    if(window_solution == NULL || previous_trace == NULL || lift == NULL ||
       (options->mode == ST_DDROM_RUN_COLLECT && snapshots == NULL))
    {
        fail_message("FOM/collect allocation failed");
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys->fe,
        sys->vals.T,
        previous_trace);

    for(int window = 0; window < num_windows; window++){
        const double window_start_time =
            (double)(window * num_window_slabs) * sys->vals.dt;
        int window_residual_violation = 0;
        double window_residual_free = 0.0;

        memset(window_solution, 0, (size_t)local_st_dof * sizeof(double));

        /* Exact verified reference-FOM solve path. */
        STSB_solve_window(
            sys,
            time_element,
            num_window_slabs,
            window_start_time,
            window + 1,
            previous_trace,
            window_solution);

        /*
         * Evaluate the exact matrix/RHS pair left by STSB_solve_window().
         * This is a read-only check of the system that was actually solved.
         */
        if(options->check_fom_residual){
            ST_DDROM_FOM_SOLVE_RESIDUAL_SUM window_residual;
            double residual_all;
            double residual_bc_all_rhs;

            memset(&window_residual, 0, sizeof(window_residual));

            if(STSB_evaluate_constrained_residual(
                sys,
                ndof,
                window_solution,
                &window_residual) != 0)
            {
                fail_message(
                    "failed to evaluate the solved strong-BC FOM residual");
            }

            cumulative_residual.residual_all2 +=
                window_residual.residual_all2;
            cumulative_residual.rhs_all2 +=
                window_residual.rhs_all2;
            cumulative_residual.residual_free2 +=
                window_residual.residual_free2;
            cumulative_residual.rhs_free2 +=
                window_residual.rhs_free2;
            cumulative_residual.residual_dirichlet2 +=
                window_residual.residual_dirichlet2;
            cumulative_residual.rhs_dirichlet2 +=
                window_residual.rhs_dirichlet2;

            residual_all = safe_relative_norm(
                window_residual.residual_all2,
                window_residual.rhs_all2);
            window_residual_free = safe_relative_norm(
                window_residual.residual_free2,
                window_residual.rhs_free2);
            residual_bc_all_rhs = sqrt(
                fmax(window_residual.residual_dirichlet2, 0.0)) /
                fmax(
                    sqrt(fmax(window_residual.rhs_all2, 0.0)),
                    1.0e-300);

            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "STSB solved-system residual window %d: "
                    "all=%.15e free=%.15e BC/all-rhs=%.15e\n",
                    window,
                    residual_all,
                    window_residual_free,
                    residual_bc_all_rhs);
            }

            if(!isfinite(window_residual_free) ||
               window_residual_free > options->fom_residual_tolerance)
            {
                window_residual_violation = 1;
                residual_violation = 1;
            }
        }

        /* Keep the same output ordering as the verified reference FOM. */
        if(options->write_output){
            output_window(
                sys,
                time_element,
                num_window_slabs,
                window,
                window_solution,
                &file_number);
        }

        /*
         * Snapshot extraction is read-only and is the only collect-specific
         * operation inserted into the reference FOM loop.
         */
        if(options->mode == ST_DDROM_RUN_COLLECT){
            if(ROM_std_stdd_stsb_build_dirichlet_lift(
                adapter,
                window,
                lift) != ROM_STDD_SUCCESS ||
               ROM_std_stdd_stsb_store_internal_homogeneous(
                rom_system,
                window_solution,
                lift,
                snapshots
                    + (size_t)window * (size_t)internal_st_dof)
                   != ROM_STDD_SUCCESS)
            {
                fail_message("snapshot lifting failed");
            }
        }

        /* Exact reference-FOM trace propagation. */
        STSB_extract_last_right_trace(
            time_element,
            num_window_slabs,
            sys->fe.total_num_nodes,
            window_solution,
            previous_trace);

        {
            double solution_l2 = 0.0;
            double trace_l2 = 0.0;

            if(STSB_compute_window_checksum(
                sys,
                ndof,
                window_solution,
                previous_trace,
                &solution_l2,
                &trace_l2) != 0)
            {
                fail_message("failed to compute FOM reference checksum");
            }

            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "STSB reference checksum window %d: "
                    "solution=%.15e trace=%.15e\n",
                    window,
                    solution_l2,
                    trace_l2);
            }
        }

        /* Strict rejection is opt-in and happens only after reference checks. */
        if(options->strict_fom_residual && window_residual_violation){
            if(monolis_mpi_get_global_my_rank() == 0){
                fprintf(
                    stderr,
                    "ERROR: FOM-consistent residual exceeds the optional "
                    "strict tolerance.\n"
                    "  window        = %d\n"
                    "  free residual = %.15e\n"
                    "  tolerance     = %.15e\n"
                    "The snapshot file has not been written.\n",
                    window,
                    window_residual_free,
                    options->fom_residual_tolerance);
                fflush(stderr);
            }
            fail_message("strict FOM-consistent residual check failed");
        }
    }

    if(options->check_fom_residual){
        const double global_all = safe_relative_norm(
            cumulative_residual.residual_all2,
            cumulative_residual.rhs_all2);
        const double global_free = safe_relative_norm(
            cumulative_residual.residual_free2,
            cumulative_residual.rhs_free2);
        const double global_bc_all_rhs = sqrt(
            fmax(cumulative_residual.residual_dirichlet2, 0.0)) /
            fmax(
                sqrt(fmax(cumulative_residual.rhs_all2, 0.0)),
                1.0e-300);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\nST-DDROM FOM-consistent residual summary\n"
                "  strong-BC residual all          : %.15e\n"
                "  strong-BC residual free         : %.15e\n"
                "  strong-BC residual BC/all-rhs   : %.15e\n"
                "  comparison tolerance            : %.15e\n"
                "  diagnostic status               : %s\n"
                "  strict rejection enabled        : %s\n",
                global_all,
                global_free,
                global_bc_all_rhs,
                options->fom_residual_tolerance,
                residual_violation ? "ABOVE TOLERANCE" : "WITHIN TOLERANCE",
                options->strict_fom_residual ? "yes" : "no");
        }
    }

    if(options->mode == ST_DDROM_RUN_COLLECT){
        make_directory_collective(options->snapshot_dir);
        if(ROM_std_stdd_write_snapshot_rank(
            rom_system,
            snapshots,
            options->snapshot_id,
            options->snapshot_dir) != ROM_STDD_SUCCESS)
        {
            fail_message("snapshot write failed");
        }
    }

    free(window_solution);
    free(previous_trace);
    free(lift);
    free(snapshots);
}

static void run_offline(
    ROM_STDD_SYSTEM* system,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ST_DDROM_OPTIONS* options)
{
    const int num_columns =
        options->num_snapshots * system->num_windows;
    double* one_case = NULL;
    double* snapshot_matrix = NULL;

    one_case = (double*)calloc(
        (size_t)system->num_windows
            * (size_t)system->internal_st_dof,
        sizeof(double));
    snapshot_matrix = (double*)calloc(
        (size_t)system->internal_st_dof * (size_t)num_columns,
        sizeof(double));

    if(one_case == NULL || snapshot_matrix == NULL){
        fail_message("offline snapshot allocation failed");
    }

    for(int snapshot = 0; snapshot < options->num_snapshots; snapshot++){
        if(ROM_std_stdd_read_snapshot_rank(
            system,
            one_case,
            snapshot,
            options->snapshot_dir) != ROM_STDD_SUCCESS)
        {
            fail_message("offline snapshot read failed");
        }

        for(int window = 0; window < system->num_windows; window++){
            const int column = snapshot * system->num_windows + window;
            for(int row = 0; row < system->internal_st_dof; row++){
                snapshot_matrix[
                    (size_t)row * (size_t)num_columns
                    + (size_t)column] =
                    one_case[
                        (size_t)window * (size_t)system->internal_st_dof
                        + (size_t)row];
            }
        }
    }

    if(ROM_std_stdd_compute_shared_local_pod(
        system,
        snapshot_matrix,
        num_columns,
        options->epsilon,
        monolis_mpi_get_global_comm()) != ROM_STDD_SUCCESS)
    {
        fail_message("local shared POD failed");
    }

    if(ROM_std_stdd_stsb_enforce_homogeneous_basis(
        adapter,
        system) != ROM_STDD_SUCCESS)
    {
        fail_message(
            "failed to enforce homogeneous Dirichlet rows in the POD basis");
    }

    make_directory_collective(options->basis_dir);

    if(ROM_std_stdd_write_basis_rank(
        system,
        num_columns,
        options->basis_dir) != ROM_STDD_SUCCESS)
    {
        fail_message("basis write failed");
    }

    if(ROM_std_stdd_stsb_prepare_operator(adapter, 0)
       != ROM_STDD_SUCCESS)
    {
        fail_message("A preparation failed");
    }

    if(ROM_std_stdd_project_operators_rank_sweep(
        system,
        ROM_std_stdd_stsb_apply_prepared_A,
        ROM_std_stdd_stsb_apply_B,
        adapter) != ROM_STDD_SUCCESS)
    {
        fail_message("distributed reduced-operator projection failed");
    }

    if(ROM_std_stdd_write_operators_rank(
        system,
        options->basis_dir) != ROM_STDD_SUCCESS)
    {
        fail_message("operator write failed");
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-DDROM offline completed: windows=%d, modes/rank=%d, "
            "training columns=%d\n",
            system->num_windows,
            system->num_modes,
            num_columns);
    }

    free(one_case);
    free(snapshot_matrix);
}

static void accumulate_validation_window(
    ST_DDROM_VALIDATION_SUM* sum,
    const ROM_STDD_SYSTEM* system,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int window_id,
    const double* q_rom,
    const double* q_reference,
    double* coefficient_reference,
    double* q_projection,
    double* trace_rom,
    double* trace_reference)
{
    const double* coefficient_rom =
        system->reduced_coef
        + (size_t)window_id * (size_t)system->num_modes_capacity;

    ROM_std_stdd_project_reference_coefficients(
        system,
        window_id,
        q_reference,
        coefficient_reference);

    memset(
        q_projection,
        0,
        (size_t)system->internal_st_dof * sizeof(double));

    for(int row = 0; row < system->internal_st_dof; row++){
        for(int mode = 0; mode < system->num_modes; mode++){
            q_projection[row] += system->basis[
                (size_t)row * (size_t)system->num_modes_capacity
                + (size_t)mode]
                * coefficient_reference[mode];
        }

        {
            const double projection_difference =
                q_projection[row] - q_reference[row];
            const double state_difference = q_rom[row] - q_reference[row];
            sum->projection_difference2 +=
                projection_difference * projection_difference;
            sum->state_difference2 += state_difference * state_difference;
            sum->reference2 += q_reference[row] * q_reference[row];
        }
    }

    for(int mode = 0; mode < system->num_modes; mode++){
        const double difference =
            coefficient_rom[mode] - coefficient_reference[mode];
        sum->coefficient_difference2 += difference * difference;
        sum->coefficient_reference2 +=
            coefficient_reference[mode] * coefficient_reference[mode];
    }

    STSB_extract_last_right_trace(
        time_element,
        num_window_slabs,
        system->num_internal_space_nodes,
        q_rom,
        trace_rom);
    STSB_extract_last_right_trace(
        time_element,
        num_window_slabs,
        system->num_internal_space_nodes,
        q_reference,
        trace_reference);

    for(int node = 0; node < system->num_internal_space_nodes; node++){
        const double difference = trace_rom[node] - trace_reference[node];
        sum->trace_difference2 += difference * difference;
        sum->trace_reference2 += trace_reference[node] * trace_reference[node];
    }
}

static void add_validation_sum(
    ST_DDROM_VALIDATION_SUM* destination,
    const ST_DDROM_VALIDATION_SUM* source)
{
    destination->projection_difference2 += source->projection_difference2;
    destination->reference2 += source->reference2;
    destination->state_difference2 += source->state_difference2;
    destination->coefficient_difference2 += source->coefficient_difference2;
    destination->coefficient_reference2 += source->coefficient_reference2;
    destination->trace_difference2 += source->trace_difference2;
    destination->trace_reference2 += source->trace_reference2;
}

static void reduce_validation_sum(
    const ST_DDROM_VALIDATION_SUM* local_sum,
    ST_DDROM_VALIDATION_SUM* global_sum)
{
    double reduced_values[7] = {
        local_sum->projection_difference2,
        local_sum->reference2,
        local_sum->state_difference2,
        local_sum->coefficient_difference2,
        local_sum->coefficient_reference2,
        local_sum->trace_difference2,
        local_sum->trace_reference2
    };

    monolis_allreduce_R(
        7,
        reduced_values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    global_sum->projection_difference2 = reduced_values[0];
    global_sum->reference2 = reduced_values[1];
    global_sum->state_difference2 = reduced_values[2];
    global_sum->coefficient_difference2 = reduced_values[3];
    global_sum->coefficient_reference2 = reduced_values[4];
    global_sum->trace_difference2 = reduced_values[5];
    global_sum->trace_reference2 = reduced_values[6];
}


static void add_residual_diagnostic_sum(
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* destination,
    const ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* source)
{
    destination->fom_residual_all2 += source->fom_residual_all2;
    destination->fom_rhs_all2 += source->fom_rhs_all2;
    destination->fom_residual_free2 += source->fom_residual_free2;
    destination->fom_rhs_free2 += source->fom_rhs_free2;
    destination->fom_residual_dirichlet2 +=
        source->fom_residual_dirichlet2;
    destination->fom_rhs_dirichlet2 += source->fom_rhs_dirichlet2;
    destination->homogeneous_dirichlet2 +=
        source->homogeneous_dirichlet2;
    destination->homogeneous_all2 += source->homogeneous_all2;
    destination->direct_projected_defect2 +=
        source->direct_projected_defect2;
    destination->assembled_defect2 += source->assembled_defect2;
    destination->identity_difference2 += source->identity_difference2;
    destination->reduced_rhs2 += source->reduced_rhs2;
    destination->self_action2 += source->self_action2;
    destination->neighbor_action2 += source->neighbor_action2;
    destination->temporal_action2 += source->temporal_action2;
}

static void pack_residual_diagnostic_sum(
    const ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* source,
    double* output)
{
    output[0] = source->fom_residual_all2;
    output[1] = source->fom_rhs_all2;
    output[2] = source->fom_residual_free2;
    output[3] = source->fom_rhs_free2;
    output[4] = source->fom_residual_dirichlet2;
    output[5] = source->fom_rhs_dirichlet2;
    output[6] = source->homogeneous_dirichlet2;
    output[7] = source->homogeneous_all2;
    output[8] = source->direct_projected_defect2;
    output[9] = source->assembled_defect2;
    output[10] = source->identity_difference2;
    output[11] = source->reduced_rhs2;
    output[12] = source->self_action2;
    output[13] = source->neighbor_action2;
    output[14] = source->temporal_action2;
}

static void unpack_residual_diagnostic_sum(
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* destination,
    const double* input)
{
    destination->fom_residual_all2 = input[0];
    destination->fom_rhs_all2 = input[1];
    destination->fom_residual_free2 = input[2];
    destination->fom_rhs_free2 = input[3];
    destination->fom_residual_dirichlet2 = input[4];
    destination->fom_rhs_dirichlet2 = input[5];
    destination->homogeneous_dirichlet2 = input[6];
    destination->homogeneous_all2 = input[7];
    destination->direct_projected_defect2 = input[8];
    destination->assembled_defect2 = input[9];
    destination->identity_difference2 = input[10];
    destination->reduced_rhs2 = input[11];
    destination->self_action2 = input[12];
    destination->neighbor_action2 = input[13];
    destination->temporal_action2 = input[14];
}

static void reduce_residual_diagnostic_sum(
    const ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* local_sum,
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* global_sum)
{
    double values[15];

    pack_residual_diagnostic_sum(local_sum, values);

    monolis_allreduce_R(
        15,
        values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    unpack_residual_diagnostic_sum(global_sum, values);
}

static double diagnostic_relative(double numerator2, double denominator2)
{
    return sqrt(fmax(numerator2, 0.0))
        / fmax(sqrt(fmax(denominator2, 0.0)), 1.0e-300);
}

static void print_residual_diagnostics_by_window(
    const ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* local_window_sums,
    int num_windows)
{
    double* values;
    const int rank = monolis_mpi_get_global_my_rank();

    if(local_window_sums == NULL || num_windows <= 0){
        return;
    }

    values = (double*)calloc(
        (size_t)(15 * num_windows),
        sizeof(double));
    if(values == NULL){
        fail_message("residual diagnostics reduction allocation failed");
    }

    for(int window = 0; window < num_windows; window++){
        pack_residual_diagnostic_sum(
            &local_window_sums[window],
            values + 15 * window);
    }

    monolis_allreduce_R(
        15 * num_windows,
        values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    if(rank == 0){
        printf("\nST-DDROM FOM-consistent lifting/reduced identity diagnostics\n");
        printf(
            "  window       raw all      raw free  raw BC/allrhs"
            "       q BC/all     direct proj       assembled"
            "        identity\n");

        for(int window = 0; window < num_windows; window++){
            ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM value;
            const double* input = values + 15 * window;

            memset(&value, 0, sizeof(value));
            unpack_residual_diagnostic_sum(&value, input);

            printf(
                "  %6d  %12.5e  %12.5e  %12.5e  %12.5e"
                "  %14.7e  %14.7e  %14.7e\n",
                window,
                diagnostic_relative(
                    value.fom_residual_all2,
                    value.fom_rhs_all2),
                diagnostic_relative(
                    value.fom_residual_free2,
                    value.fom_rhs_free2),
                diagnostic_relative(
                    value.fom_residual_dirichlet2,
                    value.fom_rhs_all2),
                diagnostic_relative(
                    value.homogeneous_dirichlet2,
                    value.homogeneous_all2),
                diagnostic_relative(
                    value.direct_projected_defect2,
                    value.reduced_rhs2),
                diagnostic_relative(
                    value.assembled_defect2,
                    value.reduced_rhs2),
                diagnostic_relative(
                    value.identity_difference2,
                    value.reduced_rhs2));
        }

        {
            ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM window0;
            memset(&window0, 0, sizeof(window0));
            unpack_residual_diagnostic_sum(&window0, values);

            printf("\nwindow 0 reduced reference decomposition\n");
            printf(
                "  ||D_self a_self|| / ||g||       : %.15e\n",
                diagnostic_relative(
                    window0.self_action2,
                    window0.reduced_rhs2));
            printf(
                "  ||D_neighbor a_neighbor||/||g|| : %.15e\n",
                diagnostic_relative(
                    window0.neighbor_action2,
                    window0.reduced_rhs2));
            printf(
                "  ||g||                            : %.15e\n",
                sqrt(fmax(window0.reduced_rhs2, 0.0)));
            printf(
                "  ||D a - g|| / ||g||             : %.15e\n",
                diagnostic_relative(
                    window0.assembled_defect2,
                    window0.reduced_rhs2));
        }
    }

    free(values);
}

static void print_validation_window_table(
    const ST_DDROM_VALIDATION_SUM* local_window_sums,
    int num_windows,
    const double* rom_window_residuals,
    const double* reference_window_defects)
{
    double* reduced_values;
    const int rank = monolis_mpi_get_global_my_rank();

    if(local_window_sums == NULL || num_windows <= 0 ||
       rom_window_residuals == NULL || reference_window_defects == NULL)
    {
        return;
    }

    reduced_values = (double*)calloc(
        (size_t)(7 * num_windows),
        sizeof(double));
    if(reduced_values == NULL){
        fail_message("window validation reduction allocation failed");
    }

    for(int window = 0; window < num_windows; window++){
        const ST_DDROM_VALIDATION_SUM* value =
            &local_window_sums[window];
        double* output = reduced_values + 7 * window;

        output[0] = value->projection_difference2;
        output[1] = value->reference2;
        output[2] = value->state_difference2;
        output[3] = value->coefficient_difference2;
        output[4] = value->coefficient_reference2;
        output[5] = value->trace_difference2;
        output[6] = value->trace_reference2;
    }

    monolis_allreduce_R(
        7 * num_windows,
        reduced_values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    if(rank == 0){
        printf("\nST-DDROM validation by window\n");
        printf(
            "  window       projection             state       coefficient"
            "             trace      ROM defect      reference defect\n");

        for(int window = 0; window < num_windows; window++){
            const double* value = reduced_values + 7 * window;
            const double reference_norm =
                fmax(sqrt(value[1]), 1.0e-300);
            const double coefficient_norm =
                fmax(sqrt(value[4]), 1.0e-300);
            const double trace_norm =
                fmax(sqrt(value[6]), 1.0e-300);

            printf(
                "  %6d  %17.9e  %17.9e  %17.9e  %17.9e"
                "  %13.5e  %17.9e\n",
                window,
                sqrt(value[0]) / reference_norm,
                sqrt(value[2]) / reference_norm,
                sqrt(value[3]) / coefficient_norm,
                sqrt(value[5]) / trace_norm,
                rom_window_residuals[window],
                reference_window_defects[window]);
        }
    }

    free(reduced_values);
}

static void print_validation_summary(
    const ST_DDROM_VALIDATION_SUM* local_sum,
    const ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* local_diagnostic_sum,
    double reduced_residual,
    double reference_reduced_defect)
{
    ST_DDROM_VALIDATION_SUM global_sum;
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM global_diagnostic_sum;
    const int rank = monolis_mpi_get_global_my_rank();

    memset(&global_sum, 0, sizeof(global_sum));
    memset(&global_diagnostic_sum, 0, sizeof(global_diagnostic_sum));

    reduce_validation_sum(local_sum, &global_sum);
    reduce_residual_diagnostic_sum(
        local_diagnostic_sum,
        &global_diagnostic_sum);

    if(rank == 0){
        printf("\nST-DDROM validation summary\n");
        printf(
            "  global projection error          : %.15e\n",
            sqrt(global_sum.projection_difference2)
                / fmax(sqrt(global_sum.reference2), 1.0e-300));
        printf(
            "  global state error               : %.15e\n",
            sqrt(global_sum.state_difference2)
                / fmax(sqrt(global_sum.reference2), 1.0e-300));
        printf(
            "  global coefficient error         : %.15e\n",
            sqrt(global_sum.coefficient_difference2)
                / fmax(sqrt(global_sum.coefficient_reference2), 1.0e-300));
        printf(
            "  global trace error               : %.15e\n",
            sqrt(global_sum.trace_difference2)
                / fmax(sqrt(global_sum.trace_reference2), 1.0e-300));
        printf(
            "  ROM reduced residual             : %.15e\n",
            reduced_residual);
        printf(
            "  assembled reference defect       : %.15e\n",
            reference_reduced_defect);
        printf(
            "  direct projected FOM defect      : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.direct_projected_defect2,
                global_diagnostic_sum.reduced_rhs2));
        printf(
            "  assembled/direct identity error  : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.identity_difference2,
                global_diagnostic_sum.reduced_rhs2));
        printf(
            "  strong-BC lifted residual all    : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.fom_residual_all2,
                global_diagnostic_sum.fom_rhs_all2));
        printf(
            "  strong-BC lifted residual free   : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.fom_residual_free2,
                global_diagnostic_sum.fom_rhs_free2));
        printf(
            "  strong-BC lifted residual BC/rhs : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.fom_residual_dirichlet2,
                global_diagnostic_sum.fom_rhs_all2));
        printf(
            "  homogeneous BC constraint error  : %.15e\n",
            diagnostic_relative(
                global_diagnostic_sum.homogeneous_dirichlet2,
                global_diagnostic_sum.homogeneous_all2));
    }
}

static void run_online(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    ROM_STDD_SYSTEM* system,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ST_DDROM_OPTIONS* options)
{
    ROM_STDD_REDUCED_SOLVER reduced_solver;
    double* initial_trace = NULL;
    double* effective_rhs = NULL;
    double* homogeneous_work = NULL;
    double* full_window_state = NULL;
    double* reference_snapshots = NULL;
    double* coefficient_reference = NULL;
    double* reference_reduced_coef = NULL;
    double* rom_window_residuals = NULL;
    double* reference_window_defects = NULL;
    double* q_projection = NULL;
    double* trace_rom = NULL;
    double* trace_reference = NULL;
    double* direct_projected_defect = NULL;
    double* assembled_defect = NULL;
    double* self_action = NULL;
    double* neighbor_action = NULL;
    double* temporal_action = NULL;
    ST_DDROM_VALIDATION_SUM validation_sum;
    ST_DDROM_VALIDATION_SUM* window_validation_sums = NULL;
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM residual_diagnostic_sum;
    ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* window_residual_diagnostics = NULL;
    int file_number = 0;
    double reduced_residual = INFINITY;
    double reference_reduced_defect = INFINITY;

    memset(&reduced_solver, 0, sizeof(reduced_solver));
    memset(&validation_sum, 0, sizeof(validation_sum));
    memset(&residual_diagnostic_sum, 0, sizeof(residual_diagnostic_sum));

    if(ROM_std_stdd_read_basis_rank(system, options->basis_dir)
       != ROM_STDD_SUCCESS ||
       ROM_std_stdd_read_operators_rank(system, options->basis_dir)
       != ROM_STDD_SUCCESS)
    {
        fail_message("online basis/operator read failed");
    }

    initial_trace = (double*)calloc(
        (size_t)sys->fe.total_num_nodes,
        sizeof(double));
    effective_rhs = (double*)calloc(
        (size_t)system->local_st_dof,
        sizeof(double));
    homogeneous_work = (double*)calloc(
        (size_t)system->local_st_dof,
        sizeof(double));
    full_window_state = (double*)calloc(
        (size_t)system->local_st_dof,
        sizeof(double));

    if(initial_trace == NULL || effective_rhs == NULL ||
       homogeneous_work == NULL || full_window_state == NULL)
    {
        fail_message("online work allocation failed");
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys->fe,
        sys->vals.T,
        initial_trace);

    for(int window = 0; window < system->num_windows; window++){
        double* reduced_rhs = system->reduced_rhs
            + (size_t)window * (size_t)system->num_modes_capacity;

        if(ROM_std_stdd_stsb_build_effective_rhs(
            adapter,
            window,
            initial_trace,
            effective_rhs) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_project_local_rhs(
            system,
            window,
            effective_rhs,
            reduced_rhs) != ROM_STDD_SUCCESS)
        {
            fail_message("online reduced RHS construction failed");
        }
    }

    if(options->reduced_solver_mode ==
       ST_DDROM_REDUCED_WINDOW_MARCHING)
    {
        if(ROM_std_stdd_solve_window_marching(
            system,
            &sys->monolis_com,
            options->reduced_max_iter,
            options->reduced_tolerance) != ROM_STDD_SUCCESS)
        {
            fail_message("distributed window-marching reduced solve failed");
        }

    }
    else{
        if(ROM_std_stdd_reduced_solver_initialize(
            &reduced_solver,
            system,
            &sys->monolis_com,
            options->reduced_max_iter,
            options->reduced_tolerance) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_reduced_solver_assemble(
            &reduced_solver,
            system) != ROM_STDD_SUCCESS)
        {
            fail_message("reduced solver initialization/assembly failed");
        }

        make_directory_collective(options->basis_dir);
        ROM_std_stdd_write_reduced_graph_rank(
            &reduced_solver,
            options->basis_dir);

        if(ROM_std_stdd_reduced_solver_solve(
            &reduced_solver,
            system) != ROM_STDD_SUCCESS)
        {
            fail_message("distributed reduced solve failed");
        }

    }

    reduced_residual =
        ROM_std_stdd_system_reduced_residual_from_coefficients(
            system,
            &sys->monolis_com,
            system->reduced_coef,
            NULL);

    if(options->validate){
        reference_snapshots = (double*)calloc(
            (size_t)system->num_windows
                * (size_t)system->internal_st_dof,
            sizeof(double));
        coefficient_reference = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        reference_reduced_coef = (double*)calloc(
            (size_t)system->num_windows
                * (size_t)system->num_modes_capacity,
            sizeof(double));
        rom_window_residuals = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        reference_window_defects = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        window_validation_sums = (ST_DDROM_VALIDATION_SUM*)calloc(
            (size_t)system->num_windows,
            sizeof(ST_DDROM_VALIDATION_SUM));
        q_projection = (double*)calloc(
            (size_t)system->internal_st_dof,
            sizeof(double));
        trace_rom = (double*)calloc(
            (size_t)system->num_internal_space_nodes,
            sizeof(double));
        trace_reference = (double*)calloc(
            (size_t)system->num_internal_space_nodes,
            sizeof(double));
        direct_projected_defect = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        assembled_defect = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        self_action = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        neighbor_action = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        temporal_action = (double*)calloc(
            (size_t)system->num_modes,
            sizeof(double));
        window_residual_diagnostics =
            (ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM*)calloc(
                (size_t)system->num_windows,
                sizeof(ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM));

        if(reference_snapshots == NULL || coefficient_reference == NULL ||
           reference_reduced_coef == NULL ||
           rom_window_residuals == NULL ||
           reference_window_defects == NULL ||
           window_validation_sums == NULL ||
           q_projection == NULL || trace_rom == NULL ||
           trace_reference == NULL ||
           direct_projected_defect == NULL ||
           assembled_defect == NULL ||
           self_action == NULL ||
           neighbor_action == NULL ||
           temporal_action == NULL ||
           window_residual_diagnostics == NULL)
        {
            fail_message("validation allocation failed");
        }

        if(ROM_std_stdd_read_snapshot_rank(
            system,
            reference_snapshots,
            options->validation_snapshot_id,
            options->snapshot_dir) != ROM_STDD_SUCCESS)
        {
            fail_message("validation snapshot read failed");
        }
    }

    for(int window = 0; window < system->num_windows; window++){
        if(ROM_std_stdd_stsb_reconstruct_full_window(
            adapter,
            system,
            window,
            homogeneous_work,
            full_window_state) != ROM_STDD_SUCCESS)
        {
            fail_message("online reconstruction failed");
        }

        if(options->validate){
            const double* reference = reference_snapshots
                + (size_t)window * (size_t)system->internal_st_dof;

            accumulate_validation_window(
                &window_validation_sums[window],
                system,
                time_element,
                num_window_slabs,
                window,
                homogeneous_work,
                reference,
                coefficient_reference,
                q_projection,
                trace_rom,
                trace_reference);

            add_validation_sum(
                &validation_sum,
                &window_validation_sums[window]);

            memcpy(
                reference_reduced_coef
                    + (size_t)window
                        * (size_t)system->num_modes_capacity,
                coefficient_reference,
                (size_t)system->num_modes * sizeof(double));
        }

        if(options->write_output){
            output_window(
                sys,
                time_element,
                num_window_slabs,
                window,
                full_window_state,
                &file_number);
        }
    }

    if(options->validate){
        for(int window = 0; window < system->num_windows; window++){
            ROM_STDD_STSB_RESIDUAL_SUM fom_residual_sum;
            ST_DDROM_RESIDUAL_DIAGNOSTIC_SUM* diagnostic =
                &window_residual_diagnostics[window];
            const double* current_reference =
                reference_snapshots
                + (size_t)window * (size_t)system->internal_st_dof;
            const double* previous_reference =
                (window > 0)
                ? reference_snapshots
                    + (size_t)(window - 1)
                        * (size_t)system->internal_st_dof
                : NULL;
            const double* reduced_rhs =
                system->reduced_rhs
                + (size_t)window
                    * (size_t)system->num_modes_capacity;

            memset(&fom_residual_sum, 0, sizeof(fom_residual_sum));
            memset(
                direct_projected_defect,
                0,
                (size_t)system->num_modes * sizeof(double));
            memset(
                assembled_defect,
                0,
                (size_t)system->num_modes * sizeof(double));
            memset(
                self_action,
                0,
                (size_t)system->num_modes * sizeof(double));
            memset(
                neighbor_action,
                0,
                (size_t)system->num_modes * sizeof(double));
            memset(
                temporal_action,
                0,
                (size_t)system->num_modes * sizeof(double));

            if(ROM_std_stdd_stsb_evaluate_lifted_reference_residual(
                adapter,
                system,
                window,
                initial_trace,
                current_reference,
                previous_reference,
                direct_projected_defect,
                &fom_residual_sum) != ROM_STDD_SUCCESS)
            {
                fail_message("direct lifted FOM residual evaluation failed");
            }

            if(ROM_std_stdd_evaluate_local_reduced_equation(
                system,
                &sys->monolis_com,
                reference_reduced_coef,
                window,
                assembled_defect,
                self_action,
                neighbor_action,
                temporal_action) != ROM_STDD_SUCCESS)
            {
                fail_message("assembled reduced equation diagnostic failed");
            }

            diagnostic->fom_residual_all2 =
                fom_residual_sum.residual_all2;
            diagnostic->fom_rhs_all2 =
                fom_residual_sum.effective_rhs_all2;
            diagnostic->fom_residual_free2 =
                fom_residual_sum.residual_free2;
            diagnostic->fom_rhs_free2 =
                fom_residual_sum.effective_rhs_free2;
            diagnostic->fom_residual_dirichlet2 =
                fom_residual_sum.residual_dirichlet2;
            diagnostic->fom_rhs_dirichlet2 =
                fom_residual_sum.effective_rhs_dirichlet2;
            diagnostic->homogeneous_dirichlet2 =
                fom_residual_sum.homogeneous_dirichlet2;
            diagnostic->homogeneous_all2 =
                fom_residual_sum.homogeneous_all2;

            for(int mode = 0; mode < system->num_modes; mode++){
                const double identity_difference =
                    assembled_defect[mode]
                    - direct_projected_defect[mode];

                diagnostic->direct_projected_defect2 +=
                    direct_projected_defect[mode]
                    * direct_projected_defect[mode];
                diagnostic->assembled_defect2 +=
                    assembled_defect[mode] * assembled_defect[mode];
                diagnostic->identity_difference2 +=
                    identity_difference * identity_difference;
                diagnostic->reduced_rhs2 +=
                    reduced_rhs[mode] * reduced_rhs[mode];
                diagnostic->self_action2 +=
                    self_action[mode] * self_action[mode];
                diagnostic->neighbor_action2 +=
                    neighbor_action[mode] * neighbor_action[mode];
                diagnostic->temporal_action2 +=
                    temporal_action[mode] * temporal_action[mode];
            }

            add_residual_diagnostic_sum(
                &residual_diagnostic_sum,
                diagnostic);
        }

        reduced_residual =
            ROM_std_stdd_system_reduced_residual_from_coefficients(
                system,
                &sys->monolis_com,
                system->reduced_coef,
                rom_window_residuals);

        reference_reduced_defect =
            ROM_std_stdd_system_reduced_residual_from_coefficients(
                system,
                &sys->monolis_com,
                reference_reduced_coef,
                reference_window_defects);

        print_validation_window_table(
            window_validation_sums,
            system->num_windows,
            rom_window_residuals,
            reference_window_defects);

        print_residual_diagnostics_by_window(
            window_residual_diagnostics,
            system->num_windows);

        print_validation_summary(
            &validation_sum,
            &residual_diagnostic_sum,
            reduced_residual,
            reference_reduced_defect);
    }
    else if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-DDROM reduced algebraic residual: %.15e\n",
            reduced_residual);
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-DDROM online completed: windows=%d, ranks=%d, "
            "modes/rank=%d, global reduced dof=%d\n",
            system->num_windows,
            monolis_mpi_get_global_comm_size(),
            system->num_modes,
            system->num_windows
                * monolis_mpi_get_global_comm_size()
                * system->num_modes);
    }

    if(options->reduced_solver_mode ==
       ST_DDROM_REDUCED_ALL_AT_ONCE)
    {
        ROM_std_stdd_reduced_solver_finalize(&reduced_solver);
    }
    free(initial_trace);
    free(effective_rhs);
    free(homogeneous_work);
    free(full_window_state);
    free(reference_snapshots);
    free(coefficient_reference);
    free(reference_reduced_coef);
    free(rom_window_residuals);
    free(reference_window_defects);
    free(window_validation_sums);
    free(q_projection);
    free(trace_rom);
    free(trace_reference);
    free(direct_projected_defect);
    free(assembled_defect);
    free(self_action);
    free(neighbor_action);
    free(temporal_action);
    free(window_residual_diagnostics);
}

int main(int argc, char* argv[])
{
    FE_SYSTEM sys = {0};
    ST_TIME_ELEMENT time_element;
    STSB_SPATIAL_GRAPH graph;
    ROM_STDD_SYSTEM rom_system;
    ROM_STDD_STSB_CONTEXT adapter;
    ST_DDROM_OPTIONS options;
    const int standard_ddrom =
        env_int("ST_DDROM_STANDARD_DDROM", 0);
    const int degree = standard_ddrom
        ? 0
        : env_int("ST_DDROM_TIME_DEGREE", ST_DDROM_DEGREE);
    const int num_window_slabs = standard_ddrom
        ? 1
        : env_int("ST_DDROM_SLABS_PER_WINDOW", ST_DDROM_WINDOW_SLABS);
    int num_windows = 0;
    int num_internal_space_nodes = 0;
    int st_dof_per_space_node;
    int rank;
    double time_start;
    double time_end;

    memset(&graph, 0, sizeof(graph));
    memset(&rom_system, 0, sizeof(rom_system));
    memset(&adapter, 0, sizeof(adapter));

    monolis_global_initialize();
    rank = monolis_mpi_get_global_my_rank();
    time_start = monolis_get_time();

    if(standard_ddrom > 1 || degree < 0 || num_window_slabs <= 0){
        fail_message("invalid ST-DDROM time discretization option");
    }

    STSB_initialize_reference_fom(
        &sys,
        &time_element,
        &graph,
        argc,
        argv,
        degree,
        num_window_slabs,
        &num_windows,
        &num_internal_space_nodes);

    read_options(&options, sys.cond.directory);

    if(options.write_output && options.mode != ST_DDROM_RUN_OFFLINE){
        initialize_output(&sys, rank);
    }

    st_dof_per_space_node = STSB_get_num_dofs_per_node(
        &time_element,
        num_window_slabs);

    /* dG(0), one slab/window is the ordinary time-marching DDROM limit. */
    if(standard_ddrom && st_dof_per_space_node != 1){
        fail_message(
            "standard DDROM mode requires exactly one time DOF per space node");
    }

    if(ROM_std_stdd_system_initialize(
        &rom_system,
        &sys.monolis_com,
        num_windows,
        sys.fe.total_num_nodes,
        num_internal_space_nodes,
        st_dof_per_space_node,
        options.max_modes) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_context_initialize(
        &adapter,
        &sys,
        &time_element,
        num_window_slabs,
        num_windows) != ROM_STDD_SUCCESS)
    {
        fail_message("ST-DDROM initialization failed");
    }

    if(rank == 0){
        printf("ST-DDROM run mode: %s\n", mode_name(options.mode));
        printf(
            "formulation: %s\n",
            standard_ddrom
                ? "ordinary spatial DDROM limit: dG(0), one slab/window"
                : "space-time DDROM");
        printf("time dG degree: %d\n", degree);
        printf("slabs/window: %d\n", num_window_slabs);
        printf("num windows: %d\n", num_windows);
        printf("spatial MPI ranks: %d\n", monolis_mpi_get_global_comm_size());
        printf("max modes/rank: %d\n", options.max_modes);
        printf(
            "reduced solver: %s\n",
            reduced_solver_mode_name(options.reduced_solver_mode));
        printf("snapshot directory: %s\n", options.snapshot_dir);
        printf("basis directory: %s\n", options.basis_dir);
        printf("ROM operator form: FOM-consistent strong-BC A plus affine lifting\n");
        printf(
            "external FOM residual check: %s, strict rejection: %s\n",
            options.check_fom_residual ? "enabled" : "disabled",
            options.strict_fom_residual ? "enabled" : "disabled");
    }

    switch(options.mode){
    case ST_DDROM_RUN_FOM:
    case ST_DDROM_RUN_COLLECT:
        run_fom_or_collect(
            &sys,
            &time_element,
            num_window_slabs,
            num_windows,
            num_internal_space_nodes,
            &options,
            &adapter,
            &rom_system);
        break;

    case ST_DDROM_RUN_OFFLINE:
        run_offline(&rom_system, &adapter, &options);
        break;

    case ST_DDROM_RUN_ONLINE:
        run_online(
            &sys,
            &time_element,
            num_window_slabs,
            &rom_system,
            &adapter,
            &options);
        break;
    }

    ROM_std_stdd_stsb_context_finalize(&adapter);
    ROM_std_stdd_system_finalize(&rom_system);
    STSB_finalize_reference_fom(&sys, &graph);

    time_end = monolis_get_time();
    if(rank == 0){
        printf("** Total time: %f\n\n", time_end - time_start);
    }

    monolis_global_finalize();
    return 0;
}
