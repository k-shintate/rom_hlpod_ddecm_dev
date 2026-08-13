#include "core_ROM.h"
#include "core_FOM_ST_spaceblock.h"
#include "st_fom_common.h"
#include "rom_std_stddrom.h"
#include "rom_std_stddrom_stsb.h"
#include "rom_std_strom.h"

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

typedef enum {
    ST_DDROM_FORMULATION_CURRENT_STATE = 0,
    ST_DDROM_FORMULATION_INCREMENTAL
} ST_DDROM_FORMULATION;

/*
 * Reference-FOM formulation contract.
 *
 * The present STSB implementation solves for the physical current state U^n,
 * so STSB_reference_solution_is_incremental() returns 0.  An incremental
 * reference implementation must return 1 from the same function.
 */
#ifdef __cplusplus
extern "C" {
#endif
int STSB_reference_solution_is_incremental(void);
#ifdef __cplusplus
}
#endif

typedef struct {
    ST_DDROM_RUN_MODE mode;
    int snapshot_id;
    int num_snapshots;
    int max_modes;
    int validate;
    int validation_snapshot_id;
    int write_output;
    int write_comparison_vtk;
    int write_mode_vtk;
    int mode_vtk_max_modes;
    int compare_dg0_standard;
    int compare_hlpod_kernels;
    int compare_strom_reference;
    double compare_tolerance;
    double strom_compare_tolerance;
    double hlpod_compare_tolerance;
    char compare_csv[4096];
    char hlpod_compare_csv[4096];
    char strom_compare_csv[4096];
    ST_DDROM_FORMULATION formulation;
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

static const char* formulation_name(ST_DDROM_FORMULATION formulation)
{
    switch(formulation){
    case ST_DDROM_FORMULATION_CURRENT_STATE:
        return "current_state";
    case ST_DDROM_FORMULATION_INCREMENTAL:
        return "incremental";
    }
    return "unknown";
}

static ST_DDROM_FORMULATION reference_fom_formulation(void)
{
    return STSB_reference_solution_is_incremental()
        ? ST_DDROM_FORMULATION_INCREMENTAL
        : ST_DDROM_FORMULATION_CURRENT_STATE;
}

static ST_DDROM_FORMULATION read_formulation(void)
{
    const char* text = getenv("ST_DDROM_FOM_FORMULATION");

    if(text == NULL || text[0] == '\0' || strcmp(text, "auto") == 0){
        return reference_fom_formulation();
    }

    if(strcmp(text, "current") == 0 ||
       strcmp(text, "current_state") == 0 ||
       strcmp(text, "absolute") == 0)
    {
        return ST_DDROM_FORMULATION_CURRENT_STATE;
    }

    if(strcmp(text, "increment") == 0 ||
       strcmp(text, "incremental") == 0 ||
       strcmp(text, "delta") == 0)
    {
        return ST_DDROM_FORMULATION_INCREMENTAL;
    }

    fprintf(
        stderr,
        "ERROR: unknown ST_DDROM_FOM_FORMULATION=%s\n"
        "       supported values: auto, current_state, incremental\n",
        text);
    exit(EXIT_FAILURE);
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
    options->write_comparison_vtk = env_int(
        "ST_DDROM_WRITE_COMPARISON_VTK",
        options->write_output);
    options->write_mode_vtk =
        env_int("ST_DDROM_WRITE_MODE_VTK", 1);
    options->mode_vtk_max_modes =
        env_int("ST_DDROM_MODE_VTK_MAX_MODES", 0);
    options->compare_dg0_standard =
        env_int("ST_DDROM_COMPARE_STANDARD_DDROM", 0);
    options->compare_hlpod_kernels =
        env_int("ST_DDROM_COMPARE_HLPOD_KERNELS", 0);
    options->compare_strom_reference =
        env_int("ST_DDROM_COMPARE_STROM", 0);
    options->compare_tolerance =
        env_double("ST_DDROM_COMPARE_TOL", 1.0e-10);
    options->hlpod_compare_tolerance =
        env_double("ST_DDROM_HLPOD_COMPARE_TOL",
                   options->compare_tolerance);
    options->strom_compare_tolerance =
        env_double("ST_DDROM_STROM_COMPARE_TOL", 1.0e-8);
    options->formulation = read_formulation();
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
       options->write_comparison_vtk > 1 ||
       options->write_mode_vtk > 1 ||
       options->compare_dg0_standard > 1 ||
       options->compare_hlpod_kernels > 1 ||
       options->compare_strom_reference > 1 ||
       options->compare_tolerance <= 0.0 ||
       options->hlpod_compare_tolerance <= 0.0 ||
       options->strom_compare_tolerance <= 0.0 ||
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

    {
        const char* compare_csv_override =
            getenv("ST_DDROM_COMPARE_CSV");

        if(compare_csv_override != NULL &&
           compare_csv_override[0] != '\0')
        {
            if(snprintf(
                options->compare_csv,
                sizeof(options->compare_csv),
                "%s",
                compare_csv_override) >=
               (int)sizeof(options->compare_csv))
            {
                fail_message("comparison CSV path is too long");
            }
        }
        else{
            join_path(
                options->compare_csv,
                sizeof(options->compare_csv),
                case_directory,
                "dg0_stddrom_vs_standard_ddrom.csv");
        }
    }

    {
        const char* hlpod_csv_override =
            getenv("ST_DDROM_HLPOD_COMPARE_CSV");

        if(hlpod_csv_override != NULL &&
           hlpod_csv_override[0] != '\0')
        {
            if(snprintf(
                options->hlpod_compare_csv,
                sizeof(options->hlpod_compare_csv),
                "%s",
                hlpod_csv_override) >=
               (int)sizeof(options->hlpod_compare_csv))
            {
                fail_message("HLPOD comparison CSV path is too long");
            }
        }
        else{
            join_path(
                options->hlpod_compare_csv,
                sizeof(options->hlpod_compare_csv),
                case_directory,
                "dg0_stddrom_vs_hlpod_kernels.csv");
        }
    }


    {
        const char* strom_csv_override =
            getenv("ST_DDROM_STROM_COMPARE_CSV");

        if(strom_csv_override != NULL &&
           strom_csv_override[0] != '\0')
        {
            if(snprintf(
                options->strom_compare_csv,
                sizeof(options->strom_compare_csv),
                "%s",
                strom_csv_override) >=
               (int)sizeof(options->strom_compare_csv))
            {
                fail_message("sequential ST-ROM comparison CSV path is too long");
            }
        }
        else{
            join_path(
                options->strom_compare_csv,
                sizeof(options->strom_compare_csv),
                case_directory,
                "stddrom_vs_sequential_strom.csv");
        }
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


static int make_directory_tree_local(const char* directory)
{
    char path[4096];
    size_t length;

    if(directory == NULL || directory[0] == '\0'){
        return 0;
    }

    if(snprintf(path, sizeof(path), "%s", directory) >= (int)sizeof(path)){
        errno = ENAMETOOLONG;
        return 1;
    }

    length = strlen(path);
    while(length > 1 && path[length - 1] == '/'){
        path[--length] = '\0';
    }

    for(char* cursor = path + 1; *cursor != '\0'; cursor++){
        if(*cursor != '/'){
            continue;
        }

        *cursor = '\0';
        if(path[0] != '\0' &&
           mkdir(path, 0775) != 0 &&
           errno != EEXIST)
        {
            *cursor = '/';
            return 1;
        }
        *cursor = '/';
    }

    if(mkdir(path, 0775) != 0 && errno != EEXIST){
        return 1;
    }

    return 0;
}

static void make_parent_directory_collective(const char* filename)
{
    char parent[4096];
    const char* slash;
    int local_status = 0;
    int global_status;

    if(filename == NULL || filename[0] == '\0'){
        local_status = 1;
    }
    else{
        slash = strrchr(filename, '/');
        if(slash != NULL){
            const size_t length = slash == filename
                ? 1u
                : (size_t)(slash - filename);

            if(length >= sizeof(parent)){
                local_status = 1;
                errno = ENAMETOOLONG;
            }
            else{
                memcpy(parent, filename, length);
                parent[length] = '\0';

                if(make_directory_tree_local(parent) != 0){
                    fprintf(
                        stderr,
                        "ERROR [rank %d]: cannot create parent directory %s: %s\n",
                        monolis_mpi_get_global_my_rank(),
                        parent,
                        strerror(errno));
                    fflush(stderr);
                    local_status = 1;
                }
            }
        }
    }

    global_status = local_status;
    monolis_allreduce_I(
        1,
        &global_status,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    if(global_status != 0){
        fail_message(
            "sequential ST-ROM comparison CSV parent directory creation failed");
    }
}

static void write_formulation_manifest(
    const char* directory,
    ST_DDROM_FORMULATION formulation)
{
    char filename[4096];
    int status = 0;

    if(directory == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/stddrom_formulation.txt",
           directory) >= (int)sizeof(filename))
    {
        fail_message("formulation manifest path is too long");
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp = fopen(filename, "w");

        if(fp == NULL){
            status = 1;
        }
        else{
            if(fprintf(
                fp,
                "version 1\nformulation %s\n",
                formulation_name(formulation)) < 0)
            {
                status = 1;
            }
            fclose(fp);
        }
    }

    monolis_allreduce_I(
        1,
        &status,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    if(status != 0){
        fail_message("failed to write formulation manifest");
    }
}

static void require_formulation_manifest(
    const char* directory,
    ST_DDROM_FORMULATION expected)
{
    char filename[4096];
    char key[64];
    char value[64];
    int version = 0;
    int status = 0;
    ST_DDROM_FORMULATION actual = ST_DDROM_FORMULATION_CURRENT_STATE;

    if(directory == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/stddrom_formulation.txt",
           directory) >= (int)sizeof(filename))
    {
        fail_message("formulation manifest path is too long");
    }

    {
        FILE* fp = fopen(filename, "r");

        if(fp == NULL ||
           fscanf(fp, "%63s %d", key, &version) != 2 ||
           strcmp(key, "version") != 0 || version != 1 ||
           fscanf(fp, "%63s %63s", key, value) != 2 ||
           strcmp(key, "formulation") != 0)
        {
            status = 1;
        }
        else if(strcmp(value, "current_state") == 0){
            actual = ST_DDROM_FORMULATION_CURRENT_STATE;
        }
        else if(strcmp(value, "incremental") == 0){
            actual = ST_DDROM_FORMULATION_INCREMENTAL;
        }
        else{
            status = 1;
        }

        if(fp != NULL){
            fclose(fp);
        }
    }

    if(status == 0 && actual != expected){
        fprintf(
            stderr,
            "ERROR [rank %d]: formulation mismatch in %s\n"
            "  file     : %s\n"
            "  requested: %s\n",
            monolis_mpi_get_global_my_rank(),
            filename,
            formulation_name(actual),
            formulation_name(expected));
        fflush(stderr);
        status = 1;
    }

    monolis_allreduce_I(
        1,
        &status,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    if(status != 0){
        fail_message(
            "snapshot/basis formulation metadata is missing or incompatible; "
            "rebuild snapshots and basis");
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

typedef struct {
    double residual_all2;
    double rhs_all2;
    double ax_all2;

    double residual_free2;
    double rhs_free2;
    double ax_free2;

    double residual_dirichlet2;
    double rhs_dirichlet2;
    double ax_dirichlet2;
} ST_DDROM_ABSOLUTE_RESIDUAL_SUM;

/*
 * Evaluate the current strong-BC MONOLIS pair without rebuilding it.
 *
 * This diagnostic intentionally records both the conventional
 *
 *     ||A x - b|| / ||b||
 *
 * and the scale-robust quantity
 *
 *     ||A x - b|| / (||A x|| + ||b||).
 *
 * The latter remains meaningful when the free-row right-hand side is small.
 */
static int evaluate_current_strong_pair_residual(
    FE_SYSTEM* sys,
    int dof_per_space_node,
    const double* solution,
    ST_DDROM_ABSOLUTE_RESIDUAL_SUM* sum)
{
    int total_nodes;
    int internal_nodes;
    int local_dof;
    int internal_dof;
    double* solution_work = NULL;
    double* product = NULL;
    double values[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if(sys == NULL || solution == NULL || sum == NULL ||
       dof_per_space_node <= 0 || sys->monolis.mat.R.B == NULL)
    {
        return -1;
    }

    total_nodes = sys->fe.total_num_nodes;
    internal_nodes = sys->monolis.mat.N;

    if(total_nodes <= 0 || internal_nodes <= 0 ||
       internal_nodes > total_nodes ||
       total_nodes > INT_MAX / dof_per_space_node ||
       internal_nodes > INT_MAX / dof_per_space_node)
    {
        return -1;
    }

    local_dof = total_nodes * dof_per_space_node;
    internal_dof = internal_nodes * dof_per_space_node;

    solution_work = (double*)calloc((size_t)local_dof, sizeof(double));
    product = (double*)calloc((size_t)local_dof, sizeof(double));

    if(solution_work == NULL || product == NULL){
        free(solution_work);
        free(product);
        return -1;
    }

    memcpy(
        solution_work,
        solution,
        (size_t)local_dof * sizeof(double));

    monolis_mpi_update_R(
        &sys->monolis_com,
        total_nodes,
        dof_per_space_node,
        solution_work);

    monolis_matvec_product_R(
        &sys->monolis,
        &sys->monolis_com,
        solution_work,
        product);

    for(int row = 0; row < internal_dof; row++){
        const int node = row / dof_per_space_node;
        const double ax = product[row];
        const double rhs = sys->monolis.mat.R.B[row];
        const double residual = ax - rhs;
        const int is_dirichlet = sys->bc.D_bc_exists[node] ? 1 : 0;

        values[0] += residual * residual;
        values[1] += rhs * rhs;
        values[2] += ax * ax;

        if(is_dirichlet){
            values[6] += residual * residual;
            values[7] += rhs * rhs;
            values[8] += ax * ax;
        }
        else{
            values[3] += residual * residual;
            values[4] += rhs * rhs;
            values[5] += ax * ax;
        }
    }

    monolis_allreduce_R(
        9,
        values,
        MONOLIS_MPI_SUM,
        sys->monolis_com.comm);

    sum->residual_all2 = values[0];
    sum->rhs_all2 = values[1];
    sum->ax_all2 = values[2];
    sum->residual_free2 = values[3];
    sum->rhs_free2 = values[4];
    sum->ax_free2 = values[5];
    sum->residual_dirichlet2 = values[6];
    sum->rhs_dirichlet2 = values[7];
    sum->ax_dirichlet2 = values[8];

    free(solution_work);
    free(product);
    return 0;
}

static double balanced_relative_residual(
    double residual2,
    double ax2,
    double rhs2)
{
    if(residual2 < 0.0 || ax2 < 0.0 || rhs2 < 0.0 ||
       !isfinite(residual2) || !isfinite(ax2) || !isfinite(rhs2))
    {
        return INFINITY;
    }

    return sqrt(residual2) /
        fmax(sqrt(ax2) + sqrt(rhs2), 1.0e-300);
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



/*
 * BB_vtk_write_point_vals_scalar() has a legacy non-const double* API,
 * although it only writes the supplied values.  Keep comparison data const
 * in this driver and isolate the compatibility cast in one wrapper.
 */
static void write_vtk_point_scalar_const(
    FILE* fp,
    const double* values,
    int nvalue,
    const char* label)
{
    BB_vtk_write_point_vals_scalar(
        fp,
        (double*)values,
        nvalue,
        label);
}



static void write_pod_mode_piece_vtk(
    FE_SYSTEM* sys,
    const char* piece_name,
    double* component,
    double* component_normalized,
    double* component_absolute,
    double* space_time_magnitude,
    double* space_time_magnitude_normalized,
    double* dirichlet_mask,
    double* owned_node_mask,
    double* source_owned_support_mask,
    double* partition_rank_field,
    double* source_rank_field,
    double* local_mode_field,
    double* temporal_dof_field,
    double* singular_value_field)
{
    FILE* fp = NULL;
    const int nnode = sys->fe.total_num_nodes;

    fp = ROM_BB_write_fopen(
        fp,
        piece_name,
        sys->cond.directory);

    if(fp == NULL){
        fail_message("POD-mode VTK piece open failed");
    }

    switch(sys->fe.local_num_nodes){
    case 4:
        BBFE_sys_write_vtk_shape(
            fp,
            &sys->fe,
            TYPE_VTK_TETRA);
        break;
    case 8:
        BBFE_sys_write_vtk_shape(
            fp,
            &sys->fe,
            TYPE_VTK_HEXAHEDRON);
        break;
    default:
        fclose(fp);
        fail_message(
            "POD-mode VTK supports tetrahedral or hexahedral meshes");
    }

    fprintf(fp, "POINT_DATA %d\n", nnode);
    BB_vtk_write_point_vals_scalar(
        fp,
        component,
        nnode,
        "POD_mode_component_raw");
    BB_vtk_write_point_vals_scalar(
        fp,
        component_normalized,
        nnode,
        "POD_mode_component_normalized");
    BB_vtk_write_point_vals_scalar(
        fp,
        component_absolute,
        nnode,
        "absolute_POD_mode_component");
    BB_vtk_write_point_vals_scalar(
        fp,
        space_time_magnitude,
        nnode,
        "POD_mode_space_time_magnitude");
    BB_vtk_write_point_vals_scalar(
        fp,
        space_time_magnitude_normalized,
        nnode,
        "POD_mode_space_time_magnitude_normalized");
    BB_vtk_write_point_vals_scalar(
        fp,
        dirichlet_mask,
        nnode,
        "Dirichlet_mask");
    BB_vtk_write_point_vals_scalar(
        fp,
        owned_node_mask,
        nnode,
        "owned_node_mask");
    BB_vtk_write_point_vals_scalar(
        fp,
        source_owned_support_mask,
        nnode,
        "source_owned_support_mask");
    BB_vtk_write_point_vals_scalar(
        fp,
        partition_rank_field,
        nnode,
        "partition_rank");
    BB_vtk_write_point_vals_scalar(
        fp,
        source_rank_field,
        nnode,
        "source_rank");
    BB_vtk_write_point_vals_scalar(
        fp,
        local_mode_field,
        nnode,
        "local_mode_index_zero_based");
    BB_vtk_write_point_vals_scalar(
        fp,
        temporal_dof_field,
        nnode,
        "temporal_dof_index");
    BB_vtk_write_point_vals_scalar(
        fp,
        singular_value_field,
        nnode,
        "singular_value");

    fclose(fp);
}


static void write_pod_mode_vtm_master(
    FE_SYSTEM* sys,
    const char* base_name,
    int comm_size)
{
    FILE* fp = NULL;
    char master_name[4096];

    if(monolis_mpi_get_global_my_rank() != 0){
        return;
    }

    if(snprintf(
        master_name,
        sizeof(master_name),
        "%s.vtm",
        base_name) >= (int)sizeof(master_name))
    {
        fail_message("POD-mode VTM filename is too long");
    }

    fp = ROM_BB_write_fopen(
        fp,
        master_name,
        sys->cond.directory);

    if(fp == NULL){
        fail_message("POD-mode VTM master open failed");
    }

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(
        fp,
        "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
        "byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <vtkMultiBlockDataSet>\n");

    for(int piece_rank = 0; piece_rank < comm_size; piece_rank++){
        fprintf(
            fp,
            "    <DataSet index=\"%d\" name=\"rank_%04d\" "
            "file=\"%s_piece_%04d.vtk\"/>\n",
            piece_rank,
            piece_rank,
            base_name,
            piece_rank);
    }

    fprintf(fp, "  </vtkMultiBlockDataSet>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}


static void output_pod_modes_vtk(
    FE_SYSTEM* sys,
    const ROM_STDD_SYSTEM* system,
    const ST_DDROM_OPTIONS* options)
{
    const int rank = monolis_mpi_get_global_my_rank();
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int nnode = (sys != NULL) ? sys->fe.total_num_nodes : 0;
    const int com_internal_nodes =
        (sys != NULL) ? sys->monolis_com.n_internal_vertex : 0;
    const int basis_internal_nodes =
        (system != NULL) ? system->num_internal_space_nodes : 0;
    FILE* summary_fp = NULL;
    double* component = NULL;
    double* component_normalized = NULL;
    double* component_absolute = NULL;
    double* magnitude = NULL;
    double* magnitude_normalized = NULL;
    double* dirichlet_mask = NULL;
    double* owned_node_mask = NULL;
    double* source_owned_support_mask = NULL;
    double* partition_rank_field = NULL;
    double* source_rank_field = NULL;
    double* local_mode_field = NULL;
    double* temporal_dof_field = NULL;
    double* singular_value_field = NULL;
    int max_modes_written = 0;

    if(sys == NULL || system == NULL || options == NULL ||
       system->basis == NULL || system->num_modes <= 0 ||
       system->num_modes > system->num_modes_capacity ||
       system->st_dof_per_space_node <= 0 || nnode <= 0 ||
       basis_internal_nodes <= 0 || basis_internal_nodes > nnode ||
       system->local_st_dof !=
           nnode * system->st_dof_per_space_node ||
       system->internal_st_dof !=
           basis_internal_nodes * system->st_dof_per_space_node)
    {
        fprintf(
            stderr,
            "ERROR [rank %d]: invalid POD-mode VTK dimensions: "
            "nnode=%d basis_internal=%d com_internal=%d "
            "ndof=%d local_st_dof=%d internal_st_dof=%d "
            "num_modes=%d capacity=%d\n",
            rank,
            nnode,
            basis_internal_nodes,
            com_internal_nodes,
            (system != NULL) ? system->st_dof_per_space_node : -1,
            (system != NULL) ? system->local_st_dof : -1,
            (system != NULL) ? system->internal_st_dof : -1,
            (system != NULL) ? system->num_modes : -1,
            (system != NULL) ? system->num_modes_capacity : -1);
        fflush(stderr);
        fail_message("invalid POD-mode VTK dimensions");
    }

    /*
     * The ROM basis is stored on the first
     * system->num_internal_space_nodes local nodes.  This count is obtained
     * from the matrix used when snapshots and basis rows are constructed and
     * is therefore authoritative for indexing system->basis.
     *
     * Some MONOLIS configurations keep a different n_internal_vertex value in
     * the communication table (for example after constructing a space-time
     * matrix or a specialized communication object).  That difference must
     * not change basis-row indexing or prevent visualization.
     */
    if(com_internal_nodes > 0 &&
       com_internal_nodes != basis_internal_nodes)
    {
        fprintf(
            stderr,
            "POD-mode VTK note [rank %d]: basis-owned nodes=%d, "
            "communication-owned nodes=%d, local nodes=%d; "
            "using basis-owned count for mode rows.\n",
            rank,
            basis_internal_nodes,
            com_internal_nodes,
            nnode);
        fflush(stderr);
    }

    component = (double*)calloc((size_t)nnode, sizeof(double));
    component_normalized =
        (double*)calloc((size_t)nnode, sizeof(double));
    component_absolute =
        (double*)calloc((size_t)nnode, sizeof(double));
    magnitude = (double*)calloc((size_t)nnode, sizeof(double));
    magnitude_normalized =
        (double*)calloc((size_t)nnode, sizeof(double));
    dirichlet_mask = (double*)calloc((size_t)nnode, sizeof(double));
    owned_node_mask = (double*)calloc((size_t)nnode, sizeof(double));
    source_owned_support_mask =
        (double*)calloc((size_t)nnode, sizeof(double));
    partition_rank_field =
        (double*)calloc((size_t)nnode, sizeof(double));
    source_rank_field =
        (double*)calloc((size_t)nnode, sizeof(double));
    local_mode_field =
        (double*)calloc((size_t)nnode, sizeof(double));
    temporal_dof_field =
        (double*)calloc((size_t)nnode, sizeof(double));
    singular_value_field =
        (double*)calloc((size_t)nnode, sizeof(double));

    if(component == NULL || component_normalized == NULL ||
       component_absolute == NULL || magnitude == NULL ||
       magnitude_normalized == NULL || dirichlet_mask == NULL ||
       owned_node_mask == NULL || source_owned_support_mask == NULL ||
       partition_rank_field == NULL || source_rank_field == NULL ||
       local_mode_field == NULL || temporal_dof_field == NULL ||
       singular_value_field == NULL)
    {
        fail_message("POD-mode VTK allocation failed");
    }

    for(int node = 0; node < nnode; node++){
        dirichlet_mask[node] =
            sys->bc.D_bc_exists[node] ? 1.0 : 0.0;
        owned_node_mask[node] =
            node < basis_internal_nodes ? 1.0 : 0.0;
        partition_rank_field[node] = (double)rank;
    }

    if(rank == 0){
        summary_fp = ROM_BB_write_fopen(
            summary_fp,
            "stddrom_mode_summary.csv",
            sys->cond.directory);
        if(summary_fp == NULL){
            fail_message("POD-mode summary open failed");
        }
        fprintf(
            summary_fp,
            "source_rank,local_mode_zero_based,temporal_dof,"
            "component_norm,total_space_time_norm,energy_fraction,"
            "component_max_abs,space_time_max_abs,singular_value\n");
    }

    for(int source_rank = 0; source_rank < comm_size; source_rank++){
        int source_modes =
            (rank == source_rank) ? system->num_modes : 0;

        monolis_allreduce_I(
            1,
            &source_modes,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);

        if(options->mode_vtk_max_modes > 0 &&
           options->mode_vtk_max_modes < source_modes)
        {
            source_modes = options->mode_vtk_max_modes;
        }

        if(source_modes > max_modes_written){
            max_modes_written = source_modes;
        }

        for(int mode = 0; mode < source_modes; mode++){
            double singular_value =
                (rank == source_rank &&
                 system->singular_values != NULL &&
                 mode < system->num_singular_values)
                    ? system->singular_values[mode]
                    : 0.0;
            double total_mode_norm2 = 0.0;
            double global_magnitude_max = 0.0;

            monolis_allreduce_R(
                1,
                &singular_value,
                MONOLIS_MPI_SUM,
                sys->monolis_com.comm);

            memset(magnitude, 0, (size_t)nnode * sizeof(double));
            memset(
                source_owned_support_mask,
                0,
                (size_t)nnode * sizeof(double));

            if(rank == source_rank){
                for(int node = 0; node < basis_internal_nodes; node++){
                    double value2 = 0.0;

                    for(int temporal_dof = 0;
                        temporal_dof < system->st_dof_per_space_node;
                        temporal_dof++)
                    {
                        const int row =
                            node * system->st_dof_per_space_node
                            + temporal_dof;
                        const double value = system->basis[
                            (size_t)row
                                * (size_t)system->num_modes_capacity
                            + (size_t)mode];

                        value2 += value * value;
                    }

                    magnitude[node] = sqrt(value2);
                    source_owned_support_mask[node] = 1.0;
                    total_mode_norm2 += value2;
                }
            }

            monolis_allreduce_R(
                1,
                &total_mode_norm2,
                MONOLIS_MPI_SUM,
                sys->monolis_com.comm);

            monolis_mpi_update_R(
                &sys->monolis_com,
                nnode,
                1,
                magnitude);

            for(int node = 0; node < nnode; node++){
                global_magnitude_max = fmax(
                    global_magnitude_max,
                    magnitude[node]);
            }

            monolis_allreduce_R(
                1,
                &global_magnitude_max,
                MONOLIS_MPI_MAX,
                sys->monolis_com.comm);

            for(int node = 0; node < nnode; node++){
                magnitude_normalized[node] =
                    global_magnitude_max > 0.0
                        ? magnitude[node] / global_magnitude_max
                        : 0.0;
                source_rank_field[node] = (double)source_rank;
                local_mode_field[node] = (double)mode;
                singular_value_field[node] = singular_value;
            }

            for(int temporal_dof = 0;
                temporal_dof < system->st_dof_per_space_node;
                temporal_dof++)
            {
                char base_name[4096];
                char piece_name[4096];
                double component_norm2 = 0.0;
                double global_component_max = 0.0;

                memset(component, 0, (size_t)nnode * sizeof(double));

                if(rank == source_rank){
                    for(int node = 0; node < basis_internal_nodes; node++){
                        const int row =
                            node * system->st_dof_per_space_node
                            + temporal_dof;
                        const double value = system->basis[
                            (size_t)row
                                * (size_t)system->num_modes_capacity
                            + (size_t)mode];

                        component[node] = value;
                        component_norm2 += value * value;
                    }
                }

                monolis_allreduce_R(
                    1,
                    &component_norm2,
                    MONOLIS_MPI_SUM,
                    sys->monolis_com.comm);

                monolis_mpi_update_R(
                    &sys->monolis_com,
                    nnode,
                    1,
                    component);

                for(int node = 0; node < nnode; node++){
                    global_component_max = fmax(
                        global_component_max,
                        fabs(component[node]));
                }

                monolis_allreduce_R(
                    1,
                    &global_component_max,
                    MONOLIS_MPI_MAX,
                    sys->monolis_com.comm);

                for(int node = 0; node < nnode; node++){
                    component_absolute[node] = fabs(component[node]);
                    component_normalized[node] =
                        global_component_max > 0.0
                            ? component[node] / global_component_max
                            : 0.0;
                    temporal_dof_field[node] = (double)temporal_dof;
                }

                if(snprintf(
                    base_name,
                    sizeof(base_name),
                    "stddrom_mode_sr_%04d_mode_%04d_tdof_%03d",
                    source_rank,
                    mode + 1,
                    temporal_dof) >= (int)sizeof(base_name) ||
                   snprintf(
                    piece_name,
                    sizeof(piece_name),
                    "%s_piece_%04d.vtk",
                    base_name,
                    rank) >= (int)sizeof(piece_name))
                {
                    fail_message("POD-mode VTK filename is too long");
                }

                write_pod_mode_piece_vtk(
                    sys,
                    piece_name,
                    component,
                    component_normalized,
                    component_absolute,
                    magnitude,
                    magnitude_normalized,
                    dirichlet_mask,
                    owned_node_mask,
                    source_owned_support_mask,
                    partition_rank_field,
                    source_rank_field,
                    local_mode_field,
                    temporal_dof_field,
                    singular_value_field);

                /* A global sync here guarantees every piece exists first. */
                (void)monolis_get_time_global_sync();
                write_pod_mode_vtm_master(
                    sys,
                    base_name,
                    comm_size);

                if(rank == 0 && summary_fp != NULL){
                    fprintf(
                        summary_fp,
                        "%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                        source_rank,
                        mode,
                        temporal_dof,
                        sqrt(fmax(component_norm2, 0.0)),
                        sqrt(fmax(total_mode_norm2, 0.0)),
                        total_mode_norm2 > 0.0
                            ? component_norm2 / total_mode_norm2
                            : 0.0,
                        global_component_max,
                        global_magnitude_max,
                        singular_value);
                }
            }

            {
                char base_name[4096];
                char piece_name[4096];

                memset(component, 0, (size_t)nnode * sizeof(double));
                memset(
                    component_normalized,
                    0,
                    (size_t)nnode * sizeof(double));
                memset(
                    component_absolute,
                    0,
                    (size_t)nnode * sizeof(double));

                for(int node = 0; node < nnode; node++){
                    temporal_dof_field[node] = -1.0;
                }

                if(snprintf(
                    base_name,
                    sizeof(base_name),
                    "stddrom_mode_sr_%04d_mode_%04d_magnitude",
                    source_rank,
                    mode + 1) >= (int)sizeof(base_name) ||
                   snprintf(
                    piece_name,
                    sizeof(piece_name),
                    "%s_piece_%04d.vtk",
                    base_name,
                    rank) >= (int)sizeof(piece_name))
                {
                    fail_message("POD-mode magnitude filename is too long");
                }

                write_pod_mode_piece_vtk(
                    sys,
                    piece_name,
                    component,
                    component_normalized,
                    component_absolute,
                    magnitude,
                    magnitude_normalized,
                    dirichlet_mask,
                    owned_node_mask,
                    source_owned_support_mask,
                    partition_rank_field,
                    source_rank_field,
                    local_mode_field,
                    temporal_dof_field,
                    singular_value_field);

                (void)monolis_get_time_global_sync();
                write_pod_mode_vtm_master(
                    sys,
                    base_name,
                    comm_size);
            }
        }
    }

    if(summary_fp != NULL){
        fclose(summary_fp);
    }

    if(rank == 0){
        printf(
            "POD mode VTK completed: open the .vtm files, not an individual "
            "piece .vtk file. source ranks=%d, modes/rank=%d, time dofs=%d\n",
            comm_size,
            max_modes_written,
            system->st_dof_per_space_node);
        printf(
            "POD mode summary: stddrom_mode_summary.csv\n");
    }

    free(component);
    free(component_normalized);
    free(component_absolute);
    free(magnitude);
    free(magnitude_normalized);
    free(dirichlet_mask);
    free(owned_node_mask);
    free(source_owned_support_mask);
    free(partition_rank_field);
    free(source_rank_field);
    free(local_mode_field);
    free(temporal_dof_field);
    free(singular_value_field);
}

static void output_standard_comparison_vtk(
    FE_SYSTEM* sys,
    int global_step,
    double time,
    const double* fom_solution,
    const double* projected_fom_solution,
    const double* rom_solution)
{
    const int nnode = (sys != NULL) ? sys->fe.total_num_nodes : 0;
    char local_name[4096];
    const char* filename;
    FILE* fp = NULL;
    double* signed_rom_error = NULL;
    double* absolute_rom_error = NULL;
    double* relative_rom_error = NULL;
    double* signed_projection_error = NULL;
    double* absolute_projection_error = NULL;
    double* dirichlet_mask = NULL;
    double* time_value = NULL;
    double global_reference_scale = 0.0;
    double relative_floor;

    if(sys == NULL || fom_solution == NULL ||
       projected_fom_solution == NULL || rom_solution == NULL ||
       nnode <= 0 || global_step <= 0)
    {
        fail_message("invalid standard-DDROM comparison VTK argument");
    }

    signed_rom_error = (double*)calloc((size_t)nnode, sizeof(double));
    absolute_rom_error = (double*)calloc((size_t)nnode, sizeof(double));
    relative_rom_error = (double*)calloc((size_t)nnode, sizeof(double));
    signed_projection_error = (double*)calloc((size_t)nnode, sizeof(double));
    absolute_projection_error = (double*)calloc((size_t)nnode, sizeof(double));
    dirichlet_mask = (double*)calloc((size_t)nnode, sizeof(double));
    time_value = (double*)calloc((size_t)nnode, sizeof(double));

    if(signed_rom_error == NULL || absolute_rom_error == NULL ||
       relative_rom_error == NULL || signed_projection_error == NULL ||
       absolute_projection_error == NULL || dirichlet_mask == NULL ||
       time_value == NULL)
    {
        fail_message("comparison VTK allocation failed");
    }

    for(int node = 0; node < nnode; node++){
        global_reference_scale = fmax(
            global_reference_scale,
            fabs(fom_solution[node]));
    }

    monolis_allreduce_R(
        1,
        &global_reference_scale,
        MONOLIS_MPI_MAX,
        sys->monolis_com.comm);

    relative_floor = fmax(
        1.0e-14,
        1.0e-12 * global_reference_scale);

    for(int node = 0; node < nnode; node++){
        const double rom_difference =
            rom_solution[node] - fom_solution[node];
        const double projection_difference =
            projected_fom_solution[node] - fom_solution[node];

        signed_rom_error[node] = rom_difference;
        absolute_rom_error[node] = fabs(rom_difference);
        relative_rom_error[node] =
            fabs(rom_difference) /
            fmax(fabs(fom_solution[node]), relative_floor);
        signed_projection_error[node] = projection_difference;
        absolute_projection_error[node] = fabs(projection_difference);
        dirichlet_mask[node] =
            sys->bc.D_bc_exists[node] ? 1.0 : 0.0;
        time_value[node] = time;
    }

    if(snprintf(
        local_name,
        sizeof(local_name),
        "stddrom_compare_%06d.vtk",
        global_step) >= (int)sizeof(local_name))
    {
        fail_message("comparison VTK filename is too long");
    }

    filename = monolis_get_global_output_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        "./",
        local_name);

    fp = ROM_BB_write_fopen(
        fp,
        filename,
        sys->cond.directory);

    if(fp == NULL){
        fail_message("comparison VTK open failed");
    }

    switch(sys->fe.local_num_nodes){
    case 4:
        BBFE_sys_write_vtk_shape(fp, &sys->fe, TYPE_VTK_TETRA);
        break;
    case 8:
        BBFE_sys_write_vtk_shape(fp, &sys->fe, TYPE_VTK_HEXAHEDRON);
        break;
    default:
        fclose(fp);
        fail_message("comparison VTK supports tetrahedral or hexahedral meshes");
    }

    fprintf(fp, "POINT_DATA %d\n", nnode);
    write_vtk_point_scalar_const(
        fp, fom_solution, nnode, "FOM_solution");
    write_vtk_point_scalar_const(
        fp, projected_fom_solution, nnode, "POD_projection_of_FOM");
    write_vtk_point_scalar_const(
        fp, rom_solution, nnode, "ROM_solution");
    BB_vtk_write_point_vals_scalar(
        fp, signed_rom_error, nnode, "ROM_minus_FOM");
    BB_vtk_write_point_vals_scalar(
        fp, absolute_rom_error, nnode, "absolute_ROM_error");
    BB_vtk_write_point_vals_scalar(
        fp, relative_rom_error, nnode, "pointwise_relative_ROM_error");
    BB_vtk_write_point_vals_scalar(
        fp, signed_projection_error, nnode, "POD_projection_minus_FOM");
    BB_vtk_write_point_vals_scalar(
        fp, absolute_projection_error, nnode, "absolute_POD_projection_error");
    BB_vtk_write_point_vals_scalar(
        fp, dirichlet_mask, nnode, "Dirichlet_mask");
    BB_vtk_write_point_vals_scalar(
        fp, time_value, nnode, "physical_time");

    fclose(fp);

    free(signed_rom_error);
    free(absolute_rom_error);
    free(relative_rom_error);
    free(signed_projection_error);
    free(absolute_projection_error);
    free(dirichlet_mask);
    free(time_value);
}

static int apply_zero_temporal_coupling(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_STSB_CONTEXT* context =
        (ROM_STDD_STSB_CONTEXT*)user_context;

    (void)x;

    if(context == NULL || y == NULL || context->local_st_dof <= 0){
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(
        y,
        0,
        (size_t)context->local_st_dof * sizeof(double));

    return ROM_STDD_SUCCESS;
}

static int store_standard_snapshot(
    const FE_SYSTEM* sys,
    const ROM_STDD_SYSTEM* system,
    const double* current_full_state,
    const double* previous_full_state,
    ST_DDROM_FORMULATION formulation,
    double* internal_snapshot)
{
    if(sys == NULL || system == NULL || current_full_state == NULL ||
       internal_snapshot == NULL || system->st_dof_per_space_node != 1 ||
       (formulation == ST_DDROM_FORMULATION_INCREMENTAL &&
        previous_full_state == NULL))
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int node = 0; node < system->num_internal_space_nodes; node++){
        if(sys->bc.D_bc_exists[node]){
            internal_snapshot[node] = 0.0;
        }
        else if(formulation == ST_DDROM_FORMULATION_INCREMENTAL){
            internal_snapshot[node] =
                current_full_state[node] - previous_full_state[node];
        }
        else{
            internal_snapshot[node] = current_full_state[node];
        }
    }

    return ROM_STDD_SUCCESS;
}

static void prepare_standard_primitive_operator(
    ROM_STDD_STSB_CONTEXT* adapter)
{
    if(adapter == NULL || adapter->sys == NULL){
        fail_message("NULL standard-DDROM operator adapter");
    }

    /*
     * Verified HLPOD-DDROM ordering:
     *
     *   D_ij = V_i^T A_primitive V_j
     *
     * is formed before the strong Dirichlet transformation is applied to
     * the current FOM right-hand side.
     */
    monolis_clear_mat_value_rhs_R(&adapter->sys->monolis);
    monolis_copy_mat_value_R(
        &adapter->sys->monolis0,
        &adapter->sys->monolis);
}


static int build_standard_formulation_rhs(
    ROM_STDD_STSB_CONTEXT* adapter,
    ST_DDROM_FORMULATION formulation,
    int window_id,
    const double* previous_full_state,
    double* full_rhs,
    double* operator_work)
{
    if(adapter == NULL || previous_full_state == NULL ||
       full_rhs == NULL || operator_work == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /*
     * Current-state form:
     *
     *     A_tilde U^n = b_tilde(U^{n-1}, t_n).
     *
     * Increment form, derived from the identical strong pair:
     *
     *     A_tilde Delta U^n
     *       = b_tilde(U^{n-1}, t_n) - A_tilde U^{n-1}.
     *
     * At a Dirichlet row A_tilde is the identity, therefore the latter RHS
     * contains g_D^n-g_D^{n-1} automatically.
     */
    if(ROM_std_stdd_stsb_build_raw_rhs(
        adapter,
        window_id,
        previous_full_state,
        full_rhs) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    if(formulation == ST_DDROM_FORMULATION_INCREMENTAL){
        if(ROM_std_stdd_stsb_apply_prepared_A(
            previous_full_state,
            operator_work,
            adapter) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        for(int row = 0; row < adapter->local_st_dof; row++){
            full_rhs[row] -= operator_work[row];
        }
    }

    return ROM_STDD_SUCCESS;
}

static void run_fom_or_collect(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int num_windows,
    int num_internal_space_nodes,
    int standard_ddrom,
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
    double* previous_standard_full_state = NULL;
    double* physical_window_solution = NULL;
    ST_DDROM_FOM_SOLVE_RESIDUAL_SUM cumulative_residual;
    int file_number = 0;
    int residual_violation = 0;

    memset(&cumulative_residual, 0, sizeof(cumulative_residual));

    window_solution = (double*)calloc((size_t)local_st_dof, sizeof(double));
    previous_trace = (double*)calloc(
        (size_t)sys->fe.total_num_nodes,
        sizeof(double));
    lift = (double*)calloc((size_t)local_st_dof, sizeof(double));

    if(standard_ddrom){
        previous_standard_full_state = (double*)calloc(
            (size_t)local_st_dof,
            sizeof(double));
        physical_window_solution = (double*)calloc(
            (size_t)local_st_dof,
            sizeof(double));
    }

    if(options->mode == ST_DDROM_RUN_COLLECT){
        snapshots = (double*)calloc(
            (size_t)num_windows * (size_t)internal_st_dof,
            sizeof(double));
    }

    if(window_solution == NULL || previous_trace == NULL || lift == NULL ||
       (standard_ddrom &&
        (previous_standard_full_state == NULL ||
         physical_window_solution == NULL)) ||
       (options->mode == ST_DDROM_RUN_COLLECT && snapshots == NULL))
    {
        fail_message("FOM/collect allocation failed");
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys->fe,
        sys->vals.T,
        previous_trace);

    if(standard_ddrom){
        memcpy(
            previous_standard_full_state,
            sys->vals.T,
            (size_t)sys->fe.total_num_nodes * sizeof(double));

        monolis_mpi_update_R(
            &sys->monolis_com,
            sys->fe.total_num_nodes,
            1,
            previous_standard_full_state);
    }

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

        if(standard_ddrom){
            if(reference_fom_formulation() ==
               ST_DDROM_FORMULATION_INCREMENTAL)
            {
                for(int row = 0; row < local_st_dof; row++){
                    physical_window_solution[row] =
                        previous_standard_full_state[row]
                        + window_solution[row];
                }
            }
            else{
                memcpy(
                    physical_window_solution,
                    window_solution,
                    (size_t)local_st_dof * sizeof(double));
            }

            monolis_mpi_update_R(
                &sys->monolis_com,
                sys->fe.total_num_nodes,
                1,
                physical_window_solution);
        }

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
                standard_ddrom
                    ? physical_window_solution
                    : window_solution,
                &file_number);
        }

        /*
         * Snapshot extraction is read-only and is the only collect-specific
         * operation inserted into the reference FOM loop.
         */
        if(options->mode == ST_DDROM_RUN_COLLECT){
            double* snapshot = snapshots
                + (size_t)window * (size_t)internal_st_dof;

            if(standard_ddrom){
                /*
                 * The verified ordinary DDROM stores the current physical
                 * state U^n and zeroes only its Dirichlet rows.  It does not
                 * store U^n-U_D^n and it does not store a time increment.
                 */
                if(store_standard_snapshot(
                    sys,
                    rom_system,
                    physical_window_solution,
                    previous_standard_full_state,
                    options->formulation,
                    snapshot) != ROM_STDD_SUCCESS)
                {
                    fail_message(
                        "standard-DDROM formulation-matched snapshot extraction failed");
                }
            }
            else{
                if(ROM_std_stdd_stsb_build_dirichlet_lift(
                    adapter,
                    window,
                    lift) != ROM_STDD_SUCCESS ||
                   ROM_std_stdd_stsb_store_internal_homogeneous(
                    rom_system,
                    window_solution,
                    lift,
                    snapshot) != ROM_STDD_SUCCESS)
                {
                    fail_message("snapshot lifting failed");
                }
            }
        }

        /* Exact reference-FOM trace propagation. */
        STSB_extract_last_right_trace(
            time_element,
            num_window_slabs,
            sys->fe.total_num_nodes,
            standard_ddrom
                ? physical_window_solution
                : window_solution,
            previous_trace);

        {
            double solution_l2 = 0.0;
            double trace_l2 = 0.0;

            if(STSB_compute_window_checksum(
                sys,
                ndof,
                standard_ddrom
                    ? physical_window_solution
                    : window_solution,
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
        if(!standard_ddrom &&
           options->strict_fom_residual &&
           window_residual_violation)
        {
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

        if(standard_ddrom){
            memcpy(
                previous_standard_full_state,
                physical_window_solution,
                (size_t)local_st_dof * sizeof(double));
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
                standard_ddrom
                    ? "ignored in verified-HLPOD standard mode"
                    : (options->strict_fom_residual ? "yes" : "no"));
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

        write_formulation_manifest(
            options->snapshot_dir,
            options->formulation);
    }

    free(window_solution);
    free(previous_trace);
    free(lift);
    free(snapshots);
    free(previous_standard_full_state);
    free(physical_window_solution);
}

static void run_offline(
    ROM_STDD_SYSTEM* system,
    ROM_STDD_STSB_CONTEXT* adapter,
    int standard_ddrom,
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

    if(standard_ddrom){
        require_formulation_manifest(
            options->snapshot_dir,
            options->formulation);
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

    if(options->write_mode_vtk){
        output_pod_modes_vtk(
            adapter->sys,
            system,
            options);
    }

    if(standard_ddrom){
        prepare_standard_primitive_operator(adapter);
    }
    else if(ROM_std_stdd_stsb_prepare_operator(adapter, 0)
            != ROM_STDD_SUCCESS)
    {
        fail_message("A preparation failed");
    }

    if(ROM_std_stdd_project_operators_rank_sweep(
        system,
        ROM_std_stdd_stsb_apply_prepared_A,
        standard_ddrom
            ? apply_zero_temporal_coupling
            : ROM_std_stdd_stsb_apply_B,
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

    write_formulation_manifest(
        options->basis_dir,
        options->formulation);

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


static void run_online_standard_ddrom(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    ROM_STDD_SYSTEM* system,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ST_DDROM_OPTIONS* options)
{
    const int nx = sys->fe.total_num_nodes;
    const int local_dof = system->local_st_dof;
    const int internal_dof = system->internal_st_dof;
    const int capacity = system->num_modes_capacity;
    const size_t all_reduced_count =
        (size_t)system->num_windows * (size_t)capacity;

    double* previous_full_state = NULL;
    double* previous_reference_full_state = NULL;
    double* full_rhs = NULL;
    double* reference_full_rhs = NULL;
    double* operator_work = NULL;
    double* reference_operator_work = NULL;
    double* homogeneous = NULL;
    double* current_lift = NULL;
    double* full_state = NULL;
    double* reference_full_state = NULL;
    double* projected_reference_full_state = NULL;

    double* reference_snapshots = NULL;
    double* reference_coefficients = NULL;
    double* reference_reduced_rhs = NULL;
    double* projected_reference = NULL;

    double* rom_window_residuals = NULL;
    double* reference_window_defects = NULL;
    double* fom_window_residual_all = NULL;
    double* fom_window_residual_free = NULL;
    double* fom_window_residual_balanced_free = NULL;
    double* fom_window_abs_residual_free = NULL;
    double* fom_window_abs_ax_free = NULL;
    double* fom_window_abs_rhs_free = NULL;

    int file_number = 0;
    double projection_difference2 = 0.0;
    double projection_reference2 = 0.0;
    double state_difference2 = 0.0;
    double state_reference2 = 0.0;
    double coefficient_difference2 = 0.0;
    double coefficient_reference2 = 0.0;

    double fom_residual_all2 = 0.0;
    double fom_rhs_all2 = 0.0;
    double fom_residual_free2 = 0.0;
    double fom_rhs_free2 = 0.0;
    double fom_ax_free2 = 0.0;
    double fom_residual_dirichlet2 = 0.0;

    if(system->st_dof_per_space_node != 1 ||
       time_element->n_dof != 1 || num_window_slabs != 1)
    {
        fail_message(
            "HLPOD-reference standard DDROM requires dG(0), one slab/window");
    }

    require_formulation_manifest(
        options->basis_dir,
        options->formulation);

    if(options->validate){
        require_formulation_manifest(
            options->snapshot_dir,
            options->formulation);
    }

    if(ROM_std_stdd_read_basis_rank(system, options->basis_dir)
       != ROM_STDD_SUCCESS ||
       ROM_std_stdd_read_operators_rank(system, options->basis_dir)
       != ROM_STDD_SUCCESS)
    {
        fail_message("standard-DDROM basis/operator read failed");
    }

    {
        const size_t temporal_count =
            (size_t)(1 + system->num_neighbor_ranks)
            * (size_t)system->num_modes
            * (size_t)system->num_modes;
        double temporal_norm2 = 0.0;

        for(size_t i = 0; i < temporal_count; i++){
            temporal_norm2 +=
                system->temporal_blocks[i] * system->temporal_blocks[i];
        }

        if(!isfinite(temporal_norm2) || temporal_norm2 > 1.0e-28){
            fail_message(
                "standard-DDROM operator file contains temporal coupling; "
                "rebuild snapshots, basis, and operators");
        }
    }

    if(options->write_mode_vtk){
        output_pod_modes_vtk(
            sys,
            system,
            options);
    }

    previous_full_state = (double*)calloc((size_t)local_dof, sizeof(double));
    full_rhs = (double*)calloc((size_t)local_dof, sizeof(double));
    operator_work = (double*)calloc((size_t)local_dof, sizeof(double));
    homogeneous = (double*)calloc((size_t)local_dof, sizeof(double));
    current_lift = (double*)calloc((size_t)local_dof, sizeof(double));
    full_state = (double*)calloc((size_t)local_dof, sizeof(double));

    if(previous_full_state == NULL || full_rhs == NULL ||
       operator_work == NULL || homogeneous == NULL ||
       current_lift == NULL || full_state == NULL)
    {
        fail_message("standard-DDROM online allocation failed");
    }

    memcpy(
        previous_full_state,
        sys->vals.T,
        (size_t)nx * sizeof(double));

    monolis_mpi_update_R(
        &sys->monolis_com,
        nx,
        1,
        previous_full_state);

    if(options->validate){
        previous_reference_full_state = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        reference_full_rhs = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        reference_operator_work = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        reference_full_state = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        projected_reference_full_state = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        reference_snapshots = (double*)calloc(
            (size_t)system->num_windows * (size_t)internal_dof,
            sizeof(double));
        reference_coefficients = (double*)calloc(
            all_reduced_count,
            sizeof(double));
        reference_reduced_rhs = (double*)calloc(
            all_reduced_count,
            sizeof(double));
        projected_reference = (double*)calloc(
            (size_t)internal_dof,
            sizeof(double));

        rom_window_residuals = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        reference_window_defects = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_residual_all = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_residual_free = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_residual_balanced_free = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_abs_residual_free = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_abs_ax_free = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));
        fom_window_abs_rhs_free = (double*)calloc(
            (size_t)system->num_windows,
            sizeof(double));

        if(previous_reference_full_state == NULL ||
           reference_full_rhs == NULL ||
           reference_operator_work == NULL ||
           reference_full_state == NULL ||
           projected_reference_full_state == NULL ||
           reference_snapshots == NULL ||
           reference_coefficients == NULL ||
           reference_reduced_rhs == NULL ||
           projected_reference == NULL ||
           rom_window_residuals == NULL ||
           reference_window_defects == NULL ||
           fom_window_residual_all == NULL ||
           fom_window_residual_free == NULL ||
           fom_window_residual_balanced_free == NULL ||
           fom_window_abs_residual_free == NULL ||
           fom_window_abs_ax_free == NULL ||
           fom_window_abs_rhs_free == NULL)
        {
            fail_message("standard-DDROM validation allocation failed");
        }

        memcpy(
            previous_reference_full_state,
            sys->vals.T,
            (size_t)nx * sizeof(double));

        monolis_mpi_update_R(
            &sys->monolis_com,
            nx,
            1,
            previous_reference_full_state);

        if(ROM_std_stdd_read_snapshot_rank(
            system,
            reference_snapshots,
            options->validation_snapshot_id,
            options->snapshot_dir) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM validation snapshot read failed");
        }
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "standard DDROM formulation: %s "
            "(reference FOM contract: %s)\n",
            formulation_name(options->formulation),
            formulation_name(reference_fom_formulation()));
        printf(
            "standard DDROM operator/RHS path: primitive A projection, "
            "strong-BC pair, explicit halo exchange + mono_com0 matvec\n");

        if(options->formulation == ST_DDROM_FORMULATION_INCREMENTAL){
            printf(
                "standard DDROM update: solve Delta U^n from "
                "b_tilde-A_tilde U^{n-1}, then U^n=U^{n-1}+Delta U^n\n");
        }
        else{
            printf(
                "standard DDROM update: solve current U^n directly and "
                "replace the previous physical state\n");
        }

        if(options->validate){
            printf(
                "standard DDROM validation: ROM and FOM trajectories use "
                "separate previous states and formulation-matched RHS vectors\n");
        }
    }

    for(int window = 0; window < system->num_windows; window++){
        ROM_STDD_SYSTEM one_window_system = *system;
        double* reduced_rhs = system->reduced_rhs
            + (size_t)window * (size_t)capacity;
        double* reduced_coef = system->reduced_coef
            + (size_t)window * (size_t)capacity;

        if(build_standard_formulation_rhs(
            adapter,
            options->formulation,
            window,
            previous_full_state,
            full_rhs,
            operator_work) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_project_local_rhs(
            system,
            window,
            full_rhs,
            reduced_rhs) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM ROM RHS assembly/projection failed");
        }

        one_window_system.num_windows = 1;
        one_window_system.reduced_rhs = reduced_rhs;
        one_window_system.reduced_coef = reduced_coef;

        if(ROM_std_stdd_solve_window_marching(
            &one_window_system,
            &sys->monolis_com,
            options->reduced_max_iter,
            options->reduced_tolerance) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM distributed spatial solve failed");
        }

        if(ROM_std_stdd_reconstruct_local_homogeneous(
            system,
            window,
            homogeneous) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_stsb_build_dirichlet_lift(
            adapter,
            window,
            current_lift) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM state reconstruction failed");
        }

        for(int row = 0; row < local_dof; row++){
            const int node = row;

            if(options->formulation ==
               ST_DDROM_FORMULATION_INCREMENTAL)
            {
                const double boundary_increment =
                    sys->bc.D_bc_exists[node]
                        ? current_lift[row] - previous_full_state[row]
                        : 0.0;

                full_state[row] =
                    previous_full_state[row]
                    + homogeneous[row]
                    + boundary_increment;
            }
            else{
                full_state[row] =
                    homogeneous[row] + current_lift[row];
            }
        }

        monolis_mpi_update_R(
            &sys->monolis_com,
            nx,
            1,
            full_state);

        if(options->write_output){
            output_window(
                sys,
                time_element,
                num_window_slabs,
                window,
                full_state,
                &file_number);
        }

        if(options->validate){
            const double* reference = reference_snapshots
                + (size_t)window * (size_t)internal_dof;
            double* reference_coef = reference_coefficients
                + (size_t)window * (size_t)capacity;
            double* reference_rhs = reference_reduced_rhs
                + (size_t)window * (size_t)capacity;
            ST_DDROM_ABSOLUTE_RESIDUAL_SUM fom_residual;

            if(ROM_std_stdd_project_reference_coefficients(
                system,
                window,
                reference,
                reference_coef) != ROM_STDD_SUCCESS)
            {
                fail_message("standard-DDROM reference projection failed");
            }

            memset(
                projected_reference,
                0,
                (size_t)internal_dof * sizeof(double));

            for(int row = 0; row < internal_dof; row++){
                double value = 0.0;

                for(int mode = 0; mode < system->num_modes; mode++){
                    value += system->basis[
                        (size_t)row * (size_t)capacity + (size_t)mode]
                        * reference_coef[mode];
                }

                projected_reference[row] = value;

                {
                    const double projection_difference =
                        value - reference[row];

                    projection_difference2 +=
                        projection_difference * projection_difference;
                    projection_reference2 += reference[row] * reference[row];
                }
            }

            for(int row = 0; row < local_dof; row++){
                const int node = row;

                if(options->formulation ==
                   ST_DDROM_FORMULATION_INCREMENTAL)
                {
                    const double boundary_increment =
                        sys->bc.D_bc_exists[node]
                            ? current_lift[row]
                                - previous_reference_full_state[row]
                            : 0.0;
                    const double reference_increment =
                        (row < internal_dof) ? reference[row] : 0.0;
                    const double projected_increment =
                        (row < internal_dof)
                            ? projected_reference[row]
                            : 0.0;

                    reference_full_state[row] =
                        previous_reference_full_state[row]
                        + reference_increment
                        + boundary_increment;
                    projected_reference_full_state[row] =
                        previous_reference_full_state[row]
                        + projected_increment
                        + boundary_increment;
                }
                else{
                    const double reference_homogeneous =
                        (row < internal_dof) ? reference[row] : 0.0;
                    const double projected_homogeneous =
                        (row < internal_dof)
                            ? projected_reference[row]
                            : 0.0;

                    reference_full_state[row] =
                        current_lift[row] + reference_homogeneous;
                    projected_reference_full_state[row] =
                        current_lift[row] + projected_homogeneous;
                }
            }

            monolis_mpi_update_R(
                &sys->monolis_com,
                nx,
                1,
                reference_full_state);
            monolis_mpi_update_R(
                &sys->monolis_com,
                nx,
                1,
                projected_reference_full_state);

            if(options->write_comparison_vtk){
                const int global_step = window + 1;

                if(global_step % sys->vals.output_interval == 0){
                    output_standard_comparison_vtk(
                        sys,
                        global_step,
                        (double)global_step * sys->vals.dt,
                        reference_full_state,
                        projected_reference_full_state,
                        full_state);
                }
            }

            for(int row = 0; row < internal_dof; row++){
                const double state_difference =
                    full_state[row] - reference_full_state[row];

                state_difference2 += state_difference * state_difference;
                state_reference2 +=
                    reference_full_state[row] * reference_full_state[row];
            }

            for(int mode = 0; mode < system->num_modes; mode++){
                const double difference =
                    reduced_coef[mode] - reference_coef[mode];

                coefficient_difference2 += difference * difference;
                coefficient_reference2 +=
                    reference_coef[mode] * reference_coef[mode];
            }

            if(build_standard_formulation_rhs(
                adapter,
                options->formulation,
                window,
                previous_reference_full_state,
                reference_full_rhs,
                reference_operator_work) != ROM_STDD_SUCCESS ||
               ROM_std_stdd_project_local_rhs(
                system,
                window,
                reference_full_rhs,
                reference_rhs) != ROM_STDD_SUCCESS)
            {
                fail_message("standard-DDROM FOM RHS assembly/projection failed");
            }

            /*
             * The physical residual is always checked against the absolute
             * current-state strong pair.  In incremental ROM mode, rebuild
             * that pair after projecting the increment RHS.
             */
            if(options->formulation ==
               ST_DDROM_FORMULATION_INCREMENTAL)
            {
                if(ROM_std_stdd_stsb_build_raw_rhs(
                    adapter,
                    window,
                    previous_reference_full_state,
                    reference_full_rhs) != ROM_STDD_SUCCESS)
                {
                    fail_message(
                        "standard-DDROM absolute FOM pair rebuild failed");
                }
            }

            memset(&fom_residual, 0, sizeof(fom_residual));
            if(evaluate_current_strong_pair_residual(
                sys,
                1,
                reference_full_state,
                &fom_residual) != 0)
            {
                fail_message("standard-DDROM FOM full residual evaluation failed");
            }

            fom_window_residual_all[window] = safe_relative_norm(
                fom_residual.residual_all2,
                fom_residual.rhs_all2);
            fom_window_residual_free[window] = safe_relative_norm(
                fom_residual.residual_free2,
                fom_residual.rhs_free2);
            fom_window_residual_balanced_free[window] =
                balanced_relative_residual(
                    fom_residual.residual_free2,
                    fom_residual.ax_free2,
                    fom_residual.rhs_free2);
            fom_window_abs_residual_free[window] =
                sqrt(fmax(fom_residual.residual_free2, 0.0));
            fom_window_abs_ax_free[window] =
                sqrt(fmax(fom_residual.ax_free2, 0.0));
            fom_window_abs_rhs_free[window] =
                sqrt(fmax(fom_residual.rhs_free2, 0.0));

            fom_residual_all2 += fom_residual.residual_all2;
            fom_rhs_all2 += fom_residual.rhs_all2;
            fom_residual_free2 += fom_residual.residual_free2;
            fom_rhs_free2 += fom_residual.rhs_free2;
            fom_ax_free2 += fom_residual.ax_free2;
            fom_residual_dirichlet2 +=
                fom_residual.residual_dirichlet2;

            memcpy(
                previous_reference_full_state,
                reference_full_state,
                (size_t)local_dof * sizeof(double));
        }

        memcpy(
            previous_full_state,
            full_state,
            (size_t)local_dof * sizeof(double));
    }

    if(options->validate){
        double global_projection_difference2 = projection_difference2;
        double global_projection_reference2 = projection_reference2;
        double global_state_difference2 = state_difference2;
        double global_state_reference2 = state_reference2;
        double global_coefficient_difference2 = coefficient_difference2;
        double global_coefficient_reference2 = coefficient_reference2;
        double rom_reduced_residual;
        double reference_reduced_defect;
        double* saved_reduced_rhs;

        monolis_allreduce_R(
            1,
            &global_projection_difference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        monolis_allreduce_R(
            1,
            &global_projection_reference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        monolis_allreduce_R(
            1,
            &global_state_difference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        monolis_allreduce_R(
            1,
            &global_state_reference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        monolis_allreduce_R(
            1,
            &global_coefficient_difference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        monolis_allreduce_R(
            1,
            &global_coefficient_reference2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);

        rom_reduced_residual =
            ROM_std_stdd_system_reduced_residual_from_coefficients(
                system,
                &sys->monolis_com,
                system->reduced_coef,
                rom_window_residuals);

        saved_reduced_rhs = system->reduced_rhs;
        system->reduced_rhs = reference_reduced_rhs;

        reference_reduced_defect =
            ROM_std_stdd_system_reduced_residual_from_coefficients(
                system,
                &sys->monolis_com,
                reference_coefficients,
                reference_window_defects);

        system->reduced_rhs = saved_reduced_rhs;

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\nstandard DDROM per-window validation "
                "(%s formulation, separate ROM/FOM trajectories)\n",
                formulation_name(options->formulation));
            printf(
                "  window      ROM red.      FOM proj.     FOM all/rhs  "
                "FOM free/rhs FOM free/bal       ||r_f||      ||Ax_f||"
                "       ||b_f||\n");

            for(int window = 0; window < system->num_windows; window++){
                printf(
                    "  %6d  %12.5e  %12.5e  %12.5e  %12.5e  "
                    "%12.5e  %12.5e  %12.5e  %12.5e\n",
                    window,
                    rom_window_residuals[window],
                    reference_window_defects[window],
                    fom_window_residual_all[window],
                    fom_window_residual_free[window],
                    fom_window_residual_balanced_free[window],
                    fom_window_abs_residual_free[window],
                    fom_window_abs_ax_free[window],
                    fom_window_abs_rhs_free[window]);
            }

            printf(
                "\nstandard DDROM validation summary "
                "(formulation-matched HLPOD path)\n"
                "  formulation                     : %s\n"
                "  global projection error         : %.15e\n"
                "  global physical-state error     : %.15e\n"
                "  global coefficient error        : %.15e\n"
                "  ROM reduced residual            : %.15e\n"
                "  projected FOM reduced defect    : %.15e\n"
                "  FOM residual all / ||b||        : %.15e\n"
                "  FOM residual free / ||b_f||     : %.15e\n"
                "  FOM residual free balanced      : %.15e\n"
                "  FOM residual BC / ||b||         : %.15e\n"
                "  accumulated ||r_free||          : %.15e\n"
                "  accumulated ||Ax_free||         : %.15e\n"
                "  accumulated ||b_free||          : %.15e\n",
                formulation_name(options->formulation),
                safe_relative_norm(
                    global_projection_difference2,
                    global_projection_reference2),
                safe_relative_norm(
                    global_state_difference2,
                    global_state_reference2),
                safe_relative_norm(
                    global_coefficient_difference2,
                    global_coefficient_reference2),
                rom_reduced_residual,
                reference_reduced_defect,
                safe_relative_norm(fom_residual_all2, fom_rhs_all2),
                safe_relative_norm(fom_residual_free2, fom_rhs_free2),
                balanced_relative_residual(
                    fom_residual_free2,
                    fom_ax_free2,
                    fom_rhs_free2),
                sqrt(fmax(fom_residual_dirichlet2, 0.0)) /
                    fmax(sqrt(fmax(fom_rhs_all2, 0.0)), 1.0e-300),
                sqrt(fmax(fom_residual_free2, 0.0)),
                sqrt(fmax(fom_ax_free2, 0.0)),
                sqrt(fmax(fom_rhs_free2, 0.0)));
        }
    }
    else if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "standard DDROM reduced residual: %.15e\n",
            ROM_std_stdd_system_reduced_residual(
                system,
                &sys->monolis_com));
    }

    free(previous_full_state);
    free(previous_reference_full_state);
    free(full_rhs);
    free(reference_full_rhs);
    free(operator_work);
    free(reference_operator_work);
    free(homogeneous);
    free(current_lift);
    free(full_state);
    free(reference_full_state);
    free(projected_reference_full_state);
    free(reference_snapshots);
    free(reference_coefficients);
    free(reference_reduced_rhs);
    free(projected_reference);
    free(rom_window_residuals);
    free(reference_window_defects);
    free(fom_window_residual_all);
    free(fom_window_residual_free);
    free(fom_window_residual_balanced_free);
    free(fom_window_abs_residual_free);
    free(fom_window_abs_ax_free);
    free(fom_window_abs_rhs_free);
}

static double distributed_relative_difference(
    const double* first,
    const double* second,
    size_t count,
    const MONOLIS_COM* communication)
{
    double values[2] = {0.0, 0.0};

    if(first == NULL || second == NULL || communication == NULL){
        return INFINITY;
    }

    for(size_t i = 0; i < count; i++){
        const double difference = first[i] - second[i];
        values[0] += difference * difference;
        values[1] += second[i] * second[i];
    }

    monolis_allreduce_R(
        2,
        values,
        MONOLIS_MPI_SUM,
        communication->comm);

    return sqrt(fmax(values[0], 0.0)) /
        fmax(sqrt(fmax(values[1], 0.0)), 1.0e-300);
}

static void output_dg0_standard_comparison_vtk(
    FE_SYSTEM* sys,
    int global_step,
    double time,
    const double* fom_solution,
    const double* stddrom_solution,
    const double* standard_ddrom_solution,
    const double* hlpod_kernel_solution)
{
    const int nnode = (sys != NULL) ? sys->fe.total_num_nodes : 0;
    char local_name[4096];
    const char* filename;
    FILE* fp = NULL;
    double* st_minus_fom = NULL;
    double* dd_minus_fom = NULL;
    double* st_minus_dd = NULL;
    double* abs_st_minus_fom = NULL;
    double* abs_dd_minus_fom = NULL;
    double* abs_st_minus_dd = NULL;
    double* hlpod_minus_dd = NULL;
    double* abs_hlpod_minus_dd = NULL;
    double* dirichlet_mask = NULL;
    double* time_value = NULL;

    if(sys == NULL || fom_solution == NULL || stddrom_solution == NULL ||
       standard_ddrom_solution == NULL || nnode <= 0 || global_step <= 0)
    {
        fail_message("invalid dG0 comparison VTK argument");
    }

    st_minus_fom = (double*)calloc((size_t)nnode, sizeof(double));
    dd_minus_fom = (double*)calloc((size_t)nnode, sizeof(double));
    st_minus_dd = (double*)calloc((size_t)nnode, sizeof(double));
    abs_st_minus_fom = (double*)calloc((size_t)nnode, sizeof(double));
    abs_dd_minus_fom = (double*)calloc((size_t)nnode, sizeof(double));
    abs_st_minus_dd = (double*)calloc((size_t)nnode, sizeof(double));
    if(hlpod_kernel_solution != NULL){
        hlpod_minus_dd = (double*)calloc((size_t)nnode, sizeof(double));
        abs_hlpod_minus_dd = (double*)calloc((size_t)nnode, sizeof(double));
    }
    dirichlet_mask = (double*)calloc((size_t)nnode, sizeof(double));
    time_value = (double*)calloc((size_t)nnode, sizeof(double));

    if(st_minus_fom == NULL || dd_minus_fom == NULL ||
       st_minus_dd == NULL || abs_st_minus_fom == NULL ||
       abs_dd_minus_fom == NULL || abs_st_minus_dd == NULL ||
       (hlpod_kernel_solution != NULL &&
        (hlpod_minus_dd == NULL || abs_hlpod_minus_dd == NULL)) ||
       dirichlet_mask == NULL || time_value == NULL)
    {
        fail_message("dG0 comparison VTK allocation failed");
    }

    for(int node = 0; node < nnode; node++){
        st_minus_fom[node] =
            stddrom_solution[node] - fom_solution[node];
        dd_minus_fom[node] =
            standard_ddrom_solution[node] - fom_solution[node];
        st_minus_dd[node] =
            stddrom_solution[node] - standard_ddrom_solution[node];
        abs_st_minus_fom[node] = fabs(st_minus_fom[node]);
        abs_dd_minus_fom[node] = fabs(dd_minus_fom[node]);
        abs_st_minus_dd[node] = fabs(st_minus_dd[node]);
        if(hlpod_kernel_solution != NULL){
            hlpod_minus_dd[node] =
                hlpod_kernel_solution[node] - standard_ddrom_solution[node];
            abs_hlpod_minus_dd[node] = fabs(hlpod_minus_dd[node]);
        }
        dirichlet_mask[node] =
            sys->bc.D_bc_exists[node] ? 1.0 : 0.0;
        time_value[node] = time;
    }

    if(snprintf(
        local_name,
        sizeof(local_name),
        "dg0_rom_compare_%06d.vtk",
        global_step) >= (int)sizeof(local_name))
    {
        fail_message("dG0 comparison VTK filename is too long");
    }

    filename = monolis_get_global_output_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        "./",
        local_name);

    fp = ROM_BB_write_fopen(
        fp,
        filename,
        sys->cond.directory);

    if(fp == NULL){
        fail_message("dG0 comparison VTK open failed");
    }

    switch(sys->fe.local_num_nodes){
    case 4:
        BBFE_sys_write_vtk_shape(fp, &sys->fe, TYPE_VTK_TETRA);
        break;
    case 8:
        BBFE_sys_write_vtk_shape(fp, &sys->fe, TYPE_VTK_HEXAHEDRON);
        break;
    default:
        fclose(fp);
        fail_message("dG0 comparison VTK supports tetra/hexa meshes");
    }

    fprintf(fp, "POINT_DATA %d\n", nnode);
    write_vtk_point_scalar_const(
        fp, fom_solution, nnode, "FOM_solution");
    write_vtk_point_scalar_const(
        fp, stddrom_solution, nnode, "dG0_ST_DDROM_solution");
    write_vtk_point_scalar_const(
        fp, standard_ddrom_solution, nnode, "standard_DDROM_solution");
    if(hlpod_kernel_solution != NULL){
        write_vtk_point_scalar_const(
            fp, hlpod_kernel_solution, nnode, "HLPOD_kernel_solution");
    }
    BB_vtk_write_point_vals_scalar(
        fp, st_minus_fom, nnode, "dG0_ST_DDROM_minus_FOM");
    BB_vtk_write_point_vals_scalar(
        fp, dd_minus_fom, nnode, "standard_DDROM_minus_FOM");
    BB_vtk_write_point_vals_scalar(
        fp, st_minus_dd, nnode, "dG0_ST_DDROM_minus_standard_DDROM");
    BB_vtk_write_point_vals_scalar(
        fp, abs_st_minus_fom, nnode, "absolute_dG0_ST_DDROM_error");
    BB_vtk_write_point_vals_scalar(
        fp, abs_dd_minus_fom, nnode, "absolute_standard_DDROM_error");
    BB_vtk_write_point_vals_scalar(
        fp, abs_st_minus_dd, nnode, "absolute_ROM_difference");
    if(hlpod_kernel_solution != NULL){
        BB_vtk_write_point_vals_scalar(
            fp, hlpod_minus_dd, nnode,
            "HLPOD_kernel_minus_standard_DDROM");
        BB_vtk_write_point_vals_scalar(
            fp, abs_hlpod_minus_dd, nnode,
            "absolute_HLPOD_kernel_difference");
    }
    BB_vtk_write_point_vals_scalar(
        fp, dirichlet_mask, nnode, "Dirichlet_mask");
    BB_vtk_write_point_vals_scalar(
        fp, time_value, nnode, "physical_time");

    fclose(fp);
    free(st_minus_fom);
    free(dd_minus_fom);
    free(st_minus_dd);
    free(abs_st_minus_fom);
    free(abs_dd_minus_fom);
    free(abs_st_minus_dd);
    free(hlpod_minus_dd);
    free(abs_hlpod_minus_dd);
    free(dirichlet_mask);
    free(time_value);
}


typedef struct {
    HLPOD_VALUES values;
    HLPOD_MAT matrix;
    ROM_STDD_SYSTEM reduced_system;
    int initialized;
    int total_num_nodes;
    int num_internal_nodes;
    int num_modes;
    int n_neib_vec;
} ST_DDROM_HLPOD_KERNEL_CONTEXT;

static void hlpod_kernel_context_finalize(
    ST_DDROM_HLPOD_KERNEL_CONTEXT* context)
{
    if(context == NULL){
        return;
    }

    if(context->matrix.VTKV != NULL){
        BB_std_free_2d_double(
            context->matrix.VTKV,
            context->num_modes,
            context->n_neib_vec);
    }
    if(context->matrix.pod_modes != NULL){
        BB_std_free_2d_double(
            context->matrix.pod_modes,
            context->total_num_nodes,
            context->num_modes);
    }

    free(context->matrix.num_modes_internal);
    free(context->matrix.n_internal_vertex_subd);
    free(context->matrix.node_id);
    free(context->matrix.num_modes_2nddd);
    free(context->matrix.VTf);
    free(context->matrix.mode_coef);
    free(context->values.sol_vec);

    ROM_std_stdd_system_finalize(&context->reduced_system);
    memset(context, 0, sizeof(*context));
}

static void hlpod_kernel_context_initialize(
    ST_DDROM_HLPOD_KERNEL_CONTEXT* context,
    FE_SYSTEM* sys,
    const ROM_STDD_SYSTEM* standard_system)
{
    int max_modes;
    int negative_modes;
    int min_modes;
    const int num_modes = standard_system->num_modes;
    const int num_neighbors = sys->monolis_com.recv_n_neib;
    const int n_internal = standard_system->num_internal_space_nodes;
    const int nnode = sys->fe.total_num_nodes;

    if(context == NULL || sys == NULL || standard_system == NULL ||
       num_modes <= 0 || n_internal <= 0 || nnode <= 0)
    {
        fail_message("invalid HLPOD-kernel comparison context");
    }

    memset(context, 0, sizeof(*context));

    /*
     * The verified save-memory HLPOD kernel explicitly assumes that its
     * first rows are exactly the communication-owned physical vertices.
     * Do not silently remap here: a mismatch is itself an incompatibility
     * between the two implementations and must be reported.
     */
    if(n_internal != sys->monolis_com.n_internal_vertex){
        fprintf(
            stderr,
            "ERROR [rank %d]: HLPOD kernel ownership mismatch: "
            "basis rows=%d, communication-owned vertices=%d, local vertices=%d\n",
            monolis_mpi_get_global_my_rank(),
            n_internal,
            sys->monolis_com.n_internal_vertex,
            nnode);
        fprintf(
            stderr,
            "The basis/snapshot set was created with a non-owned node count. "
            "Use the ownership-fixed core_FOM_ST_spaceblock.c, delete the old "
            "snapshot/basis directories, and rerun collect + offline.\n");
        fail_message(
            "cannot execute ROM_std_hlpod_calc_reduced_mat_save_memory "
            "until the POD rows use owned vertices only");
    }

    max_modes = num_modes;
    negative_modes = -num_modes;
    monolis_allreduce_I(
        1,
        &max_modes,
        MONOLIS_MPI_MAX,
        sys->monolis_com.comm);
    monolis_allreduce_I(
        1,
        &negative_modes,
        MONOLIS_MPI_MAX,
        sys->monolis_com.comm);
    min_modes = -negative_modes;

    if(max_modes != min_modes){
        fail_message(
            "HLPOD-kernel comparison currently requires an equal mode count "
            "on every spatial rank");
    }

    context->total_num_nodes = nnode;
    context->num_internal_nodes = n_internal;
    context->num_modes = num_modes;
    context->n_neib_vec = num_modes * (1 + num_neighbors);

    context->values.num_modes = num_modes;
    context->values.num_modes_max = max_modes;
    context->values.n_neib_vec = context->n_neib_vec;
    context->values.sol_vec = (double*)calloc(
        (size_t)nnode,
        sizeof(double));

    context->matrix.pod_modes = BB_std_calloc_2d_double(
        context->matrix.pod_modes,
        nnode,
        num_modes);
    context->matrix.num_modes_internal = (int*)calloc(1, sizeof(int));
    context->matrix.n_internal_vertex_subd = (int*)calloc(1, sizeof(int));
    context->matrix.node_id = (int*)calloc((size_t)n_internal, sizeof(int));
    context->matrix.num_modes_2nddd = (int*)calloc(
        (size_t)(1 + num_neighbors),
        sizeof(int));
    context->matrix.VTf = (double*)calloc(
        (size_t)context->n_neib_vec,
        sizeof(double));
    context->matrix.mode_coef = (double*)calloc(
        (size_t)context->n_neib_vec,
        sizeof(double));

    if(context->values.sol_vec == NULL ||
       context->matrix.pod_modes == NULL ||
       context->matrix.num_modes_internal == NULL ||
       context->matrix.n_internal_vertex_subd == NULL ||
       context->matrix.node_id == NULL ||
       context->matrix.num_modes_2nddd == NULL ||
       context->matrix.VTf == NULL ||
       context->matrix.mode_coef == NULL)
    {
        hlpod_kernel_context_finalize(context);
        fail_message("HLPOD-kernel comparison allocation failed");
    }

    context->matrix.num_modes_internal[0] = num_modes;
    context->matrix.n_internal_vertex_subd[0] = n_internal;

    for(int i = 0; i < n_internal; i++){
        context->matrix.node_id[i] = i;
        for(int mode = 0; mode < num_modes; mode++){
            context->matrix.pod_modes[i][mode] =
                standard_system->basis[
                    (size_t)i *
                    (size_t)standard_system->num_modes_capacity +
                    (size_t)mode];
        }
    }

    for(int slot = 0; slot < 1 + num_neighbors; slot++){
        context->matrix.num_modes_2nddd[slot] = num_modes;
    }

    if(ROM_std_stdd_system_initialize(
        &context->reduced_system,
        &sys->monolis_com,
        standard_system->num_windows,
        standard_system->num_space_nodes,
        standard_system->num_internal_space_nodes,
        1,
        standard_system->num_modes_capacity) != ROM_STDD_SUCCESS)
    {
        hlpod_kernel_context_finalize(context);
        fail_message("HLPOD-kernel reduced-system initialization failed");
    }

    context->reduced_system.num_modes = num_modes;
    context->reduced_system.diagonal_blocks = (double*)calloc(
        (size_t)(1 + num_neighbors) *
            (size_t)num_modes * (size_t)num_modes,
        sizeof(double));
    context->reduced_system.temporal_blocks = (double*)calloc(
        (size_t)(1 + num_neighbors) *
            (size_t)num_modes * (size_t)num_modes,
        sizeof(double));

    if(context->reduced_system.diagonal_blocks == NULL ||
       context->reduced_system.temporal_blocks == NULL)
    {
        hlpod_kernel_context_finalize(context);
        fail_message("HLPOD-kernel reduced-block allocation failed");
    }

    context->initialized = 1;
}

static double hlpod_kernel_build_and_compare_matrix(
    ST_DDROM_HLPOD_KERNEL_CONTEXT* context,
    FE_SYSTEM* sys,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ROM_STDD_SYSTEM* standard_system)
{
    const int num_modes = context->num_modes;
    const int block_slots = 1 + standard_system->num_neighbor_ranks;
    const size_t block_size =
        (size_t)num_modes * (size_t)num_modes;

    prepare_standard_primitive_operator(adapter);

    ROM_std_hlpod_calc_reduced_mat_save_memory(
        &sys->monolis,
        &sys->monolis_com,
        &sys->mono_com0,
        &context->values,
        &context->matrix,
        context->total_num_nodes,
        context->n_neib_vec,
        1,
        num_modes,
        1);

    if(context->matrix.VTKV == NULL){
        fail_message("ROM_std_hlpod_calc_reduced_mat_save_memory returned NULL VTKV");
    }

    for(int slot = 0; slot < block_slots; slot++){
        const int hlpod_column_offset = slot * num_modes;
        double* block = context->reduced_system.diagonal_blocks
            + (size_t)slot * block_size;

        for(int row_mode = 0; row_mode < num_modes; row_mode++){
            for(int column_mode = 0; column_mode < num_modes; column_mode++){
                block[
                    (size_t)row_mode * (size_t)num_modes +
                    (size_t)column_mode] =
                    context->matrix.VTKV[row_mode]
                        [hlpod_column_offset + column_mode];
            }
        }
    }

    return distributed_relative_difference(
        standard_system->diagonal_blocks,
        context->reduced_system.diagonal_blocks,
        (size_t)block_slots * block_size,
        &sys->monolis_com);
}

static double hlpod_kernel_project_rhs(
    ST_DDROM_HLPOD_KERNEL_CONTEXT* context,
    FE_SYSTEM* sys,
    const double* standard_rhs,
    int window)
{
    double* stored_rhs = context->reduced_system.reduced_rhs
        + (size_t)window *
            (size_t)context->reduced_system.num_modes_capacity;

    ROM_std_hlpod_calc_reduced_rhs(
        &sys->monolis,
        &context->matrix,
        context->num_modes,
        1,
        1);

    memcpy(
        stored_rhs,
        context->matrix.VTf,
        (size_t)context->num_modes * sizeof(double));

    return distributed_relative_difference(
        standard_rhs,
        context->matrix.VTf,
        (size_t)context->num_modes,
        &sys->monolis_com);
}

static void hlpod_kernel_solve_and_reconstruct(
    ST_DDROM_HLPOD_KERNEL_CONTEXT* context,
    FE_SYSTEM* sys,
    const ST_DDROM_OPTIONS* options,
    int window,
    double* homogeneous)
{
    ROM_STDD_SYSTEM one_window_system = context->reduced_system;
    double* reduced_rhs = context->reduced_system.reduced_rhs
        + (size_t)window *
            (size_t)context->reduced_system.num_modes_capacity;
    double* reduced_coef = context->reduced_system.reduced_coef
        + (size_t)window *
            (size_t)context->reduced_system.num_modes_capacity;

    one_window_system.num_windows = 1;
    one_window_system.reduced_rhs = reduced_rhs;
    one_window_system.reduced_coef = reduced_coef;

    if(ROM_std_stdd_solve_window_marching(
        &one_window_system,
        &sys->monolis_com,
        options->reduced_max_iter,
        options->reduced_tolerance) != ROM_STDD_SUCCESS)
    {
        fail_message("HLPOD-kernel generated reduced solve failed");
    }

    memset(
        context->matrix.mode_coef,
        0,
        (size_t)context->n_neib_vec * sizeof(double));
    memcpy(
        context->matrix.mode_coef,
        reduced_coef,
        (size_t)context->num_modes * sizeof(double));

    ROM_std_hlpod_calc_sol(
        &context->values,
        &context->matrix,
        context->total_num_nodes,
        context->num_modes,
        1,
        1);

    memcpy(
        homogeneous,
        context->values.sol_vec,
        (size_t)context->total_num_nodes * sizeof(double));
    monolis_mpi_update_R(
        &sys->monolis_com,
        context->total_num_nodes,
        1,
        homogeneous);
}

static void run_dg0_standard_comparison(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    ROM_STDD_SYSTEM* st_system,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ST_DDROM_OPTIONS* options)
{
    ROM_STDD_SYSTEM dd_system;
    ST_DDROM_HLPOD_KERNEL_CONTEXT hlpod_context;
    const int nx = sys->fe.total_num_nodes;
    const int internal_dof = st_system->internal_st_dof;
    const int local_dof = st_system->local_st_dof;
    const int capacity = st_system->num_modes_capacity;
    const int block_slots = 1 + st_system->num_neighbor_ranks;
    size_t block_count;

    double* initial_trace = NULL;
    double* st_effective_rhs = NULL;
    double* identity_full_rhs = NULL;
    double* dd_solve_full_rhs = NULL;
    double* operator_work = NULL;
    double* previous_st_full = NULL;
    double* previous_dd_full = NULL;
    double* st_homogeneous = NULL;
    double* dd_homogeneous = NULL;
    double* st_full = NULL;
    double* dd_full = NULL;
    double* current_lift = NULL;
    double* dd_identity_reduced = NULL;
    double* st_total_reduced = NULL;
    double* defect = NULL;
    double* self_action = NULL;
    double* neighbor_action = NULL;
    double* temporal_action = NULL;
    double* st_window_residual = NULL;
    double* dd_window_residual = NULL;
    double* reference_snapshots = NULL;
    double* reference_homogeneous = NULL;
    double* reference_full = NULL;
    double* previous_hlpod_full = NULL;
    double* hlpod_homogeneous = NULL;
    double* hlpod_full = NULL;
    double* hlpod_window_residual = NULL;

    FILE* csv = NULL;
    FILE* hlpod_csv = NULL;
    double basis_difference;
    double diagonal_difference;
    double zero_temporal_norm;
    double max_rhs_identity = 0.0;
    double max_coefficient_difference = 0.0;
    double max_homogeneous_difference = 0.0;
    double max_physical_difference = 0.0;
    double max_trace_difference = 0.0;
    double max_st_fom_error = 0.0;
    double max_dd_fom_error = 0.0;
    double st_reduced_residual;
    double dd_reduced_residual;
    double hlpod_reduced_residual = 0.0;
    double hlpod_matrix_difference = 0.0;
    double max_hlpod_rhs_difference = 0.0;
    double max_hlpod_coefficient_difference = 0.0;
    double max_hlpod_homogeneous_difference = 0.0;
    double max_hlpod_physical_difference = 0.0;

    memset(&dd_system, 0, sizeof(dd_system));
    memset(&hlpod_context, 0, sizeof(hlpod_context));

    if(time_element == NULL || time_element->degree != 0 ||
       time_element->n_dof != 1 || num_window_slabs != 1 ||
       st_system->st_dof_per_space_node != 1)
    {
        fail_message(
            "dG0/standard-DDROM comparison requires dG(0), one slab/window");
    }

    if(options->formulation != ST_DDROM_FORMULATION_CURRENT_STATE){
        fail_message(
            "dG0 comparison currently requires current_state formulation");
    }

    require_formulation_manifest(
        options->basis_dir,
        ST_DDROM_FORMULATION_CURRENT_STATE);
    if(options->validate){
        require_formulation_manifest(
            options->snapshot_dir,
            ST_DDROM_FORMULATION_CURRENT_STATE);
    }

    if(ROM_std_stdd_system_initialize(
        &dd_system,
        &sys->monolis_com,
        st_system->num_windows,
        st_system->num_space_nodes,
        st_system->num_internal_space_nodes,
        1,
        st_system->num_modes_capacity) != ROM_STDD_SUCCESS)
    {
        fail_message("standard comparison system initialization failed");
    }

    if(ROM_std_stdd_read_basis_rank(
        st_system,
        options->basis_dir) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_read_basis_rank(
        &dd_system,
        options->basis_dir) != ROM_STDD_SUCCESS)
    {
        fail_message("dG0 comparison basis read failed");
    }

    if(st_system->num_modes != dd_system.num_modes){
        fail_message("comparison systems have different mode counts");
    }

    basis_difference = distributed_relative_difference(
        st_system->basis,
        dd_system.basis,
        (size_t)internal_dof * (size_t)capacity,
        &sys->monolis_com);

    if(ROM_std_stdd_stsb_prepare_operator(adapter, 0)
       != ROM_STDD_SUCCESS ||
       ROM_std_stdd_project_operators_rank_sweep(
        st_system,
        ROM_std_stdd_stsb_apply_prepared_A,
        ROM_std_stdd_stsb_apply_B,
        adapter) != ROM_STDD_SUCCESS)
    {
        fail_message("dG0 ST-DDROM operator projection failed");
    }

    prepare_standard_primitive_operator(adapter);
    if(ROM_std_stdd_project_operators_rank_sweep(
        &dd_system,
        ROM_std_stdd_stsb_apply_prepared_A,
        apply_zero_temporal_coupling,
        adapter) != ROM_STDD_SUCCESS)
    {
        fail_message("standard-DDROM operator projection failed");
    }

    block_count = (size_t)block_slots
        * (size_t)st_system->num_modes
        * (size_t)st_system->num_modes;

    diagonal_difference = distributed_relative_difference(
        st_system->diagonal_blocks,
        dd_system.diagonal_blocks,
        block_count,
        &sys->monolis_com);

    {
        double zero = 0.0;
        double local_norm2 = 0.0;
        for(size_t i = 0; i < block_count; i++){
            local_norm2 +=
                dd_system.temporal_blocks[i]
                * dd_system.temporal_blocks[i];
        }
        monolis_allreduce_R(
            1,
            &local_norm2,
            MONOLIS_MPI_SUM,
            sys->monolis_com.comm);
        zero_temporal_norm = sqrt(fmax(local_norm2, zero));
    }

    if(options->compare_hlpod_kernels){
        hlpod_kernel_context_initialize(
            &hlpod_context,
            sys,
            &dd_system);
        hlpod_matrix_difference = hlpod_kernel_build_and_compare_matrix(
            &hlpod_context,
            sys,
            adapter,
            &dd_system);
    }

    initial_trace = (double*)calloc((size_t)nx, sizeof(double));
    st_effective_rhs = (double*)calloc((size_t)local_dof, sizeof(double));
    identity_full_rhs = (double*)calloc((size_t)local_dof, sizeof(double));
    dd_solve_full_rhs = (double*)calloc((size_t)local_dof, sizeof(double));
    operator_work = (double*)calloc((size_t)local_dof, sizeof(double));
    previous_st_full = (double*)calloc((size_t)local_dof, sizeof(double));
    previous_dd_full = (double*)calloc((size_t)local_dof, sizeof(double));
    st_homogeneous = (double*)calloc((size_t)local_dof, sizeof(double));
    dd_homogeneous = (double*)calloc((size_t)local_dof, sizeof(double));
    st_full = (double*)calloc((size_t)local_dof, sizeof(double));
    dd_full = (double*)calloc((size_t)local_dof, sizeof(double));
    current_lift = (double*)calloc((size_t)local_dof, sizeof(double));
    dd_identity_reduced = (double*)calloc((size_t)capacity, sizeof(double));
    st_total_reduced = (double*)calloc((size_t)capacity, sizeof(double));
    defect = (double*)calloc((size_t)capacity, sizeof(double));
    self_action = (double*)calloc((size_t)capacity, sizeof(double));
    neighbor_action = (double*)calloc((size_t)capacity, sizeof(double));
    temporal_action = (double*)calloc((size_t)capacity, sizeof(double));
    st_window_residual = (double*)calloc(
        (size_t)st_system->num_windows,
        sizeof(double));
    dd_window_residual = (double*)calloc(
        (size_t)st_system->num_windows,
        sizeof(double));

    if(options->compare_hlpod_kernels){
        previous_hlpod_full = (double*)calloc(
            (size_t)local_dof, sizeof(double));
        hlpod_homogeneous = (double*)calloc(
            (size_t)local_dof, sizeof(double));
        hlpod_full = (double*)calloc(
            (size_t)local_dof, sizeof(double));
        hlpod_window_residual = (double*)calloc(
            (size_t)st_system->num_windows, sizeof(double));
    }

    if(options->validate){
        reference_snapshots = (double*)calloc(
            (size_t)st_system->num_windows * (size_t)internal_dof,
            sizeof(double));
        reference_homogeneous = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
        reference_full = (double*)calloc(
            (size_t)local_dof,
            sizeof(double));
    }

    if(initial_trace == NULL || st_effective_rhs == NULL ||
       identity_full_rhs == NULL || dd_solve_full_rhs == NULL ||
       operator_work == NULL || previous_st_full == NULL ||
       previous_dd_full == NULL || st_homogeneous == NULL ||
       dd_homogeneous == NULL || st_full == NULL || dd_full == NULL ||
       current_lift == NULL || dd_identity_reduced == NULL ||
       st_total_reduced == NULL || defect == NULL ||
       self_action == NULL || neighbor_action == NULL ||
       temporal_action == NULL || st_window_residual == NULL ||
       dd_window_residual == NULL ||
       (options->compare_hlpod_kernels &&
        (previous_hlpod_full == NULL || hlpod_homogeneous == NULL ||
         hlpod_full == NULL || hlpod_window_residual == NULL)) ||
       (options->validate &&
        (reference_snapshots == NULL || reference_homogeneous == NULL ||
         reference_full == NULL)))
    {
        fail_message("dG0 comparison allocation failed");
    }

    if(options->validate &&
       ROM_std_stdd_read_snapshot_rank(
        st_system,
        reference_snapshots,
        options->validation_snapshot_id,
        options->snapshot_dir) != ROM_STDD_SUCCESS)
    {
        fail_message("dG0 comparison reference snapshot read failed");
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys->fe,
        sys->vals.T,
        initial_trace);

    for(int window = 0; window < st_system->num_windows; window++){
        double* reduced_rhs = st_system->reduced_rhs
            + (size_t)window * (size_t)capacity;

        if(ROM_std_stdd_stsb_build_effective_rhs(
            adapter,
            window,
            initial_trace,
            st_effective_rhs) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_project_local_rhs(
            st_system,
            window,
            st_effective_rhs,
            reduced_rhs) != ROM_STDD_SUCCESS)
        {
            fail_message("dG0 ST-DDROM RHS construction failed");
        }
    }

    if(ROM_std_stdd_solve_window_marching(
        st_system,
        &sys->monolis_com,
        options->reduced_max_iter,
        options->reduced_tolerance) != ROM_STDD_SUCCESS)
    {
        fail_message("dG0 ST-DDROM solve failed");
    }

    memcpy(previous_st_full, sys->vals.T, (size_t)nx * sizeof(double));
    memcpy(previous_dd_full, sys->vals.T, (size_t)nx * sizeof(double));
    monolis_mpi_update_R(&sys->monolis_com, nx, 1, previous_st_full);
    monolis_mpi_update_R(&sys->monolis_com, nx, 1, previous_dd_full);

    if(options->compare_hlpod_kernels){
        memcpy(previous_hlpod_full, sys->vals.T, (size_t)nx * sizeof(double));
        monolis_mpi_update_R(
            &sys->monolis_com, nx, 1, previous_hlpod_full);
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        csv = fopen(options->compare_csv, "w");
        if(csv == NULL){
            fprintf(
                stderr,
                "ERROR: cannot open %s: %s\n",
                options->compare_csv,
                strerror(errno));
            fail_message("dG0 comparison CSV open failed");
        }
        fprintf(
            csv,
            "window,time,rhs_identity_error,coefficient_error,"
            "homogeneous_state_error,physical_state_error,trace_error,"
            "stddrom_vs_fom,standard_ddrom_vs_fom\n");

        if(options->compare_hlpod_kernels){
            hlpod_csv = fopen(options->hlpod_compare_csv, "w");
            if(hlpod_csv == NULL){
                fprintf(
                    stderr,
                    "ERROR: cannot open %s: %s\n",
                    options->hlpod_compare_csv,
                    strerror(errno));
                fail_message("HLPOD comparison CSV open failed");
            }
            fprintf(
                hlpod_csv,
                "window,time,VTf_rhs_difference,coefficient_difference,"
                "calc_sol_homogeneous_difference,physical_difference\n");
        }
    }

    for(int window = 0; window < st_system->num_windows; window++){
        ROM_STDD_SYSTEM one_window_system = dd_system;
        double* dd_reduced_rhs = dd_system.reduced_rhs
            + (size_t)window * (size_t)capacity;
        double* dd_reduced_coef = dd_system.reduced_coef
            + (size_t)window * (size_t)capacity;
        const double* st_reduced_coef = st_system->reduced_coef
            + (size_t)window * (size_t)capacity;
        double rhs_identity_error;
        double coefficient_error;
        double homogeneous_error;
        double physical_error;
        double trace_error;
        double st_fom_error = 0.0;
        double dd_fom_error = 0.0;
        double hlpod_rhs_difference = 0.0;
        double hlpod_coefficient_difference = 0.0;
        double hlpod_homogeneous_difference = 0.0;
        double hlpod_physical_difference = 0.0;

        if(ROM_std_stdd_stsb_reconstruct_full_window(
            adapter,
            st_system,
            window,
            st_homogeneous,
            st_full) != ROM_STDD_SUCCESS)
        {
            fail_message("dG0 ST-DDROM reconstruction failed");
        }

        memset(dd_identity_reduced, 0, (size_t)capacity * sizeof(double));
        if(build_standard_formulation_rhs(
            adapter,
            ST_DDROM_FORMULATION_CURRENT_STATE,
            window,
            previous_st_full,
            identity_full_rhs,
            operator_work) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_project_local_rhs(
            &dd_system,
            window,
            identity_full_rhs,
            dd_identity_reduced) != ROM_STDD_SUCCESS)
        {
            fail_message("dG0 RHS identity construction failed");
        }

        memset(defect, 0, (size_t)capacity * sizeof(double));
        memset(self_action, 0, (size_t)capacity * sizeof(double));
        memset(neighbor_action, 0, (size_t)capacity * sizeof(double));
        memset(temporal_action, 0, (size_t)capacity * sizeof(double));

        if(ROM_std_stdd_evaluate_local_reduced_equation(
            st_system,
            &sys->monolis_com,
            st_system->reduced_coef,
            window,
            defect,
            self_action,
            neighbor_action,
            temporal_action) != ROM_STDD_SUCCESS)
        {
            fail_message("dG0 temporal reduced action evaluation failed");
        }

        for(int mode = 0; mode < st_system->num_modes; mode++){
            st_total_reduced[mode] =
                st_system->reduced_rhs[
                    (size_t)window * (size_t)capacity + (size_t)mode]
                + temporal_action[mode];
        }

        rhs_identity_error = distributed_relative_difference(
            st_total_reduced,
            dd_identity_reduced,
            (size_t)st_system->num_modes,
            &sys->monolis_com);

        if(build_standard_formulation_rhs(
            adapter,
            ST_DDROM_FORMULATION_CURRENT_STATE,
            window,
            previous_dd_full,
            dd_solve_full_rhs,
            operator_work) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_project_local_rhs(
            &dd_system,
            window,
            dd_solve_full_rhs,
            dd_reduced_rhs) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM solve RHS construction failed");
        }

        one_window_system.num_windows = 1;
        one_window_system.reduced_rhs = dd_reduced_rhs;
        one_window_system.reduced_coef = dd_reduced_coef;

        if(ROM_std_stdd_solve_window_marching(
            &one_window_system,
            &sys->monolis_com,
            options->reduced_max_iter,
            options->reduced_tolerance) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM comparison solve failed");
        }

        if(ROM_std_stdd_reconstruct_local_homogeneous(
            &dd_system,
            window,
            dd_homogeneous) != ROM_STDD_SUCCESS ||
           ROM_std_stdd_stsb_build_dirichlet_lift(
            adapter,
            window,
            current_lift) != ROM_STDD_SUCCESS)
        {
            fail_message("standard-DDROM comparison reconstruction failed");
        }

        for(int row = 0; row < local_dof; row++){
            dd_full[row] = dd_homogeneous[row] + current_lift[row];
        }
        monolis_mpi_update_R(&sys->monolis_com, nx, 1, dd_full);

        if(options->compare_hlpod_kernels){
            /*
             * Rebuild the same current-state strong pair on the independent
             * HLPOD trajectory.  This makes the HLPOD path a real online
             * calculation rather than merely reusing the standard-DDROM RHS.
             */
            if(build_standard_formulation_rhs(
                adapter,
                ST_DDROM_FORMULATION_CURRENT_STATE,
                window,
                previous_hlpod_full,
                identity_full_rhs,
                operator_work) != ROM_STDD_SUCCESS)
            {
                fail_message("HLPOD-kernel online RHS construction failed");
            }

            hlpod_rhs_difference = hlpod_kernel_project_rhs(
                &hlpod_context,
                sys,
                dd_reduced_rhs,
                window);

            hlpod_kernel_solve_and_reconstruct(
                &hlpod_context,
                sys,
                options,
                window,
                hlpod_homogeneous);

            for(int row = 0; row < local_dof; row++){
                hlpod_full[row] =
                    hlpod_homogeneous[row] + current_lift[row];
            }
            monolis_mpi_update_R(
                &sys->monolis_com, nx, 1, hlpod_full);

            hlpod_coefficient_difference = distributed_relative_difference(
                dd_reduced_coef,
                hlpod_context.reduced_system.reduced_coef
                    + (size_t)window * (size_t)capacity,
                (size_t)st_system->num_modes,
                &sys->monolis_com);
            hlpod_homogeneous_difference = distributed_relative_difference(
                dd_homogeneous,
                hlpod_homogeneous,
                (size_t)internal_dof,
                &sys->monolis_com);
            hlpod_physical_difference = distributed_relative_difference(
                dd_full,
                hlpod_full,
                (size_t)internal_dof,
                &sys->monolis_com);
        }

        coefficient_error = distributed_relative_difference(
            st_reduced_coef,
            dd_reduced_coef,
            (size_t)st_system->num_modes,
            &sys->monolis_com);
        homogeneous_error = distributed_relative_difference(
            st_homogeneous,
            dd_homogeneous,
            (size_t)internal_dof,
            &sys->monolis_com);
        physical_error = distributed_relative_difference(
            st_full,
            dd_full,
            (size_t)internal_dof,
            &sys->monolis_com);
        trace_error = physical_error;

        if(options->validate){
            memset(
                reference_homogeneous,
                0,
                (size_t)local_dof * sizeof(double));
            memcpy(
                reference_homogeneous,
                reference_snapshots
                    + (size_t)window * (size_t)internal_dof,
                (size_t)internal_dof * sizeof(double));
            monolis_mpi_update_R(
                &sys->monolis_com,
                nx,
                1,
                reference_homogeneous);

            for(int row = 0; row < local_dof; row++){
                reference_full[row] =
                    reference_homogeneous[row] + current_lift[row];
            }

            st_fom_error = distributed_relative_difference(
                st_full,
                reference_full,
                (size_t)internal_dof,
                &sys->monolis_com);
            dd_fom_error = distributed_relative_difference(
                dd_full,
                reference_full,
                (size_t)internal_dof,
                &sys->monolis_com);

            if(options->write_comparison_vtk &&
               ((window + 1) % sys->vals.output_interval == 0))
            {
                output_dg0_standard_comparison_vtk(
                    sys,
                    window + 1,
                    (double)(window + 1) * sys->vals.dt,
                    reference_full,
                    st_full,
                    dd_full,
                    options->compare_hlpod_kernels ? hlpod_full : NULL);
            }
        }

        max_rhs_identity = fmax(max_rhs_identity, rhs_identity_error);
        max_coefficient_difference = fmax(
            max_coefficient_difference,
            coefficient_error);
        max_homogeneous_difference = fmax(
            max_homogeneous_difference,
            homogeneous_error);
        max_physical_difference = fmax(
            max_physical_difference,
            physical_error);
        max_trace_difference = fmax(max_trace_difference, trace_error);
        max_st_fom_error = fmax(max_st_fom_error, st_fom_error);
        max_dd_fom_error = fmax(max_dd_fom_error, dd_fom_error);

        if(options->compare_hlpod_kernels){
            max_hlpod_rhs_difference = fmax(
                max_hlpod_rhs_difference,
                hlpod_rhs_difference);
            max_hlpod_coefficient_difference = fmax(
                max_hlpod_coefficient_difference,
                hlpod_coefficient_difference);
            max_hlpod_homogeneous_difference = fmax(
                max_hlpod_homogeneous_difference,
                hlpod_homogeneous_difference);
            max_hlpod_physical_difference = fmax(
                max_hlpod_physical_difference,
                hlpod_physical_difference);
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(
                csv,
                "%d,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,"
                "%.15e,%.15e\n",
                window,
                (double)(window + 1) * sys->vals.dt,
                rhs_identity_error,
                coefficient_error,
                homogeneous_error,
                physical_error,
                trace_error,
                st_fom_error,
                dd_fom_error);

            if(options->compare_hlpod_kernels){
                fprintf(
                    hlpod_csv,
                    "%d,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                    window,
                    (double)(window + 1) * sys->vals.dt,
                    hlpod_rhs_difference,
                    hlpod_coefficient_difference,
                    hlpod_homogeneous_difference,
                    hlpod_physical_difference);
            }
        }

        memcpy(previous_st_full, st_full, (size_t)local_dof * sizeof(double));
        memcpy(previous_dd_full, dd_full, (size_t)local_dof * sizeof(double));
        if(options->compare_hlpod_kernels){
            memcpy(
                previous_hlpod_full,
                hlpod_full,
                (size_t)local_dof * sizeof(double));
        }
    }

    st_reduced_residual =
        ROM_std_stdd_system_reduced_residual_from_coefficients(
            st_system,
            &sys->monolis_com,
            st_system->reduced_coef,
            st_window_residual);
    dd_reduced_residual =
        ROM_std_stdd_system_reduced_residual_from_coefficients(
            &dd_system,
            &sys->monolis_com,
            dd_system.reduced_coef,
            dd_window_residual);

    if(options->compare_hlpod_kernels){
        hlpod_reduced_residual =
            ROM_std_stdd_system_reduced_residual_from_coefficients(
                &hlpod_context.reduced_system,
                &sys->monolis_com,
                hlpod_context.reduced_system.reduced_coef,
                hlpod_window_residual);
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        fclose(csv);
        csv = fopen(options->compare_csv, "a");
        if(csv != NULL){
            fprintf(csv, "# per-window reduced residuals\n");
            fprintf(csv, "# window,stddrom_reduced_residual,standard_ddrom_reduced_residual\n");
            for(int window = 0; window < st_system->num_windows; window++){
                fprintf(
                    csv,
                    "# %d,%.15e,%.15e\n",
                    window,
                    st_window_residual[window],
                    dd_window_residual[window]);
            }
            fprintf(
                csv,
                "# summary,max_rhs_identity=%.15e,max_coefficient=%.15e,"
                "max_homogeneous=%.15e,max_physical=%.15e,"
                "max_trace=%.15e,max_st_fom=%.15e,max_dd_fom=%.15e,"
                "st_residual=%.15e,dd_residual=%.15e\n",
                max_rhs_identity,
                max_coefficient_difference,
                max_homogeneous_difference,
                max_physical_difference,
                max_trace_difference,
                max_st_fom_error,
                max_dd_fom_error,
                st_reduced_residual,
                dd_reduced_residual);
            fclose(csv);
            csv = NULL;
        }

        if(options->compare_hlpod_kernels){
            if(hlpod_csv != NULL){
                fprintf(
                    hlpod_csv,
                    "# summary,VTKV_matrix_difference=%.15e,"
                    "max_VTf_rhs_difference=%.15e,"
                    "max_coefficient_difference=%.15e,"
                    "max_calc_sol_homogeneous_difference=%.15e,"
                    "max_physical_difference=%.15e,"
                    "reduced_residual=%.15e\n",
                    hlpod_matrix_difference,
                    max_hlpod_rhs_difference,
                    max_hlpod_coefficient_difference,
                    max_hlpod_homogeneous_difference,
                    max_hlpod_physical_difference,
                    hlpod_reduced_residual);
                fclose(hlpod_csv);
                hlpod_csv = NULL;
            }

            printf(
                "\nVerified HLPOD kernel comparison (actual function calls)\n"
                "  ROM_std_hlpod_calc_reduced_mat_save_memory VTKV diff : %.15e\n"
                "  ROM_std_hlpod_calc_reduced_rhs max VTf diff          : %.15e\n"
                "  generated-system coefficient difference              : %.15e\n"
                "  ROM_std_hlpod_calc_sol homogeneous difference        : %.15e\n"
                "  HLPOD-kernel physical-state difference               : %.15e\n"
                "  HLPOD-kernel reduced residual                         : %.15e\n"
                "  comparison tolerance                                  : %.15e\n"
                "  equivalence status                                     : %s\n"
                "  CSV                                                    : %s\n",
                hlpod_matrix_difference,
                max_hlpod_rhs_difference,
                max_hlpod_coefficient_difference,
                max_hlpod_homogeneous_difference,
                max_hlpod_physical_difference,
                hlpod_reduced_residual,
                options->hlpod_compare_tolerance,
                (hlpod_matrix_difference <= options->hlpod_compare_tolerance &&
                 max_hlpod_rhs_difference <= options->hlpod_compare_tolerance &&
                 max_hlpod_coefficient_difference <= options->hlpod_compare_tolerance &&
                 max_hlpod_physical_difference <= options->hlpod_compare_tolerance)
                    ? "MATCH"
                    : "MISMATCH",
                options->hlpod_compare_csv);
        }

        printf(
            "\ndG(0) ST-DDROM vs standard spatial DDROM validation\n"
            "  exact-comparison conditions       : dG(0), one slab/window, "
            "same local basis, current_state\n"
            "  basis relative difference         : %.15e\n"
            "  reduced D relative difference     : %.15e\n"
            "  standard temporal block norm      : %.15e\n"
            "  max reduced RHS identity error    : %.15e\n"
            "  max coefficient difference        : %.15e\n"
            "  max homogeneous-state difference  : %.15e\n"
            "  max physical-state difference     : %.15e\n"
            "  max right-trace difference        : %.15e\n"
            "  ST-DDROM reduced residual         : %.15e\n"
            "  standard DDROM reduced residual   : %.15e\n"
            "  max ST-DDROM/FOM state error      : %.15e\n"
            "  max standard-DDROM/FOM error      : %.15e\n"
            "  comparison tolerance              : %.15e\n"
            "  equivalence status                : %s\n"
            "  CSV                               : %s\n",
            basis_difference,
            diagonal_difference,
            zero_temporal_norm,
            max_rhs_identity,
            max_coefficient_difference,
            max_homogeneous_difference,
            max_physical_difference,
            max_trace_difference,
            st_reduced_residual,
            dd_reduced_residual,
            max_st_fom_error,
            max_dd_fom_error,
            options->compare_tolerance,
            (basis_difference <= options->compare_tolerance &&
             diagonal_difference <= options->compare_tolerance &&
             max_rhs_identity <= options->compare_tolerance &&
             max_coefficient_difference <= options->compare_tolerance &&
             max_physical_difference <= options->compare_tolerance)
                ? "MATCH"
                : "MISMATCH",
            options->compare_csv);
    }

    free(initial_trace);
    free(st_effective_rhs);
    free(identity_full_rhs);
    free(dd_solve_full_rhs);
    free(operator_work);
    free(previous_st_full);
    free(previous_dd_full);
    free(st_homogeneous);
    free(dd_homogeneous);
    free(st_full);
    free(dd_full);
    free(current_lift);
    free(dd_identity_reduced);
    free(st_total_reduced);
    free(defect);
    free(self_action);
    free(neighbor_action);
    free(temporal_action);
    free(st_window_residual);
    free(dd_window_residual);
    free(reference_snapshots);
    free(reference_homogeneous);
    free(reference_full);
    free(previous_hlpod_full);
    free(hlpod_homogeneous);
    free(hlpod_full);
    free(hlpod_window_residual);
    hlpod_kernel_context_finalize(&hlpod_context);
    ROM_std_stdd_system_finalize(&dd_system);
}


typedef struct {
    FE_SYSTEM* sys;
    ROM_STDD_STSB_CONTEXT* adapter;
    const ROM_STDD_SYSTEM* stdd;
    int num_owned_nodes;
    int num_local_nodes;
    int st_dof_per_node;
    int owned_st_dof;
    int local_st_dof;
    double* local_input;
    double* local_output;
} ST_DDROM_STROM_OPERATOR_CONTEXT;

static double** strom_compare_alloc_rows(int rows, int columns)
{
    double** matrix = NULL;
    double* data = NULL;

    if(rows <= 0 || columns <= 0){
        return NULL;
    }

    matrix = (double**)malloc((size_t)rows * sizeof(double*));
    data = (double*)calloc(
        (size_t)rows * (size_t)columns,
        sizeof(double));

    if(matrix == NULL || data == NULL){
        free(matrix);
        free(data);
        return NULL;
    }

    for(int row = 0; row < rows; row++){
        matrix[row] = data + (size_t)row * (size_t)columns;
    }

    return matrix;
}

static void strom_compare_owned_serial_to_local_node_major(
    const ST_DDROM_STROM_OPERATOR_CONTEXT* context,
    const double* serial_owned,
    double* local_node_major)
{
    memset(
        local_node_major,
        0,
        (size_t)context->local_st_dof * sizeof(double));

    for(int temporal_dof = 0;
        temporal_dof < context->st_dof_per_node;
        temporal_dof++)
    {
        for(int node = 0; node < context->num_owned_nodes; node++){
            local_node_major[
                (size_t)node * (size_t)context->st_dof_per_node
                + (size_t)temporal_dof] =
                serial_owned[
                    (size_t)temporal_dof
                        * (size_t)context->num_owned_nodes
                    + (size_t)node];
        }
    }
}

static void strom_compare_local_node_major_to_owned_serial(
    const ST_DDROM_STROM_OPERATOR_CONTEXT* context,
    const double* local_node_major,
    double* serial_owned)
{
    for(int temporal_dof = 0;
        temporal_dof < context->st_dof_per_node;
        temporal_dof++)
    {
        for(int node = 0; node < context->num_owned_nodes; node++){
            serial_owned[
                (size_t)temporal_dof
                    * (size_t)context->num_owned_nodes
                + (size_t)node] =
                local_node_major[
                    (size_t)node * (size_t)context->st_dof_per_node
                    + (size_t)temporal_dof];
        }
    }
}

static int strom_compare_prepare_A(
    int window_id,
    void* user_context)
{
    ST_DDROM_STROM_OPERATOR_CONTEXT* context =
        (ST_DDROM_STROM_OPERATOR_CONTEXT*)user_context;

    if(context == NULL || context->adapter == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_std_stdd_stsb_prepare_operator(
        context->adapter,
        window_id) == ROM_STDD_SUCCESS
        ? ROM_STROM_SUCCESS
        : ROM_STROM_ERR_OPERATOR;
}

static int strom_compare_apply_prepared_A(
    const double* x,
    double* y,
    void* user_context)
{
    ST_DDROM_STROM_OPERATOR_CONTEXT* context =
        (ST_DDROM_STROM_OPERATOR_CONTEXT*)user_context;

    if(context == NULL || x == NULL || y == NULL ||
       context->local_input == NULL || context->local_output == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_compare_owned_serial_to_local_node_major(
        context,
        x,
        context->local_input);

    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->num_local_nodes,
        context->st_dof_per_node,
        context->local_input);

    memset(
        context->local_output,
        0,
        (size_t)context->local_st_dof * sizeof(double));

    if(ROM_std_stdd_stsb_apply_prepared_A(
        context->local_input,
        context->local_output,
        context->adapter) != ROM_STDD_SUCCESS)
    {
        return ROM_STROM_ERR_OPERATOR;
    }

    strom_compare_local_node_major_to_owned_serial(
        context,
        context->local_output,
        y);

    return ROM_STROM_SUCCESS;
}

static int strom_compare_apply_A(
    int window_id,
    const double* x,
    double* y,
    void* user_context)
{
    if(strom_compare_prepare_A(window_id, user_context)
       != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_OPERATOR;
    }

    return strom_compare_apply_prepared_A(x, y, user_context);
}

static int strom_compare_apply_B(
    int current_window_id,
    const double* x_previous,
    double* y_current,
    void* user_context)
{
    ST_DDROM_STROM_OPERATOR_CONTEXT* context =
        (ST_DDROM_STROM_OPERATOR_CONTEXT*)user_context;

    (void)current_window_id;

    if(context == NULL || x_previous == NULL || y_current == NULL ||
       context->local_input == NULL || context->local_output == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_compare_owned_serial_to_local_node_major(
        context,
        x_previous,
        context->local_input);

    memset(
        context->local_output,
        0,
        (size_t)context->local_st_dof * sizeof(double));

    /* ROM_std_stdd_stsb_apply_B performs the required halo update itself. */
    if(ROM_std_stdd_stsb_apply_B(
        context->local_input,
        context->local_output,
        context->adapter) != ROM_STDD_SUCCESS)
    {
        return ROM_STROM_ERR_OPERATOR;
    }

    strom_compare_local_node_major_to_owned_serial(
        context,
        context->local_output,
        y_current);

    return ROM_STROM_SUCCESS;
}

static double strom_compare_relative_global_arrays(
    const double* lhs,
    const double* rhs,
    size_t count)
{
    double difference2 = 0.0;
    double reference2 = 0.0;

    for(size_t i = 0; i < count; i++){
        const double difference = lhs[i] - rhs[i];
        difference2 += difference * difference;
        reference2 += rhs[i] * rhs[i];
    }

    if(reference2 > 0.0){
        return sqrt(difference2 / reference2);
    }
    return sqrt(difference2);
}

static double strom_compare_relative_distributed_arrays(
    const double* lhs,
    const double* rhs,
    size_t count,
    int communicator)
{
    double values[2] = {0.0, 0.0};

    for(size_t i = 0; i < count; i++){
        const double difference = lhs[i] - rhs[i];
        values[0] += difference * difference;
        values[1] += rhs[i] * rhs[i];
    }

    monolis_allreduce_R(
        2,
        values,
        MONOLIS_MPI_SUM,
        communicator);

    if(values[1] > 0.0){
        return sqrt(values[0] / values[1]);
    }
    return sqrt(values[0]);
}

static void strom_compare_collect_stdd_global_matrix(
    const ROM_STDD_SYSTEM* system,
    const double* distributed_blocks,
    double* global_matrix,
    int global_modes,
    int local_modes,
    int my_rank,
    int communicator)
{
    const size_t local_block_size =
        (size_t)local_modes * (size_t)local_modes;
    const int local_offset = my_rank * local_modes;

    memset(
        global_matrix,
        0,
        (size_t)global_modes * (size_t)global_modes * sizeof(double));

    for(int slot = 0; slot < 1 + system->num_neighbor_ranks; slot++){
        const int column_rank = slot == 0
            ? my_rank
            : system->neighbor_ranks[slot - 1];
        const int column_offset = column_rank * local_modes;
        const double* block = distributed_blocks
            + (size_t)slot * local_block_size;

        for(int i = 0; i < local_modes; i++){
            for(int j = 0; j < local_modes; j++){
                global_matrix[
                    (size_t)(local_offset + i) * (size_t)global_modes
                    + (size_t)(column_offset + j)] =
                    block[
                        (size_t)i * (size_t)local_modes
                        + (size_t)j];
            }
        }
    }

    monolis_allreduce_R(
        global_modes * global_modes,
        global_matrix,
        MONOLIS_MPI_SUM,
        communicator);
}

static void strom_compare_collect_stdd_global_vector(
    const double* local_values,
    int local_modes,
    int global_modes,
    int my_rank,
    double* global_values,
    int communicator)
{
    memset(
        global_values,
        0,
        (size_t)global_modes * sizeof(double));

    memcpy(
        global_values + (size_t)my_rank * (size_t)local_modes,
        local_values,
        (size_t)local_modes * sizeof(double));

    monolis_allreduce_R(
        global_modes,
        global_values,
        MONOLIS_MPI_SUM,
        communicator);
}

static void strom_compare_write_matrix_csv(
    const char* summary_path,
    const char* suffix,
    const double* stdd_matrix,
    const double* strom_matrix,
    int matrix_size)
{
    char filename[4096];
    FILE* fp;

    if(monolis_mpi_get_global_my_rank() != 0){
        return;
    }

    if(snprintf(
        filename,
        sizeof(filename),
        "%s.%s.csv",
        summary_path,
        suffix) >= (int)sizeof(filename))
    {
        fail_message("ST-ROM matrix comparison path is too long");
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        fail_message("cannot open ST-ROM matrix comparison CSV");
    }

    fprintf(fp, "row,column,stddrom,strom,difference\n");
    for(int row = 0; row < matrix_size; row++){
        for(int column = 0; column < matrix_size; column++){
            const size_t index =
                (size_t)row * (size_t)matrix_size + (size_t)column;
            fprintf(
                fp,
                "%d,%d,%.17e,%.17e,%.17e\n",
                row,
                column,
                stdd_matrix[index],
                strom_matrix[index],
                stdd_matrix[index] - strom_matrix[index]);
        }
    }
    fclose(fp);
}

static void run_sequential_strom_comparison(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    ROM_STDD_SYSTEM* stdd_system,
    ROM_STDD_STSB_CONTEXT* adapter,
    const ST_DDROM_OPTIONS* options)
{
    ROM_STROM_SYSTEM strom_system;
    HLPOD_VALUES* strom_values = NULL;
    HLPOD_MAT* strom_matrices = NULL;
    ST_DDROM_STROM_OPERATOR_CONTEXT operator_context;
    const int my_rank = monolis_mpi_get_global_my_rank();
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int communicator = monolis_mpi_get_global_comm();
    const int local_modes = stdd_system->num_modes;
    const int global_modes = comm_size * local_modes;
    const int num_windows = stdd_system->num_windows;
    const int num_owned = stdd_system->num_internal_space_nodes;
    const int num_local = stdd_system->num_space_nodes;
    const int st_dof_per_node = stdd_system->st_dof_per_space_node;
    const int owned_st_dof = stdd_system->internal_st_dof;
    const int local_st_dof = stdd_system->local_st_dof;
    const size_t global_matrix_count =
        (size_t)global_modes * (size_t)global_modes;

    double* stdd_D = NULL;
    double* stdd_L = NULL;
    double* strom_D = NULL;
    double* strom_L = NULL;
    double* stdd_g = NULL;
    double* strom_g = NULL;
    double* stdd_a = NULL;
    double* raw_local_rhs = NULL;
    double* raw_owned_rhs = NULL;
    double* effective_owned_rhs = NULL;
    double* stdd_effective_local = NULL;
    double* stdd_effective_owned = NULL;
    double* lift_local = NULL;
    double* initial_trace = NULL;
    double* zero_trace = NULL;
    double* stdd_q_local = NULL;
    double* stdd_q_owned = NULL;
    double* strom_q_owned = NULL;
    double* stdd_trace_local = NULL;
    double* stdd_q_right = NULL;
    double* strom_q_right = NULL;
    double* effective_rhs_window_difference = NULL;
    FILE* csv = NULL;

    double D_difference = INFINITY;
    double L_difference = INFINITY;
    double max_effective_rhs_difference = 0.0;
    double max_g_difference = 0.0;
    double max_a_difference = 0.0;
    double max_q_difference = 0.0;
    double max_q_right_difference = 0.0;
    int local_mode_max = local_modes;
    int negative_local_mode_min = -local_modes;
    int global_mode_min;

    memset(&strom_system, 0, sizeof(strom_system));
    memset(&operator_context, 0, sizeof(operator_context));

    if(options == NULL || !options->compare_strom_reference){
        return;
    }

    if(options->formulation != ST_DDROM_FORMULATION_CURRENT_STATE){
        fail_message(
            "sequential ST-ROM comparison currently requires current_state");
    }

    monolis_allreduce_I(
        1,
        &local_mode_max,
        MONOLIS_MPI_MAX,
        communicator);
    monolis_allreduce_I(
        1,
        &negative_local_mode_min,
        MONOLIS_MPI_MAX,
        communicator);
    global_mode_min = -negative_local_mode_min;

    if(local_mode_max != global_mode_min ||
       local_mode_max != local_modes)
    {
        fail_message(
            "sequential ST-ROM comparison requires one common mode count per rank");
    }

    if((long long)global_modes * (long long)global_modes > INT_MAX){
        fail_message("sequential ST-ROM comparison matrix is too large");
    }

    strom_values = (HLPOD_VALUES*)calloc(
        (size_t)num_windows,
        sizeof(HLPOD_VALUES));
    strom_matrices = (HLPOD_MAT*)calloc(
        (size_t)num_windows,
        sizeof(HLPOD_MAT));

    if(strom_values == NULL || strom_matrices == NULL ||
       ROM_std_strom_system_initialize(
           &strom_system,
           num_windows) != ROM_STROM_SUCCESS)
    {
        fail_message("sequential ST-ROM reference initialization failed");
    }

    for(int window = 0; window < num_windows; window++){
        if(ROM_std_strom_window_initialize(
            &strom_system.windows[window],
            window,
            &strom_values[window],
            &strom_matrices[window],
            num_owned,
            1,
            num_window_slabs,
            time_element->n_dof,
            1,
            global_modes,
            1) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM window initialization failed");
        }
    }

    strom_system.windows[0].hlpod_mat->pod_modes =
        strom_compare_alloc_rows(owned_st_dof, global_modes);
    if(strom_system.windows[0].hlpod_mat->pod_modes == NULL){
        fail_message("sequential ST-ROM basis allocation failed");
    }

    strom_system.windows[0].num_modes = global_modes;
    strom_system.windows[0].hlpod_vals->num_modes = global_modes;

    for(int temporal_dof = 0;
        temporal_dof < st_dof_per_node;
        temporal_dof++)
    {
        for(int node = 0; node < num_owned; node++){
            const int serial_row = temporal_dof * num_owned + node;
            const int stdd_row = node * st_dof_per_node + temporal_dof;

            for(int mode = 0; mode < local_modes; mode++){
                strom_system.windows[0].hlpod_mat->pod_modes
                    [serial_row][my_rank * local_modes + mode] =
                    stdd_system->basis[
                        (size_t)stdd_row
                            * (size_t)stdd_system->num_modes_capacity
                        + (size_t)mode];
            }
        }
    }

    for(int window = 1; window < num_windows; window++){
        if(ROM_std_strom_attach_pod_basis(
            &strom_system.windows[window],
            &strom_system.windows[0]) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM shared-basis attachment failed");
        }
    }

    if(ROM_std_strom_update_mode_offsets(&strom_system)
       != ROM_STROM_SUCCESS)
    {
        fail_message("sequential ST-ROM mode-offset construction failed");
    }

    operator_context.sys = sys;
    operator_context.adapter = adapter;
    operator_context.stdd = stdd_system;
    operator_context.num_owned_nodes = num_owned;
    operator_context.num_local_nodes = num_local;
    operator_context.st_dof_per_node = st_dof_per_node;
    operator_context.owned_st_dof = owned_st_dof;
    operator_context.local_st_dof = local_st_dof;
    operator_context.local_input = (double*)calloc(
        (size_t)local_st_dof,
        sizeof(double));
    operator_context.local_output = (double*)calloc(
        (size_t)local_st_dof,
        sizeof(double));

    stdd_D = (double*)calloc(global_matrix_count, sizeof(double));
    stdd_L = (double*)calloc(global_matrix_count, sizeof(double));
    strom_D = (double*)calloc(global_matrix_count, sizeof(double));
    strom_L = (double*)calloc(global_matrix_count, sizeof(double));
    stdd_g = (double*)calloc((size_t)global_modes, sizeof(double));
    strom_g = (double*)calloc((size_t)global_modes, sizeof(double));
    stdd_a = (double*)calloc((size_t)global_modes, sizeof(double));
    raw_local_rhs = (double*)calloc((size_t)local_st_dof, sizeof(double));
    raw_owned_rhs = (double*)calloc((size_t)owned_st_dof, sizeof(double));
    effective_owned_rhs = (double*)calloc(
        (size_t)owned_st_dof,
        sizeof(double));
    stdd_effective_local = (double*)calloc(
        (size_t)local_st_dof,
        sizeof(double));
    stdd_effective_owned = (double*)calloc(
        (size_t)owned_st_dof,
        sizeof(double));
    lift_local = (double*)calloc((size_t)local_st_dof, sizeof(double));
    initial_trace = (double*)calloc((size_t)num_local, sizeof(double));
    zero_trace = (double*)calloc((size_t)num_local, sizeof(double));
    stdd_q_local = (double*)calloc((size_t)local_st_dof, sizeof(double));
    stdd_q_owned = (double*)calloc((size_t)owned_st_dof, sizeof(double));
    strom_q_owned = (double*)calloc((size_t)owned_st_dof, sizeof(double));
    stdd_trace_local = (double*)calloc((size_t)num_local, sizeof(double));
    stdd_q_right = (double*)calloc((size_t)num_owned, sizeof(double));
    strom_q_right = (double*)calloc((size_t)num_owned, sizeof(double));
    effective_rhs_window_difference = (double*)calloc(
        (size_t)num_windows, sizeof(double));

    if(operator_context.local_input == NULL ||
       operator_context.local_output == NULL ||
       stdd_D == NULL || stdd_L == NULL || strom_D == NULL ||
       strom_L == NULL || stdd_g == NULL || strom_g == NULL ||
       stdd_a == NULL || raw_local_rhs == NULL ||
       raw_owned_rhs == NULL || effective_owned_rhs == NULL ||
       stdd_effective_local == NULL || stdd_effective_owned == NULL ||
       lift_local == NULL || initial_trace == NULL || zero_trace == NULL ||
       stdd_q_local == NULL || stdd_q_owned == NULL ||
       strom_q_owned == NULL || stdd_trace_local == NULL ||
       stdd_q_right == NULL || strom_q_right == NULL ||
       effective_rhs_window_difference == NULL)
    {
        fail_message("sequential ST-ROM comparison allocation failed");
    }

    STSB_set_previous_trace_from_nodal_value(
        &sys->fe,
        sys->vals.T,
        initial_trace);

    /* Fill the serial ST-ROM reference lift in time-major owned-row order. */
    for(int window = 0; window < num_windows; window++){
        if(ROM_std_stdd_stsb_build_dirichlet_lift(
            adapter,
            window,
            lift_local) != ROM_STDD_SUCCESS)
        {
            fail_message("sequential ST-ROM lift construction failed");
        }

        strom_compare_local_node_major_to_owned_serial(
            &operator_context,
            lift_local,
            strom_system.windows[window].reference_state);
    }

    /* D = Phi^T A Phi through the validated sequential ST-ROM kernel. */
    if(ROM_std_strom_calc_reduced_diag_prepared(
        &strom_system.windows[0],
        strom_compare_prepare_A,
        strom_compare_apply_prepared_A,
        &operator_context) != ROM_STROM_SUCCESS)
    {
        fail_message("sequential ST-ROM D construction failed");
    }

    memcpy(
        strom_D,
        strom_system.windows[0].reduced_diag[0],
        global_matrix_count * sizeof(double));
    monolis_allreduce_R(
        global_modes * global_modes,
        strom_D,
        MONOLIS_MPI_SUM,
        communicator);

    /* L = Phi^T B Phi.  All practical staged cases contain >=2 windows. */
    if(num_windows < 2 ||
       ROM_std_strom_calc_reduced_coupling(
           &strom_system.windows[1],
           &strom_system.windows[0],
           strom_compare_apply_B,
           &operator_context) != ROM_STROM_SUCCESS)
    {
        fail_message("sequential ST-ROM L construction failed");
    }

    memcpy(
        strom_L,
        strom_system.windows[1].reduced_coupling_prev[0],
        global_matrix_count * sizeof(double));
    monolis_allreduce_R(
        global_modes * global_modes,
        strom_L,
        MONOLIS_MPI_SUM,
        communicator);

    strom_compare_collect_stdd_global_matrix(
        stdd_system,
        stdd_system->diagonal_blocks,
        stdd_D,
        global_modes,
        local_modes,
        my_rank,
        communicator);
    strom_compare_collect_stdd_global_matrix(
        stdd_system,
        stdd_system->temporal_blocks,
        stdd_L,
        global_modes,
        local_modes,
        my_rank,
        communicator);

    D_difference = strom_compare_relative_global_arrays(
        stdd_D,
        strom_D,
        global_matrix_count);
    L_difference = strom_compare_relative_global_arrays(
        stdd_L,
        strom_L,
        global_matrix_count);

    for(int window = 0; window < num_windows; window++){
        if(ROM_std_strom_set_reduced_diag_data(
            &strom_system.windows[window],
            strom_D) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM D installation failed");
        }

        if(window > 0 &&
           ROM_std_strom_set_reduced_coupling_data(
               &strom_system.windows[window],
               strom_L,
               global_modes) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM L installation failed");
        }
    }

    make_parent_directory_collective(options->strom_compare_csv);

    if(my_rank == 0){
        csv = fopen(options->strom_compare_csv, "w");
        if(csv == NULL){
            fprintf(
                stderr,
                "ERROR [rank 0]: fopen(%s): %s\n",
                options->strom_compare_csv,
                strerror(errno));
            fflush(stderr);
            fail_message("cannot open sequential ST-ROM comparison CSV");
        }
        fprintf(
            csv,
            "window,time,D_relative,L_relative,effective_rhs_relative,"
            "g_relative,a_relative,q_relative,q_right_relative\n");
    }

    /*
     * Build g independently using the sequential ST-ROM affine-lift identity:
     *   g = Phi^T (f_raw - A U_D + B U_D,prev).
     */
    for(int window = 0; window < num_windows; window++){
        const double* previous_trace = window == 0
            ? initial_trace
            : zero_trace;
        ROM_STROM_WINDOW* current = &strom_system.windows[window];
        ROM_STROM_WINDOW* previous = window > 0
            ? &strom_system.windows[window - 1]
            : NULL;

        if(ROM_std_stdd_stsb_build_raw_rhs(
            adapter,
            window,
            previous_trace,
            raw_local_rhs) != ROM_STDD_SUCCESS)
        {
            fail_message("sequential ST-ROM raw RHS construction failed");
        }

        strom_compare_local_node_major_to_owned_serial(
            &operator_context,
            raw_local_rhs,
            raw_owned_rhs);

        if(ROM_std_strom_build_effective_rhs(
            current,
            previous,
            raw_owned_rhs,
            effective_owned_rhs,
            strom_compare_apply_A,
            &operator_context,
            strom_compare_apply_B,
            &operator_context) != ROM_STROM_SUCCESS ||
           ROM_std_strom_project_rhs(
            current,
            effective_owned_rhs) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM g construction failed");
        }

        memcpy(
            strom_g,
            current->reduced_rhs,
            (size_t)global_modes * sizeof(double));
        monolis_allreduce_R(
            global_modes,
            strom_g,
            MONOLIS_MPI_SUM,
            communicator);
        memcpy(
            current->reduced_rhs,
            strom_g,
            (size_t)global_modes * sizeof(double));

        strom_compare_collect_stdd_global_vector(
            stdd_system->reduced_rhs
                + (size_t)window
                    * (size_t)stdd_system->num_modes_capacity,
            local_modes,
            global_modes,
            my_rank,
            stdd_g,
            communicator);

        max_g_difference = fmax(
            max_g_difference,
            strom_compare_relative_global_arrays(
                stdd_g,
                strom_g,
                (size_t)global_modes));

        /* Also verify the full owned effective RHS before projection. */
        if(ROM_std_stdd_stsb_build_effective_rhs(
            adapter,
            window,
            initial_trace,
            stdd_effective_local) != ROM_STDD_SUCCESS)
        {
            fail_message("ST-DDROM effective RHS reconstruction failed");
        }
        strom_compare_local_node_major_to_owned_serial(
            &operator_context,
            stdd_effective_local,
            stdd_effective_owned);

        effective_rhs_window_difference[window] =
            strom_compare_relative_distributed_arrays(
                stdd_effective_owned,
                effective_owned_rhs,
                (size_t)owned_st_dof,
                communicator);
        max_effective_rhs_difference = fmax(
            max_effective_rhs_difference,
            effective_rhs_window_difference[window]);
    }

    if(ROM_std_strom_solve_sequential(
        &strom_system,
        NULL,
        NULL) != ROM_STROM_SUCCESS)
    {
        fail_message("sequential ST-ROM reduced solve failed");
    }

    for(int window = 0; window < num_windows; window++){
        ROM_STROM_WINDOW* reference_window =
            &strom_system.windows[window];
        double* saved_reference = reference_window->reference_state;
        double effective_rhs_difference;
        double g_difference;
        double a_difference;
        double q_difference;
        double q_right_difference;

        effective_rhs_difference =
            effective_rhs_window_difference[window];

        strom_compare_collect_stdd_global_vector(
            stdd_system->reduced_rhs
                + (size_t)window
                    * (size_t)stdd_system->num_modes_capacity,
            local_modes,
            global_modes,
            my_rank,
            stdd_g,
            communicator);
        memcpy(
            strom_g,
            reference_window->reduced_rhs,
            (size_t)global_modes * sizeof(double));
        g_difference = strom_compare_relative_global_arrays(
            stdd_g,
            strom_g,
            (size_t)global_modes);

        strom_compare_collect_stdd_global_vector(
            stdd_system->reduced_coef
                + (size_t)window
                    * (size_t)stdd_system->num_modes_capacity,
            local_modes,
            global_modes,
            my_rank,
            stdd_a,
            communicator);
        a_difference = strom_compare_relative_global_arrays(
            stdd_a,
            reference_window->reduced_coef,
            (size_t)global_modes);

        if(ROM_std_stdd_reconstruct_local_homogeneous(
            stdd_system,
            window,
            stdd_q_local) != ROM_STDD_SUCCESS)
        {
            fail_message("ST-DDROM q reconstruction failed");
        }
        strom_compare_local_node_major_to_owned_serial(
            &operator_context,
            stdd_q_local,
            stdd_q_owned);

        /* q excludes the Dirichlet reference/lift. */
        reference_window->reference_state = NULL;
        if(ROM_std_strom_reconstruct_window(
            reference_window,
            strom_q_owned) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM q reconstruction failed");
        }
        reference_window->reference_state = saved_reference;

        q_difference = strom_compare_relative_distributed_arrays(
            stdd_q_owned,
            strom_q_owned,
            (size_t)owned_st_dof,
            communicator);

        memset(
            stdd_trace_local,
            0,
            (size_t)num_local * sizeof(double));
        STSB_extract_last_right_trace(
            time_element,
            num_window_slabs,
            num_local,
            stdd_q_local,
            stdd_trace_local);
        memcpy(
            stdd_q_right,
            stdd_trace_local,
            (size_t)num_owned * sizeof(double));

        if(ROM_std_strom_extract_right_trace(
            strom_q_owned,
            strom_q_right,
            num_owned,
            num_window_slabs,
            time_element->n_dof,
            time_element->e_plus) != ROM_STROM_SUCCESS)
        {
            fail_message("sequential ST-ROM right-trace extraction failed");
        }

        q_right_difference = strom_compare_relative_distributed_arrays(
            stdd_q_right,
            strom_q_right,
            (size_t)num_owned,
            communicator);

        max_g_difference = fmax(max_g_difference, g_difference);
        max_a_difference = fmax(max_a_difference, a_difference);
        max_q_difference = fmax(max_q_difference, q_difference);
        max_q_right_difference = fmax(
            max_q_right_difference,
            q_right_difference);

        if(my_rank == 0){
            fprintf(
                csv,
                "%d,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                window,
                (double)(window + 1)
                    * (double)num_window_slabs * sys->vals.dt,
                D_difference,
                L_difference,
                effective_rhs_difference,
                g_difference,
                a_difference,
                q_difference,
                q_right_difference);
        }
    }

    if(my_rank == 0){
        fprintf(
            csv,
            "# summary,degree=%d,slabs_per_window=%d,global_modes=%d,"
            "D=%.15e,L=%.15e,max_effective_rhs=%.15e,max_g=%.15e,"
            "max_a=%.15e,max_q=%.15e,max_q_right=%.15e,"
            "tolerance=%.15e,status=%s\n",
            time_element->degree,
            num_window_slabs,
            global_modes,
            D_difference,
            L_difference,
            max_effective_rhs_difference,
            max_g_difference,
            max_a_difference,
            max_q_difference,
            max_q_right_difference,
            options->strom_compare_tolerance,
            (D_difference <= options->strom_compare_tolerance &&
             L_difference <= options->strom_compare_tolerance &&
             max_effective_rhs_difference <= options->strom_compare_tolerance &&
             max_g_difference <= options->strom_compare_tolerance &&
             max_a_difference <= options->strom_compare_tolerance &&
             max_q_difference <= options->strom_compare_tolerance &&
             max_q_right_difference <= options->strom_compare_tolerance)
                ? "MATCH"
                : "MISMATCH");
        fclose(csv);
        csv = NULL;

        printf(
            "\nST-DDROM vs validated sequential ST-ROM comparison\n"
            "  exact shared basis              : yes (owned rows only)\n"
            "  ST ordering                     : DDROM=node-major, ST-ROM=time-major\n"
            "  dG degree                       : %d\n"
            "  slabs/window                    : %d\n"
            "  D relative difference           : %.15e\n"
            "  L relative difference           : %.15e\n"
            "  max effective RHS difference    : %.15e\n"
            "  max g relative difference       : %.15e\n"
            "  max a relative difference       : %.15e\n"
            "  max q relative difference       : %.15e\n"
            "  max q^R relative difference     : %.15e\n"
            "  comparison tolerance            : %.15e\n"
            "  equivalence status              : %s\n"
            "  CSV                             : %s\n",
            time_element->degree,
            num_window_slabs,
            D_difference,
            L_difference,
            max_effective_rhs_difference,
            max_g_difference,
            max_a_difference,
            max_q_difference,
            max_q_right_difference,
            options->strom_compare_tolerance,
            (D_difference <= options->strom_compare_tolerance &&
             L_difference <= options->strom_compare_tolerance &&
             max_effective_rhs_difference <= options->strom_compare_tolerance &&
             max_g_difference <= options->strom_compare_tolerance &&
             max_a_difference <= options->strom_compare_tolerance &&
             max_q_difference <= options->strom_compare_tolerance &&
             max_q_right_difference <= options->strom_compare_tolerance)
                ? "MATCH"
                : "MISMATCH",
            options->strom_compare_csv);
    }

    strom_compare_write_matrix_csv(
        options->strom_compare_csv,
        "D",
        stdd_D,
        strom_D,
        global_modes);
    strom_compare_write_matrix_csv(
        options->strom_compare_csv,
        "L",
        stdd_L,
        strom_L,
        global_modes);

    free(operator_context.local_input);
    free(operator_context.local_output);
    free(stdd_D);
    free(stdd_L);
    free(strom_D);
    free(strom_L);
    free(stdd_g);
    free(strom_g);
    free(stdd_a);
    free(raw_local_rhs);
    free(raw_owned_rhs);
    free(effective_owned_rhs);
    free(stdd_effective_local);
    free(stdd_effective_owned);
    free(lift_local);
    free(initial_trace);
    free(zero_trace);
    free(stdd_q_local);
    free(stdd_q_owned);
    free(strom_q_owned);
    free(stdd_trace_local);
    free(stdd_q_right);
    free(strom_q_right);
    free(effective_rhs_window_difference);
    ROM_std_strom_system_finalize(&strom_system);
    free(strom_values);
    free(strom_matrices);
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

    if(options->compare_strom_reference){
        run_sequential_strom_comparison(
            sys,
            time_element,
            num_window_slabs,
            system,
            adapter,
            options);
    }

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
    const int compare_hlpod_kernels =
        env_int("ST_DDROM_COMPARE_HLPOD_KERNELS", 0);
    const int compare_dg0_standard =
        env_int("ST_DDROM_COMPARE_STANDARD_DDROM", 0) ||
        compare_hlpod_kernels;
    const int degree = (standard_ddrom || compare_dg0_standard)
        ? 0
        : env_int("ST_DDROM_TIME_DEGREE", ST_DDROM_DEGREE);
    const int num_window_slabs = (standard_ddrom || compare_dg0_standard)
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

    if(standard_ddrom > 1 || compare_hlpod_kernels > 1 ||
       degree < 0 || num_window_slabs <= 0)
    {
        fail_message("invalid ST-DDROM time discretization option");
    }

    if(standard_ddrom && compare_dg0_standard){
        fail_message(
            "ST_DDROM_STANDARD_DDROM and "
            "ST_DDROM_COMPARE_STANDARD_DDROM are mutually exclusive");
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

    {
        int expected_owned_space_nodes = 0;
        const int comm_size = monolis_mpi_get_global_comm_size();

        if(ST_fom_resolve_owned_space_points(
            sys.fe.total_num_nodes,
            comm_size,
            sys.monolis_com.n_internal_vertex,
            &expected_owned_space_nodes) != ST_FOM_SUCCESS)
        {
            fail_message(
                "failed to resolve the reference-FOM owned-node count");
        }

        if(num_internal_space_nodes != expected_owned_space_nodes)
        {
            fprintf(
                stderr,
                "ERROR [rank %d]: reference-FOM owned-node contract failed: "
                "returned=%d expected-owned=%d communication-owned=%d "
                "local=%d comm_size=%d\n",
                rank,
                num_internal_space_nodes,
                expected_owned_space_nodes,
                sys.monolis_com.n_internal_vertex,
                sys.fe.total_num_nodes,
                comm_size);
            fail_message(
                "reference FOM and common ST layout disagree on the owned "
                "vertex prefix");
        }
    }

    read_options(&options, sys.cond.directory);

    if((standard_ddrom || compare_dg0_standard) &&
       options.formulation != reference_fom_formulation())
    {
        fail_message(
            "ROM formulation must match the reference-FOM contract. "
            "Use ST_DDROM_FOM_FORMULATION=auto, or change "
            "STSB_reference_solution_is_incremental() in the FOM adapter.");
    }

    if(compare_dg0_standard && options.mode != ST_DDROM_RUN_ONLINE){
        fail_message(
            "ST_DDROM_COMPARE_STANDARD_DDROM=1 requires "
            "ST_DDROM_ACTION=online");
    }

    if(options.compare_strom_reference &&
       options.mode != ST_DDROM_RUN_ONLINE)
    {
        fail_message(
            "ST_DDROM_COMPARE_STROM=1 requires ST_DDROM_ACTION=online");
    }

    if(options.compare_strom_reference &&
       (standard_ddrom || compare_dg0_standard))
    {
        fail_message(
            "ST_DDROM_COMPARE_STROM must use the general ST-DDROM path; "
            "disable ST_DDROM_STANDARD_DDROM, "
            "ST_DDROM_COMPARE_STANDARD_DDROM, and "
            "ST_DDROM_COMPARE_HLPOD_KERNELS");
    }

    if(compare_dg0_standard &&
       options.formulation != ST_DDROM_FORMULATION_CURRENT_STATE)
    {
        fail_message(
            "the dG0 equivalence comparison requires current_state");
    }

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
                : (compare_dg0_standard
                    ? "dG(0) ST-DDROM / standard-DDROM equivalence test"
                    : "space-time DDROM"));
        printf("time dG degree: %d\n", degree);
        printf("slabs/window: %d\n", num_window_slabs);
        printf("num windows: %d\n", num_windows);
        printf(
            "comparison VTK: %s%s\n",
            options.write_comparison_vtk ? "enabled" : "disabled",
            options.write_comparison_vtk && !options.validate
                ? " (requires validation; no files will be written)"
                : "");
        printf(
            "POD mode VTK: %s, modes/rank limit=%s\n",
            options.write_mode_vtk ? "enabled" : "disabled",
            options.mode_vtk_max_modes > 0 ? "user-specified" : "all");
        printf(
            "FOM/ROM algebraic formulation: %s "
            "(reference contract: %s)\n",
            formulation_name(options.formulation),
            formulation_name(reference_fom_formulation()));
        printf("spatial MPI ranks: %d\n", monolis_mpi_get_global_comm_size());
        printf("max modes/rank: %d\n", options.max_modes);
        printf(
            "reduced solver: %s\n",
            reduced_solver_mode_name(options.reduced_solver_mode));
        printf("snapshot directory: %s\n", options.snapshot_dir);
        printf("basis directory: %s\n", options.basis_dir);
        printf(
            "ROM operator form: %s\n",
            standard_ddrom
                ? (options.formulation ==
                   ST_DDROM_FORMULATION_INCREMENTAL
                    ? "primitive A; increment RHS b_tilde-A_tilde U_prev; "
                      "U=U_prev+Delta U; explicit halo update + mono_com0 matvec"
                    : "primitive A; strong-BC current-state RHS; "
                      "U=V a+U_D; explicit halo update + mono_com0 matvec")
                : "matched strong-BC pair (A_tilde, b_tilde); "
                  "explicit halo update + mono_com0 matvec");
        printf(
            "external FOM residual check: %s, strict rejection: %s\n",
            options.check_fom_residual ? "enabled" : "disabled",
            standard_ddrom
                ? "ignored in verified-HLPOD standard mode"
                : (options.strict_fom_residual ? "enabled" : "disabled"));
        if(compare_dg0_standard){
            printf(
                "dG0 standard-DDROM comparison: enabled, tolerance=%.3e\n",
                options.compare_tolerance);
            printf("comparison CSV: %s\n", options.compare_csv);
            if(options.compare_hlpod_kernels){
                printf(
                    "verified HLPOD kernel comparison: enabled, "
                    "tolerance=%.3e\n",
                    options.hlpod_compare_tolerance);
                printf(
                    "HLPOD kernel comparison CSV: %s\n",
                    options.hlpod_compare_csv);
            }
        }
        if(options.compare_strom_reference){
            printf(
                "validated sequential ST-ROM comparison: enabled, "
                "tolerance=%.3e\n",
                options.strom_compare_tolerance);
            printf(
                "sequential ST-ROM comparison CSV: %s\n",
                options.strom_compare_csv);
        }
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
            standard_ddrom,
            &options,
            &adapter,
            &rom_system);
        break;

    case ST_DDROM_RUN_OFFLINE:
        run_offline(
            &rom_system,
            &adapter,
            standard_ddrom,
            &options);
        break;

    case ST_DDROM_RUN_ONLINE:
        if(compare_dg0_standard){
            run_dg0_standard_comparison(
                &sys,
                &time_element,
                num_window_slabs,
                &rom_system,
                &adapter,
                &options);
        }
        else if(standard_ddrom){
            run_online_standard_ddrom(
                &sys,
                &time_element,
                num_window_slabs,
                &rom_system,
                &adapter,
                &options);
        }
        else{
            run_online(
                &sys,
                &time_element,
                num_window_slabs,
                &rom_system,
                &adapter,
                &options);
        }
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
