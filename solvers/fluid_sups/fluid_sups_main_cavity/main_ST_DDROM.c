#include "fluid_sups_st_model.h"
#include "rom_std_stddrom.h"
#include "rom_std_stddrom_fluid.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    FLUID_ST_DDROM_OFFLINE = 0,
    FLUID_ST_DDROM_ONLINE  = 1,
    FLUID_ST_DDROM_FULL    = 2
} FLUID_ST_DDROM_MODE;

static void ddrom_fail(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static int ddrom_env_int(
    const char* name,
    int default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    long value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }

    errno = 0;
    value = strtol(text, &end, 10);

    if(errno != 0 || end == text ||
       *end != '\0' ||
       value < INT32_MIN ||
       value > INT32_MAX)
    {
        fprintf(
            stderr,
            "ERROR: invalid integer %s=%s\n",
            name,
            text);
        exit(EXIT_FAILURE);
    }

    return (int)value;
}

static double ddrom_env_double(
    const char* name,
    double default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    double value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }

    errno = 0;
    value = strtod(text, &end);

    if(errno != 0 || end == text ||
       *end != '\0' || !isfinite(value))
    {
        fprintf(
            stderr,
            "ERROR: invalid floating-point %s=%s\n",
            name,
            text);
        exit(EXIT_FAILURE);
    }

    return value;
}

static const char* ddrom_env_string(
    const char* name,
    const char* default_value)
{
    const char* value = getenv(name);
    return value != NULL && value[0] != '\0'
        ? value
        : default_value;
}

static FLUID_ST_DDROM_MODE ddrom_read_mode(void)
{
    const char* value = ddrom_env_string(
        "FLUID_ST_DDROM_MODE",
        "online");

    if(strcmp(value, "offline") == 0){
        return FLUID_ST_DDROM_OFFLINE;
    }
    if(strcmp(value, "online") == 0){
        return FLUID_ST_DDROM_ONLINE;
    }
    if(strcmp(value, "full") == 0){
        return FLUID_ST_DDROM_FULL;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_DDROM_MODE must be "
        "offline, online, or full; got %s\n",
        value);
    exit(EXIT_FAILURE);
}

static int ddrom_total_steps(
    double finish_time,
    double dt)
{
    const double value = finish_time / dt;
    const long long rounded = llround(value);

    if(dt <= 0.0 || finish_time <= 0.0 ||
       rounded <= 0 ||
       fabs(value - (double)rounded)
           > 1.0e-10 * fmax(1.0, fabs(value)) ||
       rounded > INT32_MAX)
    {
        ddrom_fail(
            "finish_time/dt must be a positive integer");
    }

    return (int)rounded;
}

static int ddrom_make_owned_window_filename(
    char* filename,
    size_t filename_size,
    const char* directory,
    int window_id)
{
    const int rank =
        monolis_mpi_get_global_my_rank();
    int written;

    if(filename == NULL || filename_size == 0 ||
       directory == NULL || directory[0] == '\0' ||
       window_id < 0)
    {
        return 1;
    }

    written = snprintf(
        filename,
        filename_size,
        "%s/fluid_st_window_%04d_rank_%06d.bin",
        directory,
        window_id,
        rank);

    if(written < 0 ||
       (size_t)written >= filename_size)
    {
        return 1;
    }

    return 0;
}

static int ddrom_local_owned_window_exists(
    const char* directory,
    int window_id)
{
    char filename[4096];
    FILE* fp;

    if(ddrom_make_owned_window_filename(
        filename,
        sizeof(filename),
        directory,
        window_id) != 0)
    {
        return 0;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        return 0;
    }

    fclose(fp);
    return 1;
}

/*
 * Return the number of consecutive FOM window files available on EVERY
 * spatial MPI rank, starting at window 0.
 *
 * The ST-DDROM online horizon and the POD training horizon are intentionally
 * independent.  A basis may be trained with fewer windows than the online
 * solve, although for production accuracy the training interval should cover
 * the intended parameter/time regime.
 */
static int ddrom_count_common_training_windows(
    const char* directory,
    int maximum_windows,
    int communicator)
{
    int local_count = 0;
    int negative_count;

    if(directory == NULL || maximum_windows <= 0){
        return 0;
    }

    while(local_count < maximum_windows &&
          ddrom_local_owned_window_exists(
              directory,
              local_count))
    {
        local_count++;
    }

    /*
     * Need MIN(local_count), while the MONOLIS interface used in this project
     * reliably exposes MAX.  Therefore:
     *
     *   min(c_r) = -max(-c_r).
     */
    negative_count = -local_count;

    monolis_allreduce_I(
        1,
        &negative_count,
        MONOLIS_MPI_MAX,
        communicator);

    return -negative_count;
}

static int ddrom_read_owned_window(
    const char* directory,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* owned_window)
{
    ST_FOM_WINDOW_FILE_HEADER header;
    char filename[4096];
    const int rank =
        monolis_mpi_get_global_my_rank();
    const size_t count =
        (size_t)layout->num_owned_space_points
        * (size_t)ST_fom_window_dof_per_space_point(
            layout);
    FILE* fp;

    if(ddrom_make_owned_window_filename(
        filename,
        sizeof(filename),
        directory,
        window_id) != 0)
    {
        return 1;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        fprintf(
            stderr,
            "ERROR [rank %d]: cannot read FOM snapshot %s: %s\n",
            rank,
            filename,
            strerror(errno));
        return 1;
    }

    if(fread(&header, sizeof(header), 1, fp) != 1 ||
       header.magic != ST_FOM_WINDOW_MAGIC ||
       header.version != ST_FOM_WINDOW_VERSION ||
       header.window_id != (uint32_t)window_id ||
       header.degree != (uint32_t)layout->time.degree ||
       header.num_slabs != (uint32_t)layout->num_slabs ||
       header.num_time_dofs !=
           (uint32_t)layout->time.n_dof ||
       header.physical_dof !=
           (uint32_t)layout->physical_dof ||
       header.num_local_space_points !=
           (uint32_t)layout->num_space_points ||
       header.num_owned_space_points !=
           (uint32_t)layout->num_owned_space_points ||
       fread(
           owned_window,
           sizeof(double),
           count,
           fp) != count)
    {
        fclose(fp);
        fprintf(
            stderr,
            "ERROR [rank %d]: incompatible FOM window snapshot %s\n",
            rank,
            filename);
        return 1;
    }

    if(fclose(fp) != 0){
        return 1;
    }

    return 0;
}

static int ddrom_build_snapshot_matrix(
    ROM_STDD_FLUID_CONTEXT* adapter,
    const char* snapshot_directory,
    int num_training_windows,
    double** snapshot_matrix_out)
{
    ROM_STDD_SYSTEM* reduced =
        adapter->reduced;
    const int columns =
        num_training_windows;
    double* snapshot_matrix = NULL;
    double* owned_window = NULL;

    if(columns <= 0 ||
       columns > reduced->num_windows ||
       snapshot_matrix_out == NULL)
    {
        return 1;
    }

    snapshot_matrix = (double*)calloc(
        (size_t)reduced->internal_st_dof
            * (size_t)columns,
        sizeof(double));

    owned_window = (double*)calloc(
        (size_t)reduced->internal_st_dof,
        sizeof(double));

    if(snapshot_matrix == NULL ||
       owned_window == NULL)
    {
        free(snapshot_matrix);
        free(owned_window);
        return 1;
    }

    if(fluid_sups_st_rom_build_window_lift(
        adapter->fom,
        adapter->layout,
        adapter->lift) != 0)
    {
        free(snapshot_matrix);
        free(owned_window);
        return 1;
    }

    for(int window = 0; window < columns; window++){
        if(ddrom_read_owned_window(
            snapshot_directory,
            window,
            adapter->layout,
            owned_window) != 0)
        {
            free(snapshot_matrix);
            free(owned_window);
            return 1;
        }

        for(int row = 0;
            row < reduced->internal_st_dof;
            row++)
        {
            snapshot_matrix[
                (size_t)row * (size_t)columns
                + (size_t)window] =
                owned_window[row]
                - adapter->lift[row];
        }
    }

    free(owned_window);
    *snapshot_matrix_out = snapshot_matrix;
    return 0;
}




static size_t ddrom_basis_offset(
    const ROM_STDD_SYSTEM* system,
    int row,
    int mode)
{
    return (size_t)row
        * (size_t)system->num_modes_capacity
        + (size_t)mode;
}

static int ddrom_classify_split_basis_modes(
    const ROM_STDD_SYSTEM* reduced,
    int* velocity_modes,
    int* pressure_modes)
{
    const double support_tolerance = 1.0e-24;
    const int temporal_slots =
        reduced->st_dof_per_space_node
        / FLUID_SUPS_PHYSICAL_DOF;
    int count_velocity = 0;
    int count_pressure = 0;

    if(reduced == NULL || reduced->basis == NULL ||
       reduced->num_modes <= 0 ||
       reduced->st_dof_per_space_node
            % FLUID_SUPS_PHYSICAL_DOF != 0 ||
       temporal_slots <= 0 ||
       velocity_modes == NULL ||
       pressure_modes == NULL)
    {
        return 1;
    }

    for(int mode = 0; mode < reduced->num_modes; mode++){
        double velocity_norm2 = 0.0;
        double pressure_norm2 = 0.0;

        for(int node = 0;
            node < reduced->num_internal_space_nodes;
            node++)
        {
            for(int temporal = 0;
                temporal < temporal_slots;
                temporal++)
            {
                const int base_row =
                    node * reduced->st_dof_per_space_node
                    + temporal * FLUID_SUPS_PHYSICAL_DOF;

                for(int component = 0; component < 3; component++){
                    const double value =
                        reduced->basis[
                            ddrom_basis_offset(
                                reduced,
                                base_row + component,
                                mode)];
                    velocity_norm2 += value * value;
                }

                {
                    const double value =
                        reduced->basis[
                            ddrom_basis_offset(
                                reduced,
                                base_row + 3,
                                mode)];
                    pressure_norm2 += value * value;
                }
            }
        }

        if(velocity_norm2 > support_tolerance &&
           pressure_norm2 <= support_tolerance)
        {
            count_velocity++;
        }
        else if(pressure_norm2 > support_tolerance &&
                velocity_norm2 <= support_tolerance)
        {
            count_pressure++;
        }
        else{
            return 1;
        }
    }

    *velocity_modes = count_velocity;
    *pressure_modes = count_pressure;
    return 0;
}

/*
 * Build the fluid space-time POD basis in the same structural form used by
 * the ordinary fluid HLPOD-DDROM:
 *
 *          Phi = diag(Phi_v, Phi_p)
 *
 * where Phi_v and Phi_p are computed independently.
 *
 * The physical ST row ordering remains
 *
 *   [node][slab][time basis][ux,uy,uz,p].
 *
 * Velocity modes therefore have exactly zero pressure support, and pressure
 * modes have exactly zero velocity support.  The reduced coordinates
 *
 *          a = [a_v ; a_p]
 *
 * are independent, preserving the velocity/pressure saddle-point freedom
 * that is lost when one coupled 4-component POD is formed.
 */
static int ddrom_compute_separated_velocity_pressure_pod(
    ROM_STDD_FLUID_CONTEXT* adapter,
    const double* full_snapshot_matrix,
    int num_snapshot_columns,
    int max_velocity_modes,
    int max_pressure_modes,
    double velocity_energy_tolerance,
    double pressure_energy_tolerance,
    int communicator,
    int* velocity_modes_out,
    int* pressure_modes_out)
{
    ROM_STDD_SYSTEM velocity_system;
    ROM_STDD_SYSTEM pressure_system;
    ROM_STDD_SYSTEM* reduced;
    const int physical_dof = FLUID_SUPS_PHYSICAL_DOF;
    int temporal_slots;
    int velocity_dof_per_node;
    int pressure_dof_per_node;
    int total_modes;
    double* velocity_snapshots = NULL;
    double* pressure_snapshots = NULL;
    size_t velocity_rows;
    size_t pressure_rows;
    int status = ROM_STDD_SUCCESS;

    memset(&velocity_system, 0, sizeof(velocity_system));
    memset(&pressure_system, 0, sizeof(pressure_system));

    if(adapter == NULL || adapter->reduced == NULL ||
       adapter->fom == NULL || full_snapshot_matrix == NULL ||
       num_snapshot_columns <= 0 ||
       max_velocity_modes <= 0 ||
       max_pressure_modes <= 0 ||
       velocity_modes_out == NULL ||
       pressure_modes_out == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    reduced = adapter->reduced;

    if(reduced->st_dof_per_space_node % physical_dof != 0){
        return ROM_STDD_ERR_DIMENSION;
    }

    temporal_slots =
        reduced->st_dof_per_space_node / physical_dof;
    velocity_dof_per_node = 3 * temporal_slots;
    pressure_dof_per_node = temporal_slots;

    if(velocity_dof_per_node <= 0 ||
       pressure_dof_per_node <= 0 ||
       max_velocity_modes > reduced->num_modes_capacity ||
       max_pressure_modes > reduced->num_modes_capacity ||
       max_velocity_modes + max_pressure_modes
            > reduced->num_modes_capacity)
    {
        return ROM_STDD_ERR_DIMENSION;
    }

    velocity_rows =
        (size_t)reduced->num_internal_space_nodes
        * (size_t)velocity_dof_per_node;
    pressure_rows =
        (size_t)reduced->num_internal_space_nodes
        * (size_t)pressure_dof_per_node;

    velocity_snapshots = (double*)calloc(
        velocity_rows * (size_t)num_snapshot_columns,
        sizeof(double));
    pressure_snapshots = (double*)calloc(
        pressure_rows * (size_t)num_snapshot_columns,
        sizeof(double));

    if(velocity_snapshots == NULL ||
       pressure_snapshots == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    for(int node = 0;
        node < reduced->num_internal_space_nodes;
        node++)
    {
        for(int temporal = 0;
            temporal < temporal_slots;
            temporal++)
        {
            const int full_base =
                node * reduced->st_dof_per_space_node
                + temporal * physical_dof;
            const int velocity_base =
                node * velocity_dof_per_node
                + temporal * 3;
            const int pressure_row =
                node * pressure_dof_per_node
                + temporal;

            for(int column = 0;
                column < num_snapshot_columns;
                column++)
            {
                for(int component = 0;
                    component < 3;
                    component++)
                {
                    velocity_snapshots[
                        ((size_t)velocity_base
                            + (size_t)component)
                            * (size_t)num_snapshot_columns
                        + (size_t)column] =
                        full_snapshot_matrix[
                            ((size_t)full_base
                                + (size_t)component)
                                * (size_t)num_snapshot_columns
                            + (size_t)column];
                }

                pressure_snapshots[
                    (size_t)pressure_row
                        * (size_t)num_snapshot_columns
                    + (size_t)column] =
                    full_snapshot_matrix[
                        ((size_t)full_base + 3u)
                            * (size_t)num_snapshot_columns
                        + (size_t)column];
            }
        }
    }

    status = ROM_std_stdd_system_initialize(
        &velocity_system,
        &adapter->fom->mono_com,
        reduced->num_windows,
        reduced->num_space_nodes,
        reduced->num_internal_space_nodes,
        velocity_dof_per_node,
        max_velocity_modes);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    status = ROM_std_stdd_system_initialize(
        &pressure_system,
        &adapter->fom->mono_com,
        reduced->num_windows,
        reduced->num_space_nodes,
        reduced->num_internal_space_nodes,
        pressure_dof_per_node,
        max_pressure_modes);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    status = ROM_std_stdd_compute_shared_local_pod(
        &velocity_system,
        velocity_snapshots,
        num_snapshot_columns,
        velocity_energy_tolerance,
        communicator);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    status = ROM_std_stdd_compute_shared_local_pod(
        &pressure_system,
        pressure_snapshots,
        num_snapshot_columns,
        pressure_energy_tolerance,
        communicator);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    total_modes =
        velocity_system.num_modes
        + pressure_system.num_modes;

    if(total_modes <= 0 ||
       total_modes > reduced->num_modes_capacity)
    {
        status = ROM_STDD_ERR_DIMENSION;
        goto cleanup;
    }

    free(reduced->basis);
    free(reduced->singular_values);
    reduced->basis = NULL;
    reduced->singular_values = NULL;

    reduced->basis = (double*)calloc(
        (size_t)reduced->internal_st_dof
            * (size_t)reduced->num_modes_capacity,
        sizeof(double));

    /*
     * Keep both singular-value spectra in the generic basis file:
     *
     *   [sigma_v(0..Ns-1), sigma_p(0..Ns-1)].
     *
     * The online solver only requires the basis itself, but preserving both
     * spectra is useful for offline diagnostics.
     */
    reduced->singular_values = (double*)calloc(
        (size_t)(2 * num_snapshot_columns),
        sizeof(double));

    if(reduced->basis == NULL ||
       reduced->singular_values == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    reduced->num_modes = total_modes;
    reduced->num_singular_values =
        2 * num_snapshot_columns;

    memcpy(
        reduced->singular_values,
        velocity_system.singular_values,
        (size_t)num_snapshot_columns * sizeof(double));
    memcpy(
        reduced->singular_values
            + num_snapshot_columns,
        pressure_system.singular_values,
        (size_t)num_snapshot_columns * sizeof(double));

    for(int node = 0;
        node < reduced->num_internal_space_nodes;
        node++)
    {
        for(int temporal = 0;
            temporal < temporal_slots;
            temporal++)
        {
            const int full_base =
                node * reduced->st_dof_per_space_node
                + temporal * physical_dof;
            const int velocity_base =
                node * velocity_dof_per_node
                + temporal * 3;
            const int pressure_row =
                node * pressure_dof_per_node
                + temporal;

            for(int velocity_mode = 0;
                velocity_mode < velocity_system.num_modes;
                velocity_mode++)
            {
                for(int component = 0;
                    component < 3;
                    component++)
                {
                    reduced->basis[
                        ddrom_basis_offset(
                            reduced,
                            full_base + component,
                            velocity_mode)] =
                        velocity_system.basis[
                            ddrom_basis_offset(
                                &velocity_system,
                                velocity_base + component,
                                velocity_mode)];
                }
            }

            for(int pressure_mode = 0;
                pressure_mode < pressure_system.num_modes;
                pressure_mode++)
            {
                reduced->basis[
                    ddrom_basis_offset(
                        reduced,
                        full_base + 3,
                        velocity_system.num_modes
                            + pressure_mode)] =
                    pressure_system.basis[
                        ddrom_basis_offset(
                            &pressure_system,
                            pressure_row,
                            pressure_mode)];
            }
        }
    }

    *velocity_modes_out =
        velocity_system.num_modes;
    *pressure_modes_out =
        pressure_system.num_modes;

cleanup:
    free(velocity_snapshots);
    free(pressure_snapshots);
    ROM_std_stdd_system_finalize(&velocity_system);
    ROM_std_stdd_system_finalize(&pressure_system);

    return status;
}

typedef struct {
    int num_validated_windows;

    double ref2;
    double projection_diff2;
    double state_diff2;
    double dynamic_diff2;

    double velocity_ref2;
    double velocity_projection_diff2;
    double velocity_state_diff2;
    double velocity_dynamic_diff2;

    double pressure_ref2;
    double pressure_projection_diff2;
    double pressure_state_diff2;
    double pressure_dynamic_diff2;

    double coefficient_ref2;
    double coefficient_diff2;

    double trace_ref2;
    double trace_diff2;
    double trace_velocity_ref2;
    double trace_velocity_diff2;
    double trace_pressure_ref2;
    double trace_pressure_diff2;

    double max_projection_error;
    double max_state_error;
    double max_dynamic_error;
    double max_velocity_state_error;
    double max_pressure_state_error;
    double max_coefficient_error;
    double max_trace_error;

    double max_reference_full_residual;
    double max_reference_projected_defect;
} DDROM_ACCURACY_SUMMARY;

typedef struct {
    double projection_error;
    double state_error;
    double dynamic_error;

    double velocity_projection_error;
    double velocity_state_error;
    double velocity_dynamic_error;

    double pressure_projection_error;
    double pressure_state_error;
    double pressure_dynamic_error;

    double coefficient_error;

    double trace_error;
    double trace_velocity_error;
    double trace_pressure_error;

    double reference_full_residual;
    double reference_momentum_residual;
    double reference_mass_residual;
    double reference_projected_defect;
} DDROM_WINDOW_ACCURACY;

static double ddrom_safe_ratio_from_squares(
    double numerator2,
    double denominator2)
{
    return sqrt(fmax(numerator2, 0.0))
        / fmax(sqrt(fmax(denominator2, 0.0)), 1.0e-300);
}

static double ddrom_reduced_vector_global_norm(
    const ROM_STDD_FLUID_CONTEXT* adapter,
    const double* vector)
{
    double norm2 = 0.0;

    for(int mode = 0;
        mode < adapter->reduced->num_modes;
        mode++)
    {
        norm2 += vector[mode] * vector[mode];
    }

    monolis_allreduce_R(
        1,
        &norm2,
        MONOLIS_MPI_SUM,
        adapter->fom->mono_com.comm);

    return sqrt(fmax(norm2, 0.0));
}

static int ddrom_reconstruct_reference_full(
    ROM_STDD_FLUID_CONTEXT* adapter,
    const double* reference_owned,
    double* reference_full)
{
    if(adapter == NULL || reference_owned == NULL ||
       reference_full == NULL)
    {
        return 1;
    }

    memset(
        reference_full,
        0,
        (size_t)adapter->reduced->local_st_dof
            * sizeof(double));

    /*
     * MONOLIS ownership contract for this distributed fluid path:
     * owned nodes are the local prefix, followed by halo nodes.
     */
    memcpy(
        reference_full,
        reference_owned,
        (size_t)adapter->reduced->internal_st_dof
            * sizeof(double));

    monolis_mpi_update_R(
        &adapter->fom->mono_com,
        adapter->fom->fe.total_num_nodes,
        adapter->fom->window_dof,
        reference_full);

    return 0;
}

static int ddrom_validate_current_window(
    ROM_STDD_FLUID_CONTEXT* adapter,
    int window_id,
    const double* reference_owned,
    double* reference_full,
    double* reference_homogeneous,
    double* projected_homogeneous,
    double* reference_coefficients,
    double* reference_previous_trace,
    double* reference_trace,
    double* rom_trace,
    double* reference_projected_rhs,
    DDROM_WINDOW_ACCURACY* accuracy,
    DDROM_ACCURACY_SUMMARY* summary)
{
    ROM_STDD_SYSTEM* reduced;
    double local[20] = {0.0};
    double reference_total_residual = 0.0;
    double reference_momentum_residual = 0.0;
    double reference_mass_residual = 0.0;
    const double* rom_coefficients;

    if(adapter == NULL || adapter->reduced == NULL ||
       adapter->layout == NULL || adapter->fom == NULL ||
       reference_owned == NULL || reference_full == NULL ||
       reference_homogeneous == NULL ||
       projected_homogeneous == NULL ||
       reference_coefficients == NULL ||
       reference_previous_trace == NULL ||
       reference_trace == NULL || rom_trace == NULL ||
       reference_projected_rhs == NULL ||
       accuracy == NULL || summary == NULL ||
       window_id < 0 ||
       window_id >= adapter->reduced->num_windows)
    {
        return 1;
    }

    reduced = adapter->reduced;
    memset(accuracy, 0, sizeof(*accuracy));

    if(fluid_sups_st_rom_build_window_lift(
           adapter->fom,
           adapter->layout,
           adapter->lift) != 0 ||
       ddrom_reconstruct_reference_full(
           adapter,
           reference_owned,
           reference_full) != 0)
    {
        return 1;
    }

    memset(
        reference_homogeneous,
        0,
        (size_t)reduced->internal_st_dof
            * sizeof(double));
    memset(
        projected_homogeneous,
        0,
        (size_t)reduced->internal_st_dof
            * sizeof(double));
    memset(
        reference_coefficients,
        0,
        (size_t)reduced->num_modes_capacity
            * sizeof(double));

    for(int row = 0; row < reduced->internal_st_dof; row++){
        reference_homogeneous[row] =
            reference_owned[row] - adapter->lift[row];
    }

    if(ROM_std_stdd_project_reference_coefficients(
           reduced,
           window_id,
           reference_homogeneous,
           reference_coefficients) != ROM_STDD_SUCCESS)
    {
        return 1;
    }

    for(int row = 0; row < reduced->internal_st_dof; row++){
        double value = 0.0;

        for(int mode = 0; mode < reduced->num_modes; mode++){
            value += reduced->basis[
                (size_t)row
                    * (size_t)reduced->num_modes_capacity
                + (size_t)mode]
                * reference_coefficients[mode];
        }

        projected_homogeneous[row] = value;
    }

    rom_coefficients =
        reduced->reduced_coef
        + (size_t)window_id
            * (size_t)reduced->num_modes_capacity;

    for(int row = 0; row < reduced->internal_st_dof; row++){
        const int physical_component =
            (row % reduced->st_dof_per_space_node)
                % FLUID_SUPS_PHYSICAL_DOF;
        const double q_ref =
            reference_homogeneous[row];
        const double q_projection =
            projected_homogeneous[row];
        const double q_rom =
            adapter->current_window[row]
            - adapter->lift[row];

        const double projection_difference =
            q_projection - q_ref;
        const double state_difference =
            q_rom - q_ref;
        const double dynamic_difference =
            q_rom - q_projection;

        local[0] += q_ref * q_ref;
        local[1] +=
            projection_difference
            * projection_difference;
        local[2] +=
            state_difference
            * state_difference;
        local[3] +=
            dynamic_difference
            * dynamic_difference;

        if(physical_component < 3){
            local[4] += q_ref * q_ref;
            local[5] +=
                projection_difference
                * projection_difference;
            local[6] +=
                state_difference
                * state_difference;
            local[7] +=
                dynamic_difference
                * dynamic_difference;
        }
        else{
            local[8] += q_ref * q_ref;
            local[9] +=
                projection_difference
                * projection_difference;
            local[10] +=
                state_difference
                * state_difference;
            local[11] +=
                dynamic_difference
                * dynamic_difference;
        }
    }

    for(int mode = 0; mode < reduced->num_modes; mode++){
        const double reference_value =
            reference_coefficients[mode];
        const double difference =
            rom_coefficients[mode]
            - reference_value;

        local[12] +=
            reference_value * reference_value;
        local[13] += difference * difference;
    }

    if(ST_fom_extract_last_right_trace(
           adapter->layout,
           reference_full,
           reference_trace) != ST_FOM_SUCCESS ||
       ST_fom_extract_last_right_trace(
           adapter->layout,
           adapter->current_window,
           rom_trace) != ST_FOM_SUCCESS)
    {
        return 1;
    }

    for(int node = 0;
        node < adapter->fom->num_owned_space_nodes;
        node++)
    {
        for(int component = 0;
            component < FLUID_SUPS_PHYSICAL_DOF;
            component++)
        {
            const int index =
                FLUID_SUPS_PHYSICAL_DOF * node
                + component;
            const double reference_value =
                reference_trace[index];
            const double difference =
                rom_trace[index]
                - reference_value;

            local[14] +=
                reference_value * reference_value;
            local[15] +=
                difference * difference;

            if(component < 3){
                local[16] +=
                    reference_value * reference_value;
                local[17] +=
                    difference * difference;
            }
            else{
                local[18] +=
                    reference_value * reference_value;
                local[19] +=
                    difference * difference;
            }
        }
    }

    monolis_allreduce_R(
        20,
        local,
        MONOLIS_MPI_SUM,
        adapter->fom->mono_com.comm);

    accuracy->projection_error =
        ddrom_safe_ratio_from_squares(
            local[1], local[0]);
    accuracy->state_error =
        ddrom_safe_ratio_from_squares(
            local[2], local[0]);
    accuracy->dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[3], local[0]);

    accuracy->velocity_projection_error =
        ddrom_safe_ratio_from_squares(
            local[5], local[4]);
    accuracy->velocity_state_error =
        ddrom_safe_ratio_from_squares(
            local[6], local[4]);
    accuracy->velocity_dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[7], local[4]);

    accuracy->pressure_projection_error =
        ddrom_safe_ratio_from_squares(
            local[9], local[8]);
    accuracy->pressure_state_error =
        ddrom_safe_ratio_from_squares(
            local[10], local[8]);
    accuracy->pressure_dynamic_error =
        ddrom_safe_ratio_from_squares(
            local[11], local[8]);

    accuracy->coefficient_error =
        ddrom_safe_ratio_from_squares(
            local[13], local[12]);

    accuracy->trace_error =
        ddrom_safe_ratio_from_squares(
            local[15], local[14]);
    accuracy->trace_velocity_error =
        ddrom_safe_ratio_from_squares(
            local[17], local[16]);
    accuracy->trace_pressure_error =
        ddrom_safe_ratio_from_squares(
            local[19], local[18]);

    /*
     * Evaluate the stored FOM reference with its own causal previous trace,
     * not with the ROM previous trace.  This verifies that the binary
     * reference itself remains a solution of the all-at-once residual.
     */
    if(fluid_sups_st_rom_assemble_window_residual(
           adapter->fom,
           adapter->layout,
           reference_previous_trace,
           reference_full) != 0 ||
       fluid_sups_st_rom_get_window_rhs_norms(
           adapter->fom,
           adapter->layout,
           &reference_total_residual,
           &reference_momentum_residual,
           &reference_mass_residual) != 0 ||
       ROM_std_stdd_project_local_rhs(
           reduced,
           window_id,
           adapter->fom->monolis_window.mat.R.B,
           reference_projected_rhs) != ROM_STDD_SUCCESS)
    {
        return 1;
    }

    accuracy->reference_full_residual =
        reference_total_residual;
    accuracy->reference_momentum_residual =
        reference_momentum_residual;
    accuracy->reference_mass_residual =
        reference_mass_residual;
    accuracy->reference_projected_defect =
        ddrom_reduced_vector_global_norm(
            adapter,
            reference_projected_rhs);

    memcpy(
        reference_previous_trace,
        reference_trace,
        (size_t)adapter->fom->fe.total_num_nodes
            * (size_t)FLUID_SUPS_PHYSICAL_DOF
            * sizeof(double));

    summary->num_validated_windows++;

    summary->ref2 += local[0];
    summary->projection_diff2 += local[1];
    summary->state_diff2 += local[2];
    summary->dynamic_diff2 += local[3];

    summary->velocity_ref2 += local[4];
    summary->velocity_projection_diff2 += local[5];
    summary->velocity_state_diff2 += local[6];
    summary->velocity_dynamic_diff2 += local[7];

    summary->pressure_ref2 += local[8];
    summary->pressure_projection_diff2 += local[9];
    summary->pressure_state_diff2 += local[10];
    summary->pressure_dynamic_diff2 += local[11];

    summary->coefficient_ref2 += local[12];
    summary->coefficient_diff2 += local[13];

    summary->trace_ref2 += local[14];
    summary->trace_diff2 += local[15];
    summary->trace_velocity_ref2 += local[16];
    summary->trace_velocity_diff2 += local[17];
    summary->trace_pressure_ref2 += local[18];
    summary->trace_pressure_diff2 += local[19];

    summary->max_projection_error =
        fmax(summary->max_projection_error,
             accuracy->projection_error);
    summary->max_state_error =
        fmax(summary->max_state_error,
             accuracy->state_error);
    summary->max_dynamic_error =
        fmax(summary->max_dynamic_error,
             accuracy->dynamic_error);
    summary->max_velocity_state_error =
        fmax(summary->max_velocity_state_error,
             accuracy->velocity_state_error);
    summary->max_pressure_state_error =
        fmax(summary->max_pressure_state_error,
             accuracy->pressure_state_error);
    summary->max_coefficient_error =
        fmax(summary->max_coefficient_error,
             accuracy->coefficient_error);
    summary->max_trace_error =
        fmax(summary->max_trace_error,
             accuracy->trace_error);
    summary->max_reference_full_residual =
        fmax(summary->max_reference_full_residual,
             accuracy->reference_full_residual);
    summary->max_reference_projected_defect =
        fmax(summary->max_reference_projected_defect,
             accuracy->reference_projected_defect);

    return 0;
}

static void ddrom_print_accuracy_summary(
    const DDROM_ACCURACY_SUMMARY* summary)
{
    if(summary == NULL ||
       monolis_mpi_get_global_my_rank() != 0)
    {
        return;
    }

    printf(
        "\n"
        "============================================================\n"
        "Fluid ST-DDROM accuracy summary\n"
        "  validated windows                    : %d\n"
        "  projection error ||qF-PhiPhiTqF||/||qF|| : %.15e\n"
        "  ROM state error  ||qF-qROM||/||qF||      : %.15e\n"
        "  dynamic error    ||PhiPhiTqF-qROM||/||qF||: %.15e\n"
        "  velocity projection error            : %.15e\n"
        "  velocity ROM state error             : %.15e\n"
        "  velocity dynamic error               : %.15e\n"
        "  pressure projection error            : %.15e\n"
        "  pressure ROM state error             : %.15e\n"
        "  pressure dynamic error               : %.15e\n"
        "  reduced coefficient error            : %.15e\n"
        "  right-trace error                    : %.15e\n"
        "  right-trace velocity error           : %.15e\n"
        "  right-trace pressure error           : %.15e\n"
        "  max window projection error          : %.15e\n"
        "  max window ROM state error           : %.15e\n"
        "  max window dynamic error             : %.15e\n"
        "  max window velocity state error      : %.15e\n"
        "  max window pressure state error      : %.15e\n"
        "  max window coefficient error         : %.15e\n"
        "  max window right-trace error         : %.15e\n"
        "  max FOM-reference full residual      : %.15e\n"
        "  max projected defect at FOM reference: %.15e\n"
        "============================================================\n",
        summary->num_validated_windows,
        ddrom_safe_ratio_from_squares(
            summary->projection_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->state_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->dynamic_diff2,
            summary->ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_projection_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_state_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->velocity_dynamic_diff2,
            summary->velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_projection_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_state_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->pressure_dynamic_diff2,
            summary->pressure_ref2),
        ddrom_safe_ratio_from_squares(
            summary->coefficient_diff2,
            summary->coefficient_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_diff2,
            summary->trace_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_velocity_diff2,
            summary->trace_velocity_ref2),
        ddrom_safe_ratio_from_squares(
            summary->trace_pressure_diff2,
            summary->trace_pressure_ref2),
        summary->max_projection_error,
        summary->max_state_error,
        summary->max_dynamic_error,
        summary->max_velocity_state_error,
        summary->max_pressure_state_error,
        summary->max_coefficient_error,
        summary->max_trace_error,
        summary->max_reference_full_residual,
        summary->max_reference_projected_defect);
}

static int ddrom_write_accuracy_header(FILE* fp)
{
    if(fp == NULL){
        return 1;
    }

    fprintf(
        fp,
        "window,"
        "projection_error,state_error,dynamic_error,"
        "velocity_projection_error,velocity_state_error,"
        "velocity_dynamic_error,"
        "pressure_projection_error,pressure_state_error,"
        "pressure_dynamic_error,"
        "coefficient_error,"
        "trace_error,trace_velocity_error,trace_pressure_error,"
        "rom_projected_residual,rom_relative_projected_residual,"
        "rom_full_residual,rom_relative_full_residual,"
        "reference_full_residual,reference_momentum_residual,"
        "reference_mass_residual,reference_projected_defect\n");

    return 0;
}

static int ddrom_write_coefficient_header(
    FILE* fp,
    int num_modes)
{
    if(fp == NULL || num_modes <= 0){
        return 1;
    }

    fprintf(fp, "window");
    for(int mode = 0; mode < num_modes; mode++){
        fprintf(fp, ",a_%04d", mode);
    }
    fprintf(
        fp,
        ",relative_full_residual,projected_residual,"
        "coefficient_correction\n");

    return 0;
}

int main(int argc, char* argv[])
{
    FLUID_SUPS_SYSTEM fom;
    ST_FOM_WINDOW_LAYOUT layout;
    ROM_STDD_SYSTEM reduced;
    ROM_STDD_FLUID_CONTEXT adapter;
    FLUID_ST_DDROM_MODE mode;

    const char* snapshot_directory;
    const char* basis_directory;
    const char* rom_window_directory;
    const char* coefficient_csv_name;
    const char* accuracy_csv_name;

    int total_steps;
    int num_windows;
    int online_windows = 0;
    int available_snapshot_windows = 0;
    int num_training_windows = 0;
    int requested_training_windows = 0;
    int max_modes;
    int max_velocity_modes;
    int max_pressure_modes;
    int total_mode_capacity;
    int velocity_modes = 0;
    int pressure_modes = 0;
    int file_number = 0;
    int write_output;
    int write_window_binary;
    int validate_accuracy;
    double pod_energy;
    double velocity_pod_energy;
    double pressure_pod_energy;

    FILE* coefficient_csv = NULL;
    FILE* accuracy_csv = NULL;
    double* snapshot_matrix = NULL;
    double* reference_owned = NULL;
    double* reference_full = NULL;
    double* reference_homogeneous = NULL;
    double* projected_homogeneous = NULL;
    double* reference_coefficients = NULL;
    double* reference_previous_trace = NULL;
    double* reference_trace = NULL;
    double* rom_trace = NULL;
    double* reference_projected_rhs = NULL;
    DDROM_ACCURACY_SUMMARY accuracy_summary;

    int status = EXIT_FAILURE;

    memset(&fom, 0, sizeof(fom));
    memset(&layout, 0, sizeof(layout));
    memset(&reduced, 0, sizeof(reduced));
    memset(&adapter, 0, sizeof(adapter));
    memset(&accuracy_summary, 0, sizeof(accuracy_summary));

    monolis_global_initialize();

    if(fluid_sups_st_initialize(
        &fom,
        argc,
        argv) != 0)
    {
        ddrom_fail("fluid ST-DDROM FOM initialization failed");
    }

    total_steps = ddrom_total_steps(
        fom.vals.finish_time,
        fom.vals.dt);

    if(total_steps % fom.vals.st_slabs_per_window != 0){
        ddrom_fail(
            "total steps must be divisible by slabs/window");
    }

    num_windows =
        total_steps / fom.vals.st_slabs_per_window;

    if(ST_fom_layout_initialize(
        &layout,
        fom.fe.total_num_nodes,
        fom.num_owned_space_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        fom.vals.st_slabs_per_window,
        fom.vals.st_time_degree) != ST_FOM_SUCCESS)
    {
        ddrom_fail("ST-DDROM layout initialization failed");
    }

    /*
     * Backward-compatible fallback:
     *
     * FLUID_ST_DDROM_MAX_MODES=N means N velocity modes AND N pressure
     * modes unless the field-specific limits below are provided.
     */
    max_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_MODES",
        10);

    max_velocity_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_VELOCITY_MODES",
        max_modes);

    max_pressure_modes = ddrom_env_int(
        "FLUID_ST_DDROM_MAX_PRESSURE_MODES",
        max_modes);

    if(max_modes <= 0 ||
       max_velocity_modes <= 0 ||
       max_pressure_modes <= 0)
    {
        ddrom_fail(
            "ST-DDROM velocity/pressure mode limits must be positive");
    }

    if(max_velocity_modes > INT32_MAX - max_pressure_modes){
        ddrom_fail("ST-DDROM total mode capacity overflow");
    }

    total_mode_capacity =
        max_velocity_modes + max_pressure_modes;

    if(ROM_std_stdd_system_initialize(
        &reduced,
        &fom.mono_com,
        num_windows,
        fom.fe.total_num_nodes,
        fom.num_owned_space_nodes,
        fom.window_dof,
        total_mode_capacity) != ROM_STDD_SUCCESS)
    {
        ddrom_fail("ROM_STDD_SYSTEM initialization failed");
    }

    if(ROM_std_stdd_fluid_context_initialize(
        &adapter,
        &fom,
        &layout,
        &reduced,
        ddrom_env_int(
            "FLUID_ST_DDROM_NR_MAX_ITER",
            fom.vals.nr_max_iter),
        ddrom_env_double(
            "FLUID_ST_DDROM_NR_EPSILON",
            fom.vals.nr_epsilon),
        ddrom_env_double(
            "FLUID_ST_DDROM_NR_ABS_EPSILON",
            1.0e-12),
        ddrom_env_int(
            "FLUID_ST_DDROM_LINEAR_MAX_ITER",
            fom.vals.mat_max_iter),
        ddrom_env_double(
            "FLUID_ST_DDROM_LINEAR_EPSILON",
            fom.vals.mat_epsilon))
       != ROM_STDD_SUCCESS)
    {
        ddrom_fail("fluid ST-DDROM adapter initialization failed");
    }

    mode = ddrom_read_mode();

    snapshot_directory = ddrom_env_string(
        "FLUID_ST_DDROM_SNAPSHOT_DIR",
        "./fluid_st_windows");

    available_snapshot_windows =
        ddrom_count_common_training_windows(
            snapshot_directory,
            num_windows,
            fom.mono_com.comm);

    requested_training_windows =
        ddrom_env_int(
            "FLUID_ST_DDROM_TRAINING_WINDOWS",
            0);

    if(requested_training_windows < 0){
        ddrom_fail(
            "FLUID_ST_DDROM_TRAINING_WINDOWS must be >= 0");
    }

    if(requested_training_windows == 0){
        num_training_windows =
            available_snapshot_windows;
    }
    else{
        num_training_windows =
            requested_training_windows;
    }

    if((mode == FLUID_ST_DDROM_OFFLINE ||
        mode == FLUID_ST_DDROM_FULL) &&
       (num_training_windows <= 0 ||
        num_training_windows > available_snapshot_windows))
    {
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(
                stderr,
                "ERROR: insufficient common FOM training windows.\n"
                "  requested training windows : %d\n"
                "  available on every rank    : %d\n"
                "  online/runtime windows     : %d\n"
                "Set FLUID_ST_DDROM_TRAINING_WINDOWS to a value "
                "between 1 and the available count, or regenerate "
                "the missing FOM window binaries.\n",
                num_training_windows,
                available_snapshot_windows,
                num_windows);
        }
        ddrom_fail(
            "invalid ST-DDROM training-window count");
    }

    online_windows = ddrom_env_int(
        "FLUID_ST_DDROM_ONLINE_WINDOWS",
        num_windows);

    if(online_windows <= 0 ||
       online_windows > num_windows)
    {
        ddrom_fail(
            "FLUID_ST_DDROM_ONLINE_WINDOWS must be in [1,num_windows]");
    }

    basis_directory = ddrom_env_string(
        "FLUID_ST_DDROM_BASIS_DIR",
        "./fluid_st_ddrom_basis");

    rom_window_directory = ddrom_env_string(
        "FLUID_ST_DDROM_WINDOW_DIR",
        "./fluid_st_ddrom_windows");

    coefficient_csv_name = ddrom_env_string(
        "FLUID_ST_DDROM_COEFFICIENT_CSV",
        "./fluid_st_ddrom_coefficients.csv");

    accuracy_csv_name = ddrom_env_string(
        "FLUID_ST_DDROM_ACCURACY_CSV",
        "./fluid_st_ddrom_accuracy.csv");

    validate_accuracy = ddrom_env_int(
        "FLUID_ST_DDROM_VALIDATE",
        1);

    if(validate_accuracy != 0 &&
       validate_accuracy != 1)
    {
        ddrom_fail(
            "FLUID_ST_DDROM_VALIDATE must be 0 or 1");
    }

    pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_POD_ENERGY",
        1.0);

    velocity_pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_VELOCITY_POD_ENERGY",
        pod_energy);

    pressure_pod_energy = ddrom_env_double(
        "FLUID_ST_DDROM_PRESSURE_POD_ENERGY",
        pod_energy);

    if(velocity_pod_energy <= 0.0 ||
       velocity_pod_energy > 1.0 ||
       pressure_pod_energy <= 0.0 ||
       pressure_pod_energy > 1.0)
    {
        ddrom_fail(
            "velocity/pressure POD energies must be in (0,1]");
    }

    write_output = ddrom_env_int(
        "FLUID_ST_DDROM_WRITE_OUTPUT",
        1);

    write_window_binary = ddrom_env_int(
        "FLUID_ST_DDROM_WRITE_WINDOW_BINARY",
        1);

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "\nfluid ST-DDROM configuration\n"
            "  mode                    : %s\n"
            "  dG degree               : %d\n"
            "  slabs/window            : %d\n"
            "  runtime window capacity : %d\n"
            "  online solve windows    : %d\n"
            "  FOM snapshots available : %d\n"
            "  POD training windows    : %d\n"
            "  physical MPI ranks      : %d\n"
            "  local/owned nodes rank0 : %d/%d\n"
            "  ST dof/node             : %d\n"
            "  basis structure         : separated velocity/pressure\n"
            "  max velocity modes/rank : %d\n"
            "  max pressure modes/rank : %d\n"
            "  total mode capacity/rank: %d\n"
            "  velocity POD energy     : %.15e\n"
            "  pressure POD energy     : %.15e\n"
            "  snapshot directory      : %s\n"
            "  basis directory         : %s\n"
            "  accuracy validation     : %s\n",
            mode == FLUID_ST_DDROM_OFFLINE
                ? "offline"
                : (mode == FLUID_ST_DDROM_ONLINE
                    ? "online"
                    : "full"),
            layout.time.degree,
            layout.num_slabs,
            num_windows,
            online_windows,
            available_snapshot_windows,
            num_training_windows,
            monolis_mpi_get_global_comm_size(),
            fom.fe.total_num_nodes,
            fom.num_owned_space_nodes,
            fom.window_dof,
            max_velocity_modes,
            max_pressure_modes,
            total_mode_capacity,
            velocity_pod_energy,
            pressure_pod_energy,
            snapshot_directory,
            basis_directory,
            validate_accuracy ? "enabled" : "disabled");
    }

    if(mode == FLUID_ST_DDROM_OFFLINE ||
       mode == FLUID_ST_DDROM_FULL)
    {
        if(ddrom_build_snapshot_matrix(
            &adapter,
            snapshot_directory,
            num_training_windows,
            &snapshot_matrix) != 0)
        {
            ddrom_fail(
                "failed to build the homogeneous ST snapshot matrix");
        }

        if(ddrom_compute_separated_velocity_pressure_pod(
            &adapter,
            snapshot_matrix,
            num_training_windows,
            max_velocity_modes,
            max_pressure_modes,
            velocity_pod_energy,
            pressure_pod_energy,
            fom.mono_com.comm,
            &velocity_modes,
            &pressure_modes) != ROM_STDD_SUCCESS)
        {
            ddrom_fail(
                "separated velocity/pressure space-time POD failed");
        }

        if(ROM_std_stdd_fluid_enforce_homogeneous_basis(
            &adapter) != ROM_STDD_SUCCESS)
        {
            ddrom_fail(
                "homogeneous fluid ST basis enforcement failed");
        }

        if(ROM_std_stdd_write_basis_rank(
            &reduced,
            num_training_windows,
            basis_directory) != ROM_STDD_SUCCESS)
        {
            ddrom_fail("ST-DDROM basis write failed");
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "fluid ST-DDROM offline completed: "
                "velocity modes/rank=%d, pressure modes/rank=%d, "
                "total modes/rank=%d, training windows=%d "
                "(runtime windows=%d)\n",
                velocity_modes,
                pressure_modes,
                reduced.num_modes,
                num_training_windows,
                num_windows);
        }
    }

    if(mode == FLUID_ST_DDROM_ONLINE ||
       mode == FLUID_ST_DDROM_FULL)
    {
        if(mode == FLUID_ST_DDROM_ONLINE){
            if(ROM_std_stdd_read_basis_rank(
                &reduced,
                basis_directory) != ROM_STDD_SUCCESS)
            {
                ddrom_fail("ST-DDROM basis read failed");
            }

            if(ROM_std_stdd_fluid_enforce_homogeneous_basis(
                &adapter) != ROM_STDD_SUCCESS)
            {
                ddrom_fail(
                    "online homogeneous basis enforcement failed");
            }
        }

        if(ddrom_classify_split_basis_modes(
               &reduced,
               &velocity_modes,
               &pressure_modes) != 0)
        {
            ddrom_fail(
                "loaded basis is not a separated velocity/pressure basis; "
                "rerun ST-DDROM offline with the split-basis implementation");
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "fluid ST-DDROM loaded split basis: "
                "velocity modes/rank=%d pressure modes/rank=%d "
                "total=%d\n",
                velocity_modes,
                pressure_modes,
                reduced.num_modes);
        }

        {
            char rank_csv[4096];
            const int rank =
                monolis_mpi_get_global_my_rank();
            int written = snprintf(
                rank_csv,
                sizeof(rank_csv),
                "%s.rank_%06d",
                coefficient_csv_name,
                rank);

            if(written < 0 ||
               (size_t)written >= sizeof(rank_csv))
            {
                ddrom_fail("coefficient CSV path is too long");
            }

            coefficient_csv = fopen(rank_csv, "w");
            if(coefficient_csv == NULL){
                ddrom_fail("cannot open coefficient CSV");
            }

            ddrom_write_coefficient_header(
                coefficient_csv,
                reduced.num_modes);
        }

        reference_owned = (double*)calloc(
            (size_t)reduced.internal_st_dof,
            sizeof(double));

        if(reference_owned == NULL){
            ddrom_fail("reference-window allocation failed");
        }

        if(validate_accuracy && available_snapshot_windows > 0){
            const size_t local_st_dof =
                (size_t)reduced.local_st_dof;
            const size_t internal_st_dof =
                (size_t)reduced.internal_st_dof;
            const size_t trace_dof =
                (size_t)fom.fe.total_num_nodes
                * (size_t)FLUID_SUPS_PHYSICAL_DOF;

            reference_full =
                (double*)calloc(
                    local_st_dof,
                    sizeof(double));
            reference_homogeneous =
                (double*)calloc(
                    internal_st_dof,
                    sizeof(double));
            projected_homogeneous =
                (double*)calloc(
                    internal_st_dof,
                    sizeof(double));
            reference_coefficients =
                (double*)calloc(
                    (size_t)reduced.num_modes_capacity,
                    sizeof(double));
            reference_previous_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            reference_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            rom_trace =
                (double*)calloc(
                    trace_dof,
                    sizeof(double));
            reference_projected_rhs =
                (double*)calloc(
                    (size_t)reduced.num_modes_capacity,
                    sizeof(double));

            if(reference_full == NULL ||
               reference_homogeneous == NULL ||
               projected_homogeneous == NULL ||
               reference_coefficients == NULL ||
               reference_previous_trace == NULL ||
               reference_trace == NULL ||
               rom_trace == NULL ||
               reference_projected_rhs == NULL)
            {
                ddrom_fail(
                    "ST-DDROM accuracy-validation allocation failed");
            }

            memcpy(
                reference_previous_trace,
                adapter.initial_trace,
                trace_dof * sizeof(double));

            if(monolis_mpi_get_global_my_rank() == 0){
                accuracy_csv = fopen(
                    accuracy_csv_name,
                    "w");
                if(accuracy_csv == NULL){
                    ddrom_fail(
                        "cannot open ST-DDROM accuracy CSV");
                }
                ddrom_write_accuracy_header(
                    accuracy_csv);
            }
        }

        for(int window = 0;
            window < online_windows;
            window++)
        {
            ROM_STDD_FLUID_WINDOW_REPORT report;
            double relative_reference_error = NAN;
            const double* coefficient;

            if(ROM_std_stdd_fluid_solve_window(
                &adapter,
                window,
                &report) != ROM_STDD_SUCCESS)
            {
                fprintf(
                    stderr,
                    "ERROR: fluid ST-DDROM failed in window %d\n",
                    window + 1);
                goto cleanup;
            }

            coefficient =
                reduced.reduced_coef
                + (size_t)window
                    * (size_t)reduced.num_modes_capacity;

            fprintf(coefficient_csv, "%d", window);
            for(int mode_id = 0;
                mode_id < reduced.num_modes;
                mode_id++)
            {
                fprintf(
                    coefficient_csv,
                    ",%.17e",
                    coefficient[mode_id]);
            }
            fprintf(
                coefficient_csv,
                ",%.17e,%.17e,%.17e\n",
                report.relative_full_residual,
                report.final_projected_residual,
                report.final_coefficient_correction);
            fflush(coefficient_csv);

            /*
             * Detailed FOM/ROM accuracy diagnostics are evaluated whenever a
             * matching reference window exists.  The reference uses its own
             * previous FOM trace for the residual check.
             */
            if(window < available_snapshot_windows &&
               ddrom_read_owned_window(
                   snapshot_directory,
                   window,
                   &layout,
                   reference_owned) == 0)
            {
                ROM_std_stdd_fluid_owned_window_error(
                    &adapter,
                    reference_owned,
                    &relative_reference_error);

                if(validate_accuracy){
                    DDROM_WINDOW_ACCURACY window_accuracy;

                    if(ddrom_validate_current_window(
                           &adapter,
                           window,
                           reference_owned,
                           reference_full,
                           reference_homogeneous,
                           projected_homogeneous,
                           reference_coefficients,
                           reference_previous_trace,
                           reference_trace,
                           rom_trace,
                           reference_projected_rhs,
                           &window_accuracy,
                           &accuracy_summary) != 0)
                    {
                        ddrom_fail(
                            "detailed ST-DDROM accuracy validation failed");
                    }

                    if(monolis_mpi_get_global_my_rank() == 0){
                        printf(
                            "fluid ST-DDROM accuracy window %d: "
                            "projection=%e state=%e dynamic=%e "
                            "velocity(state)=%e pressure(state)=%e "
                            "coef=%e trace=%e "
                            "FOM residual=%e FOM projected defect=%e\n",
                            window + 1,
                            window_accuracy.projection_error,
                            window_accuracy.state_error,
                            window_accuracy.dynamic_error,
                            window_accuracy.velocity_state_error,
                            window_accuracy.pressure_state_error,
                            window_accuracy.coefficient_error,
                            window_accuracy.trace_error,
                            window_accuracy.reference_full_residual,
                            window_accuracy.reference_projected_defect);

                        if(accuracy_csv != NULL){
                            fprintf(
                                accuracy_csv,
                                "%d,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,"
                                "%.17e,%.17e,%.17e,"
                                "%.17e,%.17e,"
                                "%.17e,%.17e,"
                                "%.17e,%.17e,%.17e,%.17e\n",
                                window,
                                window_accuracy.projection_error,
                                window_accuracy.state_error,
                                window_accuracy.dynamic_error,
                                window_accuracy.velocity_projection_error,
                                window_accuracy.velocity_state_error,
                                window_accuracy.velocity_dynamic_error,
                                window_accuracy.pressure_projection_error,
                                window_accuracy.pressure_state_error,
                                window_accuracy.pressure_dynamic_error,
                                window_accuracy.coefficient_error,
                                window_accuracy.trace_error,
                                window_accuracy.trace_velocity_error,
                                window_accuracy.trace_pressure_error,
                                report.final_projected_residual,
                                report.relative_projected_residual,
                                report.final_full_residual,
                                report.relative_full_residual,
                                window_accuracy.reference_full_residual,
                                window_accuracy.reference_momentum_residual,
                                window_accuracy.reference_mass_residual,
                                window_accuracy.reference_projected_defect);
                            fflush(accuracy_csv);
                        }
                    }
                }
            }

            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "fluid ST-DDROM window %d completed: "
                    "Newton=%d projected rel residual=%e "
                    "full rel residual=%e "
                    "FOM window error=%e\n",
                    window + 1,
                    report.nonlinear_iterations,
                    report.relative_projected_residual,
                    report.relative_full_residual,
                    relative_reference_error);
            }

            if(write_output){
                if(ROM_std_stdd_fluid_write_current_window_vtk(
                    &adapter,
                    window,
                    &file_number) != ROM_STDD_SUCCESS)
                {
                    ddrom_fail("ROM VTK output failed");
                }
            }

            if(write_window_binary){
                if(ROM_std_stdd_fluid_write_current_window_binary(
                    &adapter,
                    window,
                    rom_window_directory) != ROM_STDD_SUCCESS)
                {
                    ddrom_fail("ROM window binary output failed");
                }
            }
        }

        if(validate_accuracy &&
           accuracy_summary.num_validated_windows > 0)
        {
            ddrom_print_accuracy_summary(
                &accuracy_summary);

            if(monolis_mpi_get_global_my_rank() == 0 &&
               accuracy_summary.num_validated_windows <
                   online_windows)
            {
                printf(
                    "NOTE: accuracy was evaluated for %d/%d windows because "
                    "only %d common FOM reference windows are available.\n",
                    accuracy_summary.num_validated_windows,
                    online_windows,
                    available_snapshot_windows);
            }
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\nfluid ST-DDROM online completed successfully.\n");
        }
    }

    status = EXIT_SUCCESS;

cleanup:
    if(coefficient_csv != NULL){
        fclose(coefficient_csv);
    }
    if(accuracy_csv != NULL){
        fclose(accuracy_csv);
    }

    free(snapshot_matrix);
    free(reference_owned);
    free(reference_full);
    free(reference_homogeneous);
    free(projected_homogeneous);
    free(reference_coefficients);
    free(reference_previous_trace);
    free(reference_trace);
    free(rom_trace);
    free(reference_projected_rhs);

    ROM_std_stdd_fluid_context_finalize(&adapter);
    ROM_std_stdd_system_finalize(&reduced);
    fluid_sups_st_finalize(&fom);
    monolis_global_finalize();

    return status;
}
