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

