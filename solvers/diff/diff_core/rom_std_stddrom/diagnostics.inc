static double stdd_system_reduced_residual_from_coefficients_impl(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    const double* coefficients,
    double* window_relative_residuals)
{
    MONOLIS_COM reduced_com;
    double* current_coefficients = NULL;
    double* previous_coefficients = NULL;
    double* reduced_values = NULL;
    double global_residual2 = 0.0;
    double global_rhs2 = 0.0;
    int block_slots;
    size_t block_size2;
    double result = INFINITY;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       spatial_com == NULL || coefficients == NULL ||
       system->diagonal_blocks == NULL ||
       system->temporal_blocks == NULL || system->reduced_rhs == NULL ||
       system->num_neighbor_ranks != spatial_com->recv_n_neib)
    {
        return INFINITY;
    }

    memset(&reduced_com, 0, sizeof(reduced_com));
    if(stdd_allocate_reduced_comm(&reduced_com, spatial_com, 1)
       != ROM_STDD_SUCCESS)
    {
        return INFINITY;
    }

    block_slots = 1 + system->num_neighbor_ranks;
    block_size2 =
        (size_t)system->num_modes * (size_t)system->num_modes;

    current_coefficients = (double*)stdd_calloc(
        (size_t)block_slots * (size_t)system->num_modes,
        sizeof(double));
    previous_coefficients = (double*)stdd_calloc(
        (size_t)block_slots * (size_t)system->num_modes,
        sizeof(double));
    reduced_values = (double*)stdd_calloc(
        (size_t)(2 * system->num_windows),
        sizeof(double));

    if(current_coefficients == NULL || previous_coefficients == NULL ||
       reduced_values == NULL)
    {
        goto cleanup;
    }

    for(int window = 0; window < system->num_windows; window++){
        const double* coefficient = coefficients
            + (size_t)window * (size_t)system->num_modes_capacity;
        const double* rhs = system->reduced_rhs
            + (size_t)window * (size_t)system->num_modes_capacity;
        double local_residual2 = 0.0;
        double local_rhs2 = 0.0;

        memset(
            current_coefficients,
            0,
            (size_t)block_slots * (size_t)system->num_modes
                * sizeof(double));
        memcpy(
            current_coefficients,
            coefficient,
            (size_t)system->num_modes * sizeof(double));
        monolis_mpi_update_R(
            &reduced_com,
            block_slots,
            system->num_modes,
            current_coefficients);

        memset(
            previous_coefficients,
            0,
            (size_t)block_slots * (size_t)system->num_modes
                * sizeof(double));
        if(window > 0){
            memcpy(
                previous_coefficients,
                coefficients
                    + (size_t)(window - 1)
                        * (size_t)system->num_modes_capacity,
                (size_t)system->num_modes * sizeof(double));
            monolis_mpi_update_R(
                &reduced_com,
                block_slots,
                system->num_modes,
                previous_coefficients);
        }

        for(int row_mode = 0; row_mode < system->num_modes; row_mode++){
            double lhs = 0.0;

            for(int slot = 0; slot < block_slots; slot++){
                const double* diagonal_block = system->diagonal_blocks
                    + (size_t)slot * block_size2;
                const double* current_slot = current_coefficients
                    + (size_t)slot * (size_t)system->num_modes;

                for(int column_mode = 0;
                    column_mode < system->num_modes;
                    column_mode++)
                {
                    lhs += diagonal_block[
                        (size_t)row_mode * (size_t)system->num_modes
                        + (size_t)column_mode]
                        * current_slot[column_mode];
                }

                if(window > 0){
                    const double* temporal_block = system->temporal_blocks
                        + (size_t)slot * block_size2;
                    const double* previous_slot = previous_coefficients
                        + (size_t)slot * (size_t)system->num_modes;

                    for(int column_mode = 0;
                        column_mode < system->num_modes;
                        column_mode++)
                    {
                        lhs -= temporal_block[
                            (size_t)row_mode * (size_t)system->num_modes
                            + (size_t)column_mode]
                            * previous_slot[column_mode];
                    }
                }
            }

            {
                const double difference = lhs - rhs[row_mode];
                local_residual2 += difference * difference;
                local_rhs2 += rhs[row_mode] * rhs[row_mode];
            }
        }

        reduced_values[2 * window] = local_residual2;
        reduced_values[2 * window + 1] = local_rhs2;
    }

    monolis_allreduce_R(
        2 * system->num_windows,
        reduced_values,
        MONOLIS_MPI_SUM,
        reduced_com.comm);

    for(int window = 0; window < system->num_windows; window++){
        const double residual2 = reduced_values[2 * window];
        const double rhs2 = reduced_values[2 * window + 1];

        global_residual2 += residual2;
        global_rhs2 += rhs2;

        if(window_relative_residuals != NULL){
            window_relative_residuals[window] =
                sqrt(residual2) / fmax(sqrt(rhs2), 1.0e-300);
        }
    }

    result = sqrt(global_residual2) /
        fmax(sqrt(global_rhs2), 1.0e-300);

cleanup:
    stdd_free_reduced_comm(&reduced_com);
    free(current_coefficients);
    free(previous_coefficients);
    free(reduced_values);
    return result;
}

double ROM_std_stdd_system_reduced_residual_from_coefficients(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    const double* coefficients,
    double* window_relative_residuals)
{
    return stdd_system_reduced_residual_from_coefficients_impl(
        system,
        spatial_com,
        coefficients,
        window_relative_residuals);
}

double ROM_std_stdd_system_reduced_residual(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com)
{
    if(system == NULL){
        return INFINITY;
    }

    return stdd_system_reduced_residual_from_coefficients_impl(
        system,
        spatial_com,
        system->reduced_coef,
        NULL);
}

int ROM_std_stdd_evaluate_local_reduced_equation(
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    const double* coefficients,
    int window_id,
    double* defect,
    double* self_action,
    double* neighbor_action,
    double* temporal_action)
{
    MONOLIS_COM reduced_com;
    double* current_coefficients = NULL;
    double* previous_coefficients = NULL;
    const double* rhs;
    int block_slots;
    size_t block_size2;
    int status = ROM_STDD_SUCCESS;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       spatial_com == NULL || coefficients == NULL ||
       window_id < 0 || window_id >= system->num_windows ||
       defect == NULL || self_action == NULL ||
       neighbor_action == NULL || temporal_action == NULL ||
       system->diagonal_blocks == NULL ||
       system->temporal_blocks == NULL ||
       system->reduced_rhs == NULL ||
       system->num_neighbor_ranks != spatial_com->recv_n_neib)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(&reduced_com, 0, sizeof(reduced_com));

    if(stdd_allocate_reduced_comm(
        &reduced_com,
        spatial_com,
        1) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_ALLOCATION;
    }

    block_slots = 1 + system->num_neighbor_ranks;
    block_size2 =
        (size_t)system->num_modes * (size_t)system->num_modes;

    current_coefficients = (double*)stdd_calloc(
        (size_t)block_slots * (size_t)system->num_modes,
        sizeof(double));
    previous_coefficients = (double*)stdd_calloc(
        (size_t)block_slots * (size_t)system->num_modes,
        sizeof(double));

    if(current_coefficients == NULL || previous_coefficients == NULL){
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    memcpy(
        current_coefficients,
        coefficients
            + (size_t)window_id
                * (size_t)system->num_modes_capacity,
        (size_t)system->num_modes * sizeof(double));

    monolis_mpi_update_R(
        &reduced_com,
        block_slots,
        system->num_modes,
        current_coefficients);

    if(window_id > 0){
        memcpy(
            previous_coefficients,
            coefficients
                + (size_t)(window_id - 1)
                    * (size_t)system->num_modes_capacity,
            (size_t)system->num_modes * sizeof(double));

        monolis_mpi_update_R(
            &reduced_com,
            block_slots,
            system->num_modes,
            previous_coefficients);
    }

    rhs = system->reduced_rhs
        + (size_t)window_id * (size_t)system->num_modes_capacity;

    for(int row_mode = 0; row_mode < system->num_modes; row_mode++){
        double self_value = 0.0;
        double neighbor_value = 0.0;
        double temporal_value = 0.0;

        {
            const double* block = system->diagonal_blocks;
            const double* coefficient = current_coefficients;

            for(int column_mode = 0;
                column_mode < system->num_modes;
                column_mode++)
            {
                self_value += block[
                    (size_t)row_mode * (size_t)system->num_modes
                    + (size_t)column_mode]
                    * coefficient[column_mode];
            }
        }

        for(int slot = 1; slot < block_slots; slot++){
            const double* block = system->diagonal_blocks
                + (size_t)slot * block_size2;
            const double* coefficient = current_coefficients
                + (size_t)slot * (size_t)system->num_modes;

            for(int column_mode = 0;
                column_mode < system->num_modes;
                column_mode++)
            {
                neighbor_value += block[
                    (size_t)row_mode * (size_t)system->num_modes
                    + (size_t)column_mode]
                    * coefficient[column_mode];
            }
        }

        if(window_id > 0){
            for(int slot = 0; slot < block_slots; slot++){
                const double* block = system->temporal_blocks
                    + (size_t)slot * block_size2;
                const double* coefficient = previous_coefficients
                    + (size_t)slot * (size_t)system->num_modes;

                for(int column_mode = 0;
                    column_mode < system->num_modes;
                    column_mode++)
                {
                    temporal_value += block[
                        (size_t)row_mode * (size_t)system->num_modes
                        + (size_t)column_mode]
                        * coefficient[column_mode];
                }
            }
        }

        self_action[row_mode] = self_value;
        neighbor_action[row_mode] = neighbor_value;
        temporal_action[row_mode] = temporal_value;
        defect[row_mode] =
            self_value
            + neighbor_value
            - temporal_value
            - rhs[row_mode];
    }

cleanup:
    stdd_free_reduced_comm(&reduced_com);
    free(current_coefficients);
    free(previous_coefficients);

    return status;
}

int ROM_std_stdd_write_reduced_graph_rank(
    const ROM_STDD_REDUCED_SOLVER* solver,
    const char* directory)
{
    char filename[4096];
    FILE* fp;
    int rank;

    if(solver == NULL || solver->graph_index == NULL ||
       solver->graph_item == NULL || directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_make_directory(directory) != ROM_STDD_SUCCESS){
        return ROM_STDD_ERR_IO;
    }

    rank = monolis_mpi_get_global_my_rank();
    if(snprintf(
        filename,
        sizeof(filename),
        "%s/stddrom_metagraph.%d",
        directory,
        rank) >= (int)sizeof(filename))
    {
        return ROM_STDD_ERR_IO;
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        return ROM_STDD_ERR_IO;
    }

    fprintf(fp, "%d\n", solver->num_total_meta_nodes);
    for(int row = 0; row < solver->num_total_meta_nodes; row++){
        const int start = solver->graph_index[row];
        const int end = solver->graph_index[row + 1];
        fprintf(fp, "%d %d", row, end - start);
        for(int p = start; p < end; p++){
            fprintf(fp, " %d", solver->graph_item[p]);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

void ROM_std_stdd_reduced_solver_finalize(
    ROM_STDD_REDUCED_SOLVER* solver)
{
    if(solver == NULL){
        return;
    }

    monolis_finalize(&solver->monolis);
    stdd_free_reduced_comm(&solver->comm);
    free(solver->graph_index);
    free(solver->graph_item);
    free(solver->solution);

    memset(solver, 0, sizeof(*solver));
}

double ROM_std_stdd_reduced_algebraic_residual(
    ROM_STDD_REDUCED_SOLVER* solver)
{
    double* product;
    double local_residual2 = 0.0;
    double local_rhs2 = 0.0;
    double global_residual2 = 0.0;
    double global_rhs2 = 0.0;
    int local_dof;

    if(solver == NULL || solver->solution == NULL ||
       solver->monolis.mat.R.B == NULL)
    {
        return INFINITY;
    }

    local_dof = solver->num_total_meta_nodes * solver->block_size;
    product = (double*)stdd_calloc((size_t)local_dof, sizeof(double));
    if(product == NULL){
        return INFINITY;
    }

    monolis_matvec_product_R(
        &solver->monolis,
        &solver->comm,
        solver->solution,
        product);

    for(int i = 0;
        i < solver->num_internal_meta_nodes * solver->block_size;
        i++)
    {
        const double difference =
            product[i] - solver->monolis.mat.R.B[i];
        local_residual2 += difference * difference;
        local_rhs2 += solver->monolis.mat.R.B[i]
            * solver->monolis.mat.R.B[i];
    }

    global_residual2 = local_residual2;
    monolis_allreduce_R(
        1,
        &global_residual2,
        MONOLIS_MPI_SUM,
        solver->comm.comm);

    global_rhs2 = local_rhs2;
    monolis_allreduce_R(
        1,
        &global_rhs2,
        MONOLIS_MPI_SUM,
        solver->comm.comm);

    free(product);

    return sqrt(global_residual2) /
        fmax(sqrt(global_rhs2), 1.0e-300);
}

double ROM_std_stdd_global_relative_error(
    const double* value,
    const double* reference,
    int local_length,
    int communicator)
{
    double local_difference2 = 0.0;
    double local_reference2 = 0.0;
    double global_difference2 = 0.0;
    double global_reference2 = 0.0;

    if(value == NULL || reference == NULL || local_length < 0){
        return INFINITY;
    }

    for(int i = 0; i < local_length; i++){
        const double difference = value[i] - reference[i];
        local_difference2 += difference * difference;
        local_reference2 += reference[i] * reference[i];
    }

    global_difference2 = local_difference2;
    monolis_allreduce_R(
        1,
        &global_difference2,
        MONOLIS_MPI_SUM,
        communicator);

    global_reference2 = local_reference2;
    monolis_allreduce_R(
        1,
        &global_reference2,
        MONOLIS_MPI_SUM,
        communicator);

    return sqrt(global_difference2) /
        fmax(sqrt(global_reference2), 1.0e-300);
}

int ROM_std_stdd_project_reference_coefficients(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* reference_homogeneous,
    double* coefficients)
{
    if(system == NULL || reference_homogeneous == NULL ||
       coefficients == NULL || system->basis == NULL ||
       window_id < 0 || window_id >= system->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int mode = 0; mode < system->num_modes; mode++){
        double value = 0.0;
        for(int row = 0; row < system->internal_st_dof; row++){
            value += system->basis[
                stdd_basis_offset(system, row, mode)]
                * reference_homogeneous[row];
        }
        coefficients[mode] = value;
    }

    return ROM_STDD_SUCCESS;
}
