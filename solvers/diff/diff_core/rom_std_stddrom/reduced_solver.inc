static int stdd_allocate_reduced_comm(
    MONOLIS_COM* reduced_com,
    const MONOLIS_COM* spatial_com,
    int num_windows)
{
    int recv_items;
    int send_items;

    if(reduced_com == NULL || spatial_com == NULL || num_windows <= 0){
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(reduced_com, 0, sizeof(*reduced_com));

    reduced_com->comm = monolis_mpi_get_global_comm();
    reduced_com->comm_size = monolis_mpi_get_global_comm_size();
    reduced_com->my_rank = monolis_mpi_get_global_my_rank();
    reduced_com->n_internal_vertex = num_windows;
    reduced_com->recv_n_neib = spatial_com->recv_n_neib;
    reduced_com->send_n_neib = spatial_com->send_n_neib;

    recv_items = num_windows * reduced_com->recv_n_neib;
    send_items = num_windows * reduced_com->send_n_neib;

    reduced_com->recv_index = (int*)stdd_calloc(
        (size_t)reduced_com->recv_n_neib + 1,
        sizeof(int));
    reduced_com->send_index = (int*)stdd_calloc(
        (size_t)reduced_com->send_n_neib + 1,
        sizeof(int));
    reduced_com->recv_item = (int*)stdd_calloc(
        (size_t)(recv_items > 0 ? recv_items : 1),
        sizeof(int));
    reduced_com->send_item = (int*)stdd_calloc(
        (size_t)(send_items > 0 ? send_items : 1),
        sizeof(int));
    reduced_com->recv_neib_pe = (int*)stdd_calloc(
        (size_t)(reduced_com->recv_n_neib > 0
            ? reduced_com->recv_n_neib : 1),
        sizeof(int));
    reduced_com->send_neib_pe = (int*)stdd_calloc(
        (size_t)(reduced_com->send_n_neib > 0
            ? reduced_com->send_n_neib : 1),
        sizeof(int));

    if(reduced_com->recv_index == NULL || reduced_com->send_index == NULL ||
       reduced_com->recv_item == NULL || reduced_com->send_item == NULL ||
       reduced_com->recv_neib_pe == NULL || reduced_com->send_neib_pe == NULL)
    {
        return ROM_STDD_ERR_ALLOCATION;
    }

    if(reduced_com->recv_n_neib > 0){
        memcpy(
            reduced_com->recv_neib_pe,
            spatial_com->recv_neib_pe,
            (size_t)reduced_com->recv_n_neib * sizeof(int));
    }
    if(reduced_com->send_n_neib > 0){
        memcpy(
            reduced_com->send_neib_pe,
            spatial_com->send_neib_pe,
            (size_t)reduced_com->send_n_neib * sizeof(int));
    }

    for(int n = 0; n < reduced_com->recv_n_neib; n++){
        reduced_com->recv_index[n + 1] = (n + 1) * num_windows;
        for(int w = 0; w < num_windows; w++){
            reduced_com->recv_item[n * num_windows + w] =
                num_windows + n * num_windows + w;
        }
    }

    for(int n = 0; n < reduced_com->send_n_neib; n++){
        reduced_com->send_index[n + 1] = (n + 1) * num_windows;
        for(int w = 0; w < num_windows; w++){
            reduced_com->send_item[n * num_windows + w] = w;
        }
    }

    return ROM_STDD_SUCCESS;
}

static void stdd_free_reduced_comm(MONOLIS_COM* reduced_com)
{
    if(reduced_com == NULL){
        return;
    }

    free(reduced_com->recv_index);
    free(reduced_com->recv_item);
    free(reduced_com->send_index);
    free(reduced_com->send_item);
    free(reduced_com->recv_neib_pe);
    free(reduced_com->send_neib_pe);

    memset(reduced_com, 0, sizeof(*reduced_com));
}

static int stdd_halo_meta_node(
    const ROM_STDD_REDUCED_SOLVER* solver,
    int neighbor_slot,
    int window_id)
{
    return solver->num_windows
        + neighbor_slot * solver->num_windows
        + window_id;
}

static int stdd_dense_lu_factor(
    double* matrix,
    int dimension,
    int* pivot)
{
    double maximum_entry = 0.0;
    double pivot_tolerance;

    if(matrix == NULL || pivot == NULL || dimension <= 0){
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int i = 0; i < dimension * dimension; i++){
        const double value = fabs(matrix[i]);
        if(!isfinite(value)){
            return ROM_STDD_ERR_DIMENSION;
        }
        if(value > maximum_entry){
            maximum_entry = value;
        }
    }

    pivot_tolerance = 1.0e-13 * fmax(maximum_entry, 1.0);

    for(int column = 0; column < dimension; column++){
        int pivot_row = column;
        double pivot_value = fabs(
            matrix[(size_t)column * (size_t)dimension + (size_t)column]);

        for(int row = column + 1; row < dimension; row++){
            const double candidate = fabs(
                matrix[(size_t)row * (size_t)dimension + (size_t)column]);
            if(candidate > pivot_value){
                pivot_value = candidate;
                pivot_row = row;
            }
        }

        if(!isfinite(pivot_value) || pivot_value <= pivot_tolerance){
            fprintf(
                stderr,
                "ERROR [rank %d]: singular/ill-conditioned local reduced "
                "self block at column %d: pivot=%.15e tolerance=%.15e\n",
                monolis_mpi_get_global_my_rank(),
                column,
                pivot_value,
                pivot_tolerance);
            return ROM_STDD_ERR_DIMENSION;
        }

        pivot[column] = pivot_row;

        if(pivot_row != column){
            for(int j = 0; j < dimension; j++){
                const size_t first =
                    (size_t)column * (size_t)dimension + (size_t)j;
                const size_t second =
                    (size_t)pivot_row * (size_t)dimension + (size_t)j;
                const double temporary = matrix[first];
                matrix[first] = matrix[second];
                matrix[second] = temporary;
            }
        }

        for(int row = column + 1; row < dimension; row++){
            const size_t multiplier_index =
                (size_t)row * (size_t)dimension + (size_t)column;
            matrix[multiplier_index] /=
                matrix[(size_t)column * (size_t)dimension + (size_t)column];

            for(int j = column + 1; j < dimension; j++){
                matrix[(size_t)row * (size_t)dimension + (size_t)j] -=
                    matrix[multiplier_index]
                    * matrix[(size_t)column * (size_t)dimension + (size_t)j];
            }
        }
    }

    return ROM_STDD_SUCCESS;
}

static int stdd_dense_lu_solve(
    const double* lu,
    const int* pivot,
    int dimension,
    double* vector)
{
    if(lu == NULL || pivot == NULL || vector == NULL || dimension <= 0){
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int column = 0; column < dimension; column++){
        if(pivot[column] != column){
            const double temporary = vector[column];
            vector[column] = vector[pivot[column]];
            vector[pivot[column]] = temporary;
        }
    }

    for(int row = 0; row < dimension; row++){
        for(int column = 0; column < row; column++){
            vector[row] -=
                lu[(size_t)row * (size_t)dimension + (size_t)column]
                * vector[column];
        }
    }

    for(int row = dimension - 1; row >= 0; row--){
        for(int column = row + 1; column < dimension; column++){
            vector[row] -=
                lu[(size_t)row * (size_t)dimension + (size_t)column]
                * vector[column];
        }
        vector[row] /=
            lu[(size_t)row * (size_t)dimension + (size_t)row];

        if(!isfinite(vector[row])){
            return ROM_STDD_ERR_DIMENSION;
        }
    }

    return ROM_STDD_SUCCESS;
}

static int stdd_left_scale_dense_block(
    const double* lu,
    const int* pivot,
    int dimension,
    const double* source,
    double* destination,
    double* work)
{
    if(lu == NULL || pivot == NULL || source == NULL ||
       destination == NULL || work == NULL || dimension <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int column = 0; column < dimension; column++){
        for(int row = 0; row < dimension; row++){
            work[row] = source[
                (size_t)row * (size_t)dimension + (size_t)column];
        }

        if(stdd_dense_lu_solve(lu, pivot, dimension, work)
           != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_DIMENSION;
        }

        for(int row = 0; row < dimension; row++){
            destination[
                (size_t)row * (size_t)dimension + (size_t)column] =
                work[row];
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_reduced_solver_initialize(
    ROM_STDD_REDUCED_SOLVER* solver,
    const ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int max_iter,
    double tolerance)
{
    int position = 0;
    int status;

    if(solver == NULL || system == NULL || spatial_com == NULL ||
       system->num_modes <= 0 || max_iter <= 0 ||
       !isfinite(tolerance) || tolerance <= 0.0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(solver, 0, sizeof(*solver));

    solver->num_windows = system->num_windows;
    solver->block_size = system->num_modes;
    solver->num_internal_meta_nodes = system->num_windows;
    solver->num_neighbor_ranks = spatial_com->recv_n_neib;
    solver->num_total_meta_nodes = system->num_windows
        * (1 + spatial_com->recv_n_neib);
    solver->max_iter = max_iter;
    solver->tolerance = tolerance;

    status = stdd_allocate_reduced_comm(
        &solver->comm,
        spatial_com,
        system->num_windows);
    if(status != ROM_STDD_SUCCESS){
        ROM_std_stdd_reduced_solver_finalize(solver);
        return status;
    }

    solver->graph_index = (int*)stdd_calloc(
        (size_t)solver->num_total_meta_nodes + 1,
        sizeof(int));

    /*
     * Internal rows:
     *   self/current + neighbour/current
     *   and, for w>0, self/previous + neighbour/previous.
     * Halo rows retain a self entry only; they are communication columns,
     * not owned equations.
     */
    solver->num_graph_items =
        (2 * system->num_windows - 1)
            * (1 + solver->num_neighbor_ranks)
        + system->num_windows * solver->num_neighbor_ranks;

    solver->graph_item = (int*)stdd_calloc(
        (size_t)solver->num_graph_items,
        sizeof(int));

    if(solver->graph_index == NULL || solver->graph_item == NULL){
        ROM_std_stdd_reduced_solver_finalize(solver);
        return ROM_STDD_ERR_ALLOCATION;
    }

    solver->graph_index[0] = 0;

    for(int w = 0; w < system->num_windows; w++){
        if(w > 0){
            solver->graph_item[position++] = w - 1;
            for(int n = 0; n < solver->num_neighbor_ranks; n++){
                solver->graph_item[position++] =
                    stdd_halo_meta_node(solver, n, w - 1);
            }
        }

        solver->graph_item[position++] = w;
        for(int n = 0; n < solver->num_neighbor_ranks; n++){
            solver->graph_item[position++] =
                stdd_halo_meta_node(solver, n, w);
        }

        solver->graph_index[w + 1] = position;
    }

    for(int node = system->num_windows;
        node < solver->num_total_meta_nodes;
        node++)
    {
        solver->graph_item[position++] = node;
        solver->graph_index[node + 1] = position;
    }

    if(position != solver->num_graph_items){
        ROM_std_stdd_reduced_solver_finalize(solver);
        return ROM_STDD_ERR_DIMENSION;
    }

    monolis_initialize(&solver->monolis);
    monolis_get_nonzero_pattern_by_nodal_graph_R(
        &solver->monolis,
        solver->num_total_meta_nodes,
        solver->block_size,
        solver->graph_index,
        solver->graph_item);

    solver->monolis.mat.N = solver->num_internal_meta_nodes;
    solver->monolis.mat.NP = solver->num_total_meta_nodes;

    solver->solution = (double*)stdd_calloc(
        (size_t)solver->num_total_meta_nodes
            * (size_t)solver->block_size,
        sizeof(double));

    if(solver->solution == NULL){
        ROM_std_stdd_reduced_solver_finalize(solver);
        return ROM_STDD_ERR_ALLOCATION;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_reduced_solver_assemble(
    ROM_STDD_REDUCED_SOLVER* solver,
    const ROM_STDD_SYSTEM* system)
{
    int block_slots;

    if(solver == NULL || system == NULL ||
       system->diagonal_blocks == NULL || system->temporal_blocks == NULL ||
       system->reduced_rhs == NULL || solver->block_size != system->num_modes ||
       solver->num_windows != system->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    monolis_clear_mat_value_rhs_R(&solver->monolis);
    memset(
        solver->solution,
        0,
        (size_t)solver->num_total_meta_nodes
            * (size_t)solver->block_size
            * sizeof(double));

    block_slots = 1 + system->num_neighbor_ranks;

    for(int w = 0; w < system->num_windows; w++){
        for(int slot = 0; slot < block_slots; slot++){
            const int current_column = slot == 0
                ? w
                : stdd_halo_meta_node(solver, slot - 1, w);
            const size_t diagonal = stdd_block_offset(system, slot);

            for(int i = 0; i < system->num_modes; i++){
                for(int j = 0; j < system->num_modes; j++){
                    monolis_add_scalar_to_sparse_matrix_R(
                        &solver->monolis,
                        w,
                        current_column,
                        i,
                        j,
                        system->diagonal_blocks[
                            diagonal
                            + (size_t)i * (size_t)system->num_modes
                            + (size_t)j]);
                }
            }

            if(w > 0){
                const int previous_column = slot == 0
                    ? w - 1
                    : stdd_halo_meta_node(solver, slot - 1, w - 1);
                const size_t temporal = stdd_block_offset(system, slot);

                for(int i = 0; i < system->num_modes; i++){
                    for(int j = 0; j < system->num_modes; j++){
                        monolis_add_scalar_to_sparse_matrix_R(
                            &solver->monolis,
                            w,
                            previous_column,
                            i,
                            j,
                            -system->temporal_blocks[
                                temporal
                                + (size_t)i * (size_t)system->num_modes
                                + (size_t)j]);
                    }
                }
            }
        }

        for(int i = 0; i < system->num_modes; i++){
            solver->monolis.mat.R.B[
                (size_t)w * (size_t)system->num_modes
                + (size_t)i] =
                system->reduced_rhs[
                    (size_t)w * (size_t)system->num_modes_capacity
                    + (size_t)i];
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_reduced_solver_solve(
    ROM_STDD_REDUCED_SOLVER* solver,
    ROM_STDD_SYSTEM* system)
{
    if(solver == NULL || system == NULL ||
       solver->solution == NULL || system->reduced_coef == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    BBFE_sys_monowrap_solve(
        &solver->monolis,
        &solver->comm,
        solver->solution,
        MONOLIS_ITER_BICGSTAB,
        MONOLIS_PREC_DIAG,
        solver->max_iter,
        solver->tolerance);

    for(int w = 0; w < system->num_windows; w++){
        for(int mode = 0; mode < system->num_modes; mode++){
            system->reduced_coef[
                (size_t)w * (size_t)system->num_modes_capacity
                + (size_t)mode] =
                solver->solution[
                    (size_t)w * (size_t)solver->block_size
                    + (size_t)mode];
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_solve_window_marching(
    ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int max_iter,
    double tolerance)
{
    ROM_STDD_REDUCED_SOLVER solver;
    double* self_lu = NULL;
    double* scaled_diagonal_blocks = NULL;
    double* previous_coefficients = NULL;
    double* reduced_rhs = NULL;
    double* work = NULL;
    int* pivot = NULL;
    int graph_position = 0;
    int block_slots;
    int monolis_initialized = 0;
    int status = ROM_STDD_SUCCESS;
    size_t block_size2;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       spatial_com == NULL || system->diagonal_blocks == NULL ||
       system->temporal_blocks == NULL || system->reduced_rhs == NULL ||
       system->reduced_coef == NULL || system->num_modes <= 0 ||
       system->num_neighbor_ranks != spatial_com->recv_n_neib ||
       max_iter <= 0 || !isfinite(tolerance) || tolerance <= 0.0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(&solver, 0, sizeof(solver));

    solver.num_windows = 1;
    solver.block_size = system->num_modes;
    solver.num_internal_meta_nodes = 1;
    solver.num_neighbor_ranks = system->num_neighbor_ranks;
    solver.num_total_meta_nodes = 1 + system->num_neighbor_ranks;
    solver.max_iter = max_iter;
    solver.tolerance = tolerance;

    status = stdd_allocate_reduced_comm(&solver.comm, spatial_com, 1);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    solver.graph_index = (int*)stdd_calloc(
        (size_t)solver.num_total_meta_nodes + 1,
        sizeof(int));
    solver.num_graph_items = 1 + 2 * solver.num_neighbor_ranks;
    solver.graph_item = (int*)stdd_calloc(
        (size_t)solver.num_graph_items,
        sizeof(int));

    if(solver.graph_index == NULL || solver.graph_item == NULL){
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    solver.graph_index[0] = 0;
    solver.graph_item[graph_position++] = 0;
    for(int neighbor = 0; neighbor < solver.num_neighbor_ranks; neighbor++){
        solver.graph_item[graph_position++] = 1 + neighbor;
    }
    solver.graph_index[1] = graph_position;

    for(int node = 1; node < solver.num_total_meta_nodes; node++){
        solver.graph_item[graph_position++] = node;
        solver.graph_index[node + 1] = graph_position;
    }

    if(graph_position != solver.num_graph_items){
        status = ROM_STDD_ERR_DIMENSION;
        goto cleanup;
    }

    monolis_initialize(&solver.monolis);
    monolis_initialized = 1;
    monolis_get_nonzero_pattern_by_nodal_graph_R(
        &solver.monolis,
        solver.num_total_meta_nodes,
        solver.block_size,
        solver.graph_index,
        solver.graph_item);
    solver.monolis.mat.N = solver.num_internal_meta_nodes;
    solver.monolis.mat.NP = solver.num_total_meta_nodes;

    solver.solution = (double*)stdd_calloc(
        (size_t)solver.num_total_meta_nodes * (size_t)solver.block_size,
        sizeof(double));

    block_slots = 1 + system->num_neighbor_ranks;
    block_size2 =
        (size_t)system->num_modes * (size_t)system->num_modes;

    self_lu = (double*)stdd_calloc(block_size2, sizeof(double));
    scaled_diagonal_blocks = (double*)stdd_calloc(
        (size_t)block_slots * block_size2,
        sizeof(double));
    previous_coefficients = (double*)stdd_calloc(
        (size_t)solver.num_total_meta_nodes * (size_t)solver.block_size,
        sizeof(double));
    reduced_rhs = (double*)stdd_calloc(
        (size_t)system->num_modes,
        sizeof(double));
    work = (double*)stdd_calloc(
        (size_t)system->num_modes,
        sizeof(double));
    pivot = (int*)stdd_calloc(
        (size_t)system->num_modes,
        sizeof(int));

    if(solver.solution == NULL || self_lu == NULL ||
       scaled_diagonal_blocks == NULL || previous_coefficients == NULL ||
       reduced_rhs == NULL || work == NULL || pivot == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    memcpy(self_lu, system->diagonal_blocks, block_size2 * sizeof(double));

    status = stdd_dense_lu_factor(
        self_lu,
        system->num_modes,
        pivot);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    for(int slot = 0; slot < block_slots; slot++){
        status = stdd_left_scale_dense_block(
            self_lu,
            pivot,
            system->num_modes,
            system->diagonal_blocks + (size_t)slot * block_size2,
            scaled_diagonal_blocks + (size_t)slot * block_size2,
            work);
        if(status != ROM_STDD_SUCCESS){
            goto cleanup;
        }
    }

    monolis_clear_mat_value_rhs_R(&solver.monolis);

    for(int slot = 0; slot < block_slots; slot++){
        const int column = slot;
        const double* block =
            scaled_diagonal_blocks + (size_t)slot * block_size2;

        for(int row_mode = 0; row_mode < system->num_modes; row_mode++){
            for(int column_mode = 0;
                column_mode < system->num_modes;
                column_mode++)
            {
                monolis_add_scalar_to_sparse_matrix_R(
                    &solver.monolis,
                    0,
                    column,
                    row_mode,
                    column_mode,
                    block[
                        (size_t)row_mode * (size_t)system->num_modes
                        + (size_t)column_mode]);
            }
        }
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-DDROM reduced solver: window marching, "
            "distributed D solve with local dense block scaling\n");
        printf(
            "ST-DDROM reduced solve per window: global dof=%d, "
            "local block size=%d\n",
            monolis_mpi_get_global_comm_size() * system->num_modes,
            system->num_modes);
    }

    for(int window = 0; window < system->num_windows; window++){
        const double* stored_rhs = system->reduced_rhs
            + (size_t)window * (size_t)system->num_modes_capacity;

        memset(
            previous_coefficients,
            0,
            (size_t)solver.num_total_meta_nodes
                * (size_t)solver.block_size
                * sizeof(double));

        if(window > 0){
            const double* previous_internal = system->reduced_coef
                + (size_t)(window - 1)
                    * (size_t)system->num_modes_capacity;

            memcpy(
                previous_coefficients,
                previous_internal,
                (size_t)system->num_modes * sizeof(double));

            monolis_mpi_update_R(
                &solver.comm,
                solver.num_total_meta_nodes,
                solver.block_size,
                previous_coefficients);
        }

        for(int row_mode = 0; row_mode < system->num_modes; row_mode++){
            double value = stored_rhs[row_mode];

            if(window > 0){
                for(int slot = 0; slot < block_slots; slot++){
                    const double* temporal_block = system->temporal_blocks
                        + (size_t)slot * block_size2;
                    const double* previous_slot = previous_coefficients
                        + (size_t)slot * (size_t)system->num_modes;

                    for(int column_mode = 0;
                        column_mode < system->num_modes;
                        column_mode++)
                    {
                        value += temporal_block[
                            (size_t)row_mode * (size_t)system->num_modes
                            + (size_t)column_mode]
                            * previous_slot[column_mode];
                    }
                }
            }

            reduced_rhs[row_mode] = value;
        }

        status = stdd_dense_lu_solve(
            self_lu,
            pivot,
            system->num_modes,
            reduced_rhs);
        if(status != ROM_STDD_SUCCESS){
            goto cleanup;
        }

        memset(
            solver.monolis.mat.R.B,
            0,
            (size_t)solver.num_total_meta_nodes
                * (size_t)solver.block_size
                * sizeof(double));
        memcpy(
            solver.monolis.mat.R.B,
            reduced_rhs,
            (size_t)system->num_modes * sizeof(double));

        memset(
            solver.solution,
            0,
            (size_t)solver.num_total_meta_nodes
                * (size_t)solver.block_size
                * sizeof(double));

        if(window > 0){
            memcpy(
                solver.solution,
                previous_coefficients,
                (size_t)system->num_modes * sizeof(double));
        }

        BBFE_sys_monowrap_solve(
            &solver.monolis,
            &solver.comm,
            solver.solution,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG,
            solver.max_iter,
            solver.tolerance);

        memcpy(
            system->reduced_coef
                + (size_t)window * (size_t)system->num_modes_capacity,
            solver.solution,
            (size_t)system->num_modes * sizeof(double));
    }

cleanup:
    if(monolis_initialized){
        monolis_finalize(&solver.monolis);
    }
    stdd_free_reduced_comm(&solver.comm);
    free(solver.graph_index);
    free(solver.graph_item);
    free(solver.solution);
    free(self_lu);
    free(scaled_diagonal_blocks);
    free(previous_coefficients);
    free(reduced_rhs);
    free(work);
    free(pivot);

    return status;
}

