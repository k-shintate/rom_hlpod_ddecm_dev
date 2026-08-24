int ROM_std_stdd_project_operators_rank_sweep(
    ROM_STDD_SYSTEM* system,
    ROM_STDD_APPLY_OPERATOR apply_A,
    ROM_STDD_APPLY_COUPLING apply_B,
    void* operator_context)
{
    double* x;
    double* y;
    int rank;
    int comm_size;
    int block_slots;
    size_t block_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->basis == NULL || system->num_modes <= 0 ||
       apply_A == NULL || apply_B == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    rank = monolis_mpi_get_global_my_rank();
    comm_size = monolis_mpi_get_global_comm_size();

    block_slots = 1 + system->num_neighbor_ranks;
    block_count = (size_t)block_slots
        * (size_t)system->num_modes
        * (size_t)system->num_modes;

    free(system->diagonal_blocks);
    free(system->temporal_blocks);
    system->diagonal_blocks = (double*)stdd_calloc(block_count, sizeof(double));
    system->temporal_blocks = (double*)stdd_calloc(block_count, sizeof(double));

    x = (double*)stdd_calloc((size_t)system->local_st_dof, sizeof(double));
    y = (double*)stdd_calloc((size_t)system->local_st_dof, sizeof(double));

    if(system->diagonal_blocks == NULL || system->temporal_blocks == NULL ||
       x == NULL || y == NULL)
    {
        free(x);
        free(y);
        return ROM_STDD_ERR_ALLOCATION;
    }

    for(int source_rank = 0; source_rank < comm_size; source_rank++){
        const int slot = stdd_find_source_slot(system, source_rank, rank);

        for(int source_mode = 0;
            source_mode < system->num_modes;
            source_mode++)
        {
            memset(x, 0, (size_t)system->local_st_dof * sizeof(double));
            memset(y, 0, (size_t)system->local_st_dof * sizeof(double));

            if(rank == source_rank){
                for(int row = 0; row < system->internal_st_dof; row++){
                    x[row] = system->basis[
                        stdd_basis_offset(system, row, source_mode)];
                }
            }

            if(apply_A(x, y, operator_context) != ROM_STDD_SUCCESS){
                free(x);
                free(y);
                return ROM_STDD_ERR_OPERATOR;
            }

            if(slot >= 0){
                const size_t block = stdd_block_offset(system, slot);
                for(int row_mode = 0;
                    row_mode < system->num_modes;
                    row_mode++)
                {
                    double value = 0.0;
                    for(int row = 0; row < system->internal_st_dof; row++){
                        value += system->basis[
                            stdd_basis_offset(system, row, row_mode)]
                            * y[row];
                    }
                    system->diagonal_blocks[
                        block
                        + (size_t)row_mode * (size_t)system->num_modes
                        + (size_t)source_mode] = value;
                }
            }

            memset(y, 0, (size_t)system->local_st_dof * sizeof(double));

            if(apply_B(x, y, operator_context) != ROM_STDD_SUCCESS){
                free(x);
                free(y);
                return ROM_STDD_ERR_OPERATOR;
            }

            if(slot >= 0){
                const size_t block = stdd_block_offset(system, slot);
                for(int row_mode = 0;
                    row_mode < system->num_modes;
                    row_mode++)
                {
                    double value = 0.0;
                    for(int row = 0; row < system->internal_st_dof; row++){
                        value += system->basis[
                            stdd_basis_offset(system, row, row_mode)]
                            * y[row];
                    }
                    system->temporal_blocks[
                        block
                        + (size_t)row_mode * (size_t)system->num_modes
                        + (size_t)source_mode] = value;
                }
            }
        }
    }

    free(x);
    free(y);
    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_write_operators_rank(
    const ROM_STDD_SYSTEM* system,
    const char* directory)
{
    ROM_STDD_OPERATOR_HEADER header;
    char filename[4096];
    FILE* fp;
    size_t block_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->num_modes <= 0 || system->diagonal_blocks == NULL ||
       system->temporal_blocks == NULL || directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_make_directory(directory) != ROM_STDD_SUCCESS ||
       stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_operators",
        0,
        0) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    memset(&header, 0, sizeof(header));
    header.magic = ROM_STDD_OPERATOR_MAGIC;
    header.version = ROM_STDD_FILE_VERSION;
    header.mpi_size = monolis_mpi_get_global_comm_size();
    header.mpi_rank = monolis_mpi_get_global_my_rank();
    header.num_windows = system->num_windows;
    header.num_modes = system->num_modes;
    header.num_neighbor_ranks = system->num_neighbor_ranks;
    header.internal_st_dof = system->internal_st_dof;

    fp = fopen(filename, "wb");
    if(fp == NULL){
        fprintf(stderr, "ERROR: cannot write %s: %s\n", filename, strerror(errno));
        return ROM_STDD_ERR_IO;
    }

    block_count = (size_t)(1 + system->num_neighbor_ranks)
        * (size_t)system->num_modes
        * (size_t)system->num_modes;

    if(stdd_write_exact(&header, sizeof(header), 1, fp) != ROM_STDD_SUCCESS ||
       stdd_write_exact(
        system->neighbor_ranks,
        sizeof(int),
        (size_t)system->num_neighbor_ranks,
        fp) != ROM_STDD_SUCCESS ||
       stdd_write_exact(
        system->diagonal_blocks,
        sizeof(double),
        block_count,
        fp) != ROM_STDD_SUCCESS ||
       stdd_write_exact(
        system->temporal_blocks,
        sizeof(double),
        block_count,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_read_operators_rank(
    ROM_STDD_SYSTEM* system,
    const char* directory)
{
    ROM_STDD_OPERATOR_HEADER header;
    char filename[4096];
    FILE* fp;
    int* file_neighbors = NULL;
    size_t block_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->num_modes <= 0 || directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_operators",
        0,
        0) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        fprintf(stderr, "ERROR: cannot read %s: %s\n", filename, strerror(errno));
        return ROM_STDD_ERR_IO;
    }

    if(stdd_read_exact(&header, sizeof(header), 1, fp) != ROM_STDD_SUCCESS){
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    if(header.magic != ROM_STDD_OPERATOR_MAGIC ||
       header.version != ROM_STDD_FILE_VERSION ||
       header.mpi_size != monolis_mpi_get_global_comm_size() ||
       header.mpi_rank != monolis_mpi_get_global_my_rank() ||
       header.num_windows != system->num_windows ||
       header.num_modes != system->num_modes ||
       header.num_neighbor_ranks != system->num_neighbor_ranks ||
       header.internal_st_dof != system->internal_st_dof)
    {
        fclose(fp);
        fprintf(stderr, "ERROR: incompatible ST-DDROM operators: %s\n", filename);
        return ROM_STDD_ERR_INCOMPATIBLE;
    }

    if(system->num_neighbor_ranks > 0){
        file_neighbors = (int*)stdd_calloc(
            (size_t)system->num_neighbor_ranks,
            sizeof(int));
        if(file_neighbors == NULL){
            fclose(fp);
            return ROM_STDD_ERR_ALLOCATION;
        }

        if(stdd_read_exact(
            file_neighbors,
            sizeof(int),
            (size_t)system->num_neighbor_ranks,
            fp) != ROM_STDD_SUCCESS)
        {
            free(file_neighbors);
            fclose(fp);
            return ROM_STDD_ERR_IO;
        }

        if(memcmp(
            file_neighbors,
            system->neighbor_ranks,
            (size_t)system->num_neighbor_ranks * sizeof(int)) != 0)
        {
            free(file_neighbors);
            fclose(fp);
            fprintf(stderr, "ERROR: reduced-neighbour order changed: %s\n", filename);
            return ROM_STDD_ERR_INCOMPATIBLE;
        }
        free(file_neighbors);
    }

    block_count = (size_t)(1 + system->num_neighbor_ranks)
        * (size_t)system->num_modes
        * (size_t)system->num_modes;

    free(system->diagonal_blocks);
    free(system->temporal_blocks);
    system->diagonal_blocks = (double*)stdd_calloc(block_count, sizeof(double));
    system->temporal_blocks = (double*)stdd_calloc(block_count, sizeof(double));

    if(system->diagonal_blocks == NULL || system->temporal_blocks == NULL){
        fclose(fp);
        return ROM_STDD_ERR_ALLOCATION;
    }

    if(stdd_read_exact(
        system->diagonal_blocks,
        sizeof(double),
        block_count,
        fp) != ROM_STDD_SUCCESS ||
       stdd_read_exact(
        system->temporal_blocks,
        sizeof(double),
        block_count,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_project_local_rhs(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* effective_full_rhs,
    double* reduced_rhs_out)
{
    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->basis == NULL || system->num_modes <= 0 ||
       effective_full_rhs == NULL || reduced_rhs_out == NULL ||
       window_id < 0 || window_id >= system->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int mode = 0; mode < system->num_modes; mode++){
        double value = 0.0;
        for(int row = 0; row < system->internal_st_dof; row++){
            value += system->basis[
                stdd_basis_offset(system, row, mode)]
                * effective_full_rhs[row];
        }
        reduced_rhs_out[mode] = value;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_reconstruct_local_homogeneous(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* local_homogeneous_state)
{
    const double* coefficient;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->basis == NULL || system->num_modes <= 0 ||
       system->reduced_coef == NULL || local_homogeneous_state == NULL ||
       window_id < 0 || window_id >= system->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(
        local_homogeneous_state,
        0,
        (size_t)system->local_st_dof * sizeof(double));

    coefficient = system->reduced_coef
        + (size_t)window_id * (size_t)system->num_modes_capacity;

    for(int row = 0; row < system->internal_st_dof; row++){
        double value = 0.0;
        for(int mode = 0; mode < system->num_modes; mode++){
            value += system->basis[
                stdd_basis_offset(system, row, mode)]
                * coefficient[mode];
        }
        local_homogeneous_state[row] = value;
    }

    return ROM_STDD_SUCCESS;
}

