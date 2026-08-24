static double** stdd_alloc_matrix_rows(int rows, int columns)
{
    double** matrix;
    double* data;

    if(rows <= 0 || columns <= 0){
        return NULL;
    }

    matrix = (double**)stdd_calloc((size_t)rows, sizeof(double*));
    data = (double*)stdd_calloc(
        (size_t)rows * (size_t)columns,
        sizeof(double));

    if(matrix == NULL || data == NULL){
        free(matrix);
        free(data);
        return NULL;
    }

    for(int i = 0; i < rows; i++){
        matrix[i] = data + (size_t)i * (size_t)columns;
    }

    return matrix;
}

static void stdd_free_matrix_rows(double*** matrix)
{
    if(matrix == NULL || *matrix == NULL){
        return;
    }
    free((*matrix)[0]);
    free(*matrix);
    *matrix = NULL;
}


static int stdd_local_lapack_svd(
    const double* row_major_matrix,
    int rows,
    int columns,
    double** left_vectors,
    double* singular_values)
{
    const int minimum_dimension = rows < columns ? rows : columns;
    double* a = NULL;
    double* u = NULL;
    double* work = NULL;
    double work_query = 0.0;
    double vt_dummy = 0.0;
    char jobu = 'S';
    char jobvt = 'N';
    int lda = rows;
    int ldu = rows;
    int ldvt = 1;
    int lwork = -1;
    int info = 0;

    if(row_major_matrix == NULL || left_vectors == NULL ||
       singular_values == NULL || rows <= 0 || columns <= 0 ||
       minimum_dimension <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if((size_t)rows > SIZE_MAX / (size_t)columns ||
       (size_t)rows > SIZE_MAX / (size_t)minimum_dimension)
    {
        return ROM_STDD_ERR_DIMENSION;
    }

    a = (double*)stdd_calloc(
        (size_t)rows * (size_t)columns,
        sizeof(double));
    u = (double*)stdd_calloc(
        (size_t)rows * (size_t)minimum_dimension,
        sizeof(double));

    if(a == NULL || u == NULL){
        free(a);
        free(u);
        return ROM_STDD_ERR_ALLOCATION;
    }

    /* Convert the generic ST-DDROM row-major snapshot matrix to LAPACK's
     * column-major A(rows,columns). */
    for(int row = 0; row < rows; row++){
        for(int column = 0; column < columns; column++){
            a[(size_t)column * (size_t)rows + (size_t)row] =
                row_major_matrix[
                    (size_t)row * (size_t)columns
                    + (size_t)column];
        }
    }

    /* Workspace query. */
    dgesvd_(
        &jobu,
        &jobvt,
        &rows,
        &columns,
        a,
        &lda,
        singular_values,
        u,
        &ldu,
        &vt_dummy,
        &ldvt,
        &work_query,
        &lwork,
        &info);

    if(info != 0 || !isfinite(work_query) || work_query < 1.0 ||
       work_query > (double)INT_MAX)
    {
        fprintf(
            stderr,
            "ERROR: local LAPACK dgesvd workspace query failed: "
            "rows=%d columns=%d info=%d work=%e\\n",
            rows,
            columns,
            info,
            work_query);
        free(a);
        free(u);
        return ROM_STDD_ERR_OPERATOR;
    }

    lwork = (int)ceil(work_query);
    if(lwork < 1){
        lwork = 1;
    }

    work = (double*)stdd_calloc((size_t)lwork, sizeof(double));
    if(work == NULL){
        free(a);
        free(u);
        return ROM_STDD_ERR_ALLOCATION;
    }

    /* dgesvd destroys A, so restore it after the workspace-query call. */
    for(int row = 0; row < rows; row++){
        for(int column = 0; column < columns; column++){
            a[(size_t)column * (size_t)rows + (size_t)row] =
                row_major_matrix[
                    (size_t)row * (size_t)columns
                    + (size_t)column];
        }
    }

    info = 0;
    dgesvd_(
        &jobu,
        &jobvt,
        &rows,
        &columns,
        a,
        &lda,
        singular_values,
        u,
        &ldu,
        &vt_dummy,
        &ldvt,
        work,
        &lwork,
        &info);

    if(info != 0){
        fprintf(
            stderr,
            "ERROR: local LAPACK dgesvd failed: "
            "rows=%d columns=%d info=%d\\n",
            rows,
            columns,
            info);
        free(work);
        free(a);
        free(u);
        return ROM_STDD_ERR_OPERATOR;
    }

    /* U is column-major.  Preserve the generic basis storage convention
     * left_vectors[row][mode]. */
    for(int mode = 0; mode < minimum_dimension; mode++){
        for(int row = 0; row < rows; row++){
            left_vectors[row][mode] =
                u[(size_t)mode * (size_t)rows + (size_t)row];
        }
    }

    free(work);
    free(a);
    free(u);
    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_compute_shared_local_pod(
    ROM_STDD_SYSTEM* system,
    const double* snapshot_matrix,
    int num_snapshot_columns,
    double energy_tolerance,
    int communicator)
{
    double** left_vectors = NULL;
    double* singular_values = NULL;
    int local_modes;
    int global_modes;
    int maximum_possible_modes;
    int global_minimum_possible_modes;
    int negative_minimum_possible_modes;
    int status = ROM_STDD_SUCCESS;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       snapshot_matrix == NULL || num_snapshot_columns <= 0 ||
       !isfinite(energy_tolerance) || energy_tolerance <= 0.0 ||
       energy_tolerance > 1.0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    maximum_possible_modes = system->internal_st_dof < num_snapshot_columns
        ? system->internal_st_dof
        : num_snapshot_columns;
    if(system->num_modes_capacity < maximum_possible_modes){
        maximum_possible_modes = system->num_modes_capacity;
    }

    left_vectors = stdd_alloc_matrix_rows(
        system->internal_st_dof,
        num_snapshot_columns);
    singular_values = (double*)stdd_calloc(
        (size_t)num_snapshot_columns,
        sizeof(double));

    if(left_vectors == NULL || singular_values == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    /*
     * POD is independent on each spatial rank.  Use local LAPACK instead of
     * creating a one-rank ScaLAPACK/BLACS grid.  The resulting left singular
     * subspace is the same POD space; only the linear-algebra backend changes.
     */
    status = stdd_local_lapack_svd(
        snapshot_matrix,
        system->internal_st_dof,
        num_snapshot_columns,
        left_vectors,
        singular_values);
    if(status != ROM_STDD_SUCCESS){
        goto cleanup;
    }

    if(energy_tolerance >= 1.0){
        local_modes = maximum_possible_modes;
    }
    else{
        local_modes = ROM_BB_estimate_num_pod_modes(
            singular_values,
            num_snapshot_columns,
            energy_tolerance);
        if(local_modes < 1){
            local_modes = 1;
        }
        if(local_modes > maximum_possible_modes){
            local_modes = maximum_possible_modes;
        }
    }

    /*
     * MONOLIS fixed-block BCSR requires one common block size.  Promote the
     * local estimate to the maximum rank selected by any spatial rank.
     */
    global_modes = local_modes;
    monolis_allreduce_I(
        1,
        &global_modes,
        MONOLIS_MPI_MAX,
        communicator);

    /*
     * Use MAX on the negated value so this remains compatible with
     * MONOLIS versions that expose SUM/MAX but not MONOLIS_MPI_MIN.
     */
    negative_minimum_possible_modes = -maximum_possible_modes;
    monolis_allreduce_I(
        1,
        &negative_minimum_possible_modes,
        MONOLIS_MPI_MAX,
        communicator);
    global_minimum_possible_modes =
        -negative_minimum_possible_modes;

    if(global_modes > global_minimum_possible_modes){
        global_modes = global_minimum_possible_modes;
    }

    if(global_modes <= 0 || global_modes > maximum_possible_modes){
        status = ROM_STDD_ERR_DIMENSION;
        goto cleanup;
    }

    free(system->basis);
    free(system->singular_values);
    system->basis = NULL;
    system->singular_values = NULL;

    system->basis = (double*)stdd_calloc(
        (size_t)system->internal_st_dof
            * (size_t)system->num_modes_capacity,
        sizeof(double));
    system->singular_values = (double*)stdd_calloc(
        (size_t)num_snapshot_columns,
        sizeof(double));

    if(system->basis == NULL || system->singular_values == NULL){
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    system->num_modes = global_modes;
    system->num_singular_values = num_snapshot_columns;

    memcpy(
        system->singular_values,
        singular_values,
        (size_t)num_snapshot_columns * sizeof(double));

    for(int row = 0; row < system->internal_st_dof; row++){
        for(int mode = 0; mode < global_modes; mode++){
            system->basis[
                stdd_basis_offset(system, row, mode)] =
                left_vectors[row][mode];
        }
    }

cleanup:
    stdd_free_matrix_rows(&left_vectors);
    free(singular_values);

    return status;
}

int ROM_std_stdd_write_basis_rank(
    const ROM_STDD_SYSTEM* system,
    int num_snapshot_columns,
    const char* directory)
{
    ROM_STDD_BASIS_HEADER header;
    char filename[4096];
    FILE* fp;
    size_t basis_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       system->basis == NULL || system->singular_values == NULL ||
       system->num_modes <= 0 || num_snapshot_columns <= 0 ||
       directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_make_directory(directory) != ROM_STDD_SUCCESS ||
       stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_basis",
        0,
        0) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    memset(&header, 0, sizeof(header));
    header.magic = ROM_STDD_BASIS_MAGIC;
    header.version = ROM_STDD_FILE_VERSION;
    header.mpi_size = monolis_mpi_get_global_comm_size();
    header.mpi_rank = monolis_mpi_get_global_my_rank();
    header.num_windows = system->num_windows;
    header.num_internal_space_nodes = system->num_internal_space_nodes;
    header.st_dof_per_space_node = system->st_dof_per_space_node;
    header.internal_st_dof = system->internal_st_dof;
    header.num_modes = system->num_modes;
    header.num_snapshot_columns = num_snapshot_columns;
    header.num_singular_values = system->num_singular_values;

    fp = fopen(filename, "wb");
    if(fp == NULL){
        fprintf(stderr, "ERROR: cannot write %s: %s\n", filename, strerror(errno));
        return ROM_STDD_ERR_IO;
    }

    basis_count = (size_t)system->internal_st_dof
        * (size_t)system->num_modes;

    if(stdd_write_exact(&header, sizeof(header), 1, fp) != ROM_STDD_SUCCESS ||
       stdd_write_exact(
        system->singular_values,
        sizeof(double),
        (size_t)system->num_singular_values,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    (void)basis_count;
    for(int row = 0; row < system->internal_st_dof; row++){
        if(stdd_write_exact(
            system->basis
                + (size_t)row * (size_t)system->num_modes_capacity,
            sizeof(double),
            (size_t)system->num_modes,
            fp) != ROM_STDD_SUCCESS)
        {
            fclose(fp);
            return ROM_STDD_ERR_IO;
        }
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_read_basis_rank(
    ROM_STDD_SYSTEM* system,
    const char* directory)
{
    ROM_STDD_BASIS_HEADER header;
    char filename[4096];
    FILE* fp;
    size_t basis_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_basis",
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

    if(header.magic != ROM_STDD_BASIS_MAGIC ||
       header.version != ROM_STDD_FILE_VERSION ||
       header.mpi_size != monolis_mpi_get_global_comm_size() ||
       header.mpi_rank != monolis_mpi_get_global_my_rank() ||
       header.num_windows != system->num_windows ||
       header.num_internal_space_nodes != system->num_internal_space_nodes ||
       header.st_dof_per_space_node != system->st_dof_per_space_node ||
       header.internal_st_dof != system->internal_st_dof ||
       header.num_modes <= 0 ||
       header.num_modes > system->num_modes_capacity ||
       header.num_singular_values <= 0)
    {
        fclose(fp);
        fprintf(stderr, "ERROR: incompatible ST-DDROM basis: %s\n", filename);
        return ROM_STDD_ERR_INCOMPATIBLE;
    }

    free(system->basis);
    free(system->singular_values);
    system->basis = NULL;
    system->singular_values = NULL;

    basis_count = (size_t)system->internal_st_dof
        * (size_t)system->num_modes_capacity;

    system->basis = (double*)stdd_calloc(basis_count, sizeof(double));
    system->singular_values = (double*)stdd_calloc(
        (size_t)header.num_singular_values,
        sizeof(double));

    if(system->basis == NULL || system->singular_values == NULL){
        fclose(fp);
        return ROM_STDD_ERR_ALLOCATION;
    }

    system->num_modes = header.num_modes;
    system->num_singular_values = header.num_singular_values;

    if(stdd_read_exact(
        system->singular_values,
        sizeof(double),
        (size_t)system->num_singular_values,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    for(int row = 0; row < system->internal_st_dof; row++){
        if(stdd_read_exact(
            system->basis
                + (size_t)row * (size_t)system->num_modes_capacity,
            sizeof(double),
            (size_t)system->num_modes,
            fp) != ROM_STDD_SUCCESS)
        {
            fclose(fp);
            return ROM_STDD_ERR_IO;
        }
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

