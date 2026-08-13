#include "rom_std_stddrom.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


/*
 * Rank-local POD uses a rank-local LAPACK SVD.
 *
 * The validated ST-DDROM POD is independent on each spatial MPI rank.  The
 * previous implementation called the MONOLIS ScaLAPACK wrapper on MPI_COMM_SELF,
 * which unnecessarily pulled ScaLAPACK/BLACS into the final executable.
 * Keeping the SVD local removes that link-time dependency without changing the
 * POD algebra or any online reduced operator/solver routine.
 */
#ifdef __cplusplus
extern "C" {
#endif
void dgesvd_(
    char* jobu,
    char* jobvt,
    int* m,
    int* n,
    double* a,
    int* lda,
    double* s,
    double* u,
    int* ldu,
    double* vt,
    int* ldvt,
    double* work,
    int* lwork,
    int* info);
#ifdef __cplusplus
}
#endif

#define ROM_STDD_SNAPSHOT_MAGIC UINT64_C(0x53544444534E5031) /* STDDSNP1 */
#define ROM_STDD_BASIS_MAGIC    UINT64_C(0x5354444442415331) /* STDDBAS1 */
#define ROM_STDD_OPERATOR_MAGIC UINT64_C(0x535444444F505231) /* STDDOPR1 */
#define ROM_STDD_FILE_VERSION   UINT32_C(3)

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t reserved;
    int32_t snapshot_id;
    int32_t mpi_size;
    int32_t mpi_rank;
    int32_t num_windows;
    int32_t num_internal_space_nodes;
    int32_t st_dof_per_space_node;
    int32_t internal_st_dof;
} ROM_STDD_SNAPSHOT_HEADER;

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t reserved;
    int32_t mpi_size;
    int32_t mpi_rank;
    int32_t num_windows;
    int32_t num_internal_space_nodes;
    int32_t st_dof_per_space_node;
    int32_t internal_st_dof;
    int32_t num_modes;
    int32_t num_snapshot_columns;
    int32_t num_singular_values;
} ROM_STDD_BASIS_HEADER;

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t reserved;
    int32_t mpi_size;
    int32_t mpi_rank;
    int32_t num_windows;
    int32_t num_modes;
    int32_t num_neighbor_ranks;
    int32_t internal_st_dof;
} ROM_STDD_OPERATOR_HEADER;

static void* stdd_calloc(size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size){
        return NULL;
    }
    return calloc(count, size);
}

static int stdd_checked_product_int(int a, int b, int* result)
{
    long long value;

    if(result == NULL || a < 0 || b < 0){
        return ROM_STDD_ERR_ARGUMENT;
    }

    value = (long long)a * (long long)b;
    if(value > INT_MAX){
        return ROM_STDD_ERR_DIMENSION;
    }

    *result = (int)value;
    return ROM_STDD_SUCCESS;
}

static int stdd_make_directory(const char* directory)
{
    if(directory == NULL || directory[0] == '\0'){
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(mkdir(directory, 0775) == 0 || errno == EEXIST){
        return ROM_STDD_SUCCESS;
    }

    fprintf(stderr, "ERROR: mkdir(%s): %s\n", directory, strerror(errno));
    return ROM_STDD_ERR_IO;
}

static int stdd_build_rank_filename(
    char* filename,
    size_t filename_size,
    const char* directory,
    const char* stem,
    int identifier,
    int include_identifier)
{
    const int rank = monolis_mpi_get_global_my_rank();
    int written;

    if(filename == NULL || filename_size == 0 || directory == NULL ||
       stem == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(include_identifier){
        written = snprintf(
            filename,
            filename_size,
            "%s/%s_%04d_rank_%06d.bin",
            directory,
            stem,
            identifier,
            rank);
    }
    else{
        written = snprintf(
            filename,
            filename_size,
            "%s/%s_rank_%06d.bin",
            directory,
            stem,
            rank);
    }

    if(written < 0 || (size_t)written >= filename_size){
        return ROM_STDD_ERR_IO;
    }

    return ROM_STDD_SUCCESS;
}

static int stdd_write_exact(
    const void* data,
    size_t element_size,
    size_t count,
    FILE* fp)
{
    if(fp == NULL || (count > 0 && data == NULL)){
        return ROM_STDD_ERR_ARGUMENT;
    }
    return fwrite(data, element_size, count, fp) == count
        ? ROM_STDD_SUCCESS
        : ROM_STDD_ERR_IO;
}

static int stdd_read_exact(
    void* data,
    size_t element_size,
    size_t count,
    FILE* fp)
{
    if(fp == NULL || (count > 0 && data == NULL)){
        return ROM_STDD_ERR_ARGUMENT;
    }
    return fread(data, element_size, count, fp) == count
        ? ROM_STDD_SUCCESS
        : ROM_STDD_ERR_IO;
}

static int stdd_find_source_slot(
    const ROM_STDD_SYSTEM* system,
    int source_rank,
    int my_rank)
{
    if(source_rank == my_rank){
        return 0;
    }

    for(int n = 0; n < system->num_neighbor_ranks; n++){
        if(system->neighbor_ranks[n] == source_rank){
            return n + 1;
        }
    }

    return -1;
}

static size_t stdd_block_offset(
    const ROM_STDD_SYSTEM* system,
    int slot)
{
    return (size_t)slot
        * (size_t)system->num_modes
        * (size_t)system->num_modes;
}

static size_t stdd_basis_offset(
    const ROM_STDD_SYSTEM* system,
    int row,
    int mode)
{
    return (size_t)row * (size_t)system->num_modes_capacity
        + (size_t)mode;
}

static int stdd_validate_system_dimensions(const ROM_STDD_SYSTEM* system)
{
    if(system == NULL || system->num_windows <= 0 ||
       system->num_space_nodes <= 0 ||
       system->num_internal_space_nodes <= 0 ||
       system->num_internal_space_nodes > system->num_space_nodes ||
       system->st_dof_per_space_node <= 0 ||
       system->local_st_dof <= 0 || system->internal_st_dof <= 0 ||
       system->num_modes_capacity <= 0)
    {
        return ROM_STDD_ERR_DIMENSION;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_system_initialize(
    ROM_STDD_SYSTEM* system,
    const MONOLIS_COM* spatial_com,
    int num_windows,
    int num_space_nodes,
    int num_internal_space_nodes,
    int st_dof_per_space_node,
    int num_modes_capacity)
{
    int local_st_dof;
    int internal_st_dof;

    if(system == NULL || spatial_com == NULL || num_windows <= 0 ||
       num_space_nodes <= 0 || num_internal_space_nodes <= 0 ||
       num_internal_space_nodes > num_space_nodes ||
       st_dof_per_space_node <= 0 || num_modes_capacity <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(system, 0, sizeof(*system));

    if(stdd_checked_product_int(
        num_space_nodes,
        st_dof_per_space_node,
        &local_st_dof) != ROM_STDD_SUCCESS ||
       stdd_checked_product_int(
        num_internal_space_nodes,
        st_dof_per_space_node,
        &internal_st_dof) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_DIMENSION;
    }

    system->num_windows = num_windows;
    system->num_modes_capacity = num_modes_capacity;
    system->num_space_nodes = num_space_nodes;
    system->num_internal_space_nodes = num_internal_space_nodes;
    system->st_dof_per_space_node = st_dof_per_space_node;
    system->local_st_dof = local_st_dof;
    system->internal_st_dof = internal_st_dof;
    system->num_neighbor_ranks = spatial_com->recv_n_neib;

    if(system->num_neighbor_ranks > 0){
        system->neighbor_ranks = (int*)stdd_calloc(
            (size_t)system->num_neighbor_ranks,
            sizeof(int));
        if(system->neighbor_ranks == NULL){
            ROM_std_stdd_system_finalize(system);
            return ROM_STDD_ERR_ALLOCATION;
        }
        memcpy(
            system->neighbor_ranks,
            spatial_com->recv_neib_pe,
            (size_t)system->num_neighbor_ranks * sizeof(int));
    }

    system->reduced_rhs = (double*)stdd_calloc(
        (size_t)num_windows * (size_t)num_modes_capacity,
        sizeof(double));
    system->reduced_coef = (double*)stdd_calloc(
        (size_t)num_windows * (size_t)num_modes_capacity,
        sizeof(double));

    if(system->reduced_rhs == NULL || system->reduced_coef == NULL){
        ROM_std_stdd_system_finalize(system);
        return ROM_STDD_ERR_ALLOCATION;
    }

    return ROM_STDD_SUCCESS;
}

void ROM_std_stdd_system_finalize(ROM_STDD_SYSTEM* system)
{
    if(system == NULL){
        return;
    }

    free(system->neighbor_ranks);
    free(system->basis);
    free(system->singular_values);
    free(system->diagonal_blocks);
    free(system->temporal_blocks);
    free(system->reduced_rhs);
    free(system->reduced_coef);

    memset(system, 0, sizeof(*system));
}

int ROM_std_stdd_write_snapshot_rank(
    const ROM_STDD_SYSTEM* system,
    const double* window_snapshots,
    int snapshot_id,
    const char* directory)
{
    ROM_STDD_SNAPSHOT_HEADER header;
    char filename[4096];
    FILE* fp;
    int mpi_size;
    int rank;
    size_t value_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       window_snapshots == NULL || snapshot_id < 0 || directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_make_directory(directory) != ROM_STDD_SUCCESS){
        return ROM_STDD_ERR_IO;
    }

    if(stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_snapshot",
        snapshot_id,
        1) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    mpi_size = monolis_mpi_get_global_comm_size();
    rank = monolis_mpi_get_global_my_rank();

    memset(&header, 0, sizeof(header));
    header.magic = ROM_STDD_SNAPSHOT_MAGIC;
    header.version = ROM_STDD_FILE_VERSION;
    header.snapshot_id = snapshot_id;
    header.mpi_size = mpi_size;
    header.mpi_rank = rank;
    header.num_windows = system->num_windows;
    header.num_internal_space_nodes = system->num_internal_space_nodes;
    header.st_dof_per_space_node = system->st_dof_per_space_node;
    header.internal_st_dof = system->internal_st_dof;

    fp = fopen(filename, "wb");
    if(fp == NULL){
        fprintf(stderr, "ERROR: cannot write %s: %s\n", filename, strerror(errno));
        return ROM_STDD_ERR_IO;
    }

    value_count = (size_t)system->num_windows
        * (size_t)system->internal_st_dof;

    if(stdd_write_exact(&header, sizeof(header), 1, fp) != ROM_STDD_SUCCESS ||
       stdd_write_exact(
        window_snapshots,
        sizeof(double),
        value_count,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    if(fclose(fp) != 0){
        return ROM_STDD_ERR_IO;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_read_snapshot_rank(
    const ROM_STDD_SYSTEM* system,
    double* window_snapshots,
    int snapshot_id,
    const char* directory)
{
    ROM_STDD_SNAPSHOT_HEADER header;
    char filename[4096];
    FILE* fp;
    size_t value_count;

    if(stdd_validate_system_dimensions(system) != ROM_STDD_SUCCESS ||
       window_snapshots == NULL || snapshot_id < 0 || directory == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(stdd_build_rank_filename(
        filename,
        sizeof(filename),
        directory,
        "stddrom_snapshot",
        snapshot_id,
        1) != ROM_STDD_SUCCESS)
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

    if(header.magic != ROM_STDD_SNAPSHOT_MAGIC ||
       header.version != ROM_STDD_FILE_VERSION ||
       header.snapshot_id != snapshot_id ||
       header.mpi_size != monolis_mpi_get_global_comm_size() ||
       header.mpi_rank != monolis_mpi_get_global_my_rank() ||
       header.num_windows != system->num_windows ||
       header.num_internal_space_nodes != system->num_internal_space_nodes ||
       header.st_dof_per_space_node != system->st_dof_per_space_node ||
       header.internal_st_dof != system->internal_st_dof)
    {
        fclose(fp);
        fprintf(stderr, "ERROR: incompatible ST-DDROM snapshot: %s\n", filename);
        return ROM_STDD_ERR_INCOMPATIBLE;
    }

    value_count = (size_t)system->num_windows
        * (size_t)system->internal_st_dof;

    if(stdd_read_exact(
        window_snapshots,
        sizeof(double),
        value_count,
        fp) != ROM_STDD_SUCCESS)
    {
        fclose(fp);
        return ROM_STDD_ERR_IO;
    }

    fclose(fp);
    return ROM_STDD_SUCCESS;
}

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
