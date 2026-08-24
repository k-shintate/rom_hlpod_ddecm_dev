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

