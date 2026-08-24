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




