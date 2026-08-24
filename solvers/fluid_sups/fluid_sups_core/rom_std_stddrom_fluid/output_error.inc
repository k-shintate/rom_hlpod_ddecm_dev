int ROM_std_stdd_fluid_write_current_window_vtk(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    int* file_number)
{
    double* trace = NULL;
    int status = ROM_STDD_SUCCESS;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       file_number == NULL ||
       window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    trace = (double*)calloc(
        (size_t)context->fom->fe.total_num_nodes
            * FLUID_SUPS_PHYSICAL_DOF,
        sizeof(double));
    if(trace == NULL){
        return ROM_STDD_ERR_ALLOCATION;
    }

    for(int slab = 0;
        slab < context->layout->num_slabs;
        slab++)
    {
        const int global_step =
            window_id * context->layout->num_slabs
            + slab + 1;
        const double time =
            (double)global_step
            * context->fom->vals.dt;

        if(global_step
               % context->fom->vals.output_interval
           != 0)
        {
            continue;
        }

        if(ST_fom_extract_slab_right_trace(
            context->layout,
            slab,
            context->current_window,
            trace) != ST_FOM_SUCCESS ||
           fluid_sups_st_import_state(
            context->fom,
            trace,
            0) != 0 ||
           fluid_sups_st_synchronize_state(
            context->fom) != 0 ||
           fluid_sups_st_write_output_named(
            context->fom,
            "st_ddrom_result",
            *file_number,
            time) != 0)
        {
            status = ROM_STDD_ERR_IO;
            goto cleanup;
        }

        (*file_number)++;
    }

    /*
     * Restore the right trace as the persistent physical state for the next
     * window and for post-processing after this output loop.
     */
    if(fluid_sups_st_import_state(
        context->fom,
        context->previous_trace,
        1) != 0 ||
       fluid_sups_st_synchronize_state(
        context->fom) != 0)
    {
        status = ROM_STDD_ERR_OPERATOR;
    }

cleanup:
    free(trace);
    return status;
}

int ROM_std_stdd_fluid_write_current_window_binary(
    ROM_STDD_FLUID_CONTEXT* context,
    int window_id,
    const char* directory)
{
    char filename[4096];
    const int rank =
        monolis_mpi_get_global_my_rank();
    const double time_start =
        (double)(window_id
            * context->layout->num_slabs)
        * context->fom->vals.dt;
    const double time_end =
        time_start
        + (double)context->layout->num_slabs
            * context->fom->vals.dt;
    int written;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       directory == NULL || directory[0] == '\0' ||
       window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ST_fom_make_directory_recursive(directory)
       != ST_FOM_SUCCESS)
    {
        return ROM_STDD_ERR_IO;
    }

    written = snprintf(
        filename,
        sizeof(filename),
        "%s/fluid_st_ddrom_window_%04d_rank_%06d.bin",
        directory,
        window_id,
        rank);

    if(written < 0 ||
       (size_t)written >= sizeof(filename))
    {
        return ROM_STDD_ERR_IO;
    }

    return ST_fom_write_owned_window_binary(
        filename,
        window_id,
        time_start,
        time_end,
        context->layout,
        context->current_window) == ST_FOM_SUCCESS
        ? ROM_STDD_SUCCESS
        : ROM_STDD_ERR_IO;
}

int ROM_std_stdd_fluid_owned_window_error(
    ROM_STDD_FLUID_CONTEXT* context,
    const double* reference_owned_window,
    double* relative_error)
{
    double difference2 = 0.0;
    double reference2 = 0.0;
    double approximation2 = 0.0;

    if(fluid_rom_validate(context) != ROM_STDD_SUCCESS ||
       reference_owned_window == NULL ||
       relative_error == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ST_fom_owned_difference_squares(
        context->layout,
        reference_owned_window,
        context->current_window,
        &difference2,
        &reference2,
        &approximation2) != ST_FOM_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    monolis_allreduce_R(
        1,
        &difference2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);
    monolis_allreduce_R(
        1,
        &reference2,
        MONOLIS_MPI_SUM,
        context->fom->mono_com.comm);

    *relative_error =
        sqrt(fmax(difference2, 0.0))
        /
        fmax(
            sqrt(fmax(reference2, 0.0)),
            1.0e-300);

    (void)approximation2;
    return ROM_STDD_SUCCESS;
}
