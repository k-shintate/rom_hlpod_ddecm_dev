int ST_fom_run_causal_dg0_windows(
    const ST_FOM_CAUSAL_OPTIONS* options,
    const ST_FOM_CAUSAL_CALLBACKS* callbacks)
{
    size_t state_length;
    size_t window_length;
    double* state = NULL;
    double* window_state = NULL;
    double* right_trace = NULL;
    int window;
    int slab;
    int file_number = 0;
    int global_step = 0;
    int status = ST_FOM_SUCCESS;

    if(options == NULL || callbacks == NULL ||
       callbacks->solve_step == NULL || callbacks->export_state == NULL ||
       options->num_windows <= 0 || options->dt <= 0.0 ||
       options->layout.time.degree != 0 ||
       options->layout.time.n_dof != 1)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    state_length = (size_t)options->layout.num_space_points
        * (size_t)options->layout.physical_dof;
    window_length = ST_fom_window_vector_length(&options->layout);
    if(state_length == 0u || window_length == 0u){
        return ST_FOM_ERR_DIMENSION;
    }

    state = (double*)calloc(state_length, sizeof(double));
    window_state = (double*)calloc(window_length, sizeof(double));
    right_trace = (double*)calloc(state_length, sizeof(double));
    if(state == NULL || window_state == NULL || right_trace == NULL){
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    for(window = 0; window < options->num_windows; window++){
        const double time_start =
            (double)(window * options->layout.num_slabs) * options->dt;
        const double time_end = time_start
            + (double)options->layout.num_slabs * options->dt;

        memset(window_state, 0, window_length * sizeof(double));

        for(slab = 0; slab < options->layout.num_slabs; slab++){
            const double time_old = time_start + (double)slab * options->dt;
            const double time_new = time_old + options->dt;
            int callback_status;

            global_step++;
            callback_status = callbacks->solve_step(
                callbacks->user_context,
                global_step,
                time_old,
                time_new);
            if(callback_status != 0){
                status = ST_FOM_ERR_CALLBACK;
                goto cleanup;
            }

            callback_status = callbacks->export_state(
                callbacks->user_context,
                state);
            if(callback_status != 0){
                status = ST_FOM_ERR_CALLBACK;
                goto cleanup;
            }

            status = ST_fom_pack_dg0_slab(
                &options->layout, slab, state, window_state);
            if(status != ST_FOM_SUCCESS){
                goto cleanup;
            }

            if(options->write_output && callbacks->write_output != NULL &&
               options->output_interval > 0 &&
               global_step % options->output_interval == 0)
            {
                callback_status = callbacks->write_output(
                    callbacks->user_context,
                    file_number,
                    global_step,
                    time_new);
                if(callback_status != 0){
                    status = ST_FOM_ERR_CALLBACK;
                    goto cleanup;
                }
                file_number++;
            }
        }

        status = ST_fom_extract_last_right_trace(
            &options->layout, window_state, right_trace);
        if(status != ST_FOM_SUCCESS){
            goto cleanup;
        }

        if(callbacks->window_complete != NULL){
            const int callback_status = callbacks->window_complete(
                callbacks->user_context,
                window,
                time_start,
                time_end,
                &options->layout,
                window_state,
                right_trace);
            if(callback_status != 0){
                status = ST_FOM_ERR_CALLBACK;
                goto cleanup;
            }
        }
    }

cleanup:
    free(state);
    free(window_state);
    free(right_trace);
    return status;
}

int ST_fom_make_directory_recursive(const char* directory)
{
    char* work;
    char* p;
    size_t length;

    if(directory == NULL || directory[0] == '\0'){
        return ST_FOM_ERR_ARGUMENT;
    }

    length = strlen(directory);
    work = (char*)malloc(length + 1u);
    if(work == NULL){
        return ST_FOM_ERR_ALLOCATION;
    }
    memcpy(work, directory, length + 1u);

    if(length > 1u && work[length - 1u] == '/'){
        work[length - 1u] = '\0';
    }

    for(p = work + 1; *p != '\0'; p++){
        if(*p == '/'){
            *p = '\0';
            if(mkdir(work, 0775) != 0 && errno != EEXIST){
                free(work);
                return ST_FOM_ERR_IO;
            }
            *p = '/';
        }
    }

    if(mkdir(work, 0775) != 0 && errno != EEXIST){
        free(work);
        return ST_FOM_ERR_IO;
    }

    free(work);
    return ST_FOM_SUCCESS;
}

int ST_fom_write_owned_window_binary(
    const char* filename,
    int window_id,
    double time_start,
    double time_end,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state)
{
    ST_FOM_WINDOW_FILE_HEADER header;
    const int block_dof = ST_fom_window_dof_per_space_point(layout);
    const size_t count = (layout != NULL && block_dof > 0)
        ? (size_t)layout->num_owned_space_points * (size_t)block_dof
        : 0u;
    FILE* fp;

    if(filename == NULL || layout == NULL || window_state == NULL ||
       window_id < 0 || count == 0u)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    memset(&header, 0, sizeof(header));
    header.magic = ST_FOM_WINDOW_MAGIC;
    header.version = ST_FOM_WINDOW_VERSION;
    header.window_id = (uint32_t)window_id;
    header.degree = (uint32_t)layout->time.degree;
    header.num_slabs = (uint32_t)layout->num_slabs;
    header.num_time_dofs = (uint32_t)layout->time.n_dof;
    header.physical_dof = (uint32_t)layout->physical_dof;
    header.num_local_space_points = (uint32_t)layout->num_space_points;
    header.num_owned_space_points = (uint32_t)layout->num_owned_space_points;
    header.time_start = time_start;
    header.time_end = time_end;

    fp = fopen(filename, "wb");
    if(fp == NULL){
        return ST_FOM_ERR_IO;
    }

    if(fwrite(&header, sizeof(header), 1u, fp) != 1u ||
       fwrite(window_state, sizeof(double), count, fp) != count)
    {
        fclose(fp);
        return ST_FOM_ERR_IO;
    }

    if(fclose(fp) != 0){
        return ST_FOM_ERR_IO;
    }
    return ST_FOM_SUCCESS;
}

int ST_fom_owned_difference_squares(
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* reference,
    const double* approximation,
    double* difference_norm2,
    double* reference_norm2,
    double* approximation_norm2)
{
    const int block_dof = ST_fom_window_dof_per_space_point(layout);
    int node;
    int q;
    double diff2 = 0.0;
    double ref2 = 0.0;
    double app2 = 0.0;

    if(layout == NULL || reference == NULL || approximation == NULL ||
       difference_norm2 == NULL || reference_norm2 == NULL ||
       approximation_norm2 == NULL || block_dof <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_owned_space_points; node++){
        for(q = 0; q < block_dof; q++){
            const size_t index = (size_t)node * (size_t)block_dof
                + (size_t)q;
            const double difference = reference[index] - approximation[index];
            diff2 += difference * difference;
            ref2 += reference[index] * reference[index];
            app2 += approximation[index] * approximation[index];
        }
    }

    *difference_norm2 = diff2;
    *reference_norm2 = ref2;
    *approximation_norm2 = app2;
    return ST_FOM_SUCCESS;
}
