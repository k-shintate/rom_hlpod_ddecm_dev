/*
 * Causal dG(0) window runner and binary/diagnostic I/O helpers.
 */

static int st_fom_layout_is_valid_for_io(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    return layout != NULL &&
        layout->num_space_points > 0 &&
        layout->num_owned_space_points > 0 &&
        layout->num_owned_space_points <= layout->num_space_points &&
        layout->physical_dof > 0 &&
        layout->num_slabs > 0 &&
        layout->time.n_dof > 0 &&
        layout->time.n_dof <= ST_FOM_MAX_TIME_DOF &&
        ST_fom_window_dof_per_space_point(layout) > 0 &&
        ST_fom_window_vector_length(layout) > 0;
}

int ST_fom_run_causal_dg0_windows(
    const ST_FOM_CAUSAL_OPTIONS* options,
    const ST_FOM_CAUSAL_CALLBACKS* callbacks)
{
    const ST_FOM_WINDOW_LAYOUT* layout;
    size_t window_length;
    size_t state_length;
    double* window_state = NULL;
    double* state = NULL;
    double* right_trace = NULL;
    int global_step = 0;
    int file_number = 0;
    int status = ST_FOM_SUCCESS;
    int window;
    int slab;

    if(options == NULL || callbacks == NULL ||
       callbacks->solve_step == NULL || callbacks->export_state == NULL ||
       options->num_windows <= 0 || !isfinite(options->dt) || options->dt <= 0.0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    layout = &options->layout;
    if(!st_fom_layout_is_valid_for_io(layout) ||
       layout->time.family != ST_FOM_TIME_FAMILY_DG ||
       layout->time.degree != 0 || layout->time.n_dof != 1)
    {
        return ST_FOM_ERR_UNSUPPORTED;
    }

    if(options->write_output){
        if(callbacks->write_output == NULL || options->output_interval <= 0){
            return ST_FOM_ERR_ARGUMENT;
        }
    }

    window_length = ST_fom_window_vector_length(layout);
    state_length = (size_t)layout->num_space_points
        * (size_t)layout->physical_dof;

    window_state = (double*)calloc(window_length, sizeof(double));
    state = (double*)calloc(state_length, sizeof(double));
    right_trace = (double*)calloc(state_length, sizeof(double));
    if(window_state == NULL || state == NULL || right_trace == NULL){
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    for(window = 0; window < options->num_windows; window++){
        const double time_start =
            (double)(window * layout->num_slabs) * options->dt;
        const double time_end =
            (double)((window + 1) * layout->num_slabs) * options->dt;

        memset(window_state, 0, window_length * sizeof(double));

        for(slab = 0; slab < layout->num_slabs; slab++){
            double time_old;
            double time_new;
            int callback_status;

            global_step++;
            time_old = (double)(global_step - 1) * options->dt;
            time_new = (double)global_step * options->dt;

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
                layout,
                slab,
                state,
                window_state);
            if(status != ST_FOM_SUCCESS){
                goto cleanup;
            }

            if(options->write_output &&
               (global_step % options->output_interval) == 0)
            {
                file_number++;
                callback_status = callbacks->write_output(
                    callbacks->user_context,
                    file_number,
                    global_step,
                    time_new);
                if(callback_status != 0){
                    status = ST_FOM_ERR_CALLBACK;
                    goto cleanup;
                }
            }
        }

        status = ST_fom_extract_last_right_trace(
            layout,
            window_state,
            right_trace);
        if(status != ST_FOM_SUCCESS){
            goto cleanup;
        }

        if(callbacks->window_complete != NULL){
            const int callback_status = callbacks->window_complete(
                callbacks->user_context,
                window,
                time_start,
                time_end,
                layout,
                window_state,
                right_trace);
            if(callback_status != 0){
                status = ST_FOM_ERR_CALLBACK;
                goto cleanup;
            }
        }
    }

cleanup:
    free(window_state);
    free(state);
    free(right_trace);
    return status;
}

int ST_fom_make_directory_recursive(const char* directory)
{
    char* path = NULL;
    size_t length;
    size_t i;

    if(directory == NULL || directory[0] == '\0'){
        return ST_FOM_ERR_ARGUMENT;
    }

    length = strlen(directory);
    path = (char*)malloc(length + 1u);
    if(path == NULL){
        return ST_FOM_ERR_ALLOCATION;
    }
    memcpy(path, directory, length + 1u);

    /* Remove trailing slashes except for root. */
    while(length > 1u && path[length - 1u] == '/'){
        path[--length] = '\0';
    }

    for(i = 1u; i < length; i++){
        if(path[i] != '/'){
            continue;
        }

        path[i] = '\0';
        if(path[0] != '\0'){
            if(mkdir(path, 0777) != 0 && errno != EEXIST){
                free(path);
                return ST_FOM_ERR_IO;
            }
        }
        path[i] = '/';
    }

    if(mkdir(path, 0777) != 0 && errno != EEXIST){
        free(path);
        return ST_FOM_ERR_IO;
    }

    free(path);
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
    FILE* fp;
    const int dof_per_space_point =
        ST_fom_window_dof_per_space_point(layout);
    size_t owned_length;

    if(filename == NULL || filename[0] == '\0' ||
       window_id < 0 || !isfinite(time_start) || !isfinite(time_end) ||
       time_end < time_start || !st_fom_layout_is_valid_for_io(layout) ||
       window_state == NULL || dof_per_space_point <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    owned_length = (size_t)layout->num_owned_space_points
        * (size_t)dof_per_space_point;

    memset(&header, 0, sizeof(header));
    header.magic = (uint32_t)ST_FOM_WINDOW_MAGIC;
    header.version = (uint32_t)ST_FOM_WINDOW_VERSION;
    header.window_id = (uint32_t)window_id;
    header.degree = (uint32_t)layout->time.degree;
    header.num_slabs = (uint32_t)layout->num_slabs;
    header.num_time_dofs = (uint32_t)layout->time.n_dof;
    header.physical_dof = (uint32_t)layout->physical_dof;
    header.num_local_space_points = (uint32_t)layout->num_space_points;
    header.num_owned_space_points =
        (uint32_t)layout->num_owned_space_points;
    header.time_start = time_start;
    header.time_end = time_end;

    fp = fopen(filename, "wb");
    if(fp == NULL){
        return ST_FOM_ERR_IO;
    }

    if(fwrite(&header, sizeof(header), 1u, fp) != 1u ||
       (owned_length > 0 &&
        fwrite(window_state, sizeof(double), owned_length, fp) != owned_length))
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
    const int dof_per_space_point =
        ST_fom_window_dof_per_space_point(layout);
    size_t count;
    size_t i;
    double diff2 = 0.0;
    double ref2 = 0.0;
    double approx2 = 0.0;

    if(!st_fom_layout_is_valid_for_io(layout) || reference == NULL ||
       approximation == NULL || difference_norm2 == NULL ||
       reference_norm2 == NULL || approximation_norm2 == NULL ||
       dof_per_space_point <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    count = (size_t)layout->num_owned_space_points
        * (size_t)dof_per_space_point;

    for(i = 0; i < count; i++){
        const double ref = reference[i];
        const double approx = approximation[i];
        const double diff = approx - ref;

        diff2 += diff * diff;
        ref2 += ref * ref;
        approx2 += approx * approx;
    }

    *difference_norm2 = diff2;
    *reference_norm2 = ref2;
    *approximation_norm2 = approx2;

    return ST_FOM_SUCCESS;
}
