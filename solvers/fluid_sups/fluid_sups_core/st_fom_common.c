#include "st_fom_common.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


typedef struct {
    int* data;
    int size;
    int capacity;
} ST_FOM_INT_VECTOR;

static int ST_fom_int_compare(const void* left, const void* right)
{
    const int a = *(const int*)left;
    const int b = *(const int*)right;
    return (a > b) - (a < b);
}

static int ST_fom_int_vector_push(
    ST_FOM_INT_VECTOR* vector,
    int value)
{
    int* resized;
    int new_capacity;

    if(vector == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }

    if(vector->size < vector->capacity){
        vector->data[vector->size++] = value;
        return ST_FOM_SUCCESS;
    }

    new_capacity = vector->capacity > 0
        ? 2 * vector->capacity
        : 8;

    if(new_capacity <= vector->capacity ||
       (size_t)new_capacity > SIZE_MAX / sizeof(int))
    {
        return ST_FOM_ERR_DIMENSION;
    }

    resized = (int*)realloc(
        vector->data,
        (size_t)new_capacity * sizeof(int));
    if(resized == NULL){
        return ST_FOM_ERR_ALLOCATION;
    }

    vector->data = resized;
    vector->capacity = new_capacity;
    vector->data[vector->size++] = value;
    return ST_FOM_SUCCESS;
}

static void ST_fom_int_vector_sort_unique_without_self(
    ST_FOM_INT_VECTOR* vector,
    int self)
{
    int read;
    int write = 0;

    if(vector == NULL || vector->size <= 0){
        return;
    }

    qsort(
        vector->data,
        (size_t)vector->size,
        sizeof(int),
        ST_fom_int_compare);

    for(read = 0; read < vector->size; read++){
        const int value = vector->data[read];

        if(value == self){
            continue;
        }

        if(write == 0 || value != vector->data[write - 1]){
            vector->data[write++] = value;
        }
    }

    vector->size = write;
}

static int ST_fom_graph_row_contains(
    const ST_FOM_SPATIAL_GRAPH* graph,
    int row,
    int column)
{
    int left;
    int right;

    if(graph == NULL || graph->index == NULL || graph->item == NULL ||
       row < 0 || row >= graph->num_nodes)
    {
        return 0;
    }

    left = graph->index[row];
    right = graph->index[row + 1];

    while(left < right){
        const int middle = left + (right - left) / 2;
        const int value = graph->item[middle];

        if(value == column){
            return 1;
        }
        if(value < column){
            left = middle + 1;
        }
        else{
            right = middle;
        }
    }

    return 0;
}

int ST_fom_spatial_graph_initialize(
    ST_FOM_SPATIAL_GRAPH* graph)
{
    if(graph == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }

    memset(graph, 0, sizeof(*graph));
    return ST_FOM_SUCCESS;
}

void ST_fom_spatial_graph_finalize(
    ST_FOM_SPATIAL_GRAPH* graph)
{
    if(graph == NULL){
        return;
    }

    free(graph->index);
    free(graph->item);
    memset(graph, 0, sizeof(*graph));
}

int ST_fom_spatial_graph_build_from_connectivity(
    ST_FOM_SPATIAL_GRAPH* graph,
    int num_nodes,
    int num_elements,
    int nodes_per_element,
    int** connectivity)
{
    ST_FOM_INT_VECTOR* rows = NULL;
    long long num_items64 = 0;
    int element;
    int local_row;
    int local_column;
    int node;
    int position = 0;
    int status = ST_FOM_SUCCESS;

    if(graph == NULL || num_nodes <= 0 ||
       num_elements < 0 || nodes_per_element <= 0 ||
       (num_elements > 0 && connectivity == NULL))
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    if(graph->index != NULL || graph->item != NULL){
        return ST_FOM_ERR_ARGUMENT;
    }

    rows = (ST_FOM_INT_VECTOR*)calloc(
        (size_t)num_nodes,
        sizeof(*rows));
    if(rows == NULL){
        return ST_FOM_ERR_ALLOCATION;
    }

    for(element = 0; element < num_elements; element++){
        if(connectivity[element] == NULL){
            status = ST_FOM_ERR_ARGUMENT;
            goto cleanup;
        }

        for(local_row = 0;
            local_row < nodes_per_element;
            local_row++)
        {
            const int row = connectivity[element][local_row];

            if(row < 0 || row >= num_nodes){
                status = ST_FOM_ERR_DIMENSION;
                goto cleanup;
            }

            for(local_column = 0;
                local_column < nodes_per_element;
                local_column++)
            {
                const int column =
                    connectivity[element][local_column];

                if(column < 0 || column >= num_nodes){
                    status = ST_FOM_ERR_DIMENSION;
                    goto cleanup;
                }

                if(column != row){
                    status = ST_fom_int_vector_push(
                        &rows[row], column);
                    if(status != ST_FOM_SUCCESS){
                        goto cleanup;
                    }
                }
            }
        }
    }

    for(node = 0; node < num_nodes; node++){
        ST_fom_int_vector_sort_unique_without_self(
            &rows[node], node);

        num_items64 += rows[node].size;
        if(num_items64 > INT_MAX){
            status = ST_FOM_ERR_DIMENSION;
            goto cleanup;
        }
    }

    graph->index = (int*)calloc(
        (size_t)num_nodes + 1u,
        sizeof(int));
    graph->item = (int*)calloc(
        (size_t)num_items64,
        sizeof(int));

    if(graph->index == NULL ||
       (num_items64 > 0 && graph->item == NULL))
    {
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    graph->num_nodes = num_nodes;
    graph->num_items = (int)num_items64;

    for(node = 0; node < num_nodes; node++){
        int k;

        graph->index[node] = position;
        for(k = 0; k < rows[node].size; k++){
            graph->item[position++] = rows[node].data[k];
        }
    }
    graph->index[num_nodes] = position;

    if(position != graph->num_items){
        status = ST_FOM_ERR_DIMENSION;
        goto cleanup;
    }

cleanup:
    if(rows != NULL){
        for(node = 0; node < num_nodes; node++){
            free(rows[node].data);
        }
        free(rows);
    }

    if(status != ST_FOM_SUCCESS){
        ST_fom_spatial_graph_finalize(graph);
    }

    return status;
}

int ST_fom_spatial_graph_validate_connectivity(
    const ST_FOM_SPATIAL_GRAPH* graph,
    int num_nodes,
    int num_elements,
    int nodes_per_element,
    int** connectivity)
{
    int element;
    int local_row;
    int local_column;

    if(graph == NULL || graph->index == NULL ||
       graph->num_nodes != num_nodes ||
       num_nodes <= 0 || num_elements < 0 ||
       nodes_per_element <= 0 ||
       (num_elements > 0 && connectivity == NULL))
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(element = 0; element < num_elements; element++){
        if(connectivity[element] == NULL){
            return ST_FOM_ERR_ARGUMENT;
        }

        for(local_row = 0;
            local_row < nodes_per_element;
            local_row++)
        {
            const int row = connectivity[element][local_row];

            if(row < 0 || row >= num_nodes){
                return ST_FOM_ERR_DIMENSION;
            }

            for(local_column = 0;
                local_column < nodes_per_element;
                local_column++)
            {
                const int column =
                    connectivity[element][local_column];

                if(column < 0 || column >= num_nodes){
                    return ST_FOM_ERR_DIMENSION;
                }

                if(row != column &&
                   !ST_fom_graph_row_contains(
                       graph, row, column))
                {
                    return ST_FOM_ERR_DIMENSION;
                }
            }
        }
    }

    return ST_FOM_SUCCESS;
}

static void ST_fom_time_clear(ST_FOM_TIME_ELEMENT* time)
{
    if(time != NULL){
        memset(time, 0, sizeof(*time));
        time->family = ST_FOM_TIME_FAMILY_DG;
    }
}

static void ST_fom_set_flux_matrices(ST_FOM_TIME_ELEMENT* time)
{
    int alpha;
    int beta;

    for(alpha = 0; alpha < time->n_dof; alpha++){
        for(beta = 0; beta < time->n_dof; beta++){
            time->E_self[alpha][beta] =
                time->e_minus[alpha] * time->e_minus[beta];
            time->E_prev[alpha][beta] =
                time->e_minus[alpha] * time->e_plus[beta];
        }
    }
}

int ST_fom_time_element_init_dg(
    ST_FOM_TIME_ELEMENT* time,
    int degree)
{
    int q;

    if(time == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }

    ST_fom_time_clear(time);
    time->degree = degree;
    time->is_nodal = 1;

    if(degree == 0){
        time->n_dof = 1;
        time->dof_tau[0] = 1.0;
        time->Mt[0][0] = 1.0;
        time->Dt[0][0] = 0.0;
        time->e_minus[0] = 1.0;
        time->e_plus[0] = 1.0;
        time->num_time_ip = 1;
        time->tau[0] = 1.0;
        time->weight[0] = 1.0;
        time->psi[0][0] = 1.0;
    }
    else if(degree == 1){
        const double s = 0.57735026918962576451;

        time->n_dof = 2;
        time->dof_tau[0] = 0.0;
        time->dof_tau[1] = 1.0;

        time->Mt[0][0] = 1.0 / 3.0;
        time->Mt[0][1] = 1.0 / 6.0;
        time->Mt[1][0] = 1.0 / 6.0;
        time->Mt[1][1] = 1.0 / 3.0;

        time->Dt[0][0] = -0.5;
        time->Dt[0][1] =  0.5;
        time->Dt[1][0] = -0.5;
        time->Dt[1][1] =  0.5;

        time->e_minus[0] = 1.0;
        time->e_minus[1] = 0.0;
        time->e_plus[0] = 0.0;
        time->e_plus[1] = 1.0;

        time->num_time_ip = 2;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5 * (1.0 + s);
        time->weight[0] = 0.5;
        time->weight[1] = 0.5;
        for(q = 0; q < time->num_time_ip; q++){
            time->psi[0][q] = 1.0 - time->tau[q];
            time->psi[1][q] = time->tau[q];
        }
    }
    else if(degree == 2){
        const double s = 0.77459666924148337704;

        time->n_dof = 3;
        time->dof_tau[0] = 0.0;
        time->dof_tau[1] = 0.5;
        time->dof_tau[2] = 1.0;

        time->Mt[0][0] =  2.0 / 15.0;
        time->Mt[0][1] =  1.0 / 15.0;
        time->Mt[0][2] = -1.0 / 30.0;
        time->Mt[1][0] =  1.0 / 15.0;
        time->Mt[1][1] =  8.0 / 15.0;
        time->Mt[1][2] =  1.0 / 15.0;
        time->Mt[2][0] = -1.0 / 30.0;
        time->Mt[2][1] =  1.0 / 15.0;
        time->Mt[2][2] =  2.0 / 15.0;

        time->Dt[0][0] = -1.0 / 2.0;
        time->Dt[0][1] =  2.0 / 3.0;
        time->Dt[0][2] = -1.0 / 6.0;
        time->Dt[1][0] = -2.0 / 3.0;
        time->Dt[1][1] =  0.0;
        time->Dt[1][2] =  2.0 / 3.0;
        time->Dt[2][0] =  1.0 / 6.0;
        time->Dt[2][1] = -2.0 / 3.0;
        time->Dt[2][2] =  1.0 / 2.0;

        time->e_minus[0] = 1.0;
        time->e_minus[1] = 0.0;
        time->e_minus[2] = 0.0;
        time->e_plus[0] = 0.0;
        time->e_plus[1] = 0.0;
        time->e_plus[2] = 1.0;

        time->num_time_ip = 3;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5;
        time->tau[2] = 0.5 * (1.0 + s);
        time->weight[0] = 5.0 / 18.0;
        time->weight[1] = 4.0 / 9.0;
        time->weight[2] = 5.0 / 18.0;
        for(q = 0; q < time->num_time_ip; q++){
            const double tau = time->tau[q];
            time->psi[0][q] = 2.0 * tau * tau - 3.0 * tau + 1.0;
            time->psi[1][q] = -4.0 * tau * tau + 4.0 * tau;
            time->psi[2][q] = 2.0 * tau * tau - tau;
        }
    }
    else{
        ST_fom_time_clear(time);
        return ST_FOM_ERR_UNSUPPORTED;
    }

    ST_fom_set_flux_matrices(time);
    return ST_FOM_SUCCESS;
}

int ST_fom_resolve_owned_space_points(
    int num_local_space_points,
    int communicator_size,
    int communication_owned_space_points,
    int* num_owned_space_points)
{
    if(num_owned_space_points == NULL ||
       num_local_space_points <= 0 || communicator_size <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    if(communicator_size == 1){
        *num_owned_space_points = num_local_space_points;
        return ST_FOM_SUCCESS;
    }

    if(communication_owned_space_points <= 0 ||
       communication_owned_space_points > num_local_space_points)
    {
        return ST_FOM_ERR_DIMENSION;
    }

    *num_owned_space_points = communication_owned_space_points;
    return ST_FOM_SUCCESS;
}

int ST_fom_layout_initialize(
    ST_FOM_WINDOW_LAYOUT* layout,
    int num_space_points,
    int num_owned_space_points,
    int physical_dof,
    int num_slabs,
    int degree)
{
    int status;

    if(layout == NULL || num_space_points <= 0 ||
       num_owned_space_points <= 0 ||
       num_owned_space_points > num_space_points ||
       physical_dof <= 0 || num_slabs <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    memset(layout, 0, sizeof(*layout));
    layout->num_space_points = num_space_points;
    layout->num_owned_space_points = num_owned_space_points;
    layout->physical_dof = physical_dof;
    layout->num_slabs = num_slabs;

    status = ST_fom_time_element_init_dg(&layout->time, degree);
    if(status != ST_FOM_SUCCESS){
        memset(layout, 0, sizeof(*layout));
        return status;
    }

    if(ST_fom_window_vector_length(layout) == 0u){
        memset(layout, 0, sizeof(*layout));
        return ST_FOM_ERR_DIMENSION;
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_window_dof_per_space_point(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    long long value;

    if(layout == NULL || layout->physical_dof <= 0 ||
       layout->num_slabs <= 0 || layout->time.n_dof <= 0)
    {
        return 0;
    }

    value = (long long)layout->physical_dof
        * (long long)layout->num_slabs
        * (long long)layout->time.n_dof;
    if(value <= 0 || value > INT_MAX){
        return 0;
    }
    return (int)value;
}

size_t ST_fom_window_vector_length(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    const int block_dof = ST_fom_window_dof_per_space_point(layout);

    if(layout == NULL || block_dof <= 0 || layout->num_space_points <= 0){
        return 0u;
    }
    if((size_t)layout->num_space_points > SIZE_MAX / (size_t)block_dof){
        return 0u;
    }
    return (size_t)layout->num_space_points * (size_t)block_dof;
}

int ST_fom_block_dof(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    int time_dof,
    int component)
{
    long long value;

    if(layout == NULL || slab < 0 || slab >= layout->num_slabs ||
       time_dof < 0 || time_dof >= layout->time.n_dof ||
       component < 0 || component >= layout->physical_dof)
    {
        return -1;
    }

    value = ((long long)slab * (long long)layout->time.n_dof
        + (long long)time_dof) * (long long)layout->physical_dof
        + (long long)component;
    if(value < 0 || value > INT_MAX){
        return -1;
    }
    return (int)value;
}

size_t ST_fom_vector_index(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int space_point,
    int slab,
    int time_dof,
    int component)
{
    const int block_dof = ST_fom_window_dof_per_space_point(layout);
    const int local_dof = ST_fom_block_dof(
        layout, slab, time_dof, component);

    if(layout == NULL || block_dof <= 0 || local_dof < 0 ||
       space_point < 0 || space_point >= layout->num_space_points)
    {
        return SIZE_MAX;
    }

    if((size_t)space_point >
       (SIZE_MAX - (size_t)local_dof) / (size_t)block_dof)
    {
        return SIZE_MAX;
    }

    return (size_t)space_point * (size_t)block_dof
        + (size_t)local_dof;
}

int ST_fom_pack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* state,
    double* window_state)
{
    int node;
    int component;

    if(layout == NULL || state == NULL || window_state == NULL ||
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t index = ST_fom_vector_index(
                layout, node, slab, 0, component);
            if(index == SIZE_MAX){
                return ST_FOM_ERR_DIMENSION;
            }
            window_state[index] = state[
                (size_t)node * (size_t)layout->physical_dof
                + (size_t)component];
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_unpack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* state)
{
    int node;
    int component;

    if(layout == NULL || state == NULL || window_state == NULL ||
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t index = ST_fom_vector_index(
                layout, node, slab, 0, component);
            if(index == SIZE_MAX){
                return ST_FOM_ERR_DIMENSION;
            }
            state[(size_t)node * (size_t)layout->physical_dof
                + (size_t)component] = window_state[index];
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_extract_slab_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* trace)
{
    int node;
    int component;
    int alpha;

    if(layout == NULL || window_state == NULL || trace == NULL ||
       slab < 0 || slab >= layout->num_slabs ||
       layout->time.n_dof <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            double value = 0.0;
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                const size_t index = ST_fom_vector_index(
                    layout, node, slab, alpha, component);
                if(index == SIZE_MAX){
                    return ST_FOM_ERR_DIMENSION;
                }
                value += layout->time.e_plus[alpha] * window_state[index];
            }
            trace[(size_t)node * (size_t)layout->physical_dof
                + (size_t)component] = value;
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_extract_last_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    double* trace)
{
    if(layout == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }
    return ST_fom_extract_slab_right_trace(
        layout, layout->num_slabs - 1, window_state, trace);
}

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
