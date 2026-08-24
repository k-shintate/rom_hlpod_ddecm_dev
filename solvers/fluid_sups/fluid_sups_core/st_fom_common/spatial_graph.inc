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

