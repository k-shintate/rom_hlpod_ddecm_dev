/*
 * Spatial nodal graph utilities for st_fom_common.c.
 *
 * The graph is stored in CSR form.  For each finite element, every pair of
 * nodes in that element is inserted as a possible matrix coupling.  Self
 * entries are retained deliberately: this makes the graph represent the full
 * element sparsity pattern and is safe for sparse-pattern construction.
 */

static int st_fom_int_compare(const void* lhs, const void* rhs)
{
    const int a = *(const int*)lhs;
    const int b = *(const int*)rhs;
    return (a > b) - (a < b);
}

static int st_fom_graph_contains(
    const ST_FOM_SPATIAL_GRAPH* graph,
    int row,
    int col)
{
    int k;

    if(graph == NULL || graph->index == NULL || graph->item == NULL ||
       row < 0 || row >= graph->num_nodes)
    {
        return 0;
    }

    for(k = graph->index[row]; k < graph->index[row + 1]; k++){
        if(graph->item[k] == col){
            return 1;
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

    graph->num_nodes = 0;
    graph->num_items = 0;
    graph->index = NULL;
    graph->item = NULL;

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

    graph->num_nodes = 0;
    graph->num_items = 0;
    graph->index = NULL;
    graph->item = NULL;
}

int ST_fom_spatial_graph_build_from_connectivity(
    ST_FOM_SPATIAL_GRAPH* graph,
    int num_nodes,
    int num_elements,
    int nodes_per_element,
    int** connectivity)
{
    int** row_items = NULL;
    int* row_count = NULL;
    int* row_capacity = NULL;
    size_t total_items = 0;
    int status = ST_FOM_SUCCESS;
    int e;
    int i;
    int j;
    int row;

    if(graph == NULL || num_nodes <= 0 || num_elements < 0 ||
       nodes_per_element <= 0 ||
       (num_elements > 0 && connectivity == NULL))
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    /* Permit rebuilding an already initialized graph. */
    ST_fom_spatial_graph_finalize(graph);

    row_items = (int**)calloc((size_t)num_nodes, sizeof(int*));
    row_count = (int*)calloc((size_t)num_nodes, sizeof(int));
    row_capacity = (int*)calloc((size_t)num_nodes, sizeof(int));
    if(row_items == NULL || row_count == NULL || row_capacity == NULL){
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    for(e = 0; e < num_elements; e++){
        if(connectivity[e] == NULL){
            status = ST_FOM_ERR_ARGUMENT;
            goto cleanup;
        }

        for(i = 0; i < nodes_per_element; i++){
            const int row_node = connectivity[e][i];

            if(row_node < 0 || row_node >= num_nodes){
                status = ST_FOM_ERR_DIMENSION;
                goto cleanup;
            }

            for(j = 0; j < nodes_per_element; j++){
                const int col_node = connectivity[e][j];
                int exists = 0;
                int k;

                if(col_node < 0 || col_node >= num_nodes){
                    status = ST_FOM_ERR_DIMENSION;
                    goto cleanup;
                }

                for(k = 0; k < row_count[row_node]; k++){
                    if(row_items[row_node][k] == col_node){
                        exists = 1;
                        break;
                    }
                }

                if(!exists){
                    if(row_count[row_node] == row_capacity[row_node]){
                        int new_capacity = row_capacity[row_node] == 0
                            ? nodes_per_element
                            : 2 * row_capacity[row_node];
                        int* grown;

                        if(new_capacity <= row_capacity[row_node]){
                            status = ST_FOM_ERR_ALLOCATION;
                            goto cleanup;
                        }

                        grown = (int*)realloc(
                            row_items[row_node],
                            (size_t)new_capacity * sizeof(int));
                        if(grown == NULL){
                            status = ST_FOM_ERR_ALLOCATION;
                            goto cleanup;
                        }

                        row_items[row_node] = grown;
                        row_capacity[row_node] = new_capacity;
                    }

                    row_items[row_node][row_count[row_node]++] = col_node;
                }
            }
        }
    }

    for(row = 0; row < num_nodes; row++){
        if(row_count[row] > 1){
            qsort(
                row_items[row],
                (size_t)row_count[row],
                sizeof(int),
                st_fom_int_compare);
        }

        total_items += (size_t)row_count[row];
        if(total_items > (size_t)INT_MAX){
            status = ST_FOM_ERR_DIMENSION;
            goto cleanup;
        }
    }

    graph->index = (int*)calloc((size_t)num_nodes + 1u, sizeof(int));
    if(graph->index == NULL){
        status = ST_FOM_ERR_ALLOCATION;
        goto cleanup;
    }

    if(total_items > 0){
        graph->item = (int*)malloc(total_items * sizeof(int));
        if(graph->item == NULL){
            status = ST_FOM_ERR_ALLOCATION;
            goto cleanup;
        }
    }

    graph->num_nodes = num_nodes;
    graph->num_items = (int)total_items;

    for(row = 0; row < num_nodes; row++){
        int k;
        graph->index[row + 1] = graph->index[row] + row_count[row];
        for(k = 0; k < row_count[row]; k++){
            graph->item[graph->index[row] + k] = row_items[row][k];
        }
    }

cleanup:
    if(row_items != NULL){
        for(row = 0; row < num_nodes; row++){
            free(row_items[row]);
        }
    }
    free(row_items);
    free(row_count);
    free(row_capacity);

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
    int e;
    int i;
    int j;

    if(graph == NULL || num_nodes <= 0 || num_elements < 0 ||
       nodes_per_element <= 0 ||
       (num_elements > 0 && connectivity == NULL))
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    if(graph->num_nodes != num_nodes || graph->index == NULL ||
       graph->num_items < 0 || graph->index[0] != 0 ||
       graph->index[num_nodes] != graph->num_items ||
       (graph->num_items > 0 && graph->item == NULL))
    {
        return ST_FOM_ERR_DIMENSION;
    }

    for(i = 0; i < num_nodes; i++){
        int k;
        if(graph->index[i] > graph->index[i + 1] ||
           graph->index[i] < 0 ||
           graph->index[i + 1] > graph->num_items)
        {
            return ST_FOM_ERR_DIMENSION;
        }
        for(k = graph->index[i]; k < graph->index[i + 1]; k++){
            if(graph->item[k] < 0 || graph->item[k] >= num_nodes){
                return ST_FOM_ERR_DIMENSION;
            }
            if(k > graph->index[i] && graph->item[k - 1] >= graph->item[k]){
                return ST_FOM_ERR_DIMENSION;
            }
        }
    }

    for(e = 0; e < num_elements; e++){
        if(connectivity[e] == NULL){
            return ST_FOM_ERR_ARGUMENT;
        }

        for(i = 0; i < nodes_per_element; i++){
            const int row_node = connectivity[e][i];
            if(row_node < 0 || row_node >= num_nodes){
                return ST_FOM_ERR_DIMENSION;
            }

            for(j = 0; j < nodes_per_element; j++){
                const int col_node = connectivity[e][j];
                if(col_node < 0 || col_node >= num_nodes){
                    return ST_FOM_ERR_DIMENSION;
                }
                if(!st_fom_graph_contains(graph, row_node, col_node)){
                    return ST_FOM_ERR_DIMENSION;
                }
            }
        }
    }

    return ST_FOM_SUCCESS;
}
