#include "core_FOM_ST_spaceblock.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int* data;
    int size;
    int capacity;
} STSB_INT_VECTOR;

static void STSB_fail(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static void* STSB_calloc(size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        STSB_fail("allocation size overflow");
    }

    void* ptr = calloc(count, size);
    if(ptr == NULL) {
        STSB_fail("memory allocation failed");
    }
    return ptr;
}

static void* STSB_realloc(void* old_ptr, size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        STSB_fail("reallocation size overflow");
    }

    void* ptr = realloc(old_ptr, count * size);
    if(ptr == NULL) {
        STSB_fail("memory reallocation failed");
    }
    return ptr;
}

static int STSB_compare_int(const void* lhs, const void* rhs)
{
    const int a = *(const int*)lhs;
    const int b = *(const int*)rhs;
    return (a > b) - (a < b);
}

static void STSB_vector_push(STSB_INT_VECTOR* vector, int value)
{
    if(vector->size == vector->capacity) {
        const int next = vector->capacity == 0 ? 8 : 2 * vector->capacity;
        if(next < vector->capacity) {
            STSB_fail("graph-row capacity overflow");
        }
        vector->data = (int*)STSB_realloc(
            vector->data,
            (size_t)next,
            sizeof(int));
        vector->capacity = next;
    }
    vector->data[vector->size++] = value;
}

static void STSB_vector_sort_unique_without_self(
    STSB_INT_VECTOR* vector,
    int self)
{
    if(vector->size == 0) {
        return;
    }

    qsort(
        vector->data,
        (size_t)vector->size,
        sizeof(int),
        STSB_compare_int);

    int write = 0;
    for(int read=0; read<vector->size; read++) {
        const int value = vector->data[read];
        if(value == self) {
            continue;
        }
        if(write == 0 || value != vector->data[write - 1]) {
            vector->data[write++] = value;
        }
    }
    vector->size = write;
}

static int STSB_graph_row_contains(
    const STSB_SPATIAL_GRAPH* graph,
    int row,
    int col)
{
    int left = graph->index[row];
    int right = graph->index[row + 1];

    while(left < right) {
        const int mid = left + (right - left) / 2;
        const int value = graph->item[mid];
        if(value == col) {
            return 1;
        }
        if(value < col) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }
    return 0;
}

int STSB_get_num_dofs_per_node(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs)
{
    if(te == NULL || te->n_dof <= 0 || num_window_slabs <= 0) {
        STSB_fail("invalid temporal dimensions");
    }
    if(num_window_slabs > INT_MAX / te->n_dof) {
        STSB_fail("ST block DOF count exceeds INT_MAX");
    }
    return num_window_slabs * te->n_dof;
}

int STSB_get_block_dof(
    int slab,
    int alpha,
    const ST_TIME_ELEMENT* te)
{
    if(te == NULL || slab < 0 || alpha < 0 || alpha >= te->n_dof) {
        STSB_fail("invalid slab or temporal DOF");
    }
    if(slab > (INT_MAX - alpha) / te->n_dof) {
        STSB_fail("ST block DOF index exceeds INT_MAX");
    }
    return slab * te->n_dof + alpha;
}

size_t STSB_get_vector_index(
    int spatial_node,
    int slab,
    int alpha,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs)
{
    if(spatial_node < 0) {
        STSB_fail("negative spatial node");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);
    const int q = STSB_get_block_dof(slab, alpha, te);

    if((size_t)spatial_node > (SIZE_MAX - (size_t)q) / (size_t)ndof) {
        STSB_fail("ST vector index exceeds SIZE_MAX");
    }
    return (size_t)spatial_node * (size_t)ndof + (size_t)q;
}

void STSB_validate_configuration(
    const BBFE_DATA* fe,
    const VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs)
{
    if(fe == NULL || vals == NULL || te == NULL) {
        STSB_fail("NULL argument in STSB_validate_configuration");
    }
    if(fe->total_num_nodes <= 0 || fe->total_num_elems < 0 ||
       fe->local_num_nodes <= 0) {
        STSB_fail("invalid spatial mesh dimensions");
    }
    if(vals->dt <= 0.0) {
        STSB_fail("time spacing must be positive");
    }
    if(te->n_dof <= 0 || !te->is_nodal) {
        STSB_fail("a nodal temporal dG element is required");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);
    if((size_t)fe->total_num_nodes > SIZE_MAX / (size_t)ndof) {
        STSB_fail("window vector size exceeds SIZE_MAX");
    }
}

void STSB_spatial_graph_initialize(STSB_SPATIAL_GRAPH* graph)
{
    if(graph == NULL) {
        STSB_fail("NULL graph in STSB_spatial_graph_initialize");
    }
    memset(graph, 0, sizeof(*graph));
}


void STSB_spatial_graph_build_from_mesh(
    STSB_SPATIAL_GRAPH* graph,
    const BBFE_DATA* fe)
{
    if(graph == NULL || fe == NULL) {
        STSB_fail("NULL argument in STSB_spatial_graph_build_from_mesh");
    }
    if(graph->index != NULL || graph->item != NULL) {
        STSB_fail("graph is already allocated");
    }
    if(fe->total_num_nodes <= 0 || fe->total_num_elems < 0 ||
       fe->local_num_nodes <= 0 || fe->conn == NULL) {
        STSB_fail("invalid FE mesh in STSB_spatial_graph_build_from_mesh");
    }

    const int num_nodes = fe->total_num_nodes;
    STSB_INT_VECTOR* rows = (STSB_INT_VECTOR*)STSB_calloc(
        (size_t)num_nodes,
        sizeof(STSB_INT_VECTOR));

    for(int elem = 0; elem < fe->total_num_elems; elem++) {
        for(int a = 0; a < fe->local_num_nodes; a++) {
            const int row = fe->conn[elem][a];
            if(row < 0 || row >= num_nodes) {
                free(rows);
                STSB_fail("element row lies outside the local node range");
            }

            for(int b = 0; b < fe->local_num_nodes; b++) {
                const int col = fe->conn[elem][b];
                if(col < 0 || col >= num_nodes) {
                    free(rows);
                    STSB_fail("element column lies outside the local node range");
                }
                if(col != row) {
                    STSB_vector_push(&rows[row], col);
                }
            }
        }
    }

    int64_t total_items_64 = 0;
    for(int node = 0; node < num_nodes; node++) {
        STSB_vector_sort_unique_without_self(&rows[node], node);
        total_items_64 += rows[node].size;
    }
    if(total_items_64 > INT_MAX) {
        for(int node = 0; node < num_nodes; node++) {
            free(rows[node].data);
        }
        free(rows);
        STSB_fail("mesh-derived spatial graph exceeds INT_MAX");
    }

    graph->num_nodes = num_nodes;
    graph->num_items = (int)total_items_64;
    graph->index = (int*)STSB_calloc((size_t)num_nodes + 1, sizeof(int));
    graph->item = (int*)STSB_calloc((size_t)graph->num_items, sizeof(int));

    int position = 0;
    for(int node = 0; node < num_nodes; node++) {
        graph->index[node] = position;
        for(int p = 0; p < rows[node].size; p++) {
            graph->item[position++] = rows[node].data[p];
        }
        free(rows[node].data);
    }
    graph->index[num_nodes] = position;
    free(rows);
}

void STSB_spatial_graph_read(
    STSB_SPATIAL_GRAPH* graph,
    const char* label,
    const char* directory)
{
    if(graph == NULL || label == NULL || directory == NULL) {
        STSB_fail("NULL argument in STSB_spatial_graph_read");
    }
    if(graph->index != NULL || graph->item != NULL) {
        STSB_fail("graph is already allocated");
    }

    const char* fname = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        label);

    FILE* fp = NULL;
    fp = BBFE_sys_read_fopen(fp, fname, directory);

    int num_nodes = 0;
    if(fscanf(fp, "%d", &num_nodes) != 1 || num_nodes <= 0) {
        fclose(fp);
        STSB_fail("invalid spatial graph header");
    }

    STSB_INT_VECTOR* rows = (STSB_INT_VECTOR*)STSB_calloc(
        (size_t)num_nodes,
        sizeof(STSB_INT_VECTOR));
    unsigned char* seen = (unsigned char*)STSB_calloc(
        (size_t)num_nodes,
        sizeof(unsigned char));

    for(int input_row=0; input_row<num_nodes; input_row++) {
        int row_id = -1;
        int num_adj = -1;
        if(fscanf(fp, "%d %d", &row_id, &num_adj) != 2 ||
           row_id < 0 || row_id >= num_nodes || num_adj < 0 ||
           seen[row_id]) {
            fclose(fp);
            STSB_fail("invalid spatial graph row header");
        }
        seen[row_id] = 1;

        for(int j=0; j<num_adj; j++) {
            int col = -1;
            if(fscanf(fp, "%d", &col) != 1 || col < 0 || col >= num_nodes) {
                fclose(fp);
                STSB_fail("invalid spatial graph column");
            }
            STSB_vector_push(&rows[row_id], col);
        }
    }
    fclose(fp);

    int64_t total_items_64 = 0;
    for(int i=0; i<num_nodes; i++) {
        if(!seen[i]) {
            STSB_fail("spatial graph has a missing row");
        }
        STSB_vector_sort_unique_without_self(&rows[i], i);
        total_items_64 += rows[i].size;
    }
    if(total_items_64 > INT_MAX) {
        STSB_fail("local spatial graph exceeds INT_MAX");
    }

    graph->num_nodes = num_nodes;
    graph->num_items = (int)total_items_64;
    graph->index = (int*)STSB_calloc((size_t)num_nodes + 1, sizeof(int));
    graph->item = (int*)STSB_calloc((size_t)graph->num_items, sizeof(int));

    int position = 0;
    for(int i=0; i<num_nodes; i++) {
        graph->index[i] = position;
        for(int p=0; p<rows[i].size; p++) {
            graph->item[position++] = rows[i].data[p];
        }
        free(rows[i].data);
    }
    graph->index[num_nodes] = position;

    free(rows);
    free(seen);
}

void STSB_spatial_graph_finalize(STSB_SPATIAL_GRAPH* graph)
{
    if(graph == NULL) {
        return;
    }
    free(graph->index);
    free(graph->item);
    memset(graph, 0, sizeof(*graph));
}

void STSB_validate_spatial_graph_against_mesh(
    const STSB_SPATIAL_GRAPH* graph,
    const BBFE_DATA* fe)
{
    if(graph == NULL || fe == NULL) {
        STSB_fail("NULL argument in graph/mesh validation");
    }
    if(graph->num_nodes != fe->total_num_nodes) {
        fprintf(
            stderr,
            "ERROR: graph and mesh node counts differ.\n"
            "  graph nodes = %d\n"
            "  mesh nodes  = %d\n",
            graph->num_nodes,
            fe->total_num_nodes);
        exit(EXIT_FAILURE);
    }

    for(int e=0; e<fe->total_num_elems; e++) {
        for(int a=0; a<fe->local_num_nodes; a++) {
            const int row = fe->conn[e][a];
            if(row < 0 || row >= graph->num_nodes) {
                STSB_fail("element connectivity lies outside graph.dat");
            }
            for(int b=0; b<fe->local_num_nodes; b++) {
                const int col = fe->conn[e][b];
                if(row == col) {
                    continue;
                }
                if(!STSB_graph_row_contains(graph, row, col)) {
                    fprintf(
                        stderr,
                        "ERROR: graph.dat lacks an element coupling.\n"
                        "  element = %d\n"
                        "  row     = %d\n"
                        "  column  = %d\n",
                        e,
                        row,
                        col);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}

void STSB_set_nonzero_pattern_from_spatial_graph(
    MONOLIS* monolis,
    const STSB_SPATIAL_GRAPH* graph,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_internal_space_nodes)
{
    if(monolis == NULL || graph == NULL || te == NULL ||
       graph->index == NULL || graph->item == NULL) {
        STSB_fail("NULL argument in nonzero-pattern setup");
    }
    if(num_internal_space_nodes <= 0 ||
       num_internal_space_nodes > graph->num_nodes) {
        STSB_fail("invalid number of internal spatial nodes");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        monolis,
        graph->num_nodes,
        ndof,
        graph->index,
        graph->item);

    monolis->mat.N = num_internal_space_nodes;
    monolis->mat.NP = graph->num_nodes;
}

void STSB_compute_element_mass_diffusion(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    int elem,
    double* Me,
    double* Ke)
{
    if(fe == NULL || basis == NULL || Me == NULL || Ke == NULL ||
       elem < 0 || elem >= fe->total_num_elems) {
        STSB_fail("invalid argument in element matrix integration");
    }

    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* val_ip = BB_std_calloc_1d_double(val_ip, np);
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double** local_x = BB_std_calloc_2d_double(local_x, nl, 3);
    double** x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
    double* a_ip = BB_std_calloc_1d_double(a_ip, np);
    double* k_ip = BB_std_calloc_1d_double(k_ip, np);

    BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, elem, fe);
    BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, elem, 3);

    for(int p=0; p<np; p++) {
        BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
        a_ip[p] = manusol_get_mass_coef(x_ip[p]);
        k_ip[p] = manusol_get_diff_coef(x_ip[p]);
    }

    for(int i=0; i<nl; i++) {
        for(int j=0; j<nl; j++) {
            for(int p=0; p<np; p++) {
                val_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                    basis->N[p][i],
                    basis->N[p][j],
                    a_ip[p]);
            }
            Me[i * nl + j] = BBFE_std_integ_calc(
                np,
                val_ip,
                basis->integ_weight,
                Jacobian_ip);

            for(int p=0; p<np; p++) {
                val_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
                    fe->geo[elem][p].grad_N[i],
                    fe->geo[elem][p].grad_N[j],
                    k_ip[p]);
            }
            Ke[i * nl + j] = BBFE_std_integ_calc(
                np,
                val_ip,
                basis->integ_weight,
                Jacobian_ip);
        }
    }

    BB_std_free_1d_double(val_ip, np);
    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_x, nl, 3);
    BB_std_free_2d_double(x_ip, np, 3);
    BB_std_free_1d_double(a_ip, np);
    BB_std_free_1d_double(k_ip, np);
}

void STSB_compute_element_source(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    const ST_TIME_ELEMENT* te,
    int elem,
    double slab_start_time,
    double dt,
    double* Fe)
{
    if(fe == NULL || basis == NULL || te == NULL || Fe == NULL ||
       elem < 0 || elem >= fe->total_num_elems || dt <= 0.0) {
        STSB_fail("invalid argument in element source integration");
    }

    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int nt = te->n_dof;

    double* val_ip = BB_std_calloc_1d_double(val_ip, np);
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double** local_x = BB_std_calloc_2d_double(local_x, nl, 3);
    double** x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
    double** v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
    double* a_ip = BB_std_calloc_1d_double(a_ip, np);
    double* k_ip = BB_std_calloc_1d_double(k_ip, np);

    BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, elem, fe);
    BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, elem, 3);

    for(int p=0; p<np; p++) {
        BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
        manusol_get_conv_vel(v_ip[p], x_ip[p]);
        a_ip[p] = manusol_get_mass_coef(x_ip[p]);
        k_ip[p] = manusol_get_diff_coef(x_ip[p]);
    }

    for(int alpha=0; alpha<nt; alpha++) {
        for(int i=0; i<nl; i++) {
            for(int p=0; p<np; p++) {
                val_ip[p] = 0.0;
                for(int q=0; q<te->num_time_ip; q++) {
                    const double coef = te->weight[q] * te->psi[alpha][q];
                    if(fabs(coef) < 1.0e-30) {
                        continue;
                    }
                    const double time = slab_start_time + te->tau[q] * dt;
                    const double source = manusol_get_source(
                        x_ip[p],
                        time,
                        a_ip[p],
                        v_ip[p],
                        k_ip[p]);
                    val_ip[p] += coef * BBFE_elemmat_convdiff_vec_source(
                        basis->N[p][i],
                        source);
                }
            }
            Fe[alpha * nl + i] = BBFE_std_integ_calc(
                np,
                val_ip,
                basis->integ_weight,
                Jacobian_ip);
        }
    }

    BB_std_free_1d_double(val_ip, np);
    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_x, nl, 3);
    BB_std_free_2d_double(x_ip, np, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_1d_double(a_ip, np);
    BB_std_free_1d_double(k_ip, np);
}

void STSB_assemble_window_matrix(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs)
{
    if(monolis == NULL || fe == NULL || basis == NULL ||
       vals == NULL || te == NULL) {
        STSB_fail("NULL argument in STSB_assemble_window_matrix");
    }

    STSB_validate_configuration(fe, vals, te, num_window_slabs);

    const int nl = fe->local_num_nodes;
    const int nt = te->n_dof;
    double* Me = (double*)STSB_calloc((size_t)nl * nl, sizeof(double));
    double* Ke = (double*)STSB_calloc((size_t)nl * nl, sizeof(double));

    for(int elem=0; elem<fe->total_num_elems; elem++) {
        STSB_compute_element_mass_diffusion(fe, basis, elem, Me, Ke);

        for(int slab=0; slab<num_window_slabs; slab++) {
            for(int alpha=0; alpha<nt; alpha++) {
                const int row_dof = STSB_get_block_dof(slab, alpha, te);

                for(int beta=0; beta<nt; beta++) {
                    const int col_dof = STSB_get_block_dof(slab, beta, te);
                    const double coef_mass =
                        (te->Dt[alpha][beta] + te->E_self[alpha][beta]) /
                        vals->dt;
                    const double coef_diff = te->Mt[alpha][beta];

                    for(int i=0; i<nl; i++) {
                        const int row_node = fe->conn[elem][i];
                        for(int j=0; j<nl; j++) {
                            const int col_node = fe->conn[elem][j];
                            const double value =
                                coef_mass * Me[i * nl + j] +
                                coef_diff * Ke[i * nl + j];

                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                row_node,
                                col_node,
                                row_dof,
                                col_dof,
                                value);
                        }
                    }
                }

                if(slab > 0) {
                    for(int beta=0; beta<nt; beta++) {
                        const double coef =
                            -te->E_prev[alpha][beta] / vals->dt;
                        if(fabs(coef) < 1.0e-30) {
                            continue;
                        }

                        const int col_dof =
                            STSB_get_block_dof(slab - 1, beta, te);

                        for(int i=0; i<nl; i++) {
                            const int row_node = fe->conn[elem][i];
                            for(int j=0; j<nl; j++) {
                                const int col_node = fe->conn[elem][j];
                                monolis_add_scalar_to_sparse_matrix_R(
                                    monolis,
                                    row_node,
                                    col_node,
                                    row_dof,
                                    col_dof,
                                    coef * Me[i * nl + j]);
                            }
                        }
                    }
                }
            }
        }
    }

    free(Me);
    free(Ke);
}

void STSB_add_previous_trace_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    const double* previous_trace)
{
    if(monolis == NULL || fe == NULL || basis == NULL || vals == NULL ||
       te == NULL || previous_trace == NULL) {
        STSB_fail("NULL argument in STSB_add_previous_trace_rhs");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);
    const int nl = fe->local_num_nodes;
    const int nt = te->n_dof;
    double* Me = (double*)STSB_calloc((size_t)nl * nl, sizeof(double));
    double* Ke = (double*)STSB_calloc((size_t)nl * nl, sizeof(double));

    for(int elem=0; elem<fe->total_num_elems; elem++) {
        STSB_compute_element_mass_diffusion(fe, basis, elem, Me, Ke);

        for(int alpha=0; alpha<nt; alpha++) {
            if(fabs(te->e_minus[alpha]) < 1.0e-30) {
                continue;
            }

            const int q = STSB_get_block_dof(0, alpha, te);
            for(int i=0; i<nl; i++) {
                double inflow = 0.0;
                for(int j=0; j<nl; j++) {
                    const int node_j = fe->conn[elem][j];
                    inflow += Me[i * nl + j] * previous_trace[node_j];
                }

                const int node_i = fe->conn[elem][i];
                monolis->mat.R.B[
                    (size_t)node_i * (size_t)ndof + (size_t)q
                ] += te->e_minus[alpha] * inflow / vals->dt;
            }
        }
    }

    free(Me);
    free(Ke);
}

void STSB_add_source_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time)
{
    if(monolis == NULL || fe == NULL || basis == NULL ||
       vals == NULL || te == NULL) {
        STSB_fail("NULL argument in STSB_add_source_rhs");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);
    const int nl = fe->local_num_nodes;
    const int nt = te->n_dof;
    double* Fe = (double*)STSB_calloc((size_t)nt * nl, sizeof(double));

    for(int slab=0; slab<num_window_slabs; slab++) {
        const double slab_start_time =
            window_start_time + (double)slab * vals->dt;

        for(int elem=0; elem<fe->total_num_elems; elem++) {
            memset(Fe, 0, (size_t)nt * (size_t)nl * sizeof(double));
            STSB_compute_element_source(
                fe,
                basis,
                te,
                elem,
                slab_start_time,
                vals->dt,
                Fe);

            for(int alpha=0; alpha<nt; alpha++) {
                const int q = STSB_get_block_dof(slab, alpha, te);
                for(int i=0; i<nl; i++) {
                    const int node_i = fe->conn[elem][i];
                    monolis->mat.R.B[
                        (size_t)node_i * (size_t)ndof + (size_t)q
                    ] += Fe[alpha * nl + i];
                }
            }
        }
    }

    free(Fe);
}

void STSB_assemble_window_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    const double* previous_trace)
{
    STSB_add_previous_trace_rhs(
        monolis,
        fe,
        basis,
        vals,
        te,
        num_window_slabs,
        previous_trace);

    STSB_add_source_rhs(
        monolis,
        fe,
        basis,
        vals,
        te,
        num_window_slabs,
        window_start_time);
}

void STSB_apply_dirichlet_bc(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    double* rhs)
{
    if(monolis == NULL || fe == NULL || bc == NULL || vals == NULL ||
       te == NULL || rhs == NULL) {
        STSB_fail("NULL argument in STSB_apply_dirichlet_bc");
    }

    for(int slab=0; slab<num_window_slabs; slab++) {
        const double slab_start_time =
            window_start_time + (double)slab * vals->dt;

        for(int alpha=0; alpha<te->n_dof; alpha++) {
            const double time =
                slab_start_time + te->dof_tau[alpha] * vals->dt;

            manusol_set_theo_sol(fe, vals->theo_sol, time);
            BBFE_manusol_set_bc_scalar(fe, bc, vals->theo_sol, time);

            const int q = STSB_get_block_dof(slab, alpha, te);
            for(int node=0; node<fe->total_num_nodes; node++) {
                if(!bc->D_bc_exists[node]) {
                    continue;
                }
                monolis_set_Dirichlet_bc_R(
                    monolis,
                    rhs,
                    node,
                    q,
                    bc->imposed_D_val[node]);
            }
        }
    }
}

void STSB_restore_dirichlet_values(
    BBFE_DATA* fe,
    BBFE_BC* bc,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    double* window_solution)
{
    if(fe == NULL || bc == NULL || vals == NULL || te == NULL ||
       window_solution == NULL) {
        STSB_fail("NULL argument in STSB_restore_dirichlet_values");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);

    for(int slab=0; slab<num_window_slabs; slab++) {
        const double slab_start_time =
            window_start_time + (double)slab * vals->dt;

        for(int alpha=0; alpha<te->n_dof; alpha++) {
            const double time =
                slab_start_time + te->dof_tau[alpha] * vals->dt;

            manusol_set_theo_sol(fe, vals->theo_sol, time);
            BBFE_manusol_set_bc_scalar(fe, bc, vals->theo_sol, time);

            const int q = STSB_get_block_dof(slab, alpha, te);
            for(int node=0; node<fe->total_num_nodes; node++) {
                if(!bc->D_bc_exists[node]) {
                    continue;
                }
                window_solution[
                    (size_t)node * (size_t)ndof + (size_t)q
                ] = bc->imposed_D_val[node];
            }
        }
    }
}

void STSB_set_previous_trace_from_nodal_value(
    const BBFE_DATA* fe,
    const double* nodal_value,
    double* previous_trace)
{
    if(fe == NULL || nodal_value == NULL || previous_trace == NULL) {
        STSB_fail("NULL argument in previous-trace initialization");
    }

    memcpy(
        previous_trace,
        nodal_value,
        (size_t)fe->total_num_nodes * sizeof(double));
}

void STSB_extract_slab_right_trace(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int slab,
    int num_space_nodes,
    const double* window_solution,
    double* nodal_trace)
{
    if(te == NULL || window_solution == NULL || nodal_trace == NULL ||
       num_space_nodes <= 0 || slab < 0 || slab >= num_window_slabs) {
        STSB_fail("invalid argument in right-trace extraction");
    }

    const int ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);

    for(int node=0; node<num_space_nodes; node++) {
        double value = 0.0;
        for(int alpha=0; alpha<te->n_dof; alpha++) {
            const int q = STSB_get_block_dof(slab, alpha, te);
            value += te->e_plus[alpha] * window_solution[
                (size_t)node * (size_t)ndof + (size_t)q
            ];
        }
        nodal_trace[node] = value;
    }
}

void STSB_extract_last_right_trace(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_space_nodes,
    const double* window_solution,
    double* previous_trace)
{
    STSB_extract_slab_right_trace(
        te,
        num_window_slabs,
        num_window_slabs - 1,
        num_space_nodes,
        window_solution,
        previous_trace);
}

static void STSB_prepare_window_system(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    const double* previous_trace)
{
    if(sys == NULL || te == NULL || previous_trace == NULL) {
        STSB_fail("NULL argument in STSB_prepare_window_system");
    }

    /*
     * Recreate the exact algebraic system used by the reference FOM.
     * This helper is called both before and after the iterative solve.
     * The second call is important because the solver wrapper is allowed to
     * alter matrix/RHS work arrays while constructing a preconditioner or
     * solving the linear system.
     */
    monolis_clear_mat_value_rhs_R(&sys->monolis);
    monolis_copy_mat_value_R(&sys->monolis0, &sys->monolis);

    STSB_assemble_window_rhs(
        &sys->monolis,
        &sys->fe,
        &sys->basis,
        &sys->vals,
        te,
        num_window_slabs,
        window_start_time,
        previous_trace);

    STSB_apply_dirichlet_bc(
        &sys->monolis,
        &sys->fe,
        &sys->bc,
        &sys->vals,
        te,
        num_window_slabs,
        window_start_time,
        sys->monolis.mat.R.B);
}

void STSB_solve_window(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    int window_id,
    const double* previous_trace,
    double* window_solution)
{
    int ndof;
    const double window_end_time =
        window_start_time + (double)num_window_slabs *
            ((sys != NULL) ? sys->vals.dt : 0.0);

    if(sys == NULL || te == NULL || previous_trace == NULL ||
       window_solution == NULL) {
        STSB_fail("NULL argument in STSB_solve_window");
    }

    ndof = STSB_get_num_dofs_per_node(te, num_window_slabs);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "\n%s ----------------- ST spatial-block window %d : "
            "[%e, %e] ----------------\n",
            CODENAME,
            window_id,
            window_start_time,
            window_end_time);
    }

    STSB_prepare_window_system(
        sys,
        te,
        num_window_slabs,
        window_start_time,
        previous_trace);

    /*
     * Keep the reference solve path unchanged.  The solution is subsequently
     * synchronized and the exact pre-solve matrix/RHS pair is reconstructed
     * for residual checks and ROM projection diagnostics.
     */
    BBFE_sys_monowrap_solve(
        &sys->monolis,
        &sys->monolis_com,
        window_solution,
        MONOLIS_ITER_BICGSTAB,
        MONOLIS_PREC_DIAG,
        sys->vals.mat_max_iter,
        sys->vals.mat_epsilon);

    /* Enforce the prescribed values on owned and local halo boundary nodes. */
    STSB_restore_dirichlet_values(
        &sys->fe,
        &sys->bc,
        &sys->vals,
        te,
        num_window_slabs,
        window_start_time,
        window_solution);

    /*
     * The next-window trace and all post-solve diagnostics require current
     * halo values.  Do not rely on the solver wrapper to leave them updated.
     */
    monolis_mpi_update_R(
        &sys->monolis_com,
        sys->fe.total_num_nodes,
        ndof,
        window_solution);

    /*
     * Restore the exact strong-Dirichlet matrix/RHS after the solver.  The
     * previous implementation evaluated the residual using solver work
     * arrays left in sys->monolis, which produced a false O(1e-1) free-row
     * residual even immediately after the FOM solve.
     */
    STSB_prepare_window_system(
        sys,
        te,
        num_window_slabs,
        window_start_time,
        previous_trace);
}


void STSB_initialize_reference_fom(
    FE_SYSTEM* sys,
    ST_TIME_ELEMENT* te,
    STSB_SPATIAL_GRAPH* graph,
    int argc,
    char* argv[],
    int degree,
    int num_window_slabs,
    int* num_windows,
    int* num_internal_space_nodes)
{
    if(sys == NULL || te == NULL || graph == NULL ||
       num_windows == NULL || num_internal_space_nodes == NULL) {
        STSB_fail("NULL argument in STSB_initialize_reference_fom");
    }

    memset(sys, 0, sizeof(*sys));
    memset(te, 0, sizeof(*te));
    STSB_spatial_graph_initialize(graph);

    sys->cond.directory = BBFE_convdiff_get_directory_name(
        argc,
        argv,
        CODENAME);

    read_calc_conditions(&sys->vals, sys->cond.directory);

    BBFE_convdiff_pre(
        &sys->fe,
        &sys->basis,
        &sys->bc,
        &sys->monolis,
        &sys->monolis_com,
        argc,
        argv,
        sys->cond.directory,
        sys->vals.num_ip_each_axis,
        true);

    memory_allocation_nodal_values(
        &sys->vals,
        sys->fe.total_num_nodes);

    manusol_set_init_value(&sys->fe, sys->vals.T);

    BBFE_elemmat_set_Jacobi_mat(&sys->fe, &sys->basis);
    BBFE_elemmat_set_shapefunc_derivative(&sys->fe, &sys->basis);

    ST_TimeElement_init_dG(te, degree);
    STSB_validate_configuration(
        &sys->fe,
        &sys->vals,
        te,
        num_window_slabs);

    *num_windows = STSB_get_num_windows(
        sys->vals.finish_time,
        sys->vals.dt,
        num_window_slabs);

    *num_internal_space_nodes = sys->monolis.mat.N;
    if(*num_internal_space_nodes <= 0 ||
       *num_internal_space_nodes > sys->fe.total_num_nodes) {
        STSB_fail("invalid internal-node count from spatial MONOLIS data");
    }

    /* The verified reference FOM reconstructs adjacency from FE connectivity. */
    STSB_spatial_graph_build_from_mesh(graph, &sys->fe);
    STSB_validate_spatial_graph_against_mesh(graph, &sys->fe);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "STSB spatial graph source: FE mesh connectivity\n");
    }

    monolis_initialize(&sys->monolis0);

    STSB_set_nonzero_pattern_from_spatial_graph(
        &sys->monolis0,
        graph,
        te,
        num_window_slabs,
        *num_internal_space_nodes);

    monolis_clear_mat_value_rhs_R(&sys->monolis0);

    STSB_assemble_window_matrix(
        &sys->monolis0,
        &sys->fe,
        &sys->basis,
        &sys->vals,
        te,
        num_window_slabs);

    monolis_copy_mat_R(&sys->monolis0, &sys->monolis);
}

void STSB_finalize_reference_fom(
    FE_SYSTEM* sys,
    STSB_SPATIAL_GRAPH* graph)
{
    if(sys == NULL || graph == NULL) {
        return;
    }

    STSB_spatial_graph_finalize(graph);
    BBFE_convdiff_finalize(&sys->fe, &sys->basis, &sys->bc);
    monolis_finalize(&sys->monolis);
    monolis_finalize(&sys->monolis0);
}

int STSB_evaluate_constrained_residual(
    FE_SYSTEM* sys,
    int st_dof_per_space_node,
    const double* window_solution,
    STSB_CONSTRAINED_RESIDUAL_SUM* global_sum)
{
    const int total_space_nodes =
        (sys != NULL) ? sys->fe.total_num_nodes : 0;
    const int internal_space_nodes =
        (sys != NULL) ? sys->monolis.mat.N : 0;
    int local_st_dof;
    int internal_st_dof;
    double* solution_work = NULL;
    double* product = NULL;
    double reduced_values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if(sys == NULL || window_solution == NULL || global_sum == NULL ||
       total_space_nodes <= 0 || internal_space_nodes <= 0 ||
       internal_space_nodes > total_space_nodes ||
       st_dof_per_space_node <= 0 || sys->monolis.mat.R.B == NULL) {
        return -1;
    }

    if(total_space_nodes > INT_MAX / st_dof_per_space_node ||
       internal_space_nodes > INT_MAX / st_dof_per_space_node) {
        return -1;
    }

    local_st_dof = total_space_nodes * st_dof_per_space_node;
    internal_st_dof = internal_space_nodes * st_dof_per_space_node;

    solution_work = (double*)STSB_calloc(
        (size_t)local_st_dof,
        sizeof(double));
    product = (double*)STSB_calloc(
        (size_t)local_st_dof,
        sizeof(double));

    memcpy(
        solution_work,
        window_solution,
        (size_t)local_st_dof * sizeof(double));

    monolis_mpi_update_R(
        &sys->monolis_com,
        total_space_nodes,
        st_dof_per_space_node,
        solution_work);

    monolis_matvec_product_R(
        &sys->monolis,
        &sys->monolis_com,
        solution_work,
        product);

    for(int row = 0; row < internal_st_dof; row++) {
        const int node = row / st_dof_per_space_node;
        const double rhs = sys->monolis.mat.R.B[row];
        const double residual = product[row] - rhs;
        const int is_dirichlet = sys->bc.D_bc_exists[node] ? 1 : 0;

        reduced_values[0] += residual * residual;
        reduced_values[1] += rhs * rhs;

        if(is_dirichlet) {
            reduced_values[4] += residual * residual;
            reduced_values[5] += rhs * rhs;
        }
        else {
            reduced_values[2] += residual * residual;
            reduced_values[3] += rhs * rhs;
        }
    }

    monolis_allreduce_R(
        6,
        reduced_values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    global_sum->residual_all2 = reduced_values[0];
    global_sum->rhs_all2 = reduced_values[1];
    global_sum->residual_free2 = reduced_values[2];
    global_sum->rhs_free2 = reduced_values[3];
    global_sum->residual_dirichlet2 = reduced_values[4];
    global_sum->rhs_dirichlet2 = reduced_values[5];

    free(solution_work);
    free(product);
    return 0;
}

int STSB_compute_window_checksum(
    const FE_SYSTEM* sys,
    int st_dof_per_space_node,
    const double* window_solution,
    const double* right_trace,
    double* solution_l2,
    double* trace_l2)
{
    if(sys == NULL || window_solution == NULL || right_trace == NULL ||
       solution_l2 == NULL || trace_l2 == NULL ||
       st_dof_per_space_node <= 0 || sys->monolis.mat.N <= 0 ||
       sys->monolis.mat.N > sys->fe.total_num_nodes) {
        return -1;
    }

    const int internal_nodes = sys->monolis.mat.N;
    double values[2] = {0.0, 0.0};

    for(int node = 0; node < internal_nodes; node++) {
        for(int q = 0; q < st_dof_per_space_node; q++) {
            const double value = window_solution[
                (size_t)node * (size_t)st_dof_per_space_node + (size_t)q];
            values[0] += value * value;
        }
        values[1] += right_trace[node] * right_trace[node];
    }

    monolis_allreduce_R(
        2,
        values,
        MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    *solution_l2 = sqrt(fmax(values[0], 0.0));
    *trace_l2 = sqrt(fmax(values[1], 0.0));
    return 0;
}

int STSB_get_num_windows(
    double finish_time,
    double dt,
    int num_window_slabs)
{
    if(finish_time <= 0.0 || dt <= 0.0 || num_window_slabs <= 0) {
        STSB_fail("invalid time conditions in STSB_get_num_windows");
    }

    const int total_num_steps = (int)llround(finish_time / dt);
    const double reconstructed_time = (double)total_num_steps * dt;
    const double tolerance = 1.0e-12 * fmax(1.0, fabs(finish_time));

    if(fabs(reconstructed_time - finish_time) > tolerance) {
        fprintf(
            stderr,
            "ERROR: finish_time / dt must be an integer.\n"
            "  finish_time = %.16e\n"
            "  dt          = %.16e\n",
            finish_time,
            dt);
        exit(EXIT_FAILURE);
    }

    if(total_num_steps % num_window_slabs != 0) {
        fprintf(
            stderr,
            "ERROR: total time steps must be divisible by window slabs.\n"
            "  total steps  = %d\n"
            "  window slabs = %d\n",
            total_num_steps,
            num_window_slabs);
        exit(EXIT_FAILURE);
    }

    return total_num_steps / num_window_slabs;
}
