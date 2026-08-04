#include "rom_std_strom.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ROM_STROM_USE_CBLAS
#include <cblas.h>
#endif

#ifdef ROM_STROM_USE_LAPACK
/* Fortran LAPACK ABI; link with -llapack -lblas (or OpenBLAS). */
extern void dgesv_(
    const int* n,
    const int* nrhs,
    double* a,
    const int* lda,
    int* ipiv,
    double* b,
    const int* ldb,
    int* info);
#endif

/*
 * Replace these includes with the actual Monolis headers used by the project.
 * The implementation uses:
 *   monolis_mpi_get_self_comm
 *   monolis_scalapack_comm_initialize
 *   monolis_scalapack_gesvd_R
 *   monolis_scalapack_comm_finalize
 * and ROM_BB_estimate_num_pod_modes.
 */
#include "monolis.h"
#include "monolis_utils.h"

#define ROM_STROM_BINARY_MAGIC UINT32_C(0x5354524D) /* "STRM" */
#define ROM_STROM_BINARY_VERSION UINT32_C(2)

static double* strom_alloc_vector(int n)
{
    if(n <= 0){
        return NULL;
    }
    return (double*)calloc((size_t)n, sizeof(double));
}

static int* strom_alloc_int_vector(int n)
{
    if(n <= 0){
        return NULL;
    }
    return (int*)calloc((size_t)n, sizeof(int));
}

static double** strom_alloc_matrix(int rows, int columns)
{
    double** matrix;
    double* data;

    if(rows <= 0 || columns <= 0){
        return NULL;
    }

    matrix = (double**)malloc((size_t)rows * sizeof(double*));
    if(matrix == NULL){
        return NULL;
    }

    data = (double*)calloc(
        (size_t)rows * (size_t)columns,
        sizeof(double));

    if(data == NULL){
        free(matrix);
        return NULL;
    }

    for(int i = 0; i < rows; i++){
        matrix[i] = data + (size_t)i * (size_t)columns;
    }

    return matrix;
}

static void strom_free_matrix(double*** matrix)
{
    if(matrix == NULL || *matrix == NULL){
        return;
    }

    free((*matrix)[0]);
    free(*matrix);
    *matrix = NULL;
}

static void strom_zero_vector(double* vector, int n)
{
    if(vector != NULL && n > 0){
        memset(vector, 0, (size_t)n * sizeof(double));
    }
}

static int strom_min_int(int a, int b)
{
    return a < b ? a : b;
}

static int strom_validate_window(const ROM_STROM_WINDOW* window)
{
    if(window == NULL ||
       window->hlpod_vals == NULL ||
       window->hlpod_mat == NULL ||
       window->spatial_dof <= 0 ||
       window->st_dof <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_index(
    int node,
    int component,
    int slab,
    int time_basis,
    int total_num_nodes,
    int dof,
    int num_time_basis)
{
    const int spatial_dof = total_num_nodes * dof;
    const int spatial_id = node * dof + component;
    const int time_id = slab * num_time_basis + time_basis;

    if(node < 0 || node >= total_num_nodes ||
       component < 0 || component >= dof ||
       slab < 0 || time_basis < 0 ||
       time_basis >= num_time_basis){
        return -1;
    }

    return time_id * spatial_dof + spatial_id;
}

int ROM_std_strom_window_initialize(
    ROM_STROM_WINDOW* window,
    int window_id,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT* hlpod_mat,
    int total_num_nodes,
    int dof,
    int num_slabs,
    int num_time_basis,
    int num_snapshots,
    int num_modes_max,
    int allocate_reference_state)
{
    long long spatial_dof;
    long long st_dof;

    if(window == NULL || hlpod_vals == NULL || hlpod_mat == NULL ||
       total_num_nodes <= 0 || dof <= 0 || num_slabs <= 0 ||
       num_time_basis <= 0 || num_snapshots <= 0 || num_modes_max <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    spatial_dof = (long long)total_num_nodes * (long long)dof;
    st_dof = spatial_dof * (long long)num_slabs * (long long)num_time_basis;

    if(spatial_dof > INT32_MAX || st_dof > INT32_MAX){
        return ROM_STROM_ERR_DIMENSION;
    }

    memset(window, 0, sizeof(*window));

    window->window_id = window_id;
    window->total_num_nodes = total_num_nodes;
    window->dof = dof;
    window->num_slabs = num_slabs;
    window->num_time_basis = num_time_basis;
    window->spatial_dof = (int)spatial_dof;
    window->st_dof = (int)st_dof;
    window->num_snapshots = num_snapshots;
    window->num_modes_max = num_modes_max;
    window->hlpod_vals = hlpod_vals;
    window->hlpod_mat = hlpod_mat;
    window->owns_pod_modes = 1;
    window->owns_singular_values = 1;

    memset(window->hlpod_vals, 0, sizeof(*window->hlpod_vals));
    memset(window->hlpod_mat, 0, sizeof(*window->hlpod_mat));

    window->hlpod_mat->snapmat =
        strom_alloc_matrix(window->st_dof, num_snapshots);

    if(window->hlpod_mat->snapmat == NULL){
        ROM_std_strom_window_finalize(window);
        return ROM_STROM_ERR_ALLOCATION;
    }

    window->reduced_rhs = strom_alloc_vector(num_modes_max);
    window->reduced_coef = strom_alloc_vector(num_modes_max);

    if(window->reduced_rhs == NULL || window->reduced_coef == NULL){
        ROM_std_strom_window_finalize(window);
        return ROM_STROM_ERR_ALLOCATION;
    }

    if(allocate_reference_state){
        window->reference_state = strom_alloc_vector(window->st_dof);
        if(window->reference_state == NULL){
            ROM_std_strom_window_finalize(window);
            return ROM_STROM_ERR_ALLOCATION;
        }
    }

    return ROM_STROM_SUCCESS;
}

void ROM_std_strom_window_finalize(ROM_STROM_WINDOW* window)
{
    if(window == NULL){
        return;
    }

    if(window->hlpod_mat != NULL){
        strom_free_matrix(&window->hlpod_mat->snapmat);

        if(window->owns_pod_modes){
            strom_free_matrix(&window->hlpod_mat->pod_modes);
        }
        else{
            window->hlpod_mat->pod_modes = NULL;
        }

        free(window->hlpod_mat->node_id);
        window->hlpod_mat->node_id = NULL;

        free(window->hlpod_mat->num_modes_internal);
        window->hlpod_mat->num_modes_internal = NULL;

        free(window->hlpod_mat->n_internal_vertex_subd);
        window->hlpod_mat->n_internal_vertex_subd = NULL;
    }

    if(window->owns_singular_values){
        free(window->singular_values);
    }
    window->singular_values = NULL;

    strom_free_matrix(&window->reduced_diag);
    strom_free_matrix(&window->reduced_coupling_prev);

    free(window->reduced_rhs);
    free(window->reduced_coef);
    free(window->reference_state);

    window->reduced_rhs = NULL;
    window->reduced_coef = NULL;
    window->reference_state = NULL;

    window->num_modes = 0;
    window->num_singular_values = 0;
    window->owns_pod_modes = 1;
    window->owns_singular_values = 1;
}

int ROM_std_strom_system_initialize(
    ROM_STROM_SYSTEM* system,
    int num_windows)
{
    if(system == NULL || num_windows <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    memset(system, 0, sizeof(*system));
    system->num_windows = num_windows;

    system->windows = (ROM_STROM_WINDOW*)calloc(
        (size_t)num_windows,
        sizeof(ROM_STROM_WINDOW));

    system->mode_offset = strom_alloc_int_vector(num_windows);

    if(system->windows == NULL || system->mode_offset == NULL){
        ROM_std_strom_system_finalize(system);
        return ROM_STROM_ERR_ALLOCATION;
    }

    return ROM_STROM_SUCCESS;
}

void ROM_std_strom_system_finalize(ROM_STROM_SYSTEM* system)
{
    if(system == NULL){
        return;
    }

    if(system->windows != NULL){
        for(int w = 0; w < system->num_windows; w++){
            ROM_std_strom_window_finalize(&system->windows[w]);
        }
    }

    free(system->windows);
    free(system->mode_offset);

    system->windows = NULL;
    system->mode_offset = NULL;
    system->num_windows = 0;
    system->total_reduced_dof = 0;
}

int ROM_std_strom_update_mode_offsets(ROM_STROM_SYSTEM* system)
{
    int offset = 0;

    if(system == NULL || system->windows == NULL ||
       system->mode_offset == NULL || system->num_windows <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int w = 0; w < system->num_windows; w++){
        if(system->windows[w].num_modes <= 0){
            return ROM_STROM_ERR_DIMENSION;
        }

        system->mode_offset[w] = offset;
        offset += system->windows[w].num_modes;
    }

    system->total_reduced_dof = offset;
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_store_snapshot(
    ROM_STROM_WINDOW* window,
    const double* window_solution,
    int snapshot_id,
    int subtract_reference_state)
{
    int status = strom_validate_window(window);
    if(status != ROM_STROM_SUCCESS || window_solution == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(window->hlpod_mat->snapmat == NULL ||
       snapshot_id < 0 || snapshot_id >= window->num_snapshots){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int i = 0; i < window->st_dof; i++){
        double value = window_solution[i];

        if(subtract_reference_state && window->reference_state != NULL){
            value -= window->reference_state[i];
        }

        window->hlpod_mat->snapmat[i][snapshot_id] = value;
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_set_pod_modes(
    ROM_STROM_WINDOW* window,
    double rom_epsilon)
{
    double** left_vectors = NULL;
    double** right_vectors_t = NULL;
    double* singular_values = NULL;
    int scalapack_comm = 0;
    int comm;
    int max_available_modes;
    int selected_modes;

    int status = strom_validate_window(window);
    if(status != ROM_STROM_SUCCESS ||
       window->hlpod_mat->snapmat == NULL ||
       rom_epsilon < 0.0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    max_available_modes = strom_min_int(
        window->st_dof,
        window->num_snapshots);

    left_vectors = strom_alloc_matrix(
        window->st_dof,
        window->num_snapshots);

    right_vectors_t = strom_alloc_matrix(
        window->num_snapshots,
        window->num_snapshots);

    singular_values = strom_alloc_vector(window->num_snapshots);

    if(left_vectors == NULL || right_vectors_t == NULL ||
       singular_values == NULL){
        strom_free_matrix(&left_vectors);
        strom_free_matrix(&right_vectors_t);
        free(singular_values);
        return ROM_STROM_ERR_ALLOCATION;
    }

    comm = monolis_mpi_get_self_comm();
    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    monolis_scalapack_gesvd_R(
        window->st_dof,
        window->num_snapshots,
        window->hlpod_mat->snapmat,
        left_vectors,
        singular_values,
        right_vectors_t,
        comm,
        scalapack_comm);

    monolis_scalapack_comm_finalize(scalapack_comm);

    if(rom_epsilon >= 1.0){
        selected_modes = window->num_modes_max;
    }
    else{
        selected_modes = ROM_BB_estimate_num_pod_modes(
            singular_values,
            max_available_modes,
            rom_epsilon);
    }

    if(selected_modes < 1){
        selected_modes = 1;
    }
    if(selected_modes > window->num_modes_max){
        selected_modes = window->num_modes_max;
    }
    if(selected_modes > max_available_modes){
        selected_modes = max_available_modes;
    }

    if(window->owns_pod_modes){
        strom_free_matrix(&window->hlpod_mat->pod_modes);
    }
    else{
        window->hlpod_mat->pod_modes = NULL;
    }
    window->owns_pod_modes = 1;
    window->hlpod_mat->pod_modes = strom_alloc_matrix(
        window->st_dof,
        selected_modes);

    if(window->hlpod_mat->pod_modes == NULL){
        strom_free_matrix(&left_vectors);
        strom_free_matrix(&right_vectors_t);
        free(singular_values);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < window->st_dof; i++){
        for(int j = 0; j < selected_modes; j++){
            window->hlpod_mat->pod_modes[i][j] = left_vectors[i][j];
        }
    }

    if(window->owns_singular_values){
        free(window->singular_values);
    }
    window->singular_values = NULL;
    window->owns_singular_values = 1;
    window->singular_values = strom_alloc_vector(max_available_modes);
    if(window->singular_values == NULL){
        strom_free_matrix(&left_vectors);
        strom_free_matrix(&right_vectors_t);
        free(singular_values);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < max_available_modes; i++){
        window->singular_values[i] = singular_values[i];
    }

    free(window->hlpod_mat->node_id);
    free(window->hlpod_mat->num_modes_internal);
    free(window->hlpod_mat->n_internal_vertex_subd);

    window->hlpod_mat->node_id = strom_alloc_int_vector(window->st_dof);
    window->hlpod_mat->num_modes_internal = strom_alloc_int_vector(1);
    window->hlpod_mat->n_internal_vertex_subd = strom_alloc_int_vector(1);

    if(window->hlpod_mat->node_id == NULL ||
       window->hlpod_mat->num_modes_internal == NULL ||
       window->hlpod_mat->n_internal_vertex_subd == NULL){
        strom_free_matrix(&left_vectors);
        strom_free_matrix(&right_vectors_t);
        free(singular_values);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < window->st_dof; i++){
        window->hlpod_mat->node_id[i] = i;
    }

    window->num_modes = selected_modes;
    window->num_singular_values = max_available_modes;
    window->hlpod_vals->num_modes = selected_modes;
    window->hlpod_mat->num_modes_internal[0] = selected_modes;
    window->hlpod_mat->n_internal_vertex_subd[0] = window->st_dof;

    strom_zero_vector(window->reduced_rhs, window->num_modes_max);
    strom_zero_vector(window->reduced_coef, window->num_modes_max);

    strom_free_matrix(&left_vectors);
    strom_free_matrix(&right_vectors_t);
    free(singular_values);

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_write_pod_modes_binary(
    const ROM_STROM_WINDOW* window,
    const char* filename)
{
    FILE* fp;
    uint32_t header[5];

    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       filename == NULL || window->num_modes <= 0 ||
       window->hlpod_mat->pod_modes == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    fp = fopen(filename, "wb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    header[0] = ROM_STROM_BINARY_MAGIC;
    header[1] = ROM_STROM_BINARY_VERSION;
    header[2] = (uint32_t)window->st_dof;
    header[3] = (uint32_t)window->num_modes;
    header[4] = (uint32_t)window->num_singular_values;

    if(fwrite(header, sizeof(uint32_t), 5, fp) != 5){
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    for(int i = 0; i < window->st_dof; i++){
        if(fwrite(
            window->hlpod_mat->pod_modes[i],
            sizeof(double),
            (size_t)window->num_modes,
            fp) != (size_t)window->num_modes){
            fclose(fp);
            return ROM_STROM_ERR_IO;
        }
    }

    if(window->num_singular_values > 0 &&
       fwrite(
           window->singular_values,
           sizeof(double),
           (size_t)window->num_singular_values,
           fp) != (size_t)window->num_singular_values){
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    fclose(fp);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_read_pod_modes_binary(
    ROM_STROM_WINDOW* window,
    const char* filename)
{
    FILE* fp;
    uint32_t header[5];
    int st_dof;
    int num_modes;
    int num_singular_values;

    if(strom_validate_window(window) != ROM_STROM_SUCCESS || filename == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    if(fread(header, sizeof(uint32_t), 5, fp) != 5 ||
       header[0] != ROM_STROM_BINARY_MAGIC ||
       header[1] != ROM_STROM_BINARY_VERSION){
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    st_dof = (int)header[2];
    num_modes = (int)header[3];
    num_singular_values = (int)header[4];

    if(st_dof != window->st_dof || num_modes <= 0 ||
       num_modes > window->num_modes_max){
        fclose(fp);
        return ROM_STROM_ERR_DIMENSION;
    }

    if(window->owns_pod_modes){
        strom_free_matrix(&window->hlpod_mat->pod_modes);
    }
    else{
        window->hlpod_mat->pod_modes = NULL;
    }
    window->owns_pod_modes = 1;
    window->hlpod_mat->pod_modes = strom_alloc_matrix(st_dof, num_modes);

    if(window->owns_singular_values){
        free(window->singular_values);
    }
    window->singular_values = NULL;
    window->owns_singular_values = 1;
    window->singular_values = num_singular_values > 0
        ? strom_alloc_vector(num_singular_values)
        : NULL;

    if(window->hlpod_mat->pod_modes == NULL ||
       (num_singular_values > 0 && window->singular_values == NULL)){
        fclose(fp);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < st_dof; i++){
        if(fread(
            window->hlpod_mat->pod_modes[i],
            sizeof(double),
            (size_t)num_modes,
            fp) != (size_t)num_modes){
            fclose(fp);
            return ROM_STROM_ERR_IO;
        }
    }

    if(num_singular_values > 0 &&
       fread(
           window->singular_values,
           sizeof(double),
           (size_t)num_singular_values,
           fp) != (size_t)num_singular_values){
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    fclose(fp);

    free(window->hlpod_mat->node_id);
    free(window->hlpod_mat->num_modes_internal);
    free(window->hlpod_mat->n_internal_vertex_subd);

    window->hlpod_mat->node_id = strom_alloc_int_vector(window->st_dof);
    window->hlpod_mat->num_modes_internal = strom_alloc_int_vector(1);
    window->hlpod_mat->n_internal_vertex_subd = strom_alloc_int_vector(1);

    if(window->hlpod_mat->node_id == NULL ||
       window->hlpod_mat->num_modes_internal == NULL ||
       window->hlpod_mat->n_internal_vertex_subd == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < window->st_dof; i++){
        window->hlpod_mat->node_id[i] = i;
    }

    window->num_modes = num_modes;
    window->num_singular_values = num_singular_values;
    window->hlpod_vals->num_modes = num_modes;
    window->hlpod_mat->num_modes_internal[0] = num_modes;
    window->hlpod_mat->n_internal_vertex_subd[0] = window->st_dof;

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_copy_pod_basis(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source)
{
    double** new_pod_modes = NULL;
    double* new_singular_values = NULL;
    int* new_node_id = NULL;
    int* new_num_modes_internal = NULL;
    int* new_n_internal_vertex_subd = NULL;

    if(strom_validate_window(destination) != ROM_STROM_SUCCESS ||
       strom_validate_window(source) != ROM_STROM_SUCCESS ||
       source->num_modes <= 0 ||
       source->hlpod_mat->pod_modes == NULL ||
       (source->num_singular_values > 0 &&
        source->singular_values == NULL))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(destination->st_dof != source->st_dof ||
       destination->num_modes_max < source->num_modes)
    {
        return ROM_STROM_ERR_DIMENSION;
    }

    if(destination == source){
        return ROM_STROM_SUCCESS;
    }

    new_pod_modes = strom_alloc_matrix(
        source->st_dof,
        source->num_modes);

    if(source->num_singular_values > 0){
        new_singular_values = strom_alloc_vector(
            source->num_singular_values);
    }

    new_node_id = strom_alloc_int_vector(source->st_dof);
    new_num_modes_internal = strom_alloc_int_vector(1);
    new_n_internal_vertex_subd = strom_alloc_int_vector(1);

    if(new_pod_modes == NULL ||
       (source->num_singular_values > 0 &&
        new_singular_values == NULL) ||
       new_node_id == NULL ||
       new_num_modes_internal == NULL ||
       new_n_internal_vertex_subd == NULL)
    {
        strom_free_matrix(&new_pod_modes);
        free(new_singular_values);
        free(new_node_id);
        free(new_num_modes_internal);
        free(new_n_internal_vertex_subd);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < source->st_dof; i++){
        for(int j = 0; j < source->num_modes; j++){
            new_pod_modes[i][j] =
                source->hlpod_mat->pod_modes[i][j];
        }
        new_node_id[i] = i;
    }

    for(int i = 0; i < source->num_singular_values; i++){
        new_singular_values[i] = source->singular_values[i];
    }

    new_num_modes_internal[0] = source->num_modes;
    new_n_internal_vertex_subd[0] = source->st_dof;

    if(destination->owns_pod_modes){
        strom_free_matrix(&destination->hlpod_mat->pod_modes);
    }
    else{
        destination->hlpod_mat->pod_modes = NULL;
    }
    if(destination->owns_singular_values){
        free(destination->singular_values);
    }
    destination->singular_values = NULL;
    free(destination->hlpod_mat->node_id);
    free(destination->hlpod_mat->num_modes_internal);
    free(destination->hlpod_mat->n_internal_vertex_subd);

    destination->hlpod_mat->pod_modes = new_pod_modes;
    destination->singular_values = new_singular_values;
    destination->owns_pod_modes = 1;
    destination->owns_singular_values = 1;
    destination->hlpod_mat->node_id = new_node_id;
    destination->hlpod_mat->num_modes_internal =
        new_num_modes_internal;
    destination->hlpod_mat->n_internal_vertex_subd =
        new_n_internal_vertex_subd;

    destination->num_modes = source->num_modes;
    destination->num_singular_values =
        source->num_singular_values;
    destination->hlpod_vals->num_modes = source->num_modes;

    /* Any old reduced blocks no longer correspond to the copied basis. */
    strom_free_matrix(&destination->reduced_diag);
    strom_free_matrix(&destination->reduced_coupling_prev);
    strom_zero_vector(
        destination->reduced_rhs,
        destination->num_modes_max);
    strom_zero_vector(
        destination->reduced_coef,
        destination->num_modes_max);

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_attach_pod_basis(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source)
{
    int* new_node_id = NULL;
    int* new_num_modes_internal = NULL;
    int* new_n_internal_vertex_subd = NULL;

    if(strom_validate_window(destination) != ROM_STROM_SUCCESS ||
       strom_validate_window(source) != ROM_STROM_SUCCESS ||
       source->num_modes <= 0 ||
       source->hlpod_mat->pod_modes == NULL ||
       destination->st_dof != source->st_dof ||
       destination->num_modes_max < source->num_modes)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(destination == source){
        return ROM_STROM_SUCCESS;
    }

    new_node_id = strom_alloc_int_vector(source->st_dof);
    new_num_modes_internal = strom_alloc_int_vector(1);
    new_n_internal_vertex_subd = strom_alloc_int_vector(1);

    if(new_node_id == NULL ||
       new_num_modes_internal == NULL ||
       new_n_internal_vertex_subd == NULL)
    {
        free(new_node_id);
        free(new_num_modes_internal);
        free(new_n_internal_vertex_subd);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < source->st_dof; i++){
        new_node_id[i] = i;
    }
    new_num_modes_internal[0] = source->num_modes;
    new_n_internal_vertex_subd[0] = source->st_dof;

    if(destination->owns_pod_modes){
        strom_free_matrix(&destination->hlpod_mat->pod_modes);
    }
    if(destination->owns_singular_values){
        free(destination->singular_values);
    }

    free(destination->hlpod_mat->node_id);
    free(destination->hlpod_mat->num_modes_internal);
    free(destination->hlpod_mat->n_internal_vertex_subd);

    destination->hlpod_mat->pod_modes = source->hlpod_mat->pod_modes;
    destination->singular_values = source->singular_values;
    destination->owns_pod_modes = 0;
    destination->owns_singular_values = 0;
    destination->hlpod_mat->node_id = new_node_id;
    destination->hlpod_mat->num_modes_internal = new_num_modes_internal;
    destination->hlpod_mat->n_internal_vertex_subd =
        new_n_internal_vertex_subd;

    destination->num_modes = source->num_modes;
    destination->num_singular_values = source->num_singular_values;
    destination->hlpod_vals->num_modes = source->num_modes;

    strom_free_matrix(&destination->reduced_diag);
    strom_free_matrix(&destination->reduced_coupling_prev);
    strom_zero_vector(destination->reduced_rhs, destination->num_modes_max);
    strom_zero_vector(destination->reduced_coef, destination->num_modes_max);

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_calc_reduced_diag(
    ROM_STROM_WINDOW* window,
    ROM_STROM_APPLY_WINDOW_OPERATOR apply_A,
    void* operator_context)
{
    double* x;
    double* y;

    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       apply_A == NULL || window->num_modes <= 0 ||
       window->hlpod_mat->pod_modes == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    x = strom_alloc_vector(window->st_dof);
    y = strom_alloc_vector(window->st_dof);

    if(x == NULL || y == NULL){
        free(x);
        free(y);
        return ROM_STROM_ERR_ALLOCATION;
    }

    strom_free_matrix(&window->reduced_diag);
    window->reduced_diag = strom_alloc_matrix(
        window->num_modes,
        window->num_modes);

    if(window->reduced_diag == NULL){
        free(x);
        free(y);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int j = 0; j < window->num_modes; j++){
        for(int p = 0; p < window->st_dof; p++){
            x[p] = window->hlpod_mat->pod_modes[p][j];
        }
        strom_zero_vector(y, window->st_dof);

        if(apply_A(window->window_id, x, y, operator_context) != 0){
            free(x);
            free(y);
            return ROM_STROM_ERR_OPERATOR;
        }

        for(int i = 0; i < window->num_modes; i++){
            double value = 0.0;
            for(int p = 0; p < window->st_dof; p++){
                value += window->hlpod_mat->pod_modes[p][i] * y[p];
            }
            window->reduced_diag[i][j] = value;
        }
    }

    free(x);
    free(y);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_calc_reduced_coupling(
    ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context)
{
    double* x_prev;
    double* y_cur;

    if(strom_validate_window(current) != ROM_STROM_SUCCESS ||
       strom_validate_window(previous) != ROM_STROM_SUCCESS ||
       apply_B == NULL || current->num_modes <= 0 ||
       previous->num_modes <= 0 ||
       current->hlpod_mat->pod_modes == NULL ||
       previous->hlpod_mat->pod_modes == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    x_prev = strom_alloc_vector(previous->st_dof);
    y_cur = strom_alloc_vector(current->st_dof);

    if(x_prev == NULL || y_cur == NULL){
        free(x_prev);
        free(y_cur);
        return ROM_STROM_ERR_ALLOCATION;
    }

    strom_free_matrix(&current->reduced_coupling_prev);
    current->reduced_coupling_prev = strom_alloc_matrix(
        current->num_modes,
        previous->num_modes);

    if(current->reduced_coupling_prev == NULL){
        free(x_prev);
        free(y_cur);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int j = 0; j < previous->num_modes; j++){
        for(int p = 0; p < previous->st_dof; p++){
            x_prev[p] = previous->hlpod_mat->pod_modes[p][j];
        }
        strom_zero_vector(y_cur, current->st_dof);

        if(apply_B(
            current->window_id,
            x_prev,
            y_cur,
            coupling_context) != 0){
            free(x_prev);
            free(y_cur);
            return ROM_STROM_ERR_OPERATOR;
        }

        for(int i = 0; i < current->num_modes; i++){
            double value = 0.0;
            for(int p = 0; p < current->st_dof; p++){
                value += current->hlpod_mat->pod_modes[p][i] * y_cur[p];
            }
            current->reduced_coupling_prev[i][j] = value;
        }
    }

    free(x_prev);
    free(y_cur);
    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_calc_reduced_diag_prepared(
    ROM_STROM_WINDOW* window,
    ROM_STROM_PREPARE_WINDOW_OPERATOR prepare_A,
    ROM_STROM_APPLY_PREPARED_OPERATOR apply_prepared_A,
    void* operator_context)
{
    const int n = window != NULL ? window->st_dof : 0;
    const int r = window != NULL ? window->num_modes : 0;
    double* a_phi = NULL;
    double* x = NULL;
    double* y = NULL;
    const double* phi = NULL;
    int status;

    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       prepare_A == NULL ||
       apply_prepared_A == NULL ||
       r <= 0 ||
       window->hlpod_mat->pod_modes == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    phi = window->hlpod_mat->pod_modes[0];
    a_phi = strom_alloc_vector(n * r);
    x = strom_alloc_vector(n);
    y = strom_alloc_vector(n);

    if(a_phi == NULL || x == NULL || y == NULL){
        free(a_phi);
        free(x);
        free(y);
        return ROM_STROM_ERR_ALLOCATION;
    }

    status = prepare_A(window->window_id, operator_context);
    if(status != ROM_STROM_SUCCESS){
        free(a_phi);
        free(x);
        free(y);
        return ROM_STROM_ERR_OPERATOR;
    }

    for(int j = 0; j < r; j++){
#ifdef ROM_STROM_USE_CBLAS
        cblas_dcopy(n, phi + j, r, x, 1);
#else
        for(int p = 0; p < n; p++){
            x[p] = phi[(size_t)p * (size_t)r + (size_t)j];
        }
#endif
        status = apply_prepared_A(x, y, operator_context);
        if(status != ROM_STROM_SUCCESS){
            free(a_phi);
            free(x);
            free(y);
            return ROM_STROM_ERR_OPERATOR;
        }
#ifdef ROM_STROM_USE_CBLAS
        cblas_dcopy(n, y, 1, a_phi + j, r);
#else
        for(int p = 0; p < n; p++){
            a_phi[(size_t)p * (size_t)r + (size_t)j] = y[p];
        }
#endif
    }

    strom_free_matrix(&window->reduced_diag);
    window->reduced_diag = strom_alloc_matrix(r, r);
    if(window->reduced_diag == NULL){
        free(a_phi);
        free(x);
        free(y);
        return ROM_STROM_ERR_ALLOCATION;
    }

#ifdef ROM_STROM_USE_CBLAS
    cblas_dgemm(
        CblasRowMajor,
        CblasTrans,
        CblasNoTrans,
        r,
        r,
        n,
        1.0,
        phi,
        r,
        a_phi,
        r,
        0.0,
        window->reduced_diag[0],
        r);
#else
    for(int i = 0; i < r; i++){
        for(int j = 0; j < r; j++){
            double value = 0.0;
            for(int p = 0; p < n; p++){
                value += phi[(size_t)p * (size_t)r + (size_t)i]
                       * a_phi[(size_t)p * (size_t)r + (size_t)j];
            }
            window->reduced_diag[i][j] = value;
        }
    }
#endif

    free(a_phi);
    free(x);
    free(y);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_calc_reduced_coupling_fast(
    ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context)
{
    const int n_cur = current != NULL ? current->st_dof : 0;
    const int n_prev = previous != NULL ? previous->st_dof : 0;
    const int r_cur = current != NULL ? current->num_modes : 0;
    const int r_prev = previous != NULL ? previous->num_modes : 0;
    const double* phi_cur = NULL;
    const double* phi_prev = NULL;
    double* b_phi = NULL;
    double* x_prev = NULL;
    double* y_cur = NULL;

    if(strom_validate_window(current) != ROM_STROM_SUCCESS ||
       strom_validate_window(previous) != ROM_STROM_SUCCESS ||
       apply_B == NULL ||
       r_cur <= 0 || r_prev <= 0 ||
       current->hlpod_mat->pod_modes == NULL ||
       previous->hlpod_mat->pod_modes == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    phi_cur = current->hlpod_mat->pod_modes[0];
    phi_prev = previous->hlpod_mat->pod_modes[0];
    b_phi = strom_alloc_vector(n_cur * r_prev);
    x_prev = strom_alloc_vector(n_prev);
    y_cur = strom_alloc_vector(n_cur);

    if(b_phi == NULL || x_prev == NULL || y_cur == NULL){
        free(b_phi);
        free(x_prev);
        free(y_cur);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int j = 0; j < r_prev; j++){
#ifdef ROM_STROM_USE_CBLAS
        cblas_dcopy(n_prev, phi_prev + j, r_prev, x_prev, 1);
#else
        for(int p = 0; p < n_prev; p++){
            x_prev[p] =
                phi_prev[(size_t)p * (size_t)r_prev + (size_t)j];
        }
#endif
        strom_zero_vector(y_cur, n_cur);
        if(apply_B(
            current->window_id,
            x_prev,
            y_cur,
            coupling_context) != ROM_STROM_SUCCESS)
        {
            free(b_phi);
            free(x_prev);
            free(y_cur);
            return ROM_STROM_ERR_OPERATOR;
        }
#ifdef ROM_STROM_USE_CBLAS
        cblas_dcopy(n_cur, y_cur, 1, b_phi + j, r_prev);
#else
        for(int p = 0; p < n_cur; p++){
            b_phi[(size_t)p * (size_t)r_prev + (size_t)j] = y_cur[p];
        }
#endif
    }

    strom_free_matrix(&current->reduced_coupling_prev);
    current->reduced_coupling_prev = strom_alloc_matrix(r_cur, r_prev);
    if(current->reduced_coupling_prev == NULL){
        free(b_phi);
        free(x_prev);
        free(y_cur);
        return ROM_STROM_ERR_ALLOCATION;
    }

#ifdef ROM_STROM_USE_CBLAS
    cblas_dgemm(
        CblasRowMajor,
        CblasTrans,
        CblasNoTrans,
        r_cur,
        r_prev,
        n_cur,
        1.0,
        phi_cur,
        r_cur,
        b_phi,
        r_prev,
        0.0,
        current->reduced_coupling_prev[0],
        r_prev);
#else
    for(int i = 0; i < r_cur; i++){
        for(int j = 0; j < r_prev; j++){
            double value = 0.0;
            for(int p = 0; p < n_cur; p++){
                value += phi_cur[(size_t)p * (size_t)r_cur + (size_t)i]
                       * b_phi[(size_t)p * (size_t)r_prev + (size_t)j];
            }
            current->reduced_coupling_prev[i][j] = value;
        }
    }
#endif

    free(b_phi);
    free(x_prev);
    free(y_cur);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_copy_reduced_diag(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source)
{
    if(strom_validate_window(destination) != ROM_STROM_SUCCESS ||
       strom_validate_window(source) != ROM_STROM_SUCCESS ||
       source->num_modes <= 0 ||
       destination->num_modes != source->num_modes ||
       source->reduced_diag == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_free_matrix(&destination->reduced_diag);
    destination->reduced_diag = strom_alloc_matrix(
        source->num_modes,
        source->num_modes);
    if(destination->reduced_diag == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    memcpy(
        destination->reduced_diag[0],
        source->reduced_diag[0],
        (size_t)source->num_modes * (size_t)source->num_modes
            * sizeof(double));
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_copy_reduced_coupling(
    ROM_STROM_WINDOW* destination,
    const ROM_STROM_WINDOW* source)
{
    if(strom_validate_window(destination) != ROM_STROM_SUCCESS ||
       strom_validate_window(source) != ROM_STROM_SUCCESS ||
       destination->num_modes != source->num_modes ||
       source->reduced_coupling_prev == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_free_matrix(&destination->reduced_coupling_prev);
    destination->reduced_coupling_prev = strom_alloc_matrix(
        source->num_modes,
        source->num_modes);
    if(destination->reduced_coupling_prev == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    memcpy(
        destination->reduced_coupling_prev[0],
        source->reduced_coupling_prev[0],
        (size_t)source->num_modes * (size_t)source->num_modes
            * sizeof(double));
    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_set_reduced_diag_data(
    ROM_STROM_WINDOW* window,
    const double* row_major_data)
{
    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       window->num_modes <= 0 ||
       row_major_data == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_free_matrix(&window->reduced_diag);
    window->reduced_diag = strom_alloc_matrix(
        window->num_modes,
        window->num_modes);
    if(window->reduced_diag == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    memcpy(
        window->reduced_diag[0],
        row_major_data,
        (size_t)window->num_modes * (size_t)window->num_modes
            * sizeof(double));
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_set_reduced_coupling_data(
    ROM_STROM_WINDOW* window,
    const double* row_major_data,
    int previous_num_modes)
{
    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       window->num_modes <= 0 ||
       previous_num_modes <= 0 ||
       row_major_data == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_free_matrix(&window->reduced_coupling_prev);
    window->reduced_coupling_prev = strom_alloc_matrix(
        window->num_modes,
        previous_num_modes);
    if(window->reduced_coupling_prev == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    memcpy(
        window->reduced_coupling_prev[0],
        row_major_data,
        (size_t)window->num_modes * (size_t)previous_num_modes
            * sizeof(double));
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_project_rhs(
    ROM_STROM_WINDOW* window,
    const double* effective_full_rhs)
{
    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       effective_full_rhs == NULL || window->num_modes <= 0 ||
       window->hlpod_mat->pod_modes == NULL ||
       window->reduced_rhs == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_zero_vector(window->reduced_rhs, window->num_modes_max);

#ifdef ROM_STROM_USE_CBLAS
    cblas_dgemv(
        CblasRowMajor,
        CblasTrans,
        window->st_dof,
        window->num_modes,
        1.0,
        window->hlpod_mat->pod_modes[0],
        window->num_modes,
        effective_full_rhs,
        1,
        0.0,
        window->reduced_rhs,
        1);
#else
    for(int i = 0; i < window->num_modes; i++){
        double value = 0.0;
        for(int p = 0; p < window->st_dof; p++){
            value += window->hlpod_mat->pod_modes[p][i]
                   * effective_full_rhs[p];
        }
        window->reduced_rhs[i] = value;
    }
#endif

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_build_effective_rhs(
    const ROM_STROM_WINDOW* current,
    const ROM_STROM_WINDOW* previous,
    const double* full_rhs,
    double* effective_rhs,
    ROM_STROM_APPLY_WINDOW_OPERATOR apply_A,
    void* operator_context,
    ROM_STROM_APPLY_WINDOW_COUPLING apply_B,
    void* coupling_context)
{
    double* work;

    if(strom_validate_window(current) != ROM_STROM_SUCCESS ||
       full_rhs == NULL || effective_rhs == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int i = 0; i < current->st_dof; i++){
        effective_rhs[i] = full_rhs[i];
    }

    work = strom_alloc_vector(current->st_dof);
    if(work == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    if(current->reference_state != NULL){
        if(apply_A == NULL){
            free(work);
            return ROM_STROM_ERR_ARGUMENT;
        }

        if(apply_A(
            current->window_id,
            current->reference_state,
            work,
            operator_context) != 0){
            free(work);
            return ROM_STROM_ERR_OPERATOR;
        }

        for(int i = 0; i < current->st_dof; i++){
            effective_rhs[i] -= work[i];
        }
    }

    if(previous != NULL && previous->reference_state != NULL){
        if(apply_B == NULL){
            free(work);
            return ROM_STROM_ERR_ARGUMENT;
        }

        strom_zero_vector(work, current->st_dof);
        if(apply_B(
            current->window_id,
            previous->reference_state,
            work,
            coupling_context) != 0){
            free(work);
            return ROM_STROM_ERR_OPERATOR;
        }

        for(int i = 0; i < current->st_dof; i++){
            effective_rhs[i] += work[i];
        }
    }

    free(work);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_solve_dense_gauss(
    int n,
    double** A,
    const double* b,
    double* x,
    void* context)
{
    double** work_A;
    double* work_b;
    const double pivot_tolerance = 1.0e-14;

    (void)context;

    if(n <= 0 || A == NULL || b == NULL || x == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    work_A = strom_alloc_matrix(n, n);
    work_b = strom_alloc_vector(n);

    if(work_A == NULL || work_b == NULL){
        strom_free_matrix(&work_A);
        free(work_b);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < n; i++){
        work_b[i] = b[i];
        for(int j = 0; j < n; j++){
            work_A[i][j] = A[i][j];
        }
    }

    for(int k = 0; k < n; k++){
        int pivot_row = k;
        double pivot_abs = fabs(work_A[k][k]);

        for(int i = k + 1; i < n; i++){
            const double candidate = fabs(work_A[i][k]);
            if(candidate > pivot_abs){
                pivot_abs = candidate;
                pivot_row = i;
            }
        }

        if(pivot_abs < pivot_tolerance){
            strom_free_matrix(&work_A);
            free(work_b);
            return ROM_STROM_ERR_SOLVER;
        }

        if(pivot_row != k){
            double* row_tmp = work_A[k];
            double rhs_tmp = work_b[k];
            work_A[k] = work_A[pivot_row];
            work_A[pivot_row] = row_tmp;
            work_b[k] = work_b[pivot_row];
            work_b[pivot_row] = rhs_tmp;
        }

        for(int i = k + 1; i < n; i++){
            const double factor = work_A[i][k] / work_A[k][k];
            work_A[i][k] = 0.0;

            for(int j = k + 1; j < n; j++){
                work_A[i][j] -= factor * work_A[k][j];
            }
            work_b[i] -= factor * work_b[k];
        }
    }

    for(int i = n - 1; i >= 0; i--){
        double value = work_b[i];
        for(int j = i + 1; j < n; j++){
            value -= work_A[i][j] * x[j];
        }
        x[i] = value / work_A[i][i];
    }

    strom_free_matrix(&work_A);
    free(work_b);
    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_solve_dense_lapack(
    int n,
    double** A,
    const double* b,
    double* x,
    void* context)
{
    (void)context;

#ifndef ROM_STROM_USE_LAPACK
    return ROM_std_strom_solve_dense_gauss(n, A, b, x, NULL);
#else
    const int nrhs = 1;
    const int lda = n;
    const int ldb = n;
    int info = 0;
    int* pivots = NULL;
    double* a_column_major = NULL;

    if(n <= 0 || A == NULL || b == NULL || x == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    pivots = strom_alloc_int_vector(n);
    a_column_major = strom_alloc_vector(n * n);
    if(pivots == NULL || a_column_major == NULL){
        free(pivots);
        free(a_column_major);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int i = 0; i < n; i++){
        x[i] = b[i];
        for(int j = 0; j < n; j++){
            a_column_major[(size_t)j * (size_t)n + (size_t)i] = A[i][j];
        }
    }

    dgesv_(
        &n,
        &nrhs,
        a_column_major,
        &lda,
        pivots,
        x,
        &ldb,
        &info);

    free(pivots);
    free(a_column_major);

    return info == 0 ? ROM_STROM_SUCCESS : ROM_STROM_ERR_SOLVER;
#endif
}

static ROM_STROM_DENSE_SOLVER strom_default_dense_solver(void)
{
#ifdef ROM_STROM_USE_LAPACK
    return ROM_std_strom_solve_dense_lapack;
#else
    return ROM_std_strom_solve_dense_gauss;
#endif
}

int ROM_std_strom_solve_sequential(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_DENSE_SOLVER dense_solver,
    void* solver_context)
{
    int max_modes = 0;
    double* rhs_work;

    if(system == NULL || system->windows == NULL || system->num_windows <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(dense_solver == NULL){
        dense_solver = strom_default_dense_solver();
    }

    for(int w = 0; w < system->num_windows; w++){
        if(system->windows[w].num_modes > max_modes){
            max_modes = system->windows[w].num_modes;
        }
    }

    rhs_work = strom_alloc_vector(max_modes);
    if(rhs_work == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int w = 0; w < system->num_windows; w++){
        ROM_STROM_WINDOW* current = &system->windows[w];

        if(current->reduced_diag == NULL ||
           current->reduced_rhs == NULL ||
           current->reduced_coef == NULL ||
           current->num_modes <= 0){
            free(rhs_work);
            return ROM_STROM_ERR_ARGUMENT;
        }

        for(int i = 0; i < current->num_modes; i++){
            rhs_work[i] = current->reduced_rhs[i];
        }

        if(w > 0){
            const ROM_STROM_WINDOW* previous = &system->windows[w - 1];

            if(current->reduced_coupling_prev == NULL){
                free(rhs_work);
                return ROM_STROM_ERR_ARGUMENT;
            }

            for(int i = 0; i < current->num_modes; i++){
                for(int j = 0; j < previous->num_modes; j++){
                    rhs_work[i] +=
                        current->reduced_coupling_prev[i][j]
                        * previous->reduced_coef[j];
                }
            }
        }

        strom_zero_vector(current->reduced_coef, current->num_modes_max);

        if(dense_solver(
            current->num_modes,
            current->reduced_diag,
            rhs_work,
            current->reduced_coef,
            solver_context) != ROM_STROM_SUCCESS){
            free(rhs_work);
            return ROM_STROM_ERR_SOLVER;
        }
    }

    free(rhs_work);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_assemble_all_at_once(
    const ROM_STROM_SYSTEM* system,
    ROM_STROM_ADD_MATRIX_ENTRY add_matrix_entry,
    ROM_STROM_SET_VECTOR_ENTRY set_rhs_entry,
    void* assembly_context)
{
    if(system == NULL || system->windows == NULL ||
       system->mode_offset == NULL || system->num_windows <= 0 ||
       add_matrix_entry == NULL || set_rhs_entry == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int w = 0; w < system->num_windows; w++){
        const ROM_STROM_WINDOW* current = &system->windows[w];
        const int row_offset = system->mode_offset[w];

        if(current->reduced_diag == NULL || current->reduced_rhs == NULL){
            return ROM_STROM_ERR_ARGUMENT;
        }

        for(int i = 0; i < current->num_modes; i++){
            for(int j = 0; j < current->num_modes; j++){
                if(add_matrix_entry(
                    row_offset + i,
                    row_offset + j,
                    current->reduced_diag[i][j],
                    assembly_context) != 0){
                    return ROM_STROM_ERR_OPERATOR;
                }
            }

            if(set_rhs_entry(
                row_offset + i,
                current->reduced_rhs[i],
                assembly_context) != 0){
                return ROM_STROM_ERR_OPERATOR;
            }
        }

        if(w > 0){
            const ROM_STROM_WINDOW* previous = &system->windows[w - 1];
            const int column_offset = system->mode_offset[w - 1];

            if(current->reduced_coupling_prev == NULL){
                return ROM_STROM_ERR_ARGUMENT;
            }

            for(int i = 0; i < current->num_modes; i++){
                for(int j = 0; j < previous->num_modes; j++){
                    if(add_matrix_entry(
                        row_offset + i,
                        column_offset + j,
                        -current->reduced_coupling_prev[i][j],
                        assembly_context) != 0){
                        return ROM_STROM_ERR_OPERATOR;
                    }
                }
            }
        }
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_solve_all_at_once_dense(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_DENSE_SOLVER dense_solver,
    void* solver_context)
{
    double** matrix;
    double* rhs;
    double* solution;
    int status;

    if(system == NULL || system->windows == NULL ||
       system->mode_offset == NULL || system->total_reduced_dof <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(dense_solver == NULL){
        dense_solver = strom_default_dense_solver();
    }

    matrix = strom_alloc_matrix(
        system->total_reduced_dof,
        system->total_reduced_dof);
    rhs = strom_alloc_vector(system->total_reduced_dof);
    solution = strom_alloc_vector(system->total_reduced_dof);

    if(matrix == NULL || rhs == NULL || solution == NULL){
        strom_free_matrix(&matrix);
        free(rhs);
        free(solution);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int w = 0; w < system->num_windows; w++){
        ROM_STROM_WINDOW* current = &system->windows[w];
        const int row_offset = system->mode_offset[w];

        for(int i = 0; i < current->num_modes; i++){
            rhs[row_offset + i] = current->reduced_rhs[i];

            for(int j = 0; j < current->num_modes; j++){
                matrix[row_offset + i][row_offset + j] =
                    current->reduced_diag[i][j];
            }
        }

        if(w > 0){
            ROM_STROM_WINDOW* previous = &system->windows[w - 1];
            const int column_offset = system->mode_offset[w - 1];

            for(int i = 0; i < current->num_modes; i++){
                for(int j = 0; j < previous->num_modes; j++){
                    matrix[row_offset + i][column_offset + j] =
                        -current->reduced_coupling_prev[i][j];
                }
            }
        }
    }

    status = dense_solver(
        system->total_reduced_dof,
        matrix,
        rhs,
        solution,
        solver_context);

    if(status != ROM_STROM_SUCCESS){
        strom_free_matrix(&matrix);
        free(rhs);
        free(solution);
        return ROM_STROM_ERR_SOLVER;
    }

    for(int w = 0; w < system->num_windows; w++){
        ROM_STROM_WINDOW* current = &system->windows[w];
        const int offset = system->mode_offset[w];

        for(int i = 0; i < current->num_modes; i++){
            current->reduced_coef[i] = solution[offset + i];
        }
    }

    strom_free_matrix(&matrix);
    free(rhs);
    free(solution);
    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_project_window_state(
    const ROM_STROM_WINDOW* window,
    const double* window_state,
    double* projected_state)
{
    double* coefficients;

    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       window_state == NULL ||
       projected_state == NULL ||
       window->num_modes <= 0 ||
       window->hlpod_mat->pod_modes == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    coefficients = strom_alloc_vector(window->num_modes);
    if(coefficients == NULL){
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int j = 0; j < window->num_modes; j++){
        double coefficient = 0.0;

        for(int i = 0; i < window->st_dof; i++){
            const double centered_value =
                window_state[i]
                - (window->reference_state != NULL
                   ? window->reference_state[i]
                   : 0.0);

            coefficient +=
                window->hlpod_mat->pod_modes[i][j]
                * centered_value;
        }

        coefficients[j] = coefficient;
    }

    for(int i = 0; i < window->st_dof; i++){
        double value = window->reference_state != NULL
            ? window->reference_state[i]
            : 0.0;

        for(int j = 0; j < window->num_modes; j++){
            value +=
                window->hlpod_mat->pod_modes[i][j]
                * coefficients[j];
        }

        projected_state[i] = value;
    }

    free(coefficients);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_reconstruct_window(
    const ROM_STROM_WINDOW* window,
    double* window_solution)
{
    if(strom_validate_window(window) != ROM_STROM_SUCCESS ||
       window_solution == NULL || window->num_modes <= 0 ||
       window->hlpod_mat->pod_modes == NULL ||
       window->reduced_coef == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

#ifdef ROM_STROM_USE_CBLAS
    if(window->reference_state != NULL){
        memcpy(
            window_solution,
            window->reference_state,
            (size_t)window->st_dof * sizeof(double));
    }
    else{
        strom_zero_vector(window_solution, window->st_dof);
    }

    cblas_dgemv(
        CblasRowMajor,
        CblasNoTrans,
        window->st_dof,
        window->num_modes,
        1.0,
        window->hlpod_mat->pod_modes[0],
        window->num_modes,
        window->reduced_coef,
        1,
        1.0,
        window_solution,
        1);
#else
    for(int p = 0; p < window->st_dof; p++){
        double value = window->reference_state != NULL
            ? window->reference_state[p]
            : 0.0;

        for(int j = 0; j < window->num_modes; j++){
            value += window->hlpod_mat->pod_modes[p][j]
                   * window->reduced_coef[j];
        }

        window_solution[p] = value;
    }
#endif

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_extract_right_trace(
    const double* window_vector,
    double* trace,
    int spatial_dof,
    int num_slabs,
    int num_time_basis,
    const double* theta_plus)
{
    const int last_slab = num_slabs - 1;

    if(window_vector == NULL || trace == NULL || theta_plus == NULL ||
       spatial_dof <= 0 || num_slabs <= 0 || num_time_basis <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int p = 0; p < spatial_dof; p++){
        double value = 0.0;

        for(int alpha = 0; alpha < num_time_basis; alpha++){
            const int row =
                (last_slab * num_time_basis + alpha) * spatial_dof + p;

            value += theta_plus[alpha] * window_vector[row];
        }

        trace[p] = value;
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_inject_left_trace(
    const double* spatial_trace,
    double* current_window_vector,
    int spatial_dof,
    int num_time_basis,
    const double* jump_left_coefficients)
{
    if(spatial_trace == NULL || current_window_vector == NULL ||
       jump_left_coefficients == NULL || spatial_dof <= 0 ||
       num_time_basis <= 0){
        return ROM_STROM_ERR_ARGUMENT;
    }

    for(int alpha = 0; alpha < num_time_basis; alpha++){
        for(int p = 0; p < spatial_dof; p++){
            const int row = alpha * spatial_dof + p;
            current_window_vector[row] +=
                jump_left_coefficients[alpha] * spatial_trace[p];
        }
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_apply_dg_coupling(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* context)
{
    ROM_STROM_DG_COUPLING_CONTEXT* dg_context;
    ROM_STROM_DG_TRANSITION* transition;
    double* trace;
    double* transformed_trace;
    int status;

    if(current_window_id <= 0 || x_prev == NULL || y_cur == NULL ||
       context == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    dg_context = (ROM_STROM_DG_COUPLING_CONTEXT*)context;

    if(current_window_id >= dg_context->num_windows ||
       dg_context->transition == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    transition = &dg_context->transition[current_window_id];

    if(transition->spatial_dof <= 0 ||
       transition->theta_plus_prev == NULL ||
       transition->jump_left_cur == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    trace = strom_alloc_vector(transition->spatial_dof);
    transformed_trace = strom_alloc_vector(transition->spatial_dof);

    if(trace == NULL || transformed_trace == NULL){
        free(trace);
        free(transformed_trace);
        return ROM_STROM_ERR_ALLOCATION;
    }

    status = ROM_std_strom_extract_right_trace(
        x_prev,
        trace,
        transition->spatial_dof,
        transition->prev_num_slabs,
        transition->prev_num_time_basis,
        transition->theta_plus_prev);

    if(status != ROM_STROM_SUCCESS){
        free(trace);
        free(transformed_trace);
        return status;
    }

    if(transition->apply_spatial_trace_operator != NULL){
        if(transition->apply_spatial_trace_operator(
            trace,
            transformed_trace,
            transition->spatial_operator_context) != 0){
            free(trace);
            free(transformed_trace);
            return ROM_STROM_ERR_OPERATOR;
        }
    }
    else{
        for(int p = 0; p < transition->spatial_dof; p++){
            transformed_trace[p] = trace[p];
        }
    }

    status = ROM_std_strom_inject_left_trace(
        transformed_trace,
        y_cur,
        transition->spatial_dof,
        transition->cur_num_time_basis,
        transition->jump_left_cur);

    free(trace);
    free(transformed_trace);
    return status;
}
