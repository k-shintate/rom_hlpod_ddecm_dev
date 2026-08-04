#include "rom_std_strom_monolis.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct {
    MONOLIS* monolis;
    int block_size;
    int total_reduced_dof;
} ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT;


static void strom_monolis_zero(
    ROM_STROM_MONOLIS_SOLVER* solver)
{
    if(solver == NULL){
        return;
    }

    memset(solver, 0, sizeof(*solver));
}


int ROM_std_strom_check_uniform_rank(
    const ROM_STROM_SYSTEM* system,
    int* block_size)
{
    int rank;

    if(system == NULL ||
       system->windows == NULL ||
       system->mode_offset == NULL ||
       block_size == NULL ||
       system->num_windows <= 0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    rank = system->windows[0].num_modes;

    if(rank <= 0){
        return ROM_STROM_ERR_DIMENSION;
    }

    for(int w = 0; w < system->num_windows; w++){
        const ROM_STROM_WINDOW* window = &system->windows[w];
        const int expected_offset = w * rank;

        if(window->num_modes != rank){
            fprintf(
                stderr,
                "ERROR: MONOLIS BCSR requires a uniform window rank.\n"
                "       window 0 rank: %d\n"
                "       window %d rank: %d\n",
                rank,
                w,
                window->num_modes);

            return ROM_STROM_ERR_DIMENSION;
        }

        if(system->mode_offset[w] != expected_offset){
            fprintf(
                stderr,
                "ERROR: reduced mode offset is incompatible with BCSR.\n"
                "       window: %d\n"
                "       actual offset: %d\n"
                "       expected offset: %d\n",
                w,
                system->mode_offset[w],
                expected_offset);

            return ROM_STROM_ERR_DIMENSION;
        }
    }

    if(system->total_reduced_dof != system->num_windows * rank){
        return ROM_STROM_ERR_DIMENSION;
    }

    *block_size = rank;
    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_build_bcsr_graph(
    const int num_windows,
    int** graph_index,
    int** graph_item,
    int* num_graph_items)
{
    int* index = NULL;
    int* item = NULL;
    int position = 0;
    int num_items;

    if(num_windows <= 0 ||
       graph_index == NULL ||
       graph_item == NULL ||
       num_graph_items == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    num_items = 2 * num_windows - 1;

    index = (int*)calloc((size_t)num_windows + 1, sizeof(int));
    item = (int*)calloc((size_t)num_items, sizeof(int));

    if(index == NULL || item == NULL){
        free(index);
        free(item);
        return ROM_STROM_ERR_ALLOCATION;
    }

    index[0] = 0;

    for(int w = 0; w < num_windows; w++){
        if(w > 0){
            item[position++] = w - 1;
        }

        /* Explicit self adjacency gives the dense D_w block. */
        item[position++] = w;
        index[w + 1] = position;
    }

    if(position != num_items){
        free(index);
        free(item);
        return ROM_STROM_ERR_DIMENSION;
    }

    *graph_index = index;
    *graph_item = item;
    *num_graph_items = num_items;

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_write_bcsr_graph(
    const ROM_STROM_MONOLIS_SOLVER* solver,
    const char* filename)
{
    FILE* fp;

    if(solver == NULL ||
       solver->graph_index == NULL ||
       solver->graph_item == NULL ||
       filename == NULL ||
       filename[0] == '\0')
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    fp = fopen(filename, "w");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    fprintf(fp, "%d\n", solver->num_windows);

    for(int w = 0; w < solver->num_windows; w++){
        const int begin = solver->graph_index[w];
        const int end = solver->graph_index[w + 1];

        fprintf(fp, "%d %d", w, end - begin);

        for(int k = begin; k < end; k++){
            fprintf(fp, " %d", solver->graph_item[k]);
        }

        fputc('\n', fp);
    }

    if(fclose(fp) != 0){
        return ROM_STROM_ERR_IO;
    }

    return ROM_STROM_SUCCESS;
}


static int strom_monolis_add_entry(
    const int row,
    const int column,
    const double value,
    void* context_pointer)
{
    ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT* context =
        (ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT*)context_pointer;

    int row_window;
    int row_mode;
    int column_window;
    int column_mode;

    if(context == NULL ||
       context->monolis == NULL ||
       row < 0 ||
       column < 0 ||
       row >= context->total_reduced_dof ||
       column >= context->total_reduced_dof)
    {
        return 1;
    }

    row_window = row / context->block_size;
    row_mode = row % context->block_size;
    column_window = column / context->block_size;
    column_mode = column % context->block_size;

    monolis_add_scalar_to_sparse_matrix_R(
        context->monolis,
        row_window,
        column_window,
        row_mode,
        column_mode,
        value);

    return 0;
}


static int strom_monolis_set_rhs(
    const int row,
    const double value,
    void* context_pointer)
{
    ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT* context =
        (ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT*)context_pointer;

    if(context == NULL ||
       context->monolis == NULL ||
       row < 0 ||
       row >= context->total_reduced_dof)
    {
        return 1;
    }

    context->monolis->mat.R.B[row] = value;
    return 0;
}


int ROM_std_strom_monolis_initialize(
    ROM_STROM_MONOLIS_SOLVER* solver,
    const ROM_STROM_SYSTEM* system,
    MONOLIS_COM* monolis_com,
    const int max_iter,
    const double tolerance)
{
    int status;
    int block_size;

    if(solver == NULL ||
       system == NULL ||
       monolis_com == NULL ||
       max_iter <= 0 ||
       !isfinite(tolerance) ||
       tolerance <= 0.0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    strom_monolis_zero(solver);

    status = ROM_std_strom_check_uniform_rank(system, &block_size);
    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    status = ROM_std_strom_build_bcsr_graph(
        system->num_windows,
        &solver->graph_index,
        &solver->graph_item,
        &solver->num_graph_items);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    solver->num_windows = system->num_windows;
    solver->block_size = block_size;
    solver->total_reduced_dof = system->total_reduced_dof;
    solver->monolis_com = monolis_com;
    solver->max_iter = max_iter;
    solver->tolerance = tolerance;

    monolis_initialize(&solver->monolis);

    solver->solution = (double*)calloc(
        (size_t)solver->total_reduced_dof,
        sizeof(double));

    if(solver->solution == NULL){
        ROM_std_strom_monolis_finalize(solver);
        return ROM_STROM_ERR_ALLOCATION;
    }

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        &solver->monolis,
        solver->num_windows,
        solver->block_size,
        solver->graph_index,
        solver->graph_item);

    /* The current ST-ROM executable is serial. */
    solver->monolis.mat.N = solver->num_windows;
    solver->monolis.mat.NP = solver->num_windows;

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_monolis_assemble(
    ROM_STROM_MONOLIS_SOLVER* solver,
    const ROM_STROM_SYSTEM* system)
{
    ROM_STROM_MONOLIS_ASSEMBLY_CONTEXT context;
    int status;

    if(solver == NULL ||
       system == NULL ||
       solver->total_reduced_dof != system->total_reduced_dof ||
       solver->block_size <= 0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    monolis_clear_mat_value_rhs_R(&solver->monolis);

    memset(
        solver->solution,
        0,
        (size_t)solver->total_reduced_dof * sizeof(double));

    context.monolis = &solver->monolis;
    context.block_size = solver->block_size;
    context.total_reduced_dof = solver->total_reduced_dof;

    status = ROM_std_strom_assemble_all_at_once(
        system,
        strom_monolis_add_entry,
        strom_monolis_set_rhs,
        &context);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    return ROM_STROM_SUCCESS;
}


double ROM_std_strom_reduced_residual_norm(
    const ROM_STROM_SYSTEM* system)
{
    double residual2 = 0.0;
    double lhs2 = 0.0;
    double rhs2 = 0.0;

    if(system == NULL ||
       system->windows == NULL ||
       system->num_windows <= 0)
    {
        return -1.0;
    }

    for(int w = 0; w < system->num_windows; w++){
        const ROM_STROM_WINDOW* current = &system->windows[w];

        if(current->reduced_diag == NULL ||
           current->reduced_rhs == NULL ||
           current->reduced_coef == NULL)
        {
            return -1.0;
        }

        for(int i = 0; i < current->num_modes; i++){
            double lhs = 0.0;

            for(int j = 0; j < current->num_modes; j++){
                lhs += current->reduced_diag[i][j]
                    * current->reduced_coef[j];
            }

            if(w > 0){
                const ROM_STROM_WINDOW* previous = &system->windows[w - 1];

                if(current->reduced_coupling_prev == NULL){
                    return -1.0;
                }

                for(int j = 0; j < previous->num_modes; j++){
                    lhs -= current->reduced_coupling_prev[i][j]
                        * previous->reduced_coef[j];
                }
            }

            {
                const double rhs = current->reduced_rhs[i];
                const double residual = lhs - rhs;

                residual2 += residual * residual;
                lhs2 += lhs * lhs;
                rhs2 += rhs * rhs;
            }
        }
    }

    return sqrt(residual2) /
        fmax(sqrt(lhs2) + sqrt(rhs2), 1.0e-30);
}


int ROM_std_strom_solve_all_at_once_monolis(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_MONOLIS_SOLVER* solver)
{
    int status;

    if(system == NULL ||
       solver == NULL ||
       solver->monolis_com == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    status = ROM_std_strom_monolis_assemble(solver, system);
    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    BBFE_sys_monowrap_solve(
        &solver->monolis,
        solver->monolis_com,
        solver->solution,
        MONOLIS_ITER_BICGSTAB,
        MONOLIS_PREC_DIAG,
        solver->max_iter,
        solver->tolerance);

    for(int w = 0; w < system->num_windows; w++){
        ROM_STROM_WINDOW* window = &system->windows[w];
        const int offset = w * solver->block_size;

        for(int i = 0; i < solver->block_size; i++){
            window->reduced_coef[i] = solver->solution[offset + i];
        }
    }

    return ROM_STROM_SUCCESS;
}


void ROM_std_strom_monolis_finalize(
    ROM_STROM_MONOLIS_SOLVER* solver)
{
    if(solver == NULL){
        return;
    }

    monolis_finalize(&solver->monolis);

    free(solver->graph_index);
    free(solver->graph_item);
    free(solver->solution);

    strom_monolis_zero(solver);
}
