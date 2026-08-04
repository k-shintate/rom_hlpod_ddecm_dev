#include "rom_std_strom_fom_adapter.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ROM_STROM_SNAPSHOT_MAGIC   UINT32_C(0x5354534E) /* STSN */
#define ROM_STROM_SNAPSHOT_VERSION UINT32_C(2)
#define ROM_STROM_REDUCED_OPERATOR_MAGIC UINT32_C(0x53544F50) /* STOP */
#define ROM_STROM_REDUCED_OPERATOR_VERSION UINT32_C(1)

static double strom_window_start_time(
    const ROM_STROM_FOM_CONTEXT* context,
    int window_id)
{
    return (double)(window_id * context->num_window_slabs)
        * context->sys->vals.dt;
}

static int strom_validate_fom_context(
    const ROM_STROM_FOM_CONTEXT* context)
{
    if(context == NULL ||
       context->sys == NULL ||
       context->te == NULL ||
       context->num_windows <= 0 ||
       context->num_window_slabs <= 0 ||
       context->num_space_nodes <= 0 ||
       context->num_time_dofs <= 0 ||
       context->previous_state_dof <= 0 ||
       context->st_dof <= 0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_STROM_SUCCESS;
}


static double** strom_alloc_rows(int rows, int columns)
{
    double** matrix = NULL;
    double* data = NULL;

    if(rows <= 0 || columns <= 0){
        return NULL;
    }

    matrix = (double**)malloc((size_t)rows * sizeof(double*));
    data = (double*)calloc(
        (size_t)rows * (size_t)columns,
        sizeof(double));

    if(matrix == NULL || data == NULL){
        free(matrix);
        free(data);
        return NULL;
    }

    for(int i = 0; i < rows; i++){
        matrix[i] = data + (size_t)i * (size_t)columns;
    }
    return matrix;
}

static void strom_free_rows(double*** matrix)
{
    if(matrix == NULL || *matrix == NULL){
        return;
    }
    free((*matrix)[0]);
    free(*matrix);
    *matrix = NULL;
}

static int strom_precompute_element_mass(
    ROM_STROM_FOM_CONTEXT* context)
{
    BBFE_DATA* fe;
    BBFE_BASIS* basis;
    const int nl = context->sys->fe.local_num_nodes;
    const int np = context->sys->basis.num_integ_points;
    const int ne = context->sys->fe.total_num_elems;
    double* val_ip = NULL;
    double* jacobian_ip = NULL;
    double** local_x = NULL;
    double** x_ip = NULL;
    double* a_ip = NULL;

    if(nl <= 0 || np <= 0 || ne <= 0){
        return ROM_STROM_ERR_DIMENSION;
    }

    context->local_num_nodes = nl;
    context->element_mass = (double*)calloc(
        (size_t)ne * (size_t)nl * (size_t)nl,
        sizeof(double));

    val_ip = (double*)calloc((size_t)np, sizeof(double));
    jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    local_x = strom_alloc_rows(nl, 3);
    x_ip = strom_alloc_rows(np, 3);
    a_ip = (double*)calloc((size_t)np, sizeof(double));

    if(context->element_mass == NULL ||
       val_ip == NULL ||
       jacobian_ip == NULL ||
       local_x == NULL ||
       x_ip == NULL ||
       a_ip == NULL)
    {
        free(context->element_mass);
        context->element_mass = NULL;
        free(val_ip);
        free(jacobian_ip);
        strom_free_rows(&local_x);
        strom_free_rows(&x_ip);
        free(a_ip);
        return ROM_STROM_ERR_ALLOCATION;
    }

    fe = &context->sys->fe;
    basis = &context->sys->basis;

    for(int e = 0; e < ne; e++){
        BBFE_elemmat_set_Jacobian_array(
            jacobian_ip,
            np,
            e,
            fe);

        BBFE_elemmat_set_local_array_vector(
            local_x,
            fe,
            fe->x,
            e,
            3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(
                x_ip[p],
                nl,
                local_x,
                basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            for(int j = 0; j < nl; j++){
                for(int p = 0; p < np; p++){
                    val_ip[p] =
                        BBFE_elemmat_convdiff_mat_mass(
                            basis->N[p][i],
                            basis->N[p][j],
                            a_ip[p]);
                }

                context->element_mass[
                    ((size_t)e * (size_t)nl + (size_t)i)
                    * (size_t)nl + (size_t)j] =
                    BBFE_std_integ_calc(
                        np,
                        val_ip,
                        basis->integ_weight,
                        jacobian_ip);
            }
        }
    }

    free(val_ip);
    free(jacobian_ip);
    strom_free_rows(&local_x);
    strom_free_rows(&x_ip);
    free(a_ip);
    return ROM_STROM_SUCCESS;
}

static void strom_zero_rhs(
    ROM_STROM_FOM_CONTEXT* context)
{
    memset(
        context->sys->monolis.mat.R.B,
        0,
        (size_t)context->st_dof * sizeof(double));
}

static int strom_reset_matrix(
    ROM_STROM_FOM_CONTEXT* context)
{
    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS){
        return ROM_STROM_ERR_ARGUMENT;
    }

    /* Remove matrix/RHS modifications made for the preceding window. */
    monolis_clear_mat_value_rhs_R(&context->sys->monolis);

    /* Restore the unmodified constant window matrix. */
    monolis_copy_mat_value_R(
        &context->sys->monolis0,
        &context->sys->monolis);

    strom_zero_rhs(context);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_context_initialize(
    ROM_STROM_FOM_CONTEXT* context,
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_windows)
{
    long long previous_state_dof;
    long long st_dof;

    if(context == NULL ||
       sys == NULL ||
       te == NULL ||
       num_window_slabs <= 0 ||
       num_windows <= 0 ||
       sys->fe.total_num_nodes <= 0 ||
       te->n_dof <= 0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    previous_state_dof =
        (long long)sys->fe.total_num_nodes * te->n_dof;

    st_dof = previous_state_dof * num_window_slabs;

    if(previous_state_dof > INT_MAX || st_dof > INT_MAX){
        return ROM_STROM_ERR_DIMENSION;
    }

    memset(context, 0, sizeof(*context));

    context->sys = sys;
    context->te = te;
    context->num_windows = num_windows;
    context->num_window_slabs = num_window_slabs;
    context->num_space_nodes = sys->fe.total_num_nodes;
    context->num_time_dofs = te->n_dof;
    context->previous_state_dof = (int)previous_state_dof;
    context->st_dof = (int)st_dof;

    /*
     * The serial ST graph must contain every scalar space-time unknown.
     * Catch accidental reuse of the ordinary spatial matrix immediately.
     */
    if(sys->monolis0.mat.N != context->st_dof ||
       sys->monolis.mat.N != context->st_dof ||
       sys->monolis0.mat.NDOF != 1 ||
       sys->monolis.mat.NDOF != 1)
    {
        memset(context, 0, sizeof(*context));
        return ROM_STROM_ERR_DIMENSION;
    }

    context->zero_previous_state = (double*)calloc(
        (size_t)context->previous_state_dof,
        sizeof(double));

    context->previous_state = (double*)calloc(
        (size_t)context->previous_state_dof,
        sizeof(double));

    context->rhs_with_previous = (double*)calloc(
        (size_t)context->st_dof,
        sizeof(double));

    context->rhs_without_previous = (double*)calloc(
        (size_t)context->st_dof,
        sizeof(double));

    if(context->zero_previous_state == NULL ||
       context->previous_state == NULL ||
       context->rhs_with_previous == NULL ||
       context->rhs_without_previous == NULL)
    {
        ROM_std_strom_fom_context_finalize(context);
        return ROM_STROM_ERR_ALLOCATION;
    }

    /* Element mass matrices are built lazily only when direct B is used. */
    return ROM_STROM_SUCCESS;
}

void ROM_std_strom_fom_context_finalize(
    ROM_STROM_FOM_CONTEXT* context)
{
    if(context == NULL){
        return;
    }

    free(context->zero_previous_state);
    free(context->previous_state);
    free(context->rhs_with_previous);
    free(context->rhs_without_previous);
    free(context->element_mass);

    memset(context, 0, sizeof(*context));
}

int ROM_std_strom_fom_prepare_operator(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id)
{
    double t_window_start;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(strom_reset_matrix(context) != ROM_STROM_SUCCESS){
        return ROM_STROM_ERR_OPERATOR;
    }

    t_window_start = strom_window_start_time(context, window_id);

    /*
     * Strong Dirichlet treatment changes the matrix.  The actual boundary
     * values enter the affine RHS, but the resulting matrix is the A_w that
     * must be projected.
     */
    BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar(
        &context->sys->monolis,
        &context->sys->fe,
        &context->sys->vals,
        &context->sys->bc,
        context->te,
        context->num_window_slabs,
        t_window_start,
        context->sys->monolis.mat.R.B);

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_fom_prepare_A_callback(
    int window_id,
    void* user_context)
{
    return ROM_std_strom_fom_prepare_operator(
        (ROM_STROM_FOM_CONTEXT*)user_context,
        window_id);
}

int ROM_std_strom_fom_apply_prepared_A_callback(
    const double* x,
    double* y,
    void* user_context)
{
    return ROM_std_strom_fom_apply_prepared_A(
        (ROM_STROM_FOM_CONTEXT*)user_context,
        x,
        y);
}

int ROM_std_strom_fom_build_dirichlet_lift(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    double* lift)
{
    double t_window_start;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       lift == NULL ||
       window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    t_window_start = strom_window_start_time(context, window_id);

    memset(
        lift,
        0,
        (size_t)context->st_dof * sizeof(double));

    for(int slab = 0; slab < context->num_window_slabs; slab++){
        const double t_s =
            t_window_start
            + (double)slab * context->sys->vals.dt;

        for(int alpha = 0; alpha < context->num_time_dofs; alpha++){
            const double t_alpha =
                t_s
                + context->te->dof_tau[alpha]
                * context->sys->vals.dt;

            manusol_set_theo_sol(
                &context->sys->fe,
                context->sys->vals.theo_sol,
                t_alpha);

            BBFE_manusol_set_bc_scalar(
                &context->sys->fe,
                &context->sys->bc,
                context->sys->vals.theo_sol,
                t_alpha);

            for(int i = 0; i < context->num_space_nodes; i++){
                if(!context->sys->bc.D_bc_exists[i]){
                    continue;
                }

                lift[
                    ST_window_gid(
                        slab,
                        alpha,
                        i,
                        context->num_time_dofs,
                        context->num_space_nodes)] =
                    context->sys->bc.imposed_D_val[i];
            }
        }
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_apply_prepared_A(
    ROM_STROM_FOM_CONTEXT* context,
    const double* x,
    double* y)
{
    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       x == NULL ||
       y == NULL)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    memset(y, 0, (size_t)context->st_dof * sizeof(double));

    monolis_matvec_product_R(
        &context->sys->monolis,
        &context->sys->monolis_com,
        (double*)x,
        y);

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_fom_apply_A(
    int window_id,
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STROM_FOM_CONTEXT* context =
        (ROM_STROM_FOM_CONTEXT*)user_context;

    if(context == NULL || x == NULL || y == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(ROM_std_strom_fom_prepare_operator(context, window_id)
       != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_OPERATOR;
    }

    return ROM_std_strom_fom_apply_prepared_A(
        context,
        x,
        y);
}

int ROM_std_strom_fom_build_full_rhs(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_state,
    double* full_rhs)
{
    double t_window_start;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       previous_state == NULL ||
       full_rhs == NULL ||
       window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(strom_reset_matrix(context) != ROM_STROM_SUCCESS){
        return ROM_STROM_ERR_OPERATOR;
    }

    t_window_start = strom_window_start_time(context, window_id);

    set_element_vec_ST_dG_window(
        &context->sys->monolis,
        &context->sys->fe,
        &context->sys->basis,
        &context->sys->vals,
        context->te,
        previous_state,
        t_window_start,
        context->num_window_slabs);

    BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar(
        &context->sys->monolis,
        &context->sys->fe,
        &context->sys->vals,
        &context->sys->bc,
        context->te,
        context->num_window_slabs,
        t_window_start,
        context->sys->monolis.mat.R.B);

    memcpy(
        full_rhs,
        context->sys->monolis.mat.R.B,
        (size_t)context->st_dof * sizeof(double));

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_apply_B(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* user_context)
{
    ROM_STROM_FOM_CONTEXT* context =
        (ROM_STROM_FOM_CONTEXT*)user_context;
    int status;

    if(context == NULL ||
       x_prev == NULL ||
       y_cur == NULL ||
       current_window_id <= 0 ||
       current_window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    /*
     * C_last x_prev: extract all temporal DOFs of the previous window's
     * final slab.  This is exactly the state expected by the FOM inflow
     * assembly set_element_vec_ST_dG_window().
     */
    ST_Window_extract_last_slab_state(
        context->te,
        context->num_space_nodes,
        context->num_window_slabs,
        x_prev,
        context->previous_state);

    status = ROM_std_strom_fom_build_full_rhs(
        context,
        current_window_id,
        context->previous_state,
        context->rhs_with_previous);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    status = ROM_std_strom_fom_build_full_rhs(
        context,
        current_window_id,
        context->zero_previous_state,
        context->rhs_without_previous);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    /*
     * Source and nonhomogeneous Dirichlet terms cancel.  What remains is
     * the linear previous-window contribution B_w x_prev.
     */
    for(int i = 0; i < context->st_dof; i++){
        y_cur[i] =
            context->rhs_with_previous[i]
            - context->rhs_without_previous[i];
    }

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_fom_apply_B_direct(
    int current_window_id,
    const double* x_prev,
    double* y_cur,
    void* user_context)
{
    ROM_STROM_FOM_CONTEXT* context =
        (ROM_STROM_FOM_CONTEXT*)user_context;
    const int nt = context != NULL ? context->num_time_dofs : 0;
    const int nx = context != NULL ? context->num_space_nodes : 0;
    const int nl = context != NULL ? context->local_num_nodes : 0;
    const int ne = context != NULL
        ? context->sys->fe.total_num_elems
        : 0;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       x_prev == NULL ||
       y_cur == NULL ||
       current_window_id <= 0 ||
       current_window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(context->element_mass == NULL){
        const int status = strom_precompute_element_mass(context);
        if(status != ROM_STROM_SUCCESS){
            return status;
        }
    }

    ST_Window_extract_last_slab_state(
        context->te,
        nx,
        context->num_window_slabs,
        x_prev,
        context->previous_state);

    memset(y_cur, 0, (size_t)context->st_dof * sizeof(double));

    for(int e = 0; e < ne; e++){
        for(int alpha = 0; alpha < nt; alpha++){
            for(int beta = 0; beta < nt; beta++){
                const double coefficient =
                    context->te->E_prev[alpha][beta]
                    / context->sys->vals.dt;

                if(fabs(coefficient) < 1.0e-30){
                    continue;
                }

                for(int i = 0; i < nl; i++){
                    double value = 0.0;
                    const int row_node = context->sys->fe.conn[e][i];

                    for(int j = 0; j < nl; j++){
                        const int column_node = context->sys->fe.conn[e][j];
                        const double mass_ij = context->element_mass[
                            ((size_t)e * (size_t)nl + (size_t)i)
                            * (size_t)nl + (size_t)j];

                        value += mass_ij
                            * context->previous_state[
                                beta * nx + column_node];
                    }

                    y_cur[
                        ST_window_gid(
                            0,
                            alpha,
                            row_node,
                            nt,
                            nx)] += coefficient * value;
                }
            }
        }
    }

    /* Difference of two strongly constrained RHS vectors is zero on BC rows. */
    for(int alpha = 0; alpha < nt; alpha++){
        for(int i = 0; i < nx; i++){
            if(context->sys->bc.D_bc_exists[i]){
                y_cur[ST_window_gid(0, alpha, i, nt, nx)] = 0.0;
            }
        }
    }

    return ROM_STROM_SUCCESS;
}

static int strom_snapshot_filename(
    char* filename,
    size_t filename_size,
    const char* directory,
    int snapshot_id,
    int window_id)
{
    if(filename == NULL ||
       filename_size == 0 ||
       directory == NULL ||
       snapshot_id < 0 ||
       window_id < 0)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(snprintf(
        filename,
        filename_size,
        "%s/snapshot_%04d_window_%04d.bin",
        directory,
        snapshot_id,
        window_id) >= (int)filename_size)
    {
        return ROM_STROM_ERR_IO;
    }

    return ROM_STROM_SUCCESS;
}


int ROM_std_strom_fom_calculate_residual(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_state,
    const double* window_solution,
    double* full_rhs,
    double* matrix_times_solution,
    double* residual,
    double* residual_norm2,
    double* matrix_times_solution_norm2,
    double* full_rhs_norm2)
{
    double r2 = 0.0;
    double ax2 = 0.0;
    double rhs2 = 0.0;
    int status;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       previous_state == NULL ||
       window_solution == NULL ||
       full_rhs == NULL ||
       matrix_times_solution == NULL ||
       residual == NULL ||
       window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    /*
     * build_full_rhs() leaves context->sys->monolis containing the exact
     * boundary-modified matrix paired with full_rhs.  Apply that prepared
     * matrix directly instead of rebuilding it in a second BC pass.
     */
    status = ROM_std_strom_fom_build_full_rhs(
        context,
        window_id,
        previous_state,
        full_rhs);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    status = ROM_std_strom_fom_apply_prepared_A(
        context,
        window_solution,
        matrix_times_solution);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    for(int i = 0; i < context->st_dof; i++){
        residual[i] = matrix_times_solution[i] - full_rhs[i];

        r2 += residual[i] * residual[i];
        ax2 += matrix_times_solution[i] * matrix_times_solution[i];
        rhs2 += full_rhs[i] * full_rhs[i];
    }

    if(residual_norm2 != NULL){
        *residual_norm2 = r2;
    }
    if(matrix_times_solution_norm2 != NULL){
        *matrix_times_solution_norm2 = ax2;
    }
    if(full_rhs_norm2 != NULL){
        *full_rhs_norm2 = rhs2;
    }

    return ROM_STROM_SUCCESS;
}


static void strom_clear_residual_norms(
    ROM_STROM_RESIDUAL_NORMS* norms)
{
    if(norms == NULL){
        return;
    }

    norms->residual_norm2 = 0.0;
    norms->lhs_norm2 = 0.0;
    norms->rhs_norm2 = 0.0;
}

static void strom_add_residual_entry(
    ROM_STROM_RESIDUAL_NORMS* norms,
    double residual,
    double lhs,
    double rhs)
{
    if(norms == NULL){
        return;
    }

    norms->residual_norm2 += residual * residual;
    norms->lhs_norm2 += lhs * lhs;
    norms->rhs_norm2 += rhs * rhs;
}

double ROM_std_strom_relative_residual(
    const ROM_STROM_RESIDUAL_NORMS* norms)
{
    double scale;

    if(norms == NULL){
        return -1.0;
    }

    scale =
        sqrt(fmax(norms->lhs_norm2, 0.0))
        + sqrt(fmax(norms->rhs_norm2, 0.0));

    return
        sqrt(fmax(norms->residual_norm2, 0.0))
        / fmax(scale, 1.0e-30);
}

void ROM_std_strom_accumulate_residual_norms(
    ROM_STROM_RESIDUAL_NORMS* destination,
    const ROM_STROM_RESIDUAL_NORMS* source)
{
    if(destination == NULL || source == NULL){
        return;
    }

    destination->residual_norm2 += source->residual_norm2;
    destination->lhs_norm2 += source->lhs_norm2;
    destination->rhs_norm2 += source->rhs_norm2;
}

int ROM_std_strom_fom_calculate_lifted_residual(
    ROM_STROM_FOM_CONTEXT* context,
    int window_id,
    const double* previous_full_state,
    const double* homogeneous_solution,
    const double* current_lift,
    double* full_rhs,
    double* matrix_times_lift,
    double* homogeneous_rhs,
    double* matrix_times_homogeneous,
    double* homogeneous_residual,
    double* full_residual,
    ROM_STROM_LIFTED_RESIDUAL* result)
{
    int status;

    if(strom_validate_fom_context(context) != ROM_STROM_SUCCESS ||
       previous_full_state == NULL ||
       homogeneous_solution == NULL ||
       current_lift == NULL ||
       full_rhs == NULL ||
       matrix_times_lift == NULL ||
       homogeneous_rhs == NULL ||
       matrix_times_homogeneous == NULL ||
       homogeneous_residual == NULL ||
       full_residual == NULL ||
       result == NULL ||
       window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    memset(result, 0, sizeof(*result));
    strom_clear_residual_norms(&result->all);
    strom_clear_residual_norms(&result->free_dofs);
    strom_clear_residual_norms(&result->dirichlet_dofs);

    /*
     * This constructs a matched pair (A_tilde, b_tilde) and leaves the
     * boundary-modified A_tilde in context->sys->monolis.
     */
    status = ROM_std_strom_fom_build_full_rhs(
        context,
        window_id,
        previous_full_state,
        full_rhs);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    /* All following products deliberately use that same prepared A_tilde. */
    status = ROM_std_strom_fom_apply_prepared_A(
        context,
        current_lift,
        matrix_times_lift);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    status = ROM_std_strom_fom_apply_prepared_A(
        context,
        homogeneous_solution,
        matrix_times_homogeneous);

    if(status != ROM_STROM_SUCCESS){
        return status;
    }

    for(int gid = 0; gid < context->st_dof; gid++){
        const int space_node = gid % context->num_space_nodes;
        const int is_dirichlet =
            context->sys->bc.D_bc_exists[space_node] != 0;

        const double rhs_hom =
            full_rhs[gid] - matrix_times_lift[gid];

        const double residual_hom =
            matrix_times_homogeneous[gid] - rhs_hom;

        const double lhs_full =
            matrix_times_lift[gid]
            + matrix_times_homogeneous[gid];

        const double residual_full_value =
            lhs_full - full_rhs[gid];

        const double identity_difference =
            residual_full_value - residual_hom;

        homogeneous_rhs[gid] = rhs_hom;
        homogeneous_residual[gid] = residual_hom;
        full_residual[gid] = residual_full_value;

        strom_add_residual_entry(
            &result->all,
            residual_hom,
            matrix_times_homogeneous[gid],
            rhs_hom);

        if(is_dirichlet){
            strom_add_residual_entry(
                &result->dirichlet_dofs,
                residual_hom,
                matrix_times_homogeneous[gid],
                rhs_hom);
        }
        else{
            strom_add_residual_entry(
                &result->free_dofs,
                residual_hom,
                matrix_times_homogeneous[gid],
                rhs_hom);
        }

        result->full_hom_difference_norm2 +=
            identity_difference * identity_difference;

        result->full_hom_scale_norm2 +=
            lhs_full * lhs_full
            + full_rhs[gid] * full_rhs[gid]
            + matrix_times_homogeneous[gid]
                * matrix_times_homogeneous[gid]
            + rhs_hom * rhs_hom;
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_write_snapshot(
    const char* directory,
    int snapshot_id,
    int window_id,
    const double* values,
    int n)
{
    char filename[4096];
    uint32_t header[5];
    FILE* fp;

    if(values == NULL ||
       n <= 0 ||
       strom_snapshot_filename(
           filename,
           sizeof(filename),
           directory,
           snapshot_id,
           window_id) != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    fp = fopen(filename, "wb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    header[0] = ROM_STROM_SNAPSHOT_MAGIC;
    header[1] = ROM_STROM_SNAPSHOT_VERSION;
    header[2] = (uint32_t)snapshot_id;
    header[3] = (uint32_t)window_id;
    header[4] = (uint32_t)n;

    if(fwrite(header, sizeof(uint32_t), 5, fp) != 5 ||
       fwrite(values, sizeof(double), (size_t)n, fp) != (size_t)n)
    {
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    fclose(fp);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_read_snapshot(
    const char* directory,
    int snapshot_id,
    int window_id,
    double* values,
    int n)
{
    char filename[4096];
    uint32_t header[5];
    FILE* fp;

    if(values == NULL ||
       n <= 0 ||
       strom_snapshot_filename(
           filename,
           sizeof(filename),
           directory,
           snapshot_id,
           window_id) != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    if(fread(header, sizeof(uint32_t), 5, fp) != 5 ||
       header[0] != ROM_STROM_SNAPSHOT_MAGIC ||
       header[1] != ROM_STROM_SNAPSHOT_VERSION ||
       (int)header[2] != snapshot_id ||
       (int)header[3] != window_id ||
       (int)header[4] != n ||
       fread(values, sizeof(double), (size_t)n, fp) != (size_t)n)
    {
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    fclose(fp);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_write_basis(
    const char* directory,
    const ROM_STROM_WINDOW* window)
{
    char filename[4096];

    if(directory == NULL ||
       window == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/window_%04d_pod.bin",
           directory,
           window->window_id) >= (int)sizeof(filename))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_std_strom_write_pod_modes_binary(window, filename);
}

int ROM_std_strom_fom_read_basis(
    const char* directory,
    ROM_STROM_WINDOW* window)
{
    char filename[4096];

    if(directory == NULL ||
       window == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/window_%04d_pod.bin",
           directory,
           window->window_id) >= (int)sizeof(filename))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_std_strom_read_pod_modes_binary(window, filename);
}

int ROM_std_strom_fom_write_shared_basis(
    const char* directory,
    const ROM_STROM_WINDOW* basis_holder)
{
    char filename[4096];

    if(directory == NULL ||
       basis_holder == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/all_windows_shared_pod.bin",
           directory) >= (int)sizeof(filename))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_std_strom_write_pod_modes_binary(
        basis_holder,
        filename);
}

int ROM_std_strom_fom_read_shared_basis(
    const char* directory,
    ROM_STROM_WINDOW* window)
{
    char filename[4096];

    if(directory == NULL ||
       window == NULL ||
       snprintf(
           filename,
           sizeof(filename),
           "%s/all_windows_shared_pod.bin",
           directory) >= (int)sizeof(filename))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    return ROM_std_strom_read_pod_modes_binary(
        window,
        filename);
}


static int strom_reduced_operator_filename(
    char* filename,
    size_t filename_size,
    const char* directory)
{
    if(filename == NULL || filename_size == 0 || directory == NULL){
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(snprintf(
        filename,
        filename_size,
        "%s/all_windows_shared_reduced_operators.bin",
        directory) >= (int)filename_size)
    {
        return ROM_STROM_ERR_IO;
    }

    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_write_reduced_operators(
    const char* directory,
    const ROM_STROM_SYSTEM* system,
    int operator_reuse,
    double dt,
    int st_dof)
{
    char filename[4096];
    uint32_t header[8];
    FILE* fp = NULL;
    int num_modes;
    int num_diag_blocks;
    int num_coupling_blocks;

    if(directory == NULL ||
       system == NULL ||
       system->windows == NULL ||
       system->num_windows <= 0 ||
       st_dof <= 0 ||
       dt <= 0.0 ||
       (operator_reuse != 0 && operator_reuse != 1))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    num_modes = system->windows[0].num_modes;
    if(num_modes <= 0){
        return ROM_STROM_ERR_DIMENSION;
    }

    num_diag_blocks = operator_reuse ? 1 : system->num_windows;
    num_coupling_blocks = system->num_windows > 1
        ? (operator_reuse ? 1 : system->num_windows - 1)
        : 0;

    if(strom_reduced_operator_filename(
        filename,
        sizeof(filename),
        directory) != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_IO;
    }

    fp = fopen(filename, "wb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    header[0] = ROM_STROM_REDUCED_OPERATOR_MAGIC;
    header[1] = ROM_STROM_REDUCED_OPERATOR_VERSION;
    header[2] = (uint32_t)system->num_windows;
    header[3] = (uint32_t)num_modes;
    header[4] = (uint32_t)operator_reuse;
    header[5] = (uint32_t)num_diag_blocks;
    header[6] = (uint32_t)num_coupling_blocks;
    header[7] = (uint32_t)st_dof;

    if(fwrite(header, sizeof(uint32_t), 8, fp) != 8 ||
       fwrite(&dt, sizeof(double), 1, fp) != 1)
    {
        fclose(fp);
        return ROM_STROM_ERR_IO;
    }

    for(int block = 0; block < num_diag_blocks; block++){
        const int w = operator_reuse ? 0 : block;
        const ROM_STROM_WINDOW* window = &system->windows[w];

        if(window->num_modes != num_modes || window->reduced_diag == NULL ||
           fwrite(
               window->reduced_diag[0],
               sizeof(double),
               (size_t)num_modes * (size_t)num_modes,
               fp) != (size_t)num_modes * (size_t)num_modes)
        {
            fclose(fp);
            return ROM_STROM_ERR_IO;
        }
    }

    for(int block = 0; block < num_coupling_blocks; block++){
        const int w = operator_reuse ? 1 : block + 1;
        const ROM_STROM_WINDOW* window = &system->windows[w];

        if(window->num_modes != num_modes ||
           window->reduced_coupling_prev == NULL ||
           fwrite(
               window->reduced_coupling_prev[0],
               sizeof(double),
               (size_t)num_modes * (size_t)num_modes,
               fp) != (size_t)num_modes * (size_t)num_modes)
        {
            fclose(fp);
            return ROM_STROM_ERR_IO;
        }
    }

    fclose(fp);
    return ROM_STROM_SUCCESS;
}

int ROM_std_strom_fom_read_reduced_operators(
    const char* directory,
    ROM_STROM_SYSTEM* system,
    int operator_reuse,
    double dt,
    int st_dof)
{
    char filename[4096];
    uint32_t header[8];
    FILE* fp = NULL;
    double file_dt = 0.0;
    double* block_data = NULL;
    int num_modes;
    int num_diag_blocks;
    int num_coupling_blocks;
    int status = ROM_STROM_SUCCESS;

    if(directory == NULL ||
       system == NULL ||
       system->windows == NULL ||
       system->num_windows <= 0 ||
       st_dof <= 0 ||
       dt <= 0.0 ||
       (operator_reuse != 0 && operator_reuse != 1))
    {
        return ROM_STROM_ERR_ARGUMENT;
    }

    if(strom_reduced_operator_filename(
        filename,
        sizeof(filename),
        directory) != ROM_STROM_SUCCESS)
    {
        return ROM_STROM_ERR_IO;
    }

    fp = fopen(filename, "rb");
    if(fp == NULL){
        return ROM_STROM_ERR_IO;
    }

    if(fread(header, sizeof(uint32_t), 8, fp) != 8 ||
       fread(&file_dt, sizeof(double), 1, fp) != 1 ||
       header[0] != ROM_STROM_REDUCED_OPERATOR_MAGIC ||
       header[1] != ROM_STROM_REDUCED_OPERATOR_VERSION ||
       (int)header[2] != system->num_windows ||
       (int)header[3] != system->windows[0].num_modes ||
       (int)header[4] != operator_reuse ||
       (int)header[7] != st_dof ||
       fabs(file_dt - dt) > 1.0e-12 * fmax(1.0, fabs(dt)))
    {
        fclose(fp);
        return ROM_STROM_ERR_DIMENSION;
    }

    num_modes = (int)header[3];
    num_diag_blocks = (int)header[5];
    num_coupling_blocks = (int)header[6];

    if(num_diag_blocks != (operator_reuse ? 1 : system->num_windows) ||
       num_coupling_blocks != (system->num_windows > 1
            ? (operator_reuse ? 1 : system->num_windows - 1)
            : 0))
    {
        fclose(fp);
        return ROM_STROM_ERR_DIMENSION;
    }

    block_data = (double*)calloc(
        (size_t)num_modes * (size_t)num_modes,
        sizeof(double));
    if(block_data == NULL){
        fclose(fp);
        return ROM_STROM_ERR_ALLOCATION;
    }

    for(int block = 0; block < num_diag_blocks; block++){
        const int w = operator_reuse ? 0 : block;

        if(fread(
            block_data,
            sizeof(double),
            (size_t)num_modes * (size_t)num_modes,
            fp) != (size_t)num_modes * (size_t)num_modes)
        {
            status = ROM_STROM_ERR_IO;
            goto cleanup;
        }

        status = ROM_std_strom_set_reduced_diag_data(
            &system->windows[w],
            block_data);
        if(status != ROM_STROM_SUCCESS){
            goto cleanup;
        }
    }

    if(operator_reuse){
        for(int w = 1; w < system->num_windows; w++){
            status = ROM_std_strom_copy_reduced_diag(
                &system->windows[w],
                &system->windows[0]);
            if(status != ROM_STROM_SUCCESS){
                goto cleanup;
            }
        }
    }

    for(int block = 0; block < num_coupling_blocks; block++){
        const int w = operator_reuse ? 1 : block + 1;

        if(fread(
            block_data,
            sizeof(double),
            (size_t)num_modes * (size_t)num_modes,
            fp) != (size_t)num_modes * (size_t)num_modes)
        {
            status = ROM_STROM_ERR_IO;
            goto cleanup;
        }

        status = ROM_std_strom_set_reduced_coupling_data(
            &system->windows[w],
            block_data,
            num_modes);
        if(status != ROM_STROM_SUCCESS){
            goto cleanup;
        }
    }

    if(operator_reuse && system->num_windows > 1){
        for(int w = 2; w < system->num_windows; w++){
            status = ROM_std_strom_copy_reduced_coupling(
                &system->windows[w],
                &system->windows[1]);
            if(status != ROM_STROM_SUCCESS){
                goto cleanup;
            }
        }
    }

cleanup:
    free(block_data);
    fclose(fp);
    return status;
}
