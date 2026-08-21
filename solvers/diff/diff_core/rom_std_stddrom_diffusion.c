#include "rom_std_stddrom_diffusion.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int diffusion_stdd_validate(
    const ROM_STDD_DIFFUSION_CONTEXT* context)
{
    size_t expected_length;

    if(context == NULL || context->fom == NULL || context->layout == NULL ||
       context->reduced == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(context->layout->physical_dof != 1 ||
       context->layout->time.degree != 0 ||
       context->layout->time.n_dof != 1 ||
       context->layout->num_slabs <= 0)
    {
        return ROM_STDD_ERR_INCOMPATIBLE;
    }

    expected_length = ST_fom_window_vector_length(context->layout);

    if(context->reduced->num_space_nodes != context->fom->fe.total_num_nodes ||
       context->reduced->num_internal_space_nodes !=
           context->layout->num_owned_space_points ||
       context->reduced->st_dof_per_space_node !=
           context->layout->num_slabs ||
       (size_t)context->reduced->local_st_dof != expected_length)
    {
        return ROM_STDD_ERR_DIMENSION;
    }

    return ROM_STDD_SUCCESS;
}

static size_t diffusion_stdd_index(
    const ROM_STDD_DIFFUSION_CONTEXT* context,
    int node,
    int slab)
{
    return (size_t)node * (size_t)context->layout->num_slabs
        + (size_t)slab;
}

static int diffusion_stdd_assemble_source(
    ROM_STDD_DIFFUSION_CONTEXT* context,
    int window_id,
    double* source)
{
    BBFE_DATA* fe;
    BBFE_BASIS* basis;
    VALUES* vals;
    int nl;
    int np;
    int ns;
    double* value_ip = NULL;
    double* jacobian_ip = NULL;
    double** local_x = NULL;
    double** x_ip = NULL;
    double** v_ip = NULL;
    double* a_ip = NULL;
    double* k_ip = NULL;
    double* f_ip = NULL;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       source == NULL || window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    fe = &context->fom->fe;
    basis = &context->fom->basis;
    vals = &context->fom->vals;
    nl = fe->local_num_nodes;
    np = basis->num_integ_points;
    ns = context->layout->num_slabs;

    memset(
        source,
        0,
        (size_t)context->reduced->local_st_dof * sizeof(double));

    value_ip = (double*)calloc((size_t)np, sizeof(double));
    jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    local_x = (double**)calloc((size_t)nl, sizeof(double*));
    x_ip = (double**)calloc((size_t)np, sizeof(double*));
    v_ip = (double**)calloc((size_t)np, sizeof(double*));
    a_ip = (double*)calloc((size_t)np, sizeof(double));
    k_ip = (double*)calloc((size_t)np, sizeof(double));
    f_ip = (double*)calloc((size_t)np, sizeof(double));

    if(value_ip == NULL || jacobian_ip == NULL || local_x == NULL ||
       x_ip == NULL || v_ip == NULL || a_ip == NULL || k_ip == NULL ||
       f_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL){
            goto allocation_error;
        }
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        v_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL || v_ip[p] == NULL){
            goto allocation_error;
        }
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            manusol_get_conv_vel(v_ip[p], x_ip[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int slab = 0; slab < ns; slab++){
            const int global_step = window_id * ns + slab + 1;
            const double time = (double)global_step * vals->dt;

            for(int p = 0; p < np; p++){
                f_ip[p] = manusol_get_source(
                    x_ip[p],
                    time,
                    a_ip[p],
                    v_ip[p],
                    k_ip[p]);
            }

            for(int i = 0; i < nl; i++){
                for(int p = 0; p < np; p++){
                    value_ip[p] = BBFE_elemmat_convdiff_vec_source(
                        basis->N[p][i],
                        f_ip[p]);
                }

                source[diffusion_stdd_index(
                    context,
                    fe->conn[e][i],
                    slab)] += BBFE_std_integ_calc(
                        np,
                        value_ip,
                        basis->integ_weight,
                        jacobian_ip);
            }
        }
    }

    for(int p = 0; p < np; p++){
        free(x_ip[p]);
        free(v_ip[p]);
    }
    for(int i = 0; i < nl; i++){
        free(local_x[i]);
    }
    free(value_ip);
    free(jacobian_ip);
    free(local_x);
    free(x_ip);
    free(v_ip);
    free(a_ip);
    free(k_ip);
    free(f_ip);

    return ROM_STDD_SUCCESS;

allocation_error:
    if(x_ip != NULL){
        for(int p = 0; p < np; p++){
            free(x_ip[p]);
        }
    }
    if(v_ip != NULL){
        for(int p = 0; p < np; p++){
            free(v_ip[p]);
        }
    }
    if(local_x != NULL){
        for(int i = 0; i < nl; i++){
            free(local_x[i]);
        }
    }
    free(value_ip);
    free(jacobian_ip);
    free(local_x);
    free(x_ip);
    free(v_ip);
    free(a_ip);
    free(k_ip);
    free(f_ip);
    return ROM_STDD_ERR_ALLOCATION;
}

int ROM_std_stdd_diffusion_context_initialize(
    ROM_STDD_DIFFUSION_CONTEXT* context,
    FE_SYSTEM* fom,
    const ST_FOM_WINDOW_LAYOUT* layout,
    ROM_STDD_SYSTEM* reduced)
{
    size_t length;

    if(context == NULL || fom == NULL || layout == NULL || reduced == NULL ||
       layout->physical_dof != 1 || layout->time.degree != 0 ||
       layout->time.n_dof != 1 || layout->num_slabs <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(context, 0, sizeof(*context));
    context->fom = fom;
    context->layout = layout;
    context->reduced = reduced;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS){
        memset(context, 0, sizeof(*context));
        return ROM_STDD_ERR_DIMENSION;
    }

    length = (size_t)reduced->local_st_dof;

    context->operator_input = (double*)calloc(length, sizeof(double));
    context->lift = (double*)calloc(length, sizeof(double));
    context->previous_input = (double*)calloc(length, sizeof(double));
    context->source = (double*)calloc(length, sizeof(double));
    context->work_A = (double*)calloc(length, sizeof(double));
    context->work_B = (double*)calloc(length, sizeof(double));
    context->effective_rhs = (double*)calloc(length, sizeof(double));

    if(context->operator_input == NULL || context->lift == NULL ||
       context->previous_input == NULL || context->source == NULL ||
       context->work_A == NULL || context->work_B == NULL ||
       context->effective_rhs == NULL)
    {
        ROM_std_stdd_diffusion_context_finalize(context);
        return ROM_STDD_ERR_ALLOCATION;
    }

    monolis_com_initialize_by_self(&context->self_com);

    monolis_initialize(&context->window_operator);
    monolis_get_nonzero_pattern_by_simple_mesh_R(
        &context->window_operator,
        fom->fe.total_num_nodes,
        fom->fe.local_num_nodes,
        layout->num_slabs,
        fom->fe.total_num_elems,
        fom->fe.conn);

    monolis_initialize(&context->temporal_operator);
    monolis_get_nonzero_pattern_by_simple_mesh_R(
        &context->temporal_operator,
        fom->fe.total_num_nodes,
        fom->fe.local_num_nodes,
        layout->num_slabs,
        fom->fe.total_num_elems,
        fom->fe.conn);

    context->operators_initialized = 1;

    return ROM_std_stdd_diffusion_build_operators(context);
}

void ROM_std_stdd_diffusion_context_finalize(
    ROM_STDD_DIFFUSION_CONTEXT* context)
{
    if(context == NULL){
        return;
    }

    if(context->operators_initialized){
        monolis_finalize(&context->window_operator);
        monolis_finalize(&context->temporal_operator);
    }

    free(context->operator_input);
    free(context->lift);
    free(context->previous_input);
    free(context->source);
    free(context->work_A);
    free(context->work_B);
    free(context->effective_rhs);

    memset(context, 0, sizeof(*context));
}

int ROM_std_stdd_diffusion_build_operators(
    ROM_STDD_DIFFUSION_CONTEXT* context)
{
    BBFE_DATA* fe;
    BBFE_BASIS* basis;
    VALUES* vals;
    int nl;
    int np;
    int ns;
    double* mass_ip = NULL;
    double* diffusion_ip = NULL;
    double* jacobian_ip = NULL;
    double** local_x = NULL;
    double** x_ip = NULL;
    double* a_ip = NULL;
    double* k_ip = NULL;
    int status = ROM_STDD_SUCCESS;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       !context->operators_initialized)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    fe = &context->fom->fe;
    basis = &context->fom->basis;
    vals = &context->fom->vals;
    nl = fe->local_num_nodes;
    np = basis->num_integ_points;
    ns = context->layout->num_slabs;

    monolis_clear_mat_value_R(&context->window_operator);
    monolis_clear_mat_value_R(&context->temporal_operator);

    mass_ip = (double*)calloc((size_t)np, sizeof(double));
    diffusion_ip = (double*)calloc((size_t)np, sizeof(double));
    jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    local_x = (double**)calloc((size_t)nl, sizeof(double*));
    x_ip = (double**)calloc((size_t)np, sizeof(double*));
    a_ip = (double*)calloc((size_t)np, sizeof(double));
    k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(mass_ip == NULL || diffusion_ip == NULL || jacobian_ip == NULL ||
       local_x == NULL || x_ip == NULL || a_ip == NULL || k_ip == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL){
            status = ROM_STDD_ERR_ALLOCATION;
            goto cleanup;
        }
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL){
            status = ROM_STDD_ERR_ALLOCATION;
            goto cleanup;
        }
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int row_node = fe->conn[e][i];

            for(int j = 0; j < nl; j++){
                const int col_node = fe->conn[e][j];
                double mass_value;
                double diffusion_value;
                double diagonal_value;
                double coupling_value;

                for(int p = 0; p < np; p++){
                    mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                        basis->N[p][i],
                        basis->N[p][j],
                        a_ip[p]);

                    diffusion_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        k_ip[p]);
                }

                mass_value = BBFE_std_integ_calc(
                    np,
                    mass_ip,
                    basis->integ_weight,
                    jacobian_ip);
                diffusion_value = BBFE_std_integ_calc(
                    np,
                    diffusion_ip,
                    basis->integ_weight,
                    jacobian_ip);

                coupling_value = mass_value / vals->dt;
                diagonal_value = coupling_value + diffusion_value;

                for(int slab = 0; slab < ns; slab++){
                    monolis_add_scalar_to_sparse_matrix_R(
                        &context->window_operator,
                        row_node,
                        col_node,
                        slab,
                        slab,
                        diagonal_value);

                    if(slab > 0){
                        monolis_add_scalar_to_sparse_matrix_R(
                            &context->window_operator,
                            row_node,
                            col_node,
                            slab,
                            slab - 1,
                            -coupling_value);
                    }
                }

                monolis_add_scalar_to_sparse_matrix_R(
                    &context->temporal_operator,
                    row_node,
                    col_node,
                    0,
                    ns - 1,
                    coupling_value);
            }
        }
    }

cleanup:
    if(x_ip != NULL){
        for(int p = 0; p < np; p++){
            free(x_ip[p]);
        }
    }
    if(local_x != NULL){
        for(int i = 0; i < nl; i++){
            free(local_x[i]);
        }
    }
    free(mass_ip);
    free(diffusion_ip);
    free(jacobian_ip);
    free(local_x);
    free(x_ip);
    free(a_ip);
    free(k_ip);

    return status;
}

int ROM_std_stdd_diffusion_apply_A(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_DIFFUSION_CONTEXT* context =
        (ROM_STDD_DIFFUSION_CONTEXT*)user_context;
    size_t length;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       x == NULL || y == NULL || context->operator_input == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    length = (size_t)context->reduced->local_st_dof;
    memcpy(context->operator_input, x, length * sizeof(double));

    monolis_mpi_update_R(
        &context->fom->monolis_com,
        context->fom->fe.total_num_nodes,
        context->layout->num_slabs,
        context->operator_input);

    memset(y, 0, length * sizeof(double));
    monolis_matvec_product_R(
        &context->window_operator,
        &context->self_com,
        context->operator_input,
        y);

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_diffusion_apply_B(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_DIFFUSION_CONTEXT* context =
        (ROM_STDD_DIFFUSION_CONTEXT*)user_context;
    size_t length;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       x == NULL || y == NULL || context->operator_input == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    length = (size_t)context->reduced->local_st_dof;
    memcpy(context->operator_input, x, length * sizeof(double));

    monolis_mpi_update_R(
        &context->fom->monolis_com,
        context->fom->fe.total_num_nodes,
        context->layout->num_slabs,
        context->operator_input);

    memset(y, 0, length * sizeof(double));
    monolis_matvec_product_R(
        &context->temporal_operator,
        &context->self_com,
        context->operator_input,
        y);

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_diffusion_enforce_homogeneous_basis(
    ROM_STDD_DIFFUSION_CONTEXT* context)
{
    ROM_STDD_SYSTEM* reduced;
    const double minimum_norm = 128.0 * DBL_EPSILON;
    int ns;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       context->reduced->basis == NULL || context->reduced->num_modes <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    reduced = context->reduced;
    ns = context->layout->num_slabs;

    for(int node = 0; node < reduced->num_internal_space_nodes; node++){
        if(!context->fom->bc.D_bc_exists[node]){
            continue;
        }

        for(int slab = 0; slab < ns; slab++){
            const int row = node * ns + slab;
            for(int mode = 0; mode < reduced->num_modes; mode++){
                reduced->basis[
                    (size_t)row * (size_t)reduced->num_modes_capacity
                    + (size_t)mode] = 0.0;
            }
        }
    }

    /* Re-orthonormalize locally after exact boundary-row zeroing. */
    for(int mode = 0; mode < reduced->num_modes; mode++){
        for(int pass = 0; pass < 2; pass++){
            for(int previous = 0; previous < mode; previous++){
                double dot = 0.0;

                for(int row = 0; row < reduced->internal_st_dof; row++){
                    dot += reduced->basis[
                        (size_t)row * (size_t)reduced->num_modes_capacity
                        + (size_t)previous]
                        * reduced->basis[
                            (size_t)row * (size_t)reduced->num_modes_capacity
                            + (size_t)mode];
                }

                for(int row = 0; row < reduced->internal_st_dof; row++){
                    reduced->basis[
                        (size_t)row * (size_t)reduced->num_modes_capacity
                        + (size_t)mode] -= dot * reduced->basis[
                            (size_t)row * (size_t)reduced->num_modes_capacity
                            + (size_t)previous];
                }
            }
        }

        {
            double norm2 = 0.0;
            for(int row = 0; row < reduced->internal_st_dof; row++){
                const double value = reduced->basis[
                    (size_t)row * (size_t)reduced->num_modes_capacity
                    + (size_t)mode];
                norm2 += value * value;
            }

            if(!isfinite(norm2) ||
               norm2 <= minimum_norm * minimum_norm)
            {
                fprintf(
                    stderr,
                    "ERROR: diffusion ST-DDROM basis mode %d vanished after "
                    "Dirichlet homogenization on rank %d\n",
                    mode,
                    monolis_mpi_get_global_my_rank());
                return ROM_STDD_ERR_DIMENSION;
            }

            {
                const double inverse_norm = 1.0 / sqrt(norm2);
                for(int row = 0; row < reduced->internal_st_dof; row++){
                    reduced->basis[
                        (size_t)row * (size_t)reduced->num_modes_capacity
                        + (size_t)mode] *= inverse_norm;
                }
            }
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_diffusion_build_lift(
    const ROM_STDD_DIFFUSION_CONTEXT* context,
    int window_id,
    double* lift_out)
{
    int ns;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       lift_out == NULL || window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    ns = context->layout->num_slabs;
    memset(
        lift_out,
        0,
        (size_t)context->reduced->local_st_dof * sizeof(double));

    for(int node = 0; node < context->fom->fe.total_num_nodes; node++){
        if(!context->fom->bc.D_bc_exists[node]){
            continue;
        }

        for(int slab = 0; slab < ns; slab++){
            const int global_step = window_id * ns + slab + 1;
            const double time =
                (double)global_step * context->fom->vals.dt;

            lift_out[diffusion_stdd_index(context, node, slab)] =
                manusol_get_sol(
                    context->fom->fe.x[node][0],
                    context->fom->fe.x[node][1],
                    context->fom->fe.x[node][2],
                    time);
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_diffusion_project_all_rhs(
    ROM_STDD_DIFFUSION_CONTEXT* context)
{
    int ns;
    size_t length;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       context->reduced->basis == NULL || context->reduced->num_modes <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    ns = context->layout->num_slabs;
    length = (size_t)context->reduced->local_st_dof;

    for(int window = 0; window < context->reduced->num_windows; window++){
        int status;

        status = diffusion_stdd_assemble_source(
            context,
            window,
            context->source);
        if(status != ROM_STDD_SUCCESS){
            return status;
        }

        status = ROM_std_stdd_diffusion_build_lift(
            context,
            window,
            context->lift);
        if(status != ROM_STDD_SUCCESS){
            return status;
        }

        status = ROM_std_stdd_diffusion_apply_A(
            context->lift,
            context->work_A,
            context);
        if(status != ROM_STDD_SUCCESS){
            return status;
        }

        memset(context->previous_input, 0, length * sizeof(double));

        if(window == 0){
            /* B only reads the previous-window last slab. */
            for(int node = 0; node < context->fom->fe.total_num_nodes; node++){
                context->previous_input[
                    diffusion_stdd_index(context, node, ns - 1)] =
                    manusol_get_sol(
                        context->fom->fe.x[node][0],
                        context->fom->fe.x[node][1],
                        context->fom->fe.x[node][2],
                        0.0);
            }
        }
        else{
            const double time_start =
                (double)(window * ns) * context->fom->vals.dt;

            for(int node = 0; node < context->fom->fe.total_num_nodes; node++){
                if(context->fom->bc.D_bc_exists[node]){
                    context->previous_input[
                        diffusion_stdd_index(context, node, ns - 1)] =
                        manusol_get_sol(
                            context->fom->fe.x[node][0],
                            context->fom->fe.x[node][1],
                            context->fom->fe.x[node][2],
                            time_start);
                }
            }
        }

        status = ROM_std_stdd_diffusion_apply_B(
            context->previous_input,
            context->work_B,
            context);
        if(status != ROM_STDD_SUCCESS){
            return status;
        }

        for(size_t row = 0; row < length; row++){
            context->effective_rhs[row] =
                context->source[row]
                - context->work_A[row]
                + context->work_B[row];
        }

        status = ROM_std_stdd_project_local_rhs(
            context->reduced,
            window,
            context->effective_rhs,
            context->reduced->reduced_rhs
                + (size_t)window
                    * (size_t)context->reduced->num_modes_capacity);
        if(status != ROM_STDD_SUCCESS){
            return status;
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_diffusion_reconstruct_window(
    ROM_STDD_DIFFUSION_CONTEXT* context,
    int window_id,
    double* full_window_out)
{
    size_t length;
    int status;

    if(diffusion_stdd_validate(context) != ROM_STDD_SUCCESS ||
       full_window_out == NULL || window_id < 0 ||
       window_id >= context->reduced->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    length = (size_t)context->reduced->local_st_dof;

    status = ROM_std_stdd_reconstruct_local_homogeneous(
        context->reduced,
        window_id,
        full_window_out);
    if(status != ROM_STDD_SUCCESS){
        return status;
    }

    monolis_mpi_update_R(
        &context->fom->monolis_com,
        context->fom->fe.total_num_nodes,
        context->layout->num_slabs,
        full_window_out);

    status = ROM_std_stdd_diffusion_build_lift(
        context,
        window_id,
        context->lift);
    if(status != ROM_STDD_SUCCESS){
        return status;
    }

    for(size_t row = 0; row < length; row++){
        full_window_out[row] += context->lift[row];
    }

    return ROM_STDD_SUCCESS;
}
