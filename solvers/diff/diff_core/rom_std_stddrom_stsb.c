#include "rom_std_stddrom_stsb.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int stdd_stsb_validate(const ROM_STDD_STSB_CONTEXT* context)
{
    if(context == NULL || context->sys == NULL ||
       context->time_element == NULL || context->num_window_slabs <= 0 ||
       context->num_windows <= 0 || context->st_dof_per_space_node <= 0 ||
       context->local_st_dof <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }
    return ROM_STDD_SUCCESS;
}

static double stdd_stsb_window_start_time(
    const ROM_STDD_STSB_CONTEXT* context,
    int window_id)
{
    return (double)(window_id * context->num_window_slabs)
        * context->sys->vals.dt;
}

/*
 * Build one matched strong-Dirichlet pair
 *
 *     A_tilde(window_id), b_tilde(window_id, previous_trace)
 *
 * using exactly the same sequence as STSB_solve_window():
 *
 *     1. restore the primitive constant window matrix,
 *     2. assemble the physical window RHS,
 *     3. apply the strong Dirichlet transformation to BOTH A and b.
 *
 * On success, context->sys->monolis contains A_tilde and its R.B array
 * contains b_tilde.  When full_rhs is non-NULL, b_tilde is copied out
 * without rebuilding A_tilde.  Every subsequent lift or homogeneous-state
 * product must use ROM_std_stdd_stsb_apply_prepared_A() immediately, so all
 * products use this exact same matrix instance.
 */
static int stdd_stsb_build_strong_pair(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* previous_trace,
    double* full_rhs)
{
    double window_start_time;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       previous_trace == NULL ||
       window_id < 0 || window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    monolis_clear_mat_value_rhs_R(&context->sys->monolis);
    monolis_copy_mat_value_R(
        &context->sys->monolis0,
        &context->sys->monolis);

    window_start_time = stdd_stsb_window_start_time(context, window_id);

    STSB_assemble_window_rhs(
        &context->sys->monolis,
        &context->sys->fe,
        &context->sys->basis,
        &context->sys->vals,
        context->time_element,
        context->num_window_slabs,
        window_start_time,
        previous_trace);

    STSB_apply_dirichlet_bc(
        &context->sys->monolis,
        &context->sys->fe,
        &context->sys->bc,
        &context->sys->vals,
        context->time_element,
        context->num_window_slabs,
        window_start_time,
        context->sys->monolis.mat.R.B);

    if(full_rhs != NULL){
        memcpy(
            full_rhs,
            context->sys->monolis.mat.R.B,
            (size_t)context->local_st_dof * sizeof(double));
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_context_initialize(
    ROM_STDD_STSB_CONTEXT* context,
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* time_element,
    int num_window_slabs,
    int num_windows)
{
    int ndof;
    size_t local_count;
    size_t trace_count;

    if(context == NULL || sys == NULL || time_element == NULL ||
       num_window_slabs <= 0 || num_windows <= 0 ||
       sys->fe.total_num_nodes <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(context, 0, sizeof(*context));

    ndof = STSB_get_num_dofs_per_node(
        time_element,
        num_window_slabs);

    context->sys = sys;
    context->time_element = time_element;
    context->num_window_slabs = num_window_slabs;
    context->num_windows = num_windows;
    context->st_dof_per_space_node = ndof;
    context->local_st_dof = sys->fe.total_num_nodes * ndof;

    /*
     * The verified HLPOD-DDROM implementation separates communication from
     * matrix-vector multiplication.  Physical basis values are exchanged
     * explicitly through sys->monolis_com, while A*x is evaluated with this
     * communication-free communicator.
     */
    monolis_com_initialize_by_self(&sys->mono_com0);

    local_count = (size_t)context->local_st_dof;
    trace_count = (size_t)sys->fe.total_num_nodes;

    context->zero_previous_trace = (double*)calloc(trace_count, sizeof(double));
    context->previous_trace = (double*)calloc(trace_count, sizeof(double));
    context->previous_state_work = (double*)calloc(local_count, sizeof(double));
    context->rhs_with_previous = (double*)calloc(local_count, sizeof(double));
    context->rhs_without_previous = (double*)calloc(local_count, sizeof(double));
    context->current_lift = (double*)calloc(local_count, sizeof(double));
    context->previous_lift = (double*)calloc(local_count, sizeof(double));
    context->operator_product = (double*)calloc(local_count, sizeof(double));

    if(context->zero_previous_trace == NULL ||
       context->previous_trace == NULL ||
       context->previous_state_work == NULL ||
       context->rhs_with_previous == NULL ||
       context->rhs_without_previous == NULL ||
       context->current_lift == NULL ||
       context->previous_lift == NULL ||
       context->operator_product == NULL)
    {
        ROM_std_stdd_stsb_context_finalize(context);
        return ROM_STDD_ERR_ALLOCATION;
    }

    return ROM_STDD_SUCCESS;
}

void ROM_std_stdd_stsb_context_finalize(
    ROM_STDD_STSB_CONTEXT* context)
{
    if(context == NULL){
        return;
    }

    free(context->zero_previous_trace);
    free(context->previous_trace);
    free(context->previous_state_work);
    free(context->rhs_with_previous);
    free(context->rhs_without_previous);
    free(context->current_lift);
    free(context->previous_lift);
    free(context->operator_product);

    memset(context, 0, sizeof(*context));
}

int ROM_std_stdd_stsb_prepare_operator(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id)
{
    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       window_id < 0 || window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /*
     * Project the exact boundary-modified operator used by the reference FOM.
     * The previous trace changes only the RHS, so a zero trace is sufficient
     * when preparing A for basis sweeps.
     */
    return stdd_stsb_build_strong_pair(
        context,
        window_id,
        context->zero_previous_trace,
        NULL);
}

int ROM_std_stdd_stsb_apply_prepared_A(
    const double* x,
    double* y,
    void* user_context)
{
    ROM_STDD_STSB_CONTEXT* context =
        (ROM_STDD_STSB_CONTEXT*)user_context;
    double* x_work;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       x == NULL || y == NULL ||
       context->previous_state_work == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /*
     * Verified HLPOD-DDROM communication convention:
     *
     *   1. explicitly exchange the owned basis/state values with the physical
     *      communicator,
     *   2. evaluate the local sparse matrix product with mono_com0, which has
     *      no MPI communication.
     *
     * This is essential during the source-rank basis sweep.  For one source
     * rank, x contains that rank's owned basis values and zeros elsewhere.
     * The explicit update places those values only in neighbouring halo rows;
     * the communication-free matvec then produces the separate A*Phi_source
     * column used for V_i^T A V_source.
     */
    x_work = context->previous_state_work;

    if(x_work != x){
        memcpy(
            x_work,
            x,
            (size_t)context->local_st_dof * sizeof(double));
    }

    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->sys->fe.total_num_nodes,
        context->st_dof_per_space_node,
        x_work);

    memset(y, 0, (size_t)context->local_st_dof * sizeof(double));

    monolis_matvec_product_R(
        &context->sys->monolis,
        &context->sys->mono_com0,
        x_work,
        y);

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_build_dirichlet_lift(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    double* lift)
{
    double window_start_time;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS || lift == NULL ||
       window_id < 0 || window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(lift, 0, (size_t)context->local_st_dof * sizeof(double));
    window_start_time = stdd_stsb_window_start_time(context, window_id);

    STSB_restore_dirichlet_values(
        &context->sys->fe,
        &context->sys->bc,
        &context->sys->vals,
        context->time_element,
        context->num_window_slabs,
        window_start_time,
        lift);

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_enforce_homogeneous_basis(
    ROM_STDD_STSB_CONTEXT* context,
    ROM_STDD_SYSTEM* system)
{
    const double minimum_norm = 128.0 * DBL_EPSILON;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       system == NULL || system->basis == NULL || system->num_modes <= 0 ||
       system->num_modes > system->num_modes_capacity ||
       system->internal_st_dof !=
           context->sys->monolis.mat.N * context->st_dof_per_space_node)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /* Enforce q=0 on every owned Dirichlet space node exactly. */
    for(int row = 0; row < system->internal_st_dof; row++){
        const int node = row / context->st_dof_per_space_node;
        if(!context->sys->bc.D_bc_exists[node]){
            continue;
        }
        for(int mode = 0; mode < system->num_modes; mode++){
            system->basis[
                (size_t)row * (size_t)system->num_modes_capacity
                + (size_t)mode] = 0.0;
        }
    }

    /* Re-orthonormalize after exact boundary zeroing. */
    for(int mode = 0; mode < system->num_modes; mode++){
        for(int pass = 0; pass < 2; pass++){
            for(int previous = 0; previous < mode; previous++){
                double dot = 0.0;
                for(int row = 0; row < system->internal_st_dof; row++){
                    dot += system->basis[
                        (size_t)row * (size_t)system->num_modes_capacity
                        + (size_t)previous]
                        * system->basis[
                            (size_t)row * (size_t)system->num_modes_capacity
                            + (size_t)mode];
                }
                for(int row = 0; row < system->internal_st_dof; row++){
                    system->basis[
                        (size_t)row * (size_t)system->num_modes_capacity
                        + (size_t)mode] -=
                        dot * system->basis[
                            (size_t)row * (size_t)system->num_modes_capacity
                            + (size_t)previous];
                }
            }
        }

        double norm2 = 0.0;
        for(int row = 0; row < system->internal_st_dof; row++){
            const double value = system->basis[
                (size_t)row * (size_t)system->num_modes_capacity
                + (size_t)mode];
            norm2 += value * value;
        }

        if(!isfinite(norm2) || norm2 <= minimum_norm * minimum_norm){
            return ROM_STDD_ERR_DIMENSION;
        }

        const double inverse_norm = 1.0 / sqrt(norm2);
        for(int row = 0; row < system->internal_st_dof; row++){
            system->basis[
                (size_t)row * (size_t)system->num_modes_capacity
                + (size_t)mode] *= inverse_norm;
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_build_raw_rhs(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* previous_trace,
    double* full_rhs)
{
    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       previous_trace == NULL || full_rhs == NULL ||
       window_id < 0 || window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    /*
     * This is deliberately a full strong-BC build, not a primitive RHS
     * assembly.  On return, full_rhs is b_tilde and the working MONOLIS
     * matrix is the matching A_tilde.
     */
    return stdd_stsb_build_strong_pair(
        context,
        window_id,
        previous_trace,
        full_rhs);
}

int ROM_std_stdd_stsb_apply_B(
    const double* x_previous_window,
    double* y_current_window,
    void* user_context)
{
    ROM_STDD_STSB_CONTEXT* context =
        (ROM_STDD_STSB_CONTEXT*)user_context;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       x_previous_window == NULL || y_current_window == NULL)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memcpy(
        context->previous_state_work,
        x_previous_window,
        (size_t)context->local_st_dof * sizeof(double));

    /* Fill halo values of the source-rank basis before trace extraction. */
    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->sys->fe.total_num_nodes,
        context->st_dof_per_space_node,
        context->previous_state_work);

    STSB_extract_last_right_trace(
        context->time_element,
        context->num_window_slabs,
        context->sys->fe.total_num_nodes,
        context->previous_state_work,
        context->previous_trace);

    if(ROM_std_stdd_stsb_build_raw_rhs(
        context,
        0,
        context->previous_trace,
        context->rhs_with_previous) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_build_raw_rhs(
        context,
        0,
        context->zero_previous_trace,
        context->rhs_without_previous) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    for(int i = 0; i < context->local_st_dof; i++){
        y_current_window[i] =
            context->rhs_with_previous[i]
            - context->rhs_without_previous[i];
    }

    /* build_raw_rhs changed the working matrix; restore strong A. */
    if(ROM_std_stdd_stsb_prepare_operator(context, 0)
       != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_build_effective_rhs(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* initial_trace,
    double* effective_rhs)
{
    const double* previous_trace;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       initial_trace == NULL || effective_rhs == NULL ||
       window_id < 0 || window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ROM_std_stdd_stsb_build_dirichlet_lift(
        context,
        window_id,
        context->current_lift) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    if(window_id == 0){
        previous_trace = initial_trace;
    }
    else{
        if(ROM_std_stdd_stsb_build_dirichlet_lift(
            context,
            window_id - 1,
            context->previous_lift) != ROM_STDD_SUCCESS)
        {
            return ROM_STDD_ERR_OPERATOR;
        }

        STSB_extract_last_right_trace(
            context->time_element,
            context->num_window_slabs,
            context->sys->fe.total_num_nodes,
            context->previous_lift,
            context->previous_trace);

        previous_trace = context->previous_trace;
    }

    /*
     * build_raw_rhs() leaves the matched strong-BC matrix A_tilde in
     * context->sys->monolis.  Apply THAT SAME matrix to the lift; do not
     * rebuild the operator in a second boundary-condition pass.
     */
    if(ROM_std_stdd_stsb_build_raw_rhs(
        context,
        window_id,
        previous_trace,
        context->rhs_with_previous) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_apply_prepared_A(
        context->current_lift,
        context->operator_product,
        context) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    for(int i = 0; i < context->local_st_dof; i++){
        effective_rhs[i] =
            context->rhs_with_previous[i]
            - context->operator_product[i];
    }

    return ROM_STDD_SUCCESS;
}


int ROM_std_stdd_stsb_evaluate_raw_full_state_residual(
    ROM_STDD_STSB_CONTEXT* context,
    int window_id,
    const double* previous_trace,
    const double* full_window_state,
    ROM_STDD_STSB_RESIDUAL_SUM* residual_sum)
{
    double* state_work = NULL;
    double* raw_rhs = NULL;
    double* operator_product = NULL;
    int internal_st_dof = 0;
    int status = ROM_STDD_SUCCESS;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       previous_trace == NULL || full_window_state == NULL ||
       residual_sum == NULL || window_id < 0 ||
       window_id >= context->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(residual_sum, 0, sizeof(*residual_sum));

    state_work = (double*)calloc(
        (size_t)context->local_st_dof, sizeof(double));
    raw_rhs = (double*)calloc(
        (size_t)context->local_st_dof, sizeof(double));
    operator_product = (double*)calloc(
        (size_t)context->local_st_dof, sizeof(double));

    if(state_work == NULL || raw_rhs == NULL || operator_product == NULL){
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    memcpy(
        state_work,
        full_window_state,
        (size_t)context->local_st_dof * sizeof(double));

    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->sys->fe.total_num_nodes,
        context->st_dof_per_space_node,
        state_work);

    /*
     * build_raw_rhs() creates the matched pair (A_tilde, b_tilde).
     * Evaluate A_tilde * U with the same prepared matrix instance.
     */
    if(ROM_std_stdd_stsb_build_raw_rhs(
        context, window_id, previous_trace, raw_rhs) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_apply_prepared_A(
        state_work, operator_product, context) != ROM_STDD_SUCCESS)
    {
        status = ROM_STDD_ERR_OPERATOR;
        goto cleanup;
    }

    internal_st_dof =
        context->sys->monolis.mat.N * context->st_dof_per_space_node;

    for(int row = 0; row < internal_st_dof; row++){
        const int node = row / context->st_dof_per_space_node;
        const double residual = operator_product[row] - raw_rhs[row];
        const double rhs_value = raw_rhs[row];
        const int is_dirichlet =
            context->sys->bc.D_bc_exists[node] ? 1 : 0;

        residual_sum->residual_all2 += residual * residual;
        residual_sum->effective_rhs_all2 += rhs_value * rhs_value;

        if(is_dirichlet){
            residual_sum->residual_dirichlet2 += residual * residual;
            residual_sum->effective_rhs_dirichlet2 +=
                rhs_value * rhs_value;
        }
        else{
            residual_sum->residual_free2 += residual * residual;
            residual_sum->effective_rhs_free2 += rhs_value * rhs_value;
        }
    }

cleanup:
    free(state_work);
    free(raw_rhs);
    free(operator_product);
    return status;
}


int ROM_std_stdd_stsb_evaluate_lifted_reference_residual(
    ROM_STDD_STSB_CONTEXT* context,
    const ROM_STDD_SYSTEM* system,
    int window_id,
    const double* initial_trace,
    const double* current_reference_internal,
    const double* previous_reference_internal,
    double* direct_projected_defect,
    ROM_STDD_STSB_RESIDUAL_SUM* residual_sum)
{
    double* current_homogeneous = NULL;
    double* previous_homogeneous = NULL;
    double* previous_full_state = NULL;
    double* previous_trace = NULL;
    double* effective_rhs = NULL;
    double* operator_product = NULL;
    int status = ROM_STDD_SUCCESS;

    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       system == NULL || initial_trace == NULL ||
       current_reference_internal == NULL ||
       direct_projected_defect == NULL || residual_sum == NULL ||
       window_id < 0 || window_id >= context->num_windows ||
       system->internal_st_dof !=
           context->sys->monolis.mat.N * context->st_dof_per_space_node ||
       system->local_st_dof != context->local_st_dof ||
       system->num_modes <= 0 || system->basis == NULL ||
       (window_id > 0 && previous_reference_internal == NULL))
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memset(residual_sum, 0, sizeof(*residual_sum));
    memset(
        direct_projected_defect,
        0,
        (size_t)system->num_modes * sizeof(double));

    current_homogeneous = (double*)calloc(
        (size_t)context->local_st_dof,
        sizeof(double));
    previous_homogeneous = (double*)calloc(
        (size_t)context->local_st_dof,
        sizeof(double));
    previous_full_state = (double*)calloc(
        (size_t)context->local_st_dof,
        sizeof(double));
    previous_trace = (double*)calloc(
        (size_t)context->sys->fe.total_num_nodes,
        sizeof(double));
    effective_rhs = (double*)calloc(
        (size_t)context->local_st_dof,
        sizeof(double));
    operator_product = (double*)calloc(
        (size_t)context->local_st_dof,
        sizeof(double));

    if(current_homogeneous == NULL ||
       previous_homogeneous == NULL ||
       previous_full_state == NULL ||
       previous_trace == NULL ||
       effective_rhs == NULL ||
       operator_product == NULL)
    {
        status = ROM_STDD_ERR_ALLOCATION;
        goto cleanup;
    }

    memcpy(
        current_homogeneous,
        current_reference_internal,
        (size_t)system->internal_st_dof * sizeof(double));

    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->sys->fe.total_num_nodes,
        context->st_dof_per_space_node,
        current_homogeneous);

    if(window_id == 0){
        memcpy(
            previous_trace,
            initial_trace,
            (size_t)context->sys->fe.total_num_nodes * sizeof(double));
    }
    else{
        memcpy(
            previous_homogeneous,
            previous_reference_internal,
            (size_t)system->internal_st_dof * sizeof(double));

        monolis_mpi_update_R(
            &context->sys->monolis_com,
            context->sys->fe.total_num_nodes,
            context->st_dof_per_space_node,
            previous_homogeneous);

        if(ROM_std_stdd_stsb_build_dirichlet_lift(
            context,
            window_id - 1,
            context->previous_lift) != ROM_STDD_SUCCESS)
        {
            status = ROM_STDD_ERR_OPERATOR;
            goto cleanup;
        }

        for(int row = 0; row < context->local_st_dof; row++){
            previous_full_state[row] =
                previous_homogeneous[row] + context->previous_lift[row];
        }

        /*
         * The trace extraction reads local and halo nodes.  Synchronize the
         * combined physical state after adding q and U_D, rather than relying
         * on the two operands having independently consistent halo values.
         */
        monolis_mpi_update_R(
            &context->sys->monolis_com,
            context->sys->fe.total_num_nodes,
            context->st_dof_per_space_node,
            previous_full_state);

        STSB_extract_last_right_trace(
            context->time_element,
            context->num_window_slabs,
            context->sys->fe.total_num_nodes,
            previous_full_state,
            previous_trace);
    }

    /*
     * build_raw_rhs() leaves A_tilde paired with b_tilde.  Use the same
     * prepared A_tilde for both A_tilde*U_D and A_tilde*q.
     */
    if(ROM_std_stdd_stsb_build_dirichlet_lift(
        context,
        window_id,
        context->current_lift) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_build_raw_rhs(
        context,
        window_id,
        previous_trace,
        context->rhs_with_previous) != ROM_STDD_SUCCESS ||
       ROM_std_stdd_stsb_apply_prepared_A(
        context->current_lift,
        context->operator_product,
        context) != ROM_STDD_SUCCESS)
    {
        status = ROM_STDD_ERR_OPERATOR;
        goto cleanup;
    }

    for(int row = 0; row < context->local_st_dof; row++){
        effective_rhs[row] =
            context->rhs_with_previous[row]
            - context->operator_product[row];
    }

    if(ROM_std_stdd_stsb_apply_prepared_A(
        current_homogeneous,
        operator_product,
        context) != ROM_STDD_SUCCESS)
    {
        status = ROM_STDD_ERR_OPERATOR;
        goto cleanup;
    }

    for(int row = 0; row < system->internal_st_dof; row++){
        const int node = row / context->st_dof_per_space_node;
        const double residual =
            operator_product[row] - effective_rhs[row];
        const double rhs_value = effective_rhs[row];
        const double q_value = current_homogeneous[row];
        const int is_dirichlet =
            context->sys->bc.D_bc_exists[node] ? 1 : 0;

        residual_sum->residual_all2 += residual * residual;
        residual_sum->effective_rhs_all2 += rhs_value * rhs_value;
        residual_sum->homogeneous_all2 += q_value * q_value;

        if(is_dirichlet){
            residual_sum->residual_dirichlet2 += residual * residual;
            residual_sum->effective_rhs_dirichlet2 +=
                rhs_value * rhs_value;
            residual_sum->homogeneous_dirichlet2 += q_value * q_value;
        }
        else{
            residual_sum->residual_free2 += residual * residual;
            residual_sum->effective_rhs_free2 += rhs_value * rhs_value;
        }

        if(!is_dirichlet){
            for(int mode = 0; mode < system->num_modes; mode++){
                direct_projected_defect[mode] +=
                    system->basis[
                        (size_t)row * (size_t)system->num_modes_capacity
                        + (size_t)mode]
                    * residual;
            }
        }
    }

cleanup:
    free(current_homogeneous);
    free(previous_homogeneous);
    free(previous_full_state);
    free(previous_trace);
    free(effective_rhs);
    free(operator_product);

    return status;
}

int ROM_std_stdd_stsb_store_internal_homogeneous(
    const ROM_STDD_SYSTEM* system,
    const double* full_window_state,
    const double* lift,
    double* internal_snapshot)
{
    if(system == NULL || full_window_state == NULL || lift == NULL ||
       internal_snapshot == NULL || system->internal_st_dof <= 0)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int row = 0; row < system->internal_st_dof; row++){
        internal_snapshot[row] = full_window_state[row] - lift[row];
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_stsb_reconstruct_full_window(
    ROM_STDD_STSB_CONTEXT* context,
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* homogeneous_work,
    double* full_window_state)
{
    if(stdd_stsb_validate(context) != ROM_STDD_SUCCESS ||
       system == NULL || homogeneous_work == NULL ||
       full_window_state == NULL || window_id < 0 ||
       window_id >= system->num_windows)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(ROM_std_stdd_reconstruct_local_homogeneous(
        system,
        window_id,
        homogeneous_work) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    monolis_mpi_update_R(
        &context->sys->monolis_com,
        context->sys->fe.total_num_nodes,
        context->st_dof_per_space_node,
        homogeneous_work);

    if(ROM_std_stdd_stsb_build_dirichlet_lift(
        context,
        window_id,
        context->current_lift) != ROM_STDD_SUCCESS)
    {
        return ROM_STDD_ERR_OPERATOR;
    }

    for(int i = 0; i < context->local_st_dof; i++){
        full_window_state[i] =
            homogeneous_work[i] + context->current_lift[i];
    }

    return ROM_STDD_SUCCESS;
}
