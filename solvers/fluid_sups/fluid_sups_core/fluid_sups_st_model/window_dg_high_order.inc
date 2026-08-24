/* ========================================================================== */
/* Higher-order dG(1)/dG(2) temporal weak form                               */
/* ========================================================================== */

static double fluid_time_basis_derivative(
    const ST_FOM_TIME_ELEMENT* time,
    int beta,
    double tau)
{
    if(time == NULL || beta < 0 || beta >= time->n_dof){
        return 0.0;
    }

    if(time->degree == 0){
        return 0.0;
    }
    if(time->degree == 1){
        return beta == 0 ? -1.0 : 1.0;
    }
    if(time->degree == 2){
        if(beta == 0){
            return 4.0 * tau - 3.0;
        }
        if(beta == 1){
            return -8.0 * tau + 4.0;
        }
        return 4.0 * tau - 1.0;
    }
    return 0.0;
}

static size_t fluid_local_time_index(
    int nt,
    int nl,
    int beta,
    int local_node,
    int component)
{
    (void)nt;
    return (((size_t)beta * (size_t)nl + (size_t)local_node)
        * (size_t)FLUID_SUPS_PHYSICAL_DOF)
        + (size_t)component;
}


typedef struct {
    int nl;
    int np;
    int nq;
    int nt;

    double* Jacobian_ip;
    double (*J_inv)[3][3];

    /* [time quadrature][space quadrature] */
    double (*v)[3];
    double (*du_time)[3];
    double (*grad_v)[3][3];
    double* pressure;
    double (*grad_pressure)[3];
    double* tau;
    double* tau_c;

    /* left trace and previous right trace, at space quadrature points */
    double (*v_left)[3];
    double (*grad_v_left)[3][3];
    double* pressure_left;
    double (*grad_pressure_left)[3];
    double (*v_previous)[3];
    double* tau_left;
    double* tau_c_left;

    /* reusable spatial nodal work */
    double** local_v;
    double** local_du;
    double* local_p;
    double** local_v_left;
    double* local_p_left;
    double** local_v_previous;
} FLUID_DG_ELEMENT_FIELDS;

static void fluid_dg_element_fields_clear(
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    if(fields != NULL){
        memset(fields, 0, sizeof(*fields));
    }
}

static int fluid_dg_element_fields_allocate(
    FLUID_DG_ELEMENT_FIELDS* fields,
    int nl,
    int np,
    int nq,
    int nt)
{
    const size_t nqp = (size_t)nq * (size_t)np;

    if(fields == NULL || nl <= 0 || np <= 0 ||
       nq <= 0 || nt <= 0)
    {
        return 1;
    }

    fluid_dg_element_fields_clear(fields);
    fields->nl = nl;
    fields->np = np;
    fields->nq = nq;
    fields->nt = nt;

    fields->Jacobian_ip =
        (double*)calloc((size_t)np, sizeof(double));
    fields->J_inv =
        (double (*)[3][3])calloc(
            (size_t)np, sizeof(*fields->J_inv));

    fields->v =
        (double (*)[3])calloc(nqp, sizeof(*fields->v));
    fields->du_time =
        (double (*)[3])calloc(nqp, sizeof(*fields->du_time));
    fields->grad_v =
        (double (*)[3][3])calloc(
            nqp, sizeof(*fields->grad_v));
    fields->pressure =
        (double*)calloc(nqp, sizeof(double));
    fields->grad_pressure =
        (double (*)[3])calloc(
            nqp, sizeof(*fields->grad_pressure));
    fields->tau =
        (double*)calloc(nqp, sizeof(double));
    fields->tau_c =
        (double*)calloc(nqp, sizeof(double));

    fields->v_left =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->v_left));
    fields->grad_v_left =
        (double (*)[3][3])calloc(
            (size_t)np, sizeof(*fields->grad_v_left));
    fields->pressure_left =
        (double*)calloc((size_t)np, sizeof(double));
    fields->grad_pressure_left =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->grad_pressure_left));
    fields->v_previous =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->v_previous));
    fields->tau_left =
        (double*)calloc((size_t)np, sizeof(double));
    fields->tau_c_left =
        (double*)calloc((size_t)np, sizeof(double));

    fields->local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_du =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_p =
        BB_std_calloc_1d_double(NULL, nl);
    fields->local_v_left =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_p_left =
        BB_std_calloc_1d_double(NULL, nl);
    fields->local_v_previous =
        BB_std_calloc_2d_double(NULL, nl, 3);

    if(fields->Jacobian_ip == NULL ||
       fields->J_inv == NULL ||
       fields->v == NULL ||
       fields->du_time == NULL ||
       fields->grad_v == NULL ||
       fields->pressure == NULL ||
       fields->grad_pressure == NULL ||
       fields->tau == NULL ||
       fields->tau_c == NULL ||
       fields->v_left == NULL ||
       fields->grad_v_left == NULL ||
       fields->pressure_left == NULL ||
       fields->grad_pressure_left == NULL ||
       fields->v_previous == NULL ||
       fields->tau_left == NULL ||
       fields->tau_c_left == NULL ||
       fields->local_v == NULL ||
       fields->local_du == NULL ||
       fields->local_p == NULL ||
       fields->local_v_left == NULL ||
       fields->local_p_left == NULL ||
       fields->local_v_previous == NULL)
    {
        return 1;
    }

    return 0;
}

static void fluid_dg_element_fields_free(
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    if(fields == NULL){
        return;
    }

    free(fields->Jacobian_ip);
    free(fields->J_inv);
    free(fields->v);
    free(fields->du_time);
    free(fields->grad_v);
    free(fields->pressure);
    free(fields->grad_pressure);
    free(fields->tau);
    free(fields->tau_c);
    free(fields->v_left);
    free(fields->grad_v_left);
    free(fields->pressure_left);
    free(fields->grad_pressure_left);
    free(fields->v_previous);
    free(fields->tau_left);
    free(fields->tau_c_left);

    if(fields->local_v != NULL){
        BB_std_free_2d_double(fields->local_v, fields->nl, 3);
    }
    if(fields->local_du != NULL){
        BB_std_free_2d_double(fields->local_du, fields->nl, 3);
    }
    if(fields->local_p != NULL){
        BB_std_free_1d_double(fields->local_p, fields->nl);
    }
    if(fields->local_v_left != NULL){
        BB_std_free_2d_double(
            fields->local_v_left, fields->nl, 3);
    }
    if(fields->local_p_left != NULL){
        BB_std_free_1d_double(
            fields->local_p_left, fields->nl);
    }
    if(fields->local_v_previous != NULL){
        BB_std_free_2d_double(
            fields->local_v_previous, fields->nl, 3);
    }

    fluid_dg_element_fields_clear(fields);
}

static int fluid_dg_build_element_fields(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    const ST_FOM_TIME_ELEMENT* time,
    int elem,
    const double* coefficient,
    const double* previous_trace_local,
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    const int nl = fields->nl;
    const int np = fields->np;
    const int nq = fields->nq;
    const int nt = fields->nt;
    int q;
    int p;
    int j;
    int beta;
    int d;

    if(fe == NULL || basis == NULL || vals == NULL ||
       time == NULL || coefficient == NULL ||
       previous_trace_local == NULL || fields == NULL)
    {
        return 1;
    }

    BBFE_elemmat_set_Jacobian_array(
        fields->Jacobian_ip, np, elem, fe);

    for(p = 0; p < np; p++){
        BB_calc_mat3d_inverse(
            fe->geo[elem][p].J,
            fe->geo[elem][p].Jacobian,
            fields->J_inv[p]);
    }

    /*
     * Temporal quadrature-point state:
     *
     *   v_q       = sum_beta psi_beta v_beta
     *   du_time_q = sum_beta dpsi_beta/dtau v_beta
     *
     * Since physical dv/dt = du_time_q/dt, the verified BE kernel can be
     * evaluated at the polynomial derivative by setting
     *
     *   v_old,eff = v_q - du_time_q.
     */
    for(q = 0; q < nq; q++){
        const double tau_time = time->tau[q];

        for(j = 0; j < nl; j++){
            for(d = 0; d < 3; d++){
                double value = 0.0;
                double derivative = 0.0;

                for(beta = 0; beta < nt; beta++){
                    const double coeff =
                        coefficient[
                            fluid_local_time_index(
                                nt, nl, beta, j, d)];
                    value += time->psi[beta][q] * coeff;
                    derivative +=
                        fluid_time_basis_derivative(
                            time, beta, tau_time) * coeff;
                }

                fields->local_v[j][d] = value;
                fields->local_du[j][d] = derivative;
            }

            fields->local_p[j] = 0.0;
            for(beta = 0; beta < nt; beta++){
                fields->local_p[j] +=
                    time->psi[beta][q]
                    * coefficient[
                        fluid_local_time_index(
                            nt, nl, beta, j, 3)];
            }
        }

        for(p = 0; p < np; p++){
            const size_t qp =
                (size_t)q * (size_t)np + (size_t)p;
            double grad_storage[3][3];
            double* grad_rows[3] = {
                grad_storage[0],
                grad_storage[1],
                grad_storage[2]
            };
            double grad_p[3];

            BBFE_std_mapping_vector3d(
                fields->v[qp],
                nl,
                fields->local_v,
                basis->N[p]);
            BBFE_std_mapping_vector3d(
                fields->du_time[qp],
                nl,
                fields->local_du,
                basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_rows,
                nl,
                fields->local_v,
                fe->geo[elem][p].grad_N);

            for(d = 0; d < 3; d++){
                int k;
                fields->grad_pressure[qp][d] = 0.0;
                for(k = 0; k < 3; k++){
                    fields->grad_v[qp][d][k] =
                        grad_storage[d][k];
                }
            }

            fields->pressure[qp] =
                BBFE_std_mapping_scalar(
                    nl,
                    fields->local_p,
                    basis->N[p]);
            BBFE_std_mapping_scalar_grad(
                grad_p,
                nl,
                fields->local_p,
                fe->geo[elem][p].grad_N);
            for(d = 0; d < 3; d++){
                fields->grad_pressure[qp][d] = grad_p[d];
            }

            fields->tau[qp] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    fields->J_inv[p],
                    fe->geo[elem][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    fields->v[qp],
                    vals->dt);
            fields->tau_c[qp] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    fields->J_inv[p],
                    fe->geo[elem][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    fields->v[qp],
                    vals->dt);
        }
    }

    /* Left trace of current slab and previous right trace. */
    for(j = 0; j < nl; j++){
        for(d = 0; d < 3; d++){
            fields->local_v_left[j][d] = 0.0;
            for(beta = 0; beta < nt; beta++){
                fields->local_v_left[j][d] +=
                    time->e_minus[beta]
                    * coefficient[
                        fluid_local_time_index(
                            nt, nl, beta, j, d)];
            }
            fields->local_v_previous[j][d] =
                previous_trace_local[
                    (size_t)j * FLUID_SUPS_PHYSICAL_DOF
                    + (size_t)d];
        }

        fields->local_p_left[j] = 0.0;
        for(beta = 0; beta < nt; beta++){
            fields->local_p_left[j] +=
                time->e_minus[beta]
                * coefficient[
                    fluid_local_time_index(
                        nt, nl, beta, j, 3)];
        }
    }

    for(p = 0; p < np; p++){
        double grad_storage[3][3];
        double* grad_rows[3] = {
            grad_storage[0],
            grad_storage[1],
            grad_storage[2]
        };
        double grad_p[3];

        BBFE_std_mapping_vector3d(
            fields->v_left[p],
            nl,
            fields->local_v_left,
            basis->N[p]);
        BBFE_std_mapping_vector3d_grad(
            grad_rows,
            nl,
            fields->local_v_left,
            fe->geo[elem][p].grad_N);
        for(d = 0; d < 3; d++){
            int k;
            for(k = 0; k < 3; k++){
                fields->grad_v_left[p][d][k] =
                    grad_storage[d][k];
            }
        }

        fields->pressure_left[p] =
            BBFE_std_mapping_scalar(
                nl,
                fields->local_p_left,
                basis->N[p]);
        BBFE_std_mapping_scalar_grad(
            grad_p,
            nl,
            fields->local_p_left,
            fe->geo[elem][p].grad_N);
        for(d = 0; d < 3; d++){
            fields->grad_pressure_left[p][d] = grad_p[d];
        }

        BBFE_std_mapping_vector3d(
            fields->v_previous[p],
            nl,
            fields->local_v_previous,
            basis->N[p]);

        fields->tau_left[p] =
            BBFE_elemmat_fluid_sups_coef_metric_tensor(
                fields->J_inv[p],
                fe->geo[elem][p].Jacobian,
                vals->density,
                vals->viscosity,
                fields->v_left[p],
                vals->dt);
        fields->tau_c_left[p] =
            BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                fields->J_inv[p],
                fe->geo[elem][p].Jacobian,
                vals->density,
                vals->viscosity,
                fields->v_left[p],
                vals->dt);
    }

    return 0;
}

/*
 * Directional derivative of the verified residual with respect to an old
 * velocity nodal coefficient.  A BE old-state perturbation produces
 *
 *   delta v_old_ip = N_j delta,
 *   delta du_time  = -N_j delta.
 *
 * Hence this derivative includes every transient SUPS/PSPG contribution
 * present in the verified kernel, without replacing it by -M/dt.
 */
static void fluid_dg_old_jacobian_ip(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    const double* velocity,
    const double* velocity_old,
    double** grad_velocity,
    double pressure,
    const double* grad_pressure,
    double tau,
    double tau_c,
    const double* du_time,
    double fd_epsilon,
    double J_old[4][3])
{
    (void)fd_epsilon;

    fluid_old_jacobian_ip(
        fe,
        basis,
        vals,
        elem,
        test_node,
        trial_node,
        p,
        velocity,
        velocity_old,
        grad_velocity,
        pressure,
        grad_pressure,
        tau,
        tau_c,
        du_time,
        J_old);
}

static int fluid_dg_gather_element_coefficients(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int elem,
    int slab,
    double* coefficient,
    double* previous_trace_local)
{
    const int nl = system->fe.local_num_nodes;
    const int nt = layout->time.n_dof;
    int j;
    int beta;
    int component;

    for(j = 0; j < nl; j++){
        const int node = system->fe.conn[elem][j];

        for(beta = 0; beta < nt; beta++){
            for(component = 0;
                component < FLUID_SUPS_PHYSICAL_DOF;
                component++)
            {
                coefficient[
                    fluid_local_time_index(
                        nt, nl, beta, j, component)] =
                    window_state[
                        ST_fom_vector_index(
                            layout,
                            node,
                            slab,
                            beta,
                            component)];
            }
        }

        for(component = 0;
            component < FLUID_SUPS_PHYSICAL_DOF;
            component++)
        {
            double value = 0.0;

            if(slab == 0){
                value =
                    previous_trace[
                        (size_t)node
                            * FLUID_SUPS_PHYSICAL_DOF
                        + (size_t)component];
            }
            else{
                for(beta = 0; beta < nt; beta++){
                    value +=
                        layout->time.e_plus[beta]
                        * window_state[
                            ST_fom_vector_index(
                                layout,
                                node,
                                slab - 1,
                                beta,
                                component)];
                }
            }

            previous_trace_local[
                (size_t)j * FLUID_SUPS_PHYSICAL_DOF
                + (size_t)component] = value;
        }
    }

    return 0;
}

static int fluid_dg_add_element_slab(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int elem,
    int slab,
    int assemble_jacobian)
{
    const int nl = system->fe.local_num_nodes;
    const int np = system->basis.num_integ_points;
    const int nt = layout->time.n_dof;
    const int nq = layout->time.num_time_ip;
    const int local_dof =
        nt * nl * FLUID_SUPS_PHYSICAL_DOF;
    const size_t residual_ip_size =
        (size_t)nt * FLUID_SUPS_PHYSICAL_DOF * (size_t)np;
    const size_t jac_ip_size =
        (size_t)nt * (size_t)nt
        * FLUID_SUPS_PHYSICAL_DOF
        * FLUID_SUPS_PHYSICAL_DOF
        * (size_t)np;
    const size_t prev_jac_ip_size =
        (size_t)nt * (size_t)nt
        * FLUID_SUPS_PHYSICAL_DOF
        * 3u
        * (size_t)np;

    double* coefficient = NULL;
    double* previous_trace_local = NULL;
    double* residual_ip = NULL;
    double* jump_ip = NULL;
    double* jac_ip = NULL;
    double* prev_jac_ip = NULL;
    FLUID_DG_ELEMENT_FIELDS fields;
    int i;
    int j;
    int p;
    int q;
    int alpha;
    int beta;
    int a;
    int b;
    int status = 1;

    fluid_dg_element_fields_clear(&fields);

    coefficient =
        (double*)calloc((size_t)local_dof, sizeof(double));
    previous_trace_local =
        (double*)calloc(
            (size_t)nl * FLUID_SUPS_PHYSICAL_DOF,
            sizeof(double));
    residual_ip =
        (double*)calloc(residual_ip_size, sizeof(double));
    jump_ip =
        (double*)calloc(
            (size_t)FLUID_SUPS_PHYSICAL_DOF * (size_t)np,
            sizeof(double));

    if(assemble_jacobian){
        jac_ip =
            (double*)calloc(jac_ip_size, sizeof(double));
        if(slab > 0){
            prev_jac_ip =
                (double*)calloc(
                    prev_jac_ip_size, sizeof(double));
        }
    }

    if(coefficient == NULL ||
       previous_trace_local == NULL ||
       residual_ip == NULL ||
       jump_ip == NULL ||
       (assemble_jacobian && jac_ip == NULL) ||
       (assemble_jacobian && slab > 0 && prev_jac_ip == NULL) ||
       fluid_dg_element_fields_allocate(
           &fields, nl, np, nq, nt) != 0)
    {
        goto cleanup;
    }

    fluid_dg_gather_element_coefficients(
        system,
        layout,
        previous_trace,
        window_state,
        elem,
        slab,
        coefficient,
        previous_trace_local);

    if(fluid_dg_build_element_fields(
        &system->fe,
        &system->basis,
        &system->vals,
        &layout->time,
        elem,
        coefficient,
        previous_trace_local,
        &fields) != 0)
    {
        goto cleanup;
    }

#define FLUID_RES_IP(alpha_, a_, p_) \
    residual_ip[(((size_t)(alpha_) * 4u + (size_t)(a_)) \
        * (size_t)np) + (size_t)(p_)]
#define FLUID_JUMP_IP(a_, p_) \
    jump_ip[((size_t)(a_) * (size_t)np) + (size_t)(p_)]
#define FLUID_JAC_IP(alpha_, beta_, a_, b_, p_) \
    jac_ip[((((((size_t)(alpha_) * (size_t)nt \
        + (size_t)(beta_)) * 4u + (size_t)(a_)) * 4u \
        + (size_t)(b_)) * (size_t)np) + (size_t)(p_))]
#define FLUID_PREV_JAC_IP(alpha_, beta_, a_, b_, p_) \
    prev_jac_ip[((((((size_t)(alpha_) * (size_t)nt \
        + (size_t)(beta_)) * 4u + (size_t)(a_)) * 3u \
        + (size_t)(b_)) * (size_t)np) + (size_t)(p_))]

    for(i = 0; i < nl; i++){
        memset(
            residual_ip,
            0,
            residual_ip_size * sizeof(double));
        memset(
            jump_ip,
            0,
            (size_t)FLUID_SUPS_PHYSICAL_DOF
                * (size_t)np
                * sizeof(double));

        /*
         * Interior temporal weak residual:
         *
         *   int_0^1 psi_alpha R_SUPS(U(tau), U_t(tau)) d tau.
         */
        for(q = 0; q < nq; q++){
            for(p = 0; p < np; p++){
                const size_t qp =
                    (size_t)q * (size_t)np + (size_t)p;
                double old_effective[3] = {
                    fields.v[qp][0] - fields.du_time[qp][0],
                    fields.v[qp][1] - fields.du_time[qp][1],
                    fields.v[qp][2] - fields.du_time[qp][2]
                };
                double grad_storage[3][3];
                double* grad_rows[3] = {
                    grad_storage[0],
                    grad_storage[1],
                    grad_storage[2]
                };
                double vec[4];
                int d;

                for(a = 0; a < 3; a++){
                    for(b = 0; b < 3; b++){
                        grad_storage[a][b] =
                            fields.grad_v[qp][a][b];
                    }
                }

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    system->basis.N[p][i],
                    system->fe.geo[elem][p].grad_N[i],
                    fields.v[qp],
                    old_effective,
                    grad_rows,
                    fields.pressure[qp],
                    fields.grad_pressure[qp],
                    system->vals.density,
                    system->vals.viscosity,
                    fields.tau[qp],
                    fields.tau_c[qp],
                    system->vals.dt,
                    fields.du_time[qp]);

                for(alpha = 0; alpha < nt; alpha++){
                    const double weight =
                        layout->time.weight[q]
                        * layout->time.psi[alpha][q];

                    for(d = 0; d < 4; d++){
                        FLUID_RES_IP(alpha, d, p) +=
                            weight * vec[d];
                    }
                }
            }
        }

        /*
         * Jacobian and dG upwind jump.  The previous-state derivative J_old
         * is used in two ways:
         *
         *  1. temporal chain rule in the slab interior,
         *  2. transient jump operator at tau=0.
         *
         * For a standard mass term J_old=-M/dt, so
         *
         *   -J_old (U^- - U_prev^+) = M/dt (U^- - U_prev^+).
         *
         * Here the same operation also retains SUPS/PSPG transient pieces.
         */
        for(j = 0; j < nl; j++){
            if(assemble_jacobian){
                memset(jac_ip, 0, jac_ip_size * sizeof(double));
                if(slab > 0){
                    memset(
                        prev_jac_ip,
                        0,
                        prev_jac_ip_size * sizeof(double));
                }
            }

            for(p = 0; p < np; p++){
                double grad_left_storage[3][3];
                double* grad_left_rows[3] = {
                    grad_left_storage[0],
                    grad_left_storage[1],
                    grad_left_storage[2]
                };
                double zero_du[3] = {0.0, 0.0, 0.0};
                double J_old_left[4][3];

                for(a = 0; a < 3; a++){
                    for(b = 0; b < 3; b++){
                        grad_left_storage[a][b] =
                            fields.grad_v_left[p][a][b];
                    }
                }

                fluid_dg_old_jacobian_ip(
                    &system->fe,
                    &system->basis,
                    &system->vals,
                    elem,
                    i,
                    j,
                    p,
                    fields.v_left[p],
                    fields.v_left[p],
                    grad_left_rows,
                    fields.pressure_left[p],
                    fields.grad_pressure_left[p],
                    fields.tau_left[p],
                    fields.tau_c_left[p],
                    zero_du,
                    system->vals.st_dg_jacobian_fd_epsilon,
                    J_old_left);

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 3; b++){
                        const double jump_nodal =
                            fields.local_v_left[j][b]
                            - fields.local_v_previous[j][b];
                        FLUID_JUMP_IP(a, p) -=
                            J_old_left[a][b] * jump_nodal;
                    }
                }

                if(assemble_jacobian){
                    /*
                     * Frozen-J_old linearization of the jump operator.  This
                     * is consistent with the modified-Newton treatment already
                     * used by the verified SUPS tangent for stabilization
                     * coefficients.
                     */
                    for(alpha = 0; alpha < nt; alpha++){
                        for(beta = 0; beta < nt; beta++){
                            for(a = 0; a < 4; a++){
                                for(b = 0; b < 3; b++){
                                    FLUID_JAC_IP(
                                        alpha, beta, a, b, p) +=
                                        -layout->time.e_minus[alpha]
                                        * layout->time.e_minus[beta]
                                        * J_old_left[a][b];

                                    if(slab > 0){
                                        FLUID_PREV_JAC_IP(
                                            alpha, beta, a, b, p) +=
                                            layout->time.e_minus[alpha]
                                            * layout->time.e_plus[beta]
                                            * J_old_left[a][b];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(assemble_jacobian){
                /*
                 * Interior temporal Jacobian.
                 *
                 * Existing current-state tangent direction:
                 *   (delta U, delta U_old, delta du) = (1, 0, 1)
                 *
                 * Existing old-state tangent direction:
                 *   (0, 1, -1)
                 *
                 * A temporal coefficient beta generates:
                 *   delta U     = psi_beta
                 *   delta du    = dpsi_beta/dtau
                 *   delta U_old = psi_beta - dpsi_beta/dtau
                 *
                 * so the exact directional combination of the verified
                 * tangents is
                 *
                 *   psi_beta J_current
                 *   + (psi_beta - dpsi_beta) J_old.
                 */
                for(q = 0; q < nq; q++){
                    const double tau_time =
                        layout->time.tau[q];

                    for(p = 0; p < np; p++){
                        const size_t qp =
                            (size_t)q * (size_t)np
                            + (size_t)p;
                        double old_effective[3] = {
                            fields.v[qp][0]
                                - fields.du_time[qp][0],
                            fields.v[qp][1]
                                - fields.du_time[qp][1],
                            fields.v[qp][2]
                                - fields.du_time[qp][2]
                        };
                        double grad_storage[3][3];
                        double* grad_rows[3] = {
                            grad_storage[0],
                            grad_storage[1],
                            grad_storage[2]
                        };
                        double J_current[4][4] = {{0.0}};
                        double J_old[4][3];

                        for(a = 0; a < 3; a++){
                            for(b = 0; b < 3; b++){
                                grad_storage[a][b] =
                                    fields.grad_v[qp][a][b];
                            }
                        }

                        BBFE_elemmat_fluid_sups_mat_NR(
                            J_current,
                            fields.J_inv[p],
                            system->basis.N[p][i],
                            system->basis.N[p][j],
                            system->fe.geo[elem][p].grad_N[i],
                            system->fe.geo[elem][p].grad_N[j],
                            fields.v[qp],
                            grad_rows,
                            fields.grad_pressure[qp],
                            system->vals.density,
                            system->vals.viscosity,
                            fields.tau[qp],
                            fields.tau_c[qp],
                            system->vals.dt,
                            fields.du_time[qp]);

                        fluid_dg_old_jacobian_ip(
                            &system->fe,
                            &system->basis,
                            &system->vals,
                            elem,
                            i,
                            j,
                            p,
                            fields.v[qp],
                            old_effective,
                            grad_rows,
                            fields.pressure[qp],
                            fields.grad_pressure[qp],
                            fields.tau[qp],
                            fields.tau_c[qp],
                            fields.du_time[qp],
                            system->vals.st_dg_jacobian_fd_epsilon,
                            J_old);

                        for(alpha = 0; alpha < nt; alpha++){
                            const double test_weight =
                                layout->time.weight[q]
                                * layout->time.psi[alpha][q];

                            for(beta = 0; beta < nt; beta++){
                                const double psi_beta =
                                    layout->time.psi[beta][q];
                                const double dpsi_beta =
                                    fluid_time_basis_derivative(
                                        &layout->time,
                                        beta,
                                        tau_time);

                                for(a = 0; a < 4; a++){
                                    for(b = 0; b < 4; b++){
                                        double value =
                                            psi_beta
                                            * J_current[a][b];

                                        if(b < 3){
                                            value +=
                                                (psi_beta - dpsi_beta)
                                                * J_old[a][b];
                                        }

                                        FLUID_JAC_IP(
                                            alpha,
                                            beta,
                                            a,
                                            b,
                                            p) +=
                                            test_weight * value;
                                    }
                                }
                            }
                        }
                    }
                }

                for(alpha = 0; alpha < nt; alpha++){
                    for(beta = 0; beta < nt; beta++){
                        for(a = 0; a < 4; a++){
                            for(b = 0; b < 4; b++){
                                const double value =
                                    BBFE_std_integ_calc(
                                        np,
                                        &FLUID_JAC_IP(
                                            alpha,
                                            beta,
                                            a,
                                            b,
                                            0),
                                        system->basis.integ_weight,
                                        fields.Jacobian_ip);
                                const int row_dof =
                                    ST_fom_block_dof(
                                        layout,
                                        slab,
                                        alpha,
                                        a);
                                const int col_dof =
                                    ST_fom_block_dof(
                                        layout,
                                        slab,
                                        beta,
                                        b);

                                monolis_add_scalar_to_sparse_matrix_R(
                                    &system->monolis_window,
                                    system->fe.conn[elem][i],
                                    system->fe.conn[elem][j],
                                    row_dof,
                                    col_dof,
                                    value);
                            }

                            if(slab > 0){
                                for(b = 0; b < 3; b++){
                                    const double value =
                                        BBFE_std_integ_calc(
                                            np,
                                            &FLUID_PREV_JAC_IP(
                                                alpha,
                                                beta,
                                                a,
                                                b,
                                                0),
                                            system->basis.integ_weight,
                                            fields.Jacobian_ip);
                                    const int row_dof =
                                        ST_fom_block_dof(
                                            layout,
                                            slab,
                                            alpha,
                                            a);
                                    const int col_dof =
                                        ST_fom_block_dof(
                                            layout,
                                            slab - 1,
                                            beta,
                                            b);

                                    monolis_add_scalar_to_sparse_matrix_R(
                                        &system->monolis_window,
                                        system->fe.conn[elem][i],
                                        system->fe.conn[elem][j],
                                        row_dof,
                                        col_dof,
                                        value);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(alpha = 0; alpha < nt; alpha++){
            for(a = 0; a < 4; a++){
                double value;

                for(p = 0; p < np; p++){
                    FLUID_RES_IP(alpha, a, p) +=
                        layout->time.e_minus[alpha]
                        * FLUID_JUMP_IP(a, p);
                }

                value = BBFE_std_integ_calc(
                    np,
                    &FLUID_RES_IP(alpha, a, 0),
                    system->basis.integ_weight,
                    fields.Jacobian_ip);

                system->monolis_window.mat.R.B[
                    system->window_dof
                        * system->fe.conn[elem][i]
                    + ST_fom_block_dof(
                        layout,
                        slab,
                        alpha,
                        a)] -= value;
            }
        }
    }

#undef FLUID_RES_IP
#undef FLUID_JUMP_IP
#undef FLUID_JAC_IP
#undef FLUID_PREV_JAC_IP

    status = 0;

cleanup:
    fluid_dg_element_fields_free(&fields);
    free(coefficient);
    free(previous_trace_local);
    free(residual_ip);
    free(jump_ip);
    free(jac_ip);
    free(prev_jac_ip);
    return status;
}

static int fluid_dg_assemble_window(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int assemble_jacobian)
{
    int elem;
    int slab;

    if(system == NULL || layout == NULL ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(elem = 0; elem < system->fe.total_num_elems; elem++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            if(fluid_dg_add_element_slab(
                system,
                layout,
                previous_trace,
                window_state,
                elem,
                slab,
                assemble_jacobian) != 0)
            {
                return 1;
            }
        }
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    return 0;
}

static void fluid_dg_constant_predictor(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state)
{
    int node;
    int slab;
    int alpha;
    int component;

    memset(
        window_state,
        0,
        ST_fom_window_vector_length(layout) * sizeof(double));

    for(node = 0; node < layout->num_space_points; node++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                for(component = 0;
                    component < FLUID_SUPS_PHYSICAL_DOF;
                    component++)
                {
                    window_state[
                        ST_fom_vector_index(
                            layout,
                            node,
                            slab,
                            alpha,
                            component)] =
                        previous_trace[
                            (size_t)node
                                * FLUID_SUPS_PHYSICAL_DOF
                            + (size_t)component];
                }
            }
        }
    }

    fluid_apply_physical_bc_to_window(
        system, layout, window_state);
    fluid_apply_pressure_gauge_to_window_state(
        system, layout, window_state);
}

static int fluid_sups_st_solve_window_all_at_once_high_order(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report)
{
    const int slabs = layout != NULL ? layout->num_slabs : 0;
    const int nstate = system != NULL
        ? system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF
        : 0;
    const size_t window_length =
        layout != NULL ? ST_fom_window_vector_length(layout) : 0;
    double* initial_slab_total = NULL;
    double* initial_slab_momentum = NULL;
    double* initial_slab_mass = NULL;
    double* final_slab_total = NULL;
    double* final_slab_momentum = NULL;
    double* final_slab_mass = NULL;
    double* final_trace = NULL;
    double initial_total = 0.0;
    double initial_momentum = 0.0;
    double initial_mass = 0.0;
    double final_total = INFINITY;
    double final_momentum = INFINITY;
    double final_mass = INFINITY;
    double relative_total = INFINITY;
    double relative_momentum = INFINITY;
    double relative_mass = INFINITY;
    const int linear_solver =
        system != NULL ? fluid_window_linear_solver_type()
                       : MONOLIS_ITER_BICGSAFE;
    const double linear_epsilon =
        system != NULL ? fluid_window_linear_epsilon(system)
                       : 1.0e-8;
    const int linear_max_iter =
        system != NULL ? fluid_window_linear_max_iter(system)
                       : 10000;
    int iteration;
    int slab;
    int converged = 0;
    int status = 1;

    if(system == NULL || layout == NULL ||
       previous_trace == NULL || window_state == NULL ||
       slabs <= 0 || nstate <= 0 ||
       layout->time.degree < 1 ||
       layout->time.degree > 2 ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       layout->num_owned_space_points !=
           system->num_owned_space_nodes ||
       system->window_dof !=
           ST_fom_window_dof_per_space_point(layout))
    {
        return 1;
    }

    initial_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_trace =
        (double*)calloc((size_t)nstate, sizeof(double));

    if(initial_slab_total == NULL ||
       initial_slab_momentum == NULL ||
       initial_slab_mass == NULL ||
       final_slab_total == NULL ||
       final_slab_momentum == NULL ||
       final_slab_mass == NULL ||
       final_trace == NULL)
    {
        goto cleanup;
    }

    if(window_report != NULL){
        memset(window_report, 0, sizeof(*window_report));
        window_report->relative_total_residual = INFINITY;
        window_report->relative_momentum_residual = INFINITY;
        window_report->relative_mass_residual = INFINITY;
    }
    if(slab_reports != NULL){
        memset(
            slab_reports,
            0,
            (size_t)slabs * sizeof(*slab_reports));
    }

    fluid_dg_constant_predictor(
        system,
        layout,
        previous_trace,
        window_state);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        window_state);

    for(iteration = 0;
        iteration < system->vals.nr_max_iter;
        iteration++)
    {
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\n%s ================ ALL-AT-ONCE dG(%d) window %d "
                "Newton %d ================\n",
                CODENAME,
                layout->time.degree,
                window_id,
                iteration + 1);
        }

        if(fluid_dg_assemble_window(
            system,
            layout,
            previous_trace,
            window_state,
            1) != 0)
        {
            goto cleanup;
        }

        if(iteration == 0){
            fluid_window_rhs_norms(
                system,
                layout,
                &initial_total,
                &initial_momentum,
                &initial_mass,
                initial_slab_total,
                initial_slab_momentum,
                initial_slab_mass);
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s dG(%d) window Newton linear solve: "
                "%d dof/node = %d physical x %d time dof x %d slabs; "
                "ONE MONOLIS solve; solver=%s eps=%.3e max_iter=%d\n",
                CODENAME,
                layout->time.degree,
                system->window_dof,
                FLUID_SUPS_PHYSICAL_DOF,
                layout->time.n_dof,
                slabs,
                fluid_window_linear_solver_name(linear_solver),
                linear_epsilon,
                linear_max_iter);
        }

        /*
         * Use a zero correction initial guess every Newton iteration.
         * Reusing the previous Newton correction is especially harmful when
         * the nonlinear residual has already dropped by several orders.
         */
        fluid_zero_window_linear_solution(system);

        BBFE_sys_monowrap_solve(
            &system->monolis_window,
            &system->mono_com,
            system->monolis_window.mat.R.X,
            linear_solver,
            MONOLIS_PREC_DIAG,
            linear_max_iter,
            linear_epsilon);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            system->monolis_window.mat.R.X);

        {
            size_t row;
            for(row = 0; row < window_length; row++){
                window_state[row] +=
                    system->monolis_window.mat.R.X[row];
            }
        }

        fluid_apply_physical_bc_to_window(
            system, layout, window_state);
        fluid_apply_pressure_gauge_to_window_state(
            system, layout, window_state);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            window_state);

        if(fluid_dg_assemble_window(
            system,
            layout,
            previous_trace,
            window_state,
            0) != 0)
        {
            goto cleanup;
        }

        fluid_window_rhs_norms(
            system,
            layout,
            &final_total,
            &final_momentum,
            &final_mass,
            final_slab_total,
            final_slab_momentum,
            final_slab_mass);

        relative_total =
            fluid_window_relative_norm(
                final_total, initial_total);
        relative_momentum =
            fluid_window_relative_norm(
                final_momentum, initial_momentum);
        relative_mass =
            fluid_window_relative_norm(
                final_mass, initial_mass);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s dG(%d) window %d Newton %d residual: "
                "total=%e momentum=%e mass=%e "
                "(abs total=%.15e initial=%.15e)\n",
                CODENAME,
                layout->time.degree,
                window_id,
                iteration + 1,
                relative_total,
                relative_momentum,
                relative_mass,
                final_total,
                initial_total);
        }

        if(isfinite(relative_total) &&
           relative_total < system->vals.nr_epsilon)
        {
            converged = 1;
            break;
        }
    }

    if(!converged){
        fprintf(
            stderr,
            "ERROR: dG(%d) all-at-once window Newton did not converge: "
            "window=%d rel_total=%.15e\n",
            layout->time.degree,
            window_id,
            relative_total);
        goto cleanup;
    }

    if(slab_reports != NULL){
        for(slab = 0; slab < slabs; slab++){
            slab_reports[slab].converged = 1;
            slab_reports[slab].nonlinear_iterations =
                iteration + 1;
            slab_reports[slab].relative_total_residual =
                fluid_window_relative_norm(
                    final_slab_total[slab],
                    initial_slab_total[slab]);
            slab_reports[slab].relative_momentum_residual =
                fluid_window_relative_norm(
                    final_slab_momentum[slab],
                    initial_slab_momentum[slab]);
            slab_reports[slab].relative_mass_residual =
                fluid_window_relative_norm(
                    final_slab_mass[slab],
                    initial_slab_mass[slab]);
        }
    }

    if(window_report != NULL){
        window_report->converged = 1;
        window_report->nonlinear_iterations =
            iteration + 1;
        window_report->relative_total_residual =
            relative_total;
        window_report->relative_momentum_residual =
            relative_momentum;
        window_report->relative_mass_residual =
            relative_mass;
    }

    if(ST_fom_extract_last_right_trace(
        layout,
        window_state,
        final_trace) != ST_FOM_SUCCESS)
    {
        goto cleanup;
    }

    if(fluid_sups_st_import_state(
        system,
        final_trace,
        1) != 0 ||
       fluid_sups_st_synchronize_state(system) != 0)
    {
        goto cleanup;
    }

    status = 0;

cleanup:
    free(initial_slab_total);
    free(initial_slab_momentum);
    free(initial_slab_mass);
    free(final_slab_total);
    free(final_slab_momentum);
    free(final_slab_mass);
    free(final_trace);
    return status;
}

int fluid_sups_st_solve_window_all_at_once_dg(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report)
{
    if(layout == NULL){
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_sups_st_solve_window_all_at_once_dg0(
            system,
            window_id,
            layout,
            previous_trace,
            window_state,
            slab_reports,
            window_report);
    }

    return fluid_sups_st_solve_window_all_at_once_high_order(
        system,
        window_id,
        layout,
        previous_trace,
        window_state,
        slab_reports,
        window_report);
}


