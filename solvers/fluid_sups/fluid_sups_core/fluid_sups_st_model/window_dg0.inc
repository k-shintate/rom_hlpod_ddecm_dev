static void fluid_add_window_current_jacobian_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int offset = FLUID_SUPS_PHYSICAL_DOF * slab;

    double*** val_ip =
        BB_std_calloc_3d_double(NULL, 4, 4, np);
    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double A[4][4];
    double J_inv[3][3];
    int e;
    int i;
    int j;
    int p;
    int a;
    int b;

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(
                nl, local_p, basis->N[p]);
        }

        for(i = 0; i < nl; i++){
            for(j = 0; j < nl; j++){
                for(p = 0; p < np; p++){
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };
                    double tau;
                    double tau_c;

                    BB_calc_mat3d_inverse(
                        fe->geo[e][p].J,
                        fe->geo[e][p].Jacobian,
                        J_inv);

                    tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        vals->density,
                        vals->viscosity,
                        v_ip[p],
                        vals->dt);
                    tau_c =
                        BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                            J_inv,
                            fe->geo[e][p].Jacobian,
                            vals->density,
                            vals->viscosity,
                            v_ip[p],
                            vals->dt);

                    BBFE_elemmat_fluid_sups_mat_NR(
                        A,
                        J_inv,
                        basis->N[p][i],
                        basis->N[p][j],
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        v_ip[p],
                        grad_v_ip[p],
                        grad_p_ip[p],
                        vals->density,
                        vals->viscosity,
                        tau,
                        tau_c,
                        vals->dt,
                        du_time);

                    for(a = 0; a < 4; a++){
                        for(b = 0; b < 4; b++){
                            val_ip[a][b][p] = A[a][b];
                            A[a][b] = 0.0;
                        }
                    }
                }

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 4; b++){
                        const double value =
                            BBFE_std_integ_calc(
                                np,
                                val_ip[a][b],
                                basis->integ_weight,
                                Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            offset + a,
                            offset + b,
                            value);
                    }
                }
            }
        }
    }

    BB_std_free_3d_double(val_ip, 4, 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
}

static void fluid_add_window_residual_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab,
    int window_dof)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int offset = FLUID_SUPS_PHYSICAL_DOF * slab;

    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** val_ip_vec =
        BB_std_calloc_2d_double(NULL, 4, np);
    double* integ_val_vec =
        BB_std_calloc_1d_double(NULL, 4);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double vec[4];
    double J_inv[3][3];
    int e;
    int i;
    int p;
    int d;

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(
                nl, local_p, basis->N[p]);
        }

        for(i = 0; i < nl; i++){
            for(p = 0; p < np; p++){
                double du_time[3] = {
                    v_ip[p][0] - v_ip_old[p][0],
                    v_ip[p][1] - v_ip_old[p][1],
                    v_ip[p][2] - v_ip_old[p][2]
                };
                double tau;
                double tau_c;

                BB_calc_mat3d_inverse(
                    fe->geo[e][p].J,
                    fe->geo[e][p].Jacobian,
                    J_inv);

                tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);
                tau_c =
                    BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        vals->density,
                        vals->viscosity,
                        v_ip[p],
                        vals->dt);

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density,
                    vals->viscosity,
                    tau,
                    tau_c,
                    vals->dt,
                    du_time);

                for(d = 0; d < 4; d++){
                    val_ip_vec[d][p] = vec[d];
                }
            }

            for(d = 0; d < 4; d++){
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np,
                    val_ip_vec[d],
                    basis->integ_weight,
                    Jacobian_ip);

                monolis->mat.R.B[
                    window_dof * fe->conn[e][i]
                    + offset + d] -= integ_val_vec[d];
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
}

typedef enum {
    FLUID_JOLD_ANALYTIC = 0,
    FLUID_JOLD_VERIFY   = 1,
    FLUID_JOLD_FD       = 2
} FLUID_JOLD_MODE;

typedef struct {
    int initialized;
    FLUID_JOLD_MODE mode;
    int verify_samples;
    int verify_stride;
    long long seen_count;
    int verified_count;
    int pass_announced;
    double verify_tolerance;
} FLUID_JOLD_CONTROL;

static FLUID_JOLD_CONTROL fluid_jold_control = {0};

static FLUID_JOLD_MODE fluid_jold_mode_from_environment(void)
{
    const char* text = getenv("FLUID_ST_JOLD_MODE");

    if(text == NULL || text[0] == '\0' ||
       strcmp(text, "verify") == 0)
    {
        return FLUID_JOLD_VERIFY;
    }

    if(strcmp(text, "analytic") == 0){
        return FLUID_JOLD_ANALYTIC;
    }

    if(strcmp(text, "fd") == 0 ||
       strcmp(text, "finite_difference") == 0)
    {
        return FLUID_JOLD_FD;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_JOLD_MODE must be "
        "verify, analytic, or fd; got '%s'.\n",
        text);
    exit(EXIT_FAILURE);
}

static void fluid_jold_initialize_control(void)
{
    if(fluid_jold_control.initialized){
        return;
    }

    fluid_jold_control.mode =
        fluid_jold_mode_from_environment();

    fluid_jold_control.verify_samples =
        fluid_env_int("FLUID_ST_JOLD_VERIFY_SAMPLES", 64);

    fluid_jold_control.verify_stride =
        fluid_env_int("FLUID_ST_JOLD_VERIFY_STRIDE", 5000);

    fluid_jold_control.verify_tolerance =
        fluid_env_double("FLUID_ST_JOLD_VERIFY_TOL", 5.0e-6);

    if(fluid_jold_control.verify_samples <= 0){
        fluid_jold_control.verify_samples = 1;
    }
    if(fluid_jold_control.verify_stride <= 0){
        fluid_jold_control.verify_stride = 1;
    }
    if(fluid_jold_control.verify_tolerance <= 0.0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_JOLD_VERIFY_TOL must be positive.\n");
        exit(EXIT_FAILURE);
    }

    fluid_jold_control.seen_count = 0;
    fluid_jold_control.verified_count = 0;
    fluid_jold_control.pass_announced = 0;
    fluid_jold_control.initialized = 1;

    if(monolis_mpi_get_global_my_rank() == 0){
        const char* mode_name =
            fluid_jold_control.mode == FLUID_JOLD_ANALYTIC
                ? "analytic"
                : (fluid_jold_control.mode == FLUID_JOLD_FD
                    ? "fd"
                    : "verify");

        printf(
            "%s J_old element kernel mode: %s",
            CODENAME,
            mode_name);

        if(fluid_jold_control.mode == FLUID_JOLD_VERIFY){
            printf(
                " (samples=%d stride=%d tolerance=%.3e; "
                "FD is used for assembly)",
                fluid_jold_control.verify_samples,
                fluid_jold_control.verify_stride,
                fluid_jold_control.verify_tolerance);
        }
        printf("\n");
    }
}

static void fluid_old_jacobian_ip_fd_reference(
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
    double base[4];
    int a;
    int b;

    BBFE_elemmat_fluid_sups_vec_NR(
        base,
        basis->N[p][test_node],
        fe->geo[elem][p].grad_N[test_node],
        velocity,
        velocity_old,
        grad_velocity,
        pressure,
        grad_pressure,
        vals->density,
        vals->viscosity,
        tau,
        tau_c,
        vals->dt,
        du_time);

    for(b = 0; b < 3; b++){
        double old_perturbed[3] = {
            velocity_old[0],
            velocity_old[1],
            velocity_old[2]
        };
        double du_perturbed[3] = {
            du_time[0],
            du_time[1],
            du_time[2]
        };
        double perturbed[4];

        const double h =
            fd_epsilon
            * fmax(
                1.0,
                fmax(
                    fabs(velocity_old[b]),
                    fabs(velocity[b])));

        const double nodal_shape =
            basis->N[p][trial_node];

        old_perturbed[b] +=
            h * nodal_shape;
        du_perturbed[b] -=
            h * nodal_shape;

        BBFE_elemmat_fluid_sups_vec_NR(
            perturbed,
            basis->N[p][test_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            old_perturbed,
            grad_velocity,
            pressure,
            grad_pressure,
            vals->density,
            vals->viscosity,
            tau,
            tau_c,
            vals->dt,
            du_perturbed);

        for(a = 0; a < 4; a++){
            J_old[a][b] =
                (perturbed[a] - base[a]) / h;
        }
    }
}

static void fluid_print_jold_matrix(
    const char* label,
    const double J_old[4][3])
{
    int a;

    fprintf(stderr, "%s\n", label);
    for(a = 0; a < 4; a++){
        fprintf(
            stderr,
            "  [% .15e % .15e % .15e]\n",
            J_old[a][0],
            J_old[a][1],
            J_old[a][2]);
    }
}

static void fluid_verify_analytic_jold_sample(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    double tau,
    const double analytic[4][3],
    const double reference[4][3])
{
    double max_error = 0.0;
    double max_abs_difference = 0.0;
    int max_a = -1;
    int max_b = -1;
    int a;
    int b;

    if(fluid_jold_control.mode != FLUID_JOLD_VERIFY ||
       fluid_jold_control.verified_count >=
           fluid_jold_control.verify_samples)
    {
        return;
    }

    /*
     * Spread verification samples across the actual nonlinear run instead of
     * consuming all samples at the first element / first Newton state.
     */
    fluid_jold_control.seen_count++;

    if(((fluid_jold_control.seen_count - 1)
        % fluid_jold_control.verify_stride) != 0)
    {
        return;
    }

    for(a = 0; a < 4; a++){
        for(b = 0; b < 3; b++){
            const double difference =
                analytic[a][b] - reference[a][b];
            const double scale =
                fmax(
                    1.0,
                    fmax(
                        fabs(analytic[a][b]),
                        fabs(reference[a][b])));
            const double error =
                fabs(difference) / scale;

            if(error > max_error){
                max_error = error;
                max_abs_difference = fabs(difference);
                max_a = a;
                max_b = b;
            }
        }
    }

    if(!isfinite(max_error) ||
       max_error > fluid_jold_control.verify_tolerance)
    {
        fprintf(
            stderr,
            "\nERROR: analytic J_old does not match the linked "
            "BBFE_elemmat_fluid_sups_vec_NR provider.\n"
            "  element/test/trial/ip : %d / %d / %d / %d\n"
            "  sampled call          : %lld\n"
            "  max scaled error      : %.15e\n"
            "  max abs difference    : %.15e\n"
            "  offending entry       : (%d,%d)\n"
            "  N_i / N_j             : %.15e / %.15e\n"
            "  tau / dt              : %.15e / %.15e\n"
            "  density                : %.15e\n"
            "The Newton matrix is still using finite-difference J_old in "
            "verify mode, so use FLUID_ST_JOLD_MODE=fd to continue without "
            "analytic checks, or inspect the provider before enabling "
            "analytic mode.\n\n",
            elem,
            test_node,
            trial_node,
            p,
            fluid_jold_control.seen_count,
            max_error,
            max_abs_difference,
            max_a,
            max_b,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            tau,
            vals->dt,
            vals->density);

        fluid_print_jold_matrix(
            "analytic J_old:", analytic);
        fluid_print_jold_matrix(
            "finite-difference J_old:", reference);

        fflush(stderr);
        exit(EXIT_FAILURE);
    }

    fluid_jold_control.verified_count++;

    if(fluid_jold_control.verified_count >=
           fluid_jold_control.verify_samples &&
       !fluid_jold_control.pass_announced)
    {
        fluid_jold_control.pass_announced = 1;

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s analytic J_old verification PASS: "
                "%d sampled element kernels distributed over %lld calls, "
                "tolerance=%.3e. "
                "verify mode continued with FD assembly.\n",
                CODENAME,
                fluid_jold_control.verified_count,
                fluid_jold_control.seen_count,
                fluid_jold_control.verify_tolerance);
        }
    }
}

static void fluid_old_jacobian_ip(
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
    double J_old[4][3])
{
    fluid_jold_initialize_control();

    if(fluid_jold_control.mode == FLUID_JOLD_ANALYTIC){
        BBFE_elemmat_fluid_sups_mat_NR_old(
            J_old,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            vals->density,
            tau,
            vals->dt);
        return;
    }

    /*
     * Both `fd` and `verify` use the already validated finite-difference
     * tangent for the ACTUAL Newton matrix.  This guarantees that turning on
     * analytic verification cannot alter convergence.
     */
    fluid_old_jacobian_ip_fd_reference(
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
        vals->st_dg_jacobian_fd_epsilon,
        J_old);

    if(fluid_jold_control.mode == FLUID_JOLD_VERIFY){
        double analytic[4][3];

        BBFE_elemmat_fluid_sups_mat_NR_old(
            analytic,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            vals->density,
            tau,
            vals->dt);

        fluid_verify_analytic_jold_sample(
            fe,
            basis,
            vals,
            elem,
            test_node,
            trial_node,
            p,
            tau,
            analytic,
            J_old);
    }
}

static void fluid_add_window_previous_jacobian_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int row_offset =
        FLUID_SUPS_PHYSICAL_DOF * slab;
    const int col_offset =
        FLUID_SUPS_PHYSICAL_DOF * (slab - 1);

    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* tau =
        BB_std_calloc_1d_double(NULL, np);
    double* tau_c =
        BB_std_calloc_1d_double(NULL, np);
    double*** tangent_ip =
        BB_std_calloc_3d_double(NULL, 4, 3, np);

    int e;
    int i;
    int j;
    int p;
    int a;
    int b;

    if(slab <= 0){
        return;
    }

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            double J_inv[3][3];

            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] =
                BBFE_std_mapping_scalar(
                    nl, local_p, basis->N[p]);

            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            tau[p] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);

            tau_c[p] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);
        }

        for(i = 0; i < nl; i++){
            for(j = 0; j < nl; j++){
                for(p = 0; p < np; p++){
                    double J_old[4][3];
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };

                    fluid_old_jacobian_ip(
                        fe,
                        basis,
                        vals,
                        e,
                        i,
                        j,
                        p,
                        v_ip[p],
                        v_ip_old[p],
                        grad_v_ip[p],
                        p_ip[p],
                        grad_p_ip[p],
                        tau[p],
                        tau_c[p],
                        du_time,
                        J_old);

                    for(a = 0; a < 4; a++){
                        for(b = 0; b < 3; b++){
                            tangent_ip[a][b][p] =
                                J_old[a][b];
                        }
                    }
                }

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 3; b++){
                        const double value =
                            BBFE_std_integ_calc(
                                np,
                                tangent_ip[a][b],
                                basis->integ_weight,
                                Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            row_offset + a,
                            col_offset + b,
                            value);
                    }
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
    BB_std_free_1d_double(tau, np);
    BB_std_free_1d_double(tau_c, np);
    BB_std_free_3d_double(tangent_ip, 4, 3, np);
}

static int fluid_load_window_state_pair(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int slab,
    double* current_state,
    double* previous_state)
{
    if(fluid_unpack_window_slab_state(
        system, layout, window_state, slab, current_state) != 0)
    {
        return 1;
    }

    if(slab == 0){
        memcpy(
            previous_state,
            previous_trace,
            (size_t)system->fe.total_num_nodes
                * FLUID_SUPS_PHYSICAL_DOF
                * sizeof(double));
    }
    else{
        if(fluid_unpack_window_slab_state(
            system,
            layout,
            window_state,
            slab - 1,
            previous_state) != 0)
        {
            return 1;
        }
    }

    fluid_load_state_pair(
        system,
        current_state,
        previous_state);
    return 0;
}

static int fluid_assemble_window_residual_only(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    const int nstate =
        system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF;
    double* current_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    double* previous_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    int slab;

    if(current_state == NULL || previous_state == NULL){
        free(current_state);
        free(previous_state);
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(slab = 0; slab < layout->num_slabs; slab++){
        if(fluid_load_window_state_pair(
            system,
            layout,
            previous_trace,
            window_state,
            slab,
            current_state,
            previous_state) != 0)
        {
            free(current_state);
            free(previous_state);
            return 1;
        }

        fluid_add_window_residual_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab,
            system->window_dof);
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    free(current_state);
    free(previous_state);
    return 0;
}

static int fluid_assemble_window_newton_system(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    const int nstate =
        system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF;
    double* current_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    double* previous_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    int slab;

    if(current_state == NULL || previous_state == NULL){
        free(current_state);
        free(previous_state);
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(slab = 0; slab < layout->num_slabs; slab++){
        if(fluid_load_window_state_pair(
            system,
            layout,
            previous_trace,
            window_state,
            slab,
            current_state,
            previous_state) != 0)
        {
            free(current_state);
            free(previous_state);
            return 1;
        }

        fluid_add_window_current_jacobian_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab);

        fluid_add_window_residual_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab,
            system->window_dof);

        if(slab > 0){
            fluid_add_window_previous_jacobian_block(
                &system->monolis_window,
                &system->fe,
                &system->basis,
                &system->vals,
                slab);
        }
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    free(current_state);
    free(previous_state);
    return 0;
}

static void fluid_window_rhs_norms(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* total_norm,
    double* momentum_norm,
    double* mass_norm,
    double* slab_total_norm,
    double* slab_momentum_norm,
    double* slab_mass_norm)
{
    const int slabs = layout->num_slabs;
    const int owned = system->num_owned_space_nodes;
    double* reduced =
        (double*)calloc((size_t)3 * (size_t)slabs, sizeof(double));
    int node;
    int slab;
    int alpha;
    int component;

    for(node = 0; node < owned; node++){
        for(slab = 0; slab < slabs; slab++){
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                for(component = 0;
                    component < FLUID_SUPS_PHYSICAL_DOF;
                    component++)
                {
                    const int block_dof = ST_fom_block_dof(
                        layout, slab, alpha, component);
                    const int index =
                        system->window_dof * node + block_dof;
                    const double value =
                        system->monolis_window.mat.R.B[index];

                    reduced[3 * slab + 0] += value * value;
                    if(component < 3){
                        reduced[3 * slab + 1] += value * value;
                    }
                    else{
                        reduced[3 * slab + 2] += value * value;
                    }
                }
            }
        }
    }

    monolis_allreduce_R(
        3 * slabs,
        reduced,
        MONOLIS_MPI_SUM,
        system->mono_com.comm);

    *total_norm = 0.0;
    *momentum_norm = 0.0;
    *mass_norm = 0.0;

    for(slab = 0; slab < slabs; slab++){
        const double total2 =
            fmax(reduced[3 * slab + 0], 0.0);
        const double momentum2 =
            fmax(reduced[3 * slab + 1], 0.0);
        const double mass2 =
            fmax(reduced[3 * slab + 2], 0.0);

        *total_norm += total2;
        *momentum_norm += momentum2;
        *mass_norm += mass2;

        if(slab_total_norm != NULL){
            slab_total_norm[slab] = sqrt(total2);
        }
        if(slab_momentum_norm != NULL){
            slab_momentum_norm[slab] = sqrt(momentum2);
        }
        if(slab_mass_norm != NULL){
            slab_mass_norm[slab] = sqrt(mass2);
        }
    }

    *total_norm = sqrt(fmax(*total_norm, 0.0));
    *momentum_norm = sqrt(fmax(*momentum_norm, 0.0));
    *mass_norm = sqrt(fmax(*mass_norm, 0.0));
    free(reduced);
}

static double fluid_window_relative_norm(
    double value,
    double initial)
{
    return initial > 1.0e-300 ? value / initial : value;
}


static int fluid_window_linear_solver_type(void)
{
    const char* value = getenv("FLUID_ST_WINDOW_LINEAR_SOLVER");

    if(value == NULL || value[0] == '\0' ||
       strcmp(value, "bicgsafe") == 0)
    {
        return MONOLIS_ITER_BICGSAFE;
    }

    if(strcmp(value, "bicgstab") == 0){
        return MONOLIS_ITER_BICGSTAB;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_WINDOW_LINEAR_SOLVER must be "
        "bicgsafe or bicgstab; got '%s'.\n",
        value);
    exit(EXIT_FAILURE);
}

static const char* fluid_window_linear_solver_name(int solver)
{
    if(solver == MONOLIS_ITER_BICGSTAB){
        return "BiCGSTAB";
    }
    return "BiCGSAFE";
}

static double fluid_window_linear_epsilon(
    const FLUID_SUPS_SYSTEM* system)
{
    const double value = fluid_env_double(
        "FLUID_ST_WINDOW_LINEAR_EPSILON",
        system->vals.mat_epsilon);

    if(value <= 0.0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_WINDOW_LINEAR_EPSILON must be positive.\n");
        exit(EXIT_FAILURE);
    }

    return value;
}

static int fluid_window_linear_max_iter(
    const FLUID_SUPS_SYSTEM* system)
{
    const int value = fluid_env_int(
        "FLUID_ST_WINDOW_LINEAR_MAX_ITER",
        system->vals.mat_max_iter);

    if(value <= 0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_WINDOW_LINEAR_MAX_ITER must be positive.\n");
        exit(EXIT_FAILURE);
    }

    return value;
}

static void fluid_zero_window_linear_solution(
    FLUID_SUPS_SYSTEM* system)
{
    const size_t length =
        (size_t)system->fe.total_num_nodes
        * (size_t)system->window_dof;

    if(system != NULL &&
       system->monolis_window.mat.R.X != NULL)
    {
        memset(
            system->monolis_window.mat.R.X,
            0,
            length * sizeof(double));
    }
}

int fluid_sups_st_solve_window_all_at_once_dg0(
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
    double* final_state = NULL;
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
       layout->time.degree != 0 ||
       layout->time.n_dof != 1 ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       layout->num_owned_space_points !=
           system->num_owned_space_nodes ||
       system->window_dof !=
           FLUID_SUPS_PHYSICAL_DOF * slabs)
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
    final_state =
        (double*)calloc((size_t)nstate, sizeof(double));

    if(initial_slab_total == NULL ||
       initial_slab_momentum == NULL ||
       initial_slab_mass == NULL ||
       final_slab_total == NULL ||
       final_slab_momentum == NULL ||
       final_slab_mass == NULL ||
       final_state == NULL)
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

    /*
     * Constant predictor across the whole window.  No slab is separately
     * solved before the window Newton iteration starts.
     */
    memset(window_state, 0, window_length * sizeof(double));
    for(slab = 0; slab < slabs; slab++){
        if(ST_fom_pack_dg0_slab(
            layout,
            slab,
            previous_trace,
            window_state) != ST_FOM_SUCCESS)
        {
            goto cleanup;
        }
    }
    fluid_apply_physical_bc_to_window(
        system, layout, window_state);

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
                "\n%s ================ ALL-AT-ONCE ST window %d "
                "Newton %d ================\n",
                CODENAME,
                window_id,
                iteration + 1);
        }

        /*
         * Assemble the COMPLETE window Jacobian and residual, then call
         * MONOLIS exactly once for this Newton iteration.
         */
        if(fluid_assemble_window_newton_system(
            system,
            layout,
            previous_trace,
            window_state) != 0)
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
                "%s window Newton linear solve: "
                "%d dof/node = %d physical dof x %d slabs; "
                "ONE MONOLIS solve; solver=%s eps=%.3e max_iter=%d\n",
                CODENAME,
                system->window_dof,
                FLUID_SUPS_PHYSICAL_DOF,
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

        if(fluid_assemble_window_residual_only(
            system,
            layout,
            previous_trace,
            window_state) != 0)
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
                "%s ALL-AT-ONCE window %d Newton %d residual: "
                "total=%e momentum=%e mass=%e\n",
                CODENAME,
                window_id,
                iteration + 1,
                relative_total,
                relative_momentum,
                relative_mass);
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
            "ERROR: all-at-once window Newton did not converge: "
            "window=%d rel_total=%.15e\n",
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

    if(ST_fom_unpack_dg0_slab(
        layout,
        slabs - 1,
        window_state,
        final_state) != ST_FOM_SUCCESS)
    {
        goto cleanup;
    }

    if(fluid_sups_st_import_state(
        system,
        final_state,
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
    free(final_state);
    return status;
}


