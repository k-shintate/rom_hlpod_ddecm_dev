static void fluid_set_element_jacobian_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    FLUID_SUPS_VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip;      double* Jacobian_ip;
    val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
    Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

    double**  val_ip_vec   = BB_std_calloc_2d_double(val_ip_vec , 4 , np);
    double*   integ_val_vec= BB_std_calloc_1d_double(integ_val_vec, 4);

    double** local_v       = BB_std_calloc_2d_double(local_v, nl, 3);
    double** v_ip          = BB_std_calloc_2d_double(v_ip   , np, 3);
    double*** grad_v_ip    = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

    double** local_v_old   = BB_std_calloc_2d_double(local_v_old, nl, 3);
    double** v_ip_old      = BB_std_calloc_2d_double(v_ip_old   , np, 3);

    double*  local_p       = BB_std_calloc_1d_double(local_p, nl);
    double*  p_ip          = BB_std_calloc_1d_double(p_ip  , np);
    double** grad_p_ip     = BB_std_calloc_2d_double(grad_p_ip, np, 3);

    double A[4][4];
    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);

        BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);
        
        for (int p = 0; p < np; ++p) {
            BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
                (void)h_e;
        
        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p)
                for (int d = 0; d < 4; ++d) val_ip_vec[d][p] = 0.0;

            for (int j = 0; j < nl; ++j) {
                for (int p = 0; p < np; ++p) {
                        double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };

                    BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                    const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                    const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                    BBFE_elemmat_fluid_sups_mat_NR(
                        A, J_inv,
                        basis->N[p][i], basis->N[p][j],
                        fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                        v_ip[p], grad_v_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            val_ip[a][b][p] = A[a][b];
                            A[a][b] = 0.0;
                        }
                    }
                }

                for (int a = 0; a < 4; ++a) {
                    for (int b = 0; b < 4; ++b) {
                        const double integ_val = BBFE_std_integ_calc(
                            np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
                    }
                }
            }
        }
    }

    BB_std_free_3d_double(val_ip , 4 , 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip   , np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old,   np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip   , np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
}

static void fluid_set_element_residuals_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    FLUID_SUPS_VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip;      double* Jacobian_ip;
    val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
    Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

    double**  val_ip_vec   = BB_std_calloc_2d_double(val_ip_vec , 4 , np);
    double*   integ_val_vec= BB_std_calloc_1d_double(integ_val_vec, 4);

    double** local_v       = BB_std_calloc_2d_double(local_v, nl, 3);
    double** v_ip          = BB_std_calloc_2d_double(v_ip   , np, 3);
    double*** grad_v_ip    = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

    double** local_v_old   = BB_std_calloc_2d_double(local_v_old, nl, 3);
    double** v_ip_old      = BB_std_calloc_2d_double(v_ip_old   , np, 3);

    double*  local_p       = BB_std_calloc_1d_double(local_p, nl);
    double*  p_ip          = BB_std_calloc_1d_double(p_ip  , np);
    double** grad_p_ip     = BB_std_calloc_2d_double(grad_p_ip, np, 3);

    double vec[4];
    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);

        BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);
        
        for (int p = 0; p < np; ++p) {
            BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }
        
        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p)
                for (int d = 0; d < 4; ++d) val_ip_vec[d][p] = 0.0;

            for (int p = 0; p < np; ++p) {
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density, vals->viscosity,
                    tau, tau_c, vals->dt,
                    du_time);

                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }

    BB_std_free_3d_double(val_ip , 4 , 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip   , np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old,   np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip   , np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
}

static void fluid_set_element_residuals_mass_con(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals)
{
    /*
     * Diagnostic decomposition only.
     *
     * The repository contains the validated full NR residual kernel
     * BBFE_elemmat_fluid_sups_vec_NR(), but no definition of the historical
     * helper BBFE_elemmat_fluid_sups_vec_NR_mass_con().
     *
     * Assemble the exact same full residual used by Newton.  The subsequent
     * BBFE_fluid_calc_internal_norm_mass() reads only the continuity/pressure
     * component, so the nonlinear solve algebra is unchanged.
     */
    fluid_set_element_residuals_NR(monolis, fe, basis, vals);
}

static void fluid_set_element_residuals_momentum_eq(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals)
{
    /*
     * Diagnostic decomposition only.
     *
     * Assemble the exact full residual used by Newton.  The subsequent
     * BBFE_fluid_calc_internal_norm_momentum() reads only the three momentum
     * components.  This does not alter the Newton Jacobian, Newton RHS,
     * update, or convergence decision based on the total residual.
     */
    fluid_set_element_residuals_NR(monolis, fe, basis, vals);
}


