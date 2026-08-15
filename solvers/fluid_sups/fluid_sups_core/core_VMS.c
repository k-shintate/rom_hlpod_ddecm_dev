
static void BBFE_vms_strong_residual_qs(
    double       rM[3],
    double*      rC,
    const double u[3],
    const double u_old[3],
    double**     grad_u,
    const double grad_p[3],
    const double rho,
    const double dt)
{
    const double eps = 1.0e-30;

    double adv[3];
    for (int d = 0; d < 3; ++d) {
        adv[d] =
            u[0] * grad_u[d][0] +
            u[1] * grad_u[d][1] +
            u[2] * grad_u[d][2];
    }

    const double inv_dt = 1.0 / fmax(dt, eps);

    for (int d = 0; d < 3; ++d) {
        rM[d] =
            (u[d] - u_old[d]) * inv_dt
            + adv[d]
            + grad_p[d] / rho;
    }

    *rC = grad_u[0][0] + grad_u[1][1] + grad_u[2][2];
}


static double BBFE_vms_nu_quasistatic_residual(
    const double tau_m,
    const double h_e,
    const double rM[3],
    const double u[3],
    const double C_vms,
    const double cap_coeff)
{
    const double eps = 1.0e-30;

    const double r_norm = sqrt(
        rM[0]*rM[0] +
        rM[1]*rM[1] +
        rM[2]*rM[2] + eps);

    /*
      quasi-static subscale:
      u' = - tau_m rM
    */
    const double u_prime_norm = tau_m * r_norm;

    double nu_vms = C_vms * h_e * u_prime_norm;

    /*
      過大散逸を防ぐ cap
    */
    const double u_norm = sqrt(
        u[0]*u[0] +
        u[1]*u[1] +
        u[2]*u[2] + eps);

    const double nu_cap = cap_coeff * h_e * u_norm;

    if (nu_vms > nu_cap) nu_vms = nu_cap;
    if (nu_vms < 0.0)    nu_vms = 0.0;

    return nu_vms;
}


void BBFE_vms_mu_eff_tau(
    double*      mu_eff,
    double*      tau,
    double*      tau_c,
    const double J_inv[3][3],
    const double Jacobian,
    const double h_e,
    const double u[3],
    const double u_old[3],
    double**     grad_u,
    const double grad_p[3],
    const double rho,
    const double mu,
    const double dt,
    const double C_vms,
    const double cap_coeff)
{
    /*
      まず分子粘性で tau0 を作る
    */
    const double tau0 =
        BBFE_elemmat_fluid_sups_coef_metric_tensor(
            J_inv, Jacobian, rho, mu, u, dt);

    double rM[3];
    double rC;

    BBFE_vms_strong_residual_qs(
        rM, &rC, u, u_old, grad_u, grad_p, rho, dt);

    const double nu_vms =
        BBFE_vms_nu_quasistatic_residual(
            tau0, h_e, rM, u, C_vms, cap_coeff);

    *mu_eff = mu + rho * nu_vms;

    /*
      有効粘性で tau, tau_c を再計算する
    */
    *tau =
        BBFE_elemmat_fluid_sups_coef_metric_tensor(
            J_inv, Jacobian, rho, *mu_eff, u, dt);

    *tau_c =
        BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
            J_inv, Jacobian, rho, *mu_eff, u, dt);
}

static void BBFE_vms_mu_eff_derivative_qs(
    double       dmu_duj[3],
    double*      dmu_dpj,
    const double J_inv[3][3],
    const double Jacobian,
    const double h_e,
    const double u[3],
    const double u_old[3],
    double**     grad_u,
    const double grad_p[3],
    const double rho,
    const double mu_mol,
    const double dt,
    const double C_vms,
    const double cap_coeff,
    const double N_j,
    const double grad_N_j[3])
{
    const double eps = 1.0e-30;

    for (int k = 0; k < 3; ++k) {
        dmu_duj[k] = 0.0;
    }
    *dmu_dpj = 0.0;

    if (rho <= eps || dt <= eps || h_e <= eps || C_vms <= 0.0) {
        return;
    }

    double rM[3];
    double rC;

    BBFE_vms_strong_residual_qs(
        rM,
        &rC,
        u,
        u_old,
        grad_u,
        grad_p,
        rho,
        dt);

    const double r_norm = sqrt(
        rM[0]*rM[0] +
        rM[1]*rM[1] +
        rM[2]*rM[2] + eps);

    const double u_norm = sqrt(
        u[0]*u[0] +
        u[1]*u[1] +
        u[2]*u[2] + eps);

    /*
      tau0 は分子粘性 mu_mol で作る。
      これは BBFE_vms_mu_eff_tau() と合わせる。
    */
    const double tau0 =
        BBFE_elemmat_fluid_sups_coef_metric_tensor(
            J_inv,
            Jacobian,
            rho,
            mu_mol,
            u,
            dt);

    /*
      tau0 の速度微分。
      圧力には直接依存しない。
    */
    double dtau0_duj[3];
    double dummy_dtauc[3];

    BBFE_elemmat_fluid_sups_coef_metric_tensor_derivative(
        dtau0_duj,
        dummy_dtauc,
        J_inv,
        Jacobian,
        rho,
        mu_mol,
        u,
        dt,
        N_j);

    const double nu_raw = C_vms * h_e * tau0 * r_norm;
    const double nu_cap = cap_coeff * h_e * u_norm;

    const int use_cap = (nu_raw > nu_cap);

    const double inv_dt = 1.0 / fmax(dt, eps);

    const double u_dot_gradNj =
        u[0]*grad_N_j[0] +
        u[1]*grad_N_j[1] +
        u[2]*grad_N_j[2];

    /*
      velocity columns
    */
    for (int k = 0; k < 3; ++k) {
        double dnu = 0.0;

        if (!use_cap) {
            double drM[3];

            for (int a = 0; a < 3; ++a) {
                drM[a] =
                    ((a == k) ? N_j * inv_dt : 0.0)
                  + N_j * grad_u[a][k]
                  + ((a == k) ? u_dot_gradNj : 0.0);
            }

            const double dr_norm =
                (rM[0]*drM[0] +
                 rM[1]*drM[1] +
                 rM[2]*drM[2]) / r_norm;

            dnu =
                C_vms * h_e *
                (dtau0_duj[k] * r_norm + tau0 * dr_norm);
        }
        else {
            /*
              cap branch:
              nu_vms = cap_coeff * h_e * |u|
            */
            dnu = cap_coeff * h_e * N_j * u[k] / u_norm;
        }

        dmu_duj[k] = rho * dnu;
    }

    /*
      pressure column
    */
    if (!use_cap) {
        double drM_dp[3] = {
            grad_N_j[0] / rho,
            grad_N_j[1] / rho,
            grad_N_j[2] / rho
        };

        const double dr_norm_dp =
            (rM[0]*drM_dp[0] +
             rM[1]*drM_dp[1] +
             rM[2]*drM_dp[2]) / r_norm;

        const double dnu_dp =
            C_vms * h_e * tau0 * dr_norm_dp;

        *dmu_dpj = rho * dnu_dp;
    }
    else {
        *dmu_dpj = 0.0;
    }
}


static void BBFE_sups_tau_mu_derivative(
    double*      dtau_dmu,
    double*      dtauc_dmu,
    const double J_inv[3][3],
    const double rho,
    const double mu_eff,
    const double tau)
{
    const double eps = 1.0e-30;

    *dtau_dmu  = 0.0;
    *dtauc_dmu = 0.0;

    if (rho <= eps || tau <= eps) {
        return;
    }

    double G[3][3];
    BBFE_metric_tensor_G(G, J_inv);

    const double GG =
        BBFE_metric_tensor_G_colon_G(G);

    const double trG =
        BBFE_metric_tensor_trace(G);

    const double nu = mu_eff / rho;

    /*
      tau = denom^{-1/2}
      denom includes 36 nu^2 G:G
    */
    *dtau_dmu =
        -36.0 * tau * tau * tau * nu * GG / rho;

    /*
      tau_c = 1 / (tau * trG) を仮定。
      もし BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC()
      の定義が違うなら、ここは合わせて変更する。
    */
    if (fabs(trG) > eps) {
        *dtauc_dmu =
            -(*dtau_dmu) / (tau * tau * trG);
    }
}

static void BBFE_elemmat_fluid_mat_vms_mu_chain(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   Jacobian,
    const double   h_e,
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   u[3],
    const double   u_old[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   rho,
    const double   mu_mol,
    const double   mu_eff,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   C_vms,
    const double   cap_coeff,
    const double   du_time[3])
{
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            mat[a][b] = 0.0;
        }
    }

    double dmu_duj[3];
    double dmu_dpj;

    BBFE_vms_mu_eff_derivative_qs(
        dmu_duj,
        &dmu_dpj,
        J_inv,
        Jacobian,
        h_e,
        u,
        u_old,
        grad_u,
        grad_p,
        rho,
        mu_mol,
        dt,
        C_vms,
        cap_coeff,
        N_j,
        grad_N_j);

    double dmu_col[4] = {
        dmu_duj[0],
        dmu_duj[1],
        dmu_duj[2],
        dmu_dpj
    };

    double dtau_dmu;
    double dtauc_dmu;

    BBFE_sups_tau_mu_derivative(
        &dtau_dmu,
        &dtauc_dmu,
        J_inv,
        rho,
        mu_eff,
        tau);

    double adv[3];
    for (int d = 0; d < 3; ++d) {
        adv[d] =
            u[0] * grad_u[d][0] +
            u[1] * grad_u[d][1] +
            u[2] * grad_u[d][2];
    }

    const double div_u =
        grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    const double vdotGradNi =
        u[0]*grad_N_i[0] +
        u[1]*grad_N_i[1] +
        u[2]*grad_N_i[2];

    const double gradN_dot_gradp =
        grad_N_i[0]*grad_p[0] +
        grad_N_i[1]*grad_p[1] +
        grad_N_i[2]*grad_p[2];

    const double gradN_dot_adv =
        grad_N_i[0]*adv[0] +
        grad_N_i[1]*adv[1] +
        grad_N_i[2]*adv[2];

    const double gradN_dot_dutime =
        grad_N_i[0]*du_time[0] +
        grad_N_i[1]*du_time[1] +
        grad_N_i[2]*du_time[2];

    /*
      dR_visc / dmu_eff

      R_visc_i[d] =
        dt * mu_eff *
        {
          grad_N_i · grad_u[d]
          + grad_N_i[d] * div(u)
        }
    */
    double visc_kernel[3];

    for (int d = 0; d < 3; ++d) {
        visc_kernel[d] =
            grad_N_i[0]*grad_u[d][0] +
            grad_N_i[1]*grad_u[d][1] +
            grad_N_i[2]*grad_u[d][2]
          + grad_N_i[d] * div_u;
    }

    /*
      tau chain kernels.
      既存の BBFE_elemmat_fluid_mat_rom_nonlinear() の
      section (5), (6) と同じ residual kernel を使う。
    */
    double Rm_tau[3];

    for (int d = 0; d < 3; ++d) {
        Rm_tau[d] =
            dt * rho * vdotGradNi * adv[d]
          + dt       * vdotGradNi * grad_p[d]
          + rho      * vdotGradNi * du_time[d];
    }

    const double Rc_tau =
        dt * gradN_dot_adv
      + dt * gradN_dot_gradp / rho
      +      gradN_dot_dutime;

    for (int c = 0; c < 4; ++c) {
        const double dmu_c   = dmu_col[c];
        const double dtau_c  = dtau_dmu  * dmu_c;
        const double dtauc_c = dtauc_dmu * dmu_c;

        /*
          1) viscosity chain
        */
        for (int d = 0; d < 3; ++d) {
            mat[d][c] += dt * dmu_c * visc_kernel[d];
        }

        /*
          2) tau chain through mu_eff
        */
        for (int d = 0; d < 3; ++d) {
            mat[d][c] += dtau_c * Rm_tau[d];
        }

        mat[3][c] += dtau_c * Rc_tau;

        /*
          3) tau_c chain through mu_eff
        */
        for (int d = 0; d < 3; ++d) {
            mat[d][c] +=
                dt * rho * dtauc_c * div_u * grad_N_i[d];
        }
    }
}


void set_element_mat_NR_linear_VMS(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip      = BB_std_calloc_3d_double(NULL, 4, 4, np);
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double** val_ip_vec    = BB_std_calloc_2d_double(NULL, 4, np);
    double*  integ_val_vec = BB_std_calloc_1d_double(NULL, 4);

    double**  local_v   = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip      = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double** local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*  local_p   = BB_std_calloc_1d_double(NULL, nl);
    double*  p_ip      = BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip = BB_std_calloc_2d_double(NULL, np, 3);

    double* mu_eff_ip = BB_std_calloc_1d_double(NULL, np);
    double* tau_ip    = BB_std_calloc_1d_double(NULL, np);
    double* tau_c_ip  = BB_std_calloc_1d_double(NULL, np);

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

            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p],
                nl,
                local_v,
                fe->geo[e][p].grad_N);

            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p],
                nl,
                local_p,
                fe->geo[e][p].grad_N);

            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }

        const double vol = BBFE_std_integ_calc_volume(
            np,
            basis->integ_weight,
            Jacobian_ip);

        const double h_e = cbrt(vol);

        /*
          quasi-static residual-based VMS:
          各積分点で mu_eff, tau, tau_c を事前計算する。
        */
        for (int p = 0; p < np; ++p) {
            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            BBFE_vms_mu_eff_tau(
                &mu_eff_ip[p],
                &tau_ip[p],
                &tau_c_ip[p],
                J_inv,
                fe->geo[e][p].Jacobian,
                h_e,
                v_ip[p],
                v_ip_old[p],
                grad_v_ip[p],
                grad_p_ip[p],
                vals->density,
                vals->viscosity,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);
        }

        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p) {
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = 0.0;
                }
            }

            for (int j = 0; j < nl; ++j) {
                for (int p = 0; p < np; ++p) {
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };

                    BB_calc_mat3d_inverse(
                        fe->geo[e][p].J,
                        fe->geo[e][p].Jacobian,
                        J_inv);

                    const double mu_eff = mu_eff_ip[p];
                    const double tau    = tau_ip[p];
                    const double tau_c  = tau_c_ip[p];

                    BBFE_elemmat_fluid_mat_rom_linear(
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
                        mu_eff,
                        tau,
                        tau_c,
                        vals->dt,
                        du_time);

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
                            np,
                            val_ip[a][b],
                            basis->integ_weight,
                            Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            a,
                            b,
                            integ_val);
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

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);

    BB_std_free_1d_double(mu_eff_ip, np);
    BB_std_free_1d_double(tau_ip, np);
    BB_std_free_1d_double(tau_c_ip, np);
}


void set_element_vec_NR_linear_VMS(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip      = BB_std_calloc_3d_double(NULL, 4, 4, np);
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double** val_ip_vec    = BB_std_calloc_2d_double(NULL, 4, np);
    double*  integ_val_vec = BB_std_calloc_1d_double(NULL, 4);

    double**  local_v   = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip      = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double** local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*  local_p   = BB_std_calloc_1d_double(NULL, nl);
    double*  p_ip      = BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip = BB_std_calloc_2d_double(NULL, np, 3);

    double* mu_eff_ip = BB_std_calloc_1d_double(NULL, np);
    double* tau_ip    = BB_std_calloc_1d_double(NULL, np);
    double* tau_c_ip  = BB_std_calloc_1d_double(NULL, np);

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

            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p],
                nl,
                local_v,
                fe->geo[e][p].grad_N);

            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p],
                nl,
                local_p,
                fe->geo[e][p].grad_N);

            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }

        const double vol = BBFE_std_integ_calc_volume(
            np,
            basis->integ_weight,
            Jacobian_ip);

        const double h_e = cbrt(vol);

        for (int p = 0; p < np; ++p) {
            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            BBFE_vms_mu_eff_tau(
                &mu_eff_ip[p],
                &tau_ip[p],
                &tau_c_ip[p],
                J_inv,
                fe->geo[e][p].Jacobian,
                h_e,
                v_ip[p],
                v_ip_old[p],
                grad_v_ip[p],
                grad_p_ip[p],
                vals->density,
                vals->viscosity,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);
        }

        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p) {
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = 0.0;
                }
            }

            for (int p = 0; p < np; ++p) {
                double du_time[3] = {
                    v_ip[p][0] - v_ip_old[p][0],
                    v_ip[p][1] - v_ip_old[p][1],
                    v_ip[p][2] - v_ip_old[p][2]
                };

                BB_calc_mat3d_inverse(
                    fe->geo[e][p].J,
                    fe->geo[e][p].Jacobian,
                    J_inv);

                const double mu_eff = mu_eff_ip[p];
                const double tau    = tau_ip[p];
                const double tau_c  = tau_c_ip[p];

                BBFE_elemmat_fluid_vec_rom_linear(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density,
                    mu_eff,
                    tau,
                    tau_c,
                    vals->dt,
                    du_time);

                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                    vec[d] = 0.0;
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np,
                    val_ip_vec[d],
                    basis->integ_weight,
                    Jacobian_ip);

                monolis->mat.R.B[4 * fe->conn[e][i] + d] -= integ_val_vec[d];
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

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);

    BB_std_free_1d_double(mu_eff_ip, np);
    BB_std_free_1d_double(tau_ip, np);
    BB_std_free_1d_double(tau_c_ip, np);
}


void set_element_mat_NR_nonlinear_VMS(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip      = BB_std_calloc_3d_double(NULL, 4, 4, np);
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double** val_ip_vec    = BB_std_calloc_2d_double(NULL, 4, np);
    double*  integ_val_vec = BB_std_calloc_1d_double(NULL, 4);

    double**  local_v   = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip      = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double** local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*  local_p   = BB_std_calloc_1d_double(NULL, nl);
    double*  p_ip      = BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip = BB_std_calloc_2d_double(NULL, np, 3);

    double* mu_eff_ip = BB_std_calloc_1d_double(NULL, np);
    double* tau_ip    = BB_std_calloc_1d_double(NULL, np);
    double* tau_c_ip  = BB_std_calloc_1d_double(NULL, np);

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

            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p],
                nl,
                local_v,
                fe->geo[e][p].grad_N);

            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p],
                nl,
                local_p,
                fe->geo[e][p].grad_N);

            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }

        const double vol = BBFE_std_integ_calc_volume(
            np,
            basis->integ_weight,
            Jacobian_ip);

        const double h_e = cbrt(vol);

        for (int p = 0; p < np; ++p) {
            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            BBFE_vms_mu_eff_tau(
                &mu_eff_ip[p],
                &tau_ip[p],
                &tau_c_ip[p],
                J_inv,
                fe->geo[e][p].Jacobian,
                h_e,
                v_ip[p],
                v_ip_old[p],
                grad_v_ip[p],
                grad_p_ip[p],
                vals->density,
                vals->viscosity,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);
        }

        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p) {
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = 0.0;
                }
            }

            for (int j = 0; j < nl; ++j) {
                for (int p = 0; p < np; ++p) {
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };

                    BB_calc_mat3d_inverse(
                        fe->geo[e][p].J,
                        fe->geo[e][p].Jacobian,
                        J_inv);

                    const double mu_eff = mu_eff_ip[p];
                    const double tau    = tau_ip[p];
                    const double tau_c  = tau_c_ip[p];

                    BBFE_elemmat_fluid_mat_rom_nonlinear(
                        A,
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        basis->N[p][i],
                        basis->N[p][j],
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        v_ip[p],
                        grad_v_ip[p],
                        grad_p_ip[p],
                        vals->density,
                        mu_eff,
                        tau,
                        tau_c,
                        vals->dt,
                        du_time);

                    double A_vms_chain[4][4];

                    BBFE_elemmat_fluid_mat_vms_mu_chain(
                        A_vms_chain,
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        h_e,
                        basis->N[p][i],
                        basis->N[p][j],
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        v_ip[p],
                        v_ip_old[p],
                        grad_v_ip[p],
                        grad_p_ip[p],
                        vals->density,
                        vals->viscosity,
                        mu_eff,
                        tau,
                        tau_c,
                        vals->dt,
                        vals->C_vms,
                        vals->vms_cap_coeff,
                        du_time);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            val_ip[a][b][p] = A[a][b] + A_vms_chain[a][b];
                            A[a][b] = 0.0;
                            A_vms_chain[a][b] = 0.0;
                        }
                    }
                }

                for (int a = 0; a < 4; ++a) {
                    for (int b = 0; b < 4; ++b) {
                        const double integ_val = BBFE_std_integ_calc(
                            np,
                            val_ip[a][b],
                            basis->integ_weight,
                            Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            a,
                            b,
                            integ_val);
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

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);

    BB_std_free_1d_double(mu_eff_ip, np);
    BB_std_free_1d_double(tau_ip, np);
    BB_std_free_1d_double(tau_c_ip, np);
}


void set_element_vec_NR_nonlinear_VMS(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip      = BB_std_calloc_3d_double(NULL, 4, 4, np);
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double** val_ip_vec    = BB_std_calloc_2d_double(NULL, 4, np);
    double*  integ_val_vec = BB_std_calloc_1d_double(NULL, 4);

    double**  local_v   = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip      = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double** local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*  local_p   = BB_std_calloc_1d_double(NULL, nl);
    double*  p_ip      = BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip = BB_std_calloc_2d_double(NULL, np, 3);

    double* mu_eff_ip = BB_std_calloc_1d_double(NULL, np);
    double* tau_ip    = BB_std_calloc_1d_double(NULL, np);
    double* tau_c_ip  = BB_std_calloc_1d_double(NULL, np);

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

            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p],
                nl,
                local_v,
                fe->geo[e][p].grad_N);

            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p],
                nl,
                local_p,
                fe->geo[e][p].grad_N);

            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }

        const double vol = BBFE_std_integ_calc_volume(
            np,
            basis->integ_weight,
            Jacobian_ip);

        const double h_e = cbrt(vol);

        for (int p = 0; p < np; ++p) {
            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            BBFE_vms_mu_eff_tau(
                &mu_eff_ip[p],
                &tau_ip[p],
                &tau_c_ip[p],
                J_inv,
                fe->geo[e][p].Jacobian,
                h_e,
                v_ip[p],
                v_ip_old[p],
                grad_v_ip[p],
                grad_p_ip[p],
                vals->density,
                vals->viscosity,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);
        }

        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p) {
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = 0.0;
                }
            }

            for (int p = 0; p < np; ++p) {
                double du_time[3] = {
                    v_ip[p][0] - v_ip_old[p][0],
                    v_ip[p][1] - v_ip_old[p][1],
                    v_ip[p][2] - v_ip_old[p][2]
                };

                BB_calc_mat3d_inverse(
                    fe->geo[e][p].J,
                    fe->geo[e][p].Jacobian,
                    J_inv);

                const double mu_eff = mu_eff_ip[p];
                const double tau    = tau_ip[p];
                const double tau_c  = tau_c_ip[p];

                BBFE_elemmat_fluid_vec_rom_nonlinear(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density,
                    mu_eff,
                    tau,
                    tau_c,
                    vals->dt,
                    du_time);

                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                    vec[d] = 0.0;
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np,
                    val_ip_vec[d],
                    basis->integ_weight,
                    Jacobian_ip);

                monolis->mat.R.B[4 * fe->conn[e][i] + d] -= integ_val_vec[d];
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

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);

    BB_std_free_1d_double(mu_eff_ip, np);
    BB_std_free_1d_double(tau_ip, np);
    BB_std_free_1d_double(tau_c_ip, np);
}