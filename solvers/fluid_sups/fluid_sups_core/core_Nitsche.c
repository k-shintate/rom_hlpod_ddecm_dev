

static int inv3x3(
    const double A[3][3],
    double invA[3][3]);

static double det3x3_local(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int hex8_grad_phys_metric(
    const double xvol[8][3],
    const double dN_dxi[8][3],
    double dN_dx[8][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double J[3][3] = {{0.0}};

    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = det3x3_local(J);
    if (fabs(detJ) < 1.0e-300) return 0;

    if (!inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);

    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void hex8_shape_grad_ref(
    const int SGN[8][3],
    const double xi[3],
    double N[8],
    double dN_dxi[8][3])
{
    for (int a = 0; a < 8; ++a) {
        const double sx = (double)SGN[a][0];
        const double sy = (double)SGN[a][1];
        const double sz = (double)SGN[a][2];

        const double X = 1.0 + sx*xi[0];
        const double Y = 1.0 + sy*xi[1];
        const double Z = 1.0 + sz*xi[2];

        N[a] = 0.125 * X * Y * Z;
        dN_dxi[a][0] = 0.125 * sx * Y * Z;
        dN_dxi[a][1] = 0.125 * X * sy * Z;
        dN_dxi[a][2] = 0.125 * X * Y * sz;
    }
}

static int inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0/det;
    /* invA = adj(A)/det = C^T/det */
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static int hex8_grad_phys(
    const double xvol[8][3],
    const double dN_dxi[8][3],
    double dN_dx[8][3])
{
    double J[3][3] = {{0.0}};

    /* J[i][j] = dx_i / dxi_j */
    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    double invJ[3][3];
    if (!inv3x3(J, invJ)) return 0;

    /* dN/dx_i = sum_j dN/dxi_j * dxi_j/dx_i = dN_dxi[j] * invJ[j][i] */
    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * invJ[j][i];
            }
        }
    }
    return 1;
}

static double eps_nn_from_grad_u(
    const double grad_u[3][3], /* grad_u[a][b] = d u_a / d x_b */
    const double n[3])
{
    double val = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double eps_ab = 0.5*(grad_u[a][b] + grad_u[b][a]);
            val += n[a] * eps_ab * n[b];
        }
    }
    return val;
}


static double spalding_residual_uplus(
    double u_plus,
    double Re_y,
    double kappa,
    double A)
{
    const double x = kappa * u_plus;

    double y_plus;

    /*
      Avoid overflow of exp(x).
      x = 50 is already very large for this purpose.
    */
    double DBL_MAX = 10e13;
    if (x > 50.0) {
        return DBL_MAX * 0.25;
    }

    y_plus =
        u_plus
        + A * (
            exp(x)
            - 1.0
            - x
            - 0.5 * x * x
            - (x * x * x) / 6.0
        );

    return u_plus * y_plus - Re_y;
}


static double compute_utau_all_yplus(
    double ut_norm,
    double y_wall,
    double rho,
    double mu_eff)
{
    const double eps = 1.0e-30;

    if (ut_norm <= eps) return 0.0;
    if (y_wall  <= eps) return 0.0;
    if (rho     <= eps) return 0.0;
    if (mu_eff  <= eps) return 0.0;

    const double nu = mu_eff / rho;
    if (nu <= eps) return 0.0;

    /*
      Re_y = U y / nu
           = u+ y+
    */
    const double Re_y = ut_norm * y_wall / nu;
    if (Re_y <= eps) return 0.0;

    const double kappa = 0.41;
    const double B     = 5.2;
    const double A     = exp(-kappa * B);

    /*
      Solve:
        F(u+) = u+ y+_Spalding(u+) - Re_y = 0
    */
    double lo = 1.0e-12;
    double hi = 1.0;

    while (
        spalding_residual_uplus(hi, Re_y, kappa, A) < 0.0
        && hi < 1.0e6
    ) {
        hi *= 2.0;
    }

    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double fmid =
            spalding_residual_uplus(mid, Re_y, kappa, A);

        if (fmid > 0.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double u_plus = 0.5 * (lo + hi);

    if (u_plus <= eps) return 0.0;

    return ut_norm / u_plus;
}



static void wall_sym_nitsche_local_vecmat(
    double       vec_i[4],
    double       mat_ij[4][4],
    double       Ni,
    double       Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    double       p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double       rho,
    double       mu_eff,       /* VMS 有効粘性: 法線 Nitsche 用 */
    double       mu_wall_law,  /* 分子粘性: all-y+ wall law 用 */
    double       h_n,
    double       y_wall,
    double       dt,
    double       gamma_n)
{
    for (int a = 0; a < 4; ++a) {
        vec_i[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            mat_ij[a][b] = 0.0;
        }
    }

    const double eps = 1.0e-30;

    const double dt_eff = fmax(dt, eps);
    const double h_eff  = fmax(h_n, 1.0e-12);
    const double y_eff  = fmax(y_wall, 1.0e-12);

    double ur[3];
    for (int a = 0; a < 3; ++a) {
        ur[a] = u_q[a] - Uw[a];
    }

    const double un = dot3(ur, n);

    double ut[3];
    for (int a = 0; a < 3; ++a) {
        ut[a] = ur[a] - un * n[a];
    }

    const double ut_norm = sqrt(dot3(ut, ut));

    /*
      法線方向: no-penetration を symmetric Nitsche で課す。
      接線方向: penalty ではなく all-y+ wall law にする。
    */
    const double beta_base =
        mu_eff / h_eff + rho * h_eff / dt_eff;

    const double beta_n = gamma_n * beta_base;

    /*
      all-y+ wall law:
        tau_w = rho * u_tau^2
        tau_vec = beta_wall * u_t
        beta_wall = rho * u_tau^2 / |u_t|

      |u_t| -> 0 では Spalding の粘性底層極限として
        beta_wall ≈ mu_wall_law / y_wall
      を使う。
    */
    double beta_wall = mu_wall_law / y_eff;

    if (ut_norm > 1.0e-14) {
        const double u_tau =
            compute_utau_all_yplus(
                ut_norm,
                y_eff,
                rho,
                mu_wall_law);

        if (u_tau > 0.0) {
            beta_wall = rho * u_tau * u_tau / ut_norm;
        }
    }

    const double gni = dot3(gradNi, n);
    const double gnj = dot3(gradNj, n);

    const double epsnn_u = eps_nn_from_grad_u(grad_u_q, n);
    const double tnn_u   = -p_q + 2.0 * mu_eff * epsnn_u;

    /*
      Residual:
        法線: symmetric Nitsche
        接線: all-y+ wall law Robin
    */
    for (int a = 0; a < 3; ++a) {
        /* normal consistency */
        vec_i[a] += dt * (-Ni * n[a] * tnn_u);

        /* normal adjoint consistency */
        vec_i[a] += dt * (-(2.0 * mu_eff * n[a] * gni) * un);

        /* normal penalty */
        vec_i[a] += dt * (Ni * beta_n * un * n[a]);

        /* tangential all-y+ wall law */
        vec_i[a] += dt * (Ni * beta_wall * ut[a]);
    }

    /*
      pressure test q = Ni:
      no-penetration Nitsche の pressure consistency
    */
    vec_i[3] += dt * Ni * un;

    /*
      Tangent matrix.
      まずは beta_wall を frozen とする Picard/Newton 混合。
      max_iter_NR = 1 の現状ではこの方が安定。
    */
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double Pt_ab =
                (a == b ? 1.0 : 0.0) - n[a] * n[b];

            /* d/d u_b of normal consistency term */
            mat_ij[a][b] +=
                dt * (-Ni * n[a] *
                      (2.0 * mu_eff * n[b] * gnj));

            /* d/d u_b of normal adjoint-consistency term */
            mat_ij[a][b] +=
                dt * (-(2.0 * mu_eff * n[a] * gni) *
                      Nj * n[b]);

            /* normal penalty */
            mat_ij[a][b] +=
                dt * (Ni * Nj * beta_n * n[a] * n[b]);

            /* tangential all-y+ wall law, frozen beta_wall */
            mat_ij[a][b] +=
                dt * (Ni * Nj * beta_wall * Pt_ab);
        }

        /* d/dp of normal consistency term */
        mat_ij[a][3] +=
            dt * (Ni * Nj * n[a]);
    }

    for (int b = 0; b < 3; ++b) {
        mat_ij[3][b] +=
            dt * (Ni * Nj * n[b]);
    }
}

void set_wall_face_vecmat_symmetric_nitsche(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    /*
     * y+ visualization arrays reset.
     * この関数が呼ばれるたびに、最新の壁面 y+ を作り直す。
     */
    if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            vals->y_plus[i] = 0.0;
            vals->y_plus_count[i] = 0;
        }
    }

    int faces_map[6][4];
    build_hex8_faces_from_sgn(SGN, faces_map);

    int nVolF = 0, nSurfF = 0;
    dl_FaceEntry* volF = dl_build_vol_faces_sorted(fe, faces_map, &nVolF);
    dl_SurfFaceEntry* sF = dl_build_surf_faces_sorted(surf, &nSurfF);

    if (volF == NULL || sF == NULL) {
        fprintf(stderr, "[wall] failed to build face lists\n");
        free(volF);
        free(sF);
        return;
    }

    int* owner = (int*)malloc(sizeof(int) * nSurfF);
    int* lface = (int*)malloc(sizeof(int) * nSurfF);

    if (owner == NULL || lface == NULL) {
        fprintf(stderr, "malloc failed in set_wall_face_vecmat_symmetric_nitsche\n");
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        exit(EXIT_FAILURE);
    }

    for (int e = 0; e < nSurfF; ++e) {
        owner[e] = -1;
        lface[e] = -1;
    }

    int matched = dl_build_owner_map_mergejoin(
        volF,
        nVolF,
        sF,
        nSurfF,
        owner,
        lface);

    if (matched != nSurfF) {
        fprintf(stderr, "[wall] surface-owner mapping failed\n");
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return;
    }

    const int ns = surf->local_num_nodes;      /* normally 4 */
    const int np = basis_surf->num_integ_points;

    if (ns != 4) {
        fprintf(stderr,
            "[wall] set_wall_face_vecmat_symmetric_nitsche assumes QUAD4 surface, ns=%d\n",
            ns);
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return;
    }

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];

        if (ke < 0) {
            fprintf(stderr,
                "[wall][BUG] owner not found: es=%d owner=%d\n",
                es, ke);
            exit(EXIT_FAILURE);
        }

        {
            int common = 0;

            for (int a = 0; a < ns; ++a) {
                const int sgid = surf->conn[es][a];

                for (int b = 0; b < 8; ++b) {
                    if (fe->conn[ke][b] == sgid) {
                        ++common;
                        break;
                    }
                }
            }

            if (common != ns) {
                fprintf(stderr,
                    "[wall][BUG] surface-owner mismatch: es=%d ke=%d common=%d/%d\n",
                    es, ke, common, ns);

                fprintf(stderr, "  surf nodes:");
                for (int a = 0; a < ns; ++a) {
                    fprintf(stderr, " %d", surf->conn[es][a]);
                }

                fprintf(stderr, "\n  vol nodes:");
                for (int b = 0; b < 8; ++b) {
                    fprintf(stderr, " %d", fe->conn[ke][b]);
                }
                fprintf(stderr, "\n");

                exit(EXIT_FAILURE);
            }
        }

        double xsurf[4][3];
        for (int a = 0; a < ns; ++a) {
            int gid = surf->conn[es][a];
            xsurf[a][0] = fe->x[gid][0];
            xsurf[a][1] = fe->x[gid][1];
            xsurf[a][2] = fe->x[gid][2];
        }

        double xvol[8][3];
        for (int a = 0; a < 8; ++a) {
            int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        double vol_e = 0.0;
        for (int pp = 0; pp < basis_vol->num_integ_points; ++pp) {
            vol_e += basis_vol->integ_weight[pp] * fe->geo[ke][pp].Jacobian;
        }

        double h_e_vms = cbrt(fabs(vol_e));
        if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;

        int s2v[4];
        int ok_s2v = 1;

        for (int a = 0; a < ns; ++a) {
            s2v[a] = -1;

            int sgid = surf->conn[es][a];

            for (int b = 0; b < 8; ++b) {
                if (fe->conn[ke][b] == sgid) {
                    s2v[a] = b;
                    break;
                }
            }

            if (s2v[a] < 0) {
                fprintf(stderr, "[wall] cannot map surface node to volume node\n");
                ok_s2v = 0;
                break;
            }
        }

        if (!ok_s2v) continue;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 8; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 8.0;
        xc[1] /= 8.0;
        xc[2] /= 8.0;

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < ns; ++a) {
            xf[0] += xsurf[a][0];
            xf[1] += xsurf[a][1];
            xf[2] += xsurf[a][2];
        }
        xf[0] /= (double)ns;
        xf[1] /= (double)ns;
        xf[2] /= (double)ns;

        double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < np; ++p) {
            double nrm[3];

            dl_quad4_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nrm);

            const double Jface = sqrt(dot3(nrm, nrm));
            if (Jface <= 1.0e-300) continue;

            double n[3] = {
                nrm[0] / Jface,
                nrm[1] / Jface,
                nrm[2] / Jface
            };

            if (dot3(n, svec) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double xi3[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < ns; ++a) {
                int lv = s2v[a];
                double Na = basis_surf->N[p][a];

                xi3[0] += Na * (double)SGN[lv][0];
                xi3[1] += Na * (double)SGN[lv][1];
                xi3[2] += Na * (double)SGN[lv][2];
            }

            double Nv[8];
            double dN_dxi[8][3];
            double dN_dx[8][3];

            double J_inv_face[3][3];
            double Jacobian_face = 0.0;

            hex8_shape_grad_ref(SGN, xi3, Nv, dN_dxi);

            if (!hex8_grad_phys_metric(
                    xvol,
                    dN_dxi,
                    dN_dx,
                    J_inv_face,
                    &Jacobian_face)) {
                continue;
            }

            double u_q[3]     = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double p_q        = 0.0;

            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3]    = {0.0, 0.0, 0.0};

            for (int a = 0; a < 8; ++a) {
                int gid = fe->conn[ke][a];

                for (int c = 0; c < 3; ++c) {
                    u_q[c]     += Nv[a] * vals->v[gid][c];
                    u_old_q[c] += Nv[a] * vals->v_old[gid][c];

                    for (int d = 0; d < 3; ++d) {
                        grad_u_q[c][d] += dN_dx[a][d] * vals->v[gid][c];
                    }
                }

                p_q += Nv[a] * vals->p[gid];

                for (int d = 0; d < 3; ++d) {
                    grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                }
            }

            double* grad_u_ptr[3] = {
                grad_u_q[0],
                grad_u_q[1],
                grad_u_q[2]
            };

            double mu_eff_face = mu_molecular;
            double tau_face    = 0.0;
            double tau_c_face  = 0.0;

            BBFE_vms_mu_eff_tau(
                &mu_eff_face,
                &tau_face,
                &tau_c_face,
                J_inv_face,
                Jacobian_face,
                h_e_vms,
                u_q,
                u_old_q,
                grad_u_ptr,
                grad_p_q,
                rho,
                mu_molecular,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);

            double h_n =
                fabs((xf[0] - xc[0]) * n[0]
                   + (xf[1] - xc[1]) * n[1]
                   + (xf[2] - xc[2]) * n[2]);

            if (h_n <= 1.0e-12) h_n = 1.0e-12;

            const double y_wall = h_n;
            const double w = basis_surf->integ_weight[p] * Jface;

            /*
             * ------------------------------------------------------------
             * y+ visualization value at this wall quadrature point
             * ------------------------------------------------------------
             */
            if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
                double ur[3];
                for (int a = 0; a < 3; ++a) {
                    ur[a] = u_q[a] - Uw[a];
                }

                const double un = dot3(ur, n);

                double ut[3];
                for (int a = 0; a < 3; ++a) {
                    ut[a] = ur[a] - un * n[a];
                }

                const double ut_norm = sqrt(dot3(ut, ut));

                const double u_tau =
                    compute_utau_all_yplus(
                        ut_norm,
                        y_wall,
                        rho,
                        mu_molecular);

                double yplus_q = 0.0;

                if (u_tau > 0.0 && rho > 0.0 && mu_molecular > 0.0) {
                    yplus_q = rho * u_tau * y_wall / mu_molecular;
                }

                /*
                 * 面上の4節点へ蓄積する。
                 * 求積点重み付き平均にしたい場合は、
                 * y_plus_weight 配列を別途持たせて w で平均する。
                 */
                for (int a = 0; a < ns; ++a) {
                    int gid = surf->conn[es][a];

                    vals->y_plus[gid] += yplus_q;
                    vals->y_plus_count[gid] += 1;
                }
            }

            /*
             * Original Nitsche residual and matrix assembly
             */
            for (int i = 0; i < 8; ++i) {
                int gi = fe->conn[ke][i];

                double vec_i[4];
                double dummy_mat[4][4];
                double zero_grad[3] = {0.0, 0.0, 0.0};

                wall_sym_nitsche_local_vecmat(
                    vec_i,
                    dummy_mat,
                    Nv[i],
                    0.0,
                    dN_dx[i],
                    zero_grad,
                    n,
                    u_q,
                    p_q,
                    grad_u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu_molecular,
                    h_n,
                    y_wall,
                    vals->dt,
                    gamma_n);

                for (int a = 0; a < 4; ++a) {
                    monolis->mat.R.B[4 * gi + a] -= w * vec_i[a];
                }

                for (int j = 0; j < 8; ++j) {
                    int gj = fe->conn[ke][j];

                    double dummy_vec[4];
                    double mat_ij[4][4];

                    wall_sym_nitsche_local_vecmat(
                        dummy_vec,
                        mat_ij,
                        Nv[i],
                        Nv[j],
                        dN_dx[i],
                        dN_dx[j],
                        n,
                        u_q,
                        p_q,
                        grad_u_q,
                        Uw,
                        rho,
                        mu_eff_face,
                        mu_molecular,
                        h_n,
                        y_wall,
                        vals->dt,
                        gamma_n);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                gi,
                                gj,
                                a,
                                b,
                                w * mat_ij[a][b]);
                        }
                    }
                }
            }
        }
    }

    /*
     * Final average at wall nodes.
     * Non-wall nodes are set to 0.0.
     */
    if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            if (vals->y_plus_count[i] > 0) {
                vals->y_plus[i] /= (double)vals->y_plus_count[i];
            } else {
                vals->y_plus[i] = 0.0;
            }
        }
    }

    free(volF);
    free(sF);
    free(owner);
    free(lface);
}

/* ---- TET4/TRI3 face-key utilities ---- */
typedef long long t4_node_id_t;

static inline void t4_sort3(t4_node_id_t k[3])
{
    if (k[1] < k[0]) { t4_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
    if (k[2] < k[1]) { t4_node_id_t t = k[1]; k[1] = k[2]; k[2] = t; }
    if (k[1] < k[0]) { t4_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
}


static inline double t4_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double t4_norm3(const double a[3])
{
    return sqrt(t4_dot3(a, a));
}

static inline void t4_normalize3(double a[3])
{
    const double n = t4_norm3(a);
    if (n > 0.0) {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}

static inline void t4_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline int t4_cmp3(const t4_node_id_t a[3], const t4_node_id_t b[3])
{
    if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : +1;
    if (a[1] != b[1]) return (a[1] < b[1]) ? -1 : +1;
    if (a[2] != b[2]) return (a[2] < b[2]) ? -1 : +1;
    return 0;
}

/* Standard outward-oriented faces for a positive-orientation TET4.
 * The orientation is not required for the owner map because keys are sorted;
 * normals are still corrected by centroid direction in the integral loop.
 */
static const int T4_FACE[4][3] = {
    {1, 2, 3},  /* opposite node 0 */
    {0, 3, 2},  /* opposite node 1 */
    {0, 1, 3},  /* opposite node 2 */
    {0, 2, 1}   /* opposite node 3 */
};

typedef struct {
    t4_node_id_t key[3];
    int owner_elem;
    int owner_lface;
} t4_FaceEntry;

typedef struct {
    t4_node_id_t key[3];
    int surf_face_id;
} t4_SurfFaceEntry;

static int t4_cmp_face(const void* A, const void* B)
{
    const t4_FaceEntry* a = (const t4_FaceEntry*)A;
    const t4_FaceEntry* b = (const t4_FaceEntry*)B;
    int c = t4_cmp3(a->key, b->key);
    if (c) return c;
    if (a->owner_elem  != b->owner_elem)  return (a->owner_elem  < b->owner_elem)  ? -1 : +1;
    if (a->owner_lface != b->owner_lface) return (a->owner_lface < b->owner_lface) ? -1 : +1;
    return 0;
}

static int t4_cmp_sface(const void* A, const void* B)
{
    const t4_SurfFaceEntry* a = (const t4_SurfFaceEntry*)A;
    const t4_SurfFaceEntry* b = (const t4_SurfFaceEntry*)B;
    int c = t4_cmp3(a->key, b->key);
    if (c) return c;
    if (a->surf_face_id != b->surf_face_id) return (a->surf_face_id < b->surf_face_id) ? -1 : +1;
    return 0;
}

static t4_FaceEntry* t4_build_vol_faces_sorted(const BBFE_DATA* fe, int* nfaces_out)
{
    const int nfaces = fe->total_num_elems * 4;
    t4_FaceEntry* arr = (t4_FaceEntry*)malloc(sizeof(t4_FaceEntry) * nfaces);
    if (!arr) {
        fprintf(stderr, "[t4][ERR] malloc vol faces failed\n");
        *nfaces_out = 0;
        return NULL;
    }

    int k = 0;
    for (int e = 0; e < fe->total_num_elems; ++e) {
        for (int f = 0; f < 4; ++f) {
            t4_node_id_t key[3] = {
                (t4_node_id_t)fe->conn[e][T4_FACE[f][0]],
                (t4_node_id_t)fe->conn[e][T4_FACE[f][1]],
                (t4_node_id_t)fe->conn[e][T4_FACE[f][2]]
            };
            t4_sort3(key);
            for (int j = 0; j < 3; ++j) arr[k].key[j] = key[j];
            arr[k].owner_elem  = e;
            arr[k].owner_lface = f;
            ++k;
        }
    }

    qsort(arr, nfaces, sizeof(t4_FaceEntry), t4_cmp_face);
    *nfaces_out = nfaces;
    return arr;
}

static t4_SurfFaceEntry* t4_build_surf_faces_sorted(const BBFE_DATA* surf, int* nsurf_out)
{
    const int nsurf = surf->total_num_elems;
    t4_SurfFaceEntry* arr = (t4_SurfFaceEntry*)malloc(sizeof(t4_SurfFaceEntry) * nsurf);
    if (!arr) {
        fprintf(stderr, "[t4][ERR] malloc surf faces failed\n");
        *nsurf_out = 0;
        return NULL;
    }

    for (int e = 0; e < nsurf; ++e) {
        t4_node_id_t key[3] = {
            (t4_node_id_t)surf->conn[e][0],
            (t4_node_id_t)surf->conn[e][1],
            (t4_node_id_t)surf->conn[e][2]
        };
        t4_sort3(key);
        for (int j = 0; j < 3; ++j) arr[e].key[j] = key[j];
        arr[e].surf_face_id = e;
    }

    qsort(arr, nsurf, sizeof(t4_SurfFaceEntry), t4_cmp_sface);
    *nsurf_out = nsurf;
    return arr;
}

static int t4_build_owner_map_mergejoin(
    const t4_FaceEntry* volF,
    int nVolF,
    const t4_SurfFaceEntry* surfF,
    int nSurfF,
    int* owner_elem,
    int* owner_lface)
{
    int i = 0;
    int j = 0;
    int matches = 0;

    while (i < nVolF && j < nSurfF) {
        int c = t4_cmp3(volF[i].key, surfF[j].key);

        if (c == 0) {
            const int es = surfF[j].surf_face_id;
            owner_elem[es]  = volF[i].owner_elem;
            owner_lface[es] = volF[i].owner_lface;
            ++matches;

            t4_node_id_t key0[3] = { volF[i].key[0], volF[i].key[1], volF[i].key[2] };
            do { ++i; } while (i < nVolF && t4_cmp3(key0, volF[i].key) == 0);
            ++j;
        } else if (c < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    return matches;
}

static int t4_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner_out,
    int** lface_out)
{
    *owner_out = NULL;
    *lface_out = NULL;

    if (fe->local_num_nodes != 4 || surf->local_num_nodes != 3) {
        fprintf(stderr,
            "[t4] expected TET4/TRI3 but got fe->local_num_nodes=%d surf->local_num_nodes=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return -1;
    }

    int nVolF = 0;
    int nSurfF = 0;
    t4_FaceEntry* volF = t4_build_vol_faces_sorted(fe, &nVolF);
    t4_SurfFaceEntry* sF = t4_build_surf_faces_sorted(surf, &nSurfF);
    if (!volF || !sF) {
        free(volF);
        free(sF);
        return -1;
    }

    int* owner = (int*)malloc(sizeof(int) * nSurfF);
    int* lface = (int*)malloc(sizeof(int) * nSurfF);
    if (!owner || !lface) {
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return -1;
    }

    for (int e = 0; e < nSurfF; ++e) {
        owner[e] = -1;
        lface[e] = -1;
    }

    const int matched = t4_build_owner_map_mergejoin(volF, nVolF, sF, nSurfF, owner, lface);
    if (matched != nSurfF) {
        int cnt = 0;
        for (int e = 0; e < nSurfF; ++e) {
            if (owner[e] < 0) {
                fprintf(stderr, "[t4] surf tri %d has NO owner\n", e);
                ++cnt;
            }
        }
        fprintf(stderr,
            "[t4] %d / %d surface triangles unmapped. Check node IDs / connectivity.\n",
            cnt,
            nSurfF);
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return -1;
    }

    free(volF);
    free(sF);
    *owner_out = owner;
    *lface_out = lface;
    return 0;
}

/* normal = dx/dr x dx/ds.  Its norm is 2*area for an affine TRI3. */
static inline void t4_tri3_normal(
    const double xtri[3][3],
    const double* dN_dr,
    const double* dN_ds,
    double normal[3])
{
    double xr[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};

    for (int a = 0; a < 3; ++a) {
        xr[0] += dN_dr[a] * xtri[a][0];
        xr[1] += dN_dr[a] * xtri[a][1];
        xr[2] += dN_dr[a] * xtri[a][2];
        xs[0] += dN_ds[a] * xtri[a][0];
        xs[1] += dN_ds[a] * xtri[a][1];
        xs[2] += dN_ds[a] * xtri[a][2];
    }

    t4_cross3(xr, xs, normal);
}

/* Correct quadrature-weight convention for reference triangle.
 * If sum(weights)=0.5, returns 1.0. If sum(weights)=1.0, returns 0.5.
 */
static double t4_tri_area_factor(const BBFE_BASIS* basis_surf)
{
    double sw = 0.0;
    for (int p = 0; p < basis_surf->num_integ_points; ++p) {
        sw += basis_surf->integ_weight[p];
    }

    if (fabs(sw) <= 1.0e-300) return 1.0;
    return 0.5 / sw;
}

static int t4_inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0 / det;
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static double t4_det3x3(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

/* TET4 reference gradients for N0=1-r-s-t, N1=r, N2=s, N3=t. */
static void t4_shape_grad_ref(double dN_dxi[4][3])
{
    dN_dxi[0][0] = -1.0; dN_dxi[0][1] = -1.0; dN_dxi[0][2] = -1.0;
    dN_dxi[1][0] =  1.0; dN_dxi[1][1] =  0.0; dN_dxi[1][2] =  0.0;
    dN_dxi[2][0] =  0.0; dN_dxi[2][1] =  1.0; dN_dxi[2][2] =  0.0;
    dN_dxi[3][0] =  0.0; dN_dxi[3][1] =  0.0; dN_dxi[3][2] =  1.0;
}

static int t4_grad_phys_metric(
    const double xvol[4][3],
    double dN_dx[4][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double dN_dxi[4][3];
    t4_shape_grad_ref(dN_dxi);

    double J[3][3] = {{0.0}};
    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = t4_det3x3(J);
    if (fabs(detJ) < 1.0e-300) return 0;
    if (!t4_inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);  /* 6 * physical tet volume */

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

static void t4_grad_tet4_exact(int ke, const BBFE_DATA* fe, const VALUES* vals, double G[3][3])
{
    memset(G, 0, sizeof(double) * 9);

    if (fe->local_num_nodes != 4) return;

    double xvol[4][3];
    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];
        xvol[a][0] = fe->x[gid][0];
        xvol[a][1] = fe->x[gid][1];
        xvol[a][2] = fe->x[gid][2];
    }

    double dN_dx[4][3];
    double J_inv[3][3];
    double Jacobian = 0.0;
    if (!t4_grad_phys_metric(xvol, dN_dx, J_inv, &Jacobian)) return;

    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                G[i][j] += vals->v[gid][i] * dN_dx[a][j];
            }
        }
    }
}

static void t4_compute_all_gradients_tet4_exact(
    const BBFE_DATA* fe,
    const VALUES* vals,
    double (*G_all)[3][3])
{
    const int ne = fe->total_num_elems;

    #pragma omp parallel for if(ne > 256) schedule(static) default(none) shared(fe, vals, G_all, ne)
    for (int ke = 0; ke < ne; ++ke) {
        t4_grad_tet4_exact(ke, fe, vals, G_all[ke]);
    }
}

static int t4_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for (int a = 0; a < 3; ++a) {
        s2v[a] = -1;
        const int sgid = surf->conn[es][a];
        for (int b = 0; b < 4; ++b) {
            if (fe->conn[ke][b] == sgid) {
                s2v[a] = b;
                break;
            }
        }
        if (s2v[a] < 0) return 0;
    }
    return 1;
}

/* Pressure contribution on TRI3/TET4 surface: integrate -p n dS. */
void calc_Cd_p_tet4_tri3(
    BBFE_DATA* surf,
    BBFE_DATA* fe,
    BBFE_BASIS* basis_surf,
    VALUES* vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out)
{
    if (D_out) *D_out = 0.0;
    if (L_out) *L_out = 0.0;

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        return;
    }

    double eU[3] = { eU_in[0], eU_in[1], eU_in[2] };
    double eP[3] = { eP_in[0], eP_in[1], eP_in[2] };
    t4_normalize3(eU);

    /* Gram-Schmidt: make eP perpendicular to eU. */
    const double ep_dot_eu = t4_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];
    t4_normalize3(eP);

    const double area_factor = t4_tri_area_factor(basis_surf);
    double Fp[3] = {0.0, 0.0, 0.0};

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0; xf[1] /= 3.0; xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xc[0] += fe->x[gid][0];
            xc[1] += fe->x[gid][1];
            xc[2] += fe->x[gid][2];
        }
        xc[0] /= 4.0; xc[1] /= 4.0; xc[2] /= 4.0;

        const double svec[3] = { xf[0] - xc[0], xf[1] - xc[1], xf[2] - xc[2] };

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(xtri, basis_surf->dN_dxi[p], basis_surf->dN_det[p], nraw);

            const double Jraw = t4_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = { nraw[0]/Jraw, nraw[1]/Jraw, nraw[2]/Jraw };
            if (t4_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double p_q = 0.0;
            for (int a = 0; a < 3; ++a) {
                const int gid = surf->conn[es][a];
                p_q += basis_surf->N[p][a] * vals->p[gid];
            }
//              p_q = 1.0;
            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;
            Fp[0] += (p_q) * n[0] * w;
            Fp[1] += (p_q) * n[1] * w;
            Fp[2] += (p_q) * n[2] * w;
        }
    }

    if (D_out) *D_out = t4_dot3(Fp, eU);
    if (L_out) *L_out = t4_dot3(Fp, eP);

    free(owner);
    free(lface);
}

/* Viscous contribution on TRI3/TET4 surface.
 * Sign convention follows your current calc_Cd_v:
 *   integrate tau*n on the fluid side, then report body force = -integral.
 * The current HEX8 code stores D/L, not D/qA or L/qA, in Cd_out/Cl_out;
 * this function keeps that behavior. Change the final two assignments if you
 * want nondimensional coefficients.
 */


/* -------------------------------------------------------------------------- */
/* Small TET4/TRI3 utilities                                                  */
/* -------------------------------------------------------------------------- */

typedef long long t4n_node_id_t;

static inline double t4n_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double t4n_norm3(const double a[3])
{
    return sqrt(t4n_dot3(a, a));
}

static inline void t4n_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void t4n_sort3(t4n_node_id_t k[3])
{
    if (k[1] < k[0]) { t4n_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
    if (k[2] < k[1]) { t4n_node_id_t t = k[1]; k[1] = k[2]; k[2] = t; }
    if (k[1] < k[0]) { t4n_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
}

static inline int t4n_cmp3(const t4n_node_id_t a[3], const t4n_node_id_t b[3])
{
    if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : +1;
    if (a[1] != b[1]) return (a[1] < b[1]) ? -1 : +1;
    if (a[2] != b[2]) return (a[2] < b[2]) ? -1 : +1;
    return 0;
}

/* Faces of TET4. Orientation is not used for matching because the keys are sorted.
 * Outward normal is corrected geometrically by the vector from cell centroid to face centroid.
 */
static const int T4N_FACE[4][3] = {
    {1, 2, 3},
    {0, 3, 2},
    {0, 1, 3},
    {0, 2, 1}
};

typedef struct {
    t4n_node_id_t key[3];
    int owner_elem;
    int owner_lface;
} t4n_FaceEntry;

typedef struct {
    t4n_node_id_t key[3];
    int surf_face_id;
} t4n_SurfFaceEntry;

static int t4n_cmp_face(const void* A, const void* B)
{
    const t4n_FaceEntry* a = (const t4n_FaceEntry*)A;
    const t4n_FaceEntry* b = (const t4n_FaceEntry*)B;
    int c = t4n_cmp3(a->key, b->key);
    if (c) return c;
    if (a->owner_elem != b->owner_elem) return (a->owner_elem < b->owner_elem) ? -1 : +1;
    if (a->owner_lface != b->owner_lface) return (a->owner_lface < b->owner_lface) ? -1 : +1;
    return 0;
}

static int t4n_cmp_sface(const void* A, const void* B)
{
    const t4n_SurfFaceEntry* a = (const t4n_SurfFaceEntry*)A;
    const t4n_SurfFaceEntry* b = (const t4n_SurfFaceEntry*)B;
    int c = t4n_cmp3(a->key, b->key);
    if (c) return c;
    if (a->surf_face_id != b->surf_face_id) return (a->surf_face_id < b->surf_face_id) ? -1 : +1;
    return 0;
}

static t4n_FaceEntry* t4n_build_vol_faces_sorted(const BBFE_DATA* fe, int* nfaces_out)
{
    const int nfaces = fe->total_num_elems * 4;
    t4n_FaceEntry* arr = (t4n_FaceEntry*)malloc(sizeof(t4n_FaceEntry) * nfaces);
    if (arr == NULL) {
        fprintf(stderr, "[t4n][ERR] malloc vol faces failed\n");
        *nfaces_out = 0;
        return NULL;
    }

    int k = 0;
    for (int e = 0; e < fe->total_num_elems; ++e) {
        for (int f = 0; f < 4; ++f) {
            t4n_node_id_t key[3] = {
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][0]],
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][1]],
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][2]]
            };
            t4n_sort3(key);
            for (int j = 0; j < 3; ++j) arr[k].key[j] = key[j];
            arr[k].owner_elem  = e;
            arr[k].owner_lface = f;
            ++k;
        }
    }

    qsort(arr, nfaces, sizeof(t4n_FaceEntry), t4n_cmp_face);
    *nfaces_out = nfaces;
    return arr;
}

static t4n_SurfFaceEntry* t4n_build_surf_faces_sorted(const BBFE_DATA* surf, int* nsurf_out)
{
    const int nsurf = surf->total_num_elems;
    t4n_SurfFaceEntry* arr = (t4n_SurfFaceEntry*)malloc(sizeof(t4n_SurfFaceEntry) * nsurf);
    if (arr == NULL) {
        fprintf(stderr, "[t4n][ERR] malloc surf faces failed\n");
        *nsurf_out = 0;
        return NULL;
    }

    for (int es = 0; es < nsurf; ++es) {
        t4n_node_id_t key[3] = {
            (t4n_node_id_t)surf->conn[es][0],
            (t4n_node_id_t)surf->conn[es][1],
            (t4n_node_id_t)surf->conn[es][2]
        };
        t4n_sort3(key);
        for (int j = 0; j < 3; ++j) arr[es].key[j] = key[j];
        arr[es].surf_face_id = es;
    }

    qsort(arr, nsurf, sizeof(t4n_SurfFaceEntry), t4n_cmp_sface);
    *nsurf_out = nsurf;
    return arr;
}

static int t4n_build_owner_map_mergejoin(
    const t4n_FaceEntry* volF,
    int nVolF,
    const t4n_SurfFaceEntry* surfF,
    int nSurfF,
    int* owner_elem,
    int* owner_lface)
{
    int i = 0;
    int j = 0;
    int matches = 0;

    while (i < nVolF && j < nSurfF) {
        int c = t4n_cmp3(volF[i].key, surfF[j].key);
        if (c == 0) {
            const int es = surfF[j].surf_face_id;
            owner_elem[es]  = volF[i].owner_elem;
            owner_lface[es] = volF[i].owner_lface;
            ++matches;

            t4n_node_id_t key0[3] = { volF[i].key[0], volF[i].key[1], volF[i].key[2] };
            do { ++i; } while (i < nVolF && t4n_cmp3(key0, volF[i].key) == 0);
            ++j;
        } else if (c < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    return matches;
}

static int t4n_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner_out,
    int** lface_out)
{
    *owner_out = NULL;
    *lface_out = NULL;

    if (fe->local_num_nodes != 4 || surf->local_num_nodes != 3) {
        fprintf(stderr,
            "[t4n] expected TET4/TRI3 but got fe->local_num_nodes=%d surf->local_num_nodes=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return -1;
    }

    int nVolF = 0;
    int nSurfF = 0;
    t4n_FaceEntry* volF = t4n_build_vol_faces_sorted(fe, &nVolF);
    t4n_SurfFaceEntry* sF = t4n_build_surf_faces_sorted(surf, &nSurfF);
    if (volF == NULL || sF == NULL) {
        free(volF);
        free(sF);
        return -1;
    }

    int* owner = (int*)malloc(sizeof(int) * nSurfF);
    int* lface = (int*)malloc(sizeof(int) * nSurfF);
    if (owner == NULL || lface == NULL) {
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return -1;
    }

    for (int es = 0; es < nSurfF; ++es) {
        owner[es] = -1;
        lface[es] = -1;
    }

    const int matched = t4n_build_owner_map_mergejoin(volF, nVolF, sF, nSurfF, owner, lface);
    if (matched != nSurfF) {
        int cnt = 0;
        for (int es = 0; es < nSurfF; ++es) {
            if (owner[es] < 0) {
                fprintf(stderr, "[t4n] surf tri %d has NO owner\n", es);
                ++cnt;
            }
        }
        fprintf(stderr,
            "[t4n] %d / %d surface triangles unmapped. Check node IDs / connectivity.\n",
            cnt,
            nSurfF);
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return -1;
    }

    free(volF);
    free(sF);
    *owner_out = owner;
    *lface_out = lface;
    return 0;
}

static int t4n_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for (int a = 0; a < 3; ++a) {
        s2v[a] = -1;
        const int sgid = surf->conn[es][a];
        for (int b = 0; b < 4; ++b) {
            if (fe->conn[ke][b] == sgid) {
                s2v[a] = b;
                break;
            }
        }
        if (s2v[a] < 0) return 0;
    }
    return 1;
}

/* normal = dx/dr x dx/ds. For affine TRI3, |normal| = 2 * physical area. */
static inline void t4n_tri3_normal(
    const double xtri[3][3],
    const double* dN_dr,
    const double* dN_ds,
    double normal[3])
{
    double xr[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};

    for (int a = 0; a < 3; ++a) {
        xr[0] += dN_dr[a] * xtri[a][0];
        xr[1] += dN_dr[a] * xtri[a][1];
        xr[2] += dN_dr[a] * xtri[a][2];
        xs[0] += dN_ds[a] * xtri[a][0];
        xs[1] += dN_ds[a] * xtri[a][1];
        xs[2] += dN_ds[a] * xtri[a][2];
    }

    t4n_cross3(xr, xs, normal);
}

/* If sum(weights)=0.5, returns 1.0. If sum(weights)=1.0, returns 0.5. */
static double t4n_tri_area_factor(const BBFE_BASIS* basis_surf)
{
    double sw = 0.0;
    for (int p = 0; p < basis_surf->num_integ_points; ++p) {
        sw += basis_surf->integ_weight[p];
    }

    if (fabs(sw) <= 1.0e-300) return 1.0;
    return 0.5 / sw;
}

static double t4n_det3x3(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int t4n_inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0 / det;
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static void t4n_shape_grad_ref(double dN_dxi[4][3])
{
    dN_dxi[0][0] = -1.0; dN_dxi[0][1] = -1.0; dN_dxi[0][2] = -1.0;
    dN_dxi[1][0] =  1.0; dN_dxi[1][1] =  0.0; dN_dxi[1][2] =  0.0;
    dN_dxi[2][0] =  0.0; dN_dxi[2][1] =  1.0; dN_dxi[2][2] =  0.0;
    dN_dxi[3][0] =  0.0; dN_dxi[3][1] =  0.0; dN_dxi[3][2] =  1.0;
}

static int t4n_grad_phys_metric(
    const double xvol[4][3],
    double dN_dx[4][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double dN_dxi[4][3];
    t4n_shape_grad_ref(dN_dxi);

    double J[3][3] = {{0.0}};
    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = t4n_det3x3(J);
    if (fabs(detJ) < 1.0e-300) return 0;
    if (!t4n_inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);  /* 6 * physical TET volume */

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

/* -------------------------------------------------------------------------- */
/* Local symmetric Nitsche + all-y+ wall law, TET helper copy                 */
/* -------------------------------------------------------------------------- */

static double t4n_eps_nn_from_grad_u(
    const double grad_u[3][3],
    const double n[3])
{
    double val = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double eps_ab = 0.5 * (grad_u[a][b] + grad_u[b][a]);
            val += n[a] * eps_ab * n[b];
        }
    }
    return val;
}

static double t4n_spalding_residual_uplus(
    double u_plus,
    double Re_y,
    double kappa,
    double A)
{
    const double x = kappa * u_plus;

    /* Keep the same overflow guard style as the verified HEX-side code. */
    const double huge_value = 10e13;
    if (x > 50.0) {
        return huge_value * 0.25;
    }

    const double y_plus =
        u_plus
        + A * (
            exp(x)
            - 1.0
            - x
            - 0.5 * x * x
            - (x * x * x) / 6.0
        );

    return u_plus * y_plus - Re_y;
}

static double t4n_compute_utau_all_yplus(
    double ut_norm,
    double y_wall,
    double rho,
    double mu_eff)
{
    const double eps = 1.0e-30;

    if (ut_norm <= eps) return 0.0;
    if (y_wall  <= eps) return 0.0;
    if (rho     <= eps) return 0.0;
    if (mu_eff  <= eps) return 0.0;

    const double nu = mu_eff / rho;
    if (nu <= eps) return 0.0;

    const double Re_y = ut_norm * y_wall / nu;
    if (Re_y <= eps) return 0.0;

    const double kappa = 0.41;
    const double B     = 5.2;
    const double A     = exp(-kappa * B);

    double lo = 1.0e-12;
    double hi = 1.0;

    while (
        t4n_spalding_residual_uplus(hi, Re_y, kappa, A) < 0.0
        && hi < 1.0e6) {
        hi *= 2.0;
    }

    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double fmid = t4n_spalding_residual_uplus(mid, Re_y, kappa, A);

        if (fmid > 0.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double u_plus = 0.5 * (lo + hi);
    if (u_plus <= eps) return 0.0;

    return ut_norm / u_plus;
}

#ifndef T4_NITSCHE_DT_PENALTY_COEFF
#define T4_NITSCHE_DT_PENALTY_COEFF 0.10
#endif

#ifndef T4_NITSCHE_INV_PENALTY_COEFF
#define T4_NITSCHE_INV_PENALTY_COEFF 1.0
#endif

#ifndef T4_NITSCHE_VIS_PENALTY_COEFF
#define T4_NITSCHE_VIS_PENALTY_COEFF 1.0
#endif

static void t4n_cauchy_traction(
    double       t[3],
    const double n[3],
    const double grad_u[3][3],
    double       p,
    double       mu_eff)
{
    for (int a = 0; a < 3; ++a) {
        double s = 0.0;
        for (int b = 0; b < 3; ++b) {
            s += (grad_u[a][b] + grad_u[b][a]) * n[b];
        }
        t[a] = -p * n[a] + mu_eff * s;
    }
}

static void t4n_velocity_test_traction(
    double       t[3],
    int          a_test,
    const double gradN[3],
    const double n[3],
    double       mu_eff)
{
    const double gn = t4n_dot3(gradN, n);

    for (int c = 0; c < 3; ++c) {
        t[c] = mu_eff * ((c == a_test ? gn : 0.0) + n[a_test] * gradN[c]);
    }
}


static void t4n_wall_sym_nitsche_local_vecmat(
    double       vec_i[4],
    double       mat_ij[4][4],
    double       Ni,
    double       Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    double       p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double       rho,
    double       mu_eff,
    double       mu_wall_law,
    double       h_n,
    double       y_wall,
    double       dt,
    double       gamma_n)
{
    (void)mu_wall_law;
    (void)y_wall;

    for (int a = 0; a < 4; ++a) {
        vec_i[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            mat_ij[a][b] = 0.0;
        }
    }

    const double eps = 1.0e-30;
    const double dt_eff = fmax(dt, eps);
    const double h_eff  = fmax(h_n, 1.0e-12);

    double r[3];
    for (int a = 0; a < 3; ++a) {
        r[a] = u_q[a] - Uw[a];
    }

    const double beta_inv =
        T4_NITSCHE_INV_PENALTY_COEFF * rho * h_eff / dt_eff;
    const double beta_vis =
        T4_NITSCHE_VIS_PENALTY_COEFF * mu_eff / h_eff;
    const double tau_B = gamma_n * fmax(beta_inv, beta_vis);

    double traction_u[3];
    t4n_cauchy_traction(
        traction_u,
        n,
        grad_u_q,
        p_q,
        mu_eff);

    /*
      Velocity test rows:

        - w . sigma(u,p)n
        - sigma(w,0)n . r
        + tau_B w . r
    */
    for (int a = 0; a < 3; ++a) {
        double traction_test_a[3];
        t4n_velocity_test_traction(
            traction_test_a,
            a,
            gradNi,
            n,
            mu_eff);

        vec_i[a] += dt * (-Ni * traction_u[a]);
        vec_i[a] += dt * (-t4n_dot3(traction_test_a, r));
        vec_i[a] += dt * ( Ni * tau_B * r[a]);
    }

    /*
      Pressure test row:

        - sigma(0,q)n . r
        = -(-q n) . r
        =  q (r . n)

      This is the full-vector Nitsche pressure-row term.
    */
    vec_i[3] += dt * Ni * t4n_dot3(r, n);

    for (int a = 0; a < 3; ++a) {
        double traction_test_a[3];
        t4n_velocity_test_traction(
            traction_test_a,
            a,
            gradNi,
            n,
            mu_eff);

        for (int b = 0; b < 3; ++b) {
            double traction_trial_b[3];
            t4n_velocity_test_traction(
                traction_trial_b,
                b,
                gradNj,
                n,
                mu_eff);

            /*
              d/du_b of -Ni * sigma(u,p)n[a]
            */
            mat_ij[a][b] += dt * (-Ni * traction_trial_b[a]);

            /*
              d/du_b of -sigma(w,0)n . r
            */
            mat_ij[a][b] += dt * (-traction_test_a[b] * Nj);

            /*
              d/du_b of tau_B Ni r[a]
              tau_B is frozen.
            */
            mat_ij[a][b] += dt * (Ni * Nj * tau_B * (a == b ? 1.0 : 0.0));
        }

        /*
          d/dp of -Ni * sigma(u,p)n[a]:

            -Ni * d(-p n[a])/dp = Ni * Nj * n[a]
        */
        mat_ij[a][3] += dt * (Ni * Nj * n[a]);
    }

    /*
      d/du_b of q (r . n)
    */
    for (int b = 0; b < 3; ++b) {
        mat_ij[3][b] += dt * (Ni * Nj * n[b]);
    }
}



/* -------------------------------------------------------------------------- */
/* TET4 / TRI3 symmetric Nitsche wall assembly                                */
/* -------------------------------------------------------------------------- */
void set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    (void)basis_vol;

    if (fe->local_num_nodes != 4 || surf->local_num_nodes != 3) {
        fprintf(stderr,
            "[wall-t4] TET4/TRI3 expected, got fe->local_num_nodes=%d surf->local_num_nodes=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return;
    }

    /*
     * y+ visualization arrays reset.
     * この関数を呼ぶたびに最新の y+ を作り直す。
     */
    const int do_yplus =
        (vals->y_plus != NULL && vals->y_plus_count != NULL);

    if (do_yplus) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            vals->y_plus[i] = 0.0;
            vals->y_plus_count[i] = 0;
        }
    }

    int* owner = NULL;
    int* lface = NULL;
    if (t4n_make_owner_map(surf, fe, &owner, &lface) != 0) {
        fprintf(stderr, "[wall-t4] surface-owner mapping failed\n");
        return;
    }

    const int ns = surf->local_num_nodes;      /* TRI3 = 3 */
    const int nv = fe->local_num_nodes;        /* TET4 = 4 */
    const int np = basis_surf->num_integ_points;
    const double area_factor = t4n_tri_area_factor(basis_surf);

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];

        if (ke < 0) {
            fprintf(stderr,
                "[wall-t4][BUG] owner not found: es=%d owner=%d\n",
                es,
                ke);
            exit(EXIT_FAILURE);
        }

        /*
         * Verify that the surface triangle and owner TET share exactly 3 nodes.
         */
        {
            int common = 0;

            for (int a = 0; a < ns; ++a) {
                const int sgid = surf->conn[es][a];

                for (int b = 0; b < nv; ++b) {
                    if (fe->conn[ke][b] == sgid) {
                        ++common;
                        break;
                    }
                }
            }

            if (common != ns) {
                fprintf(stderr,
                    "[wall-t4][BUG] surface-owner mismatch: es=%d ke=%d common=%d/%d\n",
                    es,
                    ke,
                    common,
                    ns);
                exit(EXIT_FAILURE);
            }
        }

        double xsurf[3][3];
        for (int a = 0; a < ns; ++a) {
            const int gid = surf->conn[es][a];

            xsurf[a][0] = fe->x[gid][0];
            xsurf[a][1] = fe->x[gid][1];
            xsurf[a][2] = fe->x[gid][2];
        }

        double xvol[4][3];
        for (int a = 0; a < nv; ++a) {
            const int gid = fe->conn[ke][a];

            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        int s2v[3];
        if (!t4n_surface_to_volume_nodes(surf, fe, es, ke, s2v)) {
            fprintf(stderr,
                "[wall-t4] cannot map surface node to volume node: es=%d ke=%d\n",
                es,
                ke);
            continue;
        }

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;

        if (!t4n_grad_phys_metric(xvol, dN_dx, J_inv_face, &Jacobian_face)) {
            fprintf(stderr, "[wall-t4] degenerate TET4: ke=%d\n", ke);
            continue;
        }

        /*
         * TET4 volume: physical volume = |detJ| / 6.
         */
        double h_e_vms = cbrt(fabs(Jacobian_face) / 6.0);
        if (h_e_vms <= 1.0e-12) {
            h_e_vms = 1.0e-12;
        }

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < nv; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < ns; ++a) {
            xf[0] += xsurf[a][0];
            xf[1] += xsurf[a][1];
            xf[2] += xsurf[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < np; ++p) {
            double nrm[3];

            t4n_tri3_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nrm);

            const double Jface_raw = t4n_norm3(nrm);
            if (Jface_raw <= 1.0e-300) {
                continue;
            }

            double n[3] = {
                nrm[0] / Jface_raw,
                nrm[1] / Jface_raw,
                nrm[2] / Jface_raw
            };

            if (t4n_dot3(n, svec) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            /*
             * TET4 basis at the surface quadrature point.
             * The TRI3 quadrature shape value attached to surf node a becomes
             * the TET4 shape value of its matched volume local node s2v[a].
             * The opposite TET node remains zero.
             */
            double Nv[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < ns; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double u_q[3]     = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double p_q        = 0.0;

            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3]    = {0.0, 0.0, 0.0};

            for (int a = 0; a < nv; ++a) {
                const int gid = fe->conn[ke][a];

                for (int c = 0; c < 3; ++c) {
                    u_q[c]     += Nv[a] * vals->v[gid][c];
                    u_old_q[c] += Nv[a] * vals->v_old[gid][c];

                    for (int d = 0; d < 3; ++d) {
                        grad_u_q[c][d] += dN_dx[a][d] * vals->v[gid][c];
                    }
                }

                p_q += Nv[a] * vals->p[gid];

                for (int d = 0; d < 3; ++d) {
                    grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                }
            }

            double* grad_u_ptr[3] = {
                grad_u_q[0],
                grad_u_q[1],
                grad_u_q[2]
            };

            double mu_eff_face = mu_molecular;
            double tau_face    = 0.0;
            double tau_c_face  = 0.0;

            BBFE_vms_mu_eff_tau(
                &mu_eff_face,
                &tau_face,
                &tau_c_face,
                J_inv_face,
                Jacobian_face,
                h_e_vms,
                u_q,
                u_old_q,
                grad_u_ptr,
                grad_p_q,
                rho,
                mu_molecular,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);

            double h_n = fabs(
                (xf[0] - xc[0]) * n[0]
              + (xf[1] - xc[1]) * n[1]
              + (xf[2] - xc[2]) * n[2]);

            if (h_n <= 1.0e-12) {
                h_n = 1.0e-12;
            }

            const double y_wall = h_n;
            const double w =
                basis_surf->integ_weight[p] * area_factor * Jface_raw;

            /*
             * ------------------------------------------------------------
             * y+ visualization value at this wall quadrature point
             * ------------------------------------------------------------
             */
            if (do_yplus) {
                double ur[3];
                for (int a = 0; a < 3; ++a) {
                    ur[a] = u_q[a] - Uw[a];
                }

                const double un = t4n_dot3(ur, n);

                double ut[3];
                for (int a = 0; a < 3; ++a) {
                    ut[a] = ur[a] - un * n[a];
                }

                const double ut_norm = sqrt(t4n_dot3(ut, ut));

                const double u_tau =
                    compute_utau_all_yplus(
                        ut_norm,
                        y_wall,
                        rho,
                        mu_molecular);

                double yplus_q = 0.0;

                if (u_tau > 0.0 && rho > 0.0 && mu_molecular > 0.0) {
                    yplus_q = rho * u_tau * y_wall / mu_molecular;
                }

                /*
                 * 求積点 y+ を面上 TRI3 の3節点へ蓄積する。
                 * ここでは単純平均。
                 */
                for (int a = 0; a < ns; ++a) {
                    const int gid = surf->conn[es][a];

                    vals->y_plus[gid] += yplus_q;
                    vals->y_plus_count[gid] += 1;
                }
            }

            /*
             * Original Nitsche residual and matrix assembly.
             */
            for (int i = 0; i < nv; ++i) {
                const int gi = fe->conn[ke][i];

                double vec_i[4];
                double dummy_mat[4][4];
                const double zero_grad[3] = {0.0, 0.0, 0.0};

                t4n_wall_sym_nitsche_local_vecmat(
                    vec_i,
                    dummy_mat,
                    Nv[i],
                    0.0,
                    dN_dx[i],
                    zero_grad,
                    n,
                    u_q,
                    p_q,
                    grad_u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu_molecular,
                    h_n,
                    y_wall,
                    vals->dt,
                    gamma_n);

                for (int a = 0; a < 4; ++a) {
                    monolis->mat.R.B[4 * gi + a] -= w * vec_i[a];
                }

                for (int j = 0; j < nv; ++j) {
                    const int gj = fe->conn[ke][j];

                    double dummy_vec[4];
                    double mat_ij[4][4];

                    t4n_wall_sym_nitsche_local_vecmat(
                        dummy_vec,
                        mat_ij,
                        Nv[i],
                        Nv[j],
                        dN_dx[i],
                        dN_dx[j],
                        n,
                        u_q,
                        p_q,
                        grad_u_q,
                        Uw,
                        rho,
                        mu_eff_face,
                        mu_molecular,
                        h_n,
                        y_wall,
                        vals->dt,
                        gamma_n);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                gi,
                                gj,
                                a,
                                b,
                                w * mat_ij[a][b]);
                        }
                    }
                }
            }
        }
    }

    /*
     * Final average at wall nodes.
     * 非壁面節点は 0.0 にする。
     */
    if (do_yplus) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            if (vals->y_plus_count[i] > 0) {
                vals->y_plus[i] /= (double)vals->y_plus_count[i];
            } else {
                vals->y_plus[i] = 0.0;
            }
        }
    }

    free(owner);
    free(lface);
}

/* -------------------------------------------------------------------------- */
/* HEX/TET dispatch wrapper. HEX body remains unchanged.                      */
/* -------------------------------------------------------------------------- */

void set_wall_face_vecmat_symmetric_nitsche_hex_or_tet(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    if (fe->local_num_nodes == 8) {
        if (surf->local_num_nodes != 4) {
            fprintf(stderr,
                "[wall] HEX8 volume requires QUAD4 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return;
        }

        /* Verified HEX8/QUAD4 implementation. Do not change its body. */
        set_wall_face_vecmat_symmetric_nitsche(
            monolis,
            fe,
            basis_vol,
            surf,
            basis_surf,
            vals,
            rho,
            mu_molecular,
            Uw,
            gamma_n);
        return;
    }

    if (fe->local_num_nodes == 4) {
        if (surf->local_num_nodes != 3) {
            fprintf(stderr,
                "[wall] TET4 volume requires TRI3 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return;
        }

        set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
            monolis,
            fe,
            basis_vol,
            surf,
            basis_surf,
            vals,
            rho,
            mu_molecular,
            Uw,
            gamma_n);
        return;
    }

    fprintf(stderr,
        "[wall] unsupported fe->local_num_nodes=%d. Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);
}

void calc_Cd_p_hex_or_tet(
    BBFE_DATA*  surf,
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    VALUES*     vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out)
{
    if (D_out) *D_out = 0.0;
    if (L_out) *L_out = 0.0;

    if (fe->local_num_nodes == 8) {
        /* HEX8/QUAD4: verified existing implementation. Do not change its body. */
        if (surf->local_num_nodes != 4) {
            fprintf(stderr,
                "[Cd_p] HEX8 volume requires QUAD4 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return;
        }

        calc_Cd_p(
            surf,
            fe,
            basis,
            vals,
            eU_in,
            eP_in,
            D_out,
            L_out);
        return;

    } else if (fe->local_num_nodes == 4) {
        /* TET4/TRI3: new implementation. */
        if (surf->local_num_nodes != 3) {
            fprintf(stderr,
                "[Cd_p] TET4 volume requires TRI3 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return;
        }

        calc_Cd_p_tet4_tri3(
            surf,
            fe,
            basis,
            vals,
            eU_in,
            eP_in,
            D_out,
            L_out);
        return;
    }

    fprintf(stderr,
        "[Cd_p] unsupported fe->local_num_nodes=%d. Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);
}

/*
#define T4_WALL_DIAG_NUM_REGIONS 8

enum {
    T4_WALL_DIAG_FRONT_POS_DRAG = 0,
    T4_WALL_DIAG_REAR_NEG_DRAG  = 1,
    T4_WALL_DIAG_TOP_POS_LIFT   = 2,
    T4_WALL_DIAG_BOTTOM_NEG_LIFT= 3,
    T4_WALL_DIAG_SIDE_POS       = 4,
    T4_WALL_DIAG_SIDE_NEG       = 5,
    T4_WALL_DIAG_SLANT_MIXED    = 6,
    T4_WALL_DIAG_OTHER          = 7
};

static const char* t4_wall_diag_region_name[T4_WALL_DIAG_NUM_REGIONS] = {
    "front_pos_drag",
    "rear_neg_drag",
    "top_pos_lift",
    "bottom_neg_lift",
    "side_pos",
    "side_neg",
    "slant_mixed",
    "other"
};

typedef struct {
    double area;
    double int_n[3];

    double max_abs_un_over_Uref;
    double int_abs_un;
    double int_un;
    double int_un2;
    double mean_abs_un_over_Uref;
    double mean_un_over_Uref;
    double rms_un_over_Uref;
    double Uref_scale;

    double p_min;
    double p_max;
    double p_mean;
    double p_rms;
    double p_area_sum;
    double p2_area_sum;

    double cp_min;
    double cp_max;
    double cp_mean;
    double cp_rms;
    double cp_area_sum;
    double cp2_area_sum;

    // Pressure force in the drag direction: int_Gamma p (n.eU) dGamma.
    //   If eU=(1,0,0), this is int_Gamma p n_x dGamma.
    double int_pnU;
    double area_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_pnU_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_abs_un_region[T4_WALL_DIAG_NUM_REGIONS];
} T4WallDiagnostics;
*/
static double t4wd_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double t4wd_norm3(const double a[3])
{
    return sqrt(t4wd_dot3(a, a));
}

static void t4wd_cross3(double c[3], const double a[3], const double b[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static void t4wd_output_gamma_value(
    const char* prefix,
    const char* component,
    int gamma_index,
    double value,
    double t,
    const char* directory)
{
    char fname[256];

    snprintf(
        fname,
        sizeof(fname),
        "wall_nitsche_%s_%s_gamma_%02d.txt",
        prefix,
        component,
        gamma_index);

    ROM_std_hlpod_output_add_calc_time(value, t, fname, directory);
}

/* -------------------------------------------------------------------------- */
/* TET4/TRI3 viscous + optional Nitsche/wall-law force integral                */
/* -------------------------------------------------------------------------- */

#ifndef T4_NITSCHE_DT_PENALTY_COEFF
#define T4_NITSCHE_DT_PENALTY_COEFF 0.10
#endif

static void t4_cd_nitsche_wall_force_density(
    const double n[3],
    const double u_q[3],
    const double Uw_in[3],
    double rho,
    double mu_eff,
    double mu_wall_law,
    double h_n,
    double y_wall,
    double dt,
    double gamma_n,
    double f_body[3],
    double* beta_n_out,
    double* beta_wall_out)
{
    (void)mu_wall_law;
    (void)y_wall;

    const double eps = 1.0e-30;
    const double dt_eff = fmax(dt, eps);
    const double h_eff  = fmax(h_n, 1.0e-12);

    const double Uw[3] = {
        Uw_in ? Uw_in[0] : 0.0,
        Uw_in ? Uw_in[1] : 0.0,
        Uw_in ? Uw_in[2] : 0.0
    };

    double r[3];
    for (int a = 0; a < 3; ++a) {
        r[a] = u_q[a] - Uw[a];
    }

    const double beta_inv =
        T4_NITSCHE_INV_PENALTY_COEFF * rho * h_eff / dt_eff;
    const double beta_vis =
        T4_NITSCHE_VIS_PENALTY_COEFF * mu_eff / h_eff;
    const double tau_B = gamma_n * fmax(beta_inv, beta_vis);

    /*
      Conservative weak-BC penalty reaction on the body:

        F_penalty_body = - tau_B (u - Uw).

      The stress part is evaluated separately in calc_Cd_v_tet4_tri3_core().
    */
    for (int a = 0; a < 3; ++a) {
        f_body[a] = -tau_B * r[a];
    }

    if (beta_n_out)    *beta_n_out    = tau_B;
    if (beta_wall_out) *beta_wall_out = 0.0;
}
static int calc_Cd_v_tet4_tri3_core(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    int include_wall_nitsche,
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    (void)monolis_com;
    (void)Uinf;
    (void)Aref;

    if (D_out)           *D_out           = 0.0;
    if (L_out)           *L_out           = 0.0;
    if (Cd_out)          *Cd_out          = 0.0;
    if (Cl_out)          *Cl_out          = 0.0;
    if (D_nitsche_out)   *D_nitsche_out   = 0.0;
    if (L_nitsche_out)   *L_nitsche_out   = 0.0;
    if (Cd_nitsche_out)  *Cd_nitsche_out  = 0.0;
    if (Cl_nitsche_out)  *Cl_nitsche_out  = 0.0;

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        return -1;
    }

    double (*G_all)[3][3] =
        (double(*)[3][3])malloc(sizeof(double) * fe->total_num_elems * 9);
    if (!G_all) {
        free(owner);
        free(lface);
        return -1;
    }
    t4_compute_all_gradients_tet4_exact(fe, vals, G_all);

    const double area_factor = t4_tri_area_factor(basis_surf);

    double F_vis[3] = {0.0, 0.0, 0.0};
    double F_nit[3] = {0.0, 0.0, 0.0};

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xvol[4][3];
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        int s2v[3];
        if (!t4_surface_to_volume_nodes(surf, fe, es, ke, s2v)) {
            fprintf(stderr,
                "[Cd_v_t4] cannot map surface node to volume node: es=%d ke=%d\n",
                es,
                ke);
            continue;
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        const double (*G)[3] = G_all[ke];
        const double divu = G[0][0] + G[1][1] + G[2][2];

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;
        const int have_metric =
            t4_grad_phys_metric(xvol, dN_dx, J_inv_face, &Jacobian_face);

        double h_e_vms = 1.0e-12;
        if (have_metric) {
            h_e_vms = cbrt(fabs(Jacobian_face) / 6.0);
            if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;
        }

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double Nv[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double p_q = 0.0;
            double u_q[3] = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 4; ++a) {
                const int gid = fe->conn[ke][a];

                p_q += Nv[a] * vals->p[gid];

                for (int c = 0; c < 3; ++c) {
                    u_q[c] += Nv[a] * vals->v[gid][c];

                    if (vals->v_old != NULL) {
                        u_old_q[c] += Nv[a] * vals->v_old[gid][c];
                    } else {
                        u_old_q[c] += Nv[a] * vals->v[gid][c];
                    }
                }

                if (have_metric) {
                    for (int c = 0; c < 3; ++c) {
                        for (int d = 0; d < 3; ++d) {
                            grad_u_q[c][d] +=
                                dN_dx[a][d] * vals->v[gid][c];
                        }
                    }

                    for (int d = 0; d < 3; ++d) {
                        grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                    }
                }
            }

            (void)p_q;

            /*
              Viscous traction consistent with the volume bilinear form

                  mu grad(u):grad(w) + mu div(u) div(w).

              The corresponding boundary traction is

                  t_vis = mu grad(u) n + mu n div(u).

              If the physical Cauchy stress is required instead, replace this
              block by mu(grad u + grad u^T)n - 2/3 mu div(u)n.
            */
            double tau_n[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < 3; ++i) {
                double gradun_i = 0.0;
                for (int j = 0; j < 3; ++j) {
                    gradun_i += G[i][j] * n[j];
                }
                tau_n[i] = mu * gradun_i + mu * divu * n[i];
            }

            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;

            F_vis[0] += tau_n[0] * w;
            F_vis[1] += tau_n[1] * w;
            F_vis[2] += tau_n[2] * w;

            if (include_wall_nitsche && have_metric) {
                double* grad_u_ptr[3] = {
                    grad_u_q[0],
                    grad_u_q[1],
                    grad_u_q[2]
                };

                double mu_eff_face = mu;
                double tau_face = 0.0;
                double tau_c_face = 0.0;

                BBFE_vms_mu_eff_tau(
                    &mu_eff_face,
                    &tau_face,
                    &tau_c_face,
                    J_inv_face,
                    Jacobian_face,
                    h_e_vms,
                    u_q,
                    u_old_q,
                    grad_u_ptr,
                    grad_p_q,
                    rho,
                    mu,
                    vals->dt,
                    vals->C_vms,
                    vals->vms_cap_coeff);

                double h_n = fabs(
                    (xf[0] - xc[0]) * n[0]
                  + (xf[1] - xc[1]) * n[1]
                  + (xf[2] - xc[2]) * n[2]);

                if (h_n <= 1.0e-12) h_n = 1.0e-12;

                const double y_wall = h_n;

                double f_nit_q[3];
                t4_cd_nitsche_wall_force_density(
                    n,
                    u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu,
                    h_n,
                    y_wall,
                    vals->dt,
                    gamma_n,
                    f_nit_q,
                    NULL,
                    NULL);

                F_nit[0] += f_nit_q[0] * w;
                F_nit[1] += f_nit_q[1] * w;
                F_nit[2] += f_nit_q[2] * w;
            }
        }
    }

    double eU[3] = { eU_in[0], eU_in[1], eU_in[2] };
    double eP[3] = { eP_in[0], eP_in[1], eP_in[2] };
    t4_normalize3(eU);

    const double ep_dot_eu = t4_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];
    t4_normalize3(eP);

    /*
      Backward-compatible convention:
      D_out/L_out and Cd_out/Cl_out still contain dimensional forces.
      Nondimensionalization should be done once, consistently with pressure
      and Nitsche contributions, at the caller side.
    */
    const double F_body_vis[3] = {
        -F_vis[0],
        -F_vis[1],
        -F_vis[2]
    };

    const double D = t4_dot3(F_body_vis, eU);
    const double L = t4_dot3(F_body_vis, eP);
    const double D_nit = t4_dot3(F_nit, eU);
    const double L_nit = t4_dot3(F_nit, eP);

    if (D_out)          *D_out          = D;
    if (L_out)          *L_out          = L;
    if (Cd_out)         *Cd_out         = D;
    if (Cl_out)         *Cl_out         = L;
    if (D_nitsche_out)  *D_nitsche_out  = D_nit;
    if (L_nitsche_out)  *L_nitsche_out  = L_nit;
    if (Cd_nitsche_out) *Cd_nitsche_out = D_nit;
    if (Cl_nitsche_out) *Cl_nitsche_out = L_nit;

    free(G_all);
    free(owner);
    free(lface);
    return 0;
}

/*
 * Backward-compatible entry point: original viscous force only.
 * Use calc_Cd_v_tet4_tri3_nitsche() when the wall is imposed by Nitsche
 * + all-y+ wall law and you want the wall-reaction diagnostic returned
 * separately.
 */
int calc_Cd_v_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out)
{
    return calc_Cd_v_tet4_tri3_core(
        surf,
        fe,
        basis_surf,
        monolis_com,
        vals,
        rho,
        mu,
        Uinf,
        Aref,
        eU_in,
        eP_in,
        0,
        NULL,
        0.0,
        D_out,
        L_out,
        Cd_out,
        Cl_out,
        NULL,
        NULL,
        NULL,
        NULL);
}

/*
 * Extended entry point.
 *
 * D_out/L_out and Cd_out/Cl_out are the same viscous-only quantities as the
 * verified calc_Cd_v() path.  D_nitsche_out/L_nitsche_out and
 * Cd_nitsche_out/Cl_nitsche_out are separate diagnostics for the normal
 * penalty plus tangential wall-law reaction.
 */
int calc_Cd_v_tet4_tri3_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    return calc_Cd_v_tet4_tri3_core(
        surf,
        fe,
        basis_surf,
        monolis_com,
        vals,
        rho,
        mu,
        Uinf,
        Aref,
        eU_in,
        eP_in,
        1,
        Uw,
        gamma_n,
        D_out,
        L_out,
        Cd_out,
        Cl_out,
        D_nitsche_out,
        L_nitsche_out,
        Cd_nitsche_out,
        Cl_nitsche_out);
}


/*
 * Dispatch wrapper including TET4/TRI3 Nitsche wall reaction.
 * For HEX8/QUAD4 this currently falls back to the verified existing calc_Cd_v().
 * For TET4/TRI3 it calls calc_Cd_v_tet4_tri3_nitsche().
 */
int calc_Cd_v_hex_or_tet_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    if (D_out)          *D_out          = 0.0;
    if (L_out)          *L_out          = 0.0;
    if (Cd_out)         *Cd_out         = 0.0;
    if (Cl_out)         *Cl_out         = 0.0;
    if (D_nitsche_out)  *D_nitsche_out  = 0.0;
    if (L_nitsche_out)  *L_nitsche_out  = 0.0;
    if (Cd_nitsche_out) *Cd_nitsche_out = 0.0;
    if (Cl_nitsche_out) *Cl_nitsche_out = 0.0;

    if (fe->local_num_nodes == 8) {
        if (surf->local_num_nodes != 4) {
            fprintf(stderr,
                "[Cd_v_nitsche] HEX8 volume requires QUAD4 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return -1;
        }

        /* HEX Nitsche force is not added here because the verified HEX Cd_v()
         * interface has no Uw/gamma_n. Add an analogous HEX-specific extended
         * routine if you also want HEX wall reaction diagnostics.
         */
        return calc_Cd_v(
            surf,
            fe,
            basis,
            monolis_com,
            vals,
            rho,
            mu,
            Uinf,
            Aref,
            eU_in,
            eP_in,
            D_out,
            L_out,
            Cd_out,
            Cl_out);
    }

    if (fe->local_num_nodes == 4) {
        if (surf->local_num_nodes != 3) {
            fprintf(stderr,
                "[Cd_v_nitsche] TET4 volume requires TRI3 surface, but surf->local_num_nodes=%d\n",
                surf->local_num_nodes);
            return -1;
        }

        return calc_Cd_v_tet4_tri3_nitsche(
            surf,
            fe,
            basis,
            monolis_com,
            vals,
            rho,
            mu,
            Uinf,
            Aref,
            eU_in,
            eP_in,
            Uw,
            gamma_n,
            D_out,
            L_out,
            Cd_out,
            Cl_out,
            D_nitsche_out,
            L_nitsche_out,
            Cd_nitsche_out,
            Cl_nitsche_out);
    }

    fprintf(stderr,
        "[Cd_v_nitsche] unsupported fe->local_num_nodes=%d. Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);

    return -1;
}

static void t4wd_normalize3(double a[3])
{
    const double n = t4wd_norm3(a);
    if (n > 1.0e-300) {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}

static void t4wd_prepare_dirs(
    double eU[3],
    double eP[3],
    double eS[3],
    const double eU_in[3],
    const double eP_in[3])
{
    eU[0] = eU_in[0];
    eU[1] = eU_in[1];
    eU[2] = eU_in[2];
    t4wd_normalize3(eU);

    eP[0] = eP_in[0];
    eP[1] = eP_in[1];
    eP[2] = eP_in[2];

    const double ep_dot_eu = t4wd_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];

    if (t4wd_norm3(eP) <= 1.0e-14) {
        eP[0] = 0.0;
        eP[1] = 1.0;
        eP[2] = 0.0;
        const double dot_tmp = t4wd_dot3(eP, eU);
        eP[0] -= dot_tmp * eU[0];
        eP[1] -= dot_tmp * eU[1];
        eP[2] -= dot_tmp * eU[2];
    }
    t4wd_normalize3(eP);

    t4wd_cross3(eS, eU, eP);
    t4wd_normalize3(eS);
}

static int t4wd_classify_region(
    const double n[3],
    const double eU[3],
    const double eP[3],
    const double eS[3])
{
    const double du = t4wd_dot3(n, eU);
    const double dp = t4wd_dot3(n, eP);
    const double ds = t4wd_dot3(n, eS);
    const double au = fabs(du);
    const double ap = fabs(dp);
    const double as = fabs(ds);
    const double axis_threshold = 0.70;

    if (au >= ap && au >= as && au >= axis_threshold) {
        return (du >= 0.0) ? T4_WALL_DIAG_FRONT_POS_DRAG
                           : T4_WALL_DIAG_REAR_NEG_DRAG;
    }
    if (ap >= au && ap >= as && ap >= axis_threshold) {
        return (dp >= 0.0) ? T4_WALL_DIAG_TOP_POS_LIFT
                           : T4_WALL_DIAG_BOTTOM_NEG_LIFT;
    }
    if (as >= au && as >= ap && as >= axis_threshold) {
        return (ds >= 0.0) ? T4_WALL_DIAG_SIDE_POS
                           : T4_WALL_DIAG_SIDE_NEG;
    }

    return T4_WALL_DIAG_SLANT_MIXED;
}

static void t4_wall_diag_init(T4WallDiagnostics* diag)
{
    diag->area = 0.0;
    diag->int_n[0] = 0.0;
    diag->int_n[1] = 0.0;
    diag->int_n[2] = 0.0;

    diag->max_abs_un_over_Uref = 0.0;
    diag->int_abs_un = 0.0;
    diag->int_un = 0.0;
    diag->int_un2 = 0.0;
    diag->mean_abs_un_over_Uref = 0.0;
    diag->mean_un_over_Uref = 0.0;
    diag->rms_un_over_Uref = 0.0;
    diag->Uref_scale = 1.0;

    diag->p_min = 1.0e300;
    diag->p_max = -1.0e300;
    diag->p_mean = 0.0;
    diag->p_rms = 0.0;
    diag->p_area_sum = 0.0;
    diag->p2_area_sum = 0.0;

    diag->cp_min = 1.0e300;
    diag->cp_max = -1.0e300;
    diag->cp_mean = 0.0;
    diag->cp_rms = 0.0;
    diag->cp_area_sum = 0.0;
    diag->cp2_area_sum = 0.0;

    diag->int_pnU = 0.0;

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->area_region[r] = 0.0;
        diag->int_pnU_region[r] = 0.0;
        diag->int_abs_un_region[r] = 0.0;
    }
}

static void t4_wall_diag_finalize(T4WallDiagnostics* diag)
{
    if (diag->area > 0.0) {
        diag->p_mean = diag->p_area_sum / diag->area;
        diag->p_rms = sqrt(fmax(diag->p2_area_sum / diag->area, 0.0));

        diag->cp_mean = diag->cp_area_sum / diag->area;
        diag->cp_rms = sqrt(fmax(diag->cp2_area_sum / diag->area, 0.0));

        diag->mean_abs_un_over_Uref =
            diag->int_abs_un / (diag->Uref_scale * diag->area);
        diag->mean_un_over_Uref =
            diag->int_un / (diag->Uref_scale * diag->area);
        diag->rms_un_over_Uref =
            sqrt(fmax(diag->int_un2 / diag->area, 0.0)) / diag->Uref_scale;
    } else {
        diag->p_min = 0.0;
        diag->p_max = 0.0;
        diag->cp_min = 0.0;
        diag->cp_max = 0.0;
        diag->mean_abs_un_over_Uref = 0.0;
        diag->mean_un_over_Uref = 0.0;
        diag->rms_un_over_Uref = 0.0;
    }
}

static void t4_wall_penalty_scale_finalize(
    T4WallPenaltyScaleDiagnostics* diag)
{
    if (diag->area > 0.0) {
        diag->h_mean = diag->h_sum / diag->area;
    } else {
        diag->h_mean = 0.0;
        diag->h_min = 0.0;
        diag->h_max = 0.0;
    }
}

int calc_wall_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    T4WallDiagnostics* diag)
{
    t4_wall_diag_init(diag);

    if (surf->local_num_nodes != 3 || fe->local_num_nodes != 4) {
        fprintf(stderr,
            "[wall_diag] TET4/TRI3 expected, got fe=%d surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return -1;
    }

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        fprintf(stderr, "[wall_diag] surface-owner mapping failed\n");
        return -1;
    }

    double eU[3], eP[3], eS[3];
    t4wd_prepare_dirs(eU, eP, eS, eU_in, eP_in);

    const double U_scale = fmax(fabs(Uref), 1.0e-300);
    diag->Uref_scale = U_scale;

    double qref = 0.5 * rho * Uref * Uref;
    if (fabs(qref) <= 1.0e-300) qref = 1.0;

    const double area_factor = t4_tri_area_factor(basis_surf);

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xc[0] += fe->x[gid][0];
            xc[1] += fe->x[gid][1];
            xc[2] += fe->x[gid][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4wd_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4wd_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double p_q = 0.0;
            double u_q[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 3; ++a) {
                const int gid = surf->conn[es][a];
                const double Na = basis_surf->N[p][a];
                p_q += Na * vals->p[gid];
                u_q[0] += Na * vals->v[gid][0];
                u_q[1] += Na * vals->v[gid][1];
                u_q[2] += Na * vals->v[gid][2];
            }

            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;
            const double un = t4wd_dot3(u_q, n);
            const double abs_un = fabs(un);
            const double cp = (p_q - p_ref) / qref;
            const double pnU = p_q * t4wd_dot3(n, eU);
            const int region = t4wd_classify_region(n, eU, eP, eS);

            diag->area += w;
            diag->int_n[0] += n[0] * w;
            diag->int_n[1] += n[1] * w;
            diag->int_n[2] += n[2] * w;

            diag->int_abs_un += abs_un * w;
            diag->int_un += un * w;
            diag->int_un2 += un * un * w;
            diag->max_abs_un_over_Uref =
                fmax(diag->max_abs_un_over_Uref, abs_un / U_scale);

            diag->p_min = fmin(diag->p_min, p_q);
            diag->p_max = fmax(diag->p_max, p_q);
            diag->p_area_sum += p_q * w;
            diag->p2_area_sum += p_q * p_q * w;

            diag->cp_min = fmin(diag->cp_min, cp);
            diag->cp_max = fmax(diag->cp_max, cp);
            diag->cp_area_sum += cp * w;
            diag->cp2_area_sum += cp * cp * w;

            diag->int_pnU += pnU * w;
            diag->area_region[region] += w;
            diag->int_pnU_region[region] += pnU * w;
            diag->int_abs_un_region[region] += abs_un * w;
        }
    }

    free(owner);
    free(lface);

    t4_wall_diag_finalize(diag);
    return 0;
}

void wall_diagnostics_allreduce(
    T4WallDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    double sum_pack[64];
    int k = 0;

    sum_pack[k++] = diag->area;
    sum_pack[k++] = diag->int_abs_un;
    sum_pack[k++] = diag->int_un;
    sum_pack[k++] = diag->int_un2;
    sum_pack[k++] = diag->int_pnU;
    sum_pack[k++] = diag->p_area_sum;
    sum_pack[k++] = diag->p2_area_sum;
    sum_pack[k++] = diag->cp_area_sum;
    sum_pack[k++] = diag->cp2_area_sum;
    sum_pack[k++] = diag->int_n[0];
    sum_pack[k++] = diag->int_n[1];
    sum_pack[k++] = diag->int_n[2];

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->area_region[r];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->int_pnU_region[r];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->int_abs_un_region[r];
    }

    monolis_allreduce_R(k, sum_pack, MONOLIS_MPI_SUM, monolis_com->comm);

    k = 0;
    diag->area = sum_pack[k++];
    diag->int_abs_un = sum_pack[k++];
    diag->int_un = sum_pack[k++];
    diag->int_un2 = sum_pack[k++];
    diag->int_pnU = sum_pack[k++];
    diag->p_area_sum = sum_pack[k++];
    diag->p2_area_sum = sum_pack[k++];
    diag->cp_area_sum = sum_pack[k++];
    diag->cp2_area_sum = sum_pack[k++];
    diag->int_n[0] = sum_pack[k++];
    diag->int_n[1] = sum_pack[k++];
    diag->int_n[2] = sum_pack[k++];

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->area_region[r] = sum_pack[k++];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->int_pnU_region[r] = sum_pack[k++];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->int_abs_un_region[r] = sum_pack[k++];
    }

    {
        double max_pack[3] = {
            diag->max_abs_un_over_Uref,
            diag->p_max,
            diag->cp_max
        };
        monolis_allreduce_R(3, max_pack, MONOLIS_MPI_MAX, monolis_com->comm);
        diag->max_abs_un_over_Uref = max_pack[0];
        diag->p_max = max_pack[1];
        diag->cp_max = max_pack[2];
    }

    {
        double min_pack_as_negative[2] = {
            -diag->p_min,
            -diag->cp_min
        };
        monolis_allreduce_R(
            2,
            min_pack_as_negative,
            MONOLIS_MPI_MAX,
            monolis_com->comm);
        diag->p_min = -min_pack_as_negative[0];
        diag->cp_min = -min_pack_as_negative[1];
    }

    t4_wall_diag_finalize(diag);
}

void output_wall_diagnostics_tet4_tri3(
    const T4WallDiagnostics* diag,
    double t,
    const char* directory)
{
    ROM_std_hlpod_output_add_calc_time(
        diag->max_abs_un_over_Uref,
        t,
        "wall_max_abs_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_abs_un_over_Uref,
        t,
        "wall_mean_abs_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_un_over_Uref,
        t,
        "wall_mean_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->rms_un_over_Uref,
        t,
        "wall_rms_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->area,
        t,
        "wall_area_total.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_abs_un,
        t,
        "wall_int_abs_un.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_un,
        t,
        "wall_int_un_signed.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_min,
        t,
        "wall_Cp_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_max,
        t,
        "wall_Cp_max.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_rms,
        t,
        "wall_Cp_rms.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_mean,
        t,
        "wall_Cp_mean.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_min,
        t,
        "wall_p_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_max,
        t,
        "wall_p_max.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_rms,
        t,
        "wall_p_rms.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_pnU,
        t,
        "wall_pressure_drag_int_pnU.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[0],
        t,
        "wall_int_nx.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[1],
        t,
        "wall_int_ny.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[2],
        t,
        "wall_int_nz.txt",
        directory);

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        char fname[256];

        snprintf(
            fname,
            sizeof(fname),
            "wall_pressure_drag_%s.txt",
            t4_wall_diag_region_name[r]);
        ROM_std_hlpod_output_add_calc_time(
            diag->int_pnU_region[r],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_area_%s.txt",
            t4_wall_diag_region_name[r]);
        ROM_std_hlpod_output_add_calc_time(
            diag->area_region[r],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_int_abs_un_%s.txt",
            t4_wall_diag_region_name[r]);
        ROM_std_hlpod_output_add_calc_time(
            diag->int_abs_un_region[r],
            t,
            fname,
            directory);
    }
}

int calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    T4WallPenaltyScaleDiagnostics* scale_diag)
{
    for (int g = 0; g < num_gamma; ++g) {
        if (D_penalty_by_gamma != NULL) D_penalty_by_gamma[g] = 0.0;
        if (L_penalty_by_gamma != NULL) L_penalty_by_gamma[g] = 0.0;
        if (D_tangent_by_gamma != NULL) D_tangent_by_gamma[g] = 0.0;
        if (L_tangent_by_gamma != NULL) L_tangent_by_gamma[g] = 0.0;
        if (D_total_by_gamma != NULL) D_total_by_gamma[g] = 0.0;
        if (L_total_by_gamma != NULL) L_total_by_gamma[g] = 0.0;
        if (D_penalty_mu_by_gamma != NULL) D_penalty_mu_by_gamma[g] = 0.0;
        if (L_penalty_mu_by_gamma != NULL) L_penalty_mu_by_gamma[g] = 0.0;
        if (D_penalty_dt_by_gamma != NULL) D_penalty_dt_by_gamma[g] = 0.0;
        if (L_penalty_dt_by_gamma != NULL) L_penalty_dt_by_gamma[g] = 0.0;
    }

    if (scale_diag != NULL) {
        t4_wall_penalty_scale_init(scale_diag);
    }

    if (surf->local_num_nodes != 3 || fe->local_num_nodes != 4) {
        fprintf(stderr,
            "[wall_gamma] TET4/TRI3 expected, got fe=%d surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return -1;
    }

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        fprintf(stderr, "[wall_gamma] surface-owner mapping failed\n");
        return -1;
    }

    double eU[3], eP[3], eS[3];
    t4wd_prepare_dirs(eU, eP, eS, eU_in, eP_in);
    (void)eS;

    const double area_factor = t4_tri_area_factor(basis_surf);

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xvol[4][3];
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        int s2v[3];
        if (!t4_surface_to_volume_nodes(surf, fe, es, ke, s2v)) {
            fprintf(stderr,
                "[wall_gamma] cannot map surface node to volume node: es=%d ke=%d\n",
                es,
                ke);
            continue;
        }

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;
        const int have_metric =
            t4_grad_phys_metric(xvol, dN_dx, J_inv_face, &Jacobian_face);

        double h_e_vms = 1.0e-12;
        if (have_metric) {
            h_e_vms = cbrt(fabs(Jacobian_face) / 6.0);
            if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4wd_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4wd_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double Nv[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double u_q[3] = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 4; ++a) {
                const int gid = fe->conn[ke][a];

                for (int c = 0; c < 3; ++c) {
                    u_q[c] += Nv[a] * vals->v[gid][c];

                    if (vals->v_old != NULL) {
                        u_old_q[c] += Nv[a] * vals->v_old[gid][c];
                    } else {
                        u_old_q[c] += Nv[a] * vals->v[gid][c];
                    }
                }

                if (have_metric) {
                    for (int c = 0; c < 3; ++c) {
                        for (int d = 0; d < 3; ++d) {
                            grad_u_q[c][d] +=
                                dN_dx[a][d] * vals->v[gid][c];
                        }
                    }

                    for (int d = 0; d < 3; ++d) {
                        grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                    }
                }
            }

            double mu_eff_face = mu;
            if (have_metric) {
                double* grad_u_ptr[3] = {
                    grad_u_q[0],
                    grad_u_q[1],
                    grad_u_q[2]
                };

                double tau_face = 0.0;
                double tau_c_face = 0.0;

                BBFE_vms_mu_eff_tau(
                    &mu_eff_face,
                    &tau_face,
                    &tau_c_face,
                    J_inv_face,
                    Jacobian_face,
                    h_e_vms,
                    u_q,
                    u_old_q,
                    grad_u_ptr,
                    grad_p_q,
                    rho,
                    mu,
                    vals->dt,
                    vals->C_vms,
                    vals->vms_cap_coeff);
            }

            double ur[3] = {
                u_q[0] - Uw[0],
                u_q[1] - Uw[1],
                u_q[2] - Uw[2]
            };

            const double un = t4wd_dot3(ur, n);
            double ut[3] = {
                ur[0] - un * n[0],
                ur[1] - un * n[1],
                ur[2] - un * n[2]
            };
            const double ut_norm = t4wd_norm3(ut);

            double h_n = fabs(
                (xf[0] - xc[0]) * n[0]
              + (xf[1] - xc[1]) * n[1]
              + (xf[2] - xc[2]) * n[2]);
            if (h_n <= 1.0e-12) h_n = 1.0e-12;

            const double y_wall = h_n;
            const double dt_eff = fmax(vals->dt, 1.0e-30);
            const double beta_mu_base = mu_eff_face / h_n;
            const double beta_dt_base =
                T4_NITSCHE_DT_PENALTY_COEFF * rho * h_n / dt_eff;
            const double beta_base = beta_mu_base + beta_dt_base;

            double beta_wall = mu / y_wall;
            if (ut_norm > 1.0e-14) {
                const double u_tau =
                    compute_utau_all_yplus(ut_norm, y_wall, rho, mu);
                if (u_tau > 0.0) {
                    beta_wall = rho * u_tau * u_tau / ut_norm;
                }
            }

            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;

            if (scale_diag != NULL) {
                scale_diag->area += w;
                scale_diag->h_sum += h_n * w;
                scale_diag->h_min = fmin(scale_diag->h_min, h_n);
                scale_diag->h_max = fmax(scale_diag->h_max, h_n);
            }

            double f_tangent[3];
            for (int a = 0; a < 3; ++a) {
                f_tangent[a] = -beta_wall * ut[a];
            }

            const double D_tangent = t4wd_dot3(f_tangent, eU) * w;
            const double L_tangent = t4wd_dot3(f_tangent, eP) * w;

            for (int g = 0; g < num_gamma; ++g) {
                const double beta_mu = gamma_list[g] * beta_mu_base;
                const double beta_dt = gamma_list[g] * beta_dt_base;
                const double beta_n = gamma_list[g] * beta_base;
                double f_penalty[3];
                double f_penalty_mu[3];
                double f_penalty_dt[3];

                for (int a = 0; a < 3; ++a) {
                    f_penalty_mu[a] = -beta_mu * un * n[a];
                    f_penalty_dt[a] = -beta_dt * un * n[a];
                    f_penalty[a] = -beta_n * un * n[a];
                }

                const double D_penalty = t4wd_dot3(f_penalty, eU) * w;
                const double L_penalty = t4wd_dot3(f_penalty, eP) * w;
                const double D_penalty_mu =
                    t4wd_dot3(f_penalty_mu, eU) * w;
                const double L_penalty_mu =
                    t4wd_dot3(f_penalty_mu, eP) * w;
                const double D_penalty_dt =
                    t4wd_dot3(f_penalty_dt, eU) * w;
                const double L_penalty_dt =
                    t4wd_dot3(f_penalty_dt, eP) * w;

                if (D_penalty_by_gamma != NULL) {
                    D_penalty_by_gamma[g] += D_penalty;
                }
                if (L_penalty_by_gamma != NULL) {
                    L_penalty_by_gamma[g] += L_penalty;
                }
                if (D_tangent_by_gamma != NULL) {
                    D_tangent_by_gamma[g] += D_tangent;
                }
                if (L_tangent_by_gamma != NULL) {
                    L_tangent_by_gamma[g] += L_tangent;
                }
                if (D_total_by_gamma != NULL) {
                    D_total_by_gamma[g] += D_penalty + D_tangent;
                }
                if (L_total_by_gamma != NULL) {
                    L_total_by_gamma[g] += L_penalty + L_tangent;
                }
                if (D_penalty_mu_by_gamma != NULL) {
                    D_penalty_mu_by_gamma[g] += D_penalty_mu;
                }
                if (L_penalty_mu_by_gamma != NULL) {
                    L_penalty_mu_by_gamma[g] += L_penalty_mu;
                }
                if (D_penalty_dt_by_gamma != NULL) {
                    D_penalty_dt_by_gamma[g] += D_penalty_dt;
                }
                if (L_penalty_dt_by_gamma != NULL) {
                    L_penalty_dt_by_gamma[g] += L_penalty_dt;
                }
            }
        }
    }

    free(owner);
    free(lface);

    if (scale_diag != NULL) {
        t4_wall_penalty_scale_finalize(scale_diag);
    }

    return 0;
}


int calc_wall_nitsche_gamma_sweep_split_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma)
{
    return calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
        surf,
        fe,
        basis_surf,
        vals,
        rho,
        mu,
        Uw,
        eU_in,
        eP_in,
        gamma_list,
        num_gamma,
        D_penalty_by_gamma,
        L_penalty_by_gamma,
        D_tangent_by_gamma,
        L_tangent_by_gamma,
        D_total_by_gamma,
        L_total_by_gamma,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL);
}

int calc_wall_nitsche_gamma_sweep_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_nitsche_by_gamma,
    double* L_nitsche_by_gamma)
{
    return calc_wall_nitsche_gamma_sweep_split_tet4_tri3(
        surf,
        fe,
        basis_surf,
        vals,
        rho,
        mu,
        Uw,
        eU_in,
        eP_in,
        gamma_list,
        num_gamma,
        NULL,
        NULL,
        NULL,
        NULL,
        D_nitsche_by_gamma,
        L_nitsche_by_gamma);
}

void wall_gamma_sweep_allreduce(
    int num_gamma,
    double* D_nitsche_by_gamma,
    double* L_nitsche_by_gamma,
    MONOLIS_COM* monolis_com)
{
    monolis_allreduce_R(
        num_gamma,
        D_nitsche_by_gamma,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        num_gamma,
        L_nitsche_by_gamma,
        MONOLIS_MPI_SUM,
        monolis_com->comm);
}

void wall_gamma_sweep_split_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    MONOLIS_COM* monolis_com)
{
    if (D_penalty_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_tangent_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_tangent_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_tangent_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_tangent_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_total_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_total_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_total_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_total_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
}

void wall_gamma_sweep_decomposed_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    MONOLIS_COM* monolis_com)
{
    wall_gamma_sweep_split_allreduce(
        num_gamma,
        D_penalty_by_gamma,
        L_penalty_by_gamma,
        D_tangent_by_gamma,
        L_tangent_by_gamma,
        D_total_by_gamma,
        L_total_by_gamma,
        monolis_com);

    if (D_penalty_mu_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_mu_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_mu_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_mu_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_penalty_dt_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_dt_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_dt_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_dt_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
}

void wall_penalty_scale_diagnostics_allreduce(
    T4WallPenaltyScaleDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    double sum_pack[2] = {
        diag->area,
        diag->h_sum
    };

    monolis_allreduce_R(
        2,
        sum_pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    diag->area = sum_pack[0];
    diag->h_sum = sum_pack[1];

    {
        double max_pack[1] = {diag->h_max};
        monolis_allreduce_R(
            1,
            max_pack,
            MONOLIS_MPI_MAX,
            monolis_com->comm);
        diag->h_max = max_pack[0];
    }

    {
        double min_pack_as_negative[1] = {-diag->h_min};
        monolis_allreduce_R(
            1,
            min_pack_as_negative,
            MONOLIS_MPI_MAX,
            monolis_com->comm);
        diag->h_min = -min_pack_as_negative[0];
    }

    t4_wall_penalty_scale_finalize(diag);
}

void output_wall_penalty_scale_diagnostics_tet4_tri3(
    const T4WallPenaltyScaleDiagnostics* diag,
    double t,
    const char* directory)
{
    ROM_std_hlpod_output_add_calc_time(
        T4_NITSCHE_DT_PENALTY_COEFF,
        t,
        "wall_nitsche_dt_penalty_coeff.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_mean,
        t,
        "wall_h_n_mean.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_min,
        t,
        "wall_h_n_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_max,
        t,
        "wall_h_n_max.txt",
        directory);
}

void output_wall_gamma_sweep_penalty_decomposition_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_mu_by_gamma,
    const double* L_penalty_mu_by_gamma,
    const double* D_penalty_dt_by_gamma,
    const double* L_penalty_dt_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        if (D_penalty_mu_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_mu",
                "drag",
                g,
                D_penalty_mu_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_mu_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_mu",
                "lift",
                g,
                L_penalty_mu_by_gamma[g],
                t,
                directory);
        }
        if (D_penalty_dt_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_dt",
                "drag",
                g,
                D_penalty_dt_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_dt_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_dt",
                "lift",
                g,
                L_penalty_dt_by_gamma[g],
                t,
                directory);
        }

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}

void output_wall_gamma_sweep_split_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_by_gamma,
    const double* L_penalty_by_gamma,
    const double* D_tangent_by_gamma,
    const double* L_tangent_by_gamma,
    const double* D_total_by_gamma,
    const double* L_total_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        if (D_penalty_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty",
                "drag",
                g,
                D_penalty_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty",
                "lift",
                g,
                L_penalty_by_gamma[g],
                t,
                directory);
        }
        if (D_tangent_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "tangent",
                "drag",
                g,
                D_tangent_by_gamma[g],
                t,
                directory);
        }
        if (L_tangent_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "tangent",
                "lift",
                g,
                L_tangent_by_gamma[g],
                t,
                directory);
        }
        if (D_total_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "total",
                "drag",
                g,
                D_total_by_gamma[g],
                t,
                directory);
        }
        if (L_total_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "total",
                "lift",
                g,
                L_total_by_gamma[g],
                t,
                directory);
        }

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}

void output_wall_gamma_sweep_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_nitsche_by_gamma,
    const double* L_nitsche_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_drag_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            D_nitsche_by_gamma[g],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_lift_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            L_nitsche_by_gamma[g],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}





