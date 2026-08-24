/*
 * PRI6/TRI3 Karman force diagnostics.
 *
 * Include this file once at the END of core_Nitsche_prism.c so that it can
 * reuse the PRI6 helper routines already defined in that translation unit.
 */

static void pri6_force_diag_zero(PRI6ForceDiagnostics* d)
{
    if(d != NULL) memset(d, 0, sizeof(*d));
}

static int pri6_force_diag_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    PRI6ForceDiagnostics* diag)
{
    pri6_force_diag_zero(diag);

    if(surf == NULL || fe_pri == NULL || basis_surf == NULL ||
       vals == NULL || Uw == NULL || opt == NULL || diag == NULL) {
        return -1;
    }

    if(surf->total_num_elems == 0) return 0;
    if(fe_pri->local_num_nodes != 6 || surf->local_num_nodes != 3) return -1;

    int *owner = NULL, *lface = NULL, unmapped = 0;
    if(BBFE_pri6_tri3_make_owner_map(
            surf, fe_pri, &owner, &lface, &unmapped) != 0) {
        free(owner); free(lface);
        return -1;
    }

    if(unmapped != 0) {
        free(owner); free(lface);
        return -1;
    }

    const int np = basis_surf->num_integ_points;
    const double area_factor = pri6_tri_area_factor(basis_surf);
    if(np <= 0 || area_factor <= 0.0) {
        free(owner); free(lface);
        return -1;
    }

    for(int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        const int lf = lface[es];
        if(ke < 0) continue;

        int s2v[3];
        if(!pri6_surface_to_volume_nodes(surf, fe_pri, es, ke, s2v)) {
            free(owner); free(lface);
            return -1;
        }

        double xvol[6][3], xsurf[3][3];
        for(int a = 0; a < 6; ++a) {
            const int gid = fe_pri->conn[ke][a];
            for(int d = 0; d < 3; ++d) xvol[a][d] = fe_pri->x[gid][d];
        }
        for(int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            for(int d = 0; d < 3; ++d) xsurf[a][d] = fe_pri->x[gid][d];
        }

        double xc[3] = {0.0,0.0,0.0};
        double xf[3] = {0.0,0.0,0.0};
        for(int a = 0; a < 6; ++a)
            for(int d = 0; d < 3; ++d) xc[d] += xvol[a][d]/6.0;
        for(int a = 0; a < 3; ++a)
            for(int d = 0; d < 3; ++d) xf[d] += xsurf[a][d]/3.0;

        const double outdir[3] = {
            xf[0]-xc[0], xf[1]-xc[1], xf[2]-xc[2]
        };

        double e1[3] = {
            xsurf[1][0]-xsurf[0][0],
            xsurf[1][1]-xsurf[0][1],
            xsurf[1][2]-xsurf[0][2]
        };
        double e2[3] = {
            xsurf[2][0]-xsurf[0][0],
            xsurf[2][1]-xsurf[0][1],
            xsurf[2][2]-xsurf[0][2]
        };
        double cr[3]; pri6_cross3(e1,e2,cr);
        const double face_area = 0.5*pri6_norm3(cr);
        const double volume = pri6_physical_volume(xvol);
        if(face_area <= 0.0 || volume <= 0.0) {
            free(owner); free(lface);
            return -1;
        }

        const double h_n = fmax(volume/face_area, 1.0e-12);
        const double h_e = cbrt(fmax(volume,1.0e-36));

        for(int p = 0; p < np; ++p) {
            double nraw[3];
            pri6_tri3_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jface = pri6_norm3(nraw);
            if(Jface <= 1.0e-300) continue;

            double n[3] = {
                nraw[0]/Jface, nraw[1]/Jface, nraw[2]/Jface
            };
            /* outward normal of the FLUID prism cell */
            if(pri6_dot3(n,outdir) < 0.0) {
                n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
            }

            double Nv_surface[6] = {0,0,0,0,0,0};
            for(int a = 0; a < 3; ++a)
                Nv_surface[s2v[a]] = basis_surf->N[p][a];

            double Nv[6], gradN[6][3], Jinv[3][3], detJ = 0.0;
            if(!pri6_eval_on_tri_face(
                    xvol, lf, Nv_surface, Nv, gradN, Jinv, &detJ)) {
                free(owner); free(lface);
                return -1;
            }

            double u[3] = {0,0,0};
            double uold[3] = {0,0,0};
            double grad_u[3][3] = {{0}};
            double p_q = 0.0;
            double grad_p[3] = {0,0,0};

            for(int a = 0; a < 6; ++a) {
                const int gid = fe_pri->conn[ke][a];
                p_q += Nv[a]*vals->p[gid];
                for(int c = 0; c < 3; ++c) {
                    u[c] += Nv[a]*vals->v[gid][c];
                    uold[c] += Nv[a]*vals->v_old[gid][c];
                    for(int d = 0; d < 3; ++d)
                        grad_u[c][d] += gradN[a][d]*vals->v[gid][c];
                }
                for(int d = 0; d < 3; ++d)
                    grad_p[d] += gradN[a][d]*vals->p[gid];
            }

            const double w =
                basis_surf->integ_weight[p]*area_factor*Jface;

            diag->area += w;
            diag->face_qp_count += 1.0;

            /* ---------------------------------------------------------- */
            /* Pressure force on the BODY: + p n_fluid                    */
            /* This matches the existing TET4/TRI3 diagnostic convention. */
            /* ---------------------------------------------------------- */
            for(int c = 0; c < 3; ++c)
                diag->force_pressure[c] += p_q*n[c]*w;

            /* ---------------------------------------------------------- */
            /* Existing FOM viscous-force convention.                     */
            /*                                                           */
            /* Volume viscous form:                                       */
            /*   mu grad(u):grad(w) + mu div(u) div(w)                    */
            /* boundary traction on fluid:                                */
            /*   t = mu grad(u)n + mu n div(u)                            */
            /* body reaction = -t                                         */
            /* ---------------------------------------------------------- */
            const double divu =
                grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

            double gradun[3] = {0,0,0};
            double gradTun[3] = {0,0,0};
            for(int a = 0; a < 3; ++a) {
                for(int b = 0; b < 3; ++b) {
                    gradun[a]  += grad_u[a][b]*n[b];
                    gradTun[a] += grad_u[b][a]*n[b];
                }
            }

            for(int c = 0; c < 3; ++c) {
                const double t_vol =
                    mu_molecular*(gradun[c] + divu*n[c]);
                diag->force_viscous[c] += -t_vol*w;

                /* Additional physical-Cauchy diagnostic for comparison. */
                const double t_cauchy =
                    mu_molecular*(gradun[c] + gradTun[c]);
                diag->force_cauchy_viscous[c] += -t_cauchy*w;
            }

            /* ---------------------------------------------------------- */
            /* Nitsche penalty reaction, kept SEPARATE from physical force. */
            /* ---------------------------------------------------------- */
            double* grad_uptr[3] = {
                grad_u[0], grad_u[1], grad_u[2]
            };
            double mu_eff = mu_molecular, tau = 0.0, tau_c = 0.0;
            BBFE_vms_mu_eff_tau(
                &mu_eff, &tau, &tau_c,
                Jinv, detJ, h_e,
                u, uold, grad_uptr, grad_p,
                rho, mu_molecular, vals->dt,
                vals->C_vms, vals->vms_cap_coeff);
            (void)tau; (void)tau_c;
            if(!isfinite(mu_eff) || mu_eff <= 0.0) mu_eff = mu_molecular;

            double beta = 0.0;
            pri6_get_penalty(opt, rho, mu_eff, h_n, vals->dt, &beta);

            const double rwall[3] = {
                u[0]-Uw[0], u[1]-Uw[1], u[2]-Uw[2]
            };
            const double rn = pri6_dot3(rwall,n);

            if(opt->wall_mode == PRI6_WALL_NITSCHE_NOSLIP) {
                for(int c = 0; c < 3; ++c)
                    diag->force_nitsche[c] += -beta*rwall[c]*w;
            }
            else if(opt->wall_mode == PRI6_WALL_NITSCHE_FREESLIP ||
                    opt->wall_mode == PRI6_WALL_NITSCHE_SPALDING) {
                for(int c = 0; c < 3; ++c)
                    diag->force_nitsche[c] += -beta*rn*n[c]*w;
            }

            /* ---------------------------------------------------------- */
            /* Spalding tangential wall-law reaction, if selected.         */
            /* ---------------------------------------------------------- */
            if(opt->wall_mode == PRI6_WALL_NITSCHE_SPALDING ||
               opt->wall_mode == PRI6_WALL_STRONG_NORMAL_SPALDING) {

                double Nexchange[6];
                if(opt->exchange_mode == PRI6_WALL_EXCHANGE_OPPOSITE_FACE)
                    pri6_opposite_face_weights(lf,Nv,Nexchange);
                else
                    for(int a = 0; a < 6; ++a) Nexchange[a] = Nv[a];

                double u_exchange[3] = {0,0,0};
                for(int a = 0; a < 6; ++a) {
                    const int gid = fe_pri->conn[ke][a];
                    for(int c = 0; c < 3; ++c)
                        u_exchange[c] += Nexchange[a]*vals->v[gid][c];
                }

                double ut[3] = {0,0,0};
                double beta_wall = 0.0;
                pri6_spalding_coefficients(
                    opt, u_exchange, Uw, n,
                    rho, mu_molecular, h_n,
                    ut, &beta_wall);

                for(int c = 0; c < 3; ++c)
                    diag->force_walllaw[c] += beta_wall*ut[c]*w;
            }
        }

        diag->face_count += 1.0;
    }

    free(owner);
    free(lface);
    return 0;
}

int BBFE_pri6_force_diagnostics_global(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    MONOLIS_COM* mono_com,
    PRI6ForceDiagnostics* diag)
{
    PRI6ForceDiagnostics local;
    const int st = pri6_force_diag_local(
        surf, fe_pri, basis_surf, vals,
        rho, mu_molecular, Uw, opt, &local);

    int failures = (st == 0) ? 0 : 1;
    monolis_allreduce_I(
        1, &failures, MONOLIS_MPI_SUM, mono_com->comm);
    if(failures != 0) {
        pri6_force_diag_zero(diag);
        return -1;
    }

    *diag = local;

    monolis_allreduce_R(1, &diag->area, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &diag->face_count, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &diag->face_qp_count, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(3, diag->force_pressure, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(3, diag->force_viscous, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(3, diag->force_cauchy_viscous, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(3, diag->force_nitsche, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(3, diag->force_walllaw, MONOLIS_MPI_SUM, mono_com->comm);

    return 0;
}

static double pri6_force_dot(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void pri6_force_append(
    const char* filename,
    const char* directory,
    double t,
    double value)
{
    FILE* fp = NULL;
    fp = BBFE_sys_write_add_fopen(fp, filename, directory);
    if(fp != NULL) {
        fprintf(fp, "%.17e %.17e\n", t, value);
        fclose(fp);
    }
}

void BBFE_pri6_reset_karman_force_files(const char* directory)
{
    static const char* files[] = {
        "cylinder_drag_force.txt",
        "cylinder_lift_force.txt",
        "cylinder_drag_coeff.txt",
        "cylinder_lift_coeff.txt",
        "cylinder_drag_coeff_velocity.txt",
        "cylinder_lift_coeff_velocity.txt",
        "cylinder_drag_coeff_pressure.txt",
        "cylinder_lift_coeff_pressure.txt",
        "cylinder_drag_coeff_Nitsche.txt",
        "cylinder_lift_coeff_Nitsche.txt",
        "cylinder_drag_coeff_walllaw.txt",
        "cylinder_lift_coeff_walllaw.txt",
        "cylinder_drag_coeff_cauchy.txt",
        "cylinder_lift_coeff_cauchy.txt"
    };

    for(size_t i = 0; i < sizeof(files)/sizeof(files[0]); ++i) {
        FILE* fp = NULL;
        fp = BBFE_sys_write_fopen(fp, files[i], directory);
        if(fp != NULL) fclose(fp);
    }
}

void BBFE_pri6_output_karman_force_diagnostics(
    const PRI6ForceDiagnostics* d,
    double t,
    double rho,
    double U_ref,
    double A_ref,
    const double eU[3],
    const double eL[3],
    const char* directory)
{
    if(d == NULL) return;

    const double qA = 0.5*rho*U_ref*U_ref*A_ref;
    if(qA <= 0.0) return;

    const double Dp = pri6_force_dot(d->force_pressure,eU);
    const double Lp = pri6_force_dot(d->force_pressure,eL);
    const double Dv = pri6_force_dot(d->force_viscous,eU);
    const double Lv = pri6_force_dot(d->force_viscous,eL);
    const double Dc = pri6_force_dot(d->force_cauchy_viscous,eU);
    const double Lc = pri6_force_dot(d->force_cauchy_viscous,eL);
    const double Dn = pri6_force_dot(d->force_nitsche,eU);
    const double Ln = pri6_force_dot(d->force_nitsche,eL);
    const double Dw = pri6_force_dot(d->force_walllaw,eU);
    const double Lw = pri6_force_dot(d->force_walllaw,eL);

    /* Original resolved-force convention: pressure + volume-form viscous. */
    const double Dtotal = Dp + Dv;
    const double Ltotal = Lp + Lv;

    pri6_force_append("cylinder_drag_force.txt", directory, t, Dtotal);
    pri6_force_append("cylinder_lift_force.txt", directory, t, Ltotal);
    pri6_force_append("cylinder_drag_coeff.txt", directory, t, Dtotal/qA);
    pri6_force_append("cylinder_lift_coeff.txt", directory, t, Ltotal/qA);

    /* Preserve the weak-BC split diagnostic filenames. */
    pri6_force_append("cylinder_drag_coeff_velocity.txt", directory, t, Dv/qA);
    pri6_force_append("cylinder_lift_coeff_velocity.txt", directory, t, Lv/qA);
    pri6_force_append("cylinder_drag_coeff_pressure.txt", directory, t, Dp/qA);
    pri6_force_append("cylinder_lift_coeff_pressure.txt", directory, t, Lp/qA);
    pri6_force_append("cylinder_drag_coeff_Nitsche.txt", directory, t, Dn/qA);
    pri6_force_append("cylinder_lift_coeff_Nitsche.txt", directory, t, Ln/qA);
    pri6_force_append("cylinder_drag_coeff_walllaw.txt", directory, t, Dw/qA);
    pri6_force_append("cylinder_lift_coeff_walllaw.txt", directory, t, Lw/qA);

    /* Extra comparison: pressure + physical symmetric Cauchy viscous stress. */
    pri6_force_append("cylinder_drag_coeff_cauchy.txt", directory, t, (Dp+Dc)/qA);
    pri6_force_append("cylinder_lift_coeff_cauchy.txt", directory, t, (Lp+Lc)/qA);

    fprintf(stdout,
        "[PRI6 force] t=%e area=%e faces=%.0f "
        "Cd=%e Cl=%e [p=(%e,%e) visc=(%e,%e) nit=(%e,%e) wall=(%e,%e)]\n",
        t, d->area, d->face_count,
        Dtotal/qA, Ltotal/qA,
        Dp/qA, Lp/qA,
        Dv/qA, Lv/qA,
        Dn/qA, Ln/qA,
        Dw/qA, Lw/qA);
}
