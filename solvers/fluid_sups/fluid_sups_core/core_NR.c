

#include "core_FOM.h"
#include "core_NR.h"
#include <mkl.h>

static void BBFE_metric_tensor_G(
    double       G[3][3],
    const double J_inv[3][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double gij = 0.0;
            for (int k = 0; k < 3; ++k) {
                gij += J_inv[k][i] * J_inv[k][j];
            }
            G[i][j] = gij;
        }
    }
}

static double BBFE_metric_tensor_vGv(
    const double G[3][3],
    const double v[3])
{
    double vGv = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vGv += v[i] * G[i][j] * v[j];
        }
    }
    return vGv;
}

static double BBFE_metric_tensor_G_colon_G(
    const double G[3][3])
{
    double GG = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            GG += G[i][j] * G[i][j];
        }
    }
    return GG;
}

static double BBFE_metric_tensor_trace(
    const double G[3][3])
{
    return G[0][0] + G[1][1] + G[2][2];
}

static void BBFE_metric_tensor_Gv(
    double       Gv[3],
    const double G[3][3],
    const double v[3])
{
    for (int i = 0; i < 3; ++i) {
        Gv[i] = G[i][0]*v[0] + G[i][1]*v[1] + G[i][2]*v[2];
    }
}

void BBFE_elemmat_fluid_sups_coef_metric_tensor_derivative(
    double       dtau_duj[3],
    double       dtauc_duj[3],
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt,
    const double N_j)
{
    const double EPS = 1.0e-12;

    dtau_duj[0]  = 0.0;
    dtau_duj[1]  = 0.0;
    dtau_duj[2]  = 0.0;
    dtauc_duj[0] = 0.0;
    dtauc_duj[1] = 0.0;
    dtauc_duj[2] = 0.0;

    if (density <= EPS || dt <= EPS || fabs(Jacobian) < EPS) {
        return;
    }

    double G[3][3];
    BBFE_metric_tensor_G(G, J_inv);

    const double vGv = BBFE_metric_tensor_vGv(G, v);
    const double GG  = BBFE_metric_tensor_G_colon_G(G);
    const double trG = BBFE_metric_tensor_trace(G);
    const double nu  = viscosity / density;

    const double dt_term   = 4.0 / (dt * dt);
    const double adv_term  = vGv;
    const double diff_term = 36.0 * nu * nu * GG;

    double denom = dt_term + adv_term + diff_term;
    if (denom < EPS) denom = EPS;

    const double tau = 1.0 / sqrt(denom);

    double Gv[3];
    BBFE_metric_tensor_Gv(Gv, G, v);

    for (int k = 0; k < 3; ++k) {
        /* d tau / d v_k = - tau^3 * (Gv)_k */
        dtau_duj[k] = - N_j * tau * tau * tau * Gv[k];

        /* d tau_c / d v_k = tau * (Gv)_k / trG */
        if (fabs(trG) >= EPS) {
            dtauc_duj[k] = N_j * tau * Gv[k] / trG;
        }
    }
}

void BBFE_elemmat_fluid_mat_rom_linear(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            mat[r][c] = 0.0;
        }
    }

    /* mass */
    const double M = density * N_i * N_j;

    /* viscous */
    const double gNi_dot_gNj =
        grad_N_i[0]*grad_N_j[0] +
        grad_N_i[1]*grad_N_j[1] +
        grad_N_i[2]*grad_N_j[2];

    const double D_11 = dt * viscosity * (gNi_dot_gNj + grad_N_i[0]*grad_N_j[0]);
    const double D_12 = dt * viscosity * (grad_N_i[0]*grad_N_j[1]);
    const double D_13 = dt * viscosity * (grad_N_i[0]*grad_N_j[2]);

    const double D_21 = dt * viscosity * (grad_N_i[1]*grad_N_j[0]);
    const double D_22 = dt * viscosity * (gNi_dot_gNj + grad_N_i[1]*grad_N_j[1]);
    const double D_23 = dt * viscosity * (grad_N_i[1]*grad_N_j[2]);

    const double D_31 = dt * viscosity * (grad_N_i[2]*grad_N_j[0]);
    const double D_32 = dt * viscosity * (grad_N_i[2]*grad_N_j[1]);
    const double D_33 = dt * viscosity * (gNi_dot_gNj + grad_N_i[2]*grad_N_j[2]);

    /* pressure-velocity Galerkin */
    const double G1 = -dt * grad_N_i[0] * N_j;
    const double G2 = -dt * grad_N_i[1] * N_j;
    const double G3 = -dt * grad_N_i[2] * N_j;

    /* continuity Galerkin */
    const double C1 =  dt * N_i * grad_N_j[0];
    const double C2 =  dt * N_i * grad_N_j[1];
    const double C3 =  dt * N_i * grad_N_j[2];

    mat[0][0] += M + D_11;
    mat[0][1] +=     D_12;
    mat[0][2] +=     D_13;
    mat[0][3] +=     G1;

    mat[1][0] +=     D_21;
    mat[1][1] += M + D_22;
    mat[1][2] +=     D_23;
    mat[1][3] +=     G2;

    mat[2][0] +=     D_31;
    mat[2][1] +=     D_32;
    mat[2][2] += M + D_33;
    mat[2][3] +=     G3;

    mat[3][0] += C1;
    mat[3][1] += C2;
    mat[3][2] += C3;
    mat[3][3] += 0.0;
}

void BBFE_elemmat_fluid_mat_rom_nonlinear(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   Jacobian,
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            mat[r][c] = 0.0;
        }
    }

    /* derivatives of tau, tau_c wrt trial dof U_{j,k} */
    double dtau_duj[3], dtauc_duj[3];
    BBFE_elemmat_fluid_sups_coef_metric_tensor_derivative(
        dtau_duj, dtauc_duj,
        J_inv, Jacobian, density, viscosity, v, dt, N_j);

    const double vdotGradNi =
        v[0]*grad_N_i[0] + v[1]*grad_N_i[1] + v[2]*grad_N_i[2];

    const double vdotGradNj =
        v[0]*grad_N_j[0] + v[1]*grad_N_j[1] + v[2]*grad_N_j[2];

    double adv[3];
    for (int d = 0; d < 3; ++d) {
        adv[d] = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
    }

    const double div_v =
        grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

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

    /* -------------------------------------------------
       (1) pure advection (Galerkin)
       ------------------------------------------------- */
    const double A = dt * density * N_i * vdotGradNj;
    mat[0][0] += A;
    mat[1][1] += A;
    mat[2][2] += A;

    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            mat[i][k] += dt * density * N_i * N_j * grad_u[i][k];
        }
    }

    /* -------------------------------------------------
       (2) SUPG terms with tau(v)
       ------------------------------------------------- */

    /* SUPG-time: rho * tau * (v.gradNi) * N_j */
    {
        const double M_s = density * tau * vdotGradNi * N_j;
        mat[0][0] += M_s;
        mat[1][1] += M_s;
        mat[2][2] += M_s;
    }

    /* SUPG-advection: dt rho tau (v.gradNi)(v.gradNj) */
    {
        const double A_s = dt * density * tau * vdotGradNi * vdotGradNj;
        mat[0][0] += A_s;
        mat[1][1] += A_s;
        mat[2][2] += A_s;
    }

    /* SUPG pressure block: dt tau (v.gradNi) gradNj */
    mat[0][3] += dt * tau * vdotGradNi * grad_N_j[0];
    mat[1][3] += dt * tau * vdotGradNi * grad_N_j[1];
    mat[2][3] += dt * tau * vdotGradNi * grad_N_j[2];

    /* derivative wrt delta-v through adv = v.grad_u */
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            mat[i][k] += dt * density * tau * vdotGradNi * N_j * grad_u[i][k];
        }
    }

    /* derivative wrt delta-v through vdotGradNi */
    for (int d = 0; d < 3; ++d) {
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += dt * density * tau * grad_N_i[k] * N_j * adv[d];
        }
    }

    /* derivative of SUPG-pressure residual wrt delta-v */
    for (int d = 0; d < 3; ++d) {
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += dt * tau * grad_N_i[k] * N_j * grad_p[d];
        }
    }

    /* SUPG-time residual wrt delta-v */
    for (int d = 0; d < 3; ++d) {
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += density * tau * grad_N_i[k] * N_j * du_time[d];
        }
    }

    /* -------------------------------------------------
       (3) PSPG terms with tau(v)
       ------------------------------------------------- */

    /* PSPG pressure laplacian */
    mat[3][3] += dt * tau * (
        grad_N_i[0]*grad_N_j[0] +
        grad_N_i[1]*grad_N_j[1] +
        grad_N_i[2]*grad_N_j[2]) / density;

    /* PSPG time */
    for (int k = 0; k < 3; ++k) {
        mat[3][k] += tau * grad_N_i[k] * N_j;
    }

    /* PSPG advection */
    for (int k = 0; k < 3; ++k) {
        mat[3][k] += dt * tau * grad_N_i[k] * vdotGradNj;

        const double dot_gradNi_graduk =
            grad_N_i[0]*grad_u[0][k] +
            grad_N_i[1]*grad_u[1][k] +
            grad_N_i[2]*grad_u[2][k];
        mat[3][k] += dt * tau * N_j * dot_gradNi_graduk;
    }

    /* -------------------------------------------------
       (4) LSIC terms with tau_c(v)
       ------------------------------------------------- */
    {
        const double Kls = dt * density * tau_c;

        mat[0][0] += Kls * grad_N_i[0]*grad_N_j[0];
        mat[0][1] += Kls * grad_N_i[0]*grad_N_j[1];
        mat[0][2] += Kls * grad_N_i[0]*grad_N_j[2];

        mat[1][0] += Kls * grad_N_i[1]*grad_N_j[0];
        mat[1][1] += Kls * grad_N_i[1]*grad_N_j[1];
        mat[1][2] += Kls * grad_N_i[1]*grad_N_j[2];

        mat[2][0] += Kls * grad_N_i[2]*grad_N_j[0];
        mat[2][1] += Kls * grad_N_i[2]*grad_N_j[1];
        mat[2][2] += Kls * grad_N_i[2]*grad_N_j[2];
    }

    /* -------------------------------------------------
       (5) chain-rule terms from dtau
       ------------------------------------------------- */
    for (int k = 0; k < 3; ++k) {
        const double dtau_k = dtau_duj[k];

        if (dtau_k != 0.0) {
            for (int d = 0; d < 3; ++d) {
                const double Rm_tau =
                    dt * density * vdotGradNi * adv[d]
                  + dt * vdotGradNi * grad_p[d]
                  + density * vdotGradNi * du_time[d];

                mat[d][k] += dtau_k * Rm_tau;
            }

            {
                const double Rc_tau =
                    dt * gradN_dot_adv
                  + dt * gradN_dot_gradp / density
                  + gradN_dot_dutime;

                mat[3][k] += dtau_k * Rc_tau;
            }
        }
    }

    /* -------------------------------------------------
       (6) chain-rule terms from dtauc
       ------------------------------------------------- */
    for (int k = 0; k < 3; ++k) {
        const double dtauc_k = dtauc_duj[k];

        if (dtauc_k != 0.0) {
            for (int d = 0; d < 3; ++d) {
                mat[d][k] += dt * density * dtauc_k * div_v * grad_N_i[d];
            }
        }
    }
}

void BBFE_elemmat_fluid_vec_rom_linear(
    double         vec[4],
    const double   N_i,
    const double   grad_N_i[3],
    const double   v[3],
    const double   u_old[3],
    double**       grad_u,
    const double   p_cur,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int d = 0; d < 4; ++d) {
        vec[d] = 0.0;
    }

    const double div_v =
        grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    for (int d = 0; d < 3; ++d) {
        const double gradv_dot_gradNi =
            grad_u[d][0]*grad_N_i[0] +
            grad_u[d][1]*grad_N_i[1] +
            grad_u[d][2]*grad_N_i[2];

        /* viscous */
        vec[d] += dt * viscosity * gradv_dot_gradNi;
        vec[d] += dt * viscosity * div_v * grad_N_i[d];

        /* pressure-velocity Galerkin */
        vec[d] -= dt * grad_N_i[d] * p_cur;

        /* mass */
        vec[d] += density * N_i * (v[d] - u_old[d]);
        //vec[d] += density * N_i * (v[d]);
    }

    /* continuity Galerkin */
    vec[3] += dt * N_i * div_v;
}


void BBFE_elemmat_fluid_vec_rom_linear2(
    double         vec[4],
    const double   N_i,
    const double   grad_N_i[3],
    const double   v[3],
    const double   u_old[3],
    double**       grad_u,
    const double   p_cur,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int d = 0; d < 4; ++d) {
        vec[d] = 0.0;
    }

    const double div_v =
        grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    for (int d = 0; d < 3; ++d) {
        const double gradv_dot_gradNi =
            grad_u[d][0]*grad_N_i[0] +
            grad_u[d][1]*grad_N_i[1] +
            grad_u[d][2]*grad_N_i[2];

        /* viscous */
        //vec[d] += dt * viscosity * gradv_dot_gradNi;
        //vec[d] += dt * viscosity * div_v * grad_N_i[d];

        /* pressure-velocity Galerkin */
        //vec[d] -= dt * grad_N_i[d] * p_cur;

        /* mass */
        //vec[d] += density * N_i * (v[d] - u_old[d]);
        //vec[d] += density * N_i * (v[d]);
        vec[d] += density * N_i * (- u_old[d]);
    }

    /* continuity Galerkin */
    //vec[3] += dt * N_i * div_v;
}


void BBFE_elemmat_fluid_vec_rom_nonlinear(
    double         vec[4],
    const double   N_i,
    const double   grad_N_i[3],
    const double   v[3],
    const double   u_old[3],
    double**       grad_u,
    const double   p_cur,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int d = 0; d < 4; ++d) {
        vec[d] = 0.0;
    }

    const double vdotGradNi =
        v[0]*grad_N_i[0] + v[1]*grad_N_i[1] + v[2]*grad_N_i[2];

    double adv[3];
    for (int d = 0; d < 3; ++d) {
        adv[d] = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
    }

    const double div_v =
        grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    /* -------------------------------------------------
       (1) Galerkin advection
       ------------------------------------------------- */
    for (int d = 0; d < 3; ++d) {
        vec[d] += dt * density * N_i * adv[d];
    }

    /* -------------------------------------------------
       (2) SUPG
       ------------------------------------------------- */
    for (int d = 0; d < 3; ++d) {
        vec[d] += dt * density * tau * vdotGradNi * adv[d];
        vec[d] += dt * tau * vdotGradNi * grad_p[d];
        vec[d] += density * tau * vdotGradNi * (v[d] - u_old[d]);
    }

    /* -------------------------------------------------
       (3) PSPG
       ------------------------------------------------- */
    {
        const double gradN_dot_gradp =
            grad_N_i[0]*grad_p[0] +
            grad_N_i[1]*grad_p[1] +
            grad_N_i[2]*grad_p[2];

        vec[3] += dt * tau * gradN_dot_gradp / density;
    }

    {
        const double gradN_dot_adv =
            grad_N_i[0]*adv[0] +
            grad_N_i[1]*adv[1] +
            grad_N_i[2]*adv[2];

        vec[3] += dt * tau * gradN_dot_adv;
    }

    vec[3] += tau * (
        grad_N_i[0]*(v[0] - u_old[0]) +
        grad_N_i[1]*(v[1] - u_old[1]) +
        grad_N_i[2]*(v[2] - u_old[2]));

    /* -------------------------------------------------
       (4) LSIC
       ------------------------------------------------- */
    for (int d = 0; d < 3; ++d) {
        vec[d] += dt * density * tau_c * div_v * grad_N_i[d];

        //printf("integ_val_vec = %e", vec[d]);
    }
}



void set_element_mat_NR_linear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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

                    BBFE_elemmat_fluid_mat_rom_linear(
                        A, J_inv,
                        basis->N[p][i], basis->N[p][j],
                        fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                        v_ip[p], grad_v_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);
                    /*
                    BBFE_elemmat_fluid_mat_rom_nonlinear(
                        A, J_inv, fe->geo[e][p].Jacobian,
                        basis->N[p][i], basis->N[p][j],
                        fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                        v_ip[p], grad_v_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);
                    */

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


void set_element_vec_NR_linear_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
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

                BBFE_elemmat_fluid_vec_rom_linear(
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
/*
                BBFE_elemmat_fluid_vec_rom_nonlinear(
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
*/          
            
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                    vec[d] = 0.0;
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }

    
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
              
                BBFE_elemmat_fluid_vec_rom_nonlinear(
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
                    vec[d] = 0.0;
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

void set_element_vec_NR_linear_nonlinear2(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
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

                BBFE_elemmat_fluid_vec_rom_linear2(
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
/*
                BBFE_elemmat_fluid_vec_rom_nonlinear(
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
*/          
            
                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                    vec[d] = 0.0;
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }

    
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
              
                BBFE_elemmat_fluid_vec_rom_nonlinear(
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
                    vec[d] = 0.0;
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


void set_element_vec_NR_linear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
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

                BBFE_elemmat_fluid_vec_rom_linear(
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
                    vec[d] = 0.0;
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }
}


void set_element_mat_NR_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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

                    BBFE_elemmat_fluid_mat_rom_nonlinear(
                        A, J_inv, fe->geo[e][p].Jacobian,
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


void set_element_vec_NR_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
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
                
                BBFE_elemmat_fluid_vec_rom_nonlinear(
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
                    vec[d] = 0.0;
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


void HROM_ddecm_set_residuals_NR(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t)
{
    int ns = index_snap;
    double** local_matrix;  double* local_vec;
    local_matrix   = BB_std_calloc_2d_double(local_matrix, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    local_vec   = BB_std_calloc_1d_double(local_vec, hlpod_vals->n_neib_vec);


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

    for(int n=0; n < num_subdomains; n++) {
        for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];
                
            //for (int e = 0; e < fe->total_num_elems; ++e) {

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

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
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

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                    int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                    int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                    int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                    // i–j ブロックの積分・加算
                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                                if( bc->D_bc_exists[index_j*4+b]) {
                                    for(int k1 = IS; k1 < IE; k1++){
                                        double val = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val * bc->imposed_D_val[index_j*4+b];

                                        hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] += val;
                                        hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] += val;

                                        //hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] -= val;
                                        //hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] -= val;
                                    }
                                }
                                else{
                                    for(int k1 = IS; k1 < IE; k1++) {
                                        double A = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val;
                                        local_vec[k1] = 0.0;

                                        int index1 = 0;
                                        int index2 = 0;

                                        for(int ki = 0; ki < num_neib; ki++) {
                                            for(int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; kj++) {
                                                double B = hlpod_mat->neib_vec[index_j*4+b][index2];
                                                double C = hlpod_mat->mode_coef_1stdd[index1 + kj];
                                                //double C = hlpod_mat->mode_coef[index1 + kj];

                                                local_vec[k1] += A * B * C;
                                                index2++;
                                            }
                                            index1 += hlpod_mat->max_num_neib_modes[ki];
                                        }
                                    }

                                    for(int k1 = IS; k1 < IE; k1++) {
                                        int index = ns*hlpod_vals->n_neib_vec + k1;

                                        hlpod_ddhr->matrix[index][m][n] += local_vec[k1];
                                        hlpod_ddhr->RH[index][n] += local_vec[k1];

                                        if(monolis_mpi_get_global_my_rank()==0){
//                                        	printf("local_vec[k1] = %e\n", local_vec[k1]);
                                        }

                                    }
                                }



                        }
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


#include <mkl.h>
#include <string.h>

static inline int max2_int(int a, int b){ return (a > b) ? a : b; }

void HROM_ddecm_set_residuals_NR_blas(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t)
{
    (void)num_snapshot;
    (void)num_neib;
    (void)dt;
    (void)t;

    const int ns = index_snap;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int nrv = hlpod_vals->n_neib_vec;

    /* ---------------------------------
     * compact coefficient:
     * original:
     *   index2: compact j-basis index
     *   index1+kj: padded coefficient index
     * --------------------------------- */
    double* coef_compact = (double*)mkl_malloc(sizeof(double) * (size_t)nrv, 64);
    {
        int index1 = 0;
        int index2 = 0;
        for (int ki = 0; ki < num_neib; ++ki) {
            for (int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; ++kj) {
                coef_compact[index2++] = hlpod_mat->mode_coef_1stdd[index1 + kj];
            }
            index1 += hlpod_mat->max_num_neib_modes[ki];
        }
    }

    /* ---- geometry / field work ---- */
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double**  local_v     = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip        = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip   = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double**  local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*   local_p     = BB_std_calloc_1d_double(NULL, nl);
    double*   p_ip        = BB_std_calloc_1d_double(NULL, np);
    double**  grad_p_ip   = BB_std_calloc_2d_double(NULL, np, 3);

    double*   wJ          = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau         = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau_c       = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);

    /* ---- element-local metadata ---- */
    int* conn_e = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IS_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IE_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* msz_i  = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);

    /* j-side packed info */
    double* s_j = (double*)mkl_malloc(sizeof(double) * (size_t)nl * 4, 64);
    /* s_j[j*4+b] :
     *   free column -> dot(phi_j(:,b), coef_compact)
     *   D column    -> imposed_D_val
     */

    /* i-side basis pack: column-major (msz_max x 4) */
    int msz_max_global = 0;
    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];
                msz_max_global = max2_int(msz_max_global, IE - IS);
            }
        }
    }

    double* Phi_i_buf = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global * 4, 64);
    double* y_buf     = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global, 64);

    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];

            /* -----------------------------
             * element precompute
             * ----------------------------- */
            BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

            BBFE_elemmat_set_local_array_vector(local_v,     fe, vals->v,     e, 3);
            BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);
            BBFE_elemmat_set_local_array_scalar(local_p,     fe, vals->p,     e);

            for (int p = 0; p < np; ++p) {
                double J_inv[3][3];

                BBFE_std_mapping_vector3d(v_ip[p],     nl, local_v,     basis->N[p]);
                BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
                BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
                BBFE_std_mapping_scalar_grad(grad_p_ip[p],  nl, local_p, fe->geo[e][p].grad_N);
                p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                tau[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                tau_c[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                wJ[p] = basis->integ_weight[p] * Jacobian_ip[p];
            }

            /* -----------------------------
             * element node metadata
             * ----------------------------- */
            int msz_max_elem = 0;
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];

                conn_e[i] = gi;
                IS_i[i]   = IS;
                IE_i[i]   = IE;
                msz_i[i]  = IE - IS;

                msz_max_elem = max2_int(msz_max_elem, msz_i[i]);
            }

            if (msz_max_elem <= 0) continue;

            /* -----------------------------
             * j-side scalar summary:
             *   free  -> ddot(phi_j, coef_compact)
             *   D-bc  -> imposed_D_val
             * ----------------------------- */
            for (int j = 0; j < nl; ++j) {
                const int gj = conn_e[j];
                for (int b = 0; b < 4; ++b) {
                    if (bc->D_bc_exists[gj*4 + b]) {
                        s_j[j*4 + b] = bc->imposed_D_val[gj*4 + b];
                    }
                    else {
                        s_j[j*4 + b] = cblas_ddot(
                            nrv,
                            hlpod_mat->neib_vec[gj*4 + b], 1,
                            coef_compact, 1);
                    }
                }
            }

            /* -----------------------------
             * (i,j) loop
             * ----------------------------- */
            for (int i = 0; i < nl; ++i) {
                const int gi  = conn_e[i];
                const int IS  = IS_i[i];
                const int IE  = IE_i[i];
                const int msz = msz_i[i];

                if (msz <= 0) continue;

                /* Phi_i_buf: (msz x 4), column-major */
                for (int a = 0; a < 4; ++a) {
                    memcpy(Phi_i_buf + (size_t)a * (size_t)msz_max_global,
                           &hlpod_mat->neib_vec[gi*4 + a][IS],
                           sizeof(double) * (size_t)msz);
                }

                for (int j = 0; j < nl; ++j) {
                    const int gj = conn_e[j];

                    /* acc[a + 4*b] : column-major 4x4 */
                    double acc[16];
                    for (int q = 0; q < 16; ++q) acc[q] = 0.0;

                    for (int p = 0; p < np; ++p) {
                        double A[4][4];
                        double du_time[3] = {
                            v_ip[p][0] - v_ip_old[p][0],
                            v_ip[p][1] - v_ip_old[p][1],
                            v_ip[p][2] - v_ip_old[p][2]
                        };
                        double J_inv[3][3];

                        BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity,
                            tau[p], tau_c[p], vals->dt, du_time);

                        const double wp = wJ[p];
                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                acc[a + 4*b] += A[a][b] * wp;
                            }
                        }
                    }

                    /* g[a] = sum_b acc[a,b] * s_jb */
                    double g[4];
                    for (int a = 0; a < 4; ++a) {
                        g[a] =
                            acc[a + 4*0] * s_j[j*4 + 0] +
                            acc[a + 4*1] * s_j[j*4 + 1] +
                            acc[a + 4*2] * s_j[j*4 + 2] +
                            acc[a + 4*3] * s_j[j*4 + 3];
                    }

                    /* y = Phi_i(msz x 4) * g(4) */
                    cblas_dgemv(CblasColMajor, CblasNoTrans,
                                msz, 4,
                                1.0,
                                Phi_i_buf, msz_max_global,
                                g, 1,
                                0.0,
                                y_buf, 1);

                    /* scatter add */
                    for (int rr = 0; rr < msz; ++rr) {
                        const int index = ns * nrv + (IS + rr);
                        const double val = y_buf[rr];

                        hlpod_ddhr->matrix[index][m][n] += val;
                        hlpod_ddhr->RH[index][n]        += val;
                    }
                }
            }
        }
    }

    /* ---- free ---- */
    mkl_free(y_buf);
    mkl_free(Phi_i_buf);

    mkl_free(s_j);

    mkl_free(conn_e);
    mkl_free(IS_i);
    mkl_free(IE_i);
    mkl_free(msz_i);

    mkl_free(wJ);
    mkl_free(tau);
    mkl_free(tau_c);

    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    mkl_free(coef_compact);
}

void HROM_set_element_mat_NR(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*    hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

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


    int rank = monolis_mpi_get_global_my_rank();

    
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems
            : hlpod_ddhr->ovl_num_selected_elems_D_bc;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems
            : hlpod_ddhr->ovl_id_selected_elems_D_bc;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight
            : hlpod_ddhr->ovl_elem_weight_D_bc;

        for (int m = 0; m < num; m++) {
            int e = ids[m];
 
    //for (int e = 0; e < fe->total_num_elems; ++e) {
        
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

                        /*
                        BBFE_elemmat_fluid_mat_rom_linear(
                            A, J_inv,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);
                        */

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                val_ip[a][b][p] = A[a][b];

                                if (!isfinite(A[a][b])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, a, A[a][b]);
                                }

                                A[a][b] = 0.0;
                            }
                        }
                    }

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        if( bc->D_bc_exists[index_j*4+b]) {
                        }
                        else{
                            int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                            int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                            int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                            int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                            int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                            int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                            int subdomain_id = hlpod_mat->subdomain_id_in_nodes_2nddd[index_j]; 

                            if(subdomain_id_j < num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS; k2 < JE; k2++){
                                        double val = hlpod_mat->pod_basis_hr[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j*4+b][k2];
                                        //printf("val = %e", val);
                                        hlpod_ddhr->reduced_mat[k1][k2] += weight[m] *  val;
                                        //hlpod_ddhr->reduced_mat[k1][k2] += val;
                                    }
                                }
                            }

                            if(subdomain_id_j >= num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
                                        double val = hlpod_mat->pod_basis_hr[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j*4+b][k2] ;
                                        hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += weight[m] *  val;
                                        //hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += val;
                                        // /printf("val = %e", val);
                                    }
                                }
                            }
                        }
                    }
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

void HROM_set_element_mat_NR2(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*    hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            //hlpod_mat->VTKV[i][j] = 0.0;
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

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


    int rank = monolis_mpi_get_global_my_rank();

    
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems
            : hlpod_ddhr->ovl_num_selected_elems_D_bc;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems
            : hlpod_ddhr->ovl_id_selected_elems_D_bc;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight
            : hlpod_ddhr->ovl_elem_weight_D_bc;

        for (int m = 0; m < num; m++) {
            int e = ids[m];
    
    //for (int e = 0; e < fe->total_num_elems; ++e) {
        
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

                        BBFE_elemmat_fluid_mat_rom_linear(
                            A, J_inv,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                val_ip[a][b][p] = A[a][b];

                                if (!isfinite(A[a][b])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, a, A[a][b]);
                                }

                                A[a][b] = 0.0;
                            }
                        }
                    }

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);
                            
                            for(int k1 = 0; k1 < hlpod_vals->n_neib_vec; k1++){
                                for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
                                    double val = hlpod_mat->pod_basis_hr[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j*4+b][k2];

                                    hlpod_ddhr->reduced_mat[k1][k2] += val * weight[m];
                                    //hlpod_mat->VTKV[k1][k2] += val;

                                }
                            }
                        }
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



void HROM_set_element_vec_NR(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
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

    int rank = monolis_mpi_get_global_my_rank();
    
    
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems
            : hlpod_ddhr->ovl_num_selected_elems_D_bc;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems
            : hlpod_ddhr->ovl_id_selected_elems_D_bc;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight
            : hlpod_ddhr->ovl_elem_weight_D_bc;

        for (int m = 0; m < num; m++) {
            int e = ids[m];
            
    
        //for(int e = 0; e < fe->total_num_elems; e++) {

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
            
            for (int i = 0; i < nl; ++i) {

                int index = fe->conn[e][i];
                int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

                if(subdomain_id < num_subdomains){

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
                    
                        BBFE_elemmat_fluid_vec_rom_nonlinear(
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

                                if (!isfinite(vec[d])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, d, vec[d]);
                                }

                            vec[d] = 0.0;
                        }
                    }

                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                    for (int d = 0; d < 4; ++d) {
                        integ_val_vec[d] = BBFE_std_integ_calc(
                            np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);
                        
                            if (!isfinite(integ_val_vec[d])) {
                                printf("rank = %d, bad integ_val_vec: e=%d i=%d d=%d val=%e\n", rank, e, i, d, integ_val_vec[d]);
                            }      
                        
                        for(int k = IS; k < IE; k++){
                            double val = integ_val_vec[d] * hlpod_mat->pod_modes[index*4+d][k];

                            monolis->mat.R.B[k] -= weight[m] * val;
                            //monolis->mat.R.B[k] -= val;
                        
                            double mode = hlpod_mat->pod_modes[index*4+d][k];
                            double val2 = integ_val_vec[d] * mode;

                            if (!isfinite(mode) || !isfinite(val)) {
                                printf("rank = %d, bad reduction: e=%d i=%d index=%d d=%d k=%d integ=%e mode=%e val=%e\n",
                                    rank, e, i, index, d, k, integ_val_vec[d], mode, val);
                            }

                            if(val==0.0){
                            }
                            else{
                            //hlpod_ddhr->reduced_RH[k] -= val;
                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                            }

                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                        }
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


void BBFE_elemmat_fluid_mat_mass(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            mat[r][c] = 0.0;
        }
    }

    /* mass */
    const double M = density * N_i * N_j;

    mat[0][0] += M;
    mat[1][1] += M;
    mat[2][2] += M;
}

void set_element_mat_NR_mass(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
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

                    BBFE_elemmat_fluid_mat_mass(
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


void ROM_std_hlpod_reduced_rhs_to_monollis_linear2(
    MONOLIS*		monolis,
    MONOLIS*		monolis_mass,
    MONOLIS_COM*    monolis_com,
    HLPOD_MAT*      hlpod_mat,
    double*         mode_coeff,
    double*         mode_coeff_old,
    const int       num_2nd_subdomains)
{
    int index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf_tmp[index + i]  = 0.0;
        }
        index += hlpod_mat->num_modes_internal[k];
    }

    double* vec_in = BB_std_calloc_1d_double(vec_in, index);
    double* vec_out = BB_std_calloc_1d_double(vec_out, index);

    monolis_matvec_product_R(monolis_mass, monolis_com, mode_coeff_old, hlpod_mat->VTf_tmp);

    index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            vec_in[index + i] = hlpod_mat->VTf_tmp[index + i];
            hlpod_mat->VTf_tmp[index + i] = 0.0;
        }
        index += hlpod_mat->num_modes_internal[k];
    }

    monolis_matvec_product_R(monolis, monolis_com, mode_coeff, hlpod_mat->VTf_tmp);
    
    index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            vec_out[index + i] = hlpod_mat->VTf_tmp[index + i];
            hlpod_mat->VTf_tmp[index + i] = 0.0;
        }
        index += hlpod_mat->num_modes_internal[k];
    }

    index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf_tmp[index + i]  = - ( vec_out[index + i] - vec_in[index + i] );
            //printf("%e ", hlpod_mat->VTf_tmp[index + i]);
        }
        index += hlpod_mat->num_modes_internal[k];
    }
    
}


void BBFE_elemmat_fluid_mat_rom_nonlinear_PSPG(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   Jacobian,
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            mat[r][c] = 0.0;
        }
    }

    /* derivatives of tau, tau_c wrt trial dof U_{j,k} */
    double dtau_duj[3], dtauc_duj[3];
    BBFE_elemmat_fluid_sups_coef_metric_tensor_derivative(
        dtau_duj, dtauc_duj,
        J_inv, Jacobian, density, viscosity, v, dt, N_j);

    double adv[3];
    for (int d = 0; d < 3; ++d) {
        adv[d] = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
    }

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

    const double Rc_tau =
        dt * gradN_dot_adv
      + dt * gradN_dot_gradp / density
      + gradN_dot_dutime;

    for (int k = 0; k < 3; ++k) {
        const double dtau_k = dtau_duj[k];
        if (dtau_k != 0.0) {
            {
                const double Rc_tau =
                    dt * gradN_dot_adv
                  + dt * gradN_dot_gradp / density
                  + gradN_dot_dutime;

                mat[3][k] += dtau_k * Rc_tau;
            }
        }
    }
}

/*
void HROM_ddecm_set_residuals_NR_PSPG(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t)
{
    int ns = index_snap;
    double** local_matrix;  double* local_vec;
    local_matrix   = BB_std_calloc_2d_double(local_matrix, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    local_vec   = BB_std_calloc_1d_double(local_vec, hlpod_vals->n_neib_vec);


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

    for(int n=0; n < num_subdomains; n++) {
        for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];
                
            //for (int e = 0; e < fe->total_num_elems; ++e) {

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

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
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

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                    int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                    int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                    int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                    // i–j ブロックの積分・加算
                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                                if( bc->D_bc_exists[index_j*4+b]) {
                                    for(int k1 = IS; k1 < IE; k1++){
                                        double val = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val * bc->imposed_D_val[index_j*4+b];

                                        hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] += val;
                                        hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] += val;

                                        //hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] -= val;
                                        //hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] -= val;
                                    }
                                }
                                else{
                                    for(int k1 = IS; k1 < IE; k1++) {
                                        double A = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val;
                                        local_vec[k1] = 0.0;

                                        int index1 = 0;
                                        int index2 = 0;

                                        for(int ki = 0; ki < num_neib; ki++) {
                                            for(int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; kj++) {
                                                double B = hlpod_mat->neib_vec[index_j*4+b][index2];
                                                double C = hlpod_mat->mode_coef_1stdd[index1 + kj];
                                                //double C = hlpod_mat->mode_coef[index1 + kj];

                                                local_vec[k1] += A * B * C;
                                                index2++;
                                            }
                                            index1 += hlpod_mat->max_num_neib_modes[ki];
                                        }
                                    }

                                    for(int k1 = IS; k1 < IE; k1++) {
                                        int index = 2*ns*hlpod_vals->n_neib_vec + k1;

                                        hlpod_ddhr->matrix[index][m][n] += local_vec[k1];
                                        hlpod_ddhr->RH[index][n] += local_vec[k1];

                                        if(monolis_mpi_get_global_my_rank()==0){
                                        	printf("local_vec[k1] = %e\n", local_vec[k1]);
                                        }

                                    }
                                }



                        }
                    }

                }
            }
        }
    }


    for(int n=0; n < num_subdomains; n++) {
        for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];
                
            //for (int e = 0; e < fe->total_num_elems; ++e) {

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

                        BBFE_elemmat_fluid_mat_rom_nonlinear_PSPG(
                            A, J_inv, fe->geo[e][p].Jacobian,
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

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                    int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                    int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                    int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                    // i–j ブロックの積分・加算
                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                                for(int k1 = IS; k1 < IE; k1++) {
                                    double A = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val;
                                    local_vec[k1] = 0.0;

                                    int index1 = 0;
                                    int index2 = 0;

                                    for(int ki = 0; ki < num_neib; ki++) {
                                        for(int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; kj++) {
                                            double B = hlpod_mat->neib_vec[index_j*4+b][index2];
                                            double C = hlpod_mat->mode_coef_1stdd[index1 + kj];
                                            //double C = hlpod_mat->mode_coef[index1 + kj];

                                            local_vec[k1] += A * B * C;
                                            index2++;
                                        }
                                        index1 += hlpod_mat->max_num_neib_modes[ki];
                                    }
                                }

                                for(int k1 = IS; k1 < IE; k1++) {
                                    int index = 2*ns*hlpod_vals->n_neib_vec + hlpod_vals->n_neib_vec + k1;

                                    hlpod_ddhr->matrix[index][m][n] += local_vec[k1];
                                    hlpod_ddhr->RH[index][n] += local_vec[k1];
                            }



                        }
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
*/


//static inline int max2_int(int a, int b){ return (a > b) ? a : b; }


void HROM_ddecm_set_residuals_NR_blas2(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t)
{
    (void)num_snapshot;
    (void)num_neib;
    (void)dt;
    (void)t;

    const int ns = index_snap;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int nrv = hlpod_vals->n_neib_vec;

    /* ---------------------------------
     * compact coefficient:
     * original:
     *   index2: compact j-basis index
     *   index1+kj: padded coefficient index
     * --------------------------------- */
    double* coef_compact = (double*)mkl_malloc(sizeof(double) * (size_t)nrv, 64);
    {
        int index1 = 0;
        int index2 = 0;
        for (int ki = 0; ki < num_neib; ++ki) {
            for (int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; ++kj) {
                coef_compact[index2++] = hlpod_mat->mode_coef_1stdd[index1 + kj];
            }
            index1 += hlpod_mat->max_num_neib_modes[ki];
        }
    }

    /* ---- geometry / field work ---- */
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double**  local_v     = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip        = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip   = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double**  local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*   local_p     = BB_std_calloc_1d_double(NULL, nl);
    double*   p_ip        = BB_std_calloc_1d_double(NULL, np);
    double**  grad_p_ip   = BB_std_calloc_2d_double(NULL, np, 3);

    double*   wJ          = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau         = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau_c       = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);

    /* ---- element-local metadata ---- */
    int* conn_e = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IS_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IE_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* msz_i  = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);

    /* j-side packed info */
    double* s_j = (double*)mkl_malloc(sizeof(double) * (size_t)nl * 4, 64);
    /* s_j[j*4+b] :
     *   free column -> dot(phi_j(:,b), coef_compact)
     *   D column    -> imposed_D_val
     */

    /* i-side basis pack: column-major (msz_max x 4) */
    int msz_max_global = 0;
    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];
                msz_max_global = max2_int(msz_max_global, IE - IS);
            }
        }
    }

    double* Phi_i_buf = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global * 4, 64);
    double* y_buf     = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global, 64);

    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];

            /* -----------------------------
             * element precompute
             * ----------------------------- */
            BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

            BBFE_elemmat_set_local_array_vector(local_v,     fe, vals->v,     e, 3);
            BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);
            BBFE_elemmat_set_local_array_scalar(local_p,     fe, vals->p,     e);

            for (int p = 0; p < np; ++p) {
                double J_inv[3][3];

                BBFE_std_mapping_vector3d(v_ip[p],     nl, local_v,     basis->N[p]);
                BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
                BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
                BBFE_std_mapping_scalar_grad(grad_p_ip[p],  nl, local_p, fe->geo[e][p].grad_N);
                p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                tau[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                tau_c[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                wJ[p] = basis->integ_weight[p] * Jacobian_ip[p];
            }

            /* -----------------------------
             * element node metadata
             * ----------------------------- */
            int msz_max_elem = 0;
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];

                conn_e[i] = gi;
                IS_i[i]   = IS;
                IE_i[i]   = IE;
                msz_i[i]  = IE - IS;

                msz_max_elem = max2_int(msz_max_elem, msz_i[i]);
            }

            if (msz_max_elem <= 0) continue;

            /* -----------------------------
             * j-side scalar summary:
             *   free  -> ddot(phi_j, coef_compact)
             *   D-bc  -> imposed_D_val
             * ----------------------------- */
            for (int j = 0; j < nl; ++j) {
                const int gj = conn_e[j];
                for (int b = 0; b < 4; ++b) {
                    if (bc->D_bc_exists[gj*4 + b]) {
                        s_j[j*4 + b] = bc->imposed_D_val[gj*4 + b];
                    }
                    else {
                        s_j[j*4 + b] = cblas_ddot(
                            nrv,
                            hlpod_mat->neib_vec[gj*4 + b], 1,
                            coef_compact, 1);
                    }
                }
            }

            /* -----------------------------
             * (i,j) loop
             * ----------------------------- */
            for (int i = 0; i < nl; ++i) {
                const int gi  = conn_e[i];
                const int IS  = IS_i[i];
                const int IE  = IE_i[i];
                const int msz = msz_i[i];

                if (msz <= 0) continue;

                /* Phi_i_buf: (msz x 4), column-major */
                for (int a = 0; a < 4; ++a) {
                    memcpy(Phi_i_buf + (size_t)a * (size_t)msz_max_global,
                           &hlpod_mat->neib_vec[gi*4 + a][IS],
                           sizeof(double) * (size_t)msz);
                }

                for (int j = 0; j < nl; ++j) {
                    const int gj = conn_e[j];

                    /* acc[a + 4*b] : column-major 4x4 */
                    double acc[16];
                    for (int q = 0; q < 16; ++q) acc[q] = 0.0;

                    for (int p = 0; p < np; ++p) {
                        double A[4][4];
                        double du_time[3] = {
                            v_ip[p][0] - v_ip_old[p][0],
                            v_ip[p][1] - v_ip_old[p][1],
                            v_ip[p][2] - v_ip_old[p][2]
                        };
                        double J_inv[3][3];

                        BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity,
                            tau[p], tau_c[p], vals->dt, du_time);

                        const double wp = wJ[p];
                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                acc[a + 4*b] += A[a][b] * wp;
                            }
                        }
                    }

                    /* g[a] = sum_b acc[a,b] * s_jb */
                    double g[4];
                    for (int a = 0; a < 4; ++a) {
                        g[a] =
                            acc[a + 4*0] * s_j[j*4 + 0] +
                            acc[a + 4*1] * s_j[j*4 + 1] +
                            acc[a + 4*2] * s_j[j*4 + 2] +
                            acc[a + 4*3] * s_j[j*4 + 3];
                    }

                    /* y = Phi_i(msz x 4) * g(4) */
                    cblas_dgemv(CblasColMajor, CblasNoTrans,
                                msz, 4,
                                1.0,
                                Phi_i_buf, msz_max_global,
                                g, 1,
                                0.0,
                                y_buf, 1);

                    /* scatter add */
                    for (int rr = 0; rr < msz; ++rr) {
                        const int index = 2 * ns * nrv + (IS + rr);
                        const double val = y_buf[rr];

                        hlpod_ddhr->matrix[index][m][n] += val;
                        hlpod_ddhr->RH[index][n]        += val;
                    }
                }
            }
        }
    }

    /* ---- free ---- */
    mkl_free(y_buf);
    mkl_free(Phi_i_buf);

    mkl_free(s_j);

    mkl_free(conn_e);
    mkl_free(IS_i);
    mkl_free(IE_i);
    mkl_free(msz_i);

    mkl_free(wJ);
    mkl_free(tau);
    mkl_free(tau_c);

    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    mkl_free(coef_compact);
}

void HROM_ddecm_set_residuals_NR_PSPG(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t)
{
    (void)bc;
    (void)num_snapshot;
    (void)dt;
    (void)t;

    const int ns  = index_snap;
    const int nl  = fe->local_num_nodes;
    const int np  = basis->num_integ_points;
    const int nrv = hlpod_vals->n_neib_vec;

    /* ---------------------------------
     * compact coefficient:
     * original:
     *   index2: compact j-basis index
     *   index1+kj: padded coefficient index
     * --------------------------------- */
    double* coef_compact = (double*)mkl_malloc(sizeof(double) * (size_t)nrv, 64);
    {
        int index1 = 0;
        int index2 = 0;
        for (int ki = 0; ki < num_neib; ++ki) {
            for (int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; ++kj) {
                coef_compact[index2++] = hlpod_mat->mode_coef_1stdd[index1 + kj];
            }
            index1 += hlpod_mat->max_num_neib_modes[ki];
        }
    }

    /* ---- geometry / field work ---- */
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double**  local_v     = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip        = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip   = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double**  local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*   local_p     = BB_std_calloc_1d_double(NULL, nl);
    double*   p_ip        = BB_std_calloc_1d_double(NULL, np);
    double**  grad_p_ip   = BB_std_calloc_2d_double(NULL, np, 3);

    double*   wJ          = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau         = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau_c       = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);

    /* ---- element-local metadata ---- */
    int* conn_e = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IS_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IE_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* msz_i  = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);

    /* j-side packed scalar summary */
    double* s_j = (double*)mkl_malloc(sizeof(double) * (size_t)nl * 4, 64);
    /* s_j[j*4+b] = ddot(phi_j(:,b), coef_compact) */

    /* i-side basis pack: column-major (msz_max x 4) */
    int msz_max_global = 0;
    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];
                msz_max_global = max2_int(msz_max_global, IE - IS);
            }
        }
    }

    double* Phi_i_buf = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global * 4, 64);
    double* y_buf     = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global, 64);

    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];

            /* -----------------------------
             * element precompute
             * ----------------------------- */
            BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

            BBFE_elemmat_set_local_array_vector(local_v,     fe, vals->v,     e, 3);
            BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);
            BBFE_elemmat_set_local_array_scalar(local_p,     fe, vals->p,     e);

            for (int p = 0; p < np; ++p) {
                double J_inv[3][3];

                BBFE_std_mapping_vector3d(v_ip[p],     nl, local_v,     basis->N[p]);
                BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
                BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
                BBFE_std_mapping_scalar_grad(grad_p_ip[p],  nl, local_p, fe->geo[e][p].grad_N);
                p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                tau[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                tau_c[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                wJ[p] = basis->integ_weight[p] * Jacobian_ip[p];
            }

            /* -----------------------------
             * element node metadata
             * ----------------------------- */
            int msz_max_elem = 0;
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];

                conn_e[i] = gi;
                IS_i[i]   = IS;
                IE_i[i]   = IE;
                msz_i[i]  = IE - IS;

                msz_max_elem = max2_int(msz_max_elem, msz_i[i]);
            }

            if (msz_max_elem <= 0) continue;

            /* -----------------------------
             * j-side scalar summary:
             * s_j[j*4+b] = ddot(phi_j(:,b), coef_compact)
             * PSPG側は元コードどおり Dirichlet 分岐なし
             * ----------------------------- */
            for (int j = 0; j < nl; ++j) {
                const int gj = conn_e[j];
                for (int b = 0; b < 4; ++b) {
                    s_j[j*4 + b] = cblas_ddot(
                        nrv,
                        hlpod_mat->neib_vec[gj*4 + b], 1,
                        coef_compact, 1);
                }
            }

            /* -----------------------------
             * (i,j) loop
             * ----------------------------- */
            for (int i = 0; i < nl; ++i) {
                const int gi  = conn_e[i];
                const int IS  = IS_i[i];
                const int IE  = IE_i[i];
                const int msz = msz_i[i];

                if (msz <= 0) continue;

                /* Phi_i_buf: (msz x 4), column-major */
                for (int a = 0; a < 4; ++a) {
                    memcpy(Phi_i_buf + (size_t)a * (size_t)msz_max_global,
                           &hlpod_mat->neib_vec[gi*4 + a][IS],
                           sizeof(double) * (size_t)msz);
                }

                for (int j = 0; j < nl; ++j) {
                    /* acc[a + 4*b] : column-major 4x4 */
                    double acc[16];
                    for (int q = 0; q < 16; ++q) acc[q] = 0.0;

                    for (int p = 0; p < np; ++p) {
                        double A[4][4];
                        double du_time[3] = {
                            v_ip[p][0] - v_ip_old[p][0],
                            v_ip[p][1] - v_ip_old[p][1],
                            v_ip[p][2] - v_ip_old[p][2]
                        };
                        double J_inv[3][3];

                        BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                        BBFE_elemmat_fluid_mat_rom_nonlinear_PSPG(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity,
                            tau[p], tau_c[p], vals->dt, du_time);

                        const double wp = wJ[p];
                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                acc[a + 4*b] += A[a][b] * wp;
                            }
                        }
                    }

                    /* g[a] = sum_b acc[a,b] * s_jb */
                    double g[4];
                    for (int a = 0; a < 4; ++a) {
                        g[a] =
                            acc[a + 4*0] * s_j[j*4 + 0] +
                            acc[a + 4*1] * s_j[j*4 + 1] +
                            acc[a + 4*2] * s_j[j*4 + 2] +
                            acc[a + 4*3] * s_j[j*4 + 3];
                    }

                    /* y = Phi_i(msz x 4) * g(4) */
                    cblas_dgemv(CblasColMajor, CblasNoTrans,
                                msz, 4,
                                1.0,
                                Phi_i_buf, msz_max_global,
                                g, 1,
                                0.0,
                                y_buf, 1);

                    /* scatter add */
                    for (int rr = 0; rr < msz; ++rr) {
                        const int index = 2*ns*nrv + nrv + (IS + rr);
                        const double val = y_buf[rr];

                        hlpod_ddhr->matrix[index][m][n] += val;
                        hlpod_ddhr->RH[index][n]        += val;
                    }
                }
            }
        }
    }

    /* ---- free ---- */
    mkl_free(y_buf);
    mkl_free(Phi_i_buf);

    mkl_free(s_j);

    mkl_free(conn_e);
    mkl_free(IS_i);
    mkl_free(IE_i);
    mkl_free(msz_i);

    mkl_free(wJ);
    mkl_free(tau);
    mkl_free(tau_c);

    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    mkl_free(coef_compact);
}


void HROM_ddecm_set_residuals_NR_vec(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t)
{
    printf("\n\nindex_snap = %d, num_modes = %d, num_subdomains = %d\n\n", index_snap, hlpod_vals->n_neib_vec, num_subdomains);

    int ns = index_snap;

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

    //for (int e = 0; e < fe->total_num_elems; ++e) {
	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];

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

                    BBFE_elemmat_fluid_vec_rom_nonlinear(
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
                        vec[d] = 0.0;
                    }
                }

                int index = fe->conn[e][i];

                int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];
                int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
                int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                for(int d=0; d<4; d++) {
                    integ_val_vec[d] = BBFE_std_integ_calc(
                            np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                    for(int k = IS; k < IE; k++){
                        hlpod_ddhr->matrix[ns*(hlpod_vals->n_neib_vec) + k][m][n] = integ_val_vec[d] * hlpod_mat->neib_vec[index*4 + d][k];
                        hlpod_ddhr->RH[ns*(hlpod_vals->n_neib_vec) + k][n] = integ_val_vec[d] * hlpod_mat->neib_vec[index*4 + d][k];
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


void HROM_ddecm_set_residuals_NR_blas2_decoupled(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t)
{
    (void)num_snapshot;
    (void)num_neib;
    (void)dt;
    (void)t;

    const int ns = index_snap;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int nrv = hlpod_vals->n_neib_vec;

    /* ---------------------------------
     * compact coefficient:
     * original:
     *   index2: compact j-basis index
     *   index1+kj: padded coefficient index
     * --------------------------------- */
    double* coef_compact = (double*)mkl_malloc(sizeof(double) * (size_t)nrv, 64);
    {
        int index1 = 0;
        int index2 = 0;
        for (int ki = 0; ki < num_neib; ++ki) {
            for (int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; ++kj) {
                coef_compact[index2++] = hlpod_mat->mode_coef_1stdd[index1 + kj];
            }
            index1 += hlpod_mat->max_num_neib_modes[ki];
        }
    }

    /* ---- geometry / field work ---- */
    double*   Jacobian_ip = BB_std_calloc_1d_double(NULL, np);

    double**  local_v     = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip        = BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip   = BB_std_calloc_3d_double(NULL, np, 3, 3);

    double**  local_v_old = BB_std_calloc_2d_double(NULL, nl, 3);
    double**  v_ip_old    = BB_std_calloc_2d_double(NULL, np, 3);

    double*   local_p     = BB_std_calloc_1d_double(NULL, nl);
    double*   p_ip        = BB_std_calloc_1d_double(NULL, np);
    double**  grad_p_ip   = BB_std_calloc_2d_double(NULL, np, 3);

    double*   wJ          = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau         = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);
    double*   tau_c       = (double*)mkl_malloc(sizeof(double) * (size_t)np, 64);

    /* ---- element-local metadata ---- */
    int* conn_e = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IS_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* IE_i   = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);
    int* msz_i  = (int*)mkl_malloc(sizeof(int) * (size_t)nl, 64);

    /* j-side packed info */
    double* s_j = (double*)mkl_malloc(sizeof(double) * (size_t)nl * 4, 64);
    /* s_j[j*4+b] :
     *   free column -> dot(phi_j(:,b), coef_compact)
     *   D column    -> imposed_D_val
     */

    /* i-side basis pack: column-major (msz_max x 4) */
    int msz_max_global = 0;
    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];
                msz_max_global = max2_int(msz_max_global, IE - IS);
            }
        }
    }

    double* Phi_i_buf = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global * 4, 64);
    double* y_buf     = (double*)mkl_malloc(sizeof(double) * (size_t)msz_max_global, 64);

    for (int n = 0; n < num_subdomains; ++n) {
        for (int m = 0; m < hlpod_ddhr->num_elems[n]; ++m) {
            const int e = hlpod_ddhr->elem_id_local[m][n];

            /* -----------------------------
             * element precompute
             * ----------------------------- */
            BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

            BBFE_elemmat_set_local_array_vector(local_v,     fe, vals->v,     e, 3);
            BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);
            BBFE_elemmat_set_local_array_scalar(local_p,     fe, vals->p,     e);

            for (int p = 0; p < np; ++p) {
                double J_inv[3][3];

                BBFE_std_mapping_vector3d(v_ip[p],     nl, local_v,     basis->N[p]);
                BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
                BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
                BBFE_std_mapping_scalar_grad(grad_p_ip[p],  nl, local_p, fe->geo[e][p].grad_N);
                p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                tau[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                tau_c[p] = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                wJ[p] = basis->integ_weight[p] * Jacobian_ip[p];
            }

            /* -----------------------------
             * element node metadata
             * ----------------------------- */
            int msz_max_elem = 0;
            for (int i = 0; i < nl; ++i) {
                const int gi   = fe->conn[e][i];
                const int sidi = hlpod_mat->subdomain_id_in_nodes[gi];
                const int IS   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi];
                const int IE   = hlpod_ddhr->num_neib_modes_1stdd_sum[sidi + 1];

                conn_e[i] = gi;
                IS_i[i]   = IS;
                IE_i[i]   = IE;
                msz_i[i]  = IE - IS;

                msz_max_elem = max2_int(msz_max_elem, msz_i[i]);
            }

            if (msz_max_elem <= 0) continue;

            /* -----------------------------
             * j-side scalar summary:
             *   free  -> ddot(phi_j, coef_compact)
             *   D-bc  -> imposed_D_val
             * ----------------------------- */
            for (int j = 0; j < nl; ++j) {
                const int gj = conn_e[j];
                for (int b = 0; b < 4; ++b) {
                    if (bc->D_bc_exists[gj*4 + b]) {
                        s_j[j*4 + b] = bc->imposed_D_val[gj*4 + b];
                    }
                    else {
                        s_j[j*4 + b] = cblas_ddot(
                            nrv,
                            hlpod_mat->neib_vec_decoupled_v[gj*4 + b], 1,
                            coef_compact, 1);
                    }
                }
            }

            /* -----------------------------
             * (i,j) loop
             * ----------------------------- */
            for (int i = 0; i < nl; ++i) {
                const int gi  = conn_e[i];
                const int IS  = IS_i[i];
                const int IE  = IE_i[i];
                const int msz = msz_i[i];

                if (msz <= 0) continue;

                /* Phi_i_buf: (msz x 4), column-major */
                for (int a = 0; a < 4; ++a) {
                    memcpy(Phi_i_buf + (size_t)a * (size_t)msz_max_global,
                           &hlpod_mat->neib_vec_decoupled_v[gi*4 + a][IS],
                           sizeof(double) * (size_t)msz);
                }

                for (int j = 0; j < nl; ++j) {
                    const int gj = conn_e[j];

                    /* acc[a + 4*b] : column-major 4x4 */
                    double acc[16];
                    for (int q = 0; q < 16; ++q) acc[q] = 0.0;

                    for (int p = 0; p < np; ++p) {
                        double A[4][4];
                        double du_time[3] = {
                            v_ip[p][0] - v_ip_old[p][0],
                            v_ip[p][1] - v_ip_old[p][1],
                            v_ip[p][2] - v_ip_old[p][2]
                        };
                        double J_inv[3][3];

                        BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity,
                            tau[p], tau_c[p], vals->dt, du_time);

                        const double wp = wJ[p];
                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                acc[a + 4*b] += A[a][b] * wp;
                            }
                        }
                    }

                    /* g[a] = sum_b acc[a,b] * s_jb */
                    double g[4];
                    for (int a = 0; a < 4; ++a) {
                        g[a] =
                            acc[a + 4*0] * s_j[j*4 + 0] +
                            acc[a + 4*1] * s_j[j*4 + 1] +
                            acc[a + 4*2] * s_j[j*4 + 2] +
                            acc[a + 4*3] * s_j[j*4 + 3];
                    }

                    /* y = Phi_i(msz x 4) * g(4) */
                    cblas_dgemv(CblasColMajor, CblasNoTrans,
                                msz, 4,
                                1.0,
                                Phi_i_buf, msz_max_global,
                                g, 1,
                                0.0,
                                y_buf, 1);

                    /* scatter add */
                    for (int rr = 0; rr < msz; ++rr) {
                        const int index = ns * nrv + (IS + rr);
                        const double val = y_buf[rr];

                        hlpod_ddhr->matrix[index][m][n] += val;
                        hlpod_ddhr->RH[index][n]        += val;
                    }
                }
            }
        }
    }

    /* ---- free ---- */
    mkl_free(y_buf);
    mkl_free(Phi_i_buf);

    mkl_free(s_j);

    mkl_free(conn_e);
    mkl_free(IS_i);
    mkl_free(IE_i);
    mkl_free(msz_i);

    mkl_free(wJ);
    mkl_free(tau);
    mkl_free(tau_c);

    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    mkl_free(coef_compact);
}



void HROM_set_element_mat_NR_decoupled_p(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*    hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
    /*
    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }
    */

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


    int rank = monolis_mpi_get_global_my_rank();

   
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems_p
            : hlpod_ddhr->ovl_num_selected_elems_D_bc_p;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems_p
            : hlpod_ddhr->ovl_id_selected_elems_D_bc_p;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight_p
            : hlpod_ddhr->ovl_elem_weight_D_bc_p;

        for (int m = 0; m < num; m++) {
            int e = ids[m];
 
 //   for (int e = 0; e < fe->total_num_elems; ++e) {
        
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

                        /*
                        BBFE_elemmat_fluid_mat_rom_linear(
                            A, J_inv,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);
                        */

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                val_ip[a][b][p] = A[a][b];

                                if (!isfinite(A[a][b])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, a, A[a][b]);
                                }

                                A[a][b] = 0.0;
                            }
                        }
                    }

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        if( bc->D_bc_exists[index_j*4+b]) {
                        }
                        else{
                            int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                            int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                            int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                            int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                            int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                            int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                            int subdomain_id = hlpod_mat->subdomain_id_in_nodes_2nddd[index_j]; 

                            if(subdomain_id_j < num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS; k2 < JE; k2++){
                                        double val = hlpod_mat->pod_basis_hr_decoupled_p[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr_decoupled_p[index_j*4+b][k2];
                                        //printf("val = %e", val);
                                        hlpod_ddhr->reduced_mat[k1][k2] += weight[m] *  val;
                                        //hlpod_ddhr->reduced_mat[k1][k2] += val;
                                    }
                                }
                            }

                            if(subdomain_id_j >= num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
                                        double val = hlpod_mat->pod_basis_hr_decoupled_p[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr_decoupled_p[index_j*4+b][k2] ;
                                        //hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += weight[m] *  val;
                                        hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += val;
                                        // /printf("val = %e", val);
                                    }
                                }
                            }
                        }
                    }
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


void HROM_set_element_mat_NR_decoupled_v(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*   hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

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


    int rank = monolis_mpi_get_global_my_rank();

    
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems_v
            : hlpod_ddhr->ovl_num_selected_elems_D_bc_v;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems_v
            : hlpod_ddhr->ovl_id_selected_elems_D_bc_v;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight_v
            : hlpod_ddhr->ovl_elem_weight_D_bc_v;

        for (int m = 0; m < num; m++) {
            int e = ids[m];
 
    //for (int e = 0; e < fe->total_num_elems; ++e) {
        
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

                        /*
                        BBFE_elemmat_fluid_mat_rom_linear(
                            A, J_inv,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);
                        */

                        BBFE_elemmat_fluid_mat_rom_nonlinear(
                            A, J_inv, fe->geo[e][p].Jacobian,
                            basis->N[p][i], basis->N[p][j],
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                            v_ip[p], grad_v_ip[p], grad_p_ip[p],
                            vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                        for (int a = 0; a < 4; ++a) {
                            for (int b = 0; b < 4; ++b) {
                                val_ip[a][b][p] = A[a][b];

                                if (!isfinite(A[a][b])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, a, A[a][b]);
                                }

                                A[a][b] = 0.0;
                            }
                        }
                    }

                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double integ_val = BBFE_std_integ_calc(
                                np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        if( bc->D_bc_exists[index_j*4+b]) {
                        }
                        else{
                            int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                            int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                            int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                            int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                            int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                            int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                            int subdomain_id = hlpod_mat->subdomain_id_in_nodes_2nddd[index_j]; 

                            if(subdomain_id_j < num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS; k2 < JE; k2++){
                                        double val = hlpod_mat->pod_basis_hr_decoupled_v[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr_decoupled_v[index_j*4+b][k2];
                                        //printf("val = %e", val);
                                        hlpod_ddhr->reduced_mat[k1][k2] += weight[m] *  val;
                                        //hlpod_ddhr->reduced_mat[k1][k2] += val;
                                    }
                                }
                            }

                            if(subdomain_id_j >= num_subdomains){
                                for(int k1 = IS; k1 < IE; k1++){
                                    for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
                                        double val = hlpod_mat->pod_basis_hr_decoupled_v[index_i*4+a][k1] * integ_val * hlpod_mat->pod_basis_hr_decoupled_v[index_j*4+b][k2] ;
                                        hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += weight[m] *  val;
                                        //hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += val;
                                        // /printf("val = %e", val);
                                    }
                                }
                            }
                        }
                    }
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


void HROM_set_element_vec_NR_v(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
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

    int rank = monolis_mpi_get_global_my_rank();
    
    
    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems_v
            : hlpod_ddhr->ovl_num_selected_elems_D_bc_v;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems_v
            : hlpod_ddhr->ovl_id_selected_elems_D_bc_v;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight_v
            : hlpod_ddhr->ovl_elem_weight_D_bc_v;

         printf("v num = %d, %d\n", num, rank);
          

        for (int m = 0; m < num; m++) {
            int e = ids[m];
        
    
        //for(int e = 0; e < fe->total_num_elems; e++) {

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
            
            for (int i = 0; i < nl; ++i) {

                int index = fe->conn[e][i];
                int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

                if(subdomain_id < num_subdomains){

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
                    
                        BBFE_elemmat_fluid_vec_rom_nonlinear(
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

                                if (!isfinite(vec[d])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, d, vec[d]);
                                }

                            vec[d] = 0.0;
                        }
                    }

                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                    for (int d = 0; d < 4; ++d) {
                        integ_val_vec[d] = BBFE_std_integ_calc(
                            np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);
                        
                            if (!isfinite(integ_val_vec[d])) {
                                printf("rank = %d, bad integ_val_vec: e=%d i=%d d=%d val=%e\n", rank, e, i, d, integ_val_vec[d]);
                            }      
                        
                        for(int k = IS; k < IE; k++){
                            double val = integ_val_vec[d] * hlpod_mat->pod_modes_decoupled_v[index*4+d][k];

                            monolis->mat.R.B[k] -= weight[m] * val;
                            //monolis->mat.R.B[k] -= val;
                        
                            double mode = hlpod_mat->pod_modes_decoupled_v[index*4+d][k];
                            double val2 = integ_val_vec[d] * mode;

                            if (!isfinite(mode) || !isfinite(val)) {
                                printf("rank = %d, bad reduction: e=%d i=%d index=%d d=%d k=%d integ=%e mode=%e val=%e\n",
                                    rank, e, i, index, d, k, integ_val_vec[d], mode, val);
                            }

                            if(val==0.0){
                            }
                            else{
                            //hlpod_ddhr->reduced_RH[k] -= val;
                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                            }

                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                        }
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


void HROM_set_element_vec_NR_p(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
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

    int rank = monolis_mpi_get_global_my_rank();
    

    for (int kind = 0; kind < 2; kind++) {
        int num = (kind == 0)
            ? hlpod_ddhr->ovl_num_selected_elems_p
            : hlpod_ddhr->ovl_num_selected_elems_D_bc_p;

        int *ids = (kind == 0)
            ? hlpod_ddhr->ovl_id_selected_elems_p
            : hlpod_ddhr->ovl_id_selected_elems_D_bc_p;
        
        double *weight = (kind == 0)
            ? hlpod_ddhr->ovl_elem_weight_p
            : hlpod_ddhr->ovl_elem_weight_D_bc_p;

        printf("p num = %d, %d\n", num, rank);
            

        for (int m = 0; m < num; m++) {
            int e = ids[m];
    
        //for(int e = 0; e < fe->total_num_elems; e++) {

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
            
            for (int i = 0; i < nl; ++i) {

                int index = fe->conn[e][i];
                int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

                if(subdomain_id < num_subdomains){

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
                    
                        BBFE_elemmat_fluid_vec_rom_nonlinear(
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

                                if (!isfinite(vec[d])) {
                                    printf("rank = %d, bad vec: e=%d p=%d i=%d d=%d vec=%e\n", rank, e, p, i, d, vec[d]);
                                }

                            vec[d] = 0.0;
                        }
                    }

                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                    for (int d = 0; d < 4; ++d) {
                        integ_val_vec[d] = BBFE_std_integ_calc(
                            np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);
                        
                            if (!isfinite(integ_val_vec[d])) {
                                printf("rank = %d, bad integ_val_vec: e=%d i=%d d=%d val=%e\n", rank, e, i, d, integ_val_vec[d]);
                            }      
                        
                        for(int k = IS; k < IE; k++){
                            double val = integ_val_vec[d] * hlpod_mat->pod_modes_decoupled_p[index*4+d][k];

                            monolis->mat.R.B[k] -= weight[m] * val;
                            //monolis->mat.R.B[k] -= val;
                        
                            double mode = hlpod_mat->pod_modes_decoupled_p[index*4+d][k];
                            double val2 = integ_val_vec[d] * mode;

                            if (!isfinite(mode) || !isfinite(val)) {
                                printf("rank = %d, bad reduction: e=%d i=%d index=%d d=%d k=%d integ=%e mode=%e val=%e\n",
                                    rank, e, i, index, d, k, integ_val_vec[d], mode, val);
                            }

                            if(val==0.0){
                            }
                            else{
                            //hlpod_ddhr->reduced_RH[k] -= val;
                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                            }

                            //printf("val = %e integ_val_vec = %e modes = %e", val, integ_val_vec[d], hlpod_mat->pod_modes[index*4+d][k]);
                        }
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