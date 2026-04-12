
#include "team7.h"

/* TEAM Problem 7 */
const double Sigma_coil = 5.7e7;      /* impressed current regionとして扱うなら 0 でよい */
const double Sigma_al   = 3.526e7;  /* aluminum plate */
const double Sigma_air  = 0.0;

void get_sigmas_for_prop_team7(
    int prop,
    double* sigma_mass_A,
    double* sigma_cpl,
    double* sigma_phi
){
    if(prop == 1 || prop == 2){
        /* exciting coil */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;

    } else if(prop == 3){
        /* aluminum conductor with hole */
        *sigma_mass_A = Sigma_al;
        *sigma_cpl    = Sigma_al;
        *sigma_phi    = Sigma_al;

    } else if(prop == 4){
        /* hole or unused region */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;

    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


const double mu0   = 4.0*M_PI*1e-7;      // H/m
const double Nu    = 1.0 / mu0;          // 1/H/m = A/(T·m)
const double freq  = 50.0;               // Hz
const double Omega = 2.0*M_PI*freq;      // rad/s  ★これが正しい


static const double I_RMS_team7 = 10.0;   /* TEAM benchmark rated current [A rms] */
static const double FREQ_HZ_team7 = 50.0; /* [Hz] */
static const double NI_RMS_team7 = 2742.0;

static inline double get_coil_current_team7(int prop)
{
    const double omega = 2.0 * M_PI * FREQ_HZ_team7;

    if(prop == 1){
        return NI_RMS_team7;   
    }

    return 0.0;
}

void set_element_mat_nedelec_Aphi_team7(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    (void)bc;

    const int np = basis->num_integ_points;

    double* J_ip = BB_std_calloc_1d_double(J_ip, np);
    double _Complex* val_ip_C = BB_std_calloc_1d_double_C(val_ip_C, np);

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];
        double sigma_massA, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_team7(prop, &sigma_massA, &sigma_cpl, &sigma_phi);

        int nonlinear_mu = 0;

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j=0; j<ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += Nu * BBFE_elemmat_mag_mat_curl(
                        ned->curl_N_edge[e][p][i],
                        ned->curl_N_edge[e][p][j],
                        1.0
                    );
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)(si*sj);  // ★ edge_sign

                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        //if(prop==1){
            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int j=0; j<ned->local_num_edges; ++j){
                    const int gj = ned->nedelec_conn[e][j];
                    const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = 0.0 + 0.0*I;
                        val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            ned->N_edge[e][p][j],
                            1
                        );
                    }

                    double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                    v *= (double)(si*sj);
                    monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
                }
            }
        //}

        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j=0; j<ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        ned->N_edge[e][p][j],
                        sigma_massA
                    ) * Omega * I;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)(si*sj);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm =  fe->conn[e][m];

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn =  fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[m],
                        fe->geo[e][p].grad_N[n],
                        sigma_phi
                    )*I*Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
            }
        }

        for(int j=0; j<ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn =  fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[n],
                        ned->N_edge[e][p][j],
                        sigma_cpl
                    ) *I*Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)sj;

                monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
            }
        }

        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm =  fe->conn[e][m];

            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        fe->geo[e][p].grad_N[m],
                        sigma_cpl
                    )*I*Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)si;

                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gi, 0, 0, v);
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
    BB_std_free_1d_double_C(val_ip_C, np);
}

static inline int get_coil_info_team7(int elem_prop, COIL_INFO* info)
{
    if(elem_prop != 1) return 0;

    info->axis[0] = 0.0;
    info->axis[1] = 0.0;
    info->axis[2] = 1.0;

    info->turns = 1.0;  /* NI を直接使うなら 1 でよい */

    info->area = 2500.0 * (MM_TO_M * MM_TO_M);

    /* race-track の代表中心 */
    info->center[0] = 0.5 * (44.0 + 294.0) * MM_TO_M;
    info->center[1] = 0.5 * (0.0  + 250.0) * MM_TO_M;
    info->center[2] = (19.0 + 30.0 + 50.0) * MM_TO_M;

    return 1;
}

static int get_team7_rectcoil_tangent(
    const COIL_INFO* coil,
    const double x_ip[3],
    double tdir[3]
){
    const double MM = MM_TO_M;

    /* TEAM13 の coding と同じ発想:
       コイル電流密度 J を「直線部 + 4つの角部」で区分的に与える。
       ただし寸法は TEAM7 の race-track 幾何に合わせる。 */

    const double outerX = 270.0 * MM;
    const double outerY = 270.0 * MM;
    const double innerX = 200.0 * MM;
    const double innerY = 200.0 * MM;
    const double rOuter = 45.0  * MM;
    const double rInner = 10.0  * MM;

    /* centerline rounded rectangle */
    const double a  = 0.25 * (outerX + innerX);
    const double b  = 0.25 * (outerY + innerY);
    const double rc = 0.5  * (rOuter + rInner);
    const double xs = a - rc;
    const double ys = b - rc;

    const double x = x_ip[0] - coil->center[0];
    const double y = x_ip[1] - coil->center[1];
    const double eps = 1.0e-14;

    /* ---- straight bricks (TEAM13 style piecewise definition) ---- */
    if(y >=  ys && fabs(x) <= xs){
        tdir[0] = -1.0; tdir[1] =  0.0; tdir[2] = 0.0;
        return 1;
    }
    if(y <= -ys && fabs(x) <= xs){
        tdir[0] =  1.0; tdir[1] =  0.0; tdir[2] = 0.0;
        return 1;
    }
    if(x >=  xs && fabs(y) <= ys){
        tdir[0] =  0.0; tdir[1] =  1.0; tdir[2] = 0.0;
        return 1;
    }
    if(x <= -xs && fabs(y) <= ys){
        tdir[0] =  0.0; tdir[1] = -1.0; tdir[2] = 0.0;
        return 1;
    }

    /* ---- corner arcs (same idea as TEAM13: tangent from local radius) ---- */
    double cx = 0.0, cy = 0.0;
    int found_corner = 1;

    if(x >  xs && y >  ys){
        cx = +xs; cy = +ys;   /* top-right */
    } else if(x < -xs && y >  ys){
        cx = -xs; cy = +ys;   /* top-left */
    } else if(x < -xs && y < -ys){
        cx = -xs; cy = -ys;   /* bottom-left */
    } else if(x >  xs && y < -ys){
        cx = +xs; cy = -ys;   /* bottom-right */
    } else {
        found_corner = 0;
    }

    if(found_corner){
        const double xr = x - cx;
        const double yr = y - cy;
        const double r  = sqrt(xr*xr + yr*yr);

        if(r > eps){
            /* CCW tangent along rounded rectangle */
            tdir[0] = -yr / r;
            tdir[1] =  xr / r;
            tdir[2] =  0.0;
            return 1;
        }
    }

    /* フォールバック: 最近傍の中心線接線を返す */
    {
        double best_d2 = 1.0e30;
        double best_tx = 0.0;
        double best_ty = 0.0;

        #define UPDATE_BEST(px, py, tx_, ty_) do {             const double dx_ = x - (px);             const double dy_ = y - (py);             const double d2_ = dx_*dx_ + dy_*dy_;             if(d2_ < best_d2){                 best_d2 = d2_;                 best_tx = (tx_);                 best_ty = (ty_);             }         } while(0)

        UPDATE_BEST(fmax(-xs, fmin(xs, x)),  b, -1.0,  0.0);
        UPDATE_BEST(fmax(-xs, fmin(xs, x)), -b,  1.0,  0.0);
        UPDATE_BEST( a, fmax(-ys, fmin(ys, y)),  0.0,  1.0);
        UPDATE_BEST(-a, fmax(-ys, fmin(ys, y)),  0.0, -1.0);

        const double arc_centers[4][2] = {
            { +xs, +ys }, { -xs, +ys }, { -xs, -ys }, { +xs, -ys }
        };

        for(int k = 0; k < 4; ++k){
            const double ax = arc_centers[k][0];
            const double ay = arc_centers[k][1];
            const double vx = x - ax;
            const double vy = y - ay;
            const double vn = sqrt(vx*vx + vy*vy);
            if(vn <= eps) continue;

            const double ux = vx / vn;
            const double uy = vy / vn;
            UPDATE_BEST(ax + rc*ux, ay + rc*uy, -uy, ux);
        }

        #undef UPDATE_BEST

        if(best_d2 < 1.0e29){
            const double nt = sqrt(best_tx*best_tx + best_ty*best_ty);
            if(nt > eps){
                tdir[0] = best_tx / nt;
                tdir[1] = best_ty / nt;
                tdir[2] = 0.0;
                return 1;
            }
        }
    }

    tdir[0] = 0.0;
    tdir[1] = 0.0;
    tdir[2] = 0.0;
    return 0;
}


static int get_team7_rectcoil_current_density(
    const COIL_INFO* coil,
    const double x_ip[3],
    const double J_mag,
    double Js[3]
)
{
    const double MM = MM_TO_M;
    const double eps = 1.0e-14;

    /* TEAM13 coding に合わせて、J は「一定大きさ + 区分的な方向」 */
    const double outerX = 270.0 * MM;
    const double outerY = 270.0 * MM;
    const double innerX = 200.0 * MM;
    const double innerY = 200.0 * MM;
    const double rOuter = 45.0  * MM;
    const double rInner = 10.0  * MM;

    const double a  = 0.25 * (outerX + innerX);
    const double b  = 0.25 * (outerY + innerY);
    const double rc = 0.5  * (rOuter + rInner);
    const double xs = a - rc;
    const double ys = b - rc;

    const double x = x_ip[0] - coil->center[0];
    const double y = x_ip[1] - coil->center[1];

    Js[0] = 0.0;
    Js[1] = 0.0;
    Js[2] = 0.0;

    /* ---- straight bricks: same idea as TEAM13 ---- */
    if(y >=  ys && fabs(x) <= xs){
        Js[0] = -J_mag; Js[1] =  0.0;   Js[2] = 0.0;
        return 1;
    }
    if(y <= -ys && fabs(x) <= xs){
        Js[0] =  J_mag; Js[1] =  0.0;   Js[2] = 0.0;
        return 1;
    }
    if(x >=  xs && fabs(y) <= ys){
        Js[0] =  0.0;   Js[1] =  J_mag; Js[2] = 0.0;
        return 1;
    }
    if(x <= -xs && fabs(y) <= ys){
        Js[0] =  0.0;   Js[1] = -J_mag; Js[2] = 0.0;
        return 1;
    }

    /* ---- corner arcs: same formula pattern as TEAM13 ---- */
    {
        double cx = 0.0, cy = 0.0;
        int found_corner = 1;

        if(x >  xs && y >  ys){
            cx = +xs; cy = +ys;   /* top-right */
        } else if(x < -xs && y >  ys){
            cx = -xs; cy = +ys;   /* top-left */
        } else if(x < -xs && y < -ys){
            cx = -xs; cy = -ys;   /* bottom-left */
        } else if(x >  xs && y < -ys){
            cx = +xs; cy = -ys;   /* bottom-right */
        } else {
            found_corner = 0;
        }

        if(found_corner){
            const double xr = x - cx;
            const double yr = y - cy;
            const double r  = sqrt(xr*xr + yr*yr);

            if(r > eps){
                Js[0] = J_mag * (-yr / r);
                Js[1] = J_mag * ( xr / r);
                Js[2] = 0.0;
                return 1;
            }
        }
    }

    /* inside prop==1 but outside the idealized subregions: use closest centerline tangent */
    {
        double best_d2 = 1.0e100;
        double best_tx = 0.0;
        double best_ty = 0.0;

        #define UPDATE_BEST(px, py, tx_, ty_) do { \
            const double dx_ = x - (px); \
            const double dy_ = y - (py); \
            const double d2_ = dx_*dx_ + dy_*dy_; \
            if(d2_ < best_d2){ \
                best_d2 = d2_; \
                best_tx = (tx_); \
                best_ty = (ty_); \
            } \
        } while(0)

        UPDATE_BEST(fmax(-xs, fmin(xs, x)),  b, -1.0,  0.0);
        UPDATE_BEST(fmax(-xs, fmin(xs, x)), -b,  1.0,  0.0);
        UPDATE_BEST( a, fmax(-ys, fmin(ys, y)),  0.0,  1.0);
        UPDATE_BEST(-a, fmax(-ys, fmin(ys, y)),  0.0, -1.0);

        const double arc_centers[4][2] = {
            { +xs, +ys }, { -xs, +ys }, { -xs, -ys }, { +xs, -ys }
        };

        for(int k = 0; k < 4; ++k){
            const double ax = arc_centers[k][0];
            const double ay = arc_centers[k][1];
            const double vx = x - ax;
            const double vy = y - ay;
            const double vn = sqrt(vx*vx + vy*vy);
            if(vn <= eps) continue;

            const double ux = vx / vn;
            const double uy = vy / vn;

            int ok = 0;
            if(k == 0) ok = (ux >= 0.0 && uy >= 0.0);
            if(k == 1) ok = (ux <= 0.0 && uy >= 0.0);
            if(k == 2) ok = (ux <= 0.0 && uy <= 0.0);
            if(k == 3) ok = (ux >= 0.0 && uy <= 0.0);
            if(!ok) continue;

            UPDATE_BEST(ax + rc*ux, ay + rc*uy, -uy, ux);
        }

        #undef UPDATE_BEST

        if(best_d2 < 1.0e99){
            const double nt = sqrt(best_tx*best_tx + best_ty*best_ty);
            if(nt > eps){
                Js[0] = J_mag * (best_tx / nt);
                Js[1] = J_mag * (best_ty / nt);
                Js[2] = 0.0;
                return 1;
            }
        }
    }

    return 0;
}

void set_element_vec_nedelec_Aphi_team7(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double epsB   = 1.0e-14;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e = 0; e < fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_team7(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info_team7(prop, &coil);

    
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        if(is_coil){
            double NI_t  = get_coil_current_team7(prop);
            double J_mag = NI_t / coil.area;

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double Js[3] = {0.0, 0.0, 0.0};
                    get_team7_rectcoil_current_density(&coil, x_ip, J_mag, Js);

                    val_ip_C[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis->mat.C.B[gi] -= (double)si * integ;
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}


void apply_dirichlet_bc_for_A_and_phi_team7(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned)
{
    const int nen = fe->local_num_nodes;

    int is_dir_edge_n = fe->total_num_nodes;
    bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
    build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);

    /* --------------------------------------------------------
       A (edge) boundary condition: A_tan = 0 on outer boundary
       -------------------------------------------------------- */
    int n_local_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if(nen == 4){
        n_local_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if(nen == 8){
        n_local_edges = 12;
        edge_tbl = hex_edge_conn;
    }

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int i = 0; i < n_local_edges; ++i){
            int gn1 = fe->conn[e][edge_tbl[i][0]];
            int gn2 = fe->conn[e][edge_tbl[i][1]];
            int ged = ned->nedelec_conn[e][i];

            if(!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;
            if(!is_dir_edge[ged]) continue;

            monolis_set_Dirichlet_bc_C(
                monolis,
                monolis->mat.C.B,
                ged,
                0,
                0.0 + 0.0*I
            );
        }
    }

    
    int num_nodes = fe->total_num_nodes;
    int* node_is_conductor = (int*)calloc(num_nodes, sizeof(int));

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop == 3){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 3;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop==1){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 1;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 2){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 2; 
            }
        }
    }

    for (int i = 0; i < num_nodes; ++i){
        //if (node_is_conductor[i] == 1||node_is_conductor[i] == 3) {
        if (node_is_conductor[i] == 1||node_is_conductor[i] == 2||node_is_conductor[i] == 4 ||node_is_conductor[i] == 3) {
	//if (node_is_conductor[i] == 4) {


	    monolis_set_Dirichlet_bc_C(
                monolis, 
                monolis->mat.C.B, 
                i, 
                0, 
                0.0 + 0.0*I
            );
        }  
    }
    

    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);
}

