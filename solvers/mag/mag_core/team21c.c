#include "team21c.h"

/* Air stabilization factors (tune as needed) */
static const double AIR_PHI_SIGMA_FACTOR  = 0.0; /* keep phi disabled in air for P21C-EM1 */
static const double AIR_A_MASS_SIGMA_FACTOR = 0.0; /* keep A-mass disabled in air for P21C-EM1 */

const double Sigma_coil   = 5.7143e7;
const double Sigma_shield = 5.7143e7;
const double Sigma_steel  = 1.3889e6;

void get_sigmas_for_prop_team21c(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2){
        /* exciting coil conductor */
        *sigma_mass_A = Sigma_coil * 0.1256/3.2;
        //*sigma_cpl    = Sigma_coil* 0.119;
        //*sigma_phi    = Sigma_coil* 0.119;
	    *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 3){
        /* TEAM P21C-EM1 copper shielding plate */
        *sigma_mass_A = Sigma_steel;
        *sigma_cpl    = Sigma_steel;
        *sigma_phi    = Sigma_steel;
        //*sigma_cpl    = 0.0;
        //*sigma_phi    = 0.0;
    } else if(prop == 4){
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 5){
        /* air */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}

void get_sigmas_for_prop_team21a(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1  || prop == 2){
        /* exciting coil conductor */
        *sigma_mass_A = Sigma_coil*1.0e-2;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 3){
        /* TEAM P21C-EM1 copper shielding plate */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 4){
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 5){
        /* air */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


static inline int get_coil_info_team21(int elem_prop, COIL_INFO* info) {
    info->axis[0] = 0.0;
    info->axis[1] = 0.0;
    info->axis[2] = 1.0;
    info->turns   = 300.0 * 0.7;

    /* FE coil-region equivalent cross-sectional area */
    //info->area    = (32900.0) * (MM_TO_M * MM_TO_M);
    //info->area    = (32900.0) * (MM_TO_M * MM_TO_M) * 0.119;
    //info->area    = (31144) * (MM_TO_M * MM_TO_M) * 0.12570960698;
    info->area    = (3912) * (MM_TO_M * MM_TO_M);

    {
        const double cy = 0.0 * MM_TO_M;
        const double cz = 0.0 * MM_TO_M;

        if (elem_prop == 1) {
            info->center[0] = 0.0;
            info->center[1] = cy;
            info->center[2] = -120.5 * MM_TO_M; /* coil1 center z */
            return 1;
        }
        if (elem_prop == 2) {
            info->center[0] = 0.0;
            info->center[1] = cy;
            info->center[2] =  120.5 * MM_TO_M; /* coil2 center z */
            return 1;
        }
    }
    return 0;
}

#include <float.h>

/* ==========================================================
 * TEAM 21c rounded-rectangle coil tangent model in local xy-plane
 * Geometry must match the gmsh model:
 *   outer 270 x 270 mm, corner radius 45 mm
 *   inner 200 x 200 mm, corner radius 10 mm
 *
 * We use the FE coil region centerline approximation:
 *   outer/inner の中間の角丸矩形を1本の閉ループとして扱い，
 *   その接線方向を Js の方向に使う。
 * ========================================================== */
static int get_team21c_rectcoil_tangent(
    const COIL_INFO* coil,
    const double x_ip[3],
    double tdir[3]
){
    const double MM = MM_TO_M;

    /* ---- centerline geometry (midline of coil pack) ---- */
    const double outerX = 270.0 * MM;
    const double outerY = 270.0 * MM;
    const double innerX = 200.0 * MM;
    const double innerY = 200.0 * MM;
    const double rOuter = 45.0 * MM;
    const double rInner = 10.0 * MM;

    /* centerline rounded rectangle */
    const double a = 0.25 * (outerX + innerX);   /* half width  = 117.5 mm */
    const double b = 0.25 * (outerY + innerY);   /* half height = 117.5 mm */
    const double rc = 0.5  * (rOuter + rInner);  /* corner rad  = 27.5 mm */

    /* straight-part limits of centerline */
    const double xs = a - rc;  /* 90.0 mm */
    const double ys = b - rc;  /* 90.0 mm */

    /* local coordinates around coil center */
    const double x = x_ip[0] - coil->center[0];
    const double y = x_ip[1] - coil->center[1];

    /* candidates: nearest point on 4 straight segments + 4 corner arcs */
    double best_d2 = DBL_MAX;
    double best_tx = 0.0, best_ty = 0.0;

    /* utility macro */
    #define UPDATE_BEST(px, py, tx_, ty_) do { \
        double dx_ = x - (px); \
        double dy_ = y - (py); \
        double d2_ = dx_*dx_ + dy_*dy_; \
        if(d2_ < best_d2){ \
            best_d2 = d2_; \
            best_tx = (tx_); \
            best_ty = (ty_); \
        } \
    } while(0)

    /* ---------------- 4 straight segments ----------------
     * CCW tangent:
     *   top    : (-1,  0)
     *   left   : ( 0, -1)
     *   bottom : ( 1,  0)
     *   right  : ( 0,  1)
     */

    /* top: y = +b, x in [-xs, +xs] */
    {
        double px = fmax(-xs, fmin(xs, x));
        double py = b;
        UPDATE_BEST(px, py, -1.0,  0.0);
    }

    /* bottom: y = -b, x in [-xs, +xs] */
    {
        double px = fmax(-xs, fmin(xs, x));
        double py = -b;
        UPDATE_BEST(px, py,  1.0,  0.0);
    }

    /* right: x = +a, y in [-ys, +ys] */
    {
        double px = a;
        double py = fmax(-ys, fmin(ys, y));
        UPDATE_BEST(px, py,  0.0,  1.0);
    }

    /* left: x = -a, y in [-ys, +ys] */
    {
        double px = -a;
        double py = fmax(-ys, fmin(ys, y));
        UPDATE_BEST(px, py,  0.0, -1.0);
    }

    /* ---------------- 4 corner arcs ----------------
     * Arc centers:
     *   TR: (+xs, +ys)
     *   TL: (-xs, +ys)
     *   BL: (-xs, -ys)
     *   BR: (+xs, -ys)
     *
     * For a CCW loop, tangent = (-uy, ux),
     * where u = radial unit vector from arc center to closest point.
     */

    struct ArcInfo { double cx, cy; };
    const struct ArcInfo arcs[4] = {
        { +xs, +ys },  /* TR */
        { -xs, +ys },  /* TL */
        { -xs, -ys },  /* BL */
        { +xs, -ys }   /* BR */
    };

    for(int k = 0; k < 4; ++k){
        double cx = arcs[k].cx;
        double cy = arcs[k].cy;

        double vx = x - cx;
        double vy = y - cy;
        double vn = sqrt(vx*vx + vy*vy);

        if(vn < 1.0e-20){
            continue;
        }

        double ux = vx / vn;
        double uy = vy / vn;

        /* nearest point on full circle of radius rc */
        double px = cx + rc * ux;
        double py = cy + rc * uy;

        /* keep only the quarter-arc belonging to each corner */
        int ok = 0;
        if(k == 0) ok = (ux >= 0.0 && uy >= 0.0); /* TR */
        if(k == 1) ok = (ux <= 0.0 && uy >= 0.0); /* TL */
        if(k == 2) ok = (ux <= 0.0 && uy <= 0.0); /* BL */
        if(k == 3) ok = (ux >= 0.0 && uy <= 0.0); /* BR */

        if(!ok) continue;

        /* CCW tangent */
        double tx = -uy;
        double ty =  ux;

        UPDATE_BEST(px, py, tx, ty);
    }

    #undef UPDATE_BEST

    /* map local xy tangent back to global xyz
       axis is z for this geometry */
    double nt = sqrt(best_tx*best_tx + best_ty*best_ty);
    if(nt < 1.0e-20){
        tdir[0] = 0.0;
        tdir[1] = 0.0;
        tdir[2] = 0.0;
        return 0;
    }

    tdir[0] = best_tx / nt;
    tdir[1] = best_ty / nt;
    tdir[2] = 0.0;
    return 1;
}

static const int BH_A3_N = 29;

static const double BH_A3_B[BH_A3_N] = {
    0.049, 0.101, 0.150, 0.200, 0.299, 0.399, 0.499, 0.601, 0.700, 0.801,
    0.899, 1.001, 1.099, 1.200, 1.300, 1.401, 1.449, 1.500, 1.550, 1.600,
    1.639, 1.670, 1.701, 1.729, 1.760, 1.781, 1.800, 1.830, 1.850
};

static const double BH_A3_H[BH_A3_N] = {
     115,  171,  196,  214,  245,  279,  316,  359,  405,  461,
     528,  616,  732,  898, 1154, 1606, 1965, 2506, 3291, 4430,
    5599, 6698, 7926, 9251, 10792, 11930, 13106, 14949, 16290
};

static void eval_nu_and_dnudB_team21c(double Bmag, double* nu, double* dnudB)
{
    const double mu0 = 4.0e-7 * M_PI;
    const double eps = 1.0e-12;

    if (Bmag <= BH_A3_B[0]) {
        double H = BH_A3_H[0] * (Bmag / BH_A3_B[0]);   /* 原点-第1点の線形補間 */
        double dH_dB = BH_A3_H[0] / BH_A3_B[0];

        if (Bmag < eps) {
            *nu = dH_dB;
            *dnudB = 0.0;
        } else {
            *nu = H / Bmag;
            *dnudB = 0.0;   /* この区間では一定勾配 */
        }
        return;
    }

    /* measured range: up to 1.85 T */
    if (Bmag <= BH_A3_B[BH_A3_N - 1]) {
        for (int k = 0; k < BH_A3_N - 1; ++k) {
            double B0 = BH_A3_B[k];
            double B1 = BH_A3_B[k + 1];
            double H0 = BH_A3_H[k];
            double H1 = BH_A3_H[k + 1];

            if (Bmag >= B0 && Bmag <= B1) {
                double t = (Bmag - B0) / (B1 - B0);
                double H = H0 + t * (H1 - H0);
                double dH_dB = (H1 - H0) / (B1 - B0);

                *nu = H / Bmag;
                *dnudB = (dH_dB * Bmag - H) / (Bmag * Bmag);
                return;
            }
        }
    }

    /* extrapolation above 1.85 T from the benchmark document */
    if (Bmag < 2.1) {
        double H = BH_A3_H[BH_A3_N - 1];

        for (int it = 0; it < 30; ++it) {
            double f =
                mu0 * H
                - 1.9538e-10 * H * H
                + 1.9043e-5 * H
                + 1.5729
                - Bmag;

            double df =
                mu0
                - 2.0 * 1.9538e-10 * H
                + 1.9043e-5;

            double dH = -f / df;
            H += dH;
            if (fabs(dH) < 1.0e-10) break;
        }

        {
            double dB_dH =
                mu0
                - 2.0 * 1.9538e-10 * H
                + 1.9043e-5;
            double dH_dB = 1.0 / dB_dH;

            *nu = H / Bmag;
            *dnudB = (dH_dB * Bmag - H) / (Bmag * Bmag);
            return;
        }
    } else {
        double H = (Bmag - 2.0368) / mu0;
        double dH_dB = 1.0 / mu0;

        *nu = H / Bmag;
        *dnudB = (dH_dB * Bmag - H) / (Bmag * Bmag);
        return;
    }
}

double get_reluctivity_nu_team21c(double Bmag)
{
    double nu, dnudB;
    eval_nu_and_dnudB_team21c(Bmag, &nu, &dnudB);
    return nu;
}

/* ============================================================
 * Local utility functions
 * ============================================================ */
static double local_sum_sq(const double* v, int n) {
    double s = 0.0;
    for(int i=0; i<n; ++i) s += v[i]*v[i];
    return s;
}
static double local_max_abs(const double* v, int n) {
    double m = 0.0;
    for(int i=0; i<n; ++i){
        double a = fabs(v[i]);
        if(a > m) m = a;
    }
    return m;
}

/* ============================================================
 * Jacobian Assembly (Newton)
 * ============================================================ */
void set_element_mat_NR_Aphi_team21c(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0/dt;
    const double epsB   = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];
        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_team21c(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ========= [1] A-A : Curl-Curl (tangent stiffness) ========= */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    if(nonlinear_mu){
                        double B[3]; compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        double Bmag = fmax(norm3(B), epsB);
                        double nu, dnudB;
                        eval_nu_and_dnudB_team21c(Bmag, &nu, &dnudB);
                        double alpha = dnudB / Bmag;
                        val_ip_C[p] = nu * dot3(ci,cj) + alpha * dot3(ci,B) * dot3(cj,B);
                    }else{
                        val_ip_C[p] = NU_LIN * dot3(ci,cj);
                    }
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, (double)(si*sj) * v);
            }
        }

        /* ========= [2] A-A : Mass (sigma_mass_A/dt) ========= */
        if(sigma_mass_A > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0,
                        (double)(si*sj) * v * inv_dt);
                }
            }
        }

        /* ========= [3] Phi-Phi : Laplace (sigma_phi * gradN·gradN) ========= */
        if(sigma_phi > 0.0){
            for(int i=0; i<fe->local_num_nodes; ++i){
                int gi = fe->conn[e][i];
                for(int j=0; j<fe->local_num_nodes; ++j){
                    int gj = fe->conn[e][j];
                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], sigma_phi);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v * inv_dt);
                }
            }
        }

        /* ========= [4] A-Phi : Coupling C ========= */
        if(sigma_cpl > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){      /* col: node (phi) */
                int gn = fe->conn[e][n];
                for(int j=0; j<ned->local_num_edges; ++j){ /* row: edge (A) */
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], ned->N_edge[e][p][j], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v * inv_dt);
                }
            }
        }

        /* ========= [5] Phi-A : Coupling C^T / dt ========= */
        if(sigma_cpl > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){      /* col: edge (A) */
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int n=0; n<fe->local_num_nodes; ++n){  /* row: node (phi) */
                    int gn = fe->conn[e][n];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], fe->geo[e][p].grad_N[n], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gn, gi, 0, 0, (double)si * v * inv_dt );
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

void apply_dirichlet_bc_for_A_and_phi_team21c(
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

            monolis_set_Dirichlet_bc_R(
                monolis,
                monolis->mat.R.B,
                ged,
                0,
                0.0
            );
        }
    }


    
    int num_nodes = fe->total_num_nodes;
    int* node_is_conductor = (int*)calloc(num_nodes, sizeof(int));
    
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop==3||prop == 4){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 3;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 1){
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

    /*
    for (int i = 0; i < num_nodes; ++i){
        if (node_is_conductor[i] == 5) {
            monolis_set_Dirichlet_bc_R(
                monolis, 
                monolis->mat.R.B, 
                i, 
                0, 
                0.0
            );
        }  
    }
    */


    int count = 0;
    if(monolis_mpi_get_global_my_rank()==0){
        for (int i = 0; i < num_nodes; ++i){
            if (node_is_conductor[i] == 1) {
                count++;
                if(count == 1){
                    monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    i, 
                    0, 
                    0.0
                );

                printf("\n\n\n\n\nadd D_bc coil1\n\n\n\n\n");

                }
                else{
                }
            }
        }
    }

    count = 0;
    if(monolis_mpi_get_global_my_rank()==2){
        for (int i = 0; i < num_nodes; ++i){
            if (node_is_conductor[i] == 3) {
                count++;
                if(count == 1){
                    monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    i, 
                    0, 
                    0.0
                );

                printf("\n\n\n\n\nadd D_bc shield\n\n\n\n\n");

                }
                else{
                }
            }
        }
    }

    count = 0;
    if(monolis_mpi_get_global_my_rank()==2){
        for (int i = 0; i < num_nodes; ++i){
            if (node_is_conductor[i] == 2) {
                count++;
                if(count == 1){
                    monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    i, 
                    0, 
                    0.0
                );

                printf("\n\n\n\n\nadd D_bc coil2\n\n\n\n\n");

                }
                else{
                }
            }
        }
    }


    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);
}

/* ============================================================
 * Coil current excitation settings for TEAM P21C-EM1
 *   prop==1 : Coil 1
 *   prop==2 : Coil 2
 *   I1 = +Iamp*sin(wt), I2 = -Iamp*sin(wt)
 * ============================================================ */
static const double I_RMS = 10.0;   /* TEAM benchmark rated current [A rms] */
static const double FREQ_HZ_team21c = 50.0; /* [Hz] */

static inline double get_coil_current_team21c(int prop, double t)
{
    const double omega = 2.0 * M_PI * FREQ_HZ_team21c;
    const double Iamp  = sqrt(2.0) * I_RMS;  /* peak value */

    if(prop == 1){
        return  Iamp * sin(omega * t);   
    } else if(prop == 2){
        return  -Iamp * sin(omega * t); 
    }  
/*
    if(prop == 1){
        return  Iamp;   
    } else if(prop == 2){
        return  -Iamp;   
    }
*/
    return 0.0;
}

static int get_team21c_simple_azimuthal_tangent(
    const COIL_INFO* coil,
    const double x_ip[3],
    double tdir[3]
){
    const double x = x_ip[0] - coil->center[0];
    const double y = x_ip[1] - coil->center[1];
    const double r2 = x*x + y*y;
    const double eps = 1.0e-20;

    if (r2 < eps) {
        tdir[0] = 0.0;
        tdir[1] = 0.0;
        tdir[2] = 0.0;
        return 0;
    }

    /* z軸まわりの周方向 e_theta = (-y/r, x/r, 0) */
    const double r = sqrt(r2);
    tdir[0] = -y / r;
    tdir[1] =  x / r;
    tdir[2] =  0.0;
    return 1;
}

static const int USE_SIMPLE_TDIR_TEST = 0;  /* 1: 簡易周方向, 0: 元の角丸矩形接線 */

/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * Voltage excitation version
 * ============================================================ */
void set_element_vec_NR_Aphi_team21c(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0 / dt;
    const double epsB   = 1.0e-14;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e = 0; e < fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_team21c(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info_team21(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation)
         *     current source amplitude is no longer prescribed by CURRENT_AMP
         *     use updated phase current g_phase_current3[phase]
         * ========================================================== */
        if(is_coil){
            double I_t = get_coil_current_team21c(prop, current_time);
            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double tdir[3];
                    double Js[3] = {0.0, 0.0, 0.0};
                    int ok_tdir = 0;

                    if (USE_SIMPLE_TDIR_TEST) {
                        ok_tdir = get_team21c_simple_azimuthal_tangent(&coil, x_ip, tdir);
                    } else {
                        ok_tdir = get_team21c_rectcoil_tangent(&coil, x_ip, tdir);
                    }

                    if (ok_tdir) {
                        Js[0] = J_mag * tdir[0];
                        Js[1] = J_mag * tdir[1];
                        Js[2] = J_mag * tdir[2];
                    }

                    val_ip_C[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis->mat.R.B[gi] += (double)si * integ;
            }
        }

        /* ==========================================================
         * [2] A-mass term
         * ========================================================== */
        if(sigma_mass_A > 0.0){
            int use_history = (prop == 1 || prop == 2 || prop == 3 || prop == 5);

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                double acc = 0.0;
                for(int j = 0; j < ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    double coeffA;
                    if(use_history){
                        coeffA = (x_curr[gj] - x_prev[gj]);
                    }else{
                        coeffA = (x_curr[gj] - x_prev[gj]);
                    }

                    for(int p = 0; p < np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            ned->N_edge[e][p][j],
                            sigma_mass_A
                        );
                    }

                    double mij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)(si * sj) * mij * coeffA * inv_dt;
                }

                monolis->mat.R.B[gi] -= acc;
            }
        }

        /* ==========================================================
         * [3] Curl-curl stiffness term
         * ========================================================== */
        for(int i = 0; i < ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int j = 0; j < ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];
                double a_val = x_curr[gj];

                for(int p = 0; p < np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    double nu_val;
                    if(nonlinear_mu){
                        double B[3];
                        compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        nu_val = get_reluctivity_nu_team21c(fmax(norm3(B), epsB));
                    }else{
                        nu_val = NU_LIN;
                    }

                    val_ip_C[p] = nu_val * dot3(ci, cj);
                }

                double kij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)(si * sj) * kij * a_val;
            }

            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [4] A-Phi coupling in A-equation
         * ========================================================== */
        for(int j = 0; j < ned->local_num_edges; ++j){
            int gj = ned->nedelec_conn[e][j];
            int sj = ned->edge_sign[e][j];

            double acc = 0.0;
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                double phi_val = x_curr[gn] - x_prev[gn];

                for(int p = 0; p < np; ++p){
                    val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[n],
                        ned->N_edge[e][p][j],
                        sigma_cpl
                    );
                }

                double cjn = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)sj * cjn * phi_val * inv_dt;
            }

            monolis->mat.R.B[gj] -= acc;
        }
    

        /* ==========================================================
         * [5] Phi-A-equation
         * ========================================================== */
        if(sigma_cpl > 0.0){
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int i = 0; i < ned->local_num_edges; ++i){
                    int gi = ned->nedelec_conn[e][i];
                    int si = ned->edge_sign[e][i];
                    double da = x_curr[gi] - x_prev[gi];

                    for(int p = 0; p < np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            fe->geo[e][p].grad_N[n],
                            sigma_cpl
                        );
                    }

                    double gin = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)si * gin * da * inv_dt;
                }

                monolis->mat.R.B[gn] -= acc;
            }
        }

        /* ==========================================================
         * [6] Phi-equation
         * ========================================================== */
        if(sigma_phi > 0.0){
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int m = 0; m < fe->local_num_nodes; ++m){
                    int gm = fe->conn[e][m];
                    double phi_m = x_curr[gm] - x_prev[gm];

                    for(int p = 0; p < np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n],
                            fe->geo[e][p].grad_N[m],
                            sigma_phi
                        );
                    }

                    double knm = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += knm * phi_m * inv_dt;
                }

                monolis->mat.R.B[gn] -= acc;
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

double calc_copper_shield_loss_EM1(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt
){
    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    double loss_local = 0.0;

    for(int e = 0; e < fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != 3) continue; /* shield only */

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double dA_dt[3]    = {0.0, 0.0, 0.0};
            double grad_phi[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double dai = (x_curr[gi] - x_prev[gi]) * inv_dt;

                dA_dt[0] += (double)si * dai * ned->N_edge[e][p][i][0];
                dA_dt[1] += (double)si * dai * ned->N_edge[e][p][i][1];
                dA_dt[2] += (double)si * dai * ned->N_edge[e][p][i][2];
            }

            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                double phi_n = x_curr[gn];

                grad_phi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                grad_phi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                grad_phi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
            }

            /* E = -dA/dt - grad(phi) */
            double E0 = -dA_dt[0] - grad_phi[0];
            double E1 = -dA_dt[1] - grad_phi[1];
            double E2 = -dA_dt[2] - grad_phi[2];

            double e2 = E0*E0 + E1*E1 + E2*E2;

            /* p = J·E = sigma |E|^2 */
            val_ip[p] = Sigma_shield * e2;
        }

        loss_local += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    /* MPI sum */
    double loss_global = loss_local;
    monolis_allreduce_R(1, &loss_global, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    return loss_global;
}

void log_copper_shield_loss_EM1(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    double shield_loss_inst
){
    if(sys->monolis_com.my_rank != 0) return;

    FILE* fp;
    fp = BBFE_sys_write_add_fopen(fp, "team21c_em1_shield_loss.csv", sys->cond.directory);

    if(step == 0){
        fprintf(fp, "Step,Time,dt,I1,I2,ShieldLossInstant\n");
    }

    fprintf(fp, "%d,%.6e,%.6e,%.6e,%.6e,%.6e\n",
            step, t, dt,
            get_coil_current_team21c(1, t),
            get_coil_current_team21c(2, t),
            shield_loss_inst);

    fclose(fp);
}


void log_copper_shield_loss_EM1_cycle_average(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    double shield_loss_inst
){
    const double T_period = 1.0 / FREQ_HZ_team21c;

    static double cycle_start_time = 0.0;
    static double sum_loss_dt = 0.0;

    sum_loss_dt += shield_loss_inst * dt;

    if((t - cycle_start_time + dt) >= T_period - 1.0e-14){
        double cycle_len = t - cycle_start_time + dt;
        if(cycle_len < 1.0e-14) cycle_len = T_period;

        double shield_loss_avg = sum_loss_dt / cycle_len;

        if(sys->monolis_com.my_rank == 0){
            FILE* fp;
            fp = BBFE_sys_write_add_fopen(
                fp, "team21c_em1_shield_loss_cycle.csv", sys->cond.directory
            );

            if(step == 0 || fabs(cycle_start_time) < 1.0e-14){
                fprintf(fp, "CycleStart,CycleEnd,CycleLength,ShieldLossAvg\n");
            }

            fprintf(
                fp,
                "%.6e,%.6e,%.6e,%.6e\n",
                cycle_start_time, t + dt, cycle_len, shield_loss_avg
            );

            fclose(fp);
        }

        cycle_start_time = t + dt;
        sum_loss_dt = 0.0;
    }
}

SHIELD_LOSS_DIAG calc_copper_shield_loss_EM1_diag(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt
){
    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    SHIELD_LOSS_DIAG out;
    memset(&out, 0, sizeof(SHIELD_LOSS_DIAG));

    for(int e = 0; e < fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != 3) continue; /* shield only */

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        double vol_elem = 0.0;
        for(int p = 0; p < np; ++p){
            val_ip[p] = 1.0;
        }
        vol_elem = BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
        out.vol_shield += vol_elem;

        for(int p = 0; p < np; ++p){
            double dA_dt[3]    = {0.0, 0.0, 0.0};
            double grad_phi[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double dai = (x_curr[gi] - x_prev[gi]) * inv_dt;

                dA_dt[0] += (double)si * dai * ned->N_edge[e][p][i][0];
                dA_dt[1] += (double)si * dai * ned->N_edge[e][p][i][1];
                dA_dt[2] += (double)si * dai * ned->N_edge[e][p][i][2];
            }

            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                double phi_n = (x_curr[gn] - x_prev[gn]) * inv_dt;

                grad_phi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                grad_phi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                grad_phi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
            }

            /* E = -dA/dt - grad(phi) */
            double E0 = -dA_dt[0] - grad_phi[0];
            double E1 = -dA_dt[1] - grad_phi[1];
            double E2 = -dA_dt[2] - grad_phi[2];

            double a2 = dA_dt[0]*dA_dt[0] + dA_dt[1]*dA_dt[1] + dA_dt[2]*dA_dt[2];
            double p2 = grad_phi[0]*grad_phi[0] + grad_phi[1]*grad_phi[1] + grad_phi[2]*grad_phi[2];
            double ep = dA_dt[0]*grad_phi[0] + dA_dt[1]*grad_phi[1] + dA_dt[2]*grad_phi[2];
            double e2 = E0*E0 + E1*E1 + E2*E2;

            double w = basis->integ_weight[p] * Jacobian_ip[p];

            out.loss_A     += Sigma_shield * a2 * w;
            out.loss_phi   += Sigma_shield * p2 * w;
            out.loss_cross += Sigma_shield * 2.0 * ep * w;
            out.loss_total += Sigma_shield * e2 * w;

            out.rms_dA_dt    += a2 * w;
            out.rms_grad_phi += p2 * w;
            out.rms_E        += e2 * w;

            double na = sqrt(a2);
            double npg = sqrt(p2);
            double ne = sqrt(e2);

            if(na  > out.max_dA_dt)    out.max_dA_dt = na;
            if(npg > out.max_grad_phi) out.max_grad_phi = npg;
            if(ne  > out.max_E)        out.max_E = ne;
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    /* MPI sum for additive quantities */
    {
        double sendbuf[8], recvbuf[8];
        sendbuf[0] = out.loss_total;
        sendbuf[1] = out.loss_A;
        sendbuf[2] = out.loss_phi;
        sendbuf[3] = out.loss_cross;
        sendbuf[4] = out.vol_shield;
        sendbuf[5] = out.rms_dA_dt;
        sendbuf[6] = out.rms_grad_phi;
        sendbuf[7] = out.rms_E;

        memcpy(recvbuf, sendbuf, sizeof(sendbuf));
        monolis_allreduce_R(8, recvbuf, MONOLIS_MPI_SUM, sys->monolis_com.comm);

        out.loss_total   = recvbuf[0];
        out.loss_A       = recvbuf[1];
        out.loss_phi     = recvbuf[2];
        out.loss_cross   = recvbuf[3];
        out.vol_shield   = recvbuf[4];
        out.rms_dA_dt    = recvbuf[5];
        out.rms_grad_phi = recvbuf[6];
        out.rms_E        = recvbuf[7];
    }

    /* MPI max for maxima */
    {
        double sendbuf[3], recvbuf[3];
        sendbuf[0] = out.max_dA_dt;
        sendbuf[1] = out.max_grad_phi;
        sendbuf[2] = out.max_E;

        memcpy(recvbuf, sendbuf, sizeof(sendbuf));
        monolis_allreduce_R(3, recvbuf, MONOLIS_MPI_MAX, sys->monolis_com.comm);

        out.max_dA_dt    = recvbuf[0];
        out.max_grad_phi = recvbuf[1];
        out.max_E        = recvbuf[2];
    }

    if(out.vol_shield > 1.0e-30){
        out.rms_dA_dt    = sqrt(out.rms_dA_dt    / out.vol_shield);
        out.rms_grad_phi = sqrt(out.rms_grad_phi / out.vol_shield);
        out.rms_E        = sqrt(out.rms_E        / out.vol_shield);
    } else {
        out.rms_dA_dt = 0.0;
        out.rms_grad_phi = 0.0;
        out.rms_E = 0.0;
    }

    return out;
}


void log_copper_shield_loss_EM1_diag(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    const SHIELD_LOSS_DIAG* d
){
    if(sys->monolis_com.my_rank != 0) return;

    FILE* fp;
    fp = BBFE_sys_write_add_fopen(fp, "team21c_em1_shield_loss_diag.csv", sys->cond.directory);

    if(step == 0){
        fprintf(fp,
            "Step,Time,dt,I1,I2,"
            "LossTotal,LossA,LossPhi,LossCross,"
            "VolShield,RMS_dA_dt,RMS_gradPhi,RMS_E,"
            "Max_dA_dt,Max_gradPhi,Max_E\n");
    }

    fprintf(fp,
        "%d,%.6e,%.6e,%.6e,%.6e,"
        "%.6e,%.6e,%.6e,%.6e,"
        "%.6e,%.6e,%.6e,%.6e,"
        "%.6e,%.6e,%.6e\n",
        step, t, dt,
        get_coil_current_team21c(1, t),
        get_coil_current_team21c(2, t),
        d->loss_total, d->loss_A, d->loss_phi, d->loss_cross,
        d->vol_shield, d->rms_dA_dt, d->rms_grad_phi, d->rms_E,
        d->max_dA_dt, d->max_grad_phi, d->max_E
    );

    fclose(fp);
}

typedef struct {
    double vol;              /* coil volume */
    double int_J_dot_t_dV;   /* ∫_Omega J·t dV */
    double NI_equiv;         /* (1/L) ∫_Omega J·t dV */
    double NI_target;        /* turns * I(t) */
    double rel_err;          /* (NI_equiv - NI_target)/NI_target */
} COIL_AMPERE_TURN_DIAG;

COIL_AMPERE_TURN_DIAG calc_coil_ampere_turn_diag_team21c(
    FE_SYSTEM* sys,
    int coil_prop,      /* 1 or 2 */
    double current_time
){
    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    COIL_AMPERE_TURN_DIAG out;
    memset(&out, 0, sizeof(COIL_AMPERE_TURN_DIAG));

    COIL_INFO coil;
    if(!get_coil_info_team21(coil_prop, &coil)){
        BB_std_free_1d_double(Jacobian_ip, np);
        BB_std_free_1d_double(val_ip, np);
        return out;
    }

    const double I_t   = get_coil_current_team21c(coil_prop, current_time);
    const double J_mag = (coil.turns * I_t) / coil.area;

    /* gmsh geometry: coil is extruded along z by 217 mm */
    const double coil_length = 217.0 * MM_TO_M;

    double vol_local = 0.0;
    double int_local = 0.0;

    for(int e = 0; e < fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != coil_prop) continue;

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* element volume */
        for(int p = 0; p < np; ++p){
            val_ip[p] = 1.0;
        }
        vol_local += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);

        /* ∫ J·t dV */
        for(int p = 0; p < np; ++p){
            double x_ip[3];
            double tdir[3] = {0.0, 0.0, 0.0};
            double Js[3]   = {0.0, 0.0, 0.0};

            get_interp_coords(e, p, fe, basis, x_ip);

            if(get_team21c_rectcoil_tangent(&coil, x_ip, tdir)){
                Js[0] = J_mag * tdir[0];
                Js[1] = J_mag * tdir[1];
                Js[2] = J_mag * tdir[2];
            }

            val_ip[p] = Js[0]*tdir[0] + Js[1]*tdir[1] + Js[2]*tdir[2];
        }
        int_local += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    /* MPI sum */
    {
        double sendbuf[2], recvbuf[2];
        sendbuf[0] = vol_local;
        sendbuf[1] = int_local;
        memcpy(recvbuf, sendbuf, sizeof(sendbuf));
        monolis_allreduce_R(2, recvbuf, MONOLIS_MPI_SUM, sys->monolis_com.comm);

        out.vol            = recvbuf[0];
        out.int_J_dot_t_dV = recvbuf[1];
    }

    out.NI_target = coil.turns * I_t;
    out.NI_equiv  = out.int_J_dot_t_dV / coil_length;

    if(fabs(out.NI_target) > 1.0e-30){
        out.rel_err = (out.NI_equiv - out.NI_target) / out.NI_target;
    } else {
        out.rel_err = 0.0;
    }

    return out;
}

void log_coil_ampere_turn_diag_team21c(
    FE_SYSTEM* sys,
    int step,
    double current_time
){
    COIL_AMPERE_TURN_DIAG c1 =
        calc_coil_ampere_turn_diag_team21c(sys, 1, current_time);
    COIL_AMPERE_TURN_DIAG c2 =
        calc_coil_ampere_turn_diag_team21c(sys, 2, current_time);

    if(sys->monolis_com.my_rank != 0) return;

    FILE* fp;
    fp = BBFE_sys_write_add_fopen(fp, "team21c_coil_ampere_turn_diag.csv", sys->cond.directory);

    if(step == 0){
        fprintf(fp,
            "Step,Time,"
            "NI_target_1,NI_equiv_1,RelErr_1,Vol_1,IntJt_1,"
            "NI_target_2,NI_equiv_2,RelErr_2,Vol_2,IntJt_2\n");
    }

    fprintf(fp,
        "%d,%.6e,"
        "%.6e,%.6e,%.6e,%.6e,%.6e,"
        "%.6e,%.6e,%.6e,%.6e,%.6e\n",
        step, current_time,
        c1.NI_target, c1.NI_equiv, c1.rel_err, c1.vol, c1.int_J_dot_t_dV,
        c2.NI_target, c2.NI_equiv, c2.rel_err, c2.vol, c2.int_J_dot_t_dV
    );

    fclose(fp);
}

SHIELD_FIELD_INT_DIAG calc_shield_field_integrals_EM1(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt,
    int use_phi_in_shield   /* 1: grad phi を含める, 0: grad phi = 0 */
){
    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    SHIELD_FIELD_INT_DIAG out;
    memset(&out, 0, sizeof(SHIELD_FIELD_INT_DIAG));

    for(int e = 0; e < fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != 3) continue; /* shield only */

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* volume */
        for(int p = 0; p < np; ++p){
            val_ip[p] = 1.0;
        }
        out.vol_shield += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);

        /* ∫ |dA/dt|^2 dV */
        for(int p = 0; p < np; ++p){
            double dA_dt[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double dai = (x_curr[gi] - x_prev[gi]) * inv_dt;

                dA_dt[0] += (double)si * dai * ned->N_edge[e][p][i][0];
                dA_dt[1] += (double)si * dai * ned->N_edge[e][p][i][1];
                dA_dt[2] += (double)si * dai * ned->N_edge[e][p][i][2];
            }

            val_ip[p] = dA_dt[0]*dA_dt[0] + dA_dt[1]*dA_dt[1] + dA_dt[2]*dA_dt[2];
        }
        out.int_dA2 += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);

        /* ∫ |grad phi|^2 dV */
        for(int p = 0; p < np; ++p){
            double grad_phi[3] = {0.0, 0.0, 0.0};

            if(use_phi_in_shield){
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_n = x_curr[gn];

                    grad_phi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                    grad_phi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                    grad_phi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
                }
            }

            val_ip[p] = grad_phi[0]*grad_phi[0]
                      + grad_phi[1]*grad_phi[1]
                      + grad_phi[2]*grad_phi[2];
        }
        out.int_gradphi2 += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);

        /* ∫ (dA/dt · grad phi) dV */
        for(int p = 0; p < np; ++p){
            double dA_dt[3]    = {0.0, 0.0, 0.0};
            double grad_phi[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double dai = (x_curr[gi] - x_prev[gi]) * inv_dt;

                dA_dt[0] += (double)si * dai * ned->N_edge[e][p][i][0];
                dA_dt[1] += (double)si * dai * ned->N_edge[e][p][i][1];
                dA_dt[2] += (double)si * dai * ned->N_edge[e][p][i][2];
            }

            if(use_phi_in_shield){
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_n = x_curr[gn];

                    grad_phi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                    grad_phi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                    grad_phi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
                }
            }

            val_ip[p] = dA_dt[0]*grad_phi[0]
                      + dA_dt[1]*grad_phi[1]
                      + dA_dt[2]*grad_phi[2];
        }
        out.int_cross += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    /* MPI sum */
    {
        double sendbuf[4], recvbuf[4];
        sendbuf[0] = out.int_dA2;
        sendbuf[1] = out.int_gradphi2;
        sendbuf[2] = out.int_cross;
        sendbuf[3] = out.vol_shield;

        memcpy(recvbuf, sendbuf, sizeof(sendbuf));
        monolis_allreduce_R(4, recvbuf, MONOLIS_MPI_SUM, sys->monolis_com.comm);

        out.int_dA2      = recvbuf[0];
        out.int_gradphi2 = recvbuf[1];
        out.int_cross    = recvbuf[2];
        out.vol_shield   = recvbuf[3];
    }

    /* correlation coefficient */
    {
        const double denom = sqrt(out.int_dA2 * out.int_gradphi2);
        if(denom > 1.0e-30){
            out.rho = out.int_cross / denom;
        } else {
            out.rho = 0.0;
        }
    }

    return out;
}

void log_shield_field_integrals_EM1(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    const SHIELD_FIELD_INT_DIAG* d
){
    if(sys->monolis_com.my_rank != 0) return;

    FILE* fp;
    fp = BBFE_sys_write_add_fopen(fp, "team21c_em1_shield_field_integrals.csv", sys->cond.directory);

    if(step == 0){
        fprintf(fp,
            "Step,Time,dt,I1,I2,"
            "Int_dA2,Int_gradPhi2,Int_cross,Rho,VolShield\n");
    }

    fprintf(fp,
        "%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
        step, t, dt,
        get_coil_current_team21c(1, t),
        get_coil_current_team21c(2, t),
        d->int_dA2,
        d->int_gradphi2,
        d->int_cross,
        d->rho,
        d->vol_shield
    );

    fclose(fp);
}

const double mu0   = 4.0*M_PI*1e-7;      // H/m
const double Nu    = 1.0 / mu0;          // 1/H/m = A/(T·m)
const double freq  = 50.0;               // Hz
const double Omega = 2.0*M_PI*freq;      // rad/s  ★これが正しい


void set_element_mat_nedelec_Aphi_team21a2(
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
        get_sigmas_for_prop_team21c(prop, &sigma_massA, &sigma_cpl, &sigma_phi);

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


static const double I_RMS_team21a2 = 10.0;   /* TEAM benchmark rated current [A rms] */
static const double FREQ_HZ_team21a2 = 50.0; /* [Hz] */

static inline double get_coil_current_team21a2(int prop, double t)
{
    const double omega = 2.0 * M_PI * FREQ_HZ_team21a2;
    const double Iamp  = sqrt(2.0) * I_RMS_team21a2;  /* peak value */    

    if(prop == 1){
        return  Iamp;   
    } else if(prop == 2){
        return  -Iamp;   
    }

    return 0.0;
}


void set_element_vec_nedelec_Aphi_team21a2(
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
        get_sigmas_for_prop_team21c(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info_team21(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation)
         *     current source amplitude is no longer prescribed by CURRENT_AMP
         *     use updated phase current g_phase_current3[phase]
         * ========================================================== */
        if(is_coil){
            double I_t = get_coil_current_team21a2(prop, 0.0);
            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double tdir[3];
                    double Js[3] = {0.0, 0.0, 0.0};
                    int ok_tdir = 0;

                    if (USE_SIMPLE_TDIR_TEST) {
                        ok_tdir = get_team21c_simple_azimuthal_tangent(&coil, x_ip, tdir);
                    } else {
                        ok_tdir = get_team21c_rectcoil_tangent(&coil, x_ip, tdir);
                    }

                    if (ok_tdir) {
                        Js[0] = J_mag * tdir[0];
                        Js[1] = J_mag * tdir[1];
                        Js[2] = J_mag * tdir[2];
                    }

                    val_ip_C[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                double _Complex val = (double)si * integ;
                monolis->mat.C.B[gi] -= val;
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}


void apply_dirichlet_bc_for_A_and_phi_team21a2(
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

        if(prop == 4){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 4;
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
/*
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 4){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 4; 
            }
        }
    }
*/
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 3){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 3; 
            }
        }
    }


    for (int i = 0; i < num_nodes; ++i){
        if (node_is_conductor[i] == 1||node_is_conductor[i] == 2||node_is_conductor[i] == 4) {
        //if (node_is_conductor[i] == 1||node_is_conductor[i] == 2||node_is_conductor[i] == 4 ||node_is_conductor[i] == 3) {
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

double calc_copper_shield_loss_EM1_freq(
    FE_SYSTEM* sys,
    const double _Complex* x_c,
    double freq_hz,
    const char* directory
){
    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;
    const double omega = 2.0 * M_PI * freq_hz;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    double loss_local = 0.0;
    int total_num_elems = 0;


    FILE* fp;
    int BUFFER_SIZE = 1024;
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    char filename[BUFFER_SIZE];
    int graph_ndof;
    int tmp;

    snprintf(fname_n_internal_graph, BUFFER_SIZE, "parted.0/graph_elem.dat.n_internal.%d", monolis_mpi_get_global_my_rank());
    fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
    fscanf(fp, "%s %d", char_n_internal, &(tmp));
    fscanf(fp, "%d", &(total_num_elems));
    fclose(fp);

    /*
    snprintf(fname_n_internal_graph, BUFFER_SIZE, "parted.0/elem.dat.%d", filename, monolis_mpi_get_global_my_rank());
    fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
    // read the num of elements
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &(total_num_elems));
	*/

    for(int e = 0; e < total_num_elems; ++e){
        if(ned->elem_prop[e] != 3) continue; /* shield only */

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double _Complex A[3] = {0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
            double _Complex grad_phi[3] = {0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double _Complex ai = x_c[gi];

                A[0] += (double)si * ai * ned->N_edge[e][p][i][0];
                A[1] += (double)si * ai * ned->N_edge[e][p][i][1];
                A[2] += (double)si * ai * ned->N_edge[e][p][i][2];
            }

            /* phi を使う定式化ならここを有効化 */
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                double _Complex phi_n = x_c[gn];

                grad_phi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                grad_phi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                grad_phi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
            }

            /* E = -jωA - grad(phi) */
            //double _Complex E0 = -(I*omega)*A[0] - grad_phi[0];
            //double _Complex E1 = -(I*omega)*A[1] - grad_phi[1];
            //double _Complex E2 = -(I*omega)*A[2] - grad_phi[2];

            //double _Complex E0 = -A[0] - grad_phi[0];
            //double _Complex E1 = -A[1] - grad_phi[1];
            //double _Complex E2 = -A[2] - grad_phi[2];

	    double _Complex E0 = -(I*omega)*A[0] - (I*omega)*grad_phi[0];
            double _Complex E1 = -(I*omega)*A[1] - (I*omega)*grad_phi[1];
            double _Complex E2 = -(I*omega)*A[2] - (I*omega)*grad_phi[2];



            double e2 =
                creal(E0*conj(E0) + E1*conj(E1) + E2*conj(E2));

            /* peak phasor 前提なら 0.5 を掛ける */
            val_ip[p] = 0.5 * Sigma_steel * e2;
        }

        loss_local += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    double loss_global = loss_local;
    monolis_allreduce_R(1, &loss_global, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    return loss_global;
}




/* ============================================================
 * 3x3 tensor utilities for anisotropic sigma_eff
 * ============================================================ */
static inline void mat3_zero(double A[3][3])
{
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            A[i][j] = 0.0;
        }
    }
}

static inline void mat3_copy(double A[3][3], const double B[3][3])
{
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            A[i][j] = B[i][j];
        }
    }
}

static inline double dot3_local(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double norm3_local(const double a[3])
{
    return sqrt(dot3_local(a,a));
}

static inline void cross3_local(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void normalize3_local(double a[3])
{
    double n = norm3_local(a);
    if(n > 1.0e-30){
        a[0] /= n; a[1] /= n; a[2] /= n;
    }
}

static inline void matvec3_local(const double A[3][3], const double x[3], double y[3])
{
    for(int i=0;i<3;++i){
        y[i] = A[i][0]*x[0] + A[i][1]*x[1] + A[i][2]*x[2];
    }
}

static inline double bilinear3_local(
    const double a[3],
    const double S[3][3],
    const double b[3]
){
    double Sb[3];
    matvec3_local(S, b, Sb);
    return dot3_local(a, Sb);
}

/* S += alpha * u \otimes u */
static inline void add_outer_scaled_3x3(double S[3][3], double alpha, const double u[3])
{
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            S[i][j] += alpha * u[i] * u[j];
        }
    }
}

/* t から直交基底 (t, n1, n2) を作る */
static void build_orthonormal_frame_from_t(
    const double t_in[3],
    double t[3],
    double n1[3],
    double n2[3]
){
    t[0] = t_in[0];
    t[1] = t_in[1];
    t[2] = t_in[2];
    normalize3_local(t);

    /* t と平行でない適当な参照ベクトルを選ぶ */
    double ref[3] = {0.0, 0.0, 1.0};
    if(fabs(t[2]) > 0.9){
        ref[0] = 1.0; ref[1] = 0.0; ref[2] = 0.0;
    }

    cross3_local(ref, t, n1);
    normalize3_local(n1);

    cross3_local(t, n1, n2);
    normalize3_local(n2);
}

/* isotropic tensor: s * I */
static inline void set_isotropic_tensor(double S[3][3], double s)
{
    mat3_zero(S);
    S[0][0] = s;
    S[1][1] = s;
    S[2][2] = s;
}

/* ============================================================
 * sigma_eff tensor model
 *
 * coil:
 *   sigma_eff = sigma_parallel(w) * (t \otimes t)
 *             + sigma_n1(w)       * (n1 \otimes n1)
 *             + sigma_n2(w)       * (n2 \otimes n2)
 *
 * shield:
 *   isotropic sigma
 * ============================================================ */
static void get_sigma_tensors_for_prop_team21c(
    int prop,
    const double x_ip[3],
    double sigma_mass_A[3][3],
    double sigma_cpl[3][3],
    double sigma_phi[3][3]
){
    mat3_zero(sigma_mass_A);
    mat3_zero(sigma_cpl);
    mat3_zero(sigma_phi);

    if(prop == 1 || prop == 2){
        /* ---- coil region: anisotropic effective tensor ---- */
        COIL_INFO coil;
        if(!get_coil_info_team21(prop, &coil)){
            return;
        }

        double tdir[3];
        if(!get_team21c_rectcoil_tangent(&coil, x_ip, tdir)){
            return;
        }

        double t[3], n1[3], n2[3];
        build_orthonormal_frame_from_t(tdir, t, n1, n2);

        /* -----------------------------------------------
         * ここを必要に応じて周波数依存に差し替える
         * まずは実数・定数版
         *
         * 例:
         *   sigma_parallel = fill * Sigma_coil;
         *   sigma_n1 = 0;
         *   sigma_n2 = 0;
         *
         * 元コードの縮約値 Sigma_coil * 0.1256/3.2 を
         * とりあえず巻線方向成分として使うなら下記
         * ----------------------------------------------- */
        const double sigma_parallel = Sigma_coil * 0.1256;
        const double sigma_n1_val   = 0.0;
        const double sigma_n2_val   = 0.0;

        add_outer_scaled_3x3(sigma_mass_A, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_mass_A, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_mass_A, sigma_n2_val,   n2);

        /* A-phi, phi-A, phi-phi も同じ tensor を使うなら */
        add_outer_scaled_3x3(sigma_cpl, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_cpl, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_cpl, sigma_n2_val,   n2);

        add_outer_scaled_3x3(sigma_phi, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_phi, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_phi, sigma_n2_val,   n2);

        /* もし従来どおり coil では cpl/phi を切りたいなら:
         * mat3_zero(sigma_cpl);
         * mat3_zero(sigma_phi);
         */
    }
    else if(prop == 3){
        /* shield: isotropic */
        set_isotropic_tensor(sigma_mass_A, Sigma_steel);
        set_isotropic_tensor(sigma_cpl,    Sigma_steel);
        set_isotropic_tensor(sigma_phi,    Sigma_steel);
    }
    else{
        /* air / steel core / other => zero for conductivity part */
        mat3_zero(sigma_mass_A);
        mat3_zero(sigma_cpl);
        mat3_zero(sigma_phi);
    }
}

void set_element_mat_NR_Aphi_team21c_hom(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0/dt;
    const double epsB   = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ========= [1] A-A : Curl-Curl (tangent stiffness) ========= */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    if(nonlinear_mu){
                        double B[3]; compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        double Bmag = fmax(norm3(B), epsB);
                        double nu, dnudB;
                        eval_nu_and_dnudB_team21c(Bmag, &nu, &dnudB);
                        double alpha = dnudB / Bmag;
                        val_ip_C[p] = nu * dot3(ci,cj) + alpha * dot3(ci,B) * dot3(cj,B);
                    }else{
                        val_ip_C[p] = NU_LIN * dot3(ci,cj);
                    }
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, (double)(si*sj) * v);
            }
        }

        /* ========= [2] A-A : Mass ((sigma_eff/dt) N_i·N_j) ========= */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        ned->N_edge[e][p][i],
                        sigma_mass_A,
                        ned->N_edge[e][p][j]
                    );
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                if(fabs(v) > 0.0){
                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis, gi, gj, 0, 0, (double)(si*sj) * v * inv_dt);
                }
            }
        }

        /* ========= [3] Phi-Phi : gradN_i · sigma_eff · gradN_j ========= */
        for(int i=0; i<fe->local_num_nodes; ++i){
            int gi = fe->conn[e][i];

            for(int j=0; j<fe->local_num_nodes; ++j){
                int gj = fe->conn[e][j];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        fe->geo[e][p].grad_N[i],
                        sigma_phi_t,
                        fe->geo[e][p].grad_N[j]
                    );
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                if(fabs(v) > 0.0){
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v);
                }
            }
        }

        /* ========= [4] A-Phi : (N_edge_j) · sigma_eff · gradN_n ========= */
        for(int n=0; n<fe->local_num_nodes; ++n){      /* col: phi */
            int gn = fe->conn[e][n];

            for(int j=0; j<ned->local_num_edges; ++j){ /* row: A */
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        ned->N_edge[e][p][j],
                        sigma_cpl_t,
                        fe->geo[e][p].grad_N[n]
                    );
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                if(fabs(v) > 0.0){
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v);
                }
            }
        }

        /* ========= [5] Phi-A : gradN_n · sigma_eff · N_edge_i / dt ========= */
        for(int i=0; i<ned->local_num_edges; ++i){      /* col: A */
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int n=0; n<fe->local_num_nodes; ++n){   /* row: phi */
                int gn = fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        fe->geo[e][p].grad_N[n],
                        sigma_cpl_t,
                        ned->N_edge[e][p][i]
                    );
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                if(fabs(v) > 0.0){
                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis, gn, gi, 0, 0, (double)si * v * inv_dt);
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

void set_element_vec_NR_Aphi_team21c_hom(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0 / dt;
    const double epsB   = 1.0e-14;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e = 0; e < fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info_team21(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation)
         *     ここは元コードをそのまま維持
         * ========================================================== */
        if(is_coil){
            double I_t = get_coil_current_team21c(prop, current_time);
            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double tdir[3];
                    double Js[3] = {0.0, 0.0, 0.0};
                    int ok_tdir = 0;

                    if (USE_SIMPLE_TDIR_TEST) {
                        ok_tdir = get_team21c_simple_azimuthal_tangent(&coil, x_ip, tdir);
                    } else {
                        ok_tdir = get_team21c_rectcoil_tangent(&coil, x_ip, tdir);
                    }

                    if (ok_tdir) {
                        Js[0] = J_mag * tdir[0];
                        Js[1] = J_mag * tdir[1];
                        Js[2] = J_mag * tdir[2];
                    }

                    val_ip_C[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                /* residual B = -F(x) */
                monolis->mat.R.B[gi] += (double)si * integ;
            }
        }

        /* ==========================================================
         * [2] A-equation : Curl-Curl term
         *     ここも元コードのまま
         * ========================================================== */
        for(int i = 0; i < ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int j = 0; j < ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];
                double aj = x_curr[gj];

                for(int p = 0; p < np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    if(nonlinear_mu){
                        double B[3];
                        compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        double Bmag = fmax(norm3(B), epsB);

                        double nu, dnudB;
                        eval_nu_and_dnudB_team21c(Bmag, &nu, &dnudB);
                        double alpha = dnudB / Bmag;

                        val_ip_C[p] = nu * dot3(ci, cj) + alpha * dot3(ci, B) * dot3(cj, B);
                    }else{
                        val_ip_C[p] = NU_LIN * dot3(ci, cj);
                    }
                }

                double kij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)(si * sj) * kij * aj;
            }

            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [3] A-A mass term : N_i · sigma_eff · (A^{n+1}-A^n)/dt
         * ========================================================== */
        for(int i = 0; i < ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int j = 0; j < ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];
                double da = x_curr[gj] - x_prev[gj];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        ned->N_edge[e][p][i],
                        sigma_mass_A,
                        ned->N_edge[e][p][j]
                    );
                }

                double mij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)(si * sj) * mij * da * inv_dt;
            }

            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [4] A-Phi-equation
         *     A-row residual: \int N_i · sigma_eff · grad(phi) dV
         * ========================================================== */
        for(int i = 0; i < ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                double phi_val = x_curr[gn];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        ned->N_edge[e][p][i],
                        sigma_cpl_t,
                        fe->geo[e][p].grad_N[n]
                    );
                }

                double cin = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)si * cin * phi_val;
            }

            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [5] Phi-A-equation
         *     phi-row residual: \int grad(N_n) · sigma_eff · (A^{n+1}-A^n)/dt dV
         * ========================================================== */
        for(int n = 0; n < fe->local_num_nodes; ++n){
            int gn = fe->conn[e][n];

            double acc = 0.0;
            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double da = x_curr[gi] - x_prev[gi];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        fe->geo[e][p].grad_N[n],
                        sigma_cpl_t,
                        ned->N_edge[e][p][i]
                    );
                }

                double gni = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)si * gni * da * inv_dt;
            }

            monolis->mat.R.B[gn] -= acc;
        }

        /* ==========================================================
         * [6] Phi-equation
         *     phi-row residual: \int grad(N_n) · sigma_eff · grad(phi) dV
         * ========================================================== */
        for(int n = 0; n < fe->local_num_nodes; ++n){
            int gn = fe->conn[e][n];

            double acc = 0.0;
            for(int m = 0; m < fe->local_num_nodes; ++m){
                int gm = fe->conn[e][m];
                double phi_m = x_curr[gm];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl_t[3][3], sigma_phi_t[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21c(
                        prop, x_ip, sigma_mass_A, sigma_cpl_t, sigma_phi_t);

                    val_ip_C[p] = bilinear3_local(
                        fe->geo[e][p].grad_N[n],
                        sigma_phi_t,
                        fe->geo[e][p].grad_N[m]
                    );
                }

                double knm = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += knm * phi_m;
            }

            monolis->mat.R.B[gn] -= acc;
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}


static void get_sigma_tensors_for_prop_team21a2(
    int prop,
    const double x_ip[3],
    double sigma_mass_A[3][3],
    double sigma_cpl[3][3],
    double sigma_phi[3][3]
){
    mat3_zero(sigma_mass_A);
    mat3_zero(sigma_cpl);
    mat3_zero(sigma_phi);

    if(prop == 1 || prop == 2){
        COIL_INFO coil;
        if(!get_coil_info_team21(prop, &coil)) return;

        double tdir[3];
        int ok_tdir = 0;

        if (USE_SIMPLE_TDIR_TEST) {
            ok_tdir = get_team21c_simple_azimuthal_tangent(&coil, x_ip, tdir);
        } else {
            ok_tdir = get_team21c_rectcoil_tangent(&coil, x_ip, tdir);
        }
        if(!ok_tdir) return;

        double t[3], n1[3], n2[3];
        build_orthonormal_frame_from_t(tdir, t, n1, n2);

        /* まずは巻線方向のみ */
        const double sigma_parallel = Sigma_coil * 0.1256;
        const double sigma_n1_val   = 0.0;
        const double sigma_n2_val   = 0.0;

        add_outer_scaled_3x3(sigma_mass_A, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_mass_A, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_mass_A, sigma_n2_val,   n2);

        /* current-driven coil: keep phi coupling OFF */
        add_outer_scaled_3x3(sigma_cpl, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_cpl, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_cpl, sigma_n2_val,   n2);

        add_outer_scaled_3x3(sigma_phi, sigma_parallel, t);
        add_outer_scaled_3x3(sigma_phi, sigma_n1_val,   n1);
        add_outer_scaled_3x3(sigma_phi, sigma_n2_val,   n2);

        mat3_zero(sigma_cpl);
        mat3_zero(sigma_phi);
    }
    else if(prop == 3){
        set_isotropic_tensor(sigma_mass_A, Sigma_steel);
        set_isotropic_tensor(sigma_cpl,    Sigma_steel);
        set_isotropic_tensor(sigma_phi,    Sigma_steel);
    }
    else{
        mat3_zero(sigma_mass_A);
        mat3_zero(sigma_cpl);
        mat3_zero(sigma_phi);
    }
}

void set_element_mat_nedelec_Aphi_team21a2_hom(
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
        int nonlinear_mu = 0;

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        /* [1] Curl-Curl */
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
                v *= (double)(si*sj);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        /* [2] A-A : jω N_i · sigma_mass_A · N_j */
        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j=0; j<ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl[3][3], sigma_phi[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21a2(
                        prop, x_ip, sigma_mass_A, sigma_cpl, sigma_phi);

                    val_ip_C[p] =
                        bilinear3_local(
                            ned->N_edge[e][p][i],
                            sigma_mass_A,
                            ned->N_edge[e][p][j]
                        ) * I * Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)(si*sj);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        /* [3] phi-phi : jω gradN_m · sigma_phi · gradN_n */
        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl[3][3], sigma_phi[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21a2(
                        prop, x_ip, sigma_mass_A, sigma_cpl, sigma_phi);

                    val_ip_C[p] =
                        bilinear3_local(
                            fe->geo[e][p].grad_N[m],
                            sigma_phi,
                            fe->geo[e][p].grad_N[n]
                        ) * I * Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
            }
        }

        /* [4] A-phi : jω N_j · sigma_cpl · gradN_n */
        for(int j=0; j<ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl[3][3], sigma_phi[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21a2(
                        prop, x_ip, sigma_mass_A, sigma_cpl, sigma_phi);

                    val_ip_C[p] =
                        bilinear3_local(
                            ned->N_edge[e][p][j],
                            sigma_cpl,
                            fe->geo[e][p].grad_N[n]
                        ) * I * Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)sj;
                monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
            }
        }

        /* [5] phi-A : jω gradN_m · sigma_cpl · N_i */
        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    double sigma_mass_A[3][3], sigma_cpl[3][3], sigma_phi[3][3];

                    get_interp_coords(e, p, fe, basis, x_ip);
                    get_sigma_tensors_for_prop_team21a2(
                        prop, x_ip, sigma_mass_A, sigma_cpl, sigma_phi);

                    val_ip_C[p] =
                        bilinear3_local(
                            fe->geo[e][p].grad_N[m],
                            sigma_cpl,
                            ned->N_edge[e][p][i]
                        ) * I * Omega;
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