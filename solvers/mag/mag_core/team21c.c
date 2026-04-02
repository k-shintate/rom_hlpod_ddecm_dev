#include "team21c.h"

/* Air stabilization factors (tune as needed) */
static const double AIR_PHI_SIGMA_FACTOR  = 0.0; /* keep phi disabled in air for P21C-EM1 */
static const double AIR_A_MASS_SIGMA_FACTOR = 0.0; /* keep A-mass disabled in air for P21C-EM1 */

const double Sigma_coil   = 5.7143e7;
const double Sigma_shield = 5.7143e7;
const double Sigma_steel  = 6.484e6;

void get_sigmas_for_prop_team21c2(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2){
        /* exciting coil conductor */
        *sigma_mass_A = Sigma_coil;
        *sigma_cpl    = Sigma_coil;
        *sigma_phi    = Sigma_coil;
    } else if(prop == 3){
        /* TEAM P21C-EM1 copper shielding plate */
        *sigma_mass_A = Sigma_shield;
        *sigma_cpl    = Sigma_shield;
        *sigma_phi    = Sigma_shield;
    } else if(prop == 4){
        *sigma_mass_A = Sigma_steel;
        *sigma_cpl    = Sigma_steel;
        *sigma_phi    = Sigma_steel;
    } else if(prop == 5){
        /* air */
        *sigma_mass_A = Sigma_steel*1.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = Sigma_steel*1.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}

void get_sigmas_for_prop_team21c(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2){
        /* exciting coil conductor */
        *sigma_mass_A = Sigma_coil;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 3){
        /* TEAM P21C-EM1 copper shielding plate */
        *sigma_mass_A = Sigma_shield;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 4){
        *sigma_mass_A = Sigma_steel;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 5){
        /* air */
        *sigma_mass_A = Sigma_steel*1.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = Sigma_steel*1.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


static inline int get_coil_info_team21(int elem_prop, COIL_INFO* info) {
    info->axis[0] = 0.0;
    info->axis[1] = 1.0;
    info->axis[2] = 0.0;
    info->turns   = 300.0;
    info->area    = 13.04 * (MM_TO_M * MM_TO_M);

    {
        const double cy = 30.0 * MM_TO_M;
        const double cz = 20.0 * MM_TO_M;

        if (elem_prop == 1) {
            info->center[0] = 2.5 * MM_TO_M;
            info->center[1] = cy;
            info->center[2] = cz;
            return 1;
        }
        if (elem_prop == 2) {
            info->center[0] = 35.0 * MM_TO_M;
            info->center[1] = cy;
            info->center[2] = cz;
            return 1;
        }
    }
    return 0;
}

/* --- Nonlinear Material Property (Brauer Law) --- */
/* nu(B) = 100 + 10 * exp(2 * |B|^2) */
void eval_nu_and_dnudB_team21c(double Bmag, double* nu, double* dnudB)
{
    const double B2_MAX = 6.25;   /* clamp: B <= 2.5 T 相当 */
    double B2 = Bmag * Bmag;

    if (B2 > B2_MAX) {
        B2 = B2_MAX;
    }

    const double exp_val = exp(2.0 * B2);

    *nu    = 100.0 + 10.0 * exp_val;
    *dnudB = 40.0 * Bmag * exp_val;
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

        int nonlinear_mu = (prop == 4) ? 1 : 0;
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
        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;

        if(sigma_laplace > 0.0){
            for(int i=0; i<fe->local_num_nodes; ++i){
                int gi = fe->conn[e][i];
                for(int j=0; j<fe->local_num_nodes; ++j){
                    int gj = fe->conn[e][j];
                    for(int p=0; p<np; ++p){
                        /* sigma_phi -> sigma_laplace に変更 */
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], sigma_laplace);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v);
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
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v);
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
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gn, gi, 0, 0, (double)si * v * inv_dt);
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

    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);
}

/* ============================================================
 * Coil current excitation settings for TEAM P21C-EM1
 *   prop==1 : Coil 1
 *   prop==2 : Coil 2
 *   I1 = +Iamp*sin(wt), I2 = -Iamp*sin(wt)
 * ============================================================ */
static const double I_RMS = 10.0;   /* TEAM benchmark rated current [A rms] */
static const double FREQ_HZ_team21 = 50.0; /* [Hz] */

static inline double get_coil_current_team21c(int prop, double t)
{
    const double omega = 2.0 * M_PI * FREQ_HZ_team21;
    const double Iamp  = sqrt(2.0) * I_RMS;  /* peak value */

    if(prop == 1){
        return  Iamp * sin(omega * t);   /* Coil 1 */
    } else if(prop == 2){
        return -Iamp * sin(omega * t);   /* Coil 2: opposite direction */
    }
    return 0.0;
}

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

        int nonlinear_mu = (prop == 4) ? 1 : 0;

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

                    double d[3] = {
                        x_ip[0] - coil.center[0],
                        x_ip[1] - coil.center[1],
                        x_ip[2] - coil.center[2]
                    };
                    double da = dot3(d, coil.axis);

                    double r[3] = {
                        d[0] - da * coil.axis[0],
                        d[1] - da * coil.axis[1],
                        d[2] - da * coil.axis[2]
                    };

                    double tdir[3] = {
                        coil.axis[1] * r[2] - coil.axis[2] * r[1],
                        coil.axis[2] * r[0] - coil.axis[0] * r[2],
                        coil.axis[0] * r[1] - coil.axis[1] * r[0]
                    };
                    double n_tdir = norm3(tdir);

                    double Js[3] = {0.0, 0.0, 0.0};
                    if(n_tdir > eps_r){
                        double inv_n = 1.0 / n_tdir;
                        Js[0] = J_mag * tdir[0] * inv_n;
                        Js[1] = J_mag * tdir[1] * inv_n;
                        Js[2] = J_mag * tdir[2] * inv_n;
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
        if(sigma_cpl > 0.0){
            for(int j = 0; j < ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                double acc = 0.0;
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_val = x_curr[gn];

                    for(int p = 0; p < np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n],
                            ned->N_edge[e][p][j],
                            sigma_cpl
                        );
                    }

                    double cjn = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)sj * cjn * phi_val;
                }

                monolis->mat.R.B[gj] -= acc;
            }
        }

        /* ==========================================================
         * [5] Phi-equation
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

        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;
        if(sigma_laplace > 0.0){
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int m = 0; m < fe->local_num_nodes; ++m){
                    int gm = fe->conn[e][m];
                    double phi_m = x_curr[gm];

                    for(int p = 0; p < np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n],
                            fe->geo[e][p].grad_N[m],
                            sigma_laplace
                        );
                    }

                    double knm = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += knm * phi_m;
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
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    double loss = 0.0;

    for(int e = 0; e < fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != 3) continue; /* copper shield only */

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double dA_dt[3]   = {0.0, 0.0, 0.0};
            double grad_phi[3]= {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                double dai = (x_curr[gi] - x_prev[gi]) / dt;

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

            double E[3];
            E[0] = -dA_dt[0] - grad_phi[0];
            E[1] = -dA_dt[1] - grad_phi[1];
            E[2] = -dA_dt[2] - grad_phi[2];

            double e2 = E[0]*E[0] + E[1]*E[1] + E[2]*E[2];
            val_ip[p] = Sigma_shield * e2; /* J^2/sigma = sigma * E^2 */
        }

        loss += BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    return loss;
}

