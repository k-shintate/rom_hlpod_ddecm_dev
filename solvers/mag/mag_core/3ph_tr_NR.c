
#include "3ph_tr_NR.h"

/* Air stabilization factors (tune as needed) */
static const double AIR_PHI_SIGMA_FACTOR  = 1.0e-5; /* phi-phi stabilization relative to copper sigma */
static const double AIR_A_MASS_SIGMA_FACTOR = 1.0e-5; /* A-A mass stabilization relative to copper sigma */

const double Sigma_coil = 5.77e7;                  // 無次元化（典型導電率）
const double Sigma_core =  3.72e3;                  // 無次元化（典型導電率）

/* ============================================================
 * Material coefficient policy (prop: 1..3 coil, 4 core, 5 air)
 * ============================================================ */
void get_sigmas_for_prop(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2 || prop == 3){
        /* coil (conductor) */
        *sigma_mass_A = Sigma_coil;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 4){
        /* core (A-only by default; phi is not solved here) */
        *sigma_mass_A = Sigma_core;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else if(prop == 5){
        /* air: coupling cut; add small penalties for solvability */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


/* --- Nonlinear Material Property (Brauer Law) --- */
/* nu(B) = 100 + 10 * exp(2 * |B|^2) */
void eval_nu_and_dnudB(double Bmag, double* nu, double* dnudB){
    /* 真空の磁気抵抗率 (約 8.0e5) */
    const double NU0 = 1.0 / (4.0 * M_PI * 1.0e-7); 
    
    double B2 = Bmag * Bmag;
    
    /* 指数部が爆発しないようにクランプ (例: B=2.5T程度で制限) */
    /* exp(2 * 6.25) = exp(12.5) approx 2.6e5 */
    if (B2 > 6.25) B2 = 6.25; 

    double exp_val = exp(2.0 * B2);
    *nu = 100.0 + 10.0 * exp_val;
    *dnudB = 40.0 * Bmag * exp_val;

    /* オプション: nu が真空の透磁率を超えないように物理的なリミッタをかける
       (飽和後は真空の傾きに近づくのが一般的であるため) */
    if(*nu > NU0) {
        *nu = NU0;
        *dnudB = 0.0; // あるいは飽和域の傾き
    }
}


double get_reluctivity_nu(double Bmag) {
    double B2 = Bmag * Bmag;
    return 100.0 + 10.0 * exp(2.0 * B2);
}

/* ============================================================
 * B = curl(A) at integration point
 * ============================================================ */
void compute_B_ip(
    const NEDELEC* ned, int e, int p,
    const double* x_curr, int nEdge,
    double B[3]
){
    B[0]=B[1]=B[2]=0.0;
    for(int m=0; m<ned->local_num_edges; ++m){
        int gm = ned->nedelec_conn[e][m];
        int sm = ned->edge_sign[e][m];
        const double* curlNm = ned->curl_N_edge[e][p][m];
        double am = x_curr[gm] * (double)sm;
        B[0] += am * curlNm[0];
        B[1] += am * curlNm[1];
        B[2] += am * curlNm[2];
    }
}

/* integration point coordinates */
void get_interp_coords(
    int e, int p,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    double x_ip[3]
){
    x_ip[0] = 0.0; x_ip[1] = 0.0; x_ip[2] = 0.0;
    for(int n=0; n < fe->local_num_nodes; ++n){
        int global_node_id = fe->conn[e][n];
        double N_val = basis->N[p][n];
        x_ip[0] += N_val * fe->x[global_node_id][0];
        x_ip[1] += N_val * fe->x[global_node_id][1];
        x_ip[2] += N_val * fe->x[global_node_id][2];
    }
}

void update_Aphi_NR(double* x_curr, const double* delta, int n_dof_total, double relaxation){
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] += relaxation * delta[i];
    }
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
void set_element_mat_NR_Aphi(
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
        get_sigmas_for_prop(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

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
                        eval_nu_and_dnudB(Bmag, &nu, &dnudB);
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

void apply_dirichlet_bc_for_A_and_phi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned)
{
    int num_nodes = fe->total_num_nodes;

    const int nen = fe->local_num_nodes;

    static int  is_dir_edge_n = 0;
    is_dir_edge_n  = fe->total_num_nodes;
    bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
    build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);

    /* --------------------------------------------------------
       2. A (Edge) に対する境界条件: A_tan = 0
       外部境界(bc->D_bc_existsが両端でtrue)のエッジを0固定
       -------------------------------------------------------- */
    
    /* 要素タイプに応じたエッジテーブル (Tet/Hex) */
    //const int nen = fe->local_num_nodes;
    int n_local_edges = 0;
    const int (*edge_tbl)[2] = NULL;
    
    if (nen == 4) {
        n_local_edges = 6;
        edge_tbl = tet_edge_conn; /* bbfe_std内で定義されていると仮定 */
    } else if (nen == 8) {
        n_local_edges = 12;
        edge_tbl = hex_edge_conn;
    }

    for (int e = 0; e < fe->total_num_elems; ++e){
        for (int i = 0; i < n_local_edges; ++i){
            int n1_local = edge_tbl[i][0];
            int n2_local = edge_tbl[i][1];

            int gn1 = fe->conn[e][n1_local];
            int gn2 = fe->conn[e][n2_local];

            const int gid = ned->nedelec_conn[e][i];

            /* 両端点がDirichlet境界上にあるエッジは、境界エッジとみなす */
            if (!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;
            
            if(!is_dir_edge[gid]) continue; 

                int ged = ned->nedelec_conn[e][i]; /* Edge Global ID */
                
                /* A = 0.0 を設定 (Monolisの関数で行固定) */
                monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    ged,  /* 行番号 = Edge ID */
                    0,    /* DOF index (スカラーなので0) */
                    0.0   /* 指定値 */
                );
            
        }
    }

    
    double* node_is_conductor = (double*)calloc(num_nodes, sizeof(double));
    
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

        if(prop == 1 || prop == 2 || prop == 3 || prop == 4){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                //node_is_conductor[gn] = 1; 
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

/*
    int count = 0;
    if(monolis_mpi_get_global_my_rank()==0){
        for (int i = 0; i < num_nodes; ++i){
            if (node_is_conductor[i] == 4) {
                count++;
                if(count == 1){
                    monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    i, 
                    0, 
                    0.0
                );

                printf("\n\n\n\n\nadd D_bc\n\n\n\n\n");

                }
                else{
                }
            }
        }
    }

*/
    free(node_is_conductor);
}

void debug_max_B_and_nu_core(
    FE_SYSTEM* sys,
    const double* x_curr,
    int n_dof_total,
    int it, int step, double t,
    const char* directory
){
    BBFE_DATA*  fe   = &(sys->fe);
    BBFE_BASIS* basis= &(sys->basis);
    NEDELEC*    ned  = &(sys->ned);

    const int np   = basis->num_integ_points;
    const double epsB = 1.0e-14;

    /* local maxima */
    double local_max_B  = 0.0;
    double local_max_nu = 0.0;

    /* NaN/Inf counter (doubleでreduce) */
    double local_bad_nu = 0.0;

    double local_max_x = 0.0;
    for(int i=0;i<n_dof_total;++i){
        double a = fabs(x_curr[i]);
        if(a > local_max_x) local_max_x = a;
    }

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop != 4) continue; /* coreだけ */

        for(int p=0; p<np; ++p){
            double B[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B);

            double Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
            if(Bmag > local_max_B) local_max_B = Bmag;

            double nu = get_reluctivity_nu(fmax(Bmag, epsB));
            if(!isfinite(nu)){
                local_bad_nu += 1.0;
            } else {
                if(nu > local_max_nu) local_max_nu = nu;
            }
        }
    }

    /* MPI reduce (MAX/SUM) */
    double g_max_B  = local_max_B;
    double g_max_nu = local_max_nu;
    double g_bad_nu = local_bad_nu;
    double g_max_x  = local_max_x;

    monolis_allreduce_R(1, &g_max_B,  MONOLIS_MPI_MAX, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_max_nu, MONOLIS_MPI_MAX, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_bad_nu, MONOLIS_MPI_SUM, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_max_x,  MONOLIS_MPI_MAX, sys->monolis_com.comm);

    if(sys->monolis_com.my_rank == 0){
        printf("[DBG step=%d it=%d t=%.6e] core: max|x|=%.6e  max|B|=%.6e  max(nu)=%.6e  bad_nu(count)=%.0f\n",
               step, it, t, g_max_x, g_max_B, g_max_nu, g_bad_nu);
    }

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = BBFE_sys_write_add_fopen(fp, "NR_prop.txt", sys->cond.directory);
		fprintf(fp, "[DBG step=%d it=%d t=%.6e] core: max|x|=%.6e  max|B|=%.6e  max(nu)=%.6e  bad_nu(count)=%.0f\n",
               step, it, t, g_max_x, g_max_B, g_max_nu, g_bad_nu);
		fclose(fp);
	}
}

/* ============================================================
 * Voltage excitation settings
 * ============================================================ */
static const double V_AMP  = 100.0;   /* [V] phase voltage amplitude (temporary) */
static const double R_COIL = 0.1;     /* [Ohm] equivalent phase winding resistance */

/* 各相電流（時間ステップ更新） */
double g_phase_current2[3]      = {0.0, 0.0, 0.0};
double g_phase_current2_prev[3] = {0.0, 0.0, 0.0};

/* ============================================================
 * Applied phase voltage
 * ============================================================ */
static inline double get_phase_voltage(int phase_idx, double t)
{
    double omega = 2.0 * M_PI * FREQ_HZ;

    if(phase_idx == 0) return V_AMP * sin(omega * t);
    if(phase_idx == 1) return V_AMP * sin(omega * t - 2.0 * M_PI / 3.0);
    if(phase_idx == 2) return V_AMP * sin(omega * t + 2.0 * M_PI / 3.0);

    return 0.0;
}

/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * Voltage excitation version
 * ============================================================ */
void set_element_vec_NR_Aphi(
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

    if(current_time == dt){
        printf("intialize g_phase_current2");
        for(int k = 0; k < 3; k++){
            double Vapp = get_phase_voltage(k, current_time);
            g_phase_current2[k] = Vapp / R_COIL;
        }
    }

    for(int e = 0; e < fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        int nonlinear_mu = (prop == 4) ? 1 : 0;

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation)
         *     current source amplitude is no longer prescribed by CURRENT_AMP
         *     use updated phase current g_phase_current2[phase]
         * ========================================================== */
        if(is_coil){
            int phase_idx = prop - 1; /* prop 1,2,3 -> U,V,W */
            double I_t = 0.0;
            if(0 <= phase_idx && phase_idx < 3){
                I_t = g_phase_current2[phase_idx];
            }

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
                        nu_val = get_reluctivity_nu(fmax(norm3(B), epsB));
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

void log_accuracy_metrics(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    int step,
    double t,
    double dt,
    double P_input
){
    BBFE_DATA* fe     = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC* ned      = &(sys->ned);

    const int np    = basis->num_integ_points;
    const int nEdge = fe->total_num_nodes;
    const double inv_dt = 1.0 / dt;

    /* ============================================================
     * TRANSMAG comparison settings
     * ============================================================ */
    /* 飽和判定しきい値 [T]
     * 実機BHに合わせるなら別途調整してください
     */
    const double B_SAT_THRESHOLD = 1.5;

    /* 1周期 */
    const double T_period = 1.0 / FREQ_HZ;

    /* ============================================================
     * Probe Settings (Monitoring Points)
     * ============================================================ */
    const double search_radius = 2.0e-3; // 探索半径 2 mm
    const double target_pos1[3] = { 0.0025, 0.0350, 0.0200 }; // 測定点1
    const double target_pos2[3] = { 0.0350, 0.0350, 0.0200 }; // 測定点2
    const double target_pos3[3] = { 0.0675, 0.0350, 0.0200 }; // 測定点3

    /* ============================================================
     * Data Accumulators Definition
     *
     * [MPI_SUM target] (Integrals)
     *  0-2   : Reserved
     *  3-5   : V_fem (EMF) for U, V, W
     *  6     : Total Loss
     *  7     : Total Mag Energy Change-like quantity
     *  8-10  : Mag Energy Change-like quantity (Coil, Core, Air)
     * 11-13  : Volume (Coil, Core, Air)
     * 16-18  : Loss Breakdown (Coil, Core, Air)
     * 19     : Core saturated volume (|B| > B_SAT_THRESHOLD)
     * 20-22  : Probe1 (Nu, Dist, B)  <-- local only, later MAX reduce
     * 23-25  : Probe2 (Nu, Dist, B)
     * 26-28  : Probe3 (Nu, Dist, B)
     *
     * [MPI_MAX target]
     * 14     : Max B_mag (Global)
     * 15     : Max Nu (Air Region check)
     * 20-28  : Probe values
     * ============================================================ */
    const int N_LOG = 29;
    double local_vals[29];
    for(int i = 0; i < N_LOG; i++) local_vals[i] = 0.0;

    /* 真空の磁気抵抗率 */
    const double NU0 = 1.0 / (4.0 * M_PI * 1.0e-7);

    /* ============================================================
     * Source Current
     * ============================================================ */
    double I_src[3] = {0.0, 0.0, 0.0};
    double omega = 2.0 * M_PI * FREQ_HZ;
    for(int k = 0; k < 3; k++){
        COIL_INFO info;
        get_coil_info(k + 1, &info);
        I_src[k] = CURRENT_AMP * sin(omega * t + info.phase_shift);
    }

    /* ============================================================
     * Read internal element count
     * ============================================================ */
    const char* fname;
    FILE* fp_in;
    char id[128];
    int tmp;
    int num_internal_elems = 0;

    fname = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        "graph_nedelec_elem.dat.n_internal"
    );
    fp_in = BBFE_sys_read_fopen(fp_in, fname, sys->cond.directory);
    fscanf(fp_in, "%s %d", id, &(tmp));
    fscanf(fp_in, "%d", &(num_internal_elems));
    fclose(fp_in);

    /* ============================================================
     * Element Loop
     * ============================================================ */
    for(int e = 0; e < num_internal_elems; ++e){
        int prop = ned->elem_prop[e];

        /* Region Identification */
        int region_type = 0; // 0:Other, 1:Coil, 2:Core, 3:Air
        if(prop == 1 || prop == 2 || prop == 3) region_type = 1;
        else if(prop == 4) region_type = 2;
        else if(prop == 5) region_type = 3;

        double sigma_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop(prop, &sigma_A, &sigma_cpl, &sigma_phi);

        /* Coil Info */
        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);
        int phase_idx = -1;
        if(prop == 1) phase_idx = 0;
        if(prop == 2) phase_idx = 1;
        if(prop == 3) phase_idx = 2;

        /* Element centroid for probing */
        double center[3] = {0.0, 0.0, 0.0};
        for(int n = 0; n < fe->local_num_nodes; ++n){
            int node_id = fe->conn[e][n];
            center[0] += fe->x[node_id][0];
            center[1] += fe->x[node_id][1];
            center[2] += fe->x[node_id][2];
        }
        center[0] /= fe->local_num_nodes;
        center[1] /= fe->local_num_nodes;
        center[2] /= fe->local_num_nodes;

        /* Distance to probes */
        double dist1 = sqrt(
            (center[0]-target_pos1[0])*(center[0]-target_pos1[0]) +
            (center[1]-target_pos1[1])*(center[1]-target_pos1[1]) +
            (center[2]-target_pos1[2])*(center[2]-target_pos1[2])
        );
        double dist2 = sqrt(
            (center[0]-target_pos2[0])*(center[0]-target_pos2[0]) +
            (center[1]-target_pos2[1])*(center[1]-target_pos2[1]) +
            (center[2]-target_pos2[2])*(center[2]-target_pos2[2])
        );
        double dist3 = sqrt(
            (center[0]-target_pos3[0])*(center[0]-target_pos3[0]) +
            (center[1]-target_pos3[1])*(center[1]-target_pos3[1]) +
            (center[2]-target_pos3[2])*(center[2]-target_pos3[2])
        );

        /* Integration Loop */
        for(int p = 0; p < np; ++p){
            double w_detJ = basis->integ_weight[p] * fe->geo[e][p].Jacobian;

            /* Volume accumulation */
            if(region_type == 1) local_vals[11] += w_detJ;
            if(region_type == 2) local_vals[12] += w_detJ;
            if(region_type == 3) local_vals[13] += w_detJ;

            /* B-field */
            double B_curr[3], B_prev[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B_curr);
            compute_B_ip(ned, e, p, x_prev, nEdge, B_prev);

            double Bmag = norm3(B_curr);

            /* Global max B */
            if(Bmag > local_vals[14]) local_vals[14] = Bmag;

            /* Reluctivity */
            double current_nu;
            if(region_type == 2){
                current_nu = get_reluctivity_nu(fmax(Bmag, 1.0e-14));
            } else {
                current_nu = NU0;
            }

            /* Air-only max nu check */
            if(region_type == 3){
                if(current_nu > local_vals[15]) local_vals[15] = current_nu;
            }

            /* Core saturated volume */
            if(region_type == 2 && Bmag > B_SAT_THRESHOLD){
                local_vals[19] += w_detJ;
            }

            /* Probe capture (core only) */
            if(region_type == 2){
                if(dist1 < search_radius){
                    local_vals[20] = current_nu;
                    local_vals[21] = dist1;
                    local_vals[22] = Bmag;
                }
                if(dist2 < search_radius){
                    local_vals[23] = current_nu;
                    local_vals[24] = dist2;
                    local_vals[25] = Bmag;
                }
                if(dist3 < search_radius){
                    local_vals[26] = current_nu;
                    local_vals[27] = dist3;
                    local_vals[28] = Bmag;
                }
            }

            /* Reconstruct A and grad(phi) */
            double A_curr[3] = {0.0, 0.0, 0.0};
            double A_prev[3] = {0.0, 0.0, 0.0};
            double grad_phi[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                const double* N = ned->N_edge[e][p][i];
                double ac = x_curr[gi] * si;
                double ap = x_prev[gi] * si;
                for(int d = 0; d < 3; d++){
                    A_curr[d] += ac * N[d];
                    A_prev[d] += ap * N[d];
                }
            }

            if(sigma_cpl > 0.0 || sigma_phi > 0.0){
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double ph = x_curr[gn];
                    const double* dN = fe->geo[e][p].grad_N[n];
                    for(int d = 0; d < 3; d++){
                        grad_phi[d] += ph * dN[d];
                    }
                }
            }

            /* E field */
            double E_vec[3];
            for(int d = 0; d < 3; d++){
                E_vec[d] = -(A_curr[d] - A_prev[d]) * inv_dt - grad_phi[d];
            }

            /* Coil source term / EMF integration */
            if(is_coil && phase_idx >= 0){
                double x_ip[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                double d_vec[3] = {
                    x_ip[0] - coil.center[0],
                    x_ip[1] - coil.center[1],
                    x_ip[2] - coil.center[2]
                };
                double da = dot3(d_vec, coil.axis);
                double r_vec[3] = {
                    d_vec[0] - da * coil.axis[0],
                    d_vec[1] - da * coil.axis[1],
                    d_vec[2] - da * coil.axis[2]
                };
                double t_vec[3] = {
                    coil.axis[1]*r_vec[2] - coil.axis[2]*r_vec[1],
                    coil.axis[2]*r_vec[0] - coil.axis[0]*r_vec[2],
                    coil.axis[0]*r_vec[1] - coil.axis[1]*r_vec[0]
                };
                double t_norm = norm3(t_vec);

                if(t_norm > 1.0e-12){
                    double J0_mag = coil.turns / coil.area;
                    double J0[3] = {
                        J0_mag * t_vec[0] / t_norm,
                        J0_mag * t_vec[1] / t_norm,
                        J0_mag * t_vec[2] / t_norm
                    };
                    local_vals[3 + phase_idx] += dot3(E_vec, J0) * w_detJ;
                }
            }

            /* Joule loss */
            double dLoss = 0.0;
            if(sigma_A > 0.0){
                dLoss = sigma_A * dot3(E_vec, E_vec) * w_detJ;
            }

            local_vals[6] += dLoss; /* Total loss */

            if(region_type == 1) local_vals[16] += dLoss; // Coil
            if(region_type == 2) local_vals[17] += dLoss; // Core
            if(region_type == 3) local_vals[18] += dLoss; // Air

            /* Magnetic energy-change-like quantity */
            double H_vec[3] = {
                current_nu * B_curr[0],
                current_nu * B_curr[1],
                current_nu * B_curr[2]
            };
            double dB_dt[3] = {
                (B_curr[0] - B_prev[0]) * inv_dt,
                (B_curr[1] - B_prev[1]) * inv_dt,
                (B_curr[2] - B_prev[2]) * inv_dt
            };
            double dW_val = dot3(H_vec, dB_dt) * w_detJ;

            local_vals[7] += dW_val;
            if(region_type == 1) local_vals[8]  += dW_val;
            if(region_type == 2) local_vals[9]  += dW_val;
            if(region_type == 3) local_vals[10] += dW_val;
        }
    }

    /* ============================================================
     * MPI Reduce
     * ============================================================ */
    double global_vals[29];

    /* SUM-reduce targets */
    double val_sum[29];
    for(int i = 0; i < 29; i++) val_sum[i] = local_vals[i];

    /* MAX targets should not be included in SUM */
    val_sum[14] = 0.0;
    val_sum[15] = 0.0;
    for(int i = 20; i < 29; i++) val_sum[i] = 0.0;

    monolis_allreduce_R(29, val_sum, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    /* MAX-reduce targets */
    double val_max[29];
    for(int i = 0; i < 29; i++) val_max[i] = 0.0;

    val_max[14] = local_vals[14];
    val_max[15] = local_vals[15];
    for(int i = 20; i < 29; i++) val_max[i] = local_vals[i];

    monolis_allreduce_R(29, val_max, MONOLIS_MPI_MAX, sys->monolis_com.comm);

    /* Merge */
    for(int i = 0; i < 29; i++) global_vals[i] = val_sum[i];
    global_vals[14] = val_max[14];
    global_vals[15] = val_max[15];
    for(int i = 20; i < 29; i++) global_vals[i] = val_max[i];

    /* ============================================================
     * Rank 0 output
     * ============================================================ */
    if(sys->monolis_com.my_rank == 0){

        /* --------------------------------------------
         * Terminal quantities
         * -------------------------------------------- */
        double V_fem[3] = { global_vals[3], global_vals[4], global_vals[5] };

        const double R_DC_COIL = 0.1; // [Ohm]
        double V_term[3];
        for(int k = 0; k < 3; k++){
            V_term[k] = R_DC_COIL * I_src[k] - V_fem[k];
        }

        /* Probe data */
        double p1_nu = global_vals[20], p1_B = global_vals[22];
        double p2_nu = global_vals[23], p2_B = global_vals[25];
        double p3_nu = global_vals[26], p3_B = global_vals[28];

        /* Loss totals */
        double Loss_Coil  = global_vals[16];
        double Loss_Core  = global_vals[17];
        double Loss_Air   = global_vals[18];
        double Loss_Total = Loss_Coil + Loss_Core + Loss_Air;

        /* Core volume and saturation ratio */
        double CoreVol = global_vals[12];
        double CoreSatVol = global_vals[19];
        double CoreSatRatio = 0.0;
        if(CoreVol > 1.0e-30){
            CoreSatRatio = CoreSatVol / CoreVol;
        }

        /* --------------------------------------------
         * Running values for TRANSMAG comparison
         * -------------------------------------------- */
        static int initialized = 0;

        static double run_max_B = -1.0;
        static double run_max_B_time = 0.0;
        static int    run_max_B_step = -1;

        static double run_maxB_Iu = 0.0;
        static double run_maxB_Iv = 0.0;
        static double run_maxB_Iw = 0.0;

        static double run_maxB_P1_B = 0.0;
        static double run_maxB_P2_B = 0.0;
        static double run_maxB_P3_B = 0.0;

        static double run_maxB_CoreSatRatio = 0.0;
        static double run_maxB_LossCore = 0.0;
        static double run_maxB_LossTotal = 0.0;

        /* Last-cycle accumulation */
        static double cycle_start_time = 0.0;
        static double sum_Iu2_dt = 0.0;
        static double sum_Iv2_dt = 0.0;
        static double sum_Iw2_dt = 0.0;
        static double sum_Pin_dt = 0.0;
        static double sum_LossCore_dt = 0.0;
        static double sum_LossTotal_dt = 0.0;
        static double sum_CoreSatRatio_dt = 0.0;
        static double cycle_max_B = 0.0;

        if(step == 0 || !initialized){
            initialized = 1;

            run_max_B = -1.0;
            run_max_B_time = t;
            run_max_B_step = step;

            run_maxB_Iu = run_maxB_Iv = run_maxB_Iw = 0.0;
            run_maxB_P1_B = run_maxB_P2_B = run_maxB_P3_B = 0.0;
            run_maxB_CoreSatRatio = 0.0;
            run_maxB_LossCore = 0.0;
            run_maxB_LossTotal = 0.0;

            cycle_start_time = t;
            sum_Iu2_dt = sum_Iv2_dt = sum_Iw2_dt = 0.0;
            sum_Pin_dt = 0.0;
            sum_LossCore_dt = 0.0;
            sum_LossTotal_dt = 0.0;
            sum_CoreSatRatio_dt = 0.0;
            cycle_max_B = 0.0;
        }

        /* Update running max-B snapshot */
        if(global_vals[14] > run_max_B){
            run_max_B = global_vals[14];
            run_max_B_time = t;
            run_max_B_step = step;

            run_maxB_Iu = I_src[0];
            run_maxB_Iv = I_src[1];
            run_maxB_Iw = I_src[2];

            run_maxB_P1_B = p1_B;
            run_maxB_P2_B = p2_B;
            run_maxB_P3_B = p3_B;

            run_maxB_CoreSatRatio = CoreSatRatio;
            run_maxB_LossCore = Loss_Core;
            run_maxB_LossTotal = Loss_Total;
        }

        /* Accumulate cycle data */
        sum_Iu2_dt += I_src[0] * I_src[0] * dt;
        sum_Iv2_dt += I_src[1] * I_src[1] * dt;
        sum_Iw2_dt += I_src[2] * I_src[2] * dt;
        sum_Pin_dt += P_input * dt;
        sum_LossCore_dt += Loss_Core * dt;
        sum_LossTotal_dt += Loss_Total * dt;
        sum_CoreSatRatio_dt += CoreSatRatio * dt;

        if(global_vals[14] > cycle_max_B){
            cycle_max_B = global_vals[14];
        }

        /* --------------------------------------------
         * Console output
         * -------------------------------------------- */
        printf(
            "  [Monitor] Step:%d t=%.4e | MaxB: %.4f [T] | SatRatio: %.4e | "
            "P1(Core): B=%.4f Nu=%.4e | P2(Core): B=%.4f Nu=%.4e | P3(Core): B=%.4f Nu=%.4e\n",
            step, t, global_vals[14], CoreSatRatio,
            p1_B, p1_nu, p2_B, p2_nu, p3_B, p3_nu
        );

        /* --------------------------------------------
         * Instantaneous CSV
         * -------------------------------------------- */
        FILE* fp_out;
        fp_out = BBFE_sys_write_add_fopen(fp_out, "terminal_voltage_log.csv", sys->cond.directory);

        if(step == 0){
            fprintf(
                fp_out,
                "Step,Time,"
                "I_u,I_v,I_w,"
                "V_t_u,V_t_v,V_t_w,"
                "EMF_u,EMF_v,EMF_w,"
                "Power_In,"
                "Max_B,Max_Nu,"
                "Loss_Coil,Loss_Core,Loss_Air,Loss_Total,"
                "Core_Vol,Core_SatVol,Core_SatRatio,"
                "P1_Nu,P1_B,P2_Nu,P2_B,P3_Nu,P3_B,"
                "Run_Max_B,Run_Max_B_Time,Run_Max_B_Step,"
                "Run_MaxB_Iu,Run_MaxB_Iv,Run_MaxB_Iw,"
                "Run_MaxB_P1_B,Run_MaxB_P2_B,Run_MaxB_P3_B,"
                "Run_MaxB_CoreSatRatio,Run_MaxB_LossCore,Run_MaxB_LossTotal\n"
            );
        }

        fprintf(
            fp_out,
            "%d,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,"
            "%.6e,%.6e,"
            "%.6e,%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%d,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e\n",
            step, t,
            I_src[0], I_src[1], I_src[2],
            V_term[0], V_term[1], V_term[2],
            -V_fem[0], -V_fem[1], -V_fem[2],
            P_input,
            global_vals[14], global_vals[15],
            Loss_Coil, Loss_Core, Loss_Air, Loss_Total,
            CoreVol, CoreSatVol, CoreSatRatio,
            p1_nu, p1_B, p2_nu, p2_B, p3_nu, p3_B,
            run_max_B, run_max_B_time, run_max_B_step,
            run_maxB_Iu, run_maxB_Iv, run_maxB_Iw,
            run_maxB_P1_B, run_maxB_P2_B, run_maxB_P3_B,
            run_maxB_CoreSatRatio, run_maxB_LossCore, run_maxB_LossTotal
        );
        fclose(fp_out);

        /* --------------------------------------------
         * Per-cycle summary CSV
         * TRANSMAG 比較では、定常末期の1周期量が有用
         * -------------------------------------------- */
        if((t - cycle_start_time + dt) >= T_period - 1.0e-14){
            double cycle_len = t - cycle_start_time + dt;
            if(cycle_len < 1.0e-14) cycle_len = T_period;

            double Irms_u = sqrt(sum_Iu2_dt / cycle_len);
            double Irms_v = sqrt(sum_Iv2_dt / cycle_len);
            double Irms_w = sqrt(sum_Iw2_dt / cycle_len);

            double Pin_avg = sum_Pin_dt / cycle_len;
            double LossCore_avg = sum_LossCore_dt / cycle_len;
            double LossTotal_avg = sum_LossTotal_dt / cycle_len;
            double CoreSatRatio_avg = sum_CoreSatRatio_dt / cycle_len;

            FILE* fp_sum;
            fp_sum = BBFE_sys_write_add_fopen(fp_sum, "transmag_cycle_summary.csv", sys->cond.directory);

            if(step == 0 || fabs(cycle_start_time) < 1.0e-14){
                fprintf(
                    fp_sum,
                    "CycleStart,CycleEnd,CycleLength,"
                    "Irms_u,Irms_v,Irms_w,"
                    "Pin_avg,LossCore_avg,LossTotal_avg,"
                    "CoreSatRatio_avg,CycleMax_B,"
                    "Run_Max_B,Run_Max_B_Time,Run_Max_B_Step\n"
                );
            }

            fprintf(
                fp_sum,
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,"
                "%.6e,%.6e,%d\n",
                cycle_start_time, t + dt, cycle_len,
                Irms_u, Irms_v, Irms_w,
                Pin_avg, LossCore_avg, LossTotal_avg,
                CoreSatRatio_avg, cycle_max_B,
                run_max_B, run_max_B_time, run_max_B_step
            );
            fclose(fp_sum);

            /* reset cycle accumulators */
            cycle_start_time = t + dt;
            sum_Iu2_dt = 0.0;
            sum_Iv2_dt = 0.0;
            sum_Iw2_dt = 0.0;
            sum_Pin_dt = 0.0;
            sum_LossCore_dt = 0.0;
            sum_LossTotal_dt = 0.0;
            sum_CoreSatRatio_dt = 0.0;
            cycle_max_B = 0.0;
        }

        /* --------------------------------------------
         * One-line TRANSMAG summary
         * 飽和最大時刻の比較用
         * -------------------------------------------- */
        FILE* fp_tm;
        fp_tm = BBFE_sys_write_fopen(fp_tm, "transmag_peak_summary.txt", sys->cond.directory);
        fprintf(fp_tm, "Run_Max_B           = %.12e\n", run_max_B);
        fprintf(fp_tm, "Run_Max_B_Time      = %.12e\n", run_max_B_time);
        fprintf(fp_tm, "Run_Max_B_Step      = %d\n", run_max_B_step);
        fprintf(fp_tm, "Run_MaxB_Iu         = %.12e\n", run_maxB_Iu);
        fprintf(fp_tm, "Run_MaxB_Iv         = %.12e\n", run_maxB_Iv);
        fprintf(fp_tm, "Run_MaxB_Iw         = %.12e\n", run_maxB_Iw);
        fprintf(fp_tm, "Run_MaxB_P1_B       = %.12e\n", run_maxB_P1_B);
        fprintf(fp_tm, "Run_MaxB_P2_B       = %.12e\n", run_maxB_P2_B);
        fprintf(fp_tm, "Run_MaxB_P3_B       = %.12e\n", run_maxB_P3_B);
        fprintf(fp_tm, "Run_MaxB_CoreSatRatio = %.12e\n", run_maxB_CoreSatRatio);
        fprintf(fp_tm, "Run_MaxB_LossCore   = %.12e\n", run_maxB_LossCore);
        fprintf(fp_tm, "Run_MaxB_LossTotal  = %.12e\n", run_maxB_LossTotal);
        fclose(fp_tm);
    }
}





/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * ============================================================ */
void set_element_vec_NR_Aphi_current(
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

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;

        get_sigmas_for_prop(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);
        int nonlinear_mu = (prop == 4) ? 1 : 0;

        /* coil source uses prop 1..3 */
        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation) -> contributes to A residual (+)
         * ========================================================== */
        
        if(is_coil){
            double omega = 2.0 * M_PI * FREQ_HZ;
            double I_t   = CURRENT_AMP * sin(omega * current_time + coil.phase_shift);
            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double d[3] = { x_ip[0]-coil.center[0], x_ip[1]-coil.center[1], x_ip[2]-coil.center[2] };
                    double da = dot3(d, coil.axis);
                    double r[3] = { d[0]-da*coil.axis[0], d[1]-da*coil.axis[1], d[2]-da*coil.axis[2] };

                    double tdir[3] = {
                        coil.axis[1]*r[2] - coil.axis[2]*r[1],
                        coil.axis[2]*r[0] - coil.axis[0]*r[2],
                        coil.axis[0]*r[1] - coil.axis[1]*r[0]
                    };
                    double n_tdir = norm3(tdir);

                    double Js[3] = {0.0, 0.0, 0.0};
                    if(n_tdir > eps_r){
                        double inv_n = 1.0/n_tdir;
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
         * [2] A-mass term:
         *   - conductor/core:  -(sigma/dt) M (A_curr - A_prev)
         *   - air stabilization: -(sigma_stab/dt) M (A_curr)  (NO history)
         * ========================================================== */
        if(sigma_mass_A > 0.0){
            int use_history = (prop == 1 || prop == 2|| prop == 3|| prop == 5); /* air uses NO history (penalty only) */

            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                double acc = 0.0;
                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj_e = ned->nedelec_conn[e][j];
                    int gj = gj_e;
                    int sj = ned->edge_sign[e][j];

                    double coeffA;
                    if(use_history){
                        coeffA = (x_curr[gj] - x_prev[gj]); /* (A^{n+1}-A^n) */
                    }else{
                        coeffA = //x_curr[gj]; 
                        (x_curr[gj] - x_prev[gj]); 
                    }

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double mij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)(si*sj) * mij * coeffA * inv_dt;
                }
                monolis->mat.R.B[gi] -= acc; /* move to residual */
            }
        }

        /* ==========================================================
         * [3] Curl-curl stiffness term: -(K(A_curr) * A_curr)
         * ========================================================== */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];
                double a_val = x_curr[gj];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    double nu_val;
                    if(nonlinear_mu){
                        double B[3]; compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        nu_val = get_reluctivity_nu(fmax(norm3(B), epsB));
                    }else{
                        nu_val = NU_LIN;
                    }
                    val_ip_C[p] = nu_val * dot3(ci, cj);
                }
                double kij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)(si*sj) * kij * a_val;
            }
            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [4] A-Phi coupling in A-equation: -( sigma_cpl * C * phi )
         * ========================================================== */
        if(sigma_cpl > 0.0){
            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                double acc = 0.0;
                for(int n=0; n<fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_val = x_curr[gn];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], ned->N_edge[e][p][j], sigma_cpl);
                    }
                    double cjn = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)sj * cjn * phi_val;
                }
                monolis->mat.R.B[gj] -= acc;
            }
        }

        /* ==========================================================
         * [5] Phi-equation:
         *    (a) - sigma_cpl/dt * C^T * (A_curr - A_prev)  (air: sigma_cpl=0)
         *    (b) - sigma_phi * Kphi * phi_curr             (air: sigma_phi>0 stabilization)
         * ========================================================== */

        /* (a) dA/dt coupling term */
        if(sigma_cpl > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int i=0; i<ned->local_num_edges; ++i){
                    int gi = ned->nedelec_conn[e][i];
                    int si = ned->edge_sign[e][i];
                    double da = x_curr[gi] - x_prev[gi];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], fe->geo[e][p].grad_N[n], sigma_cpl);
                    }
                    double gin = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)si * gin * da * inv_dt;
                }
                monolis->mat.R.B[gn] -= acc;
            }
        }

        /* (b) phi-phi Laplace/penalty term */
        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;
        if(sigma_laplace > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int m=0; m<fe->local_num_nodes; ++m){
                    int gm = fe->conn[e][m];
                    double phi_m = x_curr[gm];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], fe->geo[e][p].grad_N[m], sigma_laplace);
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