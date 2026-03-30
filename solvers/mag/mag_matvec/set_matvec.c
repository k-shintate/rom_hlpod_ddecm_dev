

#include "set_matvec.h"

void set_reduced_mat_global_para(
    MONOLIS* monolis,
    MONOLIS_COM* monolis_com,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES*         vals,
    NEDELEC* ned,
    double**        mat,
    double**        pod_modes,
    const int 		num_modes,
    const double* x_curr,
    double dt,
    const char* directory)
{
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0/dt;
    const double epsB   = 1.0e-14;

    const int nEdge = fe->total_num_nodes;
    
    const char* fname;
    FILE* fp_in;
    char id[128];
    int tmp;
    int num_internal_elems = 0;
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "graph_nedelec_elem.dat.n_internal");
    fp_in = BBFE_sys_read_fopen(fp_in, fname, directory);
    fscanf(fp_in, "%s %d", id, &(tmp));
    fscanf(fp_in, "%d", &(num_internal_elems));
    fclose(fp_in);

    //for(int e=0; e<fe->total_num_elems; ++e){
    for(int e=0; e<num_internal_elems; ++e){

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
                //monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, (double)(si*sj) * v);

                set_element_reduced_mat(
                    mat,
                    pod_modes,
                    gi,
                    gj,
                    (double)(si*sj) * v,	
                    num_modes);
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
                    //monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0,
                    //    (double)(si*sj) * v * inv_dt);

                    set_element_reduced_mat(
                        mat,
                        pod_modes,
                        gi,
                        gj,
                        (double)(si*sj) * v * inv_dt,	
                        num_modes);
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
                    //monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v);
                    set_element_reduced_mat(mat, pod_modes, gi, gj, v, num_modes);
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
                    //monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v);
                    set_element_reduced_mat(mat, pod_modes, gj, gn, (double)sj * v, num_modes);
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
                    //monolis_add_scalar_to_sparse_matrix_R(monolis, gn, gi, 0, 0, (double)si * v * inv_dt);
                    set_element_reduced_mat(mat, pod_modes, gn, gi, (double)si * v * inv_dt, num_modes);
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * ============================================================ */
void set_reduced_vec_global_para(
    MONOLIS* monolis,
    MONOLIS_COM* monolis_com,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES*         vals,
    NEDELEC* ned,
    double*         rhs,
    double**        pod_modes,
    const int		num_modes,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time,     /* time n+1 */
    const char* directory
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0 / dt;
    const double epsB   = 1.0e-14;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    const char* fname;
    FILE* fp_in;
    char id[128];
    int tmp;
    int num_internal_elems = 0;
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "graph_nedelec_elem.dat.n_internal");
    fp_in = BBFE_sys_read_fopen(fp_in, fname, directory);
    fscanf(fp_in, "%s %d", id, &(tmp));
    fscanf(fp_in, "%d", &(num_internal_elems));
    fclose(fp_in);

    //for(int e=0; e<fe->total_num_elems; ++e){
    for(int e=0; e<num_internal_elems; ++e){

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
                //monolis->mat.R.B[gi] += (double)si * integ;

                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        gi,
                        integ * (double)si,
                        num_modes);
            }
        }
        

        /* ==========================================================
         * [2] A-mass term:
         *   - conductor/core:  -(sigma/dt) M (A_curr - A_prev)
         *   - air stabilization: -(sigma_stab/dt) M (A_curr)  (NO history)
         * ========================================================== */
        if(sigma_mass_A > 0.0){
            int use_history = (prop != 1 || prop != 2|| prop != 3|| prop != 5); /* air uses NO history (penalty only) */

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
                        coeffA = x_curr[gj]; 
                        //(x_curr[gj] - x_prev[gj]); 
                    }

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double mij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)(si*sj) * mij * coeffA * inv_dt;
                }
                //monolis->mat.R.B[gi] -= acc; /* move to residual */
                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        gi,
                        -acc,
                        num_modes);
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
            //monolis->mat.R.B[gi] -= acc;
            set_element_reduced_rhs(
                    rhs,
                    pod_modes,
                    gi,
                    -acc,
                    num_modes);
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
                //monolis->mat.R.B[gj] -= acc;
                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        gj,
                        -acc,
                        num_modes);
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
                //monolis->mat.R.B[gn] -= acc;
                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        gn,
                        -acc,
                        num_modes);
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
                //monolis->mat.R.B[gn] -= acc;
                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        gn,
                        -acc,
                        num_modes);
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}