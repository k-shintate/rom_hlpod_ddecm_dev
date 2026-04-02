
#include "core_FOM.h"
#include "3ph_tr_NR.h"
#include "team21c.h"

void ROM_sys_hlpod_fe_set_snap_mat_para_ned(
    double*       	comp_vec,
    BBFE_DATA* fe,
    HLPOD_MAT*      hlpod_mat,
    BBFE_BC*      	bc,
    NEDELEC*		ned,
    const int 		total_num_nodes,
    const int		dof,
    const int 		count)
{
    for(int i = 0; i < total_num_nodes* dof; i++){
        hlpod_mat->snapmat[i][count] = comp_vec[i];
    }
}

void solver_fom_NR_Aphi_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 1;
    const double relaxation = 1.0;

    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            dx[i] = 0.0;
        }

        /* 組み立て */
        set_element_mat_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        monowrap_solve_R(
            &(sys.monolis),
            &(sys.monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        /* Newton update */
        update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), dx, A_delta, phi_delta, sys.fe.total_num_elems);

        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        if (it == max_iter_NR-1 || (conv_v && conv_p)) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }


    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_std_hlpod_set_snapmat_nobc(
                    phi,
                    &(sys.rom_p.hlpod_mat),
                    sys.fe.total_num_nodes,
                    1,
                    (int)(step/sys.vals.snapshot_interval));
        }
        else{
            ROM_std_hlpod_set_snapmat_nobc(
                    phi,
                    &(sys.rom_p.hlpod_mat),
                    sys.monolis_com.n_internal_vertex,
                    1,
                    (int)(step/sys.vals.snapshot_interval));
        }

    }

    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_sys_hlpod_fe_set_snap_mat_para_ned(
                    A,
                    &(sys.fe),
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc),
                    &(sys.ned),
                    sys.fe.total_num_nodes,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
        }
        else{
            ROM_sys_hlpod_fe_set_snap_mat_para_ned(
                    A,
                    &(sys.fe),
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc),
                    &(sys.ned),
                    sys.monolis_com.n_internal_vertex,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
            
        }

    }

    free(dx);
    free(rvec);
    free(rvec_old);

    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
}



void solver_fom_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            dx[i] = 0.0;
        }

        /* 組み立て */
        set_element_mat_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        /* 線形解は dx に入れる（x_curr と分ける！） */
        monowrap_solve_R(
            &(sys.monolis),
            &(sys.monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        /* Newton update */
        update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), dx, A_delta, phi_delta, sys.fe.total_num_elems);

        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        //int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        //int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        if (conv_v && conv_p) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }

    free(dx);
    free(rvec);
}



void solver_fom_NR_Aphi_team21c(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            dx[i] = 0.0;
        }

        /* 組み立て */
        set_element_mat_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi_team21c(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        /* 線形解は dx に入れる（x_curr と分ける！） */
        monowrap_solve_R(
            &(sys.monolis),
            &(sys.monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        /* Newton update */
        update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);
/*
        double shield_loss_inst = calc_copper_shield_loss_EM1(&sys, x_prev, x_curr, sys.vals.dt);

        log_copper_shield_loss_EM1(&sys, step, t, sys.vals.dt, shield_loss_inst);
        log_copper_shield_loss_EM1_cycle_average(&sys, step, t, sys.vals.dt, shield_loss_inst);
*/

        SHIELD_LOSS_DIAG diag;
        diag = calc_copper_shield_loss_EM1_diag(&sys, x_prev, x_curr, sys.vals.dt);

        log_copper_shield_loss_EM1(&sys, step, t, sys.vals.dt, diag.loss_total);
        log_copper_shield_loss_EM1_diag(&sys, step, t, sys.vals.dt, &diag);
        log_copper_shield_loss_EM1_cycle_average(&sys, step, t, sys.vals.dt, diag.loss_total);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), dx, A_delta, phi_delta, sys.fe.total_num_elems);

        set_element_vec_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi_team21c(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        //int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        //int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        if (conv_v && conv_p) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }


    free(dx);
    free(rvec);
}