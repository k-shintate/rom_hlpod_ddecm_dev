
#include "hlpod_core_fe.h"
#include "write_std.h"

void ROM_sys_hlpod_fe_set_snap_mat_para(
    double*       	comp_vec,
    HLPOD_MAT*      hlpod_mat,
    BBFE_BC*      	bc,
    ROM_BC*		    rom_bc,
    const int 		total_num_nodes,
    const int		dof,
    const int 		count)
{
    for(int i = 0; i < total_num_nodes* dof; i++){
        hlpod_mat->snapmat[i][count] = comp_vec[i];
    }

    for(int i = 0; i < total_num_nodes * dof; i++) {
        if( bc->D_bc_exists[i] ) {
            hlpod_mat->snapmat[i][count] = 0.0;
        }
    }

}

void ROM_sys_hlpod_fe_add_Dbc(
    double*         ansvec,
    BBFE_BC*      	bc,
    const int 		total_num_nodes,
    const int       dof)
{
    for(int i = 0; i < total_num_nodes * dof; i++) {
        if( bc->D_bc_exists[i] ) {
            ansvec[i] = bc->imposed_D_val[i];
        }
    }
}

void ROM_sys_hlpod_fe_write_pod_modes_vtk(
        BBFE_DATA*		fe,
        ROM*			rom,
        const int 		n_internal_vertex,
        const int 		num_modes,
        const int 		ndof,
        const char*     label,
        const char*		directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){
            ROM_sys_hlpod_fe_write_pod_modes(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    fe->total_num_nodes,
                    num_modes,
                    ndof,
                    label,
                    directory);
        }
        else{
            ROM_sys_hlpod_fe_write_local_pod_modes(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    num_modes,
                    ndof,
                    label,
                    directory);
        }
    }
    else{
        if(rom->hlpod_vals.bool_global_mode==false){
            ROM_sys_hlpod_fe_write_local_pod_modes(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    num_modes,
                    ndof,
                    label,
                    directory);
        }
        else{
            ROM_sys_hlpod_fe_write_pod_modes(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    num_modes,
                    ndof,
                    label,
                    directory);
        }
    }
}

void ROM_sys_hlpod_fe_write_pod_modes_vtk_diag(
        BBFE_DATA*		fe,
        ROM*			rom,
        const int 		n_internal_vertex,
        const int 		num_modes1,
        const int 		num_modes2,
        const int		ndof1,
        const int 		ndof2,
        const char*     label_v,
        const char*     label_p,
        const char*		directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){
            ROM_sys_hlpod_fe_write_pod_modes_diag(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    num_modes1,
                    label_v,
                    label_p,
                    directory);
        }
        else{
            ROM_sys_hlpod_fe_write_local_pod_modes_diag_id(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    rom->hlpod_mat.node_id,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    num_modes1,
                    label_v,
                    label_p,
                    directory);
        }
    }
    else{
        if(rom->hlpod_vals.bool_global_mode==false){
            ROM_sys_hlpod_fe_write_local_pod_modes_diag(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    num_modes1,
                    label_v,
                    label_p,
                    directory);
        }
        else{
            ROM_sys_hlpod_fe_write_local_pod_modes_diag(
                    fe,
                    rom->hlpod_mat.pod_modes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    num_modes1,
                    label_v,
                    label_p,
                    directory);
        }
    }
}

double ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
        BBFE_DATA*    fe,
        BBFE_BASIS*   basis,
        MONOLIS_COM*  monolis_com,
        double        t,
        const double* comp_vec, // [total_num_nodes]
        const double* pod_comp_vec) // scalar function(x, y, z, t)
{
    double L2_abs_error = 0.0;
    double L2_abs_theo  = 0.0;

    double* val_ip_error;
    double* val_ip_theo;
    double* Jacobian_ip;
    bool* is_internal_elem;
    val_ip_error = BB_std_calloc_1d_double(val_ip_error, basis->num_integ_points);
    val_ip_theo  = BB_std_calloc_1d_double(val_ip_theo , basis->num_integ_points);
    Jacobian_ip  = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);
    is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

    double** local_x;
    local_x   = BB_std_calloc_2d_double(local_x  , fe->local_num_nodes, 3);

    double* fem_local_val;
    fem_local_val = BB_std_calloc_1d_double(fem_local_val, fe->local_num_nodes);

    double* pod_local_val;
    pod_local_val = BB_std_calloc_1d_double(pod_local_val, fe->local_num_nodes);

    monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
        fe->local_num_nodes, fe->conn, is_internal_elem);

    for(int e=0; e<(fe->total_num_elems); e++) {
        if (!is_internal_elem[e]) continue;

        for(int i=0; i<(fe->local_num_nodes); i++) {
            fem_local_val[i] = comp_vec[ fe->conn[e][i] ];
            pod_local_val[i] = pod_comp_vec[ fe->conn[e][i] ];

            local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
            local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
            local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
        }

        BBFE_elemmat_set_Jacobian_array(
                Jacobian_ip,
                basis->num_integ_points,
                e,
                fe);

        for(int p=0; p<(basis->num_integ_points); p++) {
            double fem_val_ip;	double pod_val_ip;
            fem_val_ip = BBFE_std_mapping_scalar(
                    fe->local_num_nodes,
                    fem_local_val,
                    basis->N[p]);

            pod_val_ip = BBFE_std_mapping_scalar(
                    fe->local_num_nodes,
                    pod_local_val,
                    basis->N[p]);

            double x_ip[3];
            BBFE_std_mapping_vector3d(
                    x_ip,
                    fe->local_num_nodes,
                    local_x,
                    basis->N[p]);

            val_ip_error[p] =
                pow( fabs( fem_val_ip - pod_val_ip), 2 );
            val_ip_theo[p] =
                pow( fabs( fem_val_ip), 2 );
        }

        double integ_val_error = BBFE_std_integ_calc(
                basis->num_integ_points,
                val_ip_error,
                basis->integ_weight,
                Jacobian_ip);
        double integ_val_theo  = BBFE_std_integ_calc(
                basis->num_integ_points,
                val_ip_theo,
                basis->integ_weight,
                Jacobian_ip);

        L2_abs_error += integ_val_error;
        L2_abs_theo  += integ_val_theo;
    }

    BB_std_free_1d_double(val_ip_error, basis->num_integ_points);
    BB_std_free_1d_double(val_ip_theo,  basis->num_integ_points);
    BB_std_free_1d_double(Jacobian_ip,  basis->num_integ_points);
    BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);

    BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);
    BB_std_free_1d_double(fem_local_val, fe->local_num_nodes);
    BB_std_free_1d_double(pod_local_val, fe->local_num_nodes);

    monolis_allreduce_R(
            1,
            &L2_abs_error,
            MONOLIS_MPI_SUM,
            monolis_com->comm);

    monolis_allreduce_R(
            1,
            &L2_abs_theo,
            MONOLIS_MPI_SUM,
            monolis_com->comm);

    return ( sqrt(L2_abs_error)/sqrt(L2_abs_theo) );
}

double ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
        BBFE_DATA*      fe,
        BBFE_BASIS*     basis,
        MONOLIS_COM*    monolis_com,
        double          t,
        const double**  fem_vec, // [total_num_nodes][3]
        const double**  rom_vec)  // [total_num_nodes][3]
{
    double L2_abs_error = 0.0;
    double L2_abs_theo  = 0.0;

    double* val_ip_error = BB_std_calloc_1d_double(NULL, basis->num_integ_points);
    double* val_ip_theo  = BB_std_calloc_1d_double(NULL, basis->num_integ_points);
    double* Jacobian_ip  = BB_std_calloc_1d_double(NULL, basis->num_integ_points);
    bool*   is_internal_elem = BB_std_calloc_1d_bool(NULL, fe->total_num_elems);

    double** local_x       = BB_std_calloc_2d_double(NULL, fe->local_num_nodes, 3);
    double** local_val     = BB_std_calloc_2d_double(NULL, fe->local_num_nodes, 3);
    double** local_ans_val = BB_std_calloc_2d_double(NULL, fe->local_num_nodes, 3);

    monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
        fe->local_num_nodes, fe->conn, is_internal_elem);

    for(int e = 0; e < fe->total_num_elems; e++) {
        if (!is_internal_elem[e]) continue;

        for(int i = 0; i < fe->local_num_nodes; i++) {
            int node_id = fe->conn[e][i];

            // 近似解の値を取得
            local_val[i][0] = fem_vec[node_id][0];
            local_val[i][1] = fem_vec[node_id][1];
            local_val[i][2] = fem_vec[node_id][2];

            // 理論解の値を取得
            local_ans_val[i][0] = rom_vec[node_id][0];
            local_ans_val[i][1] = rom_vec[node_id][1];
            local_ans_val[i][2] = rom_vec[node_id][2];

            // 座標を取得
            local_x[i][0] = fe->x[node_id][0];
            local_x[i][1] = fe->x[node_id][1];
            local_x[i][2] = fe->x[node_id][2];
        }

        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip,
            basis->num_integ_points,
            e,
            fe);

        for(int p = 0; p < basis->num_integ_points; p++) {
            double val_ip[3];
            double exact_val[3];
            double x_ip[3];

            // 積分点での近似解を計算
            BBFE_std_mapping_vector3d(
                val_ip,
                fe->local_num_nodes,
                local_val,
                basis->N[p]);

            // 積分点での理論解を計算
            BBFE_std_mapping_vector3d(
                exact_val,
                fe->local_num_nodes,
                local_ans_val,
                basis->N[p]);

            // 積分点での座標を計算
            BBFE_std_mapping_vector3d(
                x_ip,
                fe->local_num_nodes,
                local_x,
                basis->N[p]);

            double diff[3];
            diff[0] = val_ip[0] - exact_val[0];
            diff[1] = val_ip[1] - exact_val[1];
            diff[2] = val_ip[2] - exact_val[2];

            val_ip_error[p] = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
            val_ip_theo[p]  = exact_val[0]*exact_val[0] + exact_val[1]*exact_val[1] + exact_val[2]*exact_val[2];
        }

        double integ_val_error = BBFE_std_integ_calc(
            basis->num_integ_points,
            val_ip_error,
            basis->integ_weight,
            Jacobian_ip);

        double integ_val_theo  = BBFE_std_integ_calc(
            basis->num_integ_points,
            val_ip_theo,
            basis->integ_weight,
            Jacobian_ip);

        L2_abs_error += integ_val_error;
        L2_abs_theo  += integ_val_theo;
    }

    BB_std_free_1d_double(val_ip_error, basis->num_integ_points);
    BB_std_free_1d_double(val_ip_theo,  basis->num_integ_points);
    BB_std_free_1d_double(Jacobian_ip,  basis->num_integ_points);
    BB_std_free_2d_double(local_x,       fe->local_num_nodes, 3);
    BB_std_free_2d_double(local_val,     fe->local_num_nodes, 3);
    BB_std_free_2d_double(local_ans_val, fe->local_num_nodes, 3);
    BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);

    monolis_allreduce_R(
        1,
        &L2_abs_error,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &L2_abs_theo,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    return ( sqrt(L2_abs_error) / sqrt(L2_abs_theo) );
}
