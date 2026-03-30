

#include "set_matvec.h"

/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
		MONOLIS*        monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
        double**        mat,
        double**        pod_modes,
        const int 		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double A[4][4];

    bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {

					double tau = BBFE_elemmat_fluid_sups_coef(
							vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

					BBFE_elemmat_fluid_sups_mat(
							A, basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], 
							v_ip[p], vals->density, vals->viscosity, tau, vals->dt);

					for(int a=0; a<4; a++){
						for(int b=0; b<4; b++) {
							val_ip[a][b][p] = A[a][b];
							A[a][b] = 0.0;
						}
					}
				}

				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        set_element_reduced_mat(
                                mat,
                                pod_modes,
                                fe->conn[e][i]*4 + a,
                                fe->conn[e][j]*4 + b,
                                integ_val,	
                                num_modes);
                    }
				}
			}
		}
	}

	BB_std_free_3d_double(val_ip, 4 , 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}


void set_D_bc_global_para(
		MONOLIS*        monolis,
        MONOLIS_COM*    monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
    	BBFE_BC*     	bc,
        double*         rhs,
        double**        pod_modes,
        const int 		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double A[4][4];

    bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {

					double tau = BBFE_elemmat_fluid_sups_coef(
							vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

					BBFE_elemmat_fluid_sups_mat(
							A, basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], 
							v_ip[p], vals->density, vals->viscosity, tau, vals->dt);

					for(int a=0; a<4; a++){
						for(int b=0; b<4; b++) {
							val_ip[a][b][p] = A[a][b];
							A[a][b] = 0.0;
						}
					}
				}

				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

							set_element_reduced_rhs_Dbc(
                                bc,
                                rhs,
                                pod_modes,
                                fe->conn[e][i]*4 + a,
                                fe->conn[e][j]*4 + b,
                                integ_val,
                                num_modes);

					}
				}
			}
		}
	}


	BB_std_free_3d_double(val_ip     , 4 , 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}

void set_reduced_vec_global_para(
		MONOLIS*        monolis,        
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
        double*         rhs,
        double**        pod_modes,
        const int		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
        if (!is_internal_elem[e]) continue;
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			double integ_val[4];

			for(int p=0; p<np; p++) {
				double tau = BBFE_elemmat_fluid_sups_coef(
						vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

				double vec[4];
				BBFE_elemmat_fluid_sups_vec(
						vec, basis->N[p][i], fe->geo[e][p].grad_N[i],
						v_ip[p], vals->density, tau, vals->dt);

				for(int d=0; d<4; d++) {
					val_ip[d][p] = vec[d];
				}
			}

			for(int d=0; d<4; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        4 * fe->conn[e][i] + d,
                        integ_val[d],
                        num_modes);
			}
		}
	}

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}


void allreduce_global_para(
    MONOLIS_COM*	monolis_com,
    double**        mat,
    double*         rhs,
    const int		num_modes)
{
    double* vec;
	vec = BB_std_calloc_1d_double(vec, num_modes);

	for(int k2 = 0; k2 < num_modes; k2++){
        vec[k2] = rhs[k2];
	}

	monolis_allreduce_R(
		num_modes,
		vec,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k2 = 0; k2 < num_modes; k2++){
        rhs[k2] = vec[k2];
	}


	double* matrix_buffer = BB_std_calloc_1d_double(matrix_buffer, num_modes * num_modes);
	
	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            matrix_buffer[k1 * num_modes + k2] = mat[k1][k2];
		}
	}

	monolis_allreduce_R(
		num_modes * num_modes,
		matrix_buffer,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            mat[k1][k2] = matrix_buffer[k1 * num_modes + k2];
		}
	}
}