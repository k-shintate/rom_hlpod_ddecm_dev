//todo: NNLS用の関数をまとめる

#include "rom_dataset.h"
#include "core_FOM.h"

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_RH_for_NNLS_para(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    printf("\n\nindex_snap = %d, num_modes = %d, num_subdomains = %d\n\n", index_snap, hlpod_vals->n_neib_vec, num_subdomains);

    int ns = index_snap;
    int nm = num_modes;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {

			int e = hlpod_ddhr->elem_id_local[m][n];
			
			BBFE_elemmat_set_Jacobian_array(
					Jacobian_ip,
					basis->num_integ_points,
					e,
					fe);

			double vol = BBFE_std_integ_calc_volume(
					basis->num_integ_points,
					basis->integ_weight,
					Jacobian_ip);

			double h_e = cbrt(vol);

			BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
			BBFE_elemmat_set_local_array_scalar(local_T, fe, hlpod_vals->sol_vec, e);

			for(int p=0; p<np; p++) {
				BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
				manusol_get_conv_vel(v_ip[p], x_ip[p]);
				a_ip[p] = manusol_get_mass_coef(x_ip[p]);
				k_ip[p] = manusol_get_diff_coef(x_ip[p]);
				T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
				f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
			}

			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
					val_ip[p] += BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i], f_ip[p]);
	/*
					val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
	*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/(dt) * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int index = fe->conn[e][i];

				int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];
				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

				for(int k = IS; k < IE; k++){
					hlpod_ddhr->matrix[ns*(hlpod_vals->n_neib_vec) + k][m][n] -= integ_val * hlpod_mat->neib_vec[index][k];
					hlpod_ddhr->RH[ns*(hlpod_vals->n_neib_vec) + k][n] -= integ_val * hlpod_mat->neib_vec[index][k]; 
				}

			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_residuals_for_NNLS_para(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        HLPOD_MAT*    hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_neib,			//1 + monolis_com->recv_n_neib
        const double    dt,
		double       	t)
{
    int ns = index_snap;

	double** local_matrix;  double* local_vec;
	local_matrix   = BB_std_calloc_2d_double(local_matrix, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    local_vec   = BB_std_calloc_1d_double(local_vec, hlpod_vals->n_neib_vec);

 	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);


	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {

			int e = hlpod_ddhr->elem_id_local[m][n];

			BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

			double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
			double h_e = cbrt(vol);

			BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

			for(int p=0; p<np; p++) {
				BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
				a_ip[p] = manusol_get_mass_coef(x_ip[p]);
				manusol_get_conv_vel(v_ip[p], x_ip[p]);
				k_ip[p] = manusol_get_diff_coef(x_ip[p]);

			}

			for(int i=0; i<nl; i++) {       //六面体一次要素は8
				for(int j=0; j<nl; j++) {

					for(int p=0; p<np; p++) {
						val_ip[p] = 0.0;
	/*
						double tau = BBFE_elemmat_convdiff_stab_coef_ns(
								k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
	/*
						val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
								basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
	*/
						val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
								fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
	/*
						val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
								fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
	*/
						//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
						val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
									basis->N[p][i], basis->N[p][j], a_ip[p]);
	/*
						val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
									fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
	*/
					}

					double integ_val = BBFE_std_integ_calc(
							np, val_ip, basis->integ_weight, Jacobian_ip);


					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

					int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

					int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
					int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
					int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

					if( bc->D_bc_exists[index_j]) {
						for(int k1 = IS; k1 < IE; k1++){
							double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * bc->imposed_D_val[index_j];

							hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] += val;
							hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] += val;
						}
					}

					else{

						for(int k1 = IS; k1 < IE; k1++) {
							double A = hlpod_mat->neib_vec[index_i][k1] * integ_val;
							local_vec[k1] = 0.0;

							int index1 = 0;
							int index2 = 0;

							for(int ki = 0; ki < num_neib; ki++) {
								for(int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; kj++) {							
									double B = hlpod_mat->neib_vec[index_j][index2];
									double C = hlpod_mat->mode_coef_1stdd[index1 + kj];

									local_vec[k1] += A * B * C;
									index2++;
								}
								index1 += hlpod_mat->max_num_neib_modes[ki];
							}
						}

						for(int k1 = IS; k1 < IE; k1++) {
							int index = ns * hlpod_vals->n_neib_vec + k1;

							hlpod_ddhr->matrix[index][m][n] += local_vec[k1];
							hlpod_ddhr->RH[index][n] += local_vec[k1];
						}

						
						
					}

				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_RH_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    int ns = index_snap;

	double** local_matrix;  double* local_vec;
	local_matrix   = BB_std_calloc_2d_double(local_matrix, hlpod_vals->num_modes, hlpod_vals->num_modes);
    local_vec   = BB_std_calloc_1d_double(local_vec, hlpod_vals->num_modes);

 	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
			int e = hlpod_ddhr->elem_id_local[m][n];
			//	for(int e=0; e<(fe->total_num_elems); e++) {
			BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

			double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
			double h_e = cbrt(vol);

			BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

			for(int p=0; p<np; p++) {
				BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
				a_ip[p] = manusol_get_mass_coef(x_ip[p]);
				manusol_get_conv_vel(v_ip[p], x_ip[p]);
				k_ip[p] = manusol_get_diff_coef(x_ip[p]);

			}

			for(int i=0; i<nl; i++) {       //六面体一次要素は8
				for(int j=0; j<nl; j++) {

					for(int p=0; p<np; p++) {
						val_ip[p] = 0.0;
	/*
						double tau = BBFE_elemmat_convdiff_stab_coef_ns(
								k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
	/*
						val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
								basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
	*/
						val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
								fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
	/*
						val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
								fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
	*/
						//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
						val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
									basis->N[p][i], basis->N[p][j], a_ip[p]);
	/*
						val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
									fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
	*/
					}

					double integ_val = BBFE_std_integ_calc(
							np, val_ip, basis->integ_weight, Jacobian_ip);

					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

                    int index_modes = 0;
                    int index_modes2 = 0;

					//Dirichlet境界条件を含むかを判別    
					if( bc->D_bc_exists[index_j]) {
                        for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
                                double val = hlpod_mat->pod_modes[index_i][k1 + index_modes] * integ_val * bc->imposed_D_val[index_j];
                                
                                hlpod_ddhr->matrix[ns*hlpod_vals->num_modes + k1 + index_modes][m][n] += val;
                                hlpod_ddhr->RH[ns*hlpod_vals->num_modes + k1 + index_modes][n] += val;
                            }
                            index_modes += hlpod_mat->num_modes_internal[m];
                        }
					}
					else{
                        index_modes = 0;
					    for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
                                local_vec[k1 + index_modes] = 0.0;
                            }
                            index_modes += hlpod_mat->num_modes_internal[m];
                        }

                        index_modes = 0;
                        index_modes2 = 0;
						for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
                     			for(int m2 = 0; m2 < num_subdomains; m2++){
                                    for(int k2 = 0; k2 < hlpod_mat->num_modes_internal[m2]; k2++){
        								local_matrix[k1 + index_modes][k2 + index_modes2] = 0.0;
                                    }
                                    index_modes2 += hlpod_mat->num_modes_internal[m2];
                                }
                                index_modes2 = 0;
							}
                            index_modes += hlpod_mat->num_modes_internal[m];
						}

                        index_modes = 0;
                        index_modes2 = 0;
						for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
                     			for(int m2 = 0; m2 < num_subdomains; m2++){
                                    for(int k2 = 0; k2 < hlpod_mat->num_modes_internal[m2]; k2++){
								        local_matrix[k1 + index_modes][k2 + index_modes2] += hlpod_mat->pod_modes[index_i][k1 + index_modes] * integ_val * hlpod_mat->pod_modes[index_j][k2 + index_modes2];
                                    }
                                    index_modes2 += hlpod_mat->num_modes_internal[m2];
                                }
                                index_modes2 = 0;
							}
                            index_modes += hlpod_mat->num_modes_internal[m];
						}

                        index_modes = 0;
                        index_modes2 = 0;
						for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
                     			for(int m2 = 0; m2 < num_subdomains; m2++){
                                    for(int k2 = 0; k2 < hlpod_mat->num_modes_internal[m2]; k2++){
								        local_vec[k1 + index_modes] += local_matrix[k1 + index_modes][k2 + index_modes2]* hlpod_mat->mode_coef[k2 + index_modes2];
                                    }
                                    index_modes2 += hlpod_mat->num_modes_internal[m2];
                                }
                                index_modes2 = 0;
							}
                            index_modes += hlpod_mat->num_modes_internal[m];
						}

                        index_modes = 0;
						for(int m = 0; m < num_subdomains; m++){
                            for(int k1 = 0; k1 < hlpod_mat->num_modes_internal[m]; k1++){
							    hlpod_ddhr->matrix[ns*hlpod_vals->num_modes + k1 + index_modes][m][n] += local_vec[k1 + index_modes];
							    hlpod_ddhr->RH[ns*hlpod_vals->num_modes + k1 + index_modes][n] += local_vec[k1 + index_modes];
						    }
                            index_modes += hlpod_mat->num_modes_internal[m];
                        }
					}

				}
			}
		}
	}

	BB_std_free_2d_double(local_matrix, hlpod_vals->num_modes, hlpod_vals->num_modes);
    BB_std_free_1d_double(local_vec, hlpod_vals->num_modes);
	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_residuals_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    printf("\n\nindex_snap = %d, num_modes = %d, num_subdomains = %d\n\n", index_snap, num_modes, num_subdomains);

    int ns = index_snap;
    int nm = num_modes;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
			int e = hlpod_ddhr->elem_id_local[m][n];

			BBFE_elemmat_set_Jacobian_array(
					Jacobian_ip,
					basis->num_integ_points,
					e,
					fe);

			double vol = BBFE_std_integ_calc_volume(
					basis->num_integ_points,
					basis->integ_weight,
					Jacobian_ip);

			double h_e = cbrt(vol);

			BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
			BBFE_elemmat_set_local_array_scalar(local_T, fe, hlpod_vals->sol_vec, e);

			for(int p=0; p<np; p++) {
				BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
				manusol_get_conv_vel(v_ip[p], x_ip[p]);
				a_ip[p] = manusol_get_mass_coef(x_ip[p]);
				k_ip[p] = manusol_get_diff_coef(x_ip[p]);
				T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
				f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
			}

			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
					val_ip[p] += BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i], f_ip[p]);
	/*
					val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
	*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/(dt) * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int index = fe->conn[e][i];
				for(int k = 0; k < num_modes*num_subdomains; k++){
					hlpod_ddhr->matrix[ns*(num_modes*num_subdomains) + k][m][n] -= integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->RH[ns*(num_modes*num_subdomains) + k][n] -= integ_val * hlpod_mat->pod_modes[index][k]; 
				}
				
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}

void ddhr_set_matvec_RH_for_NNLS2(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    printf("\n\nindex_snap = %d, num_modes = %d, num_subdomains = %d\n\n", index_snap, num_modes, num_subdomains);

    int ns = index_snap;
    int nm = num_modes;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
			int e = hlpod_ddhr->elem_id_local[m][n];

			BBFE_elemmat_set_Jacobian_array(
					Jacobian_ip,
					basis->num_integ_points,
					e,
					fe);

			double vol = BBFE_std_integ_calc_volume(
					basis->num_integ_points,
					basis->integ_weight,
					Jacobian_ip);

			double h_e = cbrt(vol);

			BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
			BBFE_elemmat_set_local_array_scalar(local_T, fe, hlpod_vals->sol_vec, e);

			for(int p=0; p<np; p++) {
				BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
				manusol_get_conv_vel(v_ip[p], x_ip[p]);
				a_ip[p] = manusol_get_mass_coef(x_ip[p]);
				k_ip[p] = manusol_get_diff_coef(x_ip[p]);
				T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
				f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
			}

			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
					val_ip[p] += BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i], f_ip[p]);
	/*
					val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
	*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/(dt) * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int index = fe->conn[e][i];
				for(int k = 0; k < num_modes*num_subdomains; k++){
					hlpod_ddhr->matrix[ns*num_modes*num_subdomains + k][m][n] += integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->RH[ns*num_modes*num_subdomains + k][n] += integ_val * hlpod_mat->pod_modes[index][k];

					hlpod_ddhr->matrix[ns*(num_modes*num_subdomains) + k + num_snapshot*(num_modes*num_subdomains)][m][n] -= integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->RH[ns*(num_modes*num_subdomains) + k + num_snapshot*(num_modes*num_subdomains)][n] -= integ_val * hlpod_mat->pod_modes[index][k]; 
				}
				
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}

void hr_set_matvec_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    printf("\n\nindex_snap = %d, num_modes = %d\n\n", index_snap, num_modes);

    int ns = index_snap;
    int nm = num_modes;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hlpod_vals->sol_vec, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;
/*
				double tau = BBFE_elemmat_convdiff_stab_coef_ns(
						k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);
/*
				val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
*/
				//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
                val_ip[p] += 1.0/(dt) * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], a_ip[p]);
/*
				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
*/
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

            int index = fe->conn[e][i];
            for(int k = 0; k < num_modes; k++){
                hlpod_hr->matrix[ ns * num_modes + k][e] += integ_val * hlpod_mat->pod_modes[index][k];
                hlpod_hr->RH[ns * num_modes + k] += integ_val * hlpod_mat->pod_modes[index][k];
            }
            
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}


void HROM_ecm_set_RH_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    printf("\n\nindex_snap = %d, num_modes = %d\n\n", index_snap, num_modes);

    int ns = index_snap;
    int nm = num_modes;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hlpod_vals->sol_vec, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;
/*
				double tau = BBFE_elemmat_convdiff_stab_coef_ns(
						k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);
/*
				val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
*/
				//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
                val_ip[p] += 1.0/(dt) * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], a_ip[p]);
/*
				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
*/
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

            int index = fe->conn[e][i];
            for(int k = 0; k < num_modes; k++){
                hlpod_hr->matrix[ ns * num_modes + k][e] -= integ_val * hlpod_mat->pod_modes[index][k];
                hlpod_hr->RH[ns * num_modes + k] -= integ_val * hlpod_mat->pod_modes[index][k];
            }
            
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}


void HROM_ecm_set_residuals_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t)
{
    int ns = index_snap;

	double** local_matrix;  double* local_vec;
	local_matrix   = BB_std_calloc_2d_double(local_matrix, num_modes, num_modes);
    local_vec   = BB_std_calloc_1d_double(local_vec, num_modes);

 	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

		}

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
*/
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
                    val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);
/*
					val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                //Dirichlet境界条件を含むかを判別    
                if( bc->D_bc_exists[index_j]) {
                    for(int k1 = 0; k1 < num_modes; k1++){
                        double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                        hlpod_hr->matrix[ ns * num_modes + k1][e] += val;
                        hlpod_hr->RH[ns * num_modes + k1] += val;
                    }
                }
                else{
                    for(int k1 = 0; k1 < num_modes; k1++){
                        for(int k2 = 0; k2 < num_modes; k2++){
                            local_matrix[k1][k2] += hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2];
                        }
                    }

                    for(int k1 = 0; k1 < num_modes; k1++){            
                        for(int k2 = 0; k2 < num_modes; k2++){
                            local_vec[k1] += local_matrix[k1][k2]* hlpod_mat->mode_coef[k2];
                        }
                    }
                    for(int k1 = 0; k1 < num_modes; k1++){   
                        hlpod_hr->matrix[ns * num_modes + k1][e] += local_vec[k1];
                        hlpod_hr->RH[ns * num_modes + k1] += local_vec[k1];
                    }
                    for(int k2 = 0; k2 < num_modes; k2++){
                        local_vec[k2] = 0.0;
                    }
                    for(int k1 = 0; k1 < num_modes; k1++){
                        for(int k2 = 0; k2 < num_modes; k2++){
                            local_matrix[k1][k2] = 0.0;
                        }
                    }
                }

			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}
