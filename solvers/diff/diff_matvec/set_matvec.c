
//todo: set_mat, set_vec関数をまとめる

#include "set_matvec.h"

void ddhr_set_reduced_mat_para_debug(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
        HLPOD_VALUES*   hlpod_vals,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }
	
	for(int e = 0; e < fe->total_num_elems;e++){

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
					
					//基底本数ループ
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

					if( bc->D_bc_exists[index_j]){
					}
					else{
						for(int k1 = 0; k1 < hlpod_vals->num_modes; k1++){
							for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
								double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * hlpod_mat->neib_vec[index_j][k2];
								hlpod_ddhr->reduced_mat[k1][k2] += val;
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

void HROM_ecm_set_reduced_mat(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_HR*       hlpod_hr,
        const int num_modes,
		const double    dt)
{
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

    for(int i = 0;i < num_modes; i++){
        for(int j =0; j < num_modes; j++){
            hlpod_hr->reduced_mat[i][j] = 0.0;
        }
        hlpod_hr->reduced_RH[i] = 0.0;
    }

    for(int h=0; h<(hlpod_hr->num_selected_elems); h++) {
        int e = hlpod_hr->id_selected_elems[h];

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

				if( bc->D_bc_exists[index_j]) {
				}
				else{
                    for(int k1 = 0; k1 < num_modes; k1++){
                        for(int k2 = 0; k2 < num_modes; k2++){
                            double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2] ;

                            hlpod_hr->reduced_mat[k1][k2] += hlpod_hr->elem_weight[h] *  val;
                        }
                    }
				}
			}
		}
	}

	for(int h=0; h<(hlpod_hr->num_selected_elems_D_bc); h++) {
        int e = hlpod_hr->id_selected_elems_D_bc[h];

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

                if( bc->D_bc_exists[index_j]) {
                }
                else{
                    for(int k1 = 0; k1 < num_modes; k1++){
                        for(int k2 = 0; k2 < num_modes; k2++){
                            double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2] ;

                            hlpod_hr->reduced_mat[k1][k2] += hlpod_hr->elem_weight_D_bc[h] *  val;
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


void HROM_ecm_set_D_bc(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_HR*       hlpod_hr,
        const int num_modes,
		const double    dt)
{
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

    for(int h=0; h<(hlpod_hr->num_selected_elems_D_bc); h++) {
        int e = hlpod_hr->id_selected_elems_D_bc[h];

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
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
                    for(int k1 = 0; k1 < num_modes; k1++){
                        double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                        hlpod_hr->reduced_RH[ k1 ] += - hlpod_hr->elem_weight_D_bc[h] * val;
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

void HROM_ecm_set_reduced_vec(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_HR*       hlpod_hr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
        const double    dt,
		double       	t)
{
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

    for(int k = 0; k < num_modes; k++){
        hlpod_hr->reduced_RH[ k ] = 0.0;
    }

    for(int h=0; h<(hlpod_hr->num_selected_elems); h++) {
        int e = hlpod_hr->id_selected_elems[h];

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
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

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
			//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
                val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
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
                double val = integ_val * hlpod_mat->pod_modes[index][k];
                hlpod_hr->reduced_RH[ k ] += hlpod_hr->elem_weight[h] * val;
            }
            
		}
	}

    for(int h=0; h<(hlpod_hr->num_selected_elems_D_bc); h++) {
        int e = hlpod_hr->id_selected_elems_D_bc[h];

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
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

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
			//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
                val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
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
                double val = integ_val * hlpod_mat->pod_modes[index][k];
                hlpod_hr->reduced_RH[ k ] += hlpod_hr->elem_weight_D_bc[h] * val;
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

void HROM_ddecm_set_reduced_vec(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
{
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

    for(int k = 0; k < num_modes * num_subdomains; k++){
        hlpod_ddhr->reduced_RH[ k ] = 0.0;
    }

	//for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
			
			int e = hlpod_ddhr->ovl_id_selected_elems[m];

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
			BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

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
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
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
					double val = integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight[m] * val;
				}
			}
		}

		for(int m=0; m<(hlpod_ddhr->ovl_num_selected_elems_D_bc); m++) {
			int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

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
			BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

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
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
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
					double val = integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
				}
				
			}
		//}
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

void HROM_ddecm_set_reduced_mat(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int i = 0; i < num_modes*num_subdomains; i++){
        for(int j =0; j < num_modes*num_subdomains; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

	//for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
			int e = hlpod_ddhr->ovl_id_selected_elems[m];

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
					
					//基底本数ループ
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

					int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

					int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
					int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
					int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

					//printf("nummodes = %d\n", num_modes);

					if( bc->D_bc_exists[index_j]){
					}
					else{
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS; k2 < JE; k2++){
								double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2];
								hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight[m] *  val;
							}
						}
					}
				}
			}
		}
	//}

	//for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
			int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];
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
					
					//基底本数ループ
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

					int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

					int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
					int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
					int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

					if( bc->D_bc_exists[index_j]) {
					}
					else{
						/*
						for(int k1 = 0; k1 < num_modes*num_subdomains; k1++){
							for(int k2 = 0; k2 < num_modes*num_subdomains; k2++){
								double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2] ;
								hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
						}
						*/

						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS; k2 < JE; k2++){
								double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * hlpod_mat->pod_modes[index_j][k2] ;
								hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
						}
					}
				}
			}
		}
	//}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}

void HROM_ddecm_set_D_bc(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

        //printf("e = %d\n", e);

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

                int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index_i];

                if( bc->D_bc_exists[index_j]) {
                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];
                    for(int k1 = IS; k1 < IE; k1++){
                        double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                        hlpod_ddhr->reduced_RH[k1] += - hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
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

void ddhr_lb_set_reduced_mat_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
        HLPOD_VALUES*   hlpod_vals,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

	for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems[m];

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
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

				int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

				int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
				int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
				int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

				if( bc->D_bc_exists[index_j]){
				}
				else{
	                for(int k1 = IS; k1 < IE; k1++){
                        for(int k2 = JS; k2 < JE; k2++){
		                    double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * hlpod_mat->neib_vec[index_j][k2];
                            hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight[m] *  val;
                        }
                    }
				}	
    		    }
		}
	}

    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

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
                }
                else{
	                for(int k1 = IS; k1 < IE; k1++){
                        for(int k2 = JS; k2 < JE; k2++){
                            double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * hlpod_mat->neib_vec[index_j][k2] ;
                            hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
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


void HROM_ddecm_set_reduced_mat_para_save_memory(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*     hlpod_vals,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

	for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems[m];

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
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

				if( bc->D_bc_exists[index_j]){
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
							double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2];
                            hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight[m] *  val;
                        }
                    }
				}

				if(subdomain_id_j >= num_subdomains){
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += hlpod_ddhr->ovl_elem_weight[m] *  val;
							}
						}
				}

				}
			}
		}
	}


    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

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

                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
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
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
						}
					}

					if(subdomain_id_j >= num_subdomains){
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
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

void HROM_ddecm_set_D_bc_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
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

    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

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

		for(int j=0; j<nl; j++) {       //六面体一次要素は8
			int index_j = fe->conn[e][j];
			if( bc->D_bc_exists[index_j]) {

			for(int i=0; i<nl; i++) {
				
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

				int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index_i];
				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                for(int k1 = IS; k1 < IE; k1++){
                    double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                    hlpod_ddhr->reduced_RH[k1] += - hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
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


void HROM_ddecm_set_reduced_vec_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*    hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
{
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

    for(int k = 0; k < hlpod_vals->n_neib_vec; k++){
        hlpod_ddhr->reduced_RH[ k ] = 0.0;
    }


    for(int m = 0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems[m];

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
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
			int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

			if(subdomain_id < num_subdomains){
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
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

				for(int k = IS; k < IE; k++){
					double val = integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight[m] * val;
				}
			}     
		}
	}

    for(int m=0; m<(hlpod_ddhr->ovl_num_selected_elems_D_bc); m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];


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
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			int index = fe->conn[e][i];
			int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];
			if(subdomain_id < num_subdomains){

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
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
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
                    double val = integ_val * hlpod_mat->pod_modes[index][k];
                    hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
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
/********/




/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
		MONOLIS*		monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        const int 		num_modes,
		const double    dt)
{
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

	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
			hlpod_mat->VTKV[k1][k2] = 0.0;
		}
	}

	bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;
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
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

				set_element_reduced_mat(
					hlpod_mat->VTKV,
					hlpod_mat->pod_modes,
					fe->conn[e][i],
					fe->conn[e][j],
					integ_val,	
					num_modes);
					
			}
		}
	}

	double* vec;
	vec = BB_std_calloc_1d_double(vec, num_modes);

	printf("%d\n", num_modes);
	for(int k2 = 0; k2 < num_modes; k2++){
		vec[k2] = hlpod_mat->VTf[k2];
	}

	monolis_allreduce_R(
		num_modes,
		vec,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k2 = 0; k2 < num_modes; k2++){
		hlpod_mat->VTf[k2] = vec[k2];
	}


	double* matrix_buffer = BB_std_calloc_1d_double(matrix_buffer, num_modes * num_modes);
	
	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
			matrix_buffer[k1 * num_modes + k2] = hlpod_mat->VTKV[k1][k2];
		}
	}

	monolis_allreduce_R(
		num_modes * num_modes,
		matrix_buffer,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
			hlpod_mat->VTKV[k1][k2] = matrix_buffer[k1 * num_modes + k2];
		}
	}

	BB_std_free_1d_double(matrix_buffer, num_modes * num_modes);


	BB_std_free_1d_double(vec, num_modes);



	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	
	BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);
}


void set_D_bc_global_para(
		MONOLIS*     	monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        const int 		num_modes,
		const double    dt)
{
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

	
	for(int i = 0; i < num_modes; i++){
		hlpod_mat->VTf[i] = 0.0;
    }

	bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

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
                
				set_element_reduced_rhs_Dbc(
					bc,
					hlpod_mat->VTf,
					hlpod_mat->pod_modes,
					fe->conn[e][i],
					fe->conn[e][j],
					integ_val,
					num_modes);

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

	BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);
}



void set_reduced_vec_global_para(
		MONOLIS*     	monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
		HLPOD_VALUES*		hlpod_vals,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
        const double    dt,
		double       	t)
{
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
	
	bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

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
			//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
                val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], a_ip[p]);
/*
				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
*/
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			set_element_reduced_rhs(
				hlpod_mat->VTf,
				hlpod_mat->pod_modes,
				fe->conn[e][i],
				integ_val,
				num_modes);

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

	
	BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);
}
