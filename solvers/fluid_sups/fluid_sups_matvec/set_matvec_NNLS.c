//todo: NNLS用の関数をまとめる

#include "rom_dataset.h"
#include "core_FOM.h"

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_RH_for_NNLS_para(
		BBFE_DATA*     	fe,
		VALUES*         vals,
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

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];

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

                int index = fe->conn[e][i];

				int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];
				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                for(int d=0; d<4; d++) {
                    integ_val[d] = BBFE_std_integ_calc(
                            np, val_ip[d], basis->integ_weight, Jacobian_ip);
                            
                    for(int k = IS; k < IE; k++){
                        hlpod_ddhr->matrix[ns*(hlpod_vals->n_neib_vec) + k][m][n] -= integ_val[d] * hlpod_mat->neib_vec[index*4 + d][k];
                        hlpod_ddhr->RH[ns*(hlpod_vals->n_neib_vec) + k][n] -= integ_val[d] * hlpod_mat->neib_vec[index*4 + d][k];

                        //hlpod_ddhr->matrix[ns*(hlpod_vals->n_neib_vec) + k][m][n] += integ_val[d] * hlpod_mat->neib_vec[index*4 + d][k];
                        //hlpod_ddhr->RH[ns*(hlpod_vals->n_neib_vec) + k][n] += integ_val[d] * hlpod_mat->neib_vec[index*4 + d][k]; 
                    }
                }

            }
        }
    }

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);

}

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_residuals_for_NNLS_para(
		BBFE_DATA*     	fe,
		VALUES*         vals,
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
    //int ns = index_snap*2;
    int ns = index_snap;

	double** local_matrix;  double* local_vec;
	local_matrix   = BB_std_calloc_2d_double(local_matrix, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    local_vec   = BB_std_calloc_1d_double(local_vec, hlpod_vals->n_neib_vec);

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

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];

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


                    int index_i = fe->conn[e][i];
                    int index_j = fe->conn[e][j];

                    int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
                    int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
                    int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

                    int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
                    int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
                    int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

                    for(int a=0; a<4; a++){
                        for(int b=0; b<4; b++) {
                            double integ_val = BBFE_std_integ_calc(
                                    np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                            if( bc->D_bc_exists[index_j*4+b]) {
                                for(int k1 = IS; k1 < IE; k1++){
                                    double val = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val * bc->imposed_D_val[index_j*4+b];

                                    hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] += val;
                                    hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] += val;

                                    //hlpod_ddhr->matrix[ns*hlpod_vals->n_neib_vec + k1][m][n] -= val;
                                    //hlpod_ddhr->RH[ns*hlpod_vals->n_neib_vec + k1][n] -= val;
                                }
                            }
                            else{
                                for(int k1 = IS; k1 < IE; k1++) {
                                    double A = hlpod_mat->neib_vec[index_i*4+a][k1] * integ_val;
                                    local_vec[k1] = 0.0;

                                    int index1 = 0;
                                    int index2 = 0;

                                    for(int ki = 0; ki < num_neib; ki++) {
                                        for(int kj = 0; kj < hlpod_mat->num_modes_1stdd_neib[ki]; kj++) {							
                                            double B = hlpod_mat->neib_vec[index_j*4+b][index2];
                                            double C = hlpod_mat->mode_coef_1stdd[index1 + kj];

                                            local_vec[k1] += A * B * C;
                                            index2++;
                                        }
                                        index1 += hlpod_mat->max_num_neib_modes[ki];
                                    }
                                }

                                for(int k1 = IS; k1 < IE; k1++) {
                                    int index = ns*hlpod_vals->n_neib_vec + k1;

                                    hlpod_ddhr->matrix[index][m][n] += local_vec[k1];
                                    hlpod_ddhr->RH[index][n] += local_vec[k1];
                                }
                            }
                        }
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


void ddhr_set_matvec_RH_for_NNLS_para_volume_const(
		BBFE_DATA*     	fe,
		VALUES*         vals,
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

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_elems[n]; m++) {
            int e = hlpod_ddhr->elem_id_local[m][n];

            BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

            double vol = BBFE_std_integ_calc_volume(
                    np, basis->integ_weight, Jacobian_ip);
            
            hlpod_ddhr->matrix[(hlpod_vals->n_neib_vec)*num_snapshot][m][n] += vol;
            hlpod_ddhr->RH[(hlpod_vals->n_neib_vec)*num_snapshot][n] += vol; 
        }
    }
    
	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}
