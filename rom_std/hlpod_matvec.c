
#include "hlpod_matvec.h"
    
static const int BUFFER_SIZE = 10000;


void ROM_std_hlpod_set_snapmat_nobc(
    double*       	comp_vec,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int		dof,
    const int 		count)
{
    for(int i = 0; i < total_num_nodes * dof; i++){
        hlpod_mat->snapmat[i][count] = comp_vec[i];
    }
}


void ROM_std_hlpod_calc_reduced_mat_seq(
    MONOLIS*		monolis,
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int		num_base,
    const int		dof)
{
    int nl = total_num_nodes * dof;
    int k = num_base;

    hlpod_mat->KV = BB_std_calloc_2d_double(hlpod_mat->KV, nl, k);
    hlpod_mat->VTKV = BB_std_calloc_2d_double(hlpod_mat->VTKV, k, k);

    for(int i = 0; i < nl; i++){
        for(int j = 0; j < k; j++){
            hlpod_mat->KV[i][j] = 0.0;
        }
    }
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            hlpod_mat->VTKV[i][j] = 0.0;
        }
    }

    double* monolis_in;  double* monolis_out;
    monolis_in = BB_std_calloc_1d_double(monolis_in, nl);
    monolis_out = BB_std_calloc_1d_double(monolis_out, nl);

    for(int i = 0; i < k; i++){
        for(int j = 0; j < nl; j++){
            monolis_in[j] = hlpod_mat->pod_modes[j][i];
        }
        monolis_matvec_product_R(monolis, monolis_com, monolis_in, monolis_out);
        for(int j = 0; j < nl; j++){
            hlpod_mat->KV[j][i] = monolis_out[j];
        }
    }

    for(int i = 0; i < k; i++){
        for(int l = 0; l < k; l++){
            for(int j = 0; j < nl; j++){
                hlpod_mat->VTKV[i][l] += hlpod_mat->pod_modes[j][i] * hlpod_mat->KV[j][l];
            }
        }
    }
    
    BB_std_free_2d_double(hlpod_mat->KV, nl, k);
    BB_std_free_1d_double(monolis_in, nl);
    BB_std_free_1d_double(monolis_out, nl);

    hlpod_mat->VTf = BB_std_calloc_1d_double(hlpod_mat->VTf,k);
    hlpod_mat->mode_coef = BB_std_calloc_1d_double(hlpod_mat->mode_coef,k);
}



void ROM_std_hlpod_calc_reduced_mat_seq_block(
    MONOLIS*		monolis,
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int		num_base,
    const int		num_2nddd,
    const int 		dof)
{
    int nl = total_num_nodes * dof;
    int total_num_modes = num_base * num_2nddd;

    hlpod_mat->KV = BB_std_calloc_2d_double(hlpod_mat->KV, nl, total_num_modes);

    double* monolis_in;  double* monolis_out;
    monolis_in = BB_std_calloc_1d_double(monolis_in, nl);
    monolis_out = BB_std_calloc_1d_double(monolis_out, nl);

    for(int i = 0; i < total_num_modes; i++){
        for(int j = 0; j < nl; j++){
            monolis_in[j] = hlpod_mat->pod_modes[j][i];
        }
        monolis_matvec_product_R(monolis, monolis_com, monolis_in, monolis_out);
        for(int j = 0; j < nl; j++){
            hlpod_mat->KV[j][i] = monolis_out[j];
        }
    }

    BB_std_free_1d_double(monolis_in, nl);
    BB_std_free_1d_double(monolis_out, nl);

    hlpod_mat->VTKV = BB_std_calloc_2d_double(hlpod_mat->VTKV, num_base*num_2nddd, num_base*num_2nddd);

    /*localの行列積*/
    int n_neib_vec = num_base * num_2nddd;
    int globalIndexCol1 = 0;
    int globalIndexCol2 = 0;

    for(int k1 = 0; k1 < num_2nddd; k1++) {
        for(int i1 = 0; i1 < hlpod_mat->num_modes_internal[k1]; i1++) {
            int localIndexRow = 0;
            int sum = 0;
            int localIndexCol1 = 0;
            int localIndexCol2 = 0;
            
            for(int k = 0; k < num_2nddd; k++) {
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++) {
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++) {
                        for(int l = 0; l < dof; ++l) {
                            localIndexRow = hlpod_mat->node_id[j + sum] * dof + l;
                            hlpod_mat->VTKV[localIndexCol2 + i][globalIndexCol2 + i1] += 
                                hlpod_mat->pod_modes[localIndexRow][localIndexCol1 + i] * 
                                hlpod_mat->KV[localIndexRow][globalIndexCol1 + i1];
                        }
                    }
                }
                localIndexCol1 += hlpod_mat->num_modes_internal[k];
                localIndexCol2 += num_base;
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
        }
        globalIndexCol1 += hlpod_mat->num_modes_internal[k1];
        globalIndexCol2 += num_base;
    }
    /****************/
    BB_std_free_2d_double(hlpod_mat->KV, nl, total_num_modes);
    BB_std_free_2d_double(hlpod_mat->VTKV, total_num_modes, total_num_modes);
}

//省メモリver
void ROM_std_hlpod_calc_reduced_mat_save_memory(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    MONOLIS_COM*    mono_com0,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int       n_neib_vec,
    const int       num_2nddd,
    const int       num_modes,
    const int 		dof)
{
    const int NDOF  = total_num_nodes * dof;

    double* monolis_in;  double* monolis_out;  double* monolis_in2;
    monolis_in = BB_std_calloc_1d_double(monolis_in, NDOF);
    monolis_in2 = BB_std_calloc_1d_double(monolis_in2, NDOF);
    monolis_out = BB_std_calloc_1d_double(monolis_out, NDOF);

    hlpod_mat->VTKV = BB_std_calloc_2d_double(hlpod_mat->VTKV, num_modes, n_neib_vec);

	double t = monolis_get_time_global_sync();
	printf("test5");

    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
        if(l < num_modes){
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
            }

            for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
                monolis_in[j] = hlpod_mat->pod_modes[j][l];
            }

            monolis_matvec_product_R(monolis, mono_com0, monolis_in, monolis_out);

            int index_row = 0;
            int sum = 0;
            int index_column = 0;
            for(int k = 0; k < num_2nddd; k++){
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                        for(int m = 0; m < dof; m++){
                            index_row = hlpod_mat->node_id[j + sum] * dof + m;
                            hlpod_mat->VTKV[index_column + i][l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis_out[index_row];
                        }
                    }
                }
                index_column += hlpod_mat->num_modes_internal[k];
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
        }
	else{
	}
        if(l < num_modes){
            monolis_mpi_update_R( monolis_com, total_num_nodes, dof, monolis_in);
        }
        else{
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
            }
            monolis_mpi_update_R( monolis_com, total_num_nodes, dof, monolis_in);
        }

//        if(l < num_modes){
            int index_column2 = num_modes;
            for(int n = 0; n <  monolis_com->recv_n_neib; n++){
                int iS =  monolis_com->recv_index[n];
                int iE =  monolis_com->recv_index[n + 1];

                for(int k = 0; k < NDOF; k++){
                    monolis_in2[k] = 0;
                }

                for(int k = iS; k < iE; k++){
                    for(int m = 0; m < dof; m++){
                        int index =  monolis_com->recv_item[k] * dof + m;
                        monolis_in2[index] = monolis_in[index];
                    }
                }

                monolis_matvec_product_R(monolis,  mono_com0, monolis_in2, monolis_out);

                int sum = 0;
                int index_row = 0;
                int index_column = 0;

                for(int k = 0; k < num_2nddd; k++){
                    for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                        for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                            for(int m = 0; m < dof; m++){

                                index_row = hlpod_mat->node_id[j + sum] * dof + m;
                                if(l < hlpod_mat->num_modes_2nddd[n + 1]){
                                    hlpod_mat->VTKV[index_column + i][index_column2 + l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis_out[index_row];
                                }
                            }
                        }
                    }
                    index_column += hlpod_mat->num_modes_internal[k];
                    sum += hlpod_mat->n_internal_vertex_subd[k];
                }
                index_column2 += hlpod_mat->num_modes_2nddd[n + 1];
            }
        //}

    }

    BB_std_free_1d_double(monolis_in, NDOF);
    BB_std_free_1d_double(monolis_out, NDOF);
    BB_std_free_1d_double(monolis_in2, NDOF);
}


void ROM_std_hlpod_set_reduced_mat(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		total_num_nodes,
    const int		num_base,
    const int		num_2nddd,
    const int 		dof)
{
    int nl = total_num_nodes * dof;
    int total_num_modes = num_base * num_2nddd;
    int index1 = 0;
    int index2 = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int m = 0; m < hlpod_mat->num_modes_internal[k]; m++){
            for(int n = 0; n < hlpod_mat->num_modes_internal[k]; n++){
                double val = hlpod_mat->VTKV[k * num_base + m][k * num_base + n];

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis,
                    k,
                    k,
                    m,
                    n,
                    val);
            }
        }
    }
    
    for(int k = 0; k < num_2nddd; k++){
    int iS = hlpod_meta->index[k];
    int iE = hlpod_meta->index[k + 1];

    for(int i = iS; i < iE; i++){
            for(int m = 0; m < hlpod_mat->num_modes_internal[k]; m++){
                for(int n = 0; n < hlpod_mat->num_modes_internal[hlpod_meta->item[i]]; n++){
                    double val = hlpod_mat->VTKV[k * num_base + m][hlpod_meta->item[i] * num_base + n];

                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis,
                        k,
                        hlpod_meta->item[i],
                        m,
                        n,
                        val);
                }
            }
        }
    }

    hlpod_mat->VTf = BB_std_calloc_1d_double(hlpod_mat->VTf, total_num_modes);
    hlpod_mat->mode_coef = BB_std_calloc_1d_double(hlpod_mat->mode_coef, total_num_modes);
}


/*
void ROM_std_hlpod_set_reduced_mat_para(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		max_num_bases,
    const int		total_num_bases,
    const int		num_2nddd)
{
    const int M = max_num_bases;
    const int n_neib_vec = hlpod_vals->n_neib_vec;

    double** mat;
    mat = BB_std_calloc_2d_double(mat, M , M );

    int index1 = 0;
    int index2 = 0;
    int num_modes = 0;

double t7 = monolis_get_time_global_sync();
        printf("test7");


    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        int iS = hlpod_mat->num_modes_1stdd[k];
        int iE = hlpod_mat->num_modes_1stdd[k+1];

        index1 = 0; 
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                mat[index1][index2] = hlpod_mat->VTKV[m][n];
            
                monolis_add_scalar_to_sparse_matrix_R(
                    monolis,
                    k,
                    k,
                    index1,
                    index2,
                    mat[index1][index2]);

                index2++;
            }
            index1++;
        }
    }
    
double t8 = monolis_get_time_global_sync();
        printf("test8");


    for(int k = 0; k < monolis_com->n_internal_vertex; k++){

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){

                if (hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){
                    
                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k+1];

                    index1 = 0;
                    for(int m = IS; m < IE; m++){					

                        int IIS = hlpod_mat->num_modes_1stdd[j];
                        int IIE = hlpod_mat->num_modes_1stdd[j + 1];
                        num_modes += IIE - IIS;
                        
                        index2 = 0;

                        for(int n = IIS; n < IIE; n++){		
                            mat[index1][index2] = hlpod_mat->VTKV[m][n];
                            
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                k,
                                hlpod_meta->item[i],
                                index1,
                                index2,
                                mat[index1][index2]);
                            
                            index2++;
                        }
                        index1++;

                    }

                }

            }

        }

    }

    //hlpod_mat->mode_coef = BB_std_calloc_1d_double(hlpod_mat->mode_coef, num_modes);
    //hlpod_mat->VTf = BB_std_calloc_1d_double(hlpod_mat->VTf, num_modes);
   

    //BB_std_free_2d_double(hlpod_mat->VTKV, total_num_bases , n_neib_vec);
    BB_std_free_2d_double(mat, M, M);

double t6 = monolis_get_time_global_sync();
        printf("test6");


}
*/


void ROM_std_hlpod_set_reduced_mat_para(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*     hlpod_meta,
    const int       max_num_bases,
    const int       total_num_bases,
    const int       num_2nddd)
{
    const int M = max_num_bases;
    const int n_neib_vec = hlpod_vals->n_neib_vec;
    int rank = monolis_mpi_get_global_my_rank();

    double** mat;
    mat = BB_std_calloc_2d_double(mat, M, M);

    int index1 = 0;
    int index2 = 0;
    int num_modes = 0;

    /* 対角ブロック */
    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        int iS = hlpod_mat->num_modes_1stdd[k];
        int iE = hlpod_mat->num_modes_1stdd[k + 1];

        double frob_sq = 0.0;

        index1 = 0;
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                double val = hlpod_mat->VTKV[m][n];
                mat[index1][index2] = val;

                frob_sq += val * val;

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis,
                    k,
                    k,
                    index1,
                    index2,
                    val);

                index2++;
            }
            index1++;
        }

        //printf("[diag] block (%d,%d) Frobenius norm = %e\n", rank,  k, k, sqrt(frob_sq));
	
/*
        if(max_num_bases > iE - iS){
		for(int m = iE - iS; m < max_num_bases; m++){
                	monolis_add_scalar_to_sparse_matrix_R(
                    		monolis,
                    		k,
                    		k,
                    		index1,
                    		index1,
                    		1);
			index1++;
		}
        }
*/	
	
	
    }

    /* 非対角ブロック */
    for(int k = 0; k < monolis_com->n_internal_vertex; k++){

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){

                if(hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){

                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k + 1];

                    double frob_sq = 0.0;

                    index1 = 0;
                    for(int m = IS; m < IE; m++){

                        int IIS = hlpod_mat->num_modes_1stdd[j];
                        int IIE = hlpod_mat->num_modes_1stdd[j + 1];
                        num_modes += IIE - IIS;

                        index2 = 0;
                        for(int n = IIS; n < IIE; n++){
                            double val = hlpod_mat->VTKV[m][n];
                            mat[index1][index2] = val;

                            frob_sq += val * val;

                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                k,
                                hlpod_meta->item[i],
                                index1,
                                index2,
                                val);

                            index2++;
                        }
                        index1++;
                    }

                    //printf("rank = %d [offdiag] block (%d,%d) Frobenius norm = %e\n", rank,
                    //    k, hlpod_meta->item[i], sqrt(frob_sq));
                }
            }
        }
    }

    BB_std_free_2d_double(hlpod_mat->VTKV, total_num_bases, n_neib_vec);
    BB_std_free_2d_double(mat, M, M);
}


void ROM_std_hlpod_calloc_mode_coef_rhs(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		max_num_bases,
    const int		total_num_bases,
    const int		num_2nddd)
{
    const int M = max_num_bases;
    const int n_neib_vec = hlpod_vals->n_neib_vec;

    //double** mat;
    //mat = BB_std_calloc_2d_double(mat, M , M );

    int index1 = 0;
    int index2 = 0;
    int num_modes = 0;

    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        int iS = hlpod_mat->num_modes_1stdd[k];
        int iE = hlpod_mat->num_modes_1stdd[k+1];

        index1 = 0; 
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                //mat[index1][index2] = hlpod_mat->VTKV[m][n];
            /*
                monolis_add_scalar_to_sparse_matrix_R(
                    monolis,
                    k,
                    k,
                    index1,
                    index2,
                    mat[index1][index2]);
                index2++;
                */
            }
            index1++;
        }
    }
    

    for(int k = 0; k < monolis_com->n_internal_vertex; k++){

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){

                if (hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){
                    
                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k+1];

                    index1 = 0;
                    for(int m = IS; m < IE; m++){					

                        int IIS = hlpod_mat->num_modes_1stdd[j];
                        int IIE = hlpod_mat->num_modes_1stdd[j + 1];
                        num_modes += IIE - IIS;
                        
                        index2 = 0;

                        for(int n = IIS; n < IIE; n++){		
                            /*
                            mat[index1][index2] = hlpod_mat->VTKV[m][n];
                            
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                k,
                                hlpod_meta->item[i],
                                index1,
                                index2,
                                mat[index1][index2]);
*/                            
                            index2++;
                        }
                        index1++;

                    }

                }

            }

        }

    }

printf("num_modes= %d\n", num_modes);

	double tt = monolis_get_time_global_sync();

    hlpod_mat->mode_coef = BB_std_calloc_1d_double(hlpod_mat->mode_coef, num_modes);
    hlpod_mat->mode_coef_old = BB_std_calloc_1d_double(hlpod_mat->mode_coef_old, num_modes);
    hlpod_mat->VTf = BB_std_calloc_1d_double(hlpod_mat->VTf, num_modes);
    hlpod_mat->VTf_tmp = BB_std_calloc_1d_double(hlpod_mat->VTf_tmp, num_modes);

}


void ROM_std_hlpod_reduced_rhs_to_monollis(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int       num_2nd_subdomains,
    const int		num_modes)
{
    int index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            monolis->mat.R.B[index + i] = hlpod_mat->VTf[index + i];
            if(monolis_mpi_get_global_my_rank()==0){
                //printf("VTf[%d] = %e\n", index + i, hlpod_mat->VTf[index + i]);
            }
        }
        index += hlpod_mat->num_modes_internal[k];
    }

    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf[index + i] = 0.0;
        }
        index += hlpod_mat->num_modes_internal[k];
    }
    
}


void ROM_std_hlpod_calc_reduced_rhs(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		max_num_bases,
    const int		num_2nddd,
    const int 		dof)
{
    int index = 0;
    int index_row = 0;
    int sum = 0;
    int index_column = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf[index + i] = 0.0;
            //printf("%lf ", hlpod_mat->VTf[index + i]);
        }
        index += hlpod_mat->num_modes_internal[k];
    }
    
    index_row = 0;
    sum = 0;
    index_column = 0;
    index = 0;

    for(int k = 0; k < num_2nddd; k++){

        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
                    index_row = hlpod_mat->node_id[j + sum] * dof + l;
                    hlpod_mat->VTf[index + i] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis->mat.R.B[index_row];
                }
            }
        }
        index_column += hlpod_mat->num_modes_internal[k];
        index += hlpod_mat->num_modes_internal[k];
        sum += hlpod_mat->n_internal_vertex_subd[k];

    }

}

void ROM_std_hlpod_calc_reduced_rhs_add(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		max_num_bases,
    const int		num_2nddd,
    const int 		dof)
{
    int index = 0;
    int index_row = 0;
    int sum = 0;
    int index_column = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf[index + i] += hlpod_mat->VTf_tmp[index + i];
            //printf("%e ", hlpod_mat->VTf_tmp[index + i]);
        }
        index += hlpod_mat->num_modes_internal[k];
    }

}

void ROM_std_hlpod_calc_sol(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int		num_modes,
    const int		num_2nddd,
    const int 		dof)
{
    for(int j = 0; j < total_num_nodes * dof; j++){
        hlpod_vals->sol_vec[j] = 0.0;
    }

    int index_row = 0;
    int index_column = 0;
    int sum = 0;
    int index = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
                    index_row = hlpod_mat->node_id[j + sum] * dof + l;
                    hlpod_vals->sol_vec[index_row] += hlpod_mat->pod_modes[index_row][index_column + i] * hlpod_mat->mode_coef[index + i];
                }
            }
        }
        index_column += hlpod_mat->num_modes_internal[k];
        index += hlpod_mat->num_modes_internal[k];
        sum += hlpod_mat->n_internal_vertex_subd[k];
    }

}


void ROM_std_hlpod_calc_sol_global_para(
    double*         ansvec,
    double**        pod_modes,
    double*         mode_coef,
    const int 		total_num_nodes,
    const int 		num_base,
    const int		dof)
{
    int nl = total_num_nodes * dof;
    int k = num_base;

    for(int j = 0; j < nl; j++){
        ansvec[j] = 0.0;
    }

    for(int i = 0; i < k; i++){
        for(int j = 0; j < nl; j++){
            ansvec[j] += pod_modes[j][i] * mode_coef[i];
        }
    }

}

/*
void ROM_std_hlpod_update_global_modes(
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*		hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int		ndof)
{
    double* vec;
    vec = BB_std_calloc_1d_double(vec, total_num_nodes*ndof);

    for(int j = 0; j < num_modes; j++){
        for(int i = 0; i < n_internal_vertex * ndof; i++){
            vec[i] = hlpod_mat->pod_modes[i][j];
        }

        monolis_mpi_update_R(monolis_com, total_num_nodes, 4, vec);

        for(int i = 0; i < total_num_nodes*4; i++){
            hlpod_mat->pod_modes[i][j] = vec[i];
        }
    }

    BB_std_free_1d_double(vec, total_num_nodes*ndof);
}
*/

void ROM_std_hlpod_update_global_modes(
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*		hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int		ndof)
{
    double* vec;
    vec = BB_std_calloc_1d_double(vec, total_num_nodes*ndof);

    for(int j = 0; j < num_modes; j++){
        for(int i = 0; i < n_internal_vertex * ndof; i++){
            vec[i] = hlpod_mat->pod_modes[i][j];
        }

        monolis_mpi_update_R(monolis_com, total_num_nodes, 1, vec);

        for(int i = 0; i < total_num_nodes; i++){
            hlpod_mat->pod_modes[i][j] = vec[i];
        }
    }

    BB_std_free_1d_double(vec, total_num_nodes*ndof);
}
