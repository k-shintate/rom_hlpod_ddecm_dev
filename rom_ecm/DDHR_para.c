
#include "DDHR_para.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";

void HROM_ddecm_memory_allocation_sol_vec(
	    HR_VALUES*      hr_vals,
        const int       total_num_nodes,
        const int       dof)
{
    hr_vals->sol_vec = BB_std_calloc_1d_double(hr_vals->sol_vec, total_num_nodes*dof);
}

void HROM_ddecm_memory_allocation_para_online(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
	    HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes)
{
    hlpod_ddhr->reduced_mat = BB_std_calloc_2d_double(hlpod_ddhr->reduced_mat, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    hlpod_ddhr->reduced_RH = BB_std_calloc_1d_double(hlpod_ddhr->reduced_RH, hlpod_vals->n_neib_vec);
}

void HROM_ddecm_memory_allocation_para(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains)
{
	int max_num_elem = 0;
	if(num_subdomains==1){
		max_num_elem = total_num_elem;
	}
	else{
		max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	}

    //for NNLS
    hlpod_ddhr->matrix = BB_std_calloc_3d_double(hlpod_ddhr->matrix, total_num_snapshot*hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
    hlpod_ddhr->RH = BB_std_calloc_2d_double(hlpod_ddhr->RH, total_num_snapshot*hlpod_vals->n_neib_vec, num_subdomains);
}

void HROM_ddecm_memory_free_para(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains)
{
	int max_num_elem = 0;
	if(num_subdomains==1){
		max_num_elem = total_num_elem;
	}
	else{
		max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	}

    BB_std_free_3d_double(hlpod_ddhr->matrix, 2*total_num_snapshot*hlpod_vals->n_neib_vec +1, max_num_elem, num_subdomains);
    BB_std_free_2d_double(hlpod_ddhr->RH, 2*total_num_snapshot*hlpod_vals->n_neib_vec +1, num_subdomains);
}

void HROM_ddecm_set_selected_elems_para_vis(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		total_num_nodes,
		const int		num_subdomains,
		const char*     directory)
{
	int nl = fe->local_num_nodes;
	const int myrank = monolis_mpi_get_global_my_rank();
	double* ECM_elem;		//非ゼロ要素可視化用ベクトル
	double* ECM_elem_weight;//非ゼロ要素可視化用重みベクトル
	double* ECM_wireframe;	//wireframe可視化用ゼロベクトル

	ECM_elem = BB_std_calloc_1d_double(ECM_elem, total_num_nodes);
	ECM_elem_weight = BB_std_calloc_1d_double(ECM_elem_weight, total_num_nodes);
	ECM_wireframe = BB_std_calloc_1d_double(ECM_wireframe, total_num_nodes);

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_selected_elems[n]; m++) {
			int e = hlpod_ddhr->id_selected_elems[m][n];
			for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight[m][n];
			}
		}

		for(int m=0; m<(hlpod_ddhr->num_selected_elems_D_bc[n]); m++) {
			int e = hlpod_ddhr->id_selected_elems_D_bc[m][n];
			for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight_D_bc[m][n];
			}
		}
	}


    for(int m = 0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems[m];
		for(int i=0; i<nl; i++) {
			int index = fe->conn[e][i];
			ECM_elem[index] = myrank + 1;
			ECM_elem_weight[index] += hlpod_ddhr->ovl_elem_weight[m];
		}
	}

    for(int m=0; m<(hlpod_ddhr->ovl_num_selected_elems_D_bc); m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];
		for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->ovl_elem_weight_D_bc[m];
		}
	}


	const char* filename;

	char fname[BUFFER_SIZE];

	snprintf(fname, BUFFER_SIZE,"DDECM/%s", OUTPUT_FILENAME_ECM_ELEM_VTK);
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname);

	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}
	
	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);

	BB_vtk_write_point_vals_scalar(fp, ECM_elem, fe->total_num_nodes, "ECM-elem");
	BB_vtk_write_point_vals_scalar(fp, ECM_elem_weight, fe->total_num_nodes, "ECM-elem_weight");
	BB_vtk_write_point_vals_scalar(fp, ECM_wireframe, fe->total_num_nodes, "ECM-wireframe");

	BB_std_free_1d_double(ECM_elem, fe->total_num_nodes);
	BB_std_free_1d_double(ECM_elem_weight, fe->total_num_nodes);
	BB_std_free_1d_double(ECM_wireframe, fe->total_num_nodes);

	fclose(fp);

}

void HROM_ddecm_to_monollis_rhs_para(
	MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
	HLPOD_MAT*    hlpod_mat,
	const int		num_2nddd,
	const int		max_num_bases)
{
	int sum = 0;
	int index = 0;

	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			monolis->mat.R.B[index + i] = hlpod_ddrh->reduced_RH[index + i];
		}
		index += hlpod_mat->num_modes_internal[k];		
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

}

void HROM_ddecm_calc_block_solution(
	MONOLIS_COM*	monolis_com,
	BBFE_DATA* 		fe,
    HR_VALUES*      hr_vals,
	HLPOD_MAT*      hlpod_mat,
	const int		num_2nddd,
    const int 		dof)
{
	for(int j = 0; j < fe->total_num_nodes * dof; j++){
		hr_vals->sol_vec[j] = 0.0;
	}

	int index_row = 0;
	int index_column = 0;
	int sum = 0;
	int index = 0;

	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
    				index_row = hlpod_mat->node_id[j + sum]* dof + l;
    				hr_vals->sol_vec[index_row] += hlpod_mat->pod_modes[index_row][index_column + i] * hlpod_mat->mode_coef[index + i];
                }
			}
		}
		index_column += hlpod_mat->num_modes_internal[k];
		index += hlpod_mat->num_modes_internal[k];
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}
}

void HROM_ddecm_add_monollis_rhs_para(
	MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
	HLPOD_MAT*    hlpod_mat,
	const int		num_2nddd,
	const int		max_num_bases)
{
	int sum = 0;
	int index = 0;

	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			monolis->mat.R.B[index + i] += hlpod_ddrh->reduced_RH[index + i];
            //printf("val = %e\n", monolis->mat.R.B[index + i]);
		}
		index += hlpod_mat->num_modes_internal[k];		
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

}
