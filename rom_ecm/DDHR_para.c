
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

void HROM_ddecm_memory_allocation_para_pre(
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
    hlpod_ddhr->matrix = BB_std_calloc_3d_double(hlpod_ddhr->matrix, hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
    hlpod_ddhr->RH = BB_std_calloc_2d_double(hlpod_ddhr->RH, hlpod_vals->n_neib_vec, num_subdomains);
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

void HROM_ddecm_set_num_modes_para_vis(
        BBFE_DATA*      fe,
        HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int       total_num_nodes,
        const int       num_subdomains,
        const char*     directory)
{
    const int nl = fe->local_num_nodes;

    double* DD_num_modes_elem;       // 各要素に持たせるモード本数
    double* DD_num_modes_node;       // 各節点に持たせるモード本数
    int*    DD_num_modes_count_node; // 節点平均用カウンタ

    DD_num_modes_elem = BB_std_calloc_1d_double(
            DD_num_modes_elem,
            fe->total_num_elems);

    DD_num_modes_node = BB_std_calloc_1d_double(
            DD_num_modes_node,
            total_num_nodes);

    DD_num_modes_count_node = (int*)calloc(total_num_nodes, sizeof(int));

    /*
     * 全要素について、
     * その要素が単一領域内にあるか、領域をまたぐかを判定する。
     *
     * 単一領域内:
     *     その領域のモード本数を入れる
     *
     * 領域をまたぐ:
     *     -1 を入れる
     */
    for(int e = 0; e < fe->total_num_elems; e++) {
        int index0 = fe->conn[e][0];
        int sd0 = hlpod_mat->subdomain_id_in_nodes[index0];

        int is_cross_subdomain = 0;

        for(int i = 1; i < nl; i++) {
            int index = fe->conn[e][i];
            int sd = hlpod_mat->subdomain_id_in_nodes[index];

            if(sd != sd0) {
                is_cross_subdomain = 1;
                break;
            }
        }

        if(is_cross_subdomain) {
            DD_num_modes_elem[e] = -1.0;
        }
        else {
            if(0 <= sd0 && sd0 < num_subdomains) {
                int num_modes_dd =
                    hlpod_ddhr->num_neib_modes_1stdd_sum[sd0 + 1]
                  - hlpod_ddhr->num_neib_modes_1stdd_sum[sd0];

                DD_num_modes_elem[e] = (double)num_modes_dd;
            }
            else {
                /*
                 * 想定外の領域IDの場合も -1 にしておく。
                 */
                DD_num_modes_elem[e] = -1.0;
            }
        }

        /*
         * 節点データにも要素値を集約する。
         * 複数要素に接続する節点では平均値にする。
         */
        for(int i = 0; i < nl; i++) {
            int index = fe->conn[e][i];

            DD_num_modes_node[index] += DD_num_modes_elem[e];
            DD_num_modes_count_node[index] += 1;
        }
    }

    /*
     * 節点値を平均化する。
     *
     * 注意:
     * 領域境界付近では -1 と正のモード本数が平均されるため、
     * -1 そのものではなく中間値になる場合がある。
     */
    for(int i = 0; i < total_num_nodes; i++) {
        if(DD_num_modes_count_node[i] > 0) {
            DD_num_modes_node[i] /=
                (double)DD_num_modes_count_node[i];
        }
    }

    const char* filename;
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE,
            "DDECM/%s",
            "mode_distribution.vtk");

    filename = monolis_get_global_output_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            "./",
            fname);

    FILE* fp;
    fp = ROM_BB_write_fopen(fp, filename, directory);

    switch(fe->local_num_nodes) {
        case 4:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
            break;

        case 8:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
            break;

        default:
            printf("ERROR: unsupported local_num_nodes = %d\n",
                    fe->local_num_nodes);
            break;
    }

    /*
     * 節点データとしてモード本数を出力
     */
    fprintf(fp, "POINT_DATA %d\n", total_num_nodes);

    BB_vtk_write_point_vals_scalar(
            fp,
            DD_num_modes_node,
            total_num_nodes,
            "DD-num_modes-node");

    BB_std_free_1d_double(DD_num_modes_elem, fe->total_num_elems);
    BB_std_free_1d_double(DD_num_modes_node, total_num_nodes);
    free(DD_num_modes_count_node);

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
            printf("val = %e\n", monolis->mat.R.B[index + i]);
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
            printf("val = %e\n", monolis->mat.R.B[index + i]);
		}
		index += hlpod_mat->num_modes_internal[k];		
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

}
