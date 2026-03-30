
#include "HR.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";


void hr_manusol_set_bc(
		BBFE_DATA* 	fe,
		BBFE_BC*   	bc,
		HLPOD_HR*   hlpod_hr,
		double     	t,
		double      (*func)(double, double, double, double) // scalar function(x, y, z, t)
		)
{
    for(int i=0; i < fe->total_num_nodes; i++) {
		bc->imposed_D_val[i] = 0.0;
    }

    for(int i = 0; i < hlpod_hr->num_selected_nodes_D_bc; i++) {
        int index = hlpod_hr->D_bc_node_id[i];
		bc->imposed_D_val[index] = func(fe->x[index][0], fe->x[index][1], fe->x[index][2], t);
    }
}

void HROM_ecm_set_element(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        HLPOD_HR*       hlpod_hr)
{
    hlpod_hr->matrix = BB_std_calloc_2d_double(hlpod_hr->matrix, total_num_snapshot * total_num_modes, total_num_elem);
    hlpod_hr->RH = BB_std_calloc_1d_double(hlpod_hr->RH, total_num_snapshot * total_num_modes);

    hlpod_hr->reduced_mat = BB_std_calloc_2d_double(hlpod_hr->reduced_mat, total_num_modes, total_num_modes);
    hlpod_hr->reduced_RH = BB_std_calloc_1d_double(hlpod_hr->reduced_RH, total_num_modes);

}

void HROM_ecm_memory_allocation_online(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        HLPOD_HR*       hlpod_hr)
{
    hlpod_hr->reduced_mat = BB_std_calloc_2d_double(hlpod_hr->reduced_mat, total_num_modes, total_num_modes);
    hlpod_hr->reduced_RH = BB_std_calloc_1d_double(hlpod_hr->reduced_RH, total_num_modes);
}


void HROM_ecm_get_selected_elems(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
        HLPOD_HR*       hlpod_hr,
		const char*		directory)
{
    int nl = fe->local_num_nodes;

    const int max_ITER = 2000;
    const double TOL = 1.0e-9;

    double* ans_vec;
    ans_vec = BB_std_calloc_1d_double(ans_vec, total_num_elem);

    double residual;

    int NNLS_row = total_num_snapshot * total_num_modes;

    monolis_optimize_nnls_R_with_sparse_solution(
        hlpod_hr->matrix, 
        hlpod_hr->RH, 
        ans_vec, NNLS_row, total_num_elem, max_ITER, TOL, &residual);

    printf("\n\nmax_iter = %d, tol = %lf, residuals = %lf\n\n", max_ITER, TOL, residual);

    int index = 0;

    for(int i = 0; i < total_num_elem; i++){
        if(ans_vec[i] != 0.0){
            index++;
        }
    }

    int total_num_selected_elems = index;

    printf("\n\n num_selected_elems = %d\n\n", index);

	hr_write_NNLS_residual(residual, 0, 0, directory);
	hr_write_NNLS_num_elems(total_num_selected_elems, 0, 0, directory);

    int* total_id_selected_elems;
    double* total_elem_weight;

    total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems);
    total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems);

    bool* bool_elem;  
    bool_elem = BB_std_calloc_1d_bool(bool_elem, total_num_selected_elems);

	hlpod_hr->D_bc_exists = BB_std_calloc_1d_bool(hlpod_hr->D_bc_exists, fe->total_num_nodes);

    index = 0;
    for(int i = 0; i < total_num_elem; i++){
        if(ans_vec[i] != 0.0){
            total_id_selected_elems[index] = i;
            total_elem_weight[index] = ans_vec[i];
            index++;
        }
    }

    for(int h=0; h<(total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
                    bool_elem[h] = true;
					hlpod_hr->D_bc_exists[index_j] = true;
                }
            }
        }
    }

	index = 0;
	for(int i = 0; i < fe->total_num_nodes; i++) {       //六面体一次要素は8
        if( hlpod_hr->D_bc_exists[i]){
			index++;
			printf("i = %d index = %d\n", i, index);
		}
    }

	//D_bcを付与する節点
	hlpod_hr->num_selected_nodes_D_bc = index;
	hlpod_hr->D_bc_node_id = BB_std_calloc_1d_int(hlpod_hr->D_bc_node_id, hlpod_hr->num_selected_nodes_D_bc);

	index = 0;
    for(int i = 0; i < fe->total_num_nodes; i++) {
        if( hlpod_hr->D_bc_exists[i]) {
			hlpod_hr->D_bc_node_id[index] = i;
			index++;
        }
    }

    index = 0;
    for(int h=0; h<(total_num_selected_elems); h++) {
        if(bool_elem[h]){
            index++;
        } 
    }

	printf("\n\n num_elem_D_bc = %d \n\n", index);

    //index = D_bcが付与された要素数
    hlpod_hr->num_selected_elems = total_num_selected_elems - index;
    hlpod_hr->num_selected_elems_D_bc = index;

    hlpod_hr->id_selected_elems = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems, hlpod_hr->num_selected_elems);
    hlpod_hr->id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems_D_bc, hlpod_hr->num_selected_elems_D_bc);
    
    hlpod_hr->elem_weight = BB_std_calloc_1d_double(hlpod_hr->elem_weight, hlpod_hr->num_selected_elems);
    hlpod_hr->elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_hr->elem_weight_D_bc, hlpod_hr->num_selected_elems_D_bc);

    int index1 = 0;
    int index2 = 0;
    for(int h=0; h<(total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

        if(bool_elem[h]) {
            hlpod_hr->id_selected_elems_D_bc[index1] = total_id_selected_elems[h];
            hlpod_hr->elem_weight_D_bc[index1] = total_elem_weight[h];
            index1++;
        }
        else{
            hlpod_hr->id_selected_elems[index2] = total_id_selected_elems[h];
            hlpod_hr->elem_weight[index2] = total_elem_weight[h];
            index2++;
        }
    }

    FILE* fp1;
    FILE* fp2;
    char fname1[BUFFER_SIZE];
    char fname2[BUFFER_SIZE];

    snprintf(fname1, BUFFER_SIZE,"DDECM/lb_selected_elem_D_bc.%d.txt", 0);
    snprintf(fname2, BUFFER_SIZE,"DDECM/lb_selected_elem.%d.txt", 0);

    fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
    fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

    fprintf(fp1, "%d\n", index1);
    fprintf(fp2, "%d\n", index2);

    int index_1 = 0;
    int index_2 = 0;

    for(int h=0; h<(total_num_selected_elems); h++) {
        if(bool_elem[h]) {
            fprintf(fp1, "%d %.30e\n", hlpod_hr->id_selected_elems_D_bc[index_1], hlpod_hr->elem_weight_D_bc[index_1]);
            index_1++;
        }
        else{
            fprintf(fp2, "%d %.30e\n", hlpod_hr->id_selected_elems[index_2], hlpod_hr->elem_weight[index_2]);

            index_2++;
        }
    }
    
    fclose(fp1);
    fclose(fp2);

    /*input
    hlpod_hr->matrix, 
    hlpod_hr->RH,
    NNLS conditions
    */

    /*output
    hlpod_hr->num_selected_elem
    hlpod_hr->id_selected_elems
    hlpod_hr->elem_weight
    */

}


void HROM_ecm_calc_solution(
	BBFE_DATA* 		fe,
	HLPOD_MAT*      hlpod_mat,
	double*         HR_T,
	BBFE_BC*     	bc,
    int 			num_base)
{
	int nl = fe->total_num_nodes;
	int k = num_base;

	for(int j = 0; j < nl; j++){
		HR_T[j] = 0.0;
	}

	for(int i = 0; i < k; i++){
		for(int j = 0; j < nl; j++){
			HR_T[j] += hlpod_mat->pod_modes[j][i] * hlpod_mat->mode_coef[i];
		}
	}
}


void HROM_ecm_monolis_set_matrix(
	MONOLIS*     	monolis,
	HLPOD_HR*      hlpod_hr,
    const int 		num_basis)
{
	const int k = num_basis;
	int* index;
	int* item;
	int* connectivity;

	index = BB_std_calloc_1d_int(index, 1);
	item = BB_std_calloc_1d_int(item, 1);
	
	index[0] = 0;
	item[0] = 1;

	monolis_get_nonzero_pattern_by_nodal_graph_R(
		monolis,
		1,					//nnode:節点数
		k,					//ndof:節点あたりの自由度
		index,				//節点グラフを定義するindex配列
		item);				//設定グラフを定義するitem配列

	connectivity = BB_std_calloc_1d_int(connectivity, 1);
	connectivity[0] = 0;

	monolis_add_matrix_to_sparse_matrix_R(
		monolis,					//MONOLIS* mat,
  		1,							//int      n_base,
  		connectivity,				//int*     connectivity,
		hlpod_hr->reduced_mat);		//double** val

	BB_std_free_1d_int(index, 1);
	BB_std_free_1d_int(item, 1);
	BB_std_free_1d_int(connectivity, 1);
}

void HROM_ecm_to_monollis_rhs(
	MONOLIS*		monolis,
    HLPOD_HR*       hlpod_rh,
	const int 		k)
{
	for(int i = 0; i < k; i++){
		monolis->mat.R.B[i] = hlpod_rh->reduced_RH[i];
	}
}

void HROM_ecm_set_selected_elems(
		BBFE_DATA*     	fe,
        HLPOD_HR*       hlpod_hr,
        const int		total_num_nodes,
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

    for(int h=0; h<(hlpod_hr->num_selected_elems); h++) {
        int e = hlpod_hr->id_selected_elems[h];
        for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
            
            ECM_elem[index] = myrank + 1;
            ECM_elem_weight[index] += hlpod_hr->elem_weight[h];
        }
    }

	for(int h=0; h<(hlpod_hr->num_selected_elems_D_bc); h++) {
        int e = hlpod_hr->id_selected_elems_D_bc[h];
        for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
            ECM_elem[index] = myrank + 1;
            ECM_elem_weight[index] += hlpod_hr->elem_weight_D_bc[h];
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

	fclose(fp);

}



void HROM_ecm_read_selected_elems(
    HLPOD_HR*     hlpod_hr,
	const int num_subdomains,
	const char* directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;


	FILE* fp1;
	FILE* fp2;
	FILE* fp3;
	FILE* fp4;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	char fname3[BUFFER_SIZE];
	char fname4[BUFFER_SIZE];

	int Index1 = 0;
	int Index2 = 0;
	int tmp;
	double val;
	int index1 = 0;
	int index2 = 0;
	int num_selected_elems;
	int num_selected_elems_D_bc;

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", m);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", m);

		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));
		Index1 += num_selected_elems;
		Index2 += num_selected_elems_D_bc;

		fclose(fp1);
		fclose(fp2);
	}

    hlpod_hr->num_selected_elems = Index1;
    hlpod_hr->num_selected_elems_D_bc  = Index2;

    hlpod_hr->elem_weight = BB_std_calloc_1d_double(hlpod_hr->elem_weight, Index1);
    hlpod_hr->id_selected_elems = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems, Index1);
    hlpod_hr->elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_hr->elem_weight_D_bc, Index2);
    hlpod_hr->id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems_D_bc, Index2);

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", m);

		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
        for(int i = 0; i < num_selected_elems; i++){
            fscanf(fp1, "%d %lf", &(hlpod_hr->id_selected_elems[index1]), &(hlpod_hr->elem_weight[index1]));
            printf("%d %lf\n", hlpod_hr->id_selected_elems[index1], hlpod_hr->elem_weight[index1]);
            index1++;
        }

		fclose(fp1);
	}

    for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", m);

		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(num_selected_elems_D_bc));

        for(int i = 0; i < num_selected_elems_D_bc; i++){
            fscanf(fp2, "%d %lf", &(hlpod_hr->id_selected_elems_D_bc[index2]), &(hlpod_hr->elem_weight_D_bc[index2]));
            
            printf("%d %lf\n", hlpod_hr->id_selected_elems_D_bc[index2], hlpod_hr->elem_weight_D_bc[index2]);
            index2++;
        }

		fclose(fp2);
	}

}
