//DDHROMに関して、オーバーラップ要素を含んで計算する方式
//内部要素の総和が分割前の要素数の総和になることを利用
//load balancing, 2階層目のbscr形式に対応

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>


#include "DDHR_para_lb.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";

static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";

//内部要素とオーバーラップ要素の出力
void HROM_ddecm_get_selected_elems_int_ovl(
	HLPOD_DDHR*     hlpod_ddhr,
	const char*     directory)
{
	double t = monolis_get_time_global_sync();

	const int myrank = monolis_mpi_get_global_my_rank();
	int num_subdomains;

	FILE* fp;
	FILE* fp1;
	FILE* fp2;
	char fname[BUFFER_SIZE];
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	char id[BUFFER_SIZE];

	int val;
	int ndof;
	int*    ovl_selected_elems;
	int*    ovl_selected_elems_D_bc;
	double* ovl_selected_elems_weight;
	double* ovl_selected_elems_weight_D_bc;

	int num_selected_elems = 0;
	int num_selected_elems_D_bc = 0;
	int meta_n_neib;

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.n_internal.%d", myrank);
	fp = BBFE_sys_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(ndof));
	fscanf(fp, "%d", &(num_subdomains));
	fclose(fp);

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = BBFE_sys_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	for (int n = 0; n < num_subdomains; n++) {
		num_selected_elems = 0;
		num_selected_elems_D_bc = 0;

		printf("num_subdomains = %d, n = %d \n\n", num_subdomains, n);

		snprintf(fname, BUFFER_SIZE, "parted.1/%s.recv.%d", INPUT_FILENAME_NODE, subdomain_id[n]);
		fp = BBFE_sys_read_fopen(fp, fname, directory);
		fscanf(fp, "%d %d", &(meta_n_neib), &(ndof));
		int* meta_list_neib;
		meta_list_neib = BB_std_calloc_1d_int(meta_list_neib, meta_n_neib);
		for (int i = 0; i < meta_n_neib; i++) {
			fscanf(fp, "%d", &(meta_list_neib[i]));
		}
		fclose(fp);

		/*自領域*/
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[n]);
		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[n]);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(val));
		num_selected_elems_D_bc += val;
		fclose(fp2);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", meta_list_neib[m]);
			fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

			fscanf(fp1, "%d", &(val));
			num_selected_elems += val;
			fclose(fp1);
		}

		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", meta_list_neib[m]);
			fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

			fscanf(fp2, "%d", &(val));
			num_selected_elems_D_bc += val;
			fclose(fp2);
		}

		ovl_selected_elems = BB_std_calloc_1d_int(ovl_selected_elems, num_selected_elems);
		ovl_selected_elems_weight = BB_std_calloc_1d_double(ovl_selected_elems_weight, num_selected_elems);
		ovl_selected_elems_D_bc = BB_std_calloc_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
		ovl_selected_elems_weight_D_bc = BB_std_calloc_1d_double(ovl_selected_elems_weight_D_bc, num_selected_elems_D_bc);

		printf("num_subdomains = %d, n = %d \n\n", num_subdomains, n);

		int index = 0;

		/*自領域*/
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[n]);
		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);
		fscanf(fp1, "%d", &(val));
		for (int j = 0; j < val; j++) {
			fscanf(fp1, "%d %lf", &(ovl_selected_elems[j + index]), &(ovl_selected_elems_weight[j + index]));
		}
		index += val;

		fclose(fp1);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", meta_list_neib[m]);
			fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);
			fscanf(fp1, "%d", &(val));
			for (int j = 0; j < val; j++) {
				fscanf(fp1, "%d %lf", &(ovl_selected_elems[j + index]), &(ovl_selected_elems_weight[j + index]));
			}
			index += val;
			fclose(fp1);
		}

		index = 0;
		/*自領域*/
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[n]);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);
		fscanf(fp2, "%d", &(val));
		for (int j = 0; j < val; j++) {
			fscanf(fp2, "%d %lf", &(ovl_selected_elems_D_bc[j + index]), &(ovl_selected_elems_weight_D_bc[j + index]));
		}
		index += val;
		fclose(fp2);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", meta_list_neib[m]);
			fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);
			fscanf(fp2, "%d", &(val));
			for (int j = 0; j < val; j++) {
				fscanf(fp2, "%d %lf", &(ovl_selected_elems_D_bc[j + index]), &(ovl_selected_elems_weight_D_bc[j + index]));
			}
			index += val;
			fclose(fp2);
		}

		bool* bool_ovl_selected_elems;
		bool* bool_ovl_selected_elems_D_bc;

		bool_ovl_selected_elems = BB_std_calloc_1d_bool(bool_ovl_selected_elems, num_selected_elems);
		bool_ovl_selected_elems_D_bc = BB_std_calloc_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);

		int* ovl_elem_local_id;
		int* ovl_elem_local_id_D_bc;

		ovl_elem_local_id = BB_std_calloc_1d_int(ovl_elem_local_id, num_selected_elems);
		ovl_elem_local_id_D_bc = BB_std_calloc_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);

		int* ovl_elem_global_id;
		int total_num_elems;
		int tmp;

		//読み込む対象の要素のidを読み込み
		snprintf(fname, BUFFER_SIZE, "parted.1/%s.%d", INPUT_FILENAME_ELEM_ID, subdomain_id[n]);
		fp = BBFE_sys_read_fopen(fp, fname, directory);
		fscanf(fp, "%s", id);
		fscanf(fp, "%d %d", &(total_num_elems), &(tmp));
		ovl_elem_global_id = BB_std_calloc_1d_int(ovl_elem_global_id, total_num_elems);
		for (int i = 0; i < total_num_elems; i++) {
			fscanf(fp, "%d", &(ovl_elem_global_id[i]));
		}
		fclose(fp);

		printf("num_subdomains = %d, n = %d \n\n", num_subdomains, n);

		int index1 = 0;
		int index2 = 0;

		//global idのセット
		for (int i = 0; i < num_selected_elems; i++) {
			for (int j = 0; j < total_num_elems; j++) {
				if (ovl_selected_elems[i] == ovl_elem_global_id[j]) {
					bool_ovl_selected_elems[i] = true;
					ovl_elem_local_id[index1] = j;
					index1++;
				}
			}
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			for (int j = 0; j < total_num_elems; j++) {
				if (ovl_selected_elems_D_bc[i] == ovl_elem_global_id[j]) {
					bool_ovl_selected_elems_D_bc[i] = true;
					ovl_elem_local_id_D_bc[index2] = j;
					index2++;
				}
			}
		}

		snprintf(fname1, BUFFER_SIZE, "DDECM/selected_elem_overlap.%d.txt", subdomain_id[n]);
		fp1 = BBFE_sys_write_fopen(fp1, fname1, directory);

		index1 = 0;
		index2 = 0;

		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {
				index1++;
			}
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				index2++;
			}
		}

		fprintf(fp1, "%d\n", index1 + index2);
		index1 = 0;
		index2 = 0;
		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {
				fprintf(fp1, "%d %.15g\n", ovl_elem_local_id[index1], ovl_selected_elems_weight[i]);
				index1++;
			}
		}
		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				fprintf(fp1, "%d %.15g\n", ovl_elem_local_id_D_bc[index2], ovl_selected_elems_weight_D_bc[i]);
				index2++;
			}
		}

        fclose(fp1);

		//内部要素の選定
		int local_dof;
		int n_internal;
		int** conn;

		snprintf(fname, BUFFER_SIZE, "parted.1/%s.%d", INPUT_FILENAME_ELEM, subdomain_id[n]);
		fp = BBFE_sys_read_fopen(fp, fname, directory);

		fscanf(fp, "%d %d", &(total_num_elems), &(local_dof));
		conn = BB_std_calloc_2d_int(conn, total_num_elems, local_dof);
		for (int i = 0; i < total_num_elems; i++) {
			for (int j = 0; j < local_dof; j++) {
				fscanf(fp, "%d", &(conn[i][j]));
			}
		}
		fclose(fp);

		snprintf(fname, BUFFER_SIZE, "parted.1/node.dat.n_internal.%d", subdomain_id[n]);
		fp = BBFE_sys_read_fopen(fp, fname, directory);
		fscanf(fp, "%s %d", id, &(tmp));
		fscanf(fp, "%d", &(n_internal));
		fclose(fp);

		/*節点ベースの出力*/
		index1 = 0;
		index2 = 0;

		const int nl = 8; //六面体一次要素限定 今後引数にする

		int num_selected_nodes = 0;
		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {	
				int index = ovl_elem_local_id[index1];

				for(int i=0; i<nl; i++) {       //六面体一次要素は8
					if (conn[index][i] < n_internal ) {
						num_selected_nodes++;
					}
				}
				index1++;
			}
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				int index = ovl_elem_local_id_D_bc[index2];

				for(int i=0; i<nl; i++) {       //六面体一次要素は8
					if (conn[index][i] < n_internal ) {
						num_selected_nodes++;
					}
				}
				index2++;
			}
		}

		snprintf(fname1, BUFFER_SIZE, "DDECM/num_selected_node.%d.txt", subdomain_id[n]);
		fp1 = BBFE_sys_write_fopen(fp1, fname1, directory);
		fprintf(fp1, "%d\n", num_selected_nodes);
		fclose(fp1);
		/**************/

		index1 = 0;
		index2 = 0;

		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {
				int index = ovl_elem_local_id[index1];

				for (int j = 0; j < local_dof; j++) {

					if (conn[index][j] > n_internal) {

						bool_ovl_selected_elems[i] = false;
					}
				}
				index1++;
			}
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				int index = ovl_elem_local_id_D_bc[index2];

				for (int j = 0; j < local_dof; j++) {

					if (conn[index][j] > n_internal) {
						bool_ovl_selected_elems_D_bc[i] = false;
					}
				}
				index2++;
			}
		}

		snprintf(fname1, BUFFER_SIZE, "DDECM/selected_elem_internal.%d.txt", subdomain_id[n]);

		fp1 = BBFE_sys_write_fopen(fp1, fname1, directory);

		index1 = 0;
		index2 = 0;
		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {
				index1++;
			}
		}
		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				index2++;
			}
		}
		fprintf(fp1, "%d\n", index1 + index2);

		index1 = 0;
		index2 = 0;
		for (int i = 0; i < num_selected_elems; i++) {
			if (bool_ovl_selected_elems[i]) {
				fprintf(fp1, "%d %.15g\n", ovl_elem_local_id[index1], ovl_selected_elems_weight[i]);
				index1++;
			}
		}
		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			if (bool_ovl_selected_elems_D_bc[i]) {
				fprintf(fp1, "%d %.15g\n", ovl_elem_local_id_D_bc[index2], ovl_selected_elems_weight_D_bc[i]);
				index2++;
			}
		}

        fclose(fp1);

		BB_std_free_1d_int(ovl_elem_global_id, total_num_elems);

		BB_std_free_2d_int(conn, total_num_elems, local_dof);

		BB_std_free_1d_int(ovl_selected_elems, num_selected_elems);
		BB_std_free_1d_double(ovl_selected_elems_weight, num_selected_elems);
		BB_std_free_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
		BB_std_free_1d_double(ovl_selected_elems_weight_D_bc, num_selected_elems_D_bc);

		BB_std_free_1d_bool(bool_ovl_selected_elems, num_selected_elems);
		BB_std_free_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);
		BB_std_free_1d_int(ovl_elem_local_id, num_selected_elems);
		BB_std_free_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);

		BB_std_free_1d_int(meta_list_neib, meta_n_neib);
	}
    double t_tmp = monolis_get_time_global_sync();
}

void HROM_ddecm_read_selected_elems_para(
	const int num_subdomains,
	const char* directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	double t = monolis_get_time_global_sync();

	FILE* fp1;
	FILE* fp2;
	FILE* fp3;
	FILE* fp4;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	char fname3[BUFFER_SIZE];
	char fname4[BUFFER_SIZE];

	snprintf(fname3, BUFFER_SIZE, "DDECM/selected_elem_D_bc.%d.txt", monolis_mpi_get_global_my_rank());
	snprintf(fname4, BUFFER_SIZE, "DDECM/selected_elem.%d.txt", monolis_mpi_get_global_my_rank());

	fp3 = ROM_BB_write_fopen(fp3, fname3, directory);
	fp4 = ROM_BB_write_fopen(fp4, fname4, directory);

	int Index1 = 0;
	int Index2 = 0;
	int tmp;
	double val;
	int index1 = 0;
	int index2 = 0;
	int num_selected_elems;
	int num_selected_elems_D_bc;

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));
		Index1 += num_selected_elems;
		Index2 += num_selected_elems_D_bc;

		fclose(fp1);
		fclose(fp2);
	}

	fprintf(fp3, "%d\n", Index1);
	fprintf(fp4, "%d\n", Index2);

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));

		for (int i = 0; i < num_selected_elems; i++) {
			fscanf(fp1, "%d %lf", &(tmp), &(val));
			fprintf(fp3, "%d %.30e\n", tmp, val);
			index1++;
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			fscanf(fp2, "%d %lf", &(tmp), &(val));
			fprintf(fp4, "%d %.30e\n", tmp, val);
			index2++;
		}

		fclose(fp1);
		fclose(fp2);
	}

	fclose(fp3);
	fclose(fp4);

	t = monolis_get_time_global_sync();
}


double ddhr_calc_tol(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
    HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int 		num_subdomains)
{
	int NNLS_row = total_num_snapshot*hlpod_vals->n_neib_vec;	//2は残差ベクトル＋右辺ベクトルを採用しているため
    double norm = 0.0;
    double* RH;

	int index_NNLS1 = 0;
	int index_NNLS2 = 0;

   double t = monolis_get_time_global_sync();

	for (int m = 0; m < num_subdomains; m++) {
		int NNLS_row = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot;
		RH = BB_std_calloc_1d_double(RH, NNLS_row);
		index_NNLS2 = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot;

		for (int p = 0; p < total_num_snapshot; p++) {

			for (int j = hlpod_ddhr->num_internal_modes_1stdd_sum[m] + hlpod_vals->n_neib_vec * p; j < hlpod_ddhr->num_internal_modes_1stdd_sum[m + 1] + hlpod_vals->n_neib_vec * p; j++) {
				RH[index_NNLS1] = hlpod_ddhr->RH[j][m];
				index_NNLS1++;
			}

			int iS = hlpod_meta->index[m];
			int iE = hlpod_meta->index[m + 1];
			for (int n = iS; n < iE; n++) {
				for (int l = 0; l < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; l++) {
					if (hlpod_meta->my_global_id[hlpod_meta->item[n]] == hlpod_meta->global_id[l]) {
						int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[l] + hlpod_vals->n_neib_vec * p;
						int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[l + 1] + hlpod_vals->n_neib_vec * p;

						for (int j = IS; j < IE; j++) {
							RH[index_NNLS1] = hlpod_ddhr->RH[j][m];

							index_NNLS1++;
						}

					}
				}
			}

		}
        
        index_NNLS1 = 0;
		index_NNLS2 = 0;

		for(int j = 0; j < NNLS_row; j++){
			norm += RH[j]*RH[j];
		}

        BB_std_free_1d_double(RH, NNLS_row);
    }

    monolis_allreduce_R(
        1,
        &norm,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    return norm;
    
}


void HROM_ddecm_write_selected_elems_svd(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*      hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

    const int comm = monolis_mpi_get_self_comm();
    int scalapack_comm;
    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;

	int index_1 = 0;
	int index_2 = 0;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	int nl = fe->local_num_nodes;
    double t1 = monolis_get_time_global_sync();

	const int max_ITER = 1000;
	const double TOL = 1.0e-10;

	double residual;

	double* ans_vec;
	double** matrix;
	double* RH;
	bool** bool_elem;
	int* total_id_selected_elems;
	double* total_elem_weight;
	int* total_num_selected_elems;

	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);
	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
	hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);

    //double global_norm = ddhr_calc_tol(monolis_com,
    //    hlpod_vals, hlpod_ddhr,	hlpod_mat, hlpod_meta, total_num_elem, total_num_snapshot, num_subdomains);

	int Index1 = 0;
	int Index2 = 0;
	int index_NNLS1 = 0;
	int index_NNLS2 = 0;

    for (int m = 0; m < num_subdomains; m++) {
		//int NNLS_row = 2 * hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot + 1; //2は残差ベクトル＋右辺ベクトルを採用しているため
		int NNLS_row = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot + 1; //2は残差ベクトル＋右辺ベクトルを採用しているため

		ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		RH = BB_std_calloc_1d_double(RH, NNLS_row);

        for(int j = 0; j < NNLS_row; j++){
            for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
                matrix[j][i] = hlpod_ddhr->matrix[j][i][m];
            }
            RH[j] = hlpod_ddhr->RH[j][m];
        }

        double local_norm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
			local_norm += RH[j]*RH[j];
		}

        //double input_TOL = TOL * sqrt(global_norm) / (num_subdomains  * sqrt(local_norm));

        double** S = BB_std_calloc_2d_double(S, NNLS_row, hlpod_ddhr->num_elems[m]);
        double* V = BB_std_calloc_1d_double(V, hlpod_ddhr->num_elems[m]);
        double** D = BB_std_calloc_2d_double(D, hlpod_ddhr->num_elems[m], hlpod_ddhr->num_elems[m]);

        double t1 = monolis_get_time();
        monolis_scalapack_gesvd_R(
                NNLS_row,
                hlpod_ddhr->num_elems[m], 
                matrix,
                S, 
                V, 
                D, 
                comm,
                scalapack_comm);
        double t2 = monolis_get_time();

        int k = ROM_BB_estimate_num_pod_modes(
            V,
            hlpod_ddhr->num_elems[m],
            10e-10);

        k = 100;

        double** S_k = BB_std_calloc_2d_double(S_k, NNLS_row, k);

        for(int i = 0; i < NNLS_row; i++){
            for(int j = 0; j < k; j++){
                S_k[i][j] = S[i][j];
            }
        }

        double** G_k = BB_std_calloc_2d_double(G_k, k, hlpod_ddhr->num_elems[m]);
        double* b_k = BB_std_calloc_1d_double(b_k, k);

        ROM_BB_transposemat_mat(
            S_k,
            matrix,
            G_k,
            NNLS_row,
            k,
            hlpod_ddhr->num_elems[m]);
        
        ROM_BB_transposemat_vec(
            S_k,
            RH,
            b_k,
            NNLS_row,
            k);

            residual = 0.0;

        monolis_optimize_nnls_R_with_sparse_solution(
            G_k,
            b_k,
            ans_vec, k, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);

		printf("\n\nmax_iter = %d, tol = %lf, residuals = %lf\n\n", max_ITER, TOL, residual);

		int index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				index++;
			}
		}

		total_num_selected_elems[m] = index;


		hr_write_NNLS_residual(residual, myrank, m, directory);
		hr_write_NNLS_num_elems(total_num_selected_elems[m], myrank, m, directory);

		total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[m]);

		index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][m];
				total_elem_weight[index] = ans_vec[i];
				index++;
			}
		}

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			for (int i = 0; i < nl; i++) {       //六面体一次要素は8
				for (int j = 0; j < nl; j++) {
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

                    for(int k = 0; k < dof; k++) {
					    if (bc->D_bc_exists[index_j*dof + k]) {
						    bool_elem[h][m] = true;
    						hlpod_ddhr->D_bc_exists[index_j][m] = true;
	    				}
                    }
				}
			}
		}

		printf("\n\n num_elem = %d \n\n", index);

		index = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				index++;
			}
		}

		//printf("\n\n num_elem_D_bc = %d \n\n", index);

		hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems[m] - index;
		hlpod_ddhr->num_selected_elems_D_bc[m] = index;
		int index1 = 0;
		int index2 = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			if (bool_elem[h][m]) {
				hlpod_ddhr->id_selected_elems_D_bc[index1][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight_D_bc[index1][m] = total_elem_weight[h];

				index1++;
			}
			else {
				hlpod_ddhr->id_selected_elems[index2][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight[index2][m] = total_elem_weight[h];

				index2++;
			}
		}

		Index1 += index1;
		Index2 += index2;

		BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		BB_std_free_1d_double(RH, NNLS_row);

		BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[m]);

		double t = monolis_get_time_global_sync();

		FILE* fp1;
		FILE* fp2;
		char fname1[BUFFER_SIZE];
		char fname2[BUFFER_SIZE];

		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

		fprintf(fp1, "%d\n", index1);
		fprintf(fp2, "%d\n", index2);

		index_1 = 0;
		index_2 = 0;

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				fprintf(fp1, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems_D_bc[index_1][m]], hlpod_ddhr->elem_weight_D_bc[index_1][m]);
				index_1++;
			}
			else {
				fprintf(fp2, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems[index_2][m]], hlpod_ddhr->elem_weight[index_2][m]);

				index_2++;
			}
		}

		fclose(fp1);
		fclose(fp2);

	}

	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);

	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	BB_std_free_3d_double(hlpod_ddhr->matrix, total_num_snapshot * hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
	BB_std_free_2d_double(hlpod_ddhr->RH, total_num_snapshot * hlpod_vals->n_neib_vec, num_subdomains);

    monolis_scalapack_comm_finalize(scalapack_comm);

	/*input
	hlpod_ddhr->matrix,
	hlpod_ddhr->RH,
	NNLS conditions
	*/

	/*output
	hlpod_ddhr->num_selected_elem
	hlpod_ddhr->id_selected_elems
	hlpod_ddhr->elem_weight
	*/

	double t = monolis_get_time_global_sync();
}


void HROM_ddecm_write_selected_elems_para(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;

	int index_1 = 0;
	int index_2 = 0;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	int nl = fe->local_num_nodes;

	printf("\n\nmyrank = %d, num_subdomains = %d\n\n", myrank, num_subdomains);
	printf("\n\nnum_elems1 = %d\n\n", hlpod_ddhr->num_elems[0]);
    double t1 = monolis_get_time_global_sync();

	const int max_ITER = 10000;
	const double TOL = 1.0e-14;

	double residual;

	double* ans_vec;
	double** matrix;
	double* RH;
	bool** bool_elem;
	int* total_id_selected_elems;
	double* total_elem_weight;
	int* total_num_selected_elems;

	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);
	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
	hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);
  //  double global_norm = ddhr_calc_tol(monolis_com,
   //     hlpod_vals, hlpod_ddhr,	hlpod_mat, hlpod_meta, total_num_elem, total_num_snapshot, num_subdomains);

	int Index1 = 0;
	int Index2 = 0;
	int index_NNLS1 = 0;
	int index_NNLS2 = 0;


    for (int m = 0; m < num_subdomains; m++) {
		int NNLS_row = 2 * hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot + 1; //2は残差ベクトル＋右辺ベクトルを採用しているため

		ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		RH = BB_std_calloc_1d_double(RH, NNLS_row);

        for(int j = 0; j < NNLS_row; j++){
            for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
                matrix[j][i] = hlpod_ddhr->matrix[j][i][m];
            }
            RH[j] = hlpod_ddhr->RH[j][m];
        }

        double local_norm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
			local_norm += RH[j]*RH[j];
		}

//        double input_TOL = TOL * sqrt(global_norm) / (num_subdomains  * sqrt(local_norm));

		index_NNLS1 = 0;
		index_NNLS2 = 0;

		residual = 0.0;

		monolis_optimize_nnls_R_with_sparse_solution(
			matrix,
			RH,
			ans_vec, NNLS_row, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);

		printf("\n\nmax_iter = %d, tol = %lf, residuals = %lf\n\n", max_ITER, TOL, residual);

		int index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				index++;
			}
		}

		total_num_selected_elems[m] = index;

		printf("\n\nnum_selected_elems = %d\n\n", index);

		hr_write_NNLS_residual(residual, myrank, m, directory);
		hr_write_NNLS_num_elems(total_num_selected_elems[m], myrank, m, directory);

		total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[m]);

		index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][m];
				total_elem_weight[index] = ans_vec[i];
				index++;
			}
		}

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			for (int i = 0; i < nl; i++) {       //六面体一次要素は8
				for (int j = 0; j < nl; j++) {
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

                    for(int k = 0; k < dof; k++) {
					    if (bc->D_bc_exists[index_j*dof + k]) {
						    bool_elem[h][m] = true;
    						hlpod_ddhr->D_bc_exists[index_j][m] = true;
	    				}
                    }
				}
			}
		}

		printf("\n\n test_num_elem_D_bc = %d \n\n", index);

		index = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				index++;
			}
		}

		//printf("\n\n num_elem_D_bc = %d \n\n", index);

		//index = D_bcが付与された要素数
		hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems[m] - index;
		hlpod_ddhr->num_selected_elems_D_bc[m] = index;

		int index1 = 0;
		int index2 = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			if (bool_elem[h][m]) {
				hlpod_ddhr->id_selected_elems_D_bc[index1][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight_D_bc[index1][m] = total_elem_weight[h];

				index1++;
			}
			else {
				hlpod_ddhr->id_selected_elems[index2][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight[index2][m] = total_elem_weight[h];

				index2++;
			}
		}

		Index1 += index1;
		Index2 += index2;

		BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		BB_std_free_1d_double(RH, NNLS_row);

		BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[m]);

		double t = monolis_get_time_global_sync();

		FILE* fp1;
		FILE* fp2;
		char fname1[BUFFER_SIZE];
		char fname2[BUFFER_SIZE];

		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

		fprintf(fp1, "%d\n", index1);
		fprintf(fp2, "%d\n", index2);

		index_1 = 0;
		index_2 = 0;

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				fprintf(fp1, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems_D_bc[index_1][m]], hlpod_ddhr->elem_weight_D_bc[index_1][m]);
				index_1++;
			}
			else {
				fprintf(fp2, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems[index_2][m]], hlpod_ddhr->elem_weight[index_2][m]);

				index_2++;
			}
		}

		fclose(fp1);
		fclose(fp2);

	}

	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);

	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	BB_std_free_3d_double(hlpod_ddhr->matrix, total_num_snapshot * hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
	BB_std_free_2d_double(hlpod_ddhr->RH, total_num_snapshot * hlpod_vals->n_neib_vec, num_subdomains);

	/*input
	hlpod_ddhr->matrix,
	hlpod_ddhr->RH,
	NNLS conditions
	*/

	/*output
	hlpod_ddhr->num_selected_elem
	hlpod_ddhr->id_selected_elems
	hlpod_ddhr->elem_weight
	*/

	double t = monolis_get_time_global_sync();
}


void HROM_ddecm_write_selected_elems_para_arbit_subd(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;

	int index_1 = 0;
	int index_2 = 0;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	int nl = fe->local_num_nodes;
    double t1 = monolis_get_time_global_sync();

	const int max_ITER = 10000;
	const double TOL = 1.0e-10;

	double residual;

	double* ans_vec;
	double** matrix;
	double* RH;
	bool** bool_elem;
	int* total_id_selected_elems;
	double* total_elem_weight;
	int* total_num_selected_elems;

	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);
	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
	hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);

	int Index1 = 0;
	int Index2 = 0;

	int index_NNLS1 = 0;
	int index_NNLS2 = 0;

    // global_norm = ddhr_calc_tol(monolis_com,
    //    hlpod_vals, hlpod_ddhr,	hlpod_mat, hlpod_meta, total_num_elem, total_num_snapshot, num_subdomains);

	for (int m = 0; m < num_subdomains; m++) {
		int NNLS_row = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot; //2は残差ベクトル＋右辺ベクトルを採用しているため
        printf("NNLS_row = %d num_elems = %d\n", NNLS_row, hlpod_ddhr->num_elems[m]);
		ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		RH = BB_std_calloc_1d_double(RH, NNLS_row);

		index_NNLS2 = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot;

		for (int p = 0; p < total_num_snapshot; p++) {

			for (int j = hlpod_ddhr->num_internal_modes_1stdd_sum[m] + hlpod_vals->n_neib_vec * p; j < hlpod_ddhr->num_internal_modes_1stdd_sum[m + 1] + hlpod_vals->n_neib_vec * p; j++) {
				for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
					matrix[index_NNLS1][i] = hlpod_ddhr->matrix[j][i][m];
				}
				RH[index_NNLS1] = hlpod_ddhr->RH[j][m];
				index_NNLS1++;
			}

			int iS = hlpod_meta->index[m];
			int iE = hlpod_meta->index[m + 1];
			for (int n = iS; n < iE; n++) {
				for (int l = 0; l < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; l++) {
					if (hlpod_meta->my_global_id[hlpod_meta->item[n]] == hlpod_meta->global_id[l]) {

						int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[l] + hlpod_vals->n_neib_vec * p;
						int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[l + 1] + hlpod_vals->n_neib_vec * p;

						for (int j = IS; j < IE; j++) {
							for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
								matrix[index_NNLS1][i] = hlpod_ddhr->matrix[j][i][m];
							}
							RH[index_NNLS1] = hlpod_ddhr->RH[j][m];

							index_NNLS1++;
						}

					}
				}
			}
		}

        double local_norm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
			local_norm += RH[j];
		}
        double local_Frovnorm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
            for(int e = 0; e < hlpod_ddhr->num_elems[m]; e++){
			    local_Frovnorm += matrix[j][e];
            }
		}

        //double scale = max(1.0, maxval(abs(local_norm)), maxval(abs(local_Frovnorm)));
        printf("local norm = %e, local_Fnorm = %e, scale = %e ", local_norm, local_Frovnorm);

        //local_norm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
			RH[j] /= local_Frovnorm;
		}
        //local_Frovnorm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
            for(int e = 0; e < hlpod_ddhr->num_elems[m]; e++){
			    matrix[j][e] /= local_Frovnorm;
            }
		}

		index_NNLS1 = 0;
		index_NNLS2 = 0;
		residual = 0.0;

		monolis_optimize_nnls_R_with_sparse_solution(
			matrix,
			RH,
			ans_vec, NNLS_row, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);

		printf("\n\nmax_iter = %d, tol = %e, residuals = %e\n\n", max_ITER, TOL, residual);

		int index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				index++;
			}
		}

		total_num_selected_elems[m] = index;

		printf("\n\nnum_selected_elems = %d\n\n", index);

		hr_write_NNLS_residual(residual, myrank, m, directory);
		hr_write_NNLS_num_elems(total_num_selected_elems[m], myrank, m, directory);

		total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[m]);

		index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][m];
				total_elem_weight[index] = ans_vec[i];
				index++;
			}
		}

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			for (int i = 0; i < nl; i++) {       //六面体一次要素は8
				for (int j = 0; j < nl; j++) {
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

                    for(int k = 0; k < dof; k++) {
					    if (bc->D_bc_exists[index_j*dof + k]) {
						    bool_elem[h][m] = true;
    						hlpod_ddhr->D_bc_exists[index_j][m] = true;
	    				}
                    }
				}
			}
		}

		index = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				index++;
			}
		}

		//printf("\n\n num_elem_D_bc = %d \n\n", index);

		//index = D_bcが付与された要素数
		hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems[m] - index;
		hlpod_ddhr->num_selected_elems_D_bc[m] = index;

		int index1 = 0;
		int index2 = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			if (bool_elem[h][m]) {
				hlpod_ddhr->id_selected_elems_D_bc[index1][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight_D_bc[index1][m] = total_elem_weight[h];

				index1++;
			}
			else {
				hlpod_ddhr->id_selected_elems[index2][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight[index2][m] = total_elem_weight[h];

				index2++;
			}
		}

		Index1 += index1;
		Index2 += index2;

		BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		BB_std_free_1d_double(RH, NNLS_row);

		BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[m]);

		double t = monolis_get_time_global_sync();

		FILE* fp1;
		FILE* fp2;
		char fname1[BUFFER_SIZE];
		char fname2[BUFFER_SIZE];

		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

		fprintf(fp1, "%d\n", index1);
		fprintf(fp2, "%d\n", index2);

		index_1 = 0;
		index_2 = 0;

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				fprintf(fp1, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems_D_bc[index_1][m]], hlpod_ddhr->elem_weight_D_bc[index_1][m]);
				index_1++;
			}
			else {
				fprintf(fp2, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems[index_2][m]], hlpod_ddhr->elem_weight[index_2][m]);

				index_2++;
			}
		}

		fclose(fp1);
		fclose(fp2);

	}

	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);


	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	BB_std_free_3d_double(hlpod_ddhr->matrix, total_num_snapshot * hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
	BB_std_free_2d_double(hlpod_ddhr->RH, total_num_snapshot * hlpod_vals->n_neib_vec, num_subdomains);

	/*input
	hlpod_ddhr->matrix,
	hlpod_ddhr->RH,
	NNLS conditions
	*/

	/*output
	hlpod_ddhr->num_selected_elem
	hlpod_ddhr->id_selected_elems
	hlpod_ddhr->elem_weight
	*/

	double t = monolis_get_time_global_sync();
}


void HROM_ecm_get_meta_neib(
	MONOLIS_COM*  	monolis_com,
	HLPOD_META*		hlpod_meta,
	const char*     directory)
{
	/*ファイル読み込み関連*/
	int num_metagraph_nodes;
	int tmp;

	char filename[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	FILE* fp;

	int* n_internal;

	snprintf(filename, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.recv.%d", monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, filename, directory);

	fscanf(fp, "%d %d", &(hlpod_meta->num_neib) , &(tmp));
	hlpod_meta->neib_id = BB_std_calloc_1d_int(hlpod_meta->neib_id, hlpod_meta->num_neib);

	for(int i = 0; i < hlpod_meta->num_neib; i++) {
		fscanf(fp, "%d", &(hlpod_meta->neib_id[i]));
	}
	fclose(fp);

	hlpod_meta->n_internal = BB_std_calloc_1d_int(hlpod_meta->n_internal, hlpod_meta->num_neib);
	int index_internal = 0;

	hlpod_meta->n_internal_sum = 0;
	for(int m = 0; m < hlpod_meta->num_neib; m++){
		snprintf(filename, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.n_internal.%d", hlpod_meta->neib_id[m]);
		fp = ROM_BB_read_fopen(fp, filename, directory);
		fscanf(fp, "%s %d", id, &(tmp));
		fscanf(fp, "%d", &(hlpod_meta->n_internal[m]));
		hlpod_meta->n_internal_sum += hlpod_meta->n_internal[m];
		fclose(fp);
	}

	hlpod_meta->global_id = BB_std_calloc_1d_int(hlpod_meta->global_id, hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex);

	snprintf(filename, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, filename, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(num_metagraph_nodes) , &(tmp));

	for(int i = 0; i < monolis_com->n_internal_vertex; i++) {
		fscanf(fp, "%d", &(hlpod_meta->global_id[i]));
		index_internal++;
	}
	fclose(fp);

	for (int m = 0; m < hlpod_meta->num_neib; m++){
		snprintf(filename, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", hlpod_meta->neib_id[m]);
		fp = ROM_BB_read_fopen(fp, filename, directory);
		fscanf(fp, "%s", id);
		fscanf(fp, "%d %d", &(num_metagraph_nodes) , &(tmp));

		for(int i = 0; i < hlpod_meta->n_internal[m]; i++) {
			fscanf(fp, "%d", &(hlpod_meta->global_id[index_internal]));
			index_internal++;
		}
		fclose(fp);
	}

	snprintf(filename, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, filename, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(num_metagraph_nodes) , &(tmp));
	hlpod_meta->my_global_id = BB_std_calloc_1d_int(hlpod_meta->my_global_id, num_metagraph_nodes);

	for(int i = 0; i < num_metagraph_nodes; i++) {
		fscanf(fp, "%d", &(hlpod_meta->my_global_id[i]));
		printf("%d\n", hlpod_meta->my_global_id[i]);
	}
	fclose(fp);
	/******/
}


void HROM_ddecm_set_neib(
		MONOLIS_COM*  	monolis_com,
		HLPOD_MAT* 	hlpod_mat,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_subdomains,
		const int       num_snapshots,
		const char*     directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	char id[BUFFER_SIZE];
	int tmp;
	int ndof;
	int num_2nd_subdomains;
	char fname[BUFFER_SIZE];
	char char_id[BUFFER_SIZE];
	FILE* fp;

	/*隣接関係の読み込み 別の関数にした方がよい*/ 
	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.n_internal.%d",myrank);
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(ndof));
	fscanf(fp, "%d", &(num_2nd_subdomains));		//自領域を構成するpod計算領域数
	fclose(fp);

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_2nd_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d",myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for(int i = 0; i < num_2nd_subdomains; i++){
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	BB_std_free_1d_int(subdomain_id, num_2nd_subdomains);


	char filename[BUFFER_SIZE];

	//1stddの基底本数の共有
	snprintf(filename, BUFFER_SIZE,"DDECM/n_modes_internal.%d.txt", myrank);
	fp = ROM_BB_write_fopen(fp, filename, directory);

	fprintf(fp, "%d\n", monolis_com->n_internal_vertex);
	for(int j = 0; j < monolis_com->n_internal_vertex; j++){
		fprintf(fp, "%d\n", hlpod_mat->num_modes_internal[j]);
	}
	fclose(fp);
	/**/

	double t = monolis_get_time_global_sync();

	hlpod_ddhr->num_neib_modes_1stdd = BB_std_calloc_1d_int(hlpod_ddhr->num_neib_modes_1stdd, hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex);

	snprintf(filename, BUFFER_SIZE, "DDECM/n_modes_internal.%d.txt", monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, filename, directory);

	int index_internal = 0;
	fscanf(fp, "%d",&(tmp));
	for(int i = 0; i < monolis_com->n_internal_vertex; i++) {
		fscanf(fp, "%d", &(hlpod_ddhr->num_neib_modes_1stdd[i]));
		index_internal++;
	}
	fclose(fp);

	for (int m = 0; m < hlpod_meta->num_neib; m++){
		snprintf(filename, BUFFER_SIZE, "DDECM/n_modes_internal.%d.txt", hlpod_meta->neib_id[m]);
		fp = ROM_BB_read_fopen(fp, filename, directory);

		fscanf(fp, "%d",&(tmp));
		for(int i = 0; i < hlpod_meta->n_internal[m]; i++) {
			fscanf(fp, "%d", &(hlpod_ddhr->num_neib_modes_1stdd[index_internal]));
			index_internal++;
		}
		fclose(fp);
	}

	for(int i = 0; i < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; i++) {
		printf("%d\n", hlpod_ddhr->num_neib_modes_1stdd[i]);
	}

    printf("monolis_com->n_internal_vertex = %d\n", monolis_com->n_internal_vertex);
    double t1 = monolis_get_time_global_sync();
	hlpod_ddhr->num_modes_1stdd = BB_std_calloc_1d_int(hlpod_ddhr->num_modes_1stdd, monolis_com->n_internal_vertex);

	//1stddの隣接領域を含めた総基底本数の計算
	for (int k = 0; k < monolis_com->n_internal_vertex; k++) {
		int iS = hlpod_meta->index[k];
		int iE = hlpod_meta->index[k + 1];

		for (int i = iS; i < iE; i++) {
			int item_index = hlpod_meta->item[i];
			int global_id_value = hlpod_meta->my_global_id[item_index];

			printf("global_id = %d\n", global_id_value);

			for (int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++) {
				if (global_id_value == hlpod_meta->global_id[j]) {
					hlpod_ddhr->num_modes_1stdd[k] += hlpod_ddhr->num_neib_modes_1stdd[j];
				}
			}
		}
	}

	for(int k = 0; k < monolis_com->n_internal_vertex; k++){
		printf("num_modes_1stdd = %d\n", hlpod_ddhr->num_modes_1stdd[k]);
	}

	//自領域を含めた基底本数
	for(int k = 0; k < monolis_com->n_internal_vertex; k++){
		hlpod_ddhr->num_modes_1stdd[k] += hlpod_ddhr->num_neib_modes_1stdd[k];
	}

	hlpod_ddhr->num_neib_modes_1stdd_sum = BB_std_calloc_1d_int(hlpod_ddhr->num_neib_modes_1stdd_sum, hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex + 1);
	hlpod_ddhr->num_neib_modes_1stdd_sum[0] = 0;
	for(int i = 0; i < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; i++){
		hlpod_ddhr->num_neib_modes_1stdd_sum[i + 1] = hlpod_ddhr->num_neib_modes_1stdd_sum[i] + hlpod_ddhr->num_neib_modes_1stdd[i];
	}

	hlpod_ddhr->num_internal_modes_1stdd_sum = BB_std_calloc_1d_int(hlpod_ddhr->num_internal_modes_1stdd_sum, monolis_com->n_internal_vertex + 1);
	hlpod_ddhr->num_internal_modes_1stdd_sum[0] = 0;
	for(int i = 0; i < monolis_com->n_internal_vertex; i++){
		hlpod_ddhr->num_internal_modes_1stdd_sum[i + 1] = hlpod_ddhr->num_internal_modes_1stdd_sum[i] + hlpod_ddhr->num_neib_modes_1stdd[i];
	}

    double t2 = monolis_get_time_global_sync();
}


void HROM_ddecm_set_element_para(
		BBFE_DATA*     	fe,
		HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
		const char*     directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	char id[BUFFER_SIZE];
	int tmp;
	int ndof;
	int global_id;
	int num_2nd_subdomains;
	char fname[BUFFER_SIZE];
	char char_id[BUFFER_SIZE];
	FILE* fp;

    //elem_idの読み込み
	int* local_elems_id;
	local_elems_id = BB_std_calloc_1d_int(local_elems_id , fe->total_num_elems);
	int* global_elems_id;
	global_elems_id = BB_std_calloc_1d_int(global_elems_id , fe->total_num_elems);

	snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(tmp), &(tmp));
	for(int i = 0; i < fe->total_num_elems; i++){
		fscanf(fp, "%d", &(global_elems_id[i]));	//ソート対象
	}
	fclose(fp);

	for(int i = 0; i < fe->total_num_elems; i++){
		local_elems_id[i] = i;
	}

	ROM_BB_bubble_sort_with_id(global_elems_id, local_elems_id, fe->total_num_elems);
    /**/

	hlpod_ddhr->parallel_elems_id = BB_std_calloc_1d_int(hlpod_ddhr->parallel_elems_id , fe->total_num_elems);

	snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(tmp), &(tmp));
	for(int i = 0; i < fe->total_num_elems; i++){
		fscanf(fp, "%d", &(hlpod_ddhr->parallel_elems_id[i]));	//ソート対象
	}
	fclose(fp);


    /*隣接関係の読み込み 別の関数にした方がよい*/ 
	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.n_internal.%d",myrank);
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(ndof));
	fscanf(fp, "%d", &(num_2nd_subdomains));
	fclose(fp);

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_2nd_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d",myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for(int i = 0; i < num_2nd_subdomains; i++){
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);
    /**/

	int num_elems_internal;

	snprintf(fname, BUFFER_SIZE, "parted.0/elem.dat.n_internal.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", char_id, &(tmp));
	fscanf(fp, "%d", &(num_elems_internal));
	printf("num_elems_internal = %d\n", num_elems_internal);
	fclose(fp);


	hlpod_ddhr->num_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_elems, num_subdomains);
	for(int m = 0; m < num_subdomains; m++){	
		snprintf(fname, BUFFER_SIZE, "parted.1/%s.n_internal.%d", INPUT_FILENAME_ELEM, subdomain_id[m]);

		fp = ROM_BB_read_fopen(fp, fname, directory);
		fscanf(fp, "%s %d", char_id, &(tmp));
		fscanf(fp, "%d", &(hlpod_ddhr->num_elems[m]));
	}
	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);

	hlpod_ddhr->elem_id_local = BB_std_calloc_2d_int(hlpod_ddhr->elem_id_local, max_num_elem, num_subdomains);

	int index  = 0;

	for(int m = 0; m < num_subdomains; m++){	
		snprintf(fname, BUFFER_SIZE, "parted.1/%s.%d", INPUT_FILENAME_ELEM_ID, subdomain_id[m]);

		fp = ROM_BB_read_fopen(fp, fname, directory);
		fscanf(fp, "%s", char_id);
		fscanf(fp, "%d %d", &(tmp), &(tmp));

		for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			fscanf(fp, "%d", &(global_id));
			int val = ROM_BB_binarySearch(global_elems_id, global_id, fe->total_num_elems);
			hlpod_ddhr->elem_id_local[index][m] = local_elems_id[val];
			index++;

		}
		
		fclose(fp);

		index = 0;
	}

	snprintf(fname, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_ELEM, monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", char_id, &(tmp));
	fscanf(fp, "%d", &(hlpod_ddhr->num_internal_elems));
	fclose(fp);

	hlpod_ddhr->total_num_elems = BB_std_calloc_1d_int(hlpod_ddhr->total_num_elems, 1);	//ovl要素も含んだ全要素数

	for(int m = 0; m < 1; m++){	
		snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, monolis_mpi_get_global_my_rank());

		fp = ROM_BB_read_fopen(fp, fname, directory);
		fscanf(fp, "%s", char_id);
		fscanf(fp, "%d %d", &(hlpod_ddhr->total_num_elems[m]), &(tmp));

		hlpod_ddhr->ovl_elem_global_id = BB_std_calloc_2d_int(hlpod_ddhr->ovl_elem_global_id, hlpod_ddhr->total_num_elems[0], 1);

		for(int i = 0; i < hlpod_ddhr->total_num_elems[0]; i++) {
			fscanf(fp, "%d", &(hlpod_ddhr->ovl_elem_global_id[i][m]));
		}

		fclose(fp);
	}

	BB_std_free_1d_int(subdomain_id, num_2nd_subdomains);
	BB_std_free_1d_int(local_elems_id , fe->total_num_elems);
	BB_std_free_1d_int(global_elems_id , fe->total_num_elems);

}


void HROM_ddecm_get_selected_elema_add(
	HLPOD_DDHR*     hlpod_ddhr,
	const int       num_parallel_subdomains,
	const char*     directory)
{
	double t = monolis_get_time_global_sync();

	const int myrank = monolis_mpi_get_global_my_rank();
	FILE* fp1;
	FILE* fp2;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	int val;

	int*    ovl_selected_elems;
	int*    ovl_selected_elems_D_bc;
	double* ovl_selected_elems_weight;
	double* ovl_selected_elems_weight_D_bc;

	int num_selected_elems = 0;
	int num_selected_elems_D_bc = 0;

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname1, BUFFER_SIZE,"DDECM/selected_elem.%d.txt", m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
	}

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM/selected_elem_D_bc.%d.txt", m);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(val));
		num_selected_elems_D_bc += val;
		fclose(fp2);
	}

	ovl_selected_elems = BB_std_calloc_1d_int(ovl_selected_elems, num_selected_elems);
	ovl_selected_elems_weight = BB_std_calloc_1d_double(ovl_selected_elems_weight, num_selected_elems);
	ovl_selected_elems_D_bc = BB_std_calloc_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	ovl_selected_elems_weight_D_bc = BB_std_calloc_1d_double(ovl_selected_elems_weight_D_bc, num_selected_elems_D_bc);

	int index = 0;

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname1, BUFFER_SIZE,"DDECM/selected_elem.%d.txt", m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));

		for(int j = 0; j < val; j++){
			fscanf(fp1, "%d %lf", &(ovl_selected_elems[j+index]), &(ovl_selected_elems_weight[j+index]));
		}

		index += val;

		fclose(fp1);
	}

	index = 0;
	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM/selected_elem_D_bc.%d.txt", m);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(val));

		for(int j = 0; j < val; j++){
			fscanf(fp2, "%d %lf", &(ovl_selected_elems_D_bc[j+index]), &(ovl_selected_elems_weight_D_bc[j+index]));
		}

		index += val;

		fclose(fp2);
	}

	bool*   bool_ovl_selected_elems;
	bool*   bool_ovl_selected_elems_D_bc;

	bool_ovl_selected_elems = BB_std_calloc_1d_bool(bool_ovl_selected_elems, num_selected_elems);
	bool_ovl_selected_elems_D_bc = BB_std_calloc_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);

	int*    ovl_elem_local_id;
	int*    ovl_elem_local_id_D_bc;

	ovl_elem_local_id = BB_std_calloc_1d_int(ovl_elem_local_id, num_selected_elems);
	ovl_elem_local_id_D_bc = BB_std_calloc_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);


	int index1 = 0;
	int index2 = 0;

	for(int i = 0; i < num_selected_elems; i++){
		for(int j = 0; j < hlpod_ddhr->total_num_elems[0]; j++){
			if(ovl_selected_elems[i] == hlpod_ddhr->ovl_elem_global_id[j][0]){
				bool_ovl_selected_elems[i] = true;
				ovl_elem_local_id[index1] = j;
				index1++;
			}
		}
	}
	for(int i = 0; i < num_selected_elems_D_bc; i++){
		for(int j = 0; j < hlpod_ddhr->total_num_elems[0]; j++){
			if(ovl_selected_elems_D_bc[i] == hlpod_ddhr->ovl_elem_global_id[j][0]){
				bool_ovl_selected_elems_D_bc[i] = true;
				ovl_elem_local_id_D_bc[index2] = j;
				index2++;
			}
		}
	}

	hlpod_ddhr->ovl_id_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems, index1);
	hlpod_ddhr->ovl_elem_weight = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight, index1);
	hlpod_ddhr->ovl_id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_D_bc, index2);
	hlpod_ddhr->ovl_elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_D_bc, index2);

	hlpod_ddhr->ovl_num_selected_elems = index1;
	hlpod_ddhr->ovl_num_selected_elems_D_bc = index2;

	index1 = 0;
	index2 = 0;

	for(int i = 0; i < num_selected_elems; i++){
		if(bool_ovl_selected_elems[i]){

			hlpod_ddhr->ovl_id_selected_elems[index1] = ovl_elem_local_id[index1];
			hlpod_ddhr->ovl_elem_weight[index1] = ovl_selected_elems_weight[i];
			index1++;
		}
	}

	for(int i = 0; i < num_selected_elems_D_bc; i++){
		if(bool_ovl_selected_elems_D_bc[i]){
			hlpod_ddhr->ovl_id_selected_elems_D_bc[index2] = ovl_elem_local_id_D_bc[index2];
			hlpod_ddhr->ovl_elem_weight_D_bc[index2] = ovl_selected_elems_weight_D_bc[i];
			index2++;
		}
	}

	BB_std_free_1d_int(ovl_selected_elems, num_selected_elems);
	BB_std_free_1d_double(ovl_selected_elems_weight, num_selected_elems);
	BB_std_free_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	BB_std_free_1d_double(ovl_selected_elems_weight_D_bc, num_selected_elems_D_bc);

	BB_std_free_1d_bool(bool_ovl_selected_elems, num_selected_elems);
	BB_std_free_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	BB_std_free_1d_int(ovl_elem_local_id, num_selected_elems);
	BB_std_free_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);

}

/*
void HROM_ddecm_calc_block_mat_bcsr(
	MONOLIS*     	monolis,
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT* 	hlpod_mat,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_META*		hlpod_meta,
	const int 		max_num_bases,
	const int		num_2nddd,
	const char*		directory)
{
	const int M = max_num_bases;
	const int total_num_bases = hlpod_vals->num_modes;

	const int n_neib = hlpod_vals->n_neib_vec;

	double** L_in;
	L_in = BB_std_calloc_2d_double(L_in, M , M );

	int index1 = 0;
	int index2 = 0;

	for(int k = 0; k < monolis_com->n_internal_vertex; k++){

		int iS = hlpod_ddhr->num_internal_modes_1stdd_sum[k];
		int iE = hlpod_ddhr->num_internal_modes_1stdd_sum[k+1];
		int num_modes = iE - iS;

		index1 = 0;
		for(int m = iS; m < iE; m++){
			index2 = 0;
			for(int n = iS; n < iE; n++){
				L_in[index1][index2] = hlpod_ddhr->reduced_mat[m][n];
			
				monolis_add_scalar_to_sparse_matrix_R(
					monolis,
					k,
					k,
					index1,
					index2,
					L_in[index1][index2]);
				index2++;
			}
			index1++;
		}
	}
	

	for(int k = 0; k < monolis_com->n_internal_vertex; k++){

		int iS = hlpod_meta->index[k];
		int iE = hlpod_meta->index[k + 1];

		for(int i = iS; i < iE; i++){
			for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){

				if (hlpod_meta->my_global_id[hlpod_meta->item[i]] == hlpod_meta->global_id[j]){
					
					int IS = hlpod_ddhr->num_internal_modes_1stdd_sum[k];
					int IE = hlpod_ddhr->num_internal_modes_1stdd_sum[k+1];

					index1 = 0;
					for(int m = IS; m < IE; m++){					

						int IIS = hlpod_ddhr->num_neib_modes_1stdd_sum[j];
						int IIE = hlpod_ddhr->num_neib_modes_1stdd_sum[j + 1];
						
						index2 = 0;

						for(int n = IIS; n < IIE; n++){		
							L_in[index1][index2] = hlpod_ddhr->reduced_mat[m][n];

                            //printf("val = %e", L_in[index1][index2] );
							
							monolis_add_scalar_to_sparse_matrix_R(
								monolis,
								k,
								hlpod_meta->item[i],
								index1,
								index2,
								L_in[index1][index2]);
							
							index2++;
						}
						index1++;

					}

				}

			}

		}

	}

	BB_std_free_2d_double(L_in, M, M);
}
*/

void HROM_ddecm_calc_block_mat_bcsr(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_META*     hlpod_meta,
    const int       max_num_bases,
    const int       num_2nddd,
    const char*     directory)
{
    const int M = max_num_bases;
    const int total_num_bases = hlpod_vals->num_modes;
    const int n_neib = hlpod_vals->n_neib_vec;

    int rank = monolis_mpi_get_global_my_rank();

    double** L_in;
    L_in = BB_std_calloc_2d_double(L_in, M, M);

    int index1 = 0;
    int index2 = 0;

    /* 対角ブロック */
    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        int iS = hlpod_ddhr->num_internal_modes_1stdd_sum[k];
        int iE = hlpod_ddhr->num_internal_modes_1stdd_sum[k + 1];

        double frob_sq = 0.0;

        index1 = 0;
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                double val = hlpod_ddhr->reduced_mat[m][n];
                L_in[index1][index2] = val;

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

        double frob = sqrt(frob_sq);
        printf("rank = %d [diag block] k=%d, Frobenius norm = %e\n", rank, k, frob);
    }

    /* 非対角ブロック */
    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){

                if(hlpod_meta->my_global_id[hlpod_meta->item[i]] == hlpod_meta->global_id[j]){

                    int IS = hlpod_ddhr->num_internal_modes_1stdd_sum[k];
                    int IE = hlpod_ddhr->num_internal_modes_1stdd_sum[k + 1];

                    double frob_sq = 0.0;

                    index1 = 0;
                    for(int m = IS; m < IE; m++){

                        int IIS = hlpod_ddhr->num_neib_modes_1stdd_sum[j];
                        int IIE = hlpod_ddhr->num_neib_modes_1stdd_sum[j + 1];

                        index2 = 0;
                        for(int n = IIS; n < IIE; n++){
                            double val = hlpod_ddhr->reduced_mat[m][n];
                            L_in[index1][index2] = val;

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

                    double frob = sqrt(frob_sq);
                    printf("rank = %d, [offdiag block] row=%d, col=%d, Frobenius norm = %e\n", rank,
                        k, hlpod_meta->item[i], frob);
                }
            }
        }
    }

    BB_std_free_2d_double(L_in, M, M);
}



/*for visualization*/
void ddhr_lb_set_selected_elems_para(
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
	BB_std_calloc_1d_double(ECM_elem_weight, fe->total_num_nodes);
	BB_std_calloc_1d_double(ECM_wireframe, fe->total_num_nodes);

	fclose(fp);

}


//level1領域の最大基底本数の共有
void HROM_ddecm_get_neib_max_num_modes(
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT* 	    hlpod_mat,
    const int       np,
	const int       num_my_modes)
{
	hlpod_mat->max_num_neib_modes = BB_std_calloc_1d_int(hlpod_mat->max_num_neib_modes, np);

	hlpod_mat->max_num_neib_modes[0] = num_my_modes;
	monolis_mpi_update_I(monolis_com, np, 1, hlpod_mat->max_num_neib_modes);
}


//level1領域の選択された基底(p-adaptive)本数の共有
void HROM_ddecm_get_neib_num_modes(
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT* 	    hlpod_mat,
    const int       np,
	const int       num_my_modes)
{
	hlpod_mat->num_modes_1stdd_neib = BB_std_calloc_1d_int(hlpod_mat->num_modes_1stdd_neib, np);

	hlpod_mat->num_modes_1stdd_neib[0] = num_my_modes;
	monolis_mpi_update_I(monolis_com, np, 1, hlpod_mat->num_modes_1stdd_neib);

	hlpod_mat->num_neib_modes_sum = BB_std_calloc_1d_int(hlpod_mat->num_neib_modes_sum, np);	
	hlpod_mat->num_neib_modes_sum[0] = num_my_modes;
	for(int i = 1; i < np; i++){
		hlpod_mat->num_neib_modes_sum[i] = hlpod_mat->num_neib_modes_sum[i-1] + hlpod_mat->num_modes_1stdd_neib[i];
	}
}


void HROM_ddecm_get_neib_subdomain_id(
	MONOLIS_COM*  	monolis_com,
	HLPOD_MAT* 	    hlpod_mat,
	const int 		num_modes)		//num_2nd_subdomains
{
    int n_neib_vec;
    int n_vec = num_modes;    //自領域のベクトル数
    
    monolis_mpi_get_n_neib_vector(
        monolis_com,
        n_vec,
        &n_neib_vec);       //出力：自領域と隣接領域の合計ベクトル数

    const int np = monolis_com->n_internal_vertex + monolis_com->recv_index[monolis_com->recv_n_neib]; //配列サイズ
    const int n_internal_vertex = monolis_com->n_internal_vertex;

    double** my_vec;
    my_vec = BB_std_calloc_2d_double(my_vec, np, n_vec);

    for(int i = 0; i < n_vec; i++){	
        for(int j = 0; j < n_internal_vertex; j++){
			my_vec[j][i] = hlpod_mat->subdomain_id_in_nodes_internal[j][i];
        }
    }

	double t2 = monolis_get_time_global_sync();

	BB_std_free_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, monolis_com->n_internal_vertex + monolis_com->recv_index[monolis_com->recv_n_neib], num_modes);
    double** neib_vec = BB_std_calloc_2d_double(neib_vec, np, n_neib_vec);
	hlpod_mat->subdomain_id_in_nodes = BB_std_calloc_1d_int(hlpod_mat->subdomain_id_in_nodes, np);
    const int n_dof = 1;    //計算点が持つ自由度

    monolis_mpi_get_neib_vector_R(
        monolis_com,
        np,						//配列サイズ
        n_dof,					//計算点が持つ自由度
        n_vec,					//自領域のベクトル数
        n_neib_vec,				//自領域と隣接領域の合計ベクトル数
        my_vec,					//自領域のベクトル
        neib_vec);	//自領域と隣接領域が並んだベクトル

	for(int i = 0; i < n_neib_vec; i++){
		for(int j = 0; j < np; j++){
			if(neib_vec[j][i] != 0){
				hlpod_mat->subdomain_id_in_nodes[j] = i;
			}
		}
	}

	BB_std_free_2d_double(my_vec, np, n_vec);
	BB_std_free_2d_double(neib_vec, np, n_neib_vec);
}


void HROM_ddecm_get_neib_coordinates_pre(
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int       np,				//並列計算領域数
	const int       max_num_basis)	//level 1の基底本数 (並列計算領域が担当する基底本数の総和)
{
    const int num_basis = max_num_basis * np;
	hlpod_mat->mode_coef_1stdd = BB_std_calloc_1d_double(hlpod_mat->mode_coef_1stdd, num_basis);
}


void HROM_get_neib_coordinates(
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int       np,				//並列計算領域数
	const int       max_num_basis,	//level 1の基底本数 (並列計算領域が担当する基底本数の総和)
	const int 		num_subdomains,
	const int		max_num_bases)
{
    const int num_basis = max_num_basis * np;

	int index_row = 0;
	int index_column = 0;
	int index = 0;

	for(int k = 0; k < num_subdomains; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			hlpod_mat->mode_coef_1stdd[index_column + i] = hlpod_mat->mode_coef[index + i];
		}

		index_column += hlpod_mat->num_modes_internal[k];
		index += hlpod_mat->num_modes_internal[k];
	}

	monolis_mpi_update_R(monolis_com, num_basis, max_num_basis, hlpod_mat->mode_coef_1stdd);
}


void HROM_set_bc_id(
		BBFE_DATA* 	fe,
		BBFE_BC*   	bc,
		HLPOD_DDHR* hlpod_ddhr,
		const int   num_dofs_on_node,
        HLPOD_MAT*	hlpod_mat)
{
	int nl = fe->local_num_nodes;
    int j = 0;
    hlpod_mat->hr_D_bc_node_id = BB_std_calloc_1d_int(hlpod_mat->hr_D_bc_node_id, bc->num_D_bcs);

    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];
	
		for(int i=0; i<nl; i++) {
			int index_j = fe->conn[e][i];
			if( bc->D_bc_exists[index_j] && hlpod_mat->hr_D_bc_node_id[j] == 0) {
				hlpod_mat->hr_D_bc_node_id[j] = index_j;
				j++;
			}
		}
	}

	hlpod_mat->num_hr_D_bc_nodes = j;

}


void HROM_ddecm_get_neib_subdomain_id_2nddd(
	MONOLIS_COM*  	monolis_com,
	HLPOD_MAT* 	hlpod_mat,
	const int 		num_modes)		//num_2nd_subdomains
{
    int n_neib_vec;
    int n_vec = 1;    //自領域のベクトル数
    
    monolis_mpi_get_n_neib_vector(
        monolis_com,
        n_vec,
        &n_neib_vec);       //出力：自領域と隣接領域の合計ベクトル数

    const int np = monolis_com->n_internal_vertex + monolis_com->recv_index[monolis_com->recv_n_neib]; //配列サイズ
    const int n_internal_vertex = monolis_com->n_internal_vertex;

    double** my_vec;
    my_vec = BB_std_calloc_2d_double(my_vec, np, n_vec);

    for(int i = 0; i < n_vec; i++){	
        for(int j = 0; j < n_internal_vertex; j++){
			my_vec[j][i] = 1;	//非零
        }
    }

	double t2 = monolis_get_time_global_sync();

    double** neib_vec = BB_std_calloc_2d_double(neib_vec, np, n_neib_vec);

	hlpod_mat->subdomain_id_in_nodes_2nddd = BB_std_calloc_1d_int(hlpod_mat->subdomain_id_in_nodes_2nddd, np);

    const int n_dof = 1;    //計算点が持つ自由度

    monolis_mpi_get_neib_vector_R(
        monolis_com,
        np,						//配列サイズ
        n_dof,					//計算点が持つ自由度
        n_vec,					//自領域のベクトル数
        n_neib_vec,				//自領域と隣接領域の合計ベクトル数
        my_vec,					//自領域のベクトル
        neib_vec);	//自領域と隣接領域が並んだベクトル

	for(int i = 0; i < n_neib_vec; i++){
		for(int j = 0; j < np; j++){
			if(neib_vec[j][i] != 0){
				hlpod_mat->subdomain_id_in_nodes_2nddd[j] = i;
			}
		}
	}

	BB_std_free_2d_double(my_vec, np, n_vec);
	BB_std_free_2d_double(neib_vec, np, n_neib_vec);
}

void HROM_ddecm_write_selected_elems_para_arbit_subd_svd(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();
    const int comm = monolis_mpi_get_self_comm();
    int scalapack_comm;
    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;

	int index_1 = 0;
	int index_2 = 0;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	int nl = fe->local_num_nodes;
    double t1 = monolis_get_time_global_sync();

	const int max_ITER = 50000;
	//const double TOL = 1.0e-20;

	double residual;

	double* ans_vec;
	double** matrix;
	double* RH;
	bool** bool_elem;
	int* total_id_selected_elems;
	double* total_elem_weight;
	int* total_num_selected_elems;

	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);
	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
	hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);
	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);

	int Index1 = 0;
	int Index2 = 0;

	int index_NNLS1 = 0;
	int index_NNLS2 = 0;

    //double global_norm = ddhr_calc_tol(monolis_com,
    //    hlpod_vals, hlpod_ddhr,	hlpod_mat, hlpod_meta, total_num_elem, total_num_snapshot, num_subdomains);

	for (int m = 0; m < num_subdomains; m++) {
		int NNLS_row = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot; //2は残差ベクトル＋右辺ベクトルを採用しているため
        //int NNLS_row = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot +1; //2は残差ベクトル＋右辺ベクトルを採用しているため
        printf("NNLS_row = %d num_elems = %d\n", NNLS_row, hlpod_ddhr->num_elems[m]);
		ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		RH = BB_std_calloc_1d_double(RH, NNLS_row);

		index_NNLS2 = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot;

		for (int p = 0; p < total_num_snapshot; p++) {

			for (int j = hlpod_ddhr->num_internal_modes_1stdd_sum[m] + hlpod_vals->n_neib_vec * p; j < hlpod_ddhr->num_internal_modes_1stdd_sum[m + 1] + hlpod_vals->n_neib_vec * p; j++) {
				for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
					matrix[index_NNLS1][i] = hlpod_ddhr->matrix[j][i][m];
				}
				RH[index_NNLS1] = hlpod_ddhr->RH[j][m];
				index_NNLS1++;
			}

			int iS = hlpod_meta->index[m];
			int iE = hlpod_meta->index[m + 1];
			for (int n = iS; n < iE; n++) {
				for (int l = 0; l < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; l++) {
					if (hlpod_meta->my_global_id[hlpod_meta->item[n]] == hlpod_meta->global_id[l]) {

						int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[l] + hlpod_vals->n_neib_vec * p;
						int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[l + 1] + hlpod_vals->n_neib_vec * p;

						for (int j = IS; j < IE; j++) {
							for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
								matrix[index_NNLS1][i] = hlpod_ddhr->matrix[j][i][m];
							}
							RH[index_NNLS1] = hlpod_ddhr->RH[j][m];

							index_NNLS1++;
						}

					}
				}
			}
		}

        double local_norm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
			local_norm += RH[j];
		}
	

        double local_Frovnorm = 0.0;
		for(int j = 0; j < NNLS_row; j++){
            for(int e = 0; e < hlpod_ddhr->num_elems[m]; e++){
			    local_Frovnorm += matrix[j][e];
            }
		}

        printf("local norm = %e, local_Fnorm = %e ", local_norm, local_Frovnorm);



        double input_TOL = 1.0;
	       //	TOL * sqrt(global_norm) / (num_subdomains  * sqrt(local_norm));

        double** S = BB_std_calloc_2d_double(S, NNLS_row, hlpod_ddhr->num_elems[m]);
        double* V = BB_std_calloc_1d_double(V, hlpod_ddhr->num_elems[m]);
        double** D = BB_std_calloc_2d_double(D, hlpod_ddhr->num_elems[m], hlpod_ddhr->num_elems[m]);

        double t1 = monolis_get_time();
        monolis_scalapack_gesvd_R(
                NNLS_row,
                hlpod_ddhr->num_elems[m], 
                matrix,
                S, 
                V, 
                D, 
                comm,
                scalapack_comm);
        double t2 = monolis_get_time();

        int k = ROM_BB_estimate_num_pod_modes(
            V,
            hlpod_ddhr->num_elems[m],
            10e-14);

        //k = 100;

        printf("NNLS_k = %d ", k);

        double** S_k = BB_std_calloc_2d_double(S_k, NNLS_row, k);

        for(int i = 0; i < NNLS_row; i++){
            for(int j = 0; j < k; j++){
                S_k[i][j] = S[i][j];
            }
        }

        double** G_k = BB_std_calloc_2d_double(G_k, k, hlpod_ddhr->num_elems[m]);
        double* b_k = BB_std_calloc_1d_double(b_k, k);

        ROM_BB_transposemat_mat(
            S_k,
            matrix,
            G_k,
            NNLS_row,
            k,
            hlpod_ddhr->num_elems[m]);
        
        ROM_BB_transposemat_vec(
            S_k,
            RH,
            b_k,
            NNLS_row,
            k);

            residual = 0.0;

        double TOL = 1.0e-16;

        monolis_optimize_nnls_R_with_sparse_solution(
            G_k,
            b_k,
            ans_vec, k, hlpod_ddhr->num_elems[m], max_ITER,  1.0e-14, &residual);

		index_NNLS1 = 0;
		index_NNLS2 = 0;

//		residual = 0.0;
/*
		monolis_optimize_nnls_R_with_sparse_solution(
			matrix,
			RH,
			ans_vec, NNLS_row, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);
*/
		printf("\n\nmax_iter = %d, tol = %e, residuals = %e\n\n", max_ITER, TOL, residual/sqrt(local_norm));

		int index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				index++;
			}
		}

		total_num_selected_elems[m] = index;

		printf("\n\nnum_selected_elems = %d\n\n", index);

		hr_write_NNLS_residual(residual, myrank, m, directory);
		hr_write_NNLS_num_elems(total_num_selected_elems[m], myrank, m, directory);

		total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[m]);

		index = 0;
		for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			if (ans_vec[i] != 0.0) {
				total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][m];
				total_elem_weight[index] = ans_vec[i];
				index++;
			}
		}

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			for (int i = 0; i < nl; i++) {       //六面体一次要素は8
				for (int j = 0; j < nl; j++) {
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

                    for(int k = 0; k < dof; k++) {
					    if (bc->D_bc_exists[index_j*dof + k]) {
						    bool_elem[h][m] = true;
    						hlpod_ddhr->D_bc_exists[index_j][m] = true;
	    				}
                    }
				}
			}
		}

		index = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				index++;
			}
		}

		//printf("\n\n num_elem_D_bc = %d \n\n", index);

		//index = D_bcが付与された要素数
		hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems[m] - index;
		hlpod_ddhr->num_selected_elems_D_bc[m] = index;

		int index1 = 0;
		int index2 = 0;
		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			if (bool_elem[h][m]) {
				hlpod_ddhr->id_selected_elems_D_bc[index1][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight_D_bc[index1][m] = total_elem_weight[h];

				index1++;
			}
			else {
				hlpod_ddhr->id_selected_elems[index2][m] = total_id_selected_elems[h];
				hlpod_ddhr->elem_weight[index2][m] = total_elem_weight[h];

				index2++;
			}
		}

		Index1 += index1;
		Index2 += index2;

		BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		BB_std_free_1d_double(RH, NNLS_row);

		BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[m]);

		double t = monolis_get_time_global_sync();

		FILE* fp1;
		FILE* fp2;
		char fname1[BUFFER_SIZE];
		char fname2[BUFFER_SIZE];

		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
		fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

		fprintf(fp1, "%d\n", index1);
		fprintf(fp2, "%d\n", index2);

		index_1 = 0;
		index_2 = 0;

		for (int h = 0; h < (total_num_selected_elems[m]); h++) {
			if (bool_elem[h][m]) {
				fprintf(fp1, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems_D_bc[index_1][m]], hlpod_ddhr->elem_weight_D_bc[index_1][m]);
				index_1++;
			}
			else {
				fprintf(fp2, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems[index_2][m]], hlpod_ddhr->elem_weight[index_2][m]);

				index_2++;
			}
		}

		fclose(fp1);
		fclose(fp2);

	}

	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);


	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	BB_std_free_3d_double(hlpod_ddhr->matrix, total_num_snapshot * hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
	BB_std_free_2d_double(hlpod_ddhr->RH, total_num_snapshot * hlpod_vals->n_neib_vec, num_subdomains);

	/*input
	hlpod_ddhr->matrix,
	hlpod_ddhr->RH,
	NNLS conditions
	*/

	/*output
	hlpod_ddhr->num_selected_elem
	hlpod_ddhr->id_selected_elems
	hlpod_ddhr->elem_weight
	*/

	double t = monolis_get_time_global_sync();
}

/*
 * CLEAN REPLACEMENT BLOCK - Stage 3 global-row-key + exact-route version 4
 *
 * Keep exactly one copy of this block in DDHR_para_lb.c.
 * Remove older HROM_stage3_* helper definitions and the older
 * HROM_ddecm_write_selected_elems_para_arbit_subd_hierarchical()
 * before inserting this block.
 *
 * Stage-3 MPI-side output policy:
 *   - Do NOT use elem.dat.n_internal.{rank}.
 *   - Read only parted.0/elem.dat.id.{rank}.
 *   - Save this rank's Stage-2 selected ORIGINAL-global IDs/weights.
 *   - Save this rank's raw RHS contribution.
 *   - Save raw contribution columns for ALL rank-local elements.
 *
 * The later serial Stage-3 program forms the union of Stage-2 selected
 * original-global element IDs and then resolves duplicate physical-element
 * copies across rank files.  No direct MPI_* routine is used here.
 */

/*
 * CLEAN REPLACEMENT BLOCK - Stage 3 global-row-key + exact-route version 4
 *
 * Keep exactly one copy of this block in DDHR_para_lb.c.
 * Remove older HROM_stage3_* helper definitions and the older
 * HROM_ddecm_write_selected_elems_para_arbit_subd_hierarchical()
 * before inserting this block.
 *
 * Stage-3 MPI-side output policy:
 *   - Do NOT use elem.dat.n_internal.{rank}.
 *   - Read only parted.0/elem.dat.id.{rank}.
 *   - Save this rank's Stage-2 selected ORIGINAL-global IDs/weights.
 *   - Save this rank's raw RHS contribution.
 *   - Save raw contribution columns for ALL rank-local elements.
 *
 * The later serial Stage-3 program forms the union of Stage-2 selected
 * original-global element IDs and then resolves duplicate physical-element
 * copies across rank files.  No direct MPI_* routine is used here.
 */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DDHR_para_lb.h"
#include "monolis_nnls_c.h"

#ifndef HROM_STAGE3_RANK_PARTITION_DIR
#define HROM_STAGE3_RANK_PARTITION_DIR "parted.0"
#endif

static void HROM_stage3_fail(
    const char* message,
    const int rank)
{
    fprintf(
        stderr,
        "ERROR: %s (rank=%d)\n",
        message,
        rank);

    /*
     * No direct MPI_Abort.  This Stage-3 output path contains no
     * collective communication.
     */
    exit(EXIT_FAILURE);
}

static bool HROM_stage3_read_rank_element_ids(
    const int rank,
    const int expected_rank_elem_count,
    const char* directory,
    std::vector<int>& original_global_elem_id)
{
    FILE* fp = NULL;
    char fname[BUFFER_SIZE];
    char label[BUFFER_SIZE];
    int n_elem = 0;
    int tmp = 0;

    snprintf(
        fname,
        BUFFER_SIZE,
        "%s/elem.dat.id.%d",
        HROM_STAGE3_RANK_PARTITION_DIR,
        rank);

    fp = ROM_BB_read_fopen(fp, fname, directory);

    if (fp == NULL) {
        fprintf(
            stderr,
            "WARNING: Stage3 metadata: failed to open %s (rank=%d)\n",
            fname,
            rank);
        return false;
    }

    if (fscanf(fp, "%s", label) != 1 ||
        fscanf(fp, "%d %d", &n_elem, &tmp) != 2 ||
        n_elem < 0) {

        fclose(fp);
        fprintf(
            stderr,
            "WARNING: Stage3 metadata: invalid header in %s (rank=%d)\n",
            fname,
            rank);
        return false;
    }

    original_global_elem_id.assign((size_t)n_elem, -1);

    for (int e = 0; e < n_elem; e++) {
        if (fscanf(
                fp,
                "%d",
                &original_global_elem_id[(size_t)e]) != 1) {

            fclose(fp);
            fprintf(
                stderr,
                "WARNING: Stage3 metadata: failed to read element ID %d/%d "
                "from %s (rank=%d)\n",
                e,
                n_elem,
                fname,
                rank);

            original_global_elem_id.clear();
            return false;
        }
    }

    fclose(fp);

    printf(
        "[Stage3 metadata] rank=%d, partition_dir=%s, "
        "elem.dat.id count=%d, expected rank elements=%d\n",
        rank,
        HROM_STAGE3_RANK_PARTITION_DIR,
        n_elem,
        expected_rank_elem_count);

    /*
     * rank-local element IDs are used as direct indices into this array,
     * so a count mismatch is not safe to ignore.
     */
    if (expected_rank_elem_count > 0 &&
        n_elem != expected_rank_elem_count) {

        fprintf(
            stderr,
            "WARNING: Stage3 output skipped on rank %d: "
            "elem.dat.id count=%d differs from expected rank elements=%d. "
            "Verify HROM_STAGE3_RANK_PARTITION_DIR='%s'.\n",
            rank,
            n_elem,
            expected_rank_elem_count,
            HROM_STAGE3_RANK_PARTITION_DIR);

        original_global_elem_id.clear();
        return false;
    }

    return true;
}

static void HROM_stage3_write_rank_file_no_mpi(
    MONOLIS_COM* monolis_com,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_DDHR* hlpod_ddhr,
    HLPOD_META* hlpod_meta,
    const int total_num_elem,
    const int total_num_snapshot,
    const int num_subdomains,
    const int* subdomain_id,
    const std::vector<int>& stage2_selected_parallel_elem,
    const std::vector<double>& stage2_selected_weight,
    const char* directory)
{
    const int myrank = monolis_mpi_get_global_my_rank();

    if (subdomain_id == NULL) {
        HROM_stage3_fail("subdomain_id is NULL", myrank);
    }

    if (stage2_selected_parallel_elem.size() !=
        stage2_selected_weight.size()) {
        HROM_stage3_fail(
            "selected element/weight size mismatch",
            myrank);
    }

    std::vector<int> original_global_elem_id;
    if (!HROM_stage3_read_rank_element_ids(
            myrank,
            total_num_elem,
            directory,
            original_global_elem_id)) {
        fprintf(stderr,
            "WARNING: Stage3 rank file is not written on rank %d; "
            "the normal Stage1/Stage2 HROM calculation continues.\n",
            myrank);
        return;
    }

    const int num_rank_local_elem =
        (int)original_global_elem_id.size();

    /*
     * Exact route information used by the original parallel writer:
     *   original-global ID -> parallel element ID -> fine subdomain ID.
     * Save it explicitly so the serial Stage-3 program never has to guess
     * how to return a final physical element to the legacy online files.
     */
    struct RouteEntry {
        int original_global_id;
        int parallel_elem_id;
        int fine_subdomain_id;
    };

    std::vector<RouteEntry> routes;
    routes.reserve((size_t)num_rank_local_elem);

    std::unordered_map<int, int> parallel_to_rank_local;
    parallel_to_rank_local.reserve(
        (size_t)num_rank_local_elem * 2 + 1);

    std::unordered_set<unsigned long long> route_seen;
    route_seen.reserve((size_t)num_rank_local_elem * 2 + 1);

    for (int m = 0; m < num_subdomains; m++) {
        const int num_local_elems = hlpod_ddhr->num_elems[m];

        for (int i = 0; i < num_local_elems; i++) {
            const int rank_local_elem_id =
                hlpod_ddhr->elem_id_local[i][m];

            if (rank_local_elem_id < 0 ||
                rank_local_elem_id >= num_rank_local_elem) {
                HROM_stage3_fail(
                    "rank-local element ID is out of elem.dat.id range",
                    myrank);
            }

            const int parallel_elem_id =
                hlpod_ddhr->parallel_elems_id[rank_local_elem_id];

            const auto old =
                parallel_to_rank_local.find(parallel_elem_id);

            if (old == parallel_to_rank_local.end()) {
                parallel_to_rank_local.emplace(
                    parallel_elem_id,
                    rank_local_elem_id);
            }
            else {
                const int old_local = old->second;
                if (original_global_elem_id[(size_t)old_local] !=
                    original_global_elem_id[(size_t)rank_local_elem_id]) {
                    HROM_stage3_fail(
                        "same parallel element ID maps to different original-global IDs",
                        myrank);
                }
            }

            /* de-duplicate exact (rank-local element, fine subdomain) route */
            const unsigned long long key =
                ((unsigned long long)(unsigned int)rank_local_elem_id << 32) |
                (unsigned long long)(unsigned int)subdomain_id[m];

            if (route_seen.insert(key).second) {
                RouteEntry e;
                e.original_global_id =
                    original_global_elem_id[(size_t)rank_local_elem_id];
                e.parallel_elem_id = parallel_elem_id;
                e.fine_subdomain_id = subdomain_id[m];
                routes.push_back(e);
            }
        }
    }

    /* Save the exact route table as text. */
    {
        char route_name[BUFFER_SIZE];
        snprintf(route_name, BUFFER_SIZE,
            "DDECM/stage3_route.%d.txt", myrank);

        FILE* route_fp = ROM_BB_write_fopen(
            route_fp, route_name, directory);

        fprintf(route_fp, "%d\n", (int)routes.size());
        for (const RouteEntry& e : routes) {
            fprintf(route_fp, "%d %d %d\n",
                e.original_global_id,
                e.parallel_elem_id,
                e.fine_subdomain_id);
        }
        fclose(route_fp);

        printf("[STAGE3-ROUTE] rank=%d entries=%zu file=%s\n",
            myrank, routes.size(), route_name);
    }

    /* This rank's Stage-2 selected elements in ORIGINAL-global IDs. */
    std::unordered_map<int, double> local_selected_weight_by_global;
    local_selected_weight_by_global.reserve(
        stage2_selected_parallel_elem.size() * 2 + 1);

    for (size_t k = 0; k < stage2_selected_parallel_elem.size(); k++) {
        const int parallel_elem_id =
            stage2_selected_parallel_elem[k];

        const auto it =
            parallel_to_rank_local.find(parallel_elem_id);

        if (it == parallel_to_rank_local.end()) {
            HROM_stage3_fail(
                "Stage-2 selected parallel element cannot be mapped to rank-local element",
                myrank);
        }

        const int rank_local_elem_id = it->second;
        const int original_global_id =
            original_global_elem_id[(size_t)rank_local_elem_id];
        const double weight = stage2_selected_weight[k];

        const auto old =
            local_selected_weight_by_global.find(original_global_id);

        if (old == local_selected_weight_by_global.end()) {
            local_selected_weight_by_global.emplace(
                original_global_id, weight);
        }
        else if (weight > old->second) {
            old->second = weight;
        }
    }

    std::vector<std::pair<int, double>> selected_pairs;
    selected_pairs.reserve(local_selected_weight_by_global.size());
    for (const auto& kv : local_selected_weight_by_global) {
        selected_pairs.emplace_back(kv.first, kv.second);
    }
    std::sort(selected_pairs.begin(), selected_pairs.end(),
        [](const std::pair<int,double>& a,
           const std::pair<int,double>& b) {
            return a.first < b.first;
        });

    std::vector<int64_t> local_selected_global;
    std::vector<double> local_selected_weight;
    local_selected_global.reserve(selected_pairs.size());
    local_selected_weight.reserve(selected_pairs.size());
    for (const auto& p : selected_pairs) {
        local_selected_global.push_back((int64_t)p.first);
        local_selected_weight.push_back(p.second);
    }

    const int num_output_columns = num_rank_local_elem;
    std::vector<int64_t> output_original_global_elem(
        (size_t)num_output_columns, -1);
    for (int c = 0; c < num_output_columns; c++) {
        output_original_global_elem[(size_t)c] =
            (int64_t)original_global_elem_id[(size_t)c];
    }

    const int rank_row_space =
        total_num_snapshot * hlpod_vals->n_neib_vec;

    if (rank_row_space <= 0) {
        HROM_stage3_fail(
            "invalid rank Stage-3 row space",
            myrank);
    }

    /*
     * IMPORTANT: j itself is only a rank-local array index.  It must NOT be
     * used as a cross-rank row ID.  For every j we store a globally meaningful
     * key:
     *   (snapshot, global mode-owner ID, mode offset inside that owner block).
     */
    struct RowKeyLocal {
        int snapshot;
        int64_t owner_global_id;
        int mode_offset;
    };

    std::vector<RowKeyLocal> row_key((size_t)rank_row_space);
    std::vector<unsigned char> row_key_valid(
        (size_t)rank_row_space, (unsigned char)0);
    std::vector<int> row_contribution_count(
        (size_t)rank_row_space, 0);
    std::vector<double> raw_RH(
        (size_t)rank_row_space, 0.0);
    std::vector<double> raw_matrix(
        (size_t)rank_row_space * (size_t)num_output_columns,
        0.0);

    auto register_key = [&] (
        const int j,
        const int p,
        const int64_t owner_global_id,
        const int mode_offset) {

        if (j < 0 || j >= rank_row_space) {
            HROM_stage3_fail(
                "Stage-3 rank-local row is out of range",
                myrank);
        }

        if (row_key_valid[(size_t)j] == 0) {
            row_key[(size_t)j].snapshot = p;
            row_key[(size_t)j].owner_global_id = owner_global_id;
            row_key[(size_t)j].mode_offset = mode_offset;
            row_key_valid[(size_t)j] = 1;
        }
        else {
            const RowKeyLocal& old = row_key[(size_t)j];
            if (old.snapshot != p ||
                old.owner_global_id != owner_global_id ||
                old.mode_offset != mode_offset) {
                fprintf(stderr,
                    "ERROR: one rank-local row j=%d maps to inconsistent "
                    "global mode keys on rank=%d. old=(%d,%lld,%d) "
                    "new=(%d,%lld,%d)\n",
                    j, myrank,
                    old.snapshot,
                    (long long)old.owner_global_id,
                    old.mode_offset,
                    p,
                    (long long)owner_global_id,
                    mode_offset);
                exit(EXIT_FAILURE);
            }
        }
    };

    for (int m = 0; m < num_subdomains; m++) {
        const int num_local_elems = hlpod_ddhr->num_elems[m];
        std::vector<unsigned char> row_seen_in_subdomain(
            (size_t)rank_row_space, (unsigned char)0);

        for (int p = 0; p < total_num_snapshot; p++) {
            const int JS =
                hlpod_ddhr->num_internal_modes_1stdd_sum[m]
                + hlpod_vals->n_neib_vec * p;
            const int JE =
                hlpod_ddhr->num_internal_modes_1stdd_sum[m + 1]
                + hlpod_vals->n_neib_vec * p;

            const int64_t internal_owner_global_id =
                (int64_t)subdomain_id[m];

            for (int j = JS; j < JE; j++) {
                if (row_seen_in_subdomain[(size_t)j] != 0) {
                    HROM_stage3_fail(
                        "duplicate Stage-3 internal row in one subdomain",
                        myrank);
                }
                row_seen_in_subdomain[(size_t)j] = 1;

                register_key(
                    j, p, internal_owner_global_id, j - JS);

                raw_RH[(size_t)j] += hlpod_ddhr->RH[j][m];
                row_contribution_count[(size_t)j]++;

                for (int i = 0; i < num_local_elems; i++) {
                    const int rank_local_elem_id =
                        hlpod_ddhr->elem_id_local[i][m];
                    if (rank_local_elem_id < 0 ||
                        rank_local_elem_id >= num_rank_local_elem) {
                        HROM_stage3_fail(
                            "Stage-3 rank-local element ID is out of range",
                            myrank);
                    }
                    raw_matrix[
                        (size_t)j * (size_t)num_output_columns
                        + (size_t)rank_local_elem_id]
                        += hlpod_ddhr->matrix[j][i][m];
                }
            }

            const int iS = hlpod_meta->index[m];
            const int iE = hlpod_meta->index[m + 1];

            for (int n = iS; n < iE; n++) {
                const int target_global_id =
                    hlpod_meta->my_global_id[
                        hlpod_meta->item[n]];

                bool found_mode_owner = false;

                for (int l = 0;
                     l < hlpod_meta->n_internal_sum
                         + monolis_com->n_internal_vertex;
                     l++) {

                    if (target_global_id != hlpod_meta->global_id[l]) {
                        continue;
                    }

                    found_mode_owner = true;

                    const int IS =
                        hlpod_ddhr->num_neib_modes_1stdd_sum[l]
                        + hlpod_vals->n_neib_vec * p;
                    const int IE =
                        hlpod_ddhr->num_neib_modes_1stdd_sum[l + 1]
                        + hlpod_vals->n_neib_vec * p;

                    for (int j = IS; j < IE; j++) {
                        if (row_seen_in_subdomain[(size_t)j] != 0) {
                            HROM_stage3_fail(
                                "duplicate Stage-3 neighbor row in one subdomain",
                                myrank);
                        }
                        row_seen_in_subdomain[(size_t)j] = 1;

                        register_key(
                            j, p, (int64_t)target_global_id, j - IS);

                        raw_RH[(size_t)j] += hlpod_ddhr->RH[j][m];
                        row_contribution_count[(size_t)j]++;

                        for (int i = 0; i < num_local_elems; i++) {
                            const int rank_local_elem_id =
                                hlpod_ddhr->elem_id_local[i][m];
                            if (rank_local_elem_id < 0 ||
                                rank_local_elem_id >= num_rank_local_elem) {
                                HROM_stage3_fail(
                                    "Stage-3 rank-local element ID is out of range",
                                    myrank);
                            }
                            raw_matrix[
                                (size_t)j * (size_t)num_output_columns
                                + (size_t)rank_local_elem_id]
                                += hlpod_ddhr->matrix[j][i][m];
                        }
                    }
                    break;
                }

                if (!found_mode_owner) {
                    HROM_stage3_fail(
                        "neighbor global mode ID was not found",
                        myrank);
                }
            }
        }
    }

    std::vector<int> used_rows;
    used_rows.reserve((size_t)rank_row_space);
    for (int j = 0; j < rank_row_space; j++) {
        if (row_contribution_count[(size_t)j] > 0) {
            if (row_key_valid[(size_t)j] == 0) {
                HROM_stage3_fail(
                    "used Stage-3 row has no global row key",
                    myrank);
            }
            used_rows.push_back(j);
        }
    }

    const int num_used_rows = (int)used_rows.size();
    std::vector<double> used_RH((size_t)num_used_rows, 0.0);
    std::vector<double> used_matrix(
        (size_t)num_used_rows * (size_t)num_output_columns,
        0.0);

    for (int r = 0; r < num_used_rows; r++) {
        const int j = used_rows[(size_t)r];
        used_RH[(size_t)r] = raw_RH[(size_t)j];
        for (int c = 0; c < num_output_columns; c++) {
            used_matrix[
                (size_t)r * (size_t)num_output_columns + (size_t)c]
                = raw_matrix[
                    (size_t)j * (size_t)num_output_columns + (size_t)c];
        }
    }

    char out_name[BUFFER_SIZE];
    snprintf(out_name, BUFFER_SIZE,
        "DDECM/stage3_rank.%d.bin", myrank);

    FILE* out = ROM_BB_write_fopen(out, out_name, directory);

    const char magic[8] = {'H','3','N','N','L','S','4','\0'};
    const int32_t version = 4;
    const int32_t rank32 = (int32_t)myrank;
    const int32_t nrow32 = (int32_t)num_used_rows;
    const int32_t nselected32 =
        (int32_t)local_selected_global.size();
    const int32_t nlocal32 = (int32_t)num_output_columns;

    if (fwrite(magic, sizeof(char), 8, out) != 8 ||
        fwrite(&version, sizeof(int32_t), 1, out) != 1 ||
        fwrite(&rank32, sizeof(int32_t), 1, out) != 1 ||
        fwrite(&nrow32, sizeof(int32_t), 1, out) != 1 ||
        fwrite(&nselected32, sizeof(int32_t), 1, out) != 1 ||
        fwrite(&nlocal32, sizeof(int32_t), 1, out) != 1) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 header", myrank);
    }

    for (int r = 0; r < num_used_rows; r++) {
        const int j = used_rows[(size_t)r];
        const int32_t snapshot32 =
            (int32_t)row_key[(size_t)j].snapshot;
        const int64_t owner64 =
            row_key[(size_t)j].owner_global_id;
        const int32_t offset32 =
            (int32_t)row_key[(size_t)j].mode_offset;

        if (fwrite(&snapshot32, sizeof(int32_t), 1, out) != 1 ||
            fwrite(&owner64, sizeof(int64_t), 1, out) != 1 ||
            fwrite(&offset32, sizeof(int32_t), 1, out) != 1) {
            fclose(out);
            HROM_stage3_fail("failed to write Stage-3 global row keys", myrank);
        }
    }

    if (num_used_rows > 0 &&
        fwrite(used_RH.data(), sizeof(double),
            (size_t)num_used_rows, out) != (size_t)num_used_rows) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 RHS", myrank);
    }

    if (!local_selected_global.empty() &&
        fwrite(local_selected_global.data(), sizeof(int64_t),
            local_selected_global.size(), out)
            != local_selected_global.size()) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 selected IDs", myrank);
    }

    if (!local_selected_weight.empty() &&
        fwrite(local_selected_weight.data(), sizeof(double),
            local_selected_weight.size(), out)
            != local_selected_weight.size()) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 selected weights", myrank);
    }

    if (num_output_columns > 0 &&
        fwrite(output_original_global_elem.data(), sizeof(int64_t),
            (size_t)num_output_columns, out) != (size_t)num_output_columns) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 local element IDs", myrank);
    }

    const size_t matrix_count =
        (size_t)num_used_rows * (size_t)num_output_columns;
    if (matrix_count > 0 &&
        fwrite(used_matrix.data(), sizeof(double), matrix_count, out)
            != matrix_count) {
        fclose(out);
        HROM_stage3_fail("failed to write Stage-3 matrix", myrank);
    }

    fclose(out);

    printf(
        "[STAGE3-OUTPUT-V4] rank=%d stage2_selected_local=%zu "
        "rank_local_columns=%d global_key_rows=%d file=%s\n",
        myrank,
        local_selected_global.size(),
        num_output_columns,
        num_used_rows,
        out_name);
}

void HROM_ddecm_write_selected_elems_para_arbit_subd_hierarchical(
    MONOLIS_COM*  monolis_com,
    BBFE_DATA*    fe,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_DDHR*   hlpod_ddhr,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_META*   hlpod_meta,
    const int     total_num_elem,
    const int     total_num_snapshot,
    const int     total_num_modes,
    const int     num_subdomains,
    const int     max_iter,
    const double  tol,
    const int     dof,
    const char*   directory)
{
    const int myrank = monolis_mpi_get_global_my_rank();

    //exit(1);
    printf(" hierarchical\n\n");

    /*
     * max_iter と tol は引数を使用する。
     * 不正値が渡された場合のみデフォルト値を使用する。
     */
    /*
     * Avoid pathological add/remove cycling from consuming thousands
     * of iterations.  The caller may request fewer than 1000, but never
     * more than 1000 in this hierarchical routine.
     */
    const int nnls_max_iter =
        (max_iter > 0) ? std::min(max_iter, 1000) : 1000;

    const double nnls_tol =
        (tol > 0.0) ? tol : 1.0e-10;

    /*
     * NNLS 解を選択要素とみなす閾値。
     * 必要に応じて nnls_tol とは別の値にしてもよい。
     */
    const double weight_tol = nnls_tol;

    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int ndof;
    int* subdomain_id = NULL;

    double start_time =
        monolis_get_time_global_sync();

    /*
     * この値は「1要素を構成する節点数」であることを前提とする。
     *
     * 元コードでは fe->local_num_nodes が使われているが、
     * これが全局所節点数を意味する場合は、
     * 要素種別に応じた節点数へ置き換える必要がある。
     */
    const int num_nodes_per_elem =
        fe->local_num_nodes;

    /*
     * 使用しない引数に対する警告抑制。
     * 必要になれば削除する。
     */
    (void)hlpod_mat;
    (void)total_num_modes;

    /*
     * 部分領域IDを読み込む。
     */
    subdomain_id =
        (int*)calloc((size_t)num_subdomains, sizeof(int));

    if (subdomain_id == NULL) {
        fprintf(
            stderr,
            "ERROR: failed to allocate subdomain_id\n");
        return;
    }

    snprintf(
        fname,
        BUFFER_SIZE,
        "metagraph_parted.0/metagraph.dat.id.%d",
        myrank);

    fp = ROM_BB_read_fopen(fp, fname, directory);

    if (fp == NULL) {
        fprintf(
            stderr,
            "ERROR: failed to open %s\n",
            fname);

        free(subdomain_id);
        return;
    }

    if (fscanf(fp, "%s", id) != 1) {
        fprintf(
            stderr,
            "ERROR: failed to read ID from %s\n",
            fname);

        fclose(fp);
        free(subdomain_id);
        return;
    }

    /*
     * 元コードの形式を維持する。
     */
    if (fscanf(fp, "%d %d", &ndof, &ndof) != 2) {
        fprintf(
            stderr,
            "ERROR: failed to read header from %s\n",
            fname);

        fclose(fp);
        free(subdomain_id);
        return;
    }

    for (int m = 0; m < num_subdomains; m++) {
        if (fscanf(fp, "%d", &subdomain_id[m]) != 1) {
            fprintf(
                stderr,
                "ERROR: failed to read subdomain ID: m = %d\n",
                m);

            fclose(fp);
            free(subdomain_id);
            return;
        }
    }

    fclose(fp);

    /*
     * parallel_elems_id は連続した配列添字とは限らない。
     * そのため、任意のグローバル要素IDをキーとして候補集合を保持する。
     */
    std::unordered_set<int> unique_element_ids;
    std::unordered_set<int> candidate_ids;
    std::unordered_set<int> initial_active_ids;

    size_t max_local_elements = 0;

    for (int m = 0; m < num_subdomains; m++) {
        max_local_elements +=
            (size_t)hlpod_ddhr->num_elems[m];
    }

    unique_element_ids.reserve(max_local_elements);
    candidate_ids.reserve(max_local_elements);
    initial_active_ids.reserve((size_t)num_subdomains);

    /*
     * 全領域に存在する要素コピーを、グローバル要素IDで
     * 重複排除し、削減率の基準となる物理要素集合を作る。
     */
    for (int m = 0; m < num_subdomains; m++) {
        const int num_local_elems =
            hlpod_ddhr->num_elems[m];

        for (int i = 0; i < num_local_elems; i++) {
            const int local_elem_id =
                hlpod_ddhr->elem_id_local[i][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            if (global_elem_id < 0) {
                fprintf(
                    stderr,
                    "ERROR: negative global element ID while "
                    "building unique element set: "
                    "rank=%d subdomain=%d local_col=%d "
                    "local_elem=%d global_elem=%d\n",
                    myrank,
                    m,
                    i,
                    local_elem_id,
                    global_elem_id);

                exit(EXIT_FAILURE);
            }

            unique_element_ids.insert(global_elem_id);
        }
    }

    int num_domains_with_candidates = 0;
    int num_domains_with_initial_active = 0;

    /*
     * Persist Stage-1 per-fine-subdomain statistics so that the later
     * serial Stage-3 executable can print one combined report over all
     * MPI ranks and all fine subdomains.
     *
     * Candidate IDs stored here are the legacy parallel element IDs.
     * The serial program maps them to original-global IDs using the
     * exact stage3_route.<rank>.txt table.
     */
    std::vector<int> stage1_local_unique_count(
        (size_t)num_subdomains, 0);
    std::vector<int> stage1_local_candidate_count(
        (size_t)num_subdomains, 0);
    std::vector<int> stage1_local_initial_active(
        (size_t)num_subdomains, 0);
    std::vector<double> stage1_local_residual(
        (size_t)num_subdomains, 0.0);
    std::vector<std::vector<int>> stage1_local_candidate_ids(
        (size_t)num_subdomains);

    size_t sum_local_unique_elements = 0;
    size_t sum_local_candidate_elements = 0;

    printf(
        "\n"
        "============================================================\n"
        "[NNLS element set before two-stage reduction]\n"
        "rank                              = %d\n"
        "number of subdomains              = %d\n"
        "input total_num_elem              = %d\n"
        "element copies over subdomains    = %zu\n"
        "unique physical elements          = %zu\n"
        "duplicate element copies          = %zu\n"
        "============================================================\n",
        myrank,
        num_subdomains,
        total_num_elem,
        max_local_elements,
        unique_element_ids.size(),
        (max_local_elements >= unique_element_ids.size())
            ? max_local_elements - unique_element_ids.size()
            : 0);

    /*
     * 各部分領域の正規化係数。
     *
     * 第1パスの局所NNLSと、第2パスの全体NNLSで
     * 同じ正規化を再現するために保存する。
     */
    double* subdomain_scale =
        (double*)calloc(
            (size_t)num_subdomains,
            sizeof(double));

    /*
     * 各領域のNNLS行数。
     */
    int* subdomain_num_rows =
        (int*)calloc(
            (size_t)num_subdomains,
            sizeof(int));

    if (subdomain_scale == NULL ||
        subdomain_num_rows == NULL) {

        fprintf(
            stderr,
            "ERROR: failed to allocate candidate arrays\n");

        free(subdomain_scale);
        free(subdomain_num_rows);
        free(subdomain_id);

        return;
    }

    /*
     * ============================================================
     * 第1パス
     *
     * 各領域で局所 sparse NNLS を解き、
     *
     * 1. 選択された全要素を候補集合へ登録
     * 2. 最大重み要素を初期 active set へ登録
     *
     * する。
     * ============================================================
     */

    /*
     * hlpod_ddhr->matrix と RH の第1添字 j が属する、
     * スナップショット込みの共通グローバル・モード行空間。
     *
     * 第1段階の NNLS_row は各領域が使用する行数であり、
     * 第2段階の共通行番号そのものではない。
     */
    const int global_NNLS_row =
        total_num_snapshot * hlpod_vals->n_neib_vec;

    if (global_NNLS_row <= 0) {
        fprintf(
            stderr,
            "ERROR: invalid global mode row count: "
            "rank=%d snapshots=%d n_neib_vec=%d rows=%d\n",
            myrank,
            total_num_snapshot,
            hlpod_vals->n_neib_vec,
            global_NNLS_row);

        exit(EXIT_FAILURE);
    }

    for (int m = 0; m < num_subdomains; m++) {
        const int NNLS_row =
            hlpod_ddhr->num_modes_1stdd[m]
            * total_num_snapshot;

        const int num_local_elems =
            hlpod_ddhr->num_elems[m];

        std::unordered_set<int> local_unique_ids;
        std::unordered_set<int> local_candidate_ids;

        local_unique_ids.reserve((size_t)num_local_elems);
        local_candidate_ids.reserve((size_t)num_local_elems);

        for (int i = 0; i < num_local_elems; i++) {
            const int local_elem_id =
                hlpod_ddhr->elem_id_local[i][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            if (global_elem_id < 0) {
                fprintf(
                    stderr,
                    "ERROR: negative global element ID in local set: "
                    "rank=%d subdomain=%d local_col=%d "
                    "local_elem=%d global_elem=%d\n",
                    myrank,
                    m,
                    i,
                    local_elem_id,
                    global_elem_id);

                exit(EXIT_FAILURE);
            }

            local_unique_ids.insert(global_elem_id);
        }

        subdomain_num_rows[m] = NNLS_row;

        /*
         * NNLS_row は、この領域で第1段階に使用する方程式数。
         * 第2段階では各方程式を元の行番号 j へ戻して組み立てる。
         */

        printf(
            "rank = %d, subdomain = %d, "
            "NNLS_row = %d, num_elems = %d\n",
            myrank,
            m,
            NNLS_row,
            num_local_elems);

        double* ans_vec = NULL;
        double** matrix = NULL;
        double* RH = NULL;

        ans_vec =
            BB_std_calloc_1d_double(
                ans_vec,
                num_local_elems);

        matrix =
            BB_std_calloc_2d_double(
                matrix,
                NNLS_row,
                num_local_elems);

        RH =
            BB_std_calloc_1d_double(
                RH,
                NNLS_row);

        int local_row = 0;

        /*
         * 現在のモード・スナップショットループを
         * そのまま使用して局所NNLS行列を構築する。
         */
        for (int p = 0;
             p < total_num_snapshot;
             p++) {

            /*
             * 内部モード。
             */
            const int JS =
                hlpod_ddhr
                    ->num_internal_modes_1stdd_sum[m]
                + hlpod_vals->n_neib_vec * p;

            const int JE =
                hlpod_ddhr
                    ->num_internal_modes_1stdd_sum[m + 1]
                + hlpod_vals->n_neib_vec * p;

            for (int j = JS; j < JE; j++) {
                if (local_row >= NNLS_row) {
                    fprintf(
                        stderr,
                        "ERROR: NNLS row overflow: "
                        "rank=%d subdomain=%d "
                        "row=%d expected=%d\n",
                        myrank,
                        m,
                        local_row,
                        NNLS_row);

                    exit(EXIT_FAILURE);
                }

                for (int i = 0;
                     i < num_local_elems;
                     i++) {

                    matrix[local_row][i] =
                        hlpod_ddhr->matrix[j][i][m];
                }

                RH[local_row] =
                    hlpod_ddhr->RH[j][m];

                local_row++;
            }

            /*
             * 近傍モード。
             */
            const int iS =
                hlpod_meta->index[m];

            const int iE =
                hlpod_meta->index[m + 1];

            for (int n = iS; n < iE; n++) {
                const int target_global_id =
                    hlpod_meta->my_global_id[
                        hlpod_meta->item[n]];

                for (int l = 0;
                     l <
                         hlpod_meta->n_internal_sum
                         + monolis_com->n_internal_vertex;
                     l++) {

                    if (target_global_id
                        != hlpod_meta->global_id[l]) {
                        continue;
                    }

                    const int IS =
                        hlpod_ddhr
                            ->num_neib_modes_1stdd_sum[l]
                        + hlpod_vals->n_neib_vec * p;

                    const int IE =
                        hlpod_ddhr
                            ->num_neib_modes_1stdd_sum[l + 1]
                        + hlpod_vals->n_neib_vec * p;

                    for (int j = IS; j < IE; j++) {
                        if (local_row >= NNLS_row) {
                            fprintf(
                                stderr,
                                "ERROR: NNLS row overflow: "
                                "rank=%d subdomain=%d "
                                "row=%d expected=%d\n",
                                myrank,
                                m,
                                local_row,
                                NNLS_row);

                            exit(EXIT_FAILURE);
                        }

                        for (int i = 0;
                             i < num_local_elems;
                             i++) {

                            matrix[local_row][i] =
                                hlpod_ddhr
                                    ->matrix[j][i][m];
                        }

                        RH[local_row] =
                            hlpod_ddhr->RH[j][m];

                        local_row++;
                    }

                    /*
                     * target_global_id に一致する l は
                     * 一意であることを前提とする。
                     */
                    break;
                }
            }
        }

        if (local_row != NNLS_row) {
            fprintf(
                stderr,
                "ERROR: NNLS row mismatch: "
                "rank=%d subdomain=%d "
                "actual=%d expected=%d\n",
                myrank,
                m,
                local_row,
                NNLS_row);

            exit(EXIT_FAILURE);
        }

        /*
         * Frobeniusノルムを計算する。
         *
         * 元コードの単純和
         *
         *   local_Frovnorm += matrix[j][e]
         *
         * では正負が打ち消し合うため、
         * 二乗和の平方根を使用する。
         */
        double local_Frovnorm = 0.0;

        for (int j = 0; j < NNLS_row; j++) {
            for (int i = 0;
                 i < num_local_elems;
                 i++) {

                const double value =
                    matrix[j][i];

                local_Frovnorm +=
                    value * value;
            }
        }

        local_Frovnorm =
            std::sqrt(local_Frovnorm);

        if (local_Frovnorm <= DBL_MIN) {
            /*
             * ゼロ行列の場合は正規化しない。
             */
            local_Frovnorm = 1.0;
        }

        subdomain_scale[m] =
            local_Frovnorm;

        for (int j = 0; j < NNLS_row; j++) {
            RH[j] /= local_Frovnorm;

            for (int i = 0;
                 i < num_local_elems;
                 i++) {

                matrix[j][i] /=
                    local_Frovnorm;
            }
        }

        double local_residual = 0.0;

        monolis_optimize_nnls_R_with_sparse_solution(
            matrix,
            RH,
            ans_vec,
            NNLS_row,
            num_local_elems,
            nnls_max_iter,
            nnls_tol,
            &local_residual);

        printf(
            "rank = %d, subdomain = %d, "
            "local residual = %.15e\n",
            myrank,
            m,
            local_residual);

        /*
         * この領域で最大重みとなった局所列。
         */
        int max_local_col = -1;
        double max_weight = -1.0;


        for (int i = 0;
             i < num_local_elems;
             i++) {

            if (ans_vec[i] <= weight_tol) {
                continue;
            }

            const int local_elem_id =
                hlpod_ddhr
                    ->elem_id_local[i][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            if (global_elem_id < 0) {
                fprintf(
                    stderr,
                    "ERROR: negative global element ID: "
                    "rank=%d subdomain=%d local_col=%d "
                    "local_elem=%d global_elem=%d\n",
                    myrank,
                    m,
                    i,
                    local_elem_id,
                    global_elem_id);

                exit(EXIT_FAILURE);
            }

            /*
             * 各領域で選択された要素の和集合。
             * global_elem_id は配列添字ではなくキーとして扱う。
             */
            candidate_ids.insert(global_elem_id);
            local_candidate_ids.insert(global_elem_id);

            if (ans_vec[i] > max_weight) {
                max_weight = ans_vec[i];
                max_local_col = i;
            }
        }

        /*
         * 各領域で最大重みの要素を初期 active にする。
         */
        if (max_local_col >= 0) {
            const int local_elem_id =
                hlpod_ddhr
                    ->elem_id_local[max_local_col][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            candidate_ids.insert(global_elem_id);
            local_candidate_ids.insert(global_elem_id);
            initial_active_ids.insert(global_elem_id);
            num_domains_with_initial_active++;

            printf(
                "rank = %d, subdomain = %d, "
                "initial active global element = %d, "
                "weight = %.15e\n",
                myrank,
                m,
                global_elem_id,
                max_weight);
        }
        else {
            printf(
                "WARNING: no positive NNLS coefficient: "
                "rank=%d subdomain=%d\n",
                myrank,
                m);
        }

        const size_t local_unique_count =
            local_unique_ids.size();

        const size_t local_candidate_count =
            local_candidate_ids.size();

        stage1_local_unique_count[(size_t)m] =
            (int)local_unique_count;
        stage1_local_candidate_count[(size_t)m] =
            (int)local_candidate_count;
        stage1_local_initial_active[(size_t)m] =
            (max_local_col >= 0) ? 1 : 0;
        stage1_local_residual[(size_t)m] =
            local_residual;

        stage1_local_candidate_ids[(size_t)m].assign(
            local_candidate_ids.begin(),
            local_candidate_ids.end());
        std::sort(
            stage1_local_candidate_ids[(size_t)m].begin(),
            stage1_local_candidate_ids[(size_t)m].end());

        sum_local_unique_elements +=
            local_unique_count;

        sum_local_candidate_elements +=
            local_candidate_count;

        if (local_candidate_count > 0) {
            num_domains_with_candidates++;
        }

        const size_t local_removed =
            (local_unique_count >= local_candidate_count)
            ? local_unique_count - local_candidate_count
            : 0;

        const double local_keep_ratio =
            (local_unique_count > 0)
            ? 100.0
                * (double)local_candidate_count
                / (double)local_unique_count
            : 0.0;

        const double local_reduction_ratio =
            (local_unique_count > 0)
            ? 100.0
                * (double)local_removed
                / (double)local_unique_count
            : 0.0;

        printf(
            "\n"
            "[Stage 1: subdomain reduction]\n"
            "rank                       = %d\n"
            "subdomain                  = %d / %d\n"
            "subdomain ID               = %d\n"
            "local element copies       = %d\n"
            "local unique elements      = %zu\n"
            "local candidate elements   = %zu\n"
            "local initial active       = %d\n"
            "local removed              = %zu\n"
            "local retained             = %.3f %%\n"
            "local reduction            = %.3f %%\n"
            "local residual             = %.15e\n",
            myrank,
            m + 1,
            num_subdomains,
            subdomain_id[m],
            num_local_elems,
            local_unique_count,
            local_candidate_count,
            (max_local_col >= 0) ? 1 : 0,
            local_removed,
            local_keep_ratio,
            local_reduction_ratio,
            local_residual);

        BB_std_free_1d_double(
            ans_vec,
            num_local_elems);

        BB_std_free_2d_double(
            matrix,
            NNLS_row,
            num_local_elems);

        BB_std_free_1d_double(
            RH,
            NNLS_row);
    }

    /*
     * ============================================================
     * 候補集合を縮約列番号へ変換する。
     * ============================================================
     */

    std::vector<int> candidate_elem(
        candidate_ids.begin(),
        candidate_ids.end());

    /*
     * unordered_set の反復順序は不定なので、
     * グローバルID順に並べて候補列順序を再現可能にする。
     */
    std::sort(
        candidate_elem.begin(),
        candidate_elem.end());

    const int num_candidates =
        (int)candidate_elem.size();

    const int num_initial_active =
        (int)initial_active_ids.size();

    const int num_unique_elements =
        (int)unique_element_ids.size();

    const int num_removed_stage1 =
        num_unique_elements - num_candidates;

    const size_t duplicate_candidates_across_domains =
        (sum_local_candidate_elements >= (size_t)num_candidates)
        ? sum_local_candidate_elements - (size_t)num_candidates
        : 0;

    const int duplicate_initial_active =
        (num_domains_with_initial_active >= num_initial_active)
        ? num_domains_with_initial_active - num_initial_active
        : 0;

    const double stage1_keep_ratio =
        (num_unique_elements > 0)
        ? 100.0
            * (double)num_candidates
            / (double)num_unique_elements
        : 0.0;

    const double stage1_reduction_ratio =
        (num_unique_elements > 0)
        ? 100.0
            * (double)num_removed_stage1
            / (double)num_unique_elements
        : 0.0;

    const double average_candidates_per_domain =
        (num_subdomains > 0)
        ? (double)sum_local_candidate_elements
            / (double)num_subdomains
        : 0.0;

    const double active_domain_coverage =
        (num_subdomains > 0)
        ? 100.0
            * (double)num_domains_with_initial_active
            / (double)num_subdomains
        : 0.0;

    printf(
        "\n"
        "============================================================\n"
        "[Stage 1 summary]\n"
        "rank                              = %d\n"
        "number of subdomains              = %d\n"
        "domains with candidates           = %d\n"
        "domains with initial active       = %d\n"
        "active-domain coverage            = %.3f %%\n"
        "global NNLS rows                  = %d\n"
        "element copies over subdomains    = %zu\n"
        "sum of local unique counts        = %zu\n"
        "unique physical elements          = %d\n"
        "sum of local candidate counts     = %zu\n"
        "unique stage-1 candidates         = %d\n"
        "candidate duplicates              = %zu\n"
        "average candidates per domain     = %.3f\n"
        "initial-active selections         = %d\n"
        "unique initial-active elements    = %d\n"
        "duplicate initial-active picks    = %d\n"
        "stage-1 removed                   = %d\n"
        "stage-1 retained                  = %.3f %%\n"
        "stage-1 reduction                 = %.3f %%\n"
        "============================================================\n",
        myrank,
        num_subdomains,
        num_domains_with_candidates,
        num_domains_with_initial_active,
        active_domain_coverage,
        global_NNLS_row,
        max_local_elements,
        sum_local_unique_elements,
        num_unique_elements,
        sum_local_candidate_elements,
        num_candidates,
        duplicate_candidates_across_domains,
        average_candidates_per_domain,
        num_domains_with_initial_active,
        num_initial_active,
        duplicate_initial_active,
        num_removed_stage1,
        stage1_keep_ratio,
        stage1_reduction_ratio);


    /*
     * Write Stage-1 statistics for the serial/global reduction report.
     *
     * File format:
     *   HROM_STAGE1_STATS_V1
     *   <num_subdomains>
     *   <fine_subdomain_id> <local_unique> <local_candidates>
     *       <initial_active> <local_residual>
     *   <candidate_parallel_id 0>
     *   ...
     */
    {
        char stage1_stats_name[BUFFER_SIZE];
        snprintf(
            stage1_stats_name,
            BUFFER_SIZE,
            "DDECM/stage1_stats.%d.txt",
            myrank);

        FILE* fp_stage1_stats =
            ROM_BB_write_fopen(
                fp_stage1_stats,
                stage1_stats_name,
                directory);

        if (fp_stage1_stats == NULL) {
            fprintf(
                stderr,
                "ERROR: failed to open Stage-1 stats output on rank %d\n",
                myrank);
            exit(EXIT_FAILURE);
        }

        fprintf(fp_stage1_stats, "HROM_STAGE1_STATS_V1\n");
        fprintf(fp_stage1_stats, "%d\n", num_subdomains);

        for (int m = 0; m < num_subdomains; m++) {
            fprintf(
                fp_stage1_stats,
                "%d %d %d %d %.30e\n",
                subdomain_id[m],
                stage1_local_unique_count[(size_t)m],
                stage1_local_candidate_count[(size_t)m],
                stage1_local_initial_active[(size_t)m],
                stage1_local_residual[(size_t)m]);

            const std::vector<int>& ids =
                stage1_local_candidate_ids[(size_t)m];

            for (size_t k = 0; k < ids.size(); k++) {
                fprintf(fp_stage1_stats, "%d\n", ids[k]);
            }
        }

        fclose(fp_stage1_stats);

        printf(
            "[STAGE1-STATS] rank=%d subdomains=%d "
            "file=%s\n",
            myrank,
            num_subdomains,
            stage1_stats_name);
    }

    if (num_candidates == 0) {
        fprintf(
            stderr,
            "WARNING: no candidate elements were selected "
            "on rank %d\n",
            myrank);

        /*
         * 出力配列を空の状態で確保する。
         */
        const int selection_capacity =
            (nnls_max_iter > 0)
            ? nnls_max_iter
            : 1;

        hlpod_ddhr->D_bc_exists =
            BB_std_calloc_2d_bool(
                hlpod_ddhr->D_bc_exists,
                fe->total_num_nodes,
                num_subdomains);

        hlpod_ddhr->id_selected_elems =
            BB_std_calloc_2d_int(
                hlpod_ddhr->id_selected_elems,
                selection_capacity,
                num_subdomains);

        hlpod_ddhr->id_selected_elems_D_bc =
            BB_std_calloc_2d_int(
                hlpod_ddhr->id_selected_elems_D_bc,
                selection_capacity,
                num_subdomains);

        hlpod_ddhr->elem_weight =
            BB_std_calloc_2d_double(
                hlpod_ddhr->elem_weight,
                selection_capacity,
                num_subdomains);

        hlpod_ddhr->elem_weight_D_bc =
            BB_std_calloc_2d_double(
                hlpod_ddhr->elem_weight_D_bc,
                selection_capacity,
                num_subdomains);

        hlpod_ddhr->num_selected_elems =
            BB_std_calloc_1d_int(
                hlpod_ddhr->num_selected_elems,
                num_subdomains);

        hlpod_ddhr->num_selected_elems_D_bc =
            BB_std_calloc_1d_int(
                hlpod_ddhr->num_selected_elems_D_bc,
                num_subdomains);

        for (int m = 0;
             m < num_subdomains;
             m++) {

            hr_write_NNLS_residual(
                0.0,
                myrank,
                m,
                directory);

            hr_write_NNLS_num_elems(
                0,
                myrank,
                m,
                directory);

            char fname_D_bc[BUFFER_SIZE];
            char fname_no_D_bc[BUFFER_SIZE];

            snprintf(
                fname_D_bc,
                BUFFER_SIZE,
                "DDECM/lb_selected_elem_D_bc.%d.txt",
                subdomain_id[m]);

            snprintf(
                fname_no_D_bc,
                BUFFER_SIZE,
                "DDECM/lb_selected_elem.%d.txt",
                subdomain_id[m]);

            FILE* fp_D_bc =
                ROM_BB_write_fopen(
                    fp_D_bc,
                    fname_D_bc,
                    directory);

            FILE* fp_no_D_bc =
                ROM_BB_write_fopen(
                    fp_no_D_bc,
                    fname_no_D_bc,
                    directory);

            fprintf(fp_D_bc, "0\n");
            fprintf(fp_no_D_bc, "0\n");

            fclose(fp_D_bc);
            fclose(fp_no_D_bc);
        }

        /*
         * Stage 3 output must also be written by a rank with no Stage-2
         * selections, because this rank may own an element selected as a
         * ghost on another rank. No MPI collective is used here.
         */
        {
            const std::vector<int> empty_selected_elem;
            const std::vector<double> empty_selected_weight;

            HROM_stage3_write_rank_file_no_mpi(
                monolis_com,
                hlpod_vals,
                hlpod_ddhr,
                hlpod_meta,
                total_num_elem,
                total_num_snapshot,
                num_subdomains,
                subdomain_id,
                empty_selected_elem,
                empty_selected_weight,
                directory);
        }

        free(subdomain_scale);
        free(subdomain_num_rows);
        free(subdomain_id);

        const int max_num_elem_empty =
            ROM_BB_findMax(
                hlpod_ddhr->num_elems,
                num_subdomains);

        BB_std_free_3d_double(
            hlpod_ddhr->matrix,
            total_num_snapshot
                * hlpod_vals->n_neib_vec,
            max_num_elem_empty,
            num_subdomains);

        BB_std_free_2d_double(
            hlpod_ddhr->RH,
            total_num_snapshot
                * hlpod_vals->n_neib_vec,
            num_subdomains);

        return;
    }

    /*
     * 任意のグローバル要素IDを、
     * 0 ～ num_candidates - 1 の候補列番号へ変換する。
     */
    std::unordered_map<int, int> candidate_col;
    candidate_col.reserve(candidate_elem.size() * 2 + 1);

    std::vector<int> active_set_init(
        (size_t)num_candidates,
        0);

    for (int c = 0; c < num_candidates; c++) {
        const int global_elem_id =
            candidate_elem[(size_t)c];

        candidate_col.emplace(
            global_elem_id,
            c);

        if (initial_active_ids.find(global_elem_id)
            != initial_active_ids.end()) {

            active_set_init[(size_t)c] = 1;
        }
    }

    /*
     * ============================================================
     * 第2パス
     *
     * 第1段階で選択された要素の和集合だけを列として使用する。
     *
     * 行番号には各領域で0から振り直した local_row を使わず、
     * hlpod_ddhr->matrix[j][i][m] と RH[j][m] の元の行番号 j を使う。
     * これにより、異なる領域にある同じグローバル・モード行の寄与を
     * 正しい同一行へ加算する。
     *
     * full_matrix の大きさ：
     *   行数 = total_num_snapshot * hlpod_vals->n_neib_vec
     *   列数 = 第1段階 active 要素の和集合数
     *
     * 最後に、実際に寄与が存在する行だけへ圧縮してNNLSへ渡す。
     * ============================================================
     */

    double** full_matrix = NULL;
    double* full_RH = NULL;

    full_matrix =
        BB_std_calloc_2d_double(
            full_matrix,
            global_NNLS_row,
            num_candidates);

    full_RH =
        BB_std_calloc_1d_double(
            full_RH,
            global_NNLS_row);

    std::vector<int> row_contribution_count(
        (size_t)global_NNLS_row,
        0);

    std::vector<int> column_contribution_count(
        (size_t)num_candidates,
        0);

    for (int m = 0;
         m < num_subdomains;
         m++) {

        const int NNLS_row =
            subdomain_num_rows[m];

        const int num_local_elems =
            hlpod_ddhr->num_elems[m];

        /*
         * 同一領域内で同じ j を複数回登録すると、RHを二重加算する。
         * そのような行対応の重複を検出する。
         */
        std::vector<unsigned char> row_seen_in_subdomain(
            (size_t)global_NNLS_row,
            (unsigned char)0);

        int local_equation_count = 0;

        for (int p = 0;
             p < total_num_snapshot;
             p++) {

            /* 内部モード。 */
            const int JS =
                hlpod_ddhr
                    ->num_internal_modes_1stdd_sum[m]
                + hlpod_vals->n_neib_vec * p;

            const int JE =
                hlpod_ddhr
                    ->num_internal_modes_1stdd_sum[m + 1]
                + hlpod_vals->n_neib_vec * p;

            for (int j = JS; j < JE; j++) {
                if (j < 0 || j >= global_NNLS_row) {
                    fprintf(
                        stderr,
                        "ERROR: stage-2 internal mode row is out of range: "
                        "rank=%d subdomain=%d snapshot=%d j=%d rows=%d\n",
                        myrank,
                        m,
                        p,
                        j,
                        global_NNLS_row);

                    exit(EXIT_FAILURE);
                }

                if (row_seen_in_subdomain[(size_t)j] != 0) {
                    fprintf(
                        stderr,
                        "ERROR: duplicate stage-2 row in one subdomain: "
                        "rank=%d subdomain=%d snapshot=%d j=%d\n",
                        myrank,
                        m,
                        p,
                        j);

                    exit(EXIT_FAILURE);
                }

                row_seen_in_subdomain[(size_t)j] = 1;
                local_equation_count++;

                for (int i = 0;
                     i < num_local_elems;
                     i++) {

                    const int local_elem_id =
                        hlpod_ddhr
                            ->elem_id_local[i][m];

                    const int global_elem_id =
                        hlpod_ddhr
                            ->parallel_elems_id[
                                local_elem_id];

                    const auto candidate_it =
                        candidate_col.find(global_elem_id);

                    /*
                     * 第1段階でactiveでなかった要素の寄与は、
                     * 第2段階には入れない。
                     */
                    if (candidate_it == candidate_col.end()) {
                        continue;
                    }

                    const int c =
                        candidate_it->second;

                    const double contribution =
                        hlpod_ddhr->matrix[j][i][m];

                    full_matrix[j][c] +=
                        contribution;

                    if (contribution != 0.0) {
                        column_contribution_count[(size_t)c]++;
                    }
                }

                full_RH[j] +=
                    hlpod_ddhr->RH[j][m];

                row_contribution_count[(size_t)j]++;
            }

            /* 近傍モード。 */
            const int iS =
                hlpod_meta->index[m];

            const int iE =
                hlpod_meta->index[m + 1];

            for (int n = iS; n < iE; n++) {
                const int target_global_id =
                    hlpod_meta->my_global_id[
                        hlpod_meta->item[n]];

                for (int l = 0;
                     l <
                         hlpod_meta->n_internal_sum
                         + monolis_com->n_internal_vertex;
                     l++) {

                    if (target_global_id
                        != hlpod_meta->global_id[l]) {
                        continue;
                    }

                    const int IS =
                        hlpod_ddhr
                            ->num_neib_modes_1stdd_sum[l]
                        + hlpod_vals->n_neib_vec * p;

                    const int IE =
                        hlpod_ddhr
                            ->num_neib_modes_1stdd_sum[l + 1]
                        + hlpod_vals->n_neib_vec * p;

                    for (int j = IS; j < IE; j++) {
                        if (j < 0 || j >= global_NNLS_row) {
                            fprintf(
                                stderr,
                                "ERROR: stage-2 neighbor mode row is out of range: "
                                "rank=%d subdomain=%d snapshot=%d j=%d rows=%d\n",
                                myrank,
                                m,
                                p,
                                j,
                                global_NNLS_row);

                            exit(EXIT_FAILURE);
                        }

                        if (row_seen_in_subdomain[(size_t)j] != 0) {
                            fprintf(
                                stderr,
                                "ERROR: duplicate stage-2 row in one subdomain: "
                                "rank=%d subdomain=%d snapshot=%d j=%d\n",
                                myrank,
                                m,
                                p,
                                j);

                            exit(EXIT_FAILURE);
                        }

                        row_seen_in_subdomain[(size_t)j] = 1;
                        local_equation_count++;

                        for (int i = 0;
                             i < num_local_elems;
                             i++) {

                            const int local_elem_id =
                                hlpod_ddhr
                                    ->elem_id_local[i][m];

                            const int global_elem_id =
                                hlpod_ddhr
                                    ->parallel_elems_id[
                                        local_elem_id];

                            const auto candidate_it =
                                candidate_col.find(
                                    global_elem_id);

                            if (candidate_it
                                == candidate_col.end()) {
                                continue;
                            }

                            const int c =
                                candidate_it->second;

                            const double contribution =
                                hlpod_ddhr
                                    ->matrix[j][i][m];

                            full_matrix[j][c] +=
                                contribution;

                            if (contribution != 0.0) {
                                column_contribution_count[(size_t)c]++;
                            }
                        }

                        full_RH[j] +=
                            hlpod_ddhr->RH[j][m];

                        row_contribution_count[(size_t)j]++;
                    }

                    break;
                }
            }
        }

        if (local_equation_count != NNLS_row) {
            fprintf(
                stderr,
                "ERROR: stage-2 mapped equation count mismatch: "
                "rank=%d subdomain=%d mapped=%d stage1_rows=%d\n",
                myrank,
                m,
                local_equation_count,
                NNLS_row);

            exit(EXIT_FAILURE);
        }
    }

    /* 候補列に少なくとも1つ寄与があることを確認する。 */
    for (int c = 0;
         c < num_candidates;
         c++) {

        if (column_contribution_count[(size_t)c] == 0) {
            fprintf(
                stderr,
                "ERROR: candidate column has no contribution: "
                "rank=%d candidate_col=%d global_elem=%d\n",
                myrank,
                c,
                candidate_elem[(size_t)c]);

            exit(EXIT_FAILURE);
        }
    }

    /* 実際に使用されたグローバル行だけを抽出する。 */
    std::vector<int> used_global_rows;
    used_global_rows.reserve((size_t)global_NNLS_row);

    for (int row = 0;
         row < global_NNLS_row;
         row++) {

        if (row_contribution_count[(size_t)row] > 0) {
            used_global_rows.push_back(row);
        }
    }

    const int stage2_NNLS_row =
        (int)used_global_rows.size();

    if (stage2_NNLS_row <= 0) {
        fprintf(
            stderr,
            "ERROR: stage-2 system has no active mode rows: rank=%d\n",
            myrank);

        exit(EXIT_FAILURE);
    }

    double** global_matrix = NULL;
    double* global_RH = NULL;
    double* global_ans = NULL;

    global_matrix =
        BB_std_calloc_2d_double(
            global_matrix,
            stage2_NNLS_row,
            num_candidates);

    global_RH =
        BB_std_calloc_1d_double(
            global_RH,
            stage2_NNLS_row);

    global_ans =
        BB_std_calloc_1d_double(
            global_ans,
            num_candidates);

    for (int row2 = 0;
         row2 < stage2_NNLS_row;
         row2++) {

        const int global_row =
            used_global_rows[(size_t)row2];

        global_RH[row2] =
            full_RH[global_row];

        for (int c = 0;
             c < num_candidates;
             c++) {

            global_matrix[row2][c] =
                full_matrix[global_row][c];
        }
    }

    /* 元のフル行空間は圧縮後に不要。 */
    BB_std_free_2d_double(
        full_matrix,
        global_NNLS_row,
        num_candidates);

    BB_std_free_1d_double(
        full_RH,
        global_NNLS_row);

    /* 最初の数行だけ、グローバル行との対応を表示する。 */
    const int num_row_debug =
        (stage2_NNLS_row < 10) ? stage2_NNLS_row : 10;

    for (int row2 = 0;
         row2 < num_row_debug;
         row2++) {

        const int global_row =
            used_global_rows[(size_t)row2];

        printf(
            "[STAGE2-ROW] rank=%d compressed_row=%d "
            "global_mode_row=%d contributions=%d RH=%.15e\n",
            myrank,
            row2,
            global_row,
            row_contribution_count[(size_t)global_row],
            global_RH[row2]);
    }

    /*
     * 全領域のactive要素寄与を組み立て、行圧縮した後に、
     * 第2段階全体で1つの共通スケールを用いて正規化する。
     */
    double global_scale_sq = 0.0;

    for (int row = 0;
         row < stage2_NNLS_row;
         row++) {

        for (int c = 0;
             c < num_candidates;
             c++) {

            const double value =
                global_matrix[row][c];

            global_scale_sq +=
                value * value;
        }
    }

    double global_scale =
        std::sqrt(global_scale_sq);

    if (global_scale <= DBL_MIN) {
        global_scale = 1.0;
    }

    for (int row = 0;
         row < stage2_NNLS_row;
         row++) {

        global_RH[row] /=
            global_scale;

        for (int c = 0;
             c < num_candidates;
             c++) {

            global_matrix[row][c] /=
                global_scale;
        }
    }

    printf(
        "\n"
        "============================================================\n"
        "[Stage 2 reduced NNLS system]\n"
        "rank                              = %d\n"
        "number of subdomains              = %d\n"
        "full global mode rows             = %d\n"
        "used stage-2 mode rows            = %d\n"
        "stage-1 active union columns      = %d\n"
        "initial active columns            = %d\n"
        "global normalization scale        = %.15e\n"
        "============================================================\n",
        myrank,
        num_subdomains,
        global_NNLS_row,
        stage2_NNLS_row,
        num_candidates,
        num_initial_active,
        global_scale);

    /*
     * ============================================================
     * 初期 active set 付き全体 sparse NNLS。
     * ============================================================
     */

    double global_residual = 0.0;

    monolis_optimize_nnls_R_with_sparse_solution_initial_set(
        global_matrix,
        global_RH,
        global_ans,
        stage2_NNLS_row,
        num_candidates,
        nnls_max_iter,
        nnls_tol,
        active_set_init.data(),
        &global_residual);

    int num_global_selected = 0;

    std::vector<int> stage2_selected_parallel_elem;
    std::vector<double> stage2_selected_weight;

    stage2_selected_parallel_elem.reserve(
        (size_t)num_candidates);

    stage2_selected_weight.reserve(
        (size_t)num_candidates);

    for (int c = 0;
         c < num_candidates;
         c++) {

        if (global_ans[c] > weight_tol) {
            num_global_selected++;

            stage2_selected_parallel_elem.push_back(
                candidate_elem[(size_t)c]);

            stage2_selected_weight.push_back(
                global_ans[c]);
        }
    }

    /*
     * Save the raw rank contribution for the later serial Stage-3 run.
     * No MPI_* routine is called: the serial executable performs the
     * union of all Stage-2 selected original-global element IDs.
     */
    HROM_stage3_write_rank_file_no_mpi(
        monolis_com,
        hlpod_vals,
        hlpod_ddhr,
        hlpod_meta,
        total_num_elem,
        total_num_snapshot,
        num_subdomains,
        subdomain_id,
        stage2_selected_parallel_elem,
        stage2_selected_weight,
        directory);

    const int num_removed_stage2 =
        num_candidates - num_global_selected;

    const int num_removed_total =
        num_unique_elements - num_global_selected;

    const double stage2_keep_ratio =
        (num_candidates > 0)
        ? 100.0
            * (double)num_global_selected
            / (double)num_candidates
        : 0.0;

    const double stage2_reduction_ratio =
        (num_candidates > 0)
        ? 100.0
            * (double)num_removed_stage2
            / (double)num_candidates
        : 0.0;

    const double total_keep_ratio =
        (num_unique_elements > 0)
        ? 100.0
            * (double)num_global_selected
            / (double)num_unique_elements
        : 0.0;

    const double total_reduction_ratio =
        (num_unique_elements > 0)
        ? 100.0
            * (double)num_removed_total
            / (double)num_unique_elements
        : 0.0;

    printf(
        "\n"
        "============================================================\n"
        "[Two-stage NNLS reduction summary]\n"
        "rank                              = %d\n"
        "number of subdomains              = %d\n"
        "domains with candidates           = %d\n"
        "domains with initial active       = %d\n"
        "active-domain coverage            = %.3f %%\n"
        "element copies over subdomains    = %zu\n"
        "unique physical elements          = %d\n"
        "sum of local candidate counts     = %zu\n"
        "unique stage-1 candidates         = %d\n"
        "candidate duplicates              = %zu\n"
        "initial-active selections         = %d\n"
        "unique initial-active elements    = %d\n"
        "duplicate initial-active picks    = %d\n"
        "stage-2 selected elements         = %d\n"
        "stage-1 removed                   = %d\n"
        "stage-1 retained                  = %.3f %%\n"
        "stage-1 reduction                 = %.3f %%\n"
        "stage-2 removed                   = %d\n"
        "stage-2 retained                  = %.3f %%\n"
        "stage-2 reduction                 = %.3f %%\n"
        "total removed                     = %d\n"
        "total retained                    = %.3f %%\n"
        "total reduction                   = %.3f %%\n"
        "global residual                   = %.15e\n"
        "============================================================\n",
        myrank,
        num_subdomains,
        num_domains_with_candidates,
        num_domains_with_initial_active,
        active_domain_coverage,
        max_local_elements,
        num_unique_elements,
        sum_local_candidate_elements,
        num_candidates,
        duplicate_candidates_across_domains,
        num_domains_with_initial_active,
        num_initial_active,
        duplicate_initial_active,
        num_global_selected,
        num_removed_stage1,
        stage1_keep_ratio,
        stage1_reduction_ratio,
        num_removed_stage2,
        stage2_keep_ratio,
        stage2_reduction_ratio,
        num_removed_total,
        total_keep_ratio,
        total_reduction_ratio,
        global_residual);

    /*
     * ============================================================
     * 全体NNLS結果を各領域へ戻す。
     * ============================================================
     */

    /*
     * num_candidates 以上の領域を確保しておけば、
     * 最終的な選択数を格納できる。
     *
     * 既存コードとの互換性のため、
     * nnls_max_iter より小さくならないようにする。
     */
    const int selection_capacity =
        (num_candidates > nnls_max_iter)
        ? num_candidates
        : nnls_max_iter;

    hlpod_ddhr->D_bc_exists =
        BB_std_calloc_2d_bool(
            hlpod_ddhr->D_bc_exists,
            fe->total_num_nodes,
            num_subdomains);

    hlpod_ddhr->id_selected_elems =
        BB_std_calloc_2d_int(
            hlpod_ddhr->id_selected_elems,
            selection_capacity,
            num_subdomains);

    hlpod_ddhr->id_selected_elems_D_bc =
        BB_std_calloc_2d_int(
            hlpod_ddhr->id_selected_elems_D_bc,
            selection_capacity,
            num_subdomains);

    hlpod_ddhr->elem_weight =
        BB_std_calloc_2d_double(
            hlpod_ddhr->elem_weight,
            selection_capacity,
            num_subdomains);

    hlpod_ddhr->elem_weight_D_bc =
        BB_std_calloc_2d_double(
            hlpod_ddhr->elem_weight_D_bc,
            selection_capacity,
            num_subdomains);

    hlpod_ddhr->num_selected_elems =
        BB_std_calloc_1d_int(
            hlpod_ddhr->num_selected_elems,
            num_subdomains);

    hlpod_ddhr->num_selected_elems_D_bc =
        BB_std_calloc_1d_int(
            hlpod_ddhr->num_selected_elems_D_bc,
            num_subdomains);

    int* total_num_selected_elems =
        (int*)calloc(
            (size_t)num_subdomains,
            sizeof(int));

    if (total_num_selected_elems == NULL) {
        fprintf(
            stderr,
            "ERROR: failed to allocate "
            "total_num_selected_elems\n");

        exit(EXIT_FAILURE);
    }

    for (int m = 0;
         m < num_subdomains;
         m++) {

        int index_D_bc = 0;
        int index_no_D_bc = 0;

        const int num_local_elems =
            hlpod_ddhr->num_elems[m];

        for (int i = 0;
             i < num_local_elems;
             i++) {

            const int local_elem_id =
                hlpod_ddhr
                    ->elem_id_local[i][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            const auto candidate_it =
                candidate_col.find(global_elem_id);

            if (candidate_it == candidate_col.end()) {
                continue;
            }

            const int c =
                candidate_it->second;

            const double weight =
                global_ans[c];

            if (weight <= weight_tol) {
                continue;
            }

            bool has_D_bc = false;

            /*
             * 選択要素にDirichlet境界節点が含まれるか確認する。
             */
            for (int a = 0;
                 a < num_nodes_per_elem;
                 a++) {

                const int node =
                    fe->conn[local_elem_id][a];

                if (node < 0 ||
                    node >= fe->total_num_nodes) {
                    continue;
                }

                for (int k = 0;
                     k < dof;
                     k++) {

                    if (bc->D_bc_exists[
                            node * dof + k]) {

                        has_D_bc = true;

                        hlpod_ddhr
                            ->D_bc_exists[node][m] =
                            true;
                    }
                }
            }

            if (has_D_bc) {
                if (index_D_bc >= selection_capacity) {
                    fprintf(
                        stderr,
                        "ERROR: selected D_bc elements "
                        "exceed capacity\n");

                    exit(EXIT_FAILURE);
                }

                hlpod_ddhr
                    ->id_selected_elems_D_bc[
                        index_D_bc][m] =
                    local_elem_id;

                hlpod_ddhr
                    ->elem_weight_D_bc[
                        index_D_bc][m] =
                    weight;

                index_D_bc++;
            }
            else {
                if (index_no_D_bc >= selection_capacity) {
                    fprintf(
                        stderr,
                        "ERROR: selected elements "
                        "exceed capacity\n");

                    exit(EXIT_FAILURE);
                }

                hlpod_ddhr
                    ->id_selected_elems[
                        index_no_D_bc][m] =
                    local_elem_id;

                hlpod_ddhr
                    ->elem_weight[
                        index_no_D_bc][m] =
                    weight;

                index_no_D_bc++;
            }
        }

        hlpod_ddhr->num_selected_elems[m] =
            index_no_D_bc;

        hlpod_ddhr->num_selected_elems_D_bc[m] =
            index_D_bc;

        total_num_selected_elems[m] =
            index_no_D_bc + index_D_bc;

        printf(
            "rank = %d, subdomain = %d, "
            "selected = %d, "
            "D_bc = %d, no_D_bc = %d\n",
            myrank,
            m,
            total_num_selected_elems[m],
            index_D_bc,
            index_no_D_bc);

        /*
         * 全体NNLSの残差を各領域の結果として出力する。
         */
        hr_write_NNLS_residual(
            global_residual,
            myrank,
            m,
            directory);

        hr_write_NNLS_num_elems(
            total_num_selected_elems[m],
            myrank,
            m,
            directory);

        /*
         * 各領域のファイル出力。
         */
        char fname_D_bc[BUFFER_SIZE];
        char fname_no_D_bc[BUFFER_SIZE];

        snprintf(
            fname_D_bc,
            BUFFER_SIZE,
            "DDECM/lb_selected_elem_D_bc.%d.txt",
            subdomain_id[m]);

        snprintf(
            fname_no_D_bc,
            BUFFER_SIZE,
            "DDECM/lb_selected_elem.%d.txt",
            subdomain_id[m]);

        FILE* fp_D_bc =
            ROM_BB_write_fopen(
                fp_D_bc,
                fname_D_bc,
                directory);

        FILE* fp_no_D_bc =
            ROM_BB_write_fopen(
                fp_no_D_bc,
                fname_no_D_bc,
                directory);

        fprintf(
            fp_D_bc,
            "%d\n",
            index_D_bc);

        fprintf(
            fp_no_D_bc,
            "%d\n",
            index_no_D_bc);

        for (int h = 0;
             h < index_D_bc;
             h++) {

            const int local_elem_id =
                hlpod_ddhr
                    ->id_selected_elems_D_bc[h][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            const double weight =
                hlpod_ddhr
                    ->elem_weight_D_bc[h][m];

            fprintf(
                fp_D_bc,
                "%d %.30e\n",
                global_elem_id,
                weight);
        }

        for (int h = 0;
             h < index_no_D_bc;
             h++) {

            const int local_elem_id =
                hlpod_ddhr
                    ->id_selected_elems[h][m];

            const int global_elem_id =
                hlpod_ddhr
                    ->parallel_elems_id[local_elem_id];

            const double weight =
                hlpod_ddhr
                    ->elem_weight[h][m];

            fprintf(
                fp_no_D_bc,
                "%d %.30e\n",
                global_elem_id,
                weight);
        }

        fclose(fp_D_bc);
        fclose(fp_no_D_bc);
    }

    /*
     * ============================================================
     * メモリ解放。
     * ============================================================
     */

    BB_std_free_2d_double(
        global_matrix,
        stage2_NNLS_row,
        num_candidates);

    BB_std_free_1d_double(
        global_RH,
        stage2_NNLS_row);

    BB_std_free_1d_double(
        global_ans,
        num_candidates);

    free(total_num_selected_elems);

    free(subdomain_scale);
    free(subdomain_num_rows);
    free(subdomain_id);

    /*
     * 元コードと同様に、
     * hlpod_ddhr 内の入力用行列を解放する。
     */
    const int max_num_elem =
        ROM_BB_findMax(
            hlpod_ddhr->num_elems,
            num_subdomains);

    BB_std_free_3d_double(
        hlpod_ddhr->matrix,
        total_num_snapshot
            * hlpod_vals->n_neib_vec,
        max_num_elem,
        num_subdomains);

    BB_std_free_2d_double(
        hlpod_ddhr->RH,
        total_num_snapshot
            * hlpod_vals->n_neib_vec,
        num_subdomains);

    const double end_time =
        monolis_get_time_global_sync();

    printf(
        "rank = %d, HROM global candidate NNLS time = %.15e\n",
        myrank,
        end_time - start_time);
}