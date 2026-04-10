//DDHROMに関して、オーバーラップ要素を含んで計算する方式
//内部要素の総和が分割前の要素数の総和になることを利用
//load balancing, 2階層目のbscr形式に対応

#include "DDHR_para_lb.h"
#include "DDHR_para_lb_decoupled.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";

static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";

//内部要素とオーバーラップ要素の出力
void HROM_ddecm_get_selected_elems_int_ovl_decoupled(
	HLPOD_DDHR*     hlpod_ddhr,
    const char*     fphs,
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
		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[n]);
		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, subdomain_id[n]);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(val));
		num_selected_elems_D_bc += val;
		fclose(fp2);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, meta_list_neib[m]);
			fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

			fscanf(fp1, "%d", &(val));
			num_selected_elems += val;
			fclose(fp1);
		}

		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, meta_list_neib[m]);
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
		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[n]);
		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);
		fscanf(fp1, "%d", &(val));
		for (int j = 0; j < val; j++) {
			fscanf(fp1, "%d %lf", &(ovl_selected_elems[j + index]), &(ovl_selected_elems_weight[j + index]));
		}
		index += val;

		fclose(fp1);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, meta_list_neib[m]);
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
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, subdomain_id[n]);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);
		fscanf(fp2, "%d", &(val));
		for (int j = 0; j < val; j++) {
			fscanf(fp2, "%d %lf", &(ovl_selected_elems_D_bc[j + index]), &(ovl_selected_elems_weight_D_bc[j + index]));
		}
		index += val;
		fclose(fp2);

		/*隣接領域*/
		for (int m = 0; m < meta_n_neib; m++) {
			snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, meta_list_neib[m]);
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

		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/selected_elem_overlap.%d.txt", fphs, subdomain_id[n]);
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

		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/num_selected_node.%d.txt", fphs, subdomain_id[n]);
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

		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/selected_elem_internal.%d.txt", fphs, subdomain_id[n]);

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

void HROM_ddecm_read_selected_elems_para_decoupled(
	const int num_subdomains,
    const char*     fphs,
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

	snprintf(fname3, BUFFER_SIZE, "DDECM_%s/selected_elem_D_bc.%d.txt", fphs, monolis_mpi_get_global_my_rank());
	snprintf(fname4, BUFFER_SIZE, "DDECM_%s/selected_elem.%d.txt", fphs, monolis_mpi_get_global_my_rank());

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
		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs,  subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[m]);

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
		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[m]);

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

void HROM_ddecm_write_selected_elems_para_arbit_subd_decoupled(
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
    const char*     fphs,
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
	const double TOL = 1.0e-12;

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

        printf("local norm = %e, local_Fnorm = %e ", local_norm, local_Frovnorm);


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

		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[m]);

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

void HROM_ddecm_set_neib_decoupled(
		MONOLIS_COM*  	monolis_com,
		HLPOD_MAT* 	hlpod_mat,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_subdomains,
		const int       num_snapshots,
        const char*     fphs,
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
	snprintf(filename, BUFFER_SIZE,"DDECM_%s/n_modes_internal.%d.txt", fphs, myrank);
	fp = ROM_BB_write_fopen(fp, filename, directory);

	fprintf(fp, "%d\n", monolis_com->n_internal_vertex);
	for(int j = 0; j < monolis_com->n_internal_vertex; j++){
		fprintf(fp, "%d\n", hlpod_mat->num_modes_internal[j]);
	}
	fclose(fp);
	/**/

	double t = monolis_get_time_global_sync();

	hlpod_ddhr->num_neib_modes_1stdd = BB_std_calloc_1d_int(hlpod_ddhr->num_neib_modes_1stdd, hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex);

	snprintf(filename, BUFFER_SIZE, "DDECM_%s/n_modes_internal.%d.txt", fphs, monolis_mpi_get_global_my_rank());
	fp = ROM_BB_read_fopen(fp, filename, directory);

	int index_internal = 0;
	fscanf(fp, "%d",&(tmp));
	for(int i = 0; i < monolis_com->n_internal_vertex; i++) {
		fscanf(fp, "%d", &(hlpod_ddhr->num_neib_modes_1stdd[i]));
		index_internal++;
	}
	fclose(fp);

	for (int m = 0; m < hlpod_meta->num_neib; m++){
		snprintf(filename, BUFFER_SIZE, "DDECM_%s/n_modes_internal.%d.txt", fphs, hlpod_meta->neib_id[m]);
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

void HROM_ddecm_get_selected_elema_add_decoupled(
	HLPOD_DDHR*     hlpod_ddhr,
	const int       num_parallel_subdomains,
    const char*     fphs,
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
	}

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
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
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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

void HROM_ddecm_get_selected_elema_add_decoupled_v(
	HLPOD_DDHR*     hlpod_ddhr,
	const int       num_parallel_subdomains,
    const char*     fphs,
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
	}

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
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
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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

	hlpod_ddhr->ovl_id_selected_elems_v = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_v, index1);
	hlpod_ddhr->ovl_elem_weight_v = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_v, index1);
	hlpod_ddhr->ovl_id_selected_elems_D_bc_v = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_D_bc_v, index2);
	hlpod_ddhr->ovl_elem_weight_D_bc_v = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_D_bc_v, index2);

	hlpod_ddhr->ovl_num_selected_elems_v = index1;
	hlpod_ddhr->ovl_num_selected_elems_D_bc_v = index2;

	index1 = 0;
	index2 = 0;

	for(int i = 0; i < num_selected_elems; i++){
		if(bool_ovl_selected_elems[i]){

			hlpod_ddhr->ovl_id_selected_elems_v[index1] = ovl_elem_local_id[index1];
			hlpod_ddhr->ovl_elem_weight_v[index1] = ovl_selected_elems_weight[i];
			index1++;
		}
	}

	for(int i = 0; i < num_selected_elems_D_bc; i++){
		if(bool_ovl_selected_elems_D_bc[i]){
			hlpod_ddhr->ovl_id_selected_elems_D_bc_v[index2] = ovl_elem_local_id_D_bc[index2];
			hlpod_ddhr->ovl_elem_weight_D_bc_v[index2] = ovl_selected_elems_weight_D_bc[i];
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


void HROM_ddecm_get_selected_elema_add_decoupled_p(
	HLPOD_DDHR*     hlpod_ddhr,
	const int       num_parallel_subdomains,
    const char*     fphs,
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
	}

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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
		snprintf(fname1, BUFFER_SIZE,"DDECM_%s/selected_elem.%d.txt", fphs, m);
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
		snprintf(fname2, BUFFER_SIZE,"DDECM_%s/selected_elem_D_bc.%d.txt", fphs, m);
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

	hlpod_ddhr->ovl_id_selected_elems_p = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_p, index1);
	hlpod_ddhr->ovl_elem_weight_p = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_p, index1);
	hlpod_ddhr->ovl_id_selected_elems_D_bc_p = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_D_bc_p, index2);
	hlpod_ddhr->ovl_elem_weight_D_bc_p = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_D_bc_p, index2);

	hlpod_ddhr->ovl_num_selected_elems_p = index1;
	hlpod_ddhr->ovl_num_selected_elems_D_bc_p = index2;

	index1 = 0;
	index2 = 0;

	for(int i = 0; i < num_selected_elems; i++){
		if(bool_ovl_selected_elems[i]){

			hlpod_ddhr->ovl_id_selected_elems_p[index1] = ovl_elem_local_id[index1];
			hlpod_ddhr->ovl_elem_weight_p[index1] = ovl_selected_elems_weight[i];
			index1++;
		}
	}

	for(int i = 0; i < num_selected_elems_D_bc; i++){
		if(bool_ovl_selected_elems_D_bc[i]){
			hlpod_ddhr->ovl_id_selected_elems_D_bc_p[index2] = ovl_elem_local_id_D_bc[index2];
			hlpod_ddhr->ovl_elem_weight_D_bc_p[index2] = ovl_selected_elems_weight_D_bc[i];
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


/*for visualization*/
void ddhr_lb_set_selected_elems_para_decoupled(
		BBFE_DATA*     	fe,
		HLPOD_DDHR*     hlpod_ddhr,
		const int		total_num_nodes,
		const int		num_subdomains,
        const char*     fphs,
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

	snprintf(fname, BUFFER_SIZE,"DDECM_%s/%s", fphs, OUTPUT_FILENAME_ECM_ELEM_VTK);
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


void HROM_ddecm_write_selected_elems_para_arbit_subd_svd_decoupled(
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
    const char*     fphs,
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
            10e-10);

        //k = 100;

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

		snprintf(fname1, BUFFER_SIZE, "DDECM_%s/lb_selected_elem_D_bc.%d.txt", fphs, subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM_%s/lb_selected_elem.%d.txt", fphs, subdomain_id[m]);

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
