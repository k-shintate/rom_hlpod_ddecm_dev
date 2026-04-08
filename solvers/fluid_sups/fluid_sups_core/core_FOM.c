
#include "core_FOM.h"
#include "core_NR.h"
#include <mkl.h>


const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 3;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-7;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 30000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;
const char*         ID_DENSITY  = "#density";
const double      DVAL_DENSITY  = 1.0;
const char*       ID_VISCOSITY  = "#viscosity";
const double    DVAL_VISCOSITY  = 0.01;

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_CAVITY = "cavity_Re%e.txt";

static const char* OUTPUT_FILENAME_KARMAN_VORTEX_P = "karman_vortex_p_val.dat";
static const char* OUTPUT_FILENAME_KARMAN_VORTEX_N = "karman_vortex_n_val.dat";

static const char* OUTPUT_FILENAME_KARMAN_VORTEX_CP = "karman_vortex_Cp_%d.dat";
static const char* OUTPUT_FILENAME_KARMAN_VORTEX_PINF = "karman_vortex_pinf_%d.dat";

static const char* INPUT_FILENAME_SURF = "surf_graph.dat";

double epsilon = 1.0e-5;

// メッシュは各方向41節点40要素 or 51節点50要素 or 101節点100要素 の六面体一次要素限定
void output_cavity_center_vx(
		VALUES*        vals,
		const char*    method,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_CAVITY, method, reynolds);

	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	for(int i=0; i<51; i++){
	 	double z = 0.02*i;
	 	int n = i*51*51 + 25*51 + 25;
	 	double vx = vals->v[n][0];

		fprintf(fp, "%lf %lf\n", z, vx);
	}

	/*
	for(int i=0; i<101; i++){
		double z = 0.01*i;
		int n = i*101*101 + 50*101 + 50;
		double vx = vals->v[n][0];

		fprintf(fp, "%lf %lf\n", z, vx);
	}
	*/

	fclose(fp);
}


void initialize_velocity_pressure_karman_vortex(
	double** v,
	double* p,
	const int total_num_nodes)
{
    // Initialize velocity array
    for (int i = 0; i < total_num_nodes; i++) {
        v[i][0] = 1.0;
    }

    // Initialize pressure array
    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void initialize_velocity_pressure_cavity(
	double** v,
    double** v_old,
	double* p,
	const int total_num_nodes)
{
    // Initialize velocity array
    for (int i = 0; i < total_num_nodes; i++) {
        v[i][0] = 1.0;
        v_old[i][0] = 1.0;
    }

    // Initialize pressure array
    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{
	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	srand((unsigned)time(NULL));

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}

	int tmp;
	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(tmp));
	printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		

        if (val == 1.0){
            int r = rand() % 9 + 1;  // 1〜9 の整数を生成
            bc->imposed_D_val[index] = val * (1 + (double)r * epsilon);
        }
    	else{
            bc->imposed_D_val[ index ] = val;
        }

	}

	fclose(fp);
}


void output_result_file_karman_vortex(
        BBFE_DATA*     fe,
		VALUES*        vals,
        double         t,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_P);

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.02;

    double lowerx = 1.25 - val;
    double upperx = 1.25 + val;

    double lowery = 0.5 - val;
    double uppery = 0.5 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        if (fe->x[i][0] > lowerx && fe->x[i][0] < upperx 
            && fe->x[i][1] > lowery && fe->x[i][1] < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->v[i][0],
                    vals->v[i][1],
                    vals->v[i][2]);
        }
    }

	fclose(fp);

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_N);
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    lowery = -0.5 - val;
    uppery = -0.5 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        if (fe->x[i][0] > lowerx && fe->x[i][0] < upperx 
            && fe->x[i][1] > lowery && fe->x[i][1] < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->v[i][0],
                    vals->v[i][1],
                    vals->v[i][2]);
        }
    }

fclose(fp);


}

void output_result_file_karman_vortex_pressure(
        BBFE_DATA*     fe,
		VALUES*        vals,
        double         t,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_CP, monolis_mpi_get_global_my_rank());

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.000001;

    double lowerr = 0.5 - val;
    double upperr = 0.5 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        double X = fe->x[i][0];
        double Y = fe->x[i][1];
        double r = sqrt(X*X + Y*Y);

        if (r > lowerr && r < upperr) {
            fprintf(fp, "%lf %d %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->p[i]);
        }
    }

	fclose(fp);
}


void output_result_file_karman_vortex_pressure_inf(
        BBFE_DATA*     fe,
                VALUES*        vals,
        double         t,
                const char*    directory)
{
        double reynolds = vals->density / vals->viscosity;
        char filename[BUFFER_SIZE];

        snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_PINF, monolis_mpi_get_global_my_rank());

        FILE* fp;
        fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.1;

    double lowerx = 21.2 - val;
    double upperx = 21.2 + val;

    double lowery = 8.4 - val;
    double uppery = 8.4 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        double X = fe->x[i][0];
        double Y = fe->x[i][1];

        if (X > lowerx && X < upperx &&
			Y > lowery && Y < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->p[i]);
        }
    }

        fclose(fp);
}

void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
    vals->v_old = BB_std_calloc_2d_double(vals->v_old, total_num_nodes, 3);
	vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);

    vals->delta_v = BB_std_calloc_2d_double(vals->delta_v, total_num_nodes, 3);
	vals->delta_p = BB_std_calloc_1d_double(vals->delta_p, total_num_nodes);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;

	vals->density          = DVAL_DENSITY;
	vals->viscosity        = DVAL_VISCOSITY;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s %s: %e\n", CODENAME, ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ID_VISCOSITY,        vals->viscosity);
	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->density), filename, ID_DENSITY, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->viscosity), filename, ID_VISCOSITY, BUFFER_SIZE, CODENAME);


		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t)
{
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
	BB_vtk_write_point_vals_vector(fp, vals->v, fe->total_num_nodes, "Velocity");
	BB_vtk_write_point_vals_scalar(fp, vals->p, fe->total_num_nodes, "Pressure");

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		const int file_num,
		double t)
{
	int myrank = monolis_mpi_get_global_my_rank();
	char fname_vtk[BUFFER_SIZE];
    
	const char* filename;
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num, myrank);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	output_result_file_vtk(
			&(sys->fe), &(sys->vals), filename, sys->cond.directory, t);

}

void BBFE_fluid_sups_read_Dirichlet_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{
	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = ROM_BB_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}

	int tmp;
	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(tmp));
	printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}

	fclose(fp);
}


void BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
                BBFE_BC*     bc,
                VALUES*      vals,
                const char*  filename,
                const char*  directory,
                const int    total_num_nodes,
                const int    block_size)
{
        srand((unsigned)time(NULL));

        bc->total_num_nodes = total_num_nodes;
        bc->block_size      = block_size;

        BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
        int n = total_num_nodes * bc->block_size;

        for(int i=0; i<n; i++) {
                bc->D_bc_exists[i]   = false;
                bc->imposed_D_val[i] = 0.0;
        }

        FILE* fp;
        fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
        if( fp == NULL ) {
                printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
                                CODENAME, filename);
                return;
        }

        int tmp;
        BB_std_scan_line(&fp, BUFFER_SIZE,
                        "%d %d", &(bc->num_D_bcs), &(tmp));
        printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

        for(int i=0; i<(bc->num_D_bcs); i++) {
                int node_id;  int block_id;  double val;
                BB_std_scan_line(&fp, BUFFER_SIZE,
                                "%d %d %lf", &node_id, &block_id, &val);

                int index = (bc->block_size)*node_id + block_id;
                bc->D_bc_exists[ index ]   = true;
                //bc->imposed_D_val[ index ] = val;
                double r = ((double)rand() / RAND_MAX) * 2.0 - 1.0; // -1.0 ~ 1.0 の乱数
                bc->imposed_D_val[index] = val * (1.0 + r * epsilon);

                if(block_id < 3){
                        vals->v[node_id][block_id] = val;
                }
                else if(block_id == 3){
                        vals->p[node_id] = val;
                }
                else{
                        printf("%s ERROR: block_id %d is not supported.\n", CODENAME, block_id);
                        exit(EXIT_FAILURE);
                }
        }

        fclose(fp);
}


void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
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

	for(int e=0; e<(fe->total_num_elems); e++) {
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

						monolis_add_scalar_to_sparse_matrix_R(
								monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
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


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
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

	for(int e=0; e<(fe->total_num_elems); e++) {
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

				monolis->mat.R.B[ 4*fe->conn[e][i] + d ] += integ_val[d];
			}
		}
	}

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step)
{
		monolis_clear_mat_value_R(&(sys.monolis));

		set_element_mat(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				4,
				&(sys.bc),
				sys.monolis.mat.R.B);

		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB_N128,
				MONOLIS_PREC_DIAG,
				sys.fe.total_num_nodes,
				sys.vals.mat_epsilon);

		BBFE_fluid_sups_renew_velocity(
				sys.vals.v, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);
		
		BBFE_fluid_sups_renew_pressure(
				sys.vals.p, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);
}



//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    const int step)
{
		monolis_clear_mat_value_R(&(sys.monolis));

		printf("%s --- prediction step ---\n", CODENAME);

		set_element_mat(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				4,
				&(sys.bc),
				sys.monolis.mat.R.B);

		monolis_show_timelog (&(sys.monolis), true);
		monolis_show_iterlog (&(sys.monolis), true);
		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB_N128,
				MONOLIS_PREC_DIAG,
				sys.fe.total_num_nodes,
				sys.vals.mat_epsilon);

		BBFE_fluid_sups_renew_velocity(
				sys.vals.v, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);

		BBFE_fluid_sups_renew_pressure(
				sys.vals.p, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);

		double* vec;
		vec = BB_std_calloc_1d_double(vec, 3*sys.fe.total_num_nodes);

		ROM_BB_vec_copy_2d_to_1d(
			sys.vals.v,
			vec,
			sys.fe.total_num_nodes);

		if(step%sys.vals.snapshot_interval == 0) {
			printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));

			if(monolis_mpi_get_global_comm_size() == 1){
				ROM_std_hlpod_set_snapmat_nobc(
						sys.vals.p,
						&(sys.rom_p.hlpod_mat),
						sys.fe.total_num_nodes,
                        1,
						(int)(step/sys.vals.snapshot_interval));
			}
			else{
				ROM_std_hlpod_set_snapmat_nobc(
						sys.vals.p,
						&(sys.rom_p.hlpod_mat),
						sys.mono_com.n_internal_vertex,
                        1,
						(int)(step/sys.vals.snapshot_interval));
			}

		}

		if(step%sys.vals.snapshot_interval == 0) {
			printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));

			if(monolis_mpi_get_global_comm_size() == 1){
				ROM_sys_hlpod_fe_set_snap_mat_para(
						vec,
						&(sys.rom_v.hlpod_mat),
						&(sys.bc),
						&(sys.rom_sups.rom_bc),	//要変更
						sys.fe.total_num_nodes,
						3,
						((int)step/sys.vals.snapshot_interval));
			}
			else{
				ROM_sys_hlpod_fe_set_snap_mat_para(
						vec,
						&(sys.rom_v.hlpod_mat),
						&(sys.bc),
						&(sys.rom_sups.rom_bc),
                        sys.mono_com.n_internal_vertex,
						3,
						((int)step/sys.vals.snapshot_interval));
			}

		}

		BB_std_free_1d_double(vec, 3*sys.fe.total_num_nodes);

}

/* supg + pspg + NR method */

void update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] += delta_v[i][j];
        }
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] += delta_p[i];
    }
}

void ROM_BB_vec_copy_2d(
        double**  in,   //input
        double**  out,  //output
        const int num1,
        const int num2)
{
    for (int i = 0; i < num1; i++) {
        for(int j = 0; j < num2; j++){
            out[i][j] = in[i][j];
        }
    }
}

double calc_internal_norm_2d(
        double**  in,       //input
        const int num1,
        const int num2)
{
    double norm = 0.0;

    for (int i = 0; i < num1; i++) {
        for(int j = 0; j < num2; j++){
            norm += in[i][j] * in[i][j];
        }
    }

    return norm;
}

void BBFE_fluid_sups_add_velocity_pressure(
                double**  v,
                double*   p,
                double*   sol_vec,
                const int total_num_nodes)
{
        for(int i=0; i<total_num_nodes; i++) {
                for(int d=0; d<3; d++) {
                        sol_vec[ 4*i + d ] = v[i][d];
                }
                sol_vec[ 4*i + 3 ] = p[i];
        }
}

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
                BBFE_BC*     bc,
                const char*  filename,
                const char*  directory,
                const int    total_num_nodes,
                const int    block_size)
{
        bc->total_num_nodes = total_num_nodes;
        bc->block_size      = block_size;

        BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
        int n = total_num_nodes * bc->block_size;

        for(int i=0; i<n; i++) {
                bc->D_bc_exists[i]   = false;
                bc->imposed_D_val[i] = 0.0;
        }

        FILE* fp;
        fp = ROM_BB_read_fopen_without_error(fp, filename, directory);
        if( fp == NULL ) {
                printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
                                CODENAME, filename);
                return;
        }

        int tmp;
        BB_std_scan_line(&fp, BUFFER_SIZE,
                        "%d %d", &(bc->num_D_bcs), &(tmp));
        printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

        for(int i=0; i<(bc->num_D_bcs); i++) {
                int node_id;  int block_id;  double val;
                BB_std_scan_line(&fp, BUFFER_SIZE,
                                "%d %d %lf", &node_id, &block_id, &val);

                int index = (bc->block_size)*node_id + block_id;
                bc->D_bc_exists[ index ]   = true;
                bc->imposed_D_val[ index ] = 0.0;
        }

        fclose(fp);
}

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt)
{
    double G[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double gij = 0.0;
            for (int k = 0; k < 3; ++k) {
                gij += J_inv[k][i] * J_inv[k][j];
            }
            G[i][j] = gij;
        }
    }

    double vGv = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vGv += v[i] * G[i][j] * v[j];
        }
    }

    double trG = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            trG += G[i][j] * G[i][j];
        }
    }

    const double nu              = viscosity / density;
    const double dt_term         = 4.0 / (dt * dt);
    const double advection_term  = vGv;
    const double diffusion_term  = 36.0 * nu * nu * trG;

    double denom = dt_term + advection_term + diffusion_term;

    if (fabs(Jacobian) < 1e-12) {
        denom = 1e-12;
    }
    return 1.0 / sqrt(denom);
}

double BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt)
{
    double G[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double gij = 0.0;
            for (int k = 0; k < 3; ++k) {
                gij += J_inv[k][i] * J_inv[k][j];
            }
            G[i][j] = gij;
        }
    }

    double vGv = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vGv += v[i] * G[i][j] * v[j];
        }
    }

    double trG2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            trG2 += G[i][j] * G[i][j];
        }
    }

    const double nu              = viscosity / density;
    const double dt_term         = 4.0 / (dt * dt);
    const double advection_term  = vGv;
    const double diffusion_term  = 36.0 * nu * nu * trG2;

    double denom = dt_term + advection_term + diffusion_term;

    double trG = 0.0;
    for (int i = 0; i < 3; ++i) {
        trG += G[i][i];
    }

    if (fabs(Jacobian) < 1e-12) {
        denom = 1e-12;
    }

    return sqrt(denom) / trG;
}

static double BBFE_elem_tau_LSIC_matrix_norm(
    const BBFE_DATA*   fe,
    const BBFE_BASIS*  basis,
    int e,                        /* element id */
    int nl,                       /* local_num_nodes */
    int np,                       /* num_integ_points */
    const double* const* v_ip,    /* [np][3]  各積分点の速度ベクトル */
    const double* Jacobian_ip,    /* [np]     |J| at qp */
    double density)
{
    const int nsd = 3;
    const int ndof = nl * nsd;
    const double EPS = 1e-30;

    if (density <= 0.0) {
        /* ρ<=0 の場合の発散回避。方針に応じて 0.0 を返すなどに調整可 */
        density = EPS;
    }

    /* S_ab = ∫ ρ N_a (u·∇N_b) dΩ  （nl×nl）*/
    double* S = (double*)calloc((size_t)nl * nl, sizeof(double));
    /* E_(aα,bβ) = ∫ ρ (∂αN_a)(∂βN_b) dΩ （ndof×ndof）*/
    double* E = (double*)calloc((size_t)ndof * ndof, sizeof(double));
    double* u_dot_gradN = (double*)malloc((size_t)nl * sizeof(double));

    if (!S || !E || !u_dot_gradN) {
        free(S); free(E); free(u_dot_gradN);
        return 0.0;
    }

    /* ループ用一時変数（C89 対応の場合は関数冒頭に宣言を移してください） */
    int p, a, b, alpha, beta;

    for (p = 0; p < np; ++p) {
        const double w = basis->integ_weight[p] * Jacobian_ip[p];
        const double* up = v_ip[p];                       /* u at qp (長さ3) */
        const double* const* dNa = fe->geo[e][p].grad_N;  /* [nl][3] */

        /* u·∇N_b を前計算 */
        for (b = 0; b < nl; ++b) {
            u_dot_gradN[b] = up[0]*dNa[b][0] + up[1]*dNa[b][1] + up[2]*dNa[b][2];
        }

        /* S_ab の積分加算 */
        for (a = 0; a < nl; ++a) {
            const double Na = basis->N[p][a];
            for (b = 0; b < nl; ++b) {
                S[a*nl + b] += density * Na * u_dot_gradN[b] * w;
            }
        }

        /* E_(aα,bβ) += ρ (∂αN_a)(∂βN_b) w */
        for (a = 0; a < nl; ++a) {
            for (b = 0; b < nl; ++b) {
                for (alpha = 0; alpha < nsd; ++alpha) {
                    const int ii_base = a*nsd + alpha;
                    for (beta = 0; beta < nsd; ++beta) {
                        const int jj = b*nsd + beta;
                        const int ii = ii_base;
                        E[ii*ndof + jj] += density * dNa[a][alpha] * dNa[b][beta] * w;
                    }
                }
            }
        }
    }

    /* Frobenius ノルム */
    double normC2 = 0.0;
    for (a = 0; a < nl; ++a) {
        for (b = 0; b < nl; ++b) {
            /* δ_{αβ} により成分方向に同じスカラが入る → nsd 倍 */
            const double Sab = S[a*nl + b];
            normC2 += (double)nsd * Sab * Sab;
        }
    }

    double normE2 = 0.0;
    {
        const int NN = ndof * ndof;
        int k;
        for (k = 0; k < NN; ++k) normE2 += E[k] * E[k];
    }

    free(u_dot_gradN);
    free(S);
    free(E);

    {
        const double nC = sqrt(normC2);
        const double nE = sqrt(normE2);
        return nC / ((nE > EPS) ? nE : EPS);
    }
}


void BBFE_elemmat_fluid_sups_mat_NR(
    double         mat[4][4],
    const double   J_inv[3][3],
    const double   N_i,
    const double   N_j,
    const double   grad_N_i[3],
    const double   grad_N_j[3],
    const double   v[3],
    double**       grad_u,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    const double vdotGradNi = BB_calc_vec3d_dot(v, grad_N_i);
    const double vdotGradNj = BB_calc_vec3d_dot(v, grad_N_j);
    const double div_v = grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    double adv[3] = {0,0,0};
    for (int d = 0; d < 3; ++d) {
        adv[d] = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
    }

    //時間項（Galerkin + SUPG）
    const double M   = density * N_i * N_j;
    const double M_s = density * tau * vdotGradNi * N_j;

    //移流（∇δu 側）
    const double A   = dt * density * N_i * vdotGradNj;
    const double A_s = dt * density * tau * vdotGradNi * vdotGradNj;

    //粘性
    const double D_11 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[0]*grad_N_j[0] );
    const double D_12 = dt * viscosity * grad_N_i[0] * grad_N_j[1];
    const double D_13 = dt * viscosity * grad_N_i[0] * grad_N_j[2];
    const double D_21 = dt * viscosity * grad_N_i[1] * grad_N_j[0];
    const double D_22 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[1]*grad_N_j[1] );
    const double D_23 = dt * viscosity * grad_N_i[1] * grad_N_j[2];
    const double D_31 = dt * viscosity * grad_N_i[2] * grad_N_j[0];
    const double D_32 = dt * viscosity * grad_N_i[2] * grad_N_j[1];
    const double D_33 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[2]*grad_N_j[2] );

    //LSIC (grad-div penalization) 
    const double Kls = dt * density * tau_c;
    const double LS_00 = Kls * grad_N_i[0]*grad_N_j[0];
    const double LS_01 = Kls * grad_N_i[0]*grad_N_j[1];
    const double LS_02 = Kls * grad_N_i[0]*grad_N_j[2];
    const double LS_10 = Kls * grad_N_i[1]*grad_N_j[0];
    const double LS_11 = Kls * grad_N_i[1]*grad_N_j[1];
    const double LS_12 = Kls * grad_N_i[1]*grad_N_j[2];
    const double LS_20 = Kls * grad_N_i[2]*grad_N_j[0];
    const double LS_21 = Kls * grad_N_i[2]*grad_N_j[1];
    const double LS_22 = Kls * grad_N_i[2]*grad_N_j[2];

    //圧力（運動量式の ∂R/∂p）---
    // Galerkin
    const double G1   = - dt * grad_N_i[0] * N_j;
    const double G2   = - dt * grad_N_i[1] * N_j;
    const double G3   = - dt * grad_N_i[2] * N_j;

    //SUPG
    const double G_s1 = + dt * tau * vdotGradNi * grad_N_j[0];
    const double G_s2 = + dt * tau * vdotGradNi * grad_N_j[1];
    const double G_s3 = + dt * tau * vdotGradNi * grad_N_j[2];

    //連続式の ∂R/∂u（Galerkin 部）
    const double C1 = dt * N_i * grad_N_j[0];
    const double C2 = dt * N_i * grad_N_j[1];
    const double C3 = dt * N_i * grad_N_j[2];

    //連続式の ∂R/∂p（PSPG 圧力）
    const double G_p = dt * tau * (grad_N_i[0]*grad_N_j[0] + grad_N_i[1]*grad_N_j[1] + grad_N_i[2]*grad_N_j[2]) / density;

    mat[0][0] = M + M_s + A + A_s + D_11 + LS_00;
    mat[0][1] = D_12 + LS_01;
    mat[0][2] = D_13 + LS_02;
    mat[0][3] = G1 + G_s1;

    mat[1][0] = D_21 + LS_10;
    mat[1][1] = M + M_s + A + A_s + D_22 + LS_11;
    mat[1][2] = D_23 + LS_12;
    mat[1][3] = G2 + G_s2;

    mat[2][0] = D_31 + LS_20;
    mat[2][1] = D_32 + LS_21;
    mat[2][2] = M + M_s + A + A_s + D_33 + LS_22;
    mat[2][3] = G3 + G_s3;

    mat[3][0] = C1;
    mat[3][1] = C2;
    mat[3][2] = C3;
    mat[3][3] = G_p;

    //移流の δv 由来（Galerkin + SUPG）
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            // Galerkin
            mat[i][k] += dt * density * N_i * N_j * grad_u[i][k];
            // SUPG
            mat[i][k] += dt * density * tau * vdotGradNi * N_j * grad_u[i][k];
        }
    }

    //SUPG： (δv·∇N_i) による運動量寄与
    for (int d = 0; d < 3; ++d) {
        const double adv_d = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += dt * density * tau * grad_N_i[k] * N_j * adv_d;
        }
    }

    //圧力(SUPG) の δv 由来
    for (int d = 0; d < 3; ++d) {
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += dt * tau * grad_N_i[k] * N_j * grad_p[d];
        }
    }

    //連続式 PSPG-移流
    for (int k = 0; k < 3; ++k) {
        // ∇δu 側
        mat[3][k] += dt * tau * grad_N_i[k] * vdotGradNj;
        // δv 側
        const double dot_gradNi_graduk =
            grad_N_i[0]*grad_u[0][k] + grad_N_i[1]*grad_u[1][k] + grad_N_i[2]*grad_u[2][k];
        mat[3][k] += dt * tau * N_j * dot_gradNi_graduk;
    }

    //時間(SUPG) の δv 側
    for (int d = 0; d < 3; ++d) {
        const double du_d = du_time[d];
        for (int k = 0; k < 3; ++k) {
            mat[d][k] += density * tau * grad_N_i[k] * N_j * du_d;
        }
    }

    //連続式 PSPG-時間
    for (int k = 0; k < 3; ++k) {
        mat[3][k] += tau  * grad_N_i[k] * N_j;
    }
}


void BBFE_elemmat_fluid_sups_vec_NR(
    double         vec[4],
    const double   N_i,
    const double   grad_N_i[3],
    const double   v[3],
    const double   u_old[3],
    double**       grad_u,
    const double   p_cur,
    const double   grad_p[3],
    const double   density,
    const double   viscosity,
    const double   tau,
    const double   tau_c,
    const double   dt,
    const double   du_time[3])
{
    for (int d = 0; d < 4; ++d) vec[d] = 0.0;
    const double vdotGradNi = BB_calc_vec3d_dot(v, grad_N_i);
    const double div_v = grad_u[0][0] + grad_u[1][1] + grad_u[2][2];

    //移流項
    double adv[3] = {0,0,0};
    for (int d = 0; d < 3; ++d) {
        adv[d] = v[0]*grad_u[d][0] + v[1]*grad_u[d][1] + v[2]*grad_u[d][2];
    }

    // 運動量方程式
    for (int d = 0; d < 3; ++d) {
        //Galerkin：移流
        vec[d] += dt * density * N_i * adv[d];

        //SUPG（移流）
        vec[d] += dt * density * tau * vdotGradNi * adv[d];

        //粘性
        const double gradv_dot_gradNi =
            grad_u[d][0]*grad_N_i[0] + grad_u[d][1]*grad_N_i[1] + grad_u[d][2]*grad_N_i[2];
        vec[d] += dt * viscosity * gradv_dot_gradNi;
        vec[d] += dt * viscosity * div_v * grad_N_i[d];

        //LSIC
        vec[d] += dt * density * tau_c * div_v * grad_N_i[d];

        //圧力：Galerkin（部分積分後）
        vec[d] -= dt * grad_N_i[d] * p_cur;

        //圧力：SUPG
        vec[d] += dt * tau * vdotGradNi * grad_p[d];
    }

    //連続式
    //Galerkin
    vec[3] += dt * N_i * div_v;

    //PSPG
    const double gradN_dot_gradp = grad_N_i[0]*grad_p[0] + grad_N_i[1]*grad_p[1] + grad_N_i[2]*grad_p[2];
    vec[3] += dt * tau * (gradN_dot_gradp) / density;

    //PSPG
    const double gradN_dot_adv = grad_N_i[0]*adv[0] + grad_N_i[1]*adv[1] + grad_N_i[2]*adv[2];
    vec[3] += dt * tau * gradN_dot_adv;

    //Galerkin時間項
    for (int d = 0; d < 3; ++d) {
        double r = 0.0;
        r += density * N_i * (v[d] - u_old[d]);
        r += density * tau * BB_calc_vec3d_dot(v, grad_N_i) * (v[d] - u_old[d]);
        vec[d] += r;
    }

    //PSPG時間項
    vec[3] += tau  * ( BB_calc_vec3d_dot(grad_N_i, v) - BB_calc_vec3d_dot(grad_N_i, u_old) );
}





void set_element_mat_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip;      double* Jacobian_ip;
    val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
    Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

    double**  val_ip_vec   = BB_std_calloc_2d_double(val_ip_vec , 4 , np);
    double*   integ_val_vec= BB_std_calloc_1d_double(integ_val_vec, 4);

    double** local_v       = BB_std_calloc_2d_double(local_v, nl, 3);
    double** v_ip          = BB_std_calloc_2d_double(v_ip   , np, 3);
    double*** grad_v_ip    = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

    double** local_v_old   = BB_std_calloc_2d_double(local_v_old, nl, 3);
    double** v_ip_old      = BB_std_calloc_2d_double(v_ip_old   , np, 3);

    double*  local_p       = BB_std_calloc_1d_double(local_p, nl);
    double*  p_ip          = BB_std_calloc_1d_double(p_ip  , np);
    double** grad_p_ip     = BB_std_calloc_2d_double(grad_p_ip, np, 3);

    double A[4][4];
    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);

        BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);
        
        for (int p = 0; p < np; ++p) {
            BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p)
                for (int d = 0; d < 4; ++d) val_ip_vec[d][p] = 0.0;

            for (int j = 0; j < nl; ++j) {
                for (int p = 0; p < np; ++p) {
                        double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };

                    BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                    const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                    const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                    BBFE_elemmat_fluid_sups_mat_NR(
                        A, J_inv,
                        basis->N[p][i], basis->N[p][j],
                        fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                        v_ip[p], grad_v_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            val_ip[a][b][p] = A[a][b];
                            A[a][b] = 0.0;
                        }
                    }
                }

                for (int a = 0; a < 4; ++a) {
                    for (int b = 0; b < 4; ++b) {
                        const double integ_val = BBFE_std_integ_calc(
                            np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
                    }
                }
            }
        }
    }

    BB_std_free_3d_double(val_ip , 4 , 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip   , np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old,   np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip   , np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
}


void set_element_vec_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip;      double* Jacobian_ip;
    val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
    Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

    double**  val_ip_vec   = BB_std_calloc_2d_double(val_ip_vec , 4 , np);
    double*   integ_val_vec= BB_std_calloc_1d_double(integ_val_vec, 4);

    double** local_v       = BB_std_calloc_2d_double(local_v, nl, 3);
    double** v_ip          = BB_std_calloc_2d_double(v_ip   , np, 3);
    double*** grad_v_ip    = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

    double** local_v_old   = BB_std_calloc_2d_double(local_v_old, nl, 3);
    double** v_ip_old      = BB_std_calloc_2d_double(v_ip_old   , np, 3);

    double*  local_p       = BB_std_calloc_1d_double(local_p, nl);
    double*  p_ip          = BB_std_calloc_1d_double(p_ip  , np);
    double** grad_p_ip     = BB_std_calloc_2d_double(grad_p_ip, np, 3);

    double vec[4];
    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);

        BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);
        
        for (int p = 0; p < np; ++p) {
            BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);
        }
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
        
        for (int i = 0; i < nl; ++i) {
            for (int p = 0; p < np; ++p)
                for (int d = 0; d < 4; ++d) val_ip_vec[d][p] = 0.0;

            for (int p = 0; p < np; ++p) {
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

                const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density, vals->viscosity,
                    tau, tau_c, vals->dt,
                    du_time);

                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                }
            }

            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }

    BB_std_free_3d_double(val_ip , 4 , 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip   , np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old,   np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip   , np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
}

double BBFE_elemmat_fluid_sups_coef_Tezduyar2003(
		BBFE_DATA*   fe,
		const int    e,
		const int    p,
		const double density,
		const double viscosity,
		const double v_ip[3],
		const double grad_v[3],
		const double dt)
{
	double nu = viscosity/density;

	double tau1_inv = 0.0;
	for(int i =0; i<fe->local_num_nodes; i++){
		double v_dot_grad = BB_calc_vec3d_dot(v_ip, fe->geo[e][p].grad_N[i]);
		tau1_inv += fabs(v_dot_grad);
	}

	double tau2_inv = 2.0 / dt;

	double r[3];
	double len_grad_v = BB_calc_vec3d_length(grad_v);
	if(len_grad_v > 1.0e-12){
		r[0] = grad_v[0] / len_grad_v;
		r[1] = grad_v[1] / len_grad_v;
		r[2] = grad_v[2] / len_grad_v;
	}
	else{
		r[0] = 0.0; r[1] = 0.0; r[2] = 0.0;
	}

	double tmp = 0.0;
	for(int i =0; i<fe->local_num_nodes; i++){
		double r_dot_grad = BB_calc_vec3d_dot(r, fe->geo[e][p].grad_N[i]);
		tmp += fabs(r_dot_grad);
	}
	double tau3_inv = nu * tmp * tmp;

	double denom = tau1_inv * tau1_inv + tau2_inv * tau2_inv + tau3_inv * tau3_inv;

	double val = sqrt(1.0/denom);

	return (val);
}


void set_element_mat_NR_Tezuer(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double*** val_ip;      double* Jacobian_ip;
    val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
    Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

    double**  val_ip_vec   = BB_std_calloc_2d_double(val_ip_vec , 4 , np);
    double*   integ_val_vec= BB_std_calloc_1d_double(integ_val_vec, 4);

    double** local_v       = BB_std_calloc_2d_double(local_v, nl, 3);
    double** v_ip          = BB_std_calloc_2d_double(v_ip   , np, 3);
    double*** grad_v_ip    = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);
    double* local_lv = BB_std_calloc_1d_double(local_lv, nl);

    double** local_v_old   = BB_std_calloc_2d_double(local_v_old, nl, 3);
    double** v_ip_old      = BB_std_calloc_2d_double(v_ip_old   , np, 3);

    double*  local_p       = BB_std_calloc_1d_double(local_p, nl);
    double*  p_ip          = BB_std_calloc_1d_double(p_ip  , np);
    double** grad_p_ip     = BB_std_calloc_2d_double(grad_p_ip, np, 3);

	double** grad_v;
	grad_v = BB_std_calloc_2d_double(grad_v, np, 3);

	double* tau_ip;
	tau_ip = BB_std_calloc_1d_double(tau_ip, np);


    double A[4][4];
    double vec[4];
    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        // 局所場取り出し
        BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(local_v_old, fe, vals->v_old, e, 3);

        BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);

        // 積分点へ写像
        for (int p = 0; p < np; ++p) {
            BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(nl, local_p, basis->N[p]);

// 要素 e の積分点 p における速度勾配の計算
			for(int i=0; i<nl; i++){
				local_lv[i] = BB_calc_vec3d_length(local_v[i]);
			}
			BBFE_std_mapping_scalar_grad(grad_v[p], nl, local_lv, fe->geo[e][p].grad_N);

			tau_ip[p] = BBFE_elemmat_fluid_sups_coef_Tezduyar2003(
					fe, e, p, vals->density, vals->viscosity, v_ip[p], grad_v[p], vals->dt);


        }
                double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
                double h_e = cbrt(vol);
/*
const double tau_e = BBFE_elemmat_fluid_sups_coef_Tezduyar2003(
                        fe, e, p,
                        vals->density, vals->viscosity, v_ip[p],
			grad_v_ip[p],
			vals->dt);
*/

const double tau_c_e = BBFE_elem_tau_LSIC_matrix_norm(
        fe, basis, e, nl, np, (const double* const*)v_ip, Jacobian_ip, vals->density);

        // i-行の組み立て
        for (int i = 0; i < nl; ++i) {

            // RHS 保管配列をクリア（p ごと）
            for (int p = 0; p < np; ++p)
                for (int d = 0; d < 4; ++d) val_ip_vec[d][p] = 0.0;

            // === 行列（i,j ブロック）：従来通り i–j–p の順でOK ===
            for (int j = 0; j < nl; ++j) {
                for (int p = 0; p < np; ++p) {
                        double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };
/*
                        double tau = BBFE_elemmat_fluid_sups_coef(
                                                vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);
*/
                    BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);
/*
const double tau_e = BBFE_elemmat_fluid_sups_coef_Tezduyar2003(
                        fe, e, p,
                        vals->density, vals->viscosity, v_ip[p],
                        grad_v_ip[p],
                        vals->dt);
*/
//		  const double tau = tau_e;
const double tau_c = tau_c_e;

		    /*
                    const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                    const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                        J_inv, fe->geo[e][p].Jacobian,
                        vals->density, vals->viscosity, v_ip[p], vals->dt);

                double rel = check_element_jacobian_at_qp(
                //        N_i,grad_N_i, N_j,grad_N_j,
                //        v,u_old, grad_u, p_cur,grad_p,
                //        rho,mu, tau,tau_c, dt, du_time, eps);
                        J_inv,
                        basis->N[p][i], fe->geo[e][p].grad_N[i],
                        basis->N[p][j], fe->geo[e][p].grad_N[j],
                        v_ip[p], v_ip_old[p], grad_v_ip[p], p_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau, tau_c, vals->dt, du_time, 10e-7);
                if(monolis_mpi_get_global_my_rank()==0){
                        printf("rel = %e\n", rel);
                }
*/
                    BBFE_elemmat_fluid_sups_mat_NR(
                        A, J_inv,
                        basis->N[p][i], basis->N[p][j],
                        fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
                        v_ip[p], grad_v_ip[p], grad_p_ip[p],
                        vals->density, vals->viscosity, tau_ip[p], tau_c, vals->dt, du_time);

                    // 行列の一時保存
                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            val_ip[a][b][p] = A[a][b];
                            A[a][b] = 0.0;
                        }
                    }
                }

                // i–j ブロックの積分・加算
                for (int a = 0; a < 4; ++a) {
                    for (int b = 0; b < 4; ++b) {
                        const double integ_val = BBFE_std_integ_calc(
                            np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
                    }
                }
            }

            // === 残差（RHS）：★ j ループの外で、i–p で一回だけ ===
            for (int p = 0; p < np; ++p) {
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                        };

                BB_calc_mat3d_inverse(fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);
/*
                const double tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);

                double tau = BBFE_elemmat_fluid_sups_coef(
                                                vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

                const double tau_c = BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv, fe->geo[e][p].Jacobian,
                    vals->density, vals->viscosity, v_ip[p], vals->dt);
*/
		/*
const double tau_e = BBFE_elemmat_fluid_sups_coef_Tezduyar2003(
                        fe, e, p,
                        vals->density, vals->viscosity, v_ip[p],
                        grad_v_ip[p],
                        vals->dt);
*/
		//const double tau = tau_e;
		const double tau_c = tau_c_e;

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density, vals->viscosity,
                    tau_ip[p], tau_c, vals->dt,
                    du_time);

                for (int d = 0; d < 4; ++d) {
                    val_ip_vec[d][p] = vec[d];
                }
            }

            // 残差の積分・加算
            for (int d = 0; d < 4; ++d) {
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np, val_ip_vec[d], basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[ 4*fe->conn[e][i] + d ] -= integ_val_vec[d];
            }
        }
    }

    BB_std_free_3d_double(val_ip , 4 , 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);

    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip   , np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old,   np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);

    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip   , np);
    BB_std_free_2d_double(grad_p_ip, np, 3);

    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
}


void solver_fom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step)
{
    printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);

    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 1;

    for(int i = 0; i < max_iter_NR; i++){
        printf("\n%s ----------------- TIme step %d : NR step %d ----------------\n", CODENAME, step, i);

        monolis_clear_mat_value_R(&(sys.monolis));

        set_element_mat_NR_Tezuer(
                &(sys.monolis),
                &(sys.fe),
                &(sys.basis),
                &(sys.vals));
        /*
        set_element_vec_NR(
                &(sys.monolis),
                &(sys.fe),
                &(sys.basis),
                &(sys.vals));
*/
        BBFE_sys_monowrap_set_Dirichlet_bc(
                &(sys.monolis),
                sys.fe.total_num_nodes,
                4,
                &(sys.bc_NR),
                sys.monolis.mat.R.B);

        ROM_monowrap_solve(
                &(sys.monolis),
                &(sys.mono_com),
                sys.monolis.mat.R.X,
                MONOLIS_ITER_BICGSAFE,
                MONOLIS_PREC_DIAG,
                sys.fe.total_num_nodes,
                sys.vals.mat_epsilon);

        BBFE_fluid_sups_renew_velocity(
                sys.vals.delta_v,
                sys.monolis.mat.R.X,
                sys.fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
                sys.vals.delta_p,
                sys.monolis.mat.R.X,
                sys.fe.total_num_nodes);

        update_velocity_pressure_NR(
                sys.vals.v,
                sys.vals.delta_v,
                sys.vals.p,
                sys.vals.delta_p,
                sys.fe.total_num_nodes);

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_2d(
            sys.vals.v,
            sys.mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys.vals.delta_v,
            sys.mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.mono_com.n_internal_vertex; ++ii) {
            double pv  = sys.vals.p[ii];
            double dpv = sys.vals.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys.vals.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys.vals.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.mono_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.mono_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
            i, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
        printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e\n",
            i, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p);

        if (conv_v && conv_p) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys.vals.v[ii][d] - sys.vals.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);

            ROM_BB_vec_copy_2d(
                sys.vals.v,
                sys.vals.v_old,
                sys.fe.total_num_nodes,
                3);
            break;
        }

    }
/*
        double* vec;
    vec = BB_std_calloc_1d_double(vec, 3*sys.fe.total_num_nodes);

    ROM_BB_vec_copy_2d_to_1d(
            sys.vals.v,
            vec,
            sys.fe.total_num_nodes);

    //if(step%sys.vals.snapshot_interval == 0) {
    if(step%2 == 0) {
            if(monolis_mpi_get_global_my_rank()==0){
    printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));
}
            if(monolis_mpi_get_global_comm_size() == 1){
                    ROM_std_hlpod_set_snapmat_nobc(
                                    sys.vals.p,
                                    &(sys.rom_p.hlpod_mat),
                                    sys.fe.total_num_nodes,
            1,
                                    (int)(step/sys.vals.snapshot_interval));
            }
            else{
                    ROM_std_hlpod_set_snapmat_nobc(
                                    sys.vals.p,
                                    &(sys.rom_p.hlpod_mat),
                                    sys.mono_com.n_internal_vertex,
            1,
                                    (int)(step/sys.vals.snapshot_interval));
            }

    }

    if(step%sys.vals.snapshot_interval == 0) {
            if(monolis_mpi_get_global_my_rank()==0){
    printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));
}
            if(monolis_mpi_get_global_comm_size() == 1){
                    //ROM_sys_hlpod_fe_set_snap_mat_para(
                    ROM_std_hlpod_set_snapmat_nobc(
                                    vec,
                                    &(sys.rom_v.hlpod_mat),
                                    //&(sys.bc),
                                    //&(sys.rom_sups.rom_bc), //要変更
                                    sys.fe.total_num_nodes,
                                    3,
                                    ((int)step/sys.vals.snapshot_interval));
            }
            else{
                    //ROM_sys_hlpod_fe_set_snap_mat_para(
                    ROM_std_hlpod_set_snapmat_nobc(
                                    vec,
                                    &(sys.rom_v.hlpod_mat),
                                    //&(sys.bc),
                                    //&(sys.rom_sups.rom_bc),
            sys.mono_com.n_internal_vertex,
                                    3,
                                    ((int)step/sys.vals.snapshot_interval));
            }

    }

    BB_std_free_1d_double(vec, 3*sys.fe.total_num_nodes);
*/
}

void memory_allocation_nodal_values_AB2(
                VALUES*         vals,
                const int       total_num_nodes)
{
        vals->v     = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
        vals->p     = BB_std_calloc_1d_double(vals->p, total_num_nodes);

        vals->delta_v = BB_std_calloc_2d_double(vals->delta_v, total_num_nodes, 3);
        vals->delta_p     = BB_std_calloc_1d_double(vals->delta_p, total_num_nodes);

        vals->v_old = BB_std_calloc_2d_double(vals->v_old, total_num_nodes, 3);
}


void read_connectivity_graph_surf(
    BBFE_DATA *fe,
    const char *filename,
    const char *directory,
    int num_integ_points)
{
    FILE* fp;
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    int graph_ndof;
    int tmp;

    snprintf(fname_n_internal_graph, BUFFER_SIZE, "parted.0/%s.n_internal.%d", filename, monolis_mpi_get_global_my_rank());
    fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
    fscanf(fp, "%s %d", char_n_internal, &(tmp));
    fscanf(fp, "%d", &(graph_ndof));
    fclose(fp);

    snprintf(fname_n_internal_graph, BUFFER_SIZE, "parted.0/%s.%d", filename, monolis_mpi_get_global_my_rank());
    fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
    // read the num of elements
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &(fe->total_num_elems));

    fe->total_num_elems = graph_ndof;
    printf("%s Num. elements: %d\n", CODENAME, fe->total_num_elems);
    fe->local_num_nodes = 4; // 1st-order rectangle
    BBFE_sys_memory_allocation_elem(fe, num_integ_points, 3);
    if (fe->total_num_elems < 0)
    {
        exit(EXIT_FAILURE);
    }
    if (fe->total_num_elems == 0)
    {
        return;
    }
    fe->conn[0][0] = 1;
    int id = 0; // 　elem id (not used)
    for (int e = 0; e < graph_ndof; e++)
    {
        if (fscanf(fp, "%d ", &id) != 1)
        {
            exit(EXIT_FAILURE);
        }
        if (fscanf(fp, "%d ", &(fe->local_num_nodes)) != 1)
        {
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < (fe->local_num_nodes); i++)
        {
            if(fscanf(fp, "%d", &(fe->conn[e][i])) != 1)
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    fclose(fp);
}

void BBFE_convdiff_set_basis_surface(
                BBFE_BASIS*   basis,
                int           local_num_nodes,
                int           num_integ_points_each_axis)
{
        switch( local_num_nodes ) {
                case 3:
                        basis->num_integ_points =
                                BBFE_std_integ_tri_set_arbitrary_points(
                                                num_integ_points_each_axis,
                                                basis->integ_point,
                                                basis->integ_weight);

                        for(int i=0; i<(basis->num_integ_points); i++) {
                                BBFE_std_shapefunc_tri1st_get_val(
                                                basis->integ_point[i],
                                                basis->N[i]);

                                BBFE_std_shapefunc_tri1st_get_derivative(
                                                basis->integ_point[i],
                                                basis->dN_dxi[i],
                                                basis->dN_det[i]);
                        }
                        printf("%s Surface element type: 1st-order triangle.\n", CODENAME);
                        break;

                case 4:
                        basis->num_integ_points =
                                BBFE_std_integ_rec_set_arbitrary_points(
                                                num_integ_points_each_axis,
                                                basis->integ_point,
                                                basis->integ_weight);

                        for(int i=0; i<(basis->num_integ_points); i++) {
                                BBFE_std_shapefunc_rec1st_get_val(
                                                basis->integ_point[i],
                                                basis->N[i]);

                                BBFE_std_shapefunc_rec1st_get_derivative(
                                                basis->integ_point[i],
                                                basis->dN_dxi[i],
                                                basis->dN_det[i]);
                        }
                        printf("%s Surface element type: 1st-order rectangle.\n", CODENAME);
                        break;

                case 6:
                case 9:
                        // should be implemented for higher order elements
                        break;
        }
}

void pre_surface(
                //MONOLIS_COM*  monolis_com_surf,
                BBFE_DATA*    surf,
                BBFE_BASIS*   basis,
                const char*   directory,
		const char*   fname,
                int           num_integ_points_each_axis)
{
        int n_axis = num_integ_points_each_axis;
        const char* filename;

        filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, fname);

        read_connectivity_graph_surf(
                        surf,
                        fname,
                        directory,
                        n_axis * n_axis);

        BBFE_sys_memory_allocation_integ(
                        basis,
                        n_axis*n_axis*n_axis,
                        2);
        BBFE_sys_memory_allocation_shapefunc(
                        basis,
                        surf->local_num_nodes,
                        1,
                        n_axis*n_axis*n_axis);

        BBFE_convdiff_set_basis_surface(
                        basis,
                        surf->local_num_nodes,
                        n_axis);
}

/* 以下、from Chat GPT */

/* 内積・正規化のユーティリティ（環境に同等関数があるならそちらでOK） */
static inline double vdot3(const double a[3], const double b[3]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline double vnorm3(const double a[3]){
    return sqrt(vdot3(a,a));
}
static inline void vscale3(double a[3], double s){
    a[0]*=s; a[1]*=s; a[2]*=s;
}
static inline void vadd3(double r[3], const double a[3]){
    r[0]+=a[0]; r[1]+=a[1]; r[2]+=a[2];
}
static inline void vsub3(double r[3], const double a[3], const double b[3]){
    r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2];
}
static inline void vcopy3(double r[3], const double a[3]){
    r[0]=a[0]; r[1]=a[1]; r[2]=a[2];
}

void calc_Cd_p(
    BBFE_DATA*  surf,
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    VALUES*     vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out)
{
    /* 出力初期化 */
    if (D_out) *D_out = 0.0;
    if (L_out) *L_out = 0.0;

    /* 方向ベクトルを正規化＋直交化（Gram–Schmidt） */
    double eU[3]; vcopy3(eU, eU_in);
    double eP0[3]; vcopy3(eP0, eP_in);
    double nrm = vnorm3(eU); if(nrm>0) vscale3(eU, 1.0/nrm);
    /* eP を eU に直交化して正規化 */
    double eP[3] = { eP0[0] - vdot3(eP0,eU)*eU[0],
                     eP0[1] - vdot3(eP0,eU)*eU[1],
                     eP0[2] - vdot3(eP0,eU)*eU[2] };
    nrm = vnorm3(eP); if(nrm>0) vscale3(eP, 1.0/nrm);

    /* 作業バッファ（面の局所節点値） */
    double** local_x = BB_std_calloc_2d_double(NULL, surf->local_num_nodes, 3);
    double*  local_p = BB_std_calloc_1d_double(NULL, surf->local_num_nodes);

    /* 圧力による合力ベクトルとチェック用量 */
    double Fp[3] = {0,0,0};
    double nsum[3] = {0,0,0};   /* Σ n dS （閉曲面なら 0 が理想）*/
    double area_sum = 0.0;      /* Σ dS (= Σ w*J) */

    /* 形状の“中心”を近似（面の全頂点平均）。円柱の側面なら十分外向き判定に使える */
    double xcenter[3] = {0,0,0};
    {
        double sum[3] = {0,0,0};
        long long cnt = 0;
        for(int e=0; e<surf->total_num_elems; ++e){
            BBFE_elemmat_set_local_array_vector(local_x, surf, fe->x, e, 3);
            for(int a=0; a<surf->local_num_nodes; ++a){
                sum[0] += local_x[a][0];
                sum[1] += local_x[a][1];
                sum[2] += local_x[a][2];
                ++cnt;
            }
        }
        if (cnt > 0){
            xcenter[0] = sum[0]/(double)cnt;
            xcenter[1] = sum[1]/(double)cnt;
            xcenter[2] = sum[2]/(double)cnt;
        }
    }

    /* 本積分ループ */
    for(int e=0; e<surf->total_num_elems; ++e){
        BBFE_elemmat_set_local_array_vector(local_x, surf, fe->x,  e, 3);
        BBFE_elemmat_set_local_array_scalar(local_p, surf, vals->p, e); /* 実圧力 */
        //for(int k=0; k<surf->local_num_nodes; ++k) local_p[k] = 1.0;      /* 一様圧力テスト */

        /* 面重心（外向き判定用） */
        double xf[3]={0,0,0};
        for(int a=0; a<surf->local_num_nodes; ++a){
            xf[0]+=local_x[a][0]; xf[1]+=local_x[a][1]; xf[2]+=local_x[a][2];
        }
        double inv = 1.0 / (double)surf->local_num_nodes;
        xf[0]*=inv; xf[1]*=inv; xf[2]*=inv;

        for(int p=0; p<basis->num_integ_points; ++p){
            /* 外向き法線（長さは任意） */
            double nvec[3];
            BBFE_set_surface_get_outward_normal_vector(
                surf->local_num_nodes, local_x,
                basis->dN_dxi[p], basis->dN_det[p], nvec);

            /* J と単位法線を分離 */
            double J = vnorm3(nvec);
            if (J <= 1e-300) continue;
            double nh[3] = { nvec[0]/J, nvec[1]/J, nvec[2]/J };

            /* 外向き保証：中心から面重心へのベクトルと同じ向きに */
            double svec[3]; vsub3(svec, xf, xcenter);
            if (vdot3(svec, nh) < 0.0){ nh[0]*=-1.0; nh[1]*=-1.0; nh[2]*=-1.0; }

            /* 求積点の圧力（線形/双線形補間） */
            double p_q = BBFE_std_mapping_scalar(
                surf->local_num_nodes, local_p, basis->N[p]);

            /* 面トラクション -p * n を積分（重み = w_q * J） */
            double w = basis->integ_weight[p] * J;
            Fp[0] += (-p_q) * nh[0] * w;
            Fp[1] += (-p_q) * nh[1] * w;
            Fp[2] += (-p_q) * nh[2] * w;
        }
    }

    *D_out = vdot3(Fp, eU);
    *L_out = vdot3(Fp, eP);

    BB_std_free_2d_double(local_x, surf->local_num_nodes, 3);
    BB_std_free_1d_double(local_p, surf->local_num_nodes);
}


/* ---- このファイル内で使うキー型（プロジェクトの node_id_t に依存しない） ---- */
typedef long long dl_node_id_t;

/* --- 参照要素の符号テーブル（固定） --- */
static const int SGN[8][3] = {
  {-1,-1,-1}, {+1,-1,-1}, {+1,+1,-1}, {-1,+1,-1},
  {-1,-1,+1}, {+1,-1,+1}, {+1,+1,+1}, {-1,+1,+1}
};

static inline int find_node_by_signs(const int sgn[8][3], int sx, int sy, int sz){
    for(int a=0;a<8;++a)
        if(sgn[a][0]==sx && sgn[a][1]==sy && sgn[a][2]==sz) return a;
    return -1;
}

/* faces[6][4] を「参照要素で外向き法線が出る CCW 順」で生成
 * 面の並び: z- , z+ , y- , x+ , y+ , x-  */
static inline void build_hex8_faces_from_sgn(const int sgn[8][3], int faces[6][4]){
    faces[0][0] = find_node_by_signs(sgn,-1,-1,-1);
    faces[0][1] = find_node_by_signs(sgn,+1,-1,-1);
    faces[0][2] = find_node_by_signs(sgn,+1,+1,-1);
    faces[0][3] = find_node_by_signs(sgn,-1,+1,-1);

    faces[1][0] = find_node_by_signs(sgn,-1,-1,+1);
    faces[1][1] = find_node_by_signs(sgn,+1,-1,+1);
    faces[1][2] = find_node_by_signs(sgn,+1,+1,+1);
    faces[1][3] = find_node_by_signs(sgn,-1,+1,+1);

    faces[2][0] = find_node_by_signs(sgn,-1,-1,-1);
    faces[2][1] = find_node_by_signs(sgn,+1,-1,-1);
    faces[2][2] = find_node_by_signs(sgn,+1,-1,+1);
    faces[2][3] = find_node_by_signs(sgn,-1,-1,+1);

    faces[3][0] = find_node_by_signs(sgn,+1,-1,-1);
    faces[3][1] = find_node_by_signs(sgn,+1,+1,-1);
    faces[3][2] = find_node_by_signs(sgn,+1,+1,+1);
    faces[3][3] = find_node_by_signs(sgn,+1,-1,+1);

    faces[4][0] = find_node_by_signs(sgn,-1,+1,-1);
    faces[4][1] = find_node_by_signs(sgn,+1,+1,-1);
    faces[4][2] = find_node_by_signs(sgn,+1,+1,+1);
    faces[4][3] = find_node_by_signs(sgn,-1,+1,+1);

    faces[5][0] = find_node_by_signs(sgn,-1,-1,-1);
    faces[5][1] = find_node_by_signs(sgn,-1,+1,-1);
    faces[5][2] = find_node_by_signs(sgn,-1,+1,+1);
    faces[5][3] = find_node_by_signs(sgn,-1,-1,+1);
}

/* ====== FaceEntry / SurfFaceEntry + merge join ====== */

static inline void dl_sort4(dl_node_id_t k[4]){
    if(k[1]<k[0]){dl_node_id_t t=k[0];k[0]=k[1];k[1]=t;}
    if(k[3]<k[2]){dl_node_id_t t=k[2];k[2]=k[3];k[3]=t;}
    if(k[2]<k[0]){dl_node_id_t t=k[0];k[0]=k[2];k[2]=t;}
    if(k[3]<k[1]){dl_node_id_t t=k[1];k[1]=k[3];k[3]=t;}
    if(k[2]<k[1]){dl_node_id_t t=k[1];k[1]=k[2];k[2]=t;}
}
static inline int dl_cmp4(const dl_node_id_t a[4], const dl_node_id_t b[4]){
    if(a[0]!=b[0]) return (a[0]<b[0])?-1:+1;
    if(a[1]!=b[1]) return (a[1]<b[1])?-1:+1;
    if(a[2]!=b[2]) return (a[2]<b[2])?-1:+1;
    if(a[3]!=b[3]) return (a[3]<b[3])?-1:+1;
    return 0;
}

typedef struct {
    dl_node_id_t key[4];
    int owner_elem;
    int owner_lface;
} dl_FaceEntry;

typedef struct {
    dl_node_id_t key[4];
    int surf_face_id;
} dl_SurfFaceEntry;

/* ---- タイブレーク付き比較（全順序化 → qsort の不安定性回避） ---- */
static int dl_cmp_face(const void* A, const void* B){
    const dl_FaceEntry* a=(const dl_FaceEntry*)A, *b=(const dl_FaceEntry*)B;
    int c = dl_cmp4(a->key, b->key);
    if (c) return c;
    if (a->owner_elem != b->owner_elem) return (a->owner_elem < b->owner_elem) ? -1 : +1;
    if (a->owner_lface != b->owner_lface) return (a->owner_lface < b->owner_lface) ? -1 : +1;
    return 0;
}
static int dl_cmp_sface(const void* A, const void* B){
    const dl_SurfFaceEntry* a=(const dl_SurfFaceEntry*)A, *b=(const dl_SurfFaceEntry*)B;
    int c = dl_cmp4(a->key, b->key);
    if (c) return c;
    if (a->surf_face_id != b->surf_face_id) return (a->surf_face_id < b->surf_face_id) ? -1 : +1;
    return 0;
}

/* faces_map[6][4] は呼び出し側で build_hex8_faces_from_sgn して渡す（グローバル排除） */
static dl_FaceEntry*
dl_build_vol_faces_sorted(const BBFE_DATA* fe, const int faces_map[6][4], int* nfaces_out)
{
    const int nfaces = fe->total_num_elems * 6;
    dl_FaceEntry* arr = (dl_FaceEntry*)malloc(sizeof(dl_FaceEntry)*nfaces);
    if(!arr){ fprintf(stderr,"[dl][ERR] malloc vol faces failed\n"); *nfaces_out=0; return NULL; }

    int k=0;
    for(int e=0;e<fe->total_num_elems;++e){
    for(int f=0;f<6;++f){
        dl_node_id_t key[4] = {
        (dl_node_id_t)fe->conn[e][faces_map[f][0]],
        (dl_node_id_t)fe->conn[e][faces_map[f][1]],
        (dl_node_id_t)fe->conn[e][faces_map[f][2]],
        (dl_node_id_t)fe->conn[e][faces_map[f][3]]
        };
        dl_sort4(key);
        for(int j=0;j<4;++j) arr[k].key[j]=key[j];
        arr[k].owner_elem=e; arr[k].owner_lface=f; ++k;
    }
    }
    qsort(arr, nfaces, sizeof(dl_FaceEntry), dl_cmp_face);
    *nfaces_out = nfaces;
    return arr;
    }

static dl_SurfFaceEntry* dl_build_surf_faces_sorted(
    const BBFE_DATA* surf, int* nsurf_out){
    const int nsurf = surf->total_num_elems;
    dl_SurfFaceEntry* arr = (dl_SurfFaceEntry*)malloc(sizeof(dl_SurfFaceEntry)*nsurf);
    if(!arr){ fprintf(stderr,"[dl][ERR] malloc surf faces failed\n"); *nsurf_out=0; return NULL; }

    for(int e=0;e<nsurf;++e){
    dl_node_id_t key[4] = {
        (dl_node_id_t)surf->conn[e][0],
        (dl_node_id_t)surf->conn[e][1],
        (dl_node_id_t)surf->conn[e][2],
        (dl_node_id_t)surf->conn[e][3]
    };
    dl_sort4(key);
    for(int j=0;j<4;++j) arr[e].key[j]=key[j];
    arr[e].surf_face_id=e;
    }
    qsort(arr, nsurf, sizeof(dl_SurfFaceEntry), dl_cmp_sface);
    *nsurf_out = nsurf;
    return arr;
}

/* owner_elem/owner_lface は [nSurfF] 確保・初期化済み（-1）前提 */
static int dl_build_owner_map_mergejoin(const dl_FaceEntry* volF, int nVolF,
                                        const dl_SurfFaceEntry* surfF, int nSurfF,
                                        int *owner_elem, int *owner_lface)
{
  int i=0, j=0, matches=0;

  while (i<nVolF && j<nSurfF){
    int c = dl_cmp4(volF[i].key, surfF[j].key);
    if (c==0){
      const int esu = surfF[j].surf_face_id;
      owner_elem[esu]  = volF[i].owner_elem;
      owner_lface[esu] = volF[i].owner_lface;
      ++matches;

      dl_node_id_t key0[4] = {volF[i].key[0],volF[i].key[1],volF[i].key[2],volF[i].key[3]};
      do { ++i; } while (i<nVolF && dl_cmp4(key0, volF[i].key)==0);
      ++j;
    } else if (c < 0){
      ++i;
    } else {
      ++j;
    }
  }
  return matches;
}

/* ====== ∇u（HEX8, セル一定LSQ） ====== */
static int dl_inv3x3(const double A[3][3], double Ai[3][3]){
  const double a=A[0][0], b=A[0][1], c=A[0][2];
  const double d=A[1][0], e=A[1][1], f=A[1][2];
  const double g=A[2][0], h=A[2][1], k=A[2][2];
  const double det = a*(e*k-f*h) - b*(d*k-f*g) + c*(d*h-e*g);

  if (fabs(det) < 1e-30) return 0;
  Ai[0][0]=(e*k-f*h)/det; Ai[0][1]=(c*h-b*k)/det; Ai[0][2]=(b*f-c*e)/det;
  Ai[1][0]=(f*g-d*k)/det; Ai[1][1]=(a*k-c*g)/det; Ai[1][2]=(c*d-a*f)/det;
  Ai[2][0]=(d*h-e*g)/det; Ai[2][1]=(b*g-a*h)/det; Ai[2][2]=(a*e-b*d)/det;

  return 1;
}

static void dl_grad_hex8_lsq(int ke, const BBFE_DATA* fe, const VALUES* vals, double G[3][3]){
    if (fe->local_num_nodes != 8){
        memset(G,0,sizeof(double)*9);
    return;
    }
    double xc[3]={0}, ubar[3]={0};
    for(int a=0;a<8;++a){
        int n=fe->conn[ke][a];
    for(int i=0;i<3;++i){ xc[i]+=fe->x[n][i]; ubar[i]+=vals->v[n][i]; }
    }
    for(int i=0;i<3;++i){ xc[i]/=8.0; ubar[i]/=8.0; }

    double A[3][3]={{0}}, b[3][3]={{0}};

    for(int a=0;a<8;++a){
        int n=fe->conn[ke][a];

        double r[3] = { fe->x[n][0]-xc[0], fe->x[n][1]-xc[1], fe->x[n][2]-xc[2] };

        A[0][0]+=r[0]*r[0]; A[0][1]+=r[0]*r[1]; A[0][2]+=r[0]*r[2];
        A[1][0]+=r[1]*r[0]; A[1][1]+=r[1]*r[1]; A[1][2]+=r[1]*r[2];
        A[2][0]+=r[2]*r[0]; A[2][1]+=r[2]*r[1]; A[2][2]+=r[2]*r[2];

        double du0=vals->v[n][0]-ubar[0], du1=vals->v[n][1]-ubar[1], du2=vals->v[n][2]-ubar[2];

        b[0][0]+=r[0]*du0; b[0][1]+=r[1]*du0; b[0][2]+=r[2]*du0;
        b[1][0]+=r[0]*du1; b[1][1]+=r[1]*du1; b[1][2]+=r[2]*du1;
        b[2][0]+=r[0]*du2; b[2][1]+=r[1]*du2; b[2][2]+=r[2]*du2;
    }
    const double eps_reg=1e-20; A[0][0]+=eps_reg; A[1][1]+=eps_reg; A[2][2]+=eps_reg;
    double Ai[3][3];

    if(!dl_inv3x3(A,Ai)){ memset(G,0,sizeof(double)*9); return; }

    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        G[i][j] = Ai[j][0]*b[i][0] + Ai[j][1]*b[i][1] + Ai[j][2]*b[i][2];

}

static void dl_compute_all_gradients_LSQ_hex8(
    const BBFE_DATA* fe, const VALUES* vals, double (*G_all)[3][3]){

    const int ne = fe->total_num_elems;

    #pragma omp parallel for if(ne>256) schedule(static) default(none) shared(fe, vals, G_all, ne)

    for(int ke=0; ke<ne; ++ke){
        dl_grad_hex8_lsq(ke, fe, vals, G_all[ke]);
    }
}

/* ====== QUAD4: 法線と J（dN/dxi, dN/deta から） ====== */
static inline void dl_quad4_normal(
    const double local_x[4][3], const double* dN_dxi, const double* dN_det, double normal[3])
{
    double dx_dxi[3]={0}, dx_det[3]={0};
    for(int a=0;a<4;++a){
        dx_dxi[0]+=dN_dxi[a]*local_x[a][0];
        dx_dxi[1]+=dN_dxi[a]*local_x[a][1];
        dx_dxi[2]+=dN_dxi[a]*local_x[a][2];
        dx_det[0]+=dN_det[a]*local_x[a][0];
        dx_det[1]+=dN_det[a]*local_x[a][1];
        dx_det[2]+=dN_det[a]*local_x[a][2];
    }
    normal[0] = dx_dxi[1]*dx_det[2] - dx_dxi[2]*dx_det[1];
    normal[1] = dx_dxi[2]*dx_det[0] - dx_dxi[0]*dx_det[2];
    normal[2] = dx_dxi[0]*dx_det[1] - dx_dxi[1]*dx_det[0];
}

/* ====== ベクトルユーティリティ ====== */
static inline double dl_dot3(const double a[3], const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static inline double dl_norm3(const double v[3]){return sqrt(dl_dot3(v,v));}
static inline void dl_normalize3(double v[3]){ double L=dl_norm3(v); if(L>0){v[0]/=L;v[1]/=L;v[2]/=L;} }

int calc_Cd_v(const BBFE_DATA* surf, const BBFE_DATA* fe, const BBFE_BASIS* basis, MONOLIS_COM* monolis_com,
                       const VALUES* vals,
                       double rho, double mu, double Uinf, double Aref,
                       const double eU_in[3], const double eP_in[3],
                       double* D_out, double* L_out, double* Cd_out, double* Cl_out)
{
    (void)monolis_com;
    double A_total = 0;

    int faces_map[6][4]; build_hex8_faces_from_sgn(SGN, faces_map);

    int nVolF=0, nSurfF=0;
    dl_FaceEntry*  volF  = dl_build_vol_faces_sorted(fe, faces_map, &nVolF);
    dl_SurfFaceEntry* sF = dl_build_surf_faces_sorted(surf, &nSurfF);
    if(!volF || !sF){ free(volF); free(sF); return -1; }

    int *owner = (int*)malloc(sizeof(int)*nSurfF);
    int *lface = (int*)malloc(sizeof(int)*nSurfF);
    if(!owner || !lface){ free(volF); free(sF); free(owner); free(lface); return -1; }
        for (int e=0; e<nSurfF; ++e){ owner[e] = -1; lface[e] = -1; }

    const int matched = dl_build_owner_map_mergejoin(volF,nVolF,sF,nSurfF,owner,lface);
    if (matched != nSurfF){
        int cnt=0;
        for (int e=0; e<nSurfF; ++e){ if (owner[e] < 0){ fprintf(stderr, "[dl] surf face %d has NO owner\n", e); ++cnt; } }
        fprintf(stderr, "[dl] %d / %d surf faces unmapped. Check node IDs / sorting / orientation.\n", cnt, nSurfF);
        free(volF); free(sF); free(owner); free(lface);
        return -1;
    }

    double (*G_all)[3][3] = (double(*)[3][3])malloc(sizeof(double)*fe->total_num_elems*9);
    if(!G_all){ free(volF); free(sF); free(owner); free(lface); return -1; }
    dl_compute_all_gradients_LSQ_hex8(fe, vals, G_all);

    double F_vis[3]={0,0,0};
    double F_pre[3]={0,0,0};
    const int nl = surf->local_num_nodes; /* QUAD4=4 */
    const int np = basis->num_integ_points;

    for(int e=0; e<surf->total_num_elems; ++e){
    double xloc[4][3];
    for(int a=0;a<nl;++a){
        int gid = surf->conn[e][a];
        xloc[a][0]=fe->x[gid][0]; xloc[a][1]=fe->x[gid][1]; xloc[a][2]=fe->x[gid][2];
    }
    double xf[3]={0,0,0};
    for(int a=0;a<nl;++a){ xf[0]+=xloc[a][0]; xf[1]+=xloc[a][1]; xf[2]+=xloc[a][2]; }
    xf[0]/=nl; xf[1]/=nl; xf[2]/=nl;

    const int ke = owner[e];
    double xc[3]={0,0,0};
    for(int a=0;a<8;++a){ int n=fe->conn[ke][a]; xc[0]+=fe->x[n][0]; xc[1]+=fe->x[n][1]; xc[2]+=fe->x[n][2]; }
    xc[0]/=8.0; xc[1]/=8.0; xc[2]/=8.0;
    double svec[3]={ xf[0]-xc[0], xf[1]-xc[1], xf[2]-xc[2] };

    const double (*G)[3] = G_all[ke];

    for(int p=0;p<np;++p){
        double nrm[3]; dl_quad4_normal(xloc, basis->dN_dxi[p], basis->dN_det[p], nrm);
        double J = dl_norm3(nrm);
        if(J<=1e-300) continue;
        double nh[3] = { nrm[0]/J, nrm[1]/J, nrm[2]/J };
        if(dl_dot3(svec,nh)<0){ nh[0]*=-1; nh[1]*=-1; nh[2]*=-1; }

        double p_ip=0.0;
        for(int a=0;a<nl;++a){ int gid=surf->conn[e][a]; p_ip += basis->N[p][a]*vals->p[gid]; }

        const double divu = G[0][0] + G[1][1] + G[2][2];
        double tau_n[3] = {0,0,0};
        for(int i=0;i<3;++i){
        double s = 0.0;
        for(int j=0;j<3;++j) s += (G[i][j] + G[j][i]) * nh[j];
        tau_n[i] = mu*s - (2.0/3.0)*mu*divu*nh[i];
        }

        /* τ = μ(∇u + ∇u^T) - (2/3) μ (∇·u) I */
        for(int i=0;i<3;++i){
            double s = 0.0;
            for(int j=0;j<3;++j) s += (G[i][j] + G[j][i]) * nh[j];
                tau_n[i] = mu*s - (2.0/3.0)*mu*divu*nh[i];
        }

        const double t_vis[3] = { tau_n[0], tau_n[1], tau_n[2] };

        const double t_pre[3] = { -p_ip*nh[0], -p_ip*nh[1], -p_ip*nh[2] };
        const double w = basis->integ_weight[p]*J;
        A_total += w;

        F_vis[0]+=t_vis[0]*w;  F_vis[1]+=t_vis[1]*w;  F_vis[2]+=t_vis[2]*w;
        F_pre[0]+=t_pre[0]*w;  F_pre[1]+=t_pre[1]*w;  F_pre[2]+=t_pre[2]*w;
        }
    }

    double eU[3]={eU_in[0],eU_in[1],eU_in[2]}, eP[3]={eP_in[0],eP_in[1],eP_in[2]};
    dl_normalize3(eU); dl_normalize3(eP);

    const double F_tot_fluid[3] = { F_vis[0], F_vis[1], F_vis[2] };
    const double F_body[3]      = { -F_tot_fluid[0], -F_tot_fluid[1], -F_tot_fluid[2] };

    const double D = dl_dot3(F_body, eU);
    const double L = dl_dot3(F_body, eP);
    const double qA = 0.5*rho*Uinf*Uinf * Aref;

    if (D_out)  *D_out  = D;
    if (L_out)  *L_out  = L;
    if (Cd_out) *Cd_out = D;
    if (Cl_out) *Cl_out = L;

    free(volF); free(sF); free(owner); free(lface); free(G_all);
    return 0;
}


#ifndef CHECK_ALLOC
#define CHECK_ALLOC(p) do { \
  if((p) == NULL){ \
    printf("alloc failed: %s (%s:%d)\n", #p, __FILE__, __LINE__); \
    abort(); \
  } \
} while(0)
#endif

static inline int max2(int a, int b){ return a > b ? a : b; }

double calc_internal_norm_1d(
        double*  in,       //input
        const int num1,
        const int num2)
{
    double norm = 0.0;

    for (int i = 0; i < num1; i++) {    
        norm += in[i] * in[i];
    }

    return norm;
}

void solver_fom_NR_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step)
{
	if(monolis_mpi_get_global_my_rank()==0){
    		printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
	}


    double* rvec = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));


        const double rel_tol_v = 1.0e-6;
const double abs_tol_v = 1.0e-12;
const double rel_tol_p = 1.0e-6;
const double abs_tol_p = 1.0e-12;
const double tiny      = 1.0e-30;
int max_iter_NR = 1;

    for(int it = 0; it < max_iter_NR; it++){
	    	if(monolis_mpi_get_global_my_rank()==0){
                	printf("\n%s ----------------- Time step %d : NR step %d ----------------\n", CODENAME, step, it);
		}
                monolis_clear_mat_value_R(&(sys.monolis));


        for(int i=0; i<sys.fe.total_num_nodes*4; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

                set_element_mat_NR_linear(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals));
                
                set_element_vec_NR_linear(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals));

		set_element_mat_NR_nonlinear(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals));

                set_element_vec_NR_nonlinear(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals));

                BBFE_sys_monowrap_set_Dirichlet_bc(
                                &(sys.monolis),
                                sys.fe.total_num_nodes,
                                4,
                                &(sys.bc_NR),
                                sys.monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys.fe.total_num_nodes*4; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }


                 ROM_monowrap_solve(
                                &(sys.monolis),
                                &(sys.mono_com),
                                sys.monolis.mat.R.X,
                                MONOLIS_ITER_BICGSAFE,
                                MONOLIS_PREC_DIAG,
                                20000,
                                sys.vals.mat_epsilon);

                BBFE_fluid_sups_renew_velocity(
                                sys.vals.delta_v,
                                sys.monolis.mat.R.X,
                                sys.fe.total_num_nodes);

                BBFE_fluid_sups_renew_pressure(
                                sys.vals.delta_p,
                                sys.monolis.mat.R.X,
                                sys.fe.total_num_nodes);

        update_velocity_pressure_NR(
                sys.vals.v,
                sys.vals.delta_v,
                sys.vals.p,
                sys.vals.delta_p,
                                sys.fe.total_num_nodes);

monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<sys.fe.total_num_nodes*4; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            //dx[i] = 0.0;
        }

	
set_element_mat_NR_Tezuer(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals));

                BBFE_sys_monowrap_set_Dirichlet_bc(
                                &(sys.monolis),
                                sys.fe.total_num_nodes,
                                4,
                                &(sys.bc_NR),
                                sys.monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        //if(i==0){
            for(int i=0; i<sys.fe.total_num_nodes*4; ++i) rvec[i] = sys.monolis.mat.R.B[i];
        //}


double norm_v = calc_internal_norm_2d(
    sys.vals.v,
    sys.mono_com.n_internal_vertex,
    3);

double norm_delta_v = calc_internal_norm_2d(
    sys.vals.delta_v,
    sys.mono_com.n_internal_vertex,
    3);

        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.mono_com.n_internal_vertex*4,
            1);

        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.mono_com.n_internal_vertex*4,
            1);



/* 圧力の L2 ノルム（内部自由度のみ） */
double norm_p = 0.0, norm_delta_p = 0.0;
for (int ii = 0; ii < sys.mono_com.n_internal_vertex; ++ii) {
    double pv  = sys.vals.p[ii];
    double dpv = sys.vals.delta_p[ii];
    norm_p       += pv  * pv;
    norm_delta_p += dpv * dpv;
}

/* L∞（最大変化量）：速度は3成分、圧力は1成分 */
double linf_delta_v_local = 0.0;
double linf_delta_p_local = 0.0;
for (int i_node = 0; i_node < sys.mono_com.n_internal_vertex; ++i_node) {
    for (int d = 0; d < 3; ++d) {
        double av = fabs(sys.vals.delta_v[i_node][d]);
        if (av > linf_delta_v_local) linf_delta_v_local = av;
    }
    double ap = fabs(sys.vals.delta_p[i_node]);
    if (ap > linf_delta_p_local) linf_delta_p_local = ap;
}

/* MPI で集約（L2 和：SUM、L∞：MAX） */
monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.mono_com.comm);



/* L∞は MAX で集約（Monolis に MAX が無ければ、自前で rank0 に gather→max でもOK） */
double linf_delta_v = linf_delta_v_local;
double linf_delta_p = linf_delta_p_local;
monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.mono_com.comm);
monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.mono_com.comm);

/* ルートを取って実ノルムに */
double nrm_v        = sqrt(norm_v);
double nrm_dv       = sqrt(norm_delta_v);
double nrm_p        = sqrt(norm_p);
double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);


/* 相対＋絶対の複合判定（ゼロ割回避） */
double denom_v = fmax(nrm_v,  tiny);
double denom_p = fmax(nrm_p,  tiny);

int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

/* ログ出力を見やすく */
if(monolis_mpi_get_global_my_rank()==0){
printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
       it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e   |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n",
       it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
}

/* 収束したら後処理 */
if (conv_v && conv_p) {
    /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
    double max_du = 0.0;
    for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {
        for (int d = 0; d < 3; ++d) {
            double du = fabs(sys.vals.v[ii][d] - sys.vals.v_old[ii][d]);
            if (du > max_du) max_du = du;
        }
    }
    if(monolis_mpi_get_global_my_rank()==0){
    	printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
    }
    ROM_BB_vec_copy_2d(
        sys.vals.v,
        sys.vals.v_old,
        sys.fe.total_num_nodes,
        3);
    break;
}

    }


    double* vec;
    vec = BB_std_calloc_1d_double(vec, 3*sys.fe.total_num_nodes);

    ROM_BB_vec_copy_2d_to_1d(
            sys.vals.v,
            vec,
            sys.fe.total_num_nodes);

    if(step%sys.vals.snapshot_interval == 0) {
            if(monolis_mpi_get_global_my_rank()==0){
    printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));
}
            if(monolis_mpi_get_global_comm_size() == 1){
                    ROM_std_hlpod_set_snapmat_nobc(
                                    sys.vals.p,
                                    &(sys.rom_p.hlpod_mat),
                                    sys.fe.total_num_nodes,
            1,
                                    (int)(step/sys.vals.snapshot_interval));
            }
            else{
                    ROM_std_hlpod_set_snapmat_nobc(
                                    sys.vals.p,
                                    &(sys.rom_p.hlpod_mat),
                                    sys.mono_com.n_internal_vertex,
            1,
                                    (int)(step/sys.vals.snapshot_interval));
            }

    }

    if(step%sys.vals.snapshot_interval == 0) {
            if(monolis_mpi_get_global_my_rank()==0){
    printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));
}
            if(monolis_mpi_get_global_comm_size() == 1){
                    //ROM_sys_hlpod_fe_set_snap_mat_para(
                    ROM_std_hlpod_set_snapmat_nobc(
                                    vec,
                                    &(sys.rom_v.hlpod_mat),
                                    //&(sys.bc),
                                    //&(sys.rom_sups.rom_bc), //要変更
                                    sys.fe.total_num_nodes,
                                    3,
                                    ((int)step/sys.vals.snapshot_interval));
            }
            else{
                    //ROM_sys_hlpod_fe_set_snap_mat_para(
                    ROM_std_hlpod_set_snapmat_nobc(
                                    vec,
                                    &(sys.rom_v.hlpod_mat),
                                    //&(sys.bc),
                                    //&(sys.rom_sups.rom_bc),
            sys.mono_com.n_internal_vertex,
                                    3,
                                    ((int)step/sys.vals.snapshot_interval));
            }

    }

    BB_std_free_1d_double(vec, 3*sys.fe.total_num_nodes);



}

