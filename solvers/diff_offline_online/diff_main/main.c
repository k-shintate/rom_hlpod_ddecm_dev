
#include "convdiff_core.h"

#include "hlpod_dataset.h"  //for ROM

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 3;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;


const double DELTA    = 1.0E-06;
const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

/*for ROM*/
const char*     ROM_ID_FINISH_TIME      = "#rom_finish_time";
const double  POD_DVAL_FINISH_TIME      = 1.0;
const char* POD_ID_OUTPUT_INTERVAL      = "#rom_output_interval";
const int POD_DVAL_OUTPUT_INTERVAL      = 1;
const char* POD_ID_SNAPSHOT_INTERVAL    = "#snapshot_interval";
const int POD_DVAL_SNAPSHOT_INTERVAL    = 1;


static const char* ROM_OUTPUT_FILENAME_VTK  = "rom_result_%06d.vtk";
static const char* POD_INPUT_FILENAME_COND  = "rom_cond.dat";

static const char* INPUT_DIRECTORYNAME_1STDD        = "parted.0/";
static const char* INPUT_DIRECTORYNAME_2NDD         = "parted.1/";
static const char* INPUT_DIRECTORYNAME_METAGRAPH    = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH         = "metagraph.dat";

static const char* OPTION_NUM_MODES = "-nm";
static const char* OPTION_NUM_1STDD = "-nd";
static const char* OPTION_PADAPTIVE = "-pa";
static const char* OPTION_SOLVER_TYPE = "-st";
/********/

typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double dt;
	double finish_time;
	int    output_interval;

	double* T;
	double* error;
	double* theo_sol;

    /*for ROM input data*/
    double  rom_finish_time;
    int     rom_output_interval;
    int     snapshot_interval;

    int     num_cases;
    double* density_cases;
    double* viscosity_cases;

} VALUES;


typedef struct
{
	const char* directory;

} CONDITIONS;

/*for ROM input data*/
typedef struct
{
    int     num_modes;		    //1領域当たりの最大基底本数
    double  rom_epsilon;	    //基底取得閾値
    int     num_subdomains;     //POD計算領域数
    int     solver_type;
} ROM_PRM;

typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;
	BBFE_BC      bc;
	MONOLIS      monolis;
	MONOLIS_COM  monolis_com;

	CONDITIONS   cond;
	VALUES       vals;

	MONOLIS      monolis0;          // for nonsteady analysis


    /*for ROM*/
	MONOLIS      monolis_rom;
	MONOLIS      monolis_rom0;
	
	VALUES       vals_rom;

	MONOLIS_COM  mono_com_rom;
	MONOLIS_COM  mono_com0;             //初期化関数
	MONOLIS_COM  mono_com_rom_solv;	    //線形ソルバ用コミュニケータ

	ROM_PRM		 rom_prm;		        //input データ用
	ROM		  	 rom;

} FE_SYSTEM;


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z ) * (0.5*t + sin( 1.0*t ));

	return val;
}


void manusol_get_conv_vel(
		double v[3],
		double x[3])
{
	v[0] = 1.0 + x[0]*x[0];
	v[1] = 1.0 + x[1]*x[1];
	v[2] = 1.0 + x[2]*x[2];
}


double manusol_get_mass_coef(
		double x[3])
{
	double val = 1.0;

	return val;
}


double manusol_get_diff_coef(
		double x[3])
{
	double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));

	return val;
}

double manusol_get_source(
		double x[3],
		double t,
		double a,
		double v[3],
		double k)
{
	double val = 0.0;

	return val;
}

/*
double manusol_get_source(
		double x[3],
		double t,
		double a,
		double v[3],
		double k)
{
	double invD  = 1.0/DELTA;
	double invD2 = invD*invD;

	double dT_dt = ( manusol_get_sol(x[0], x[1], x[2], t+DELTA) -
		manusol_get_sol(x[0], x[1], x[2], t-DELTA) ) * (invD/2.0);

	double dk_dx[3];
	double x_p[3];  double x_m[3];
	x_p[0] = x[0]+DELTA;  x_p[1] = x[1];  x_p[2] = x[2];
	x_m[0] = x[0]-DELTA;  x_m[1] = x[1];  x_m[2] = x[2];
	dk_dx[0] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) ) * (invD/2.0);
	x_p[0] = x[0];  x_p[1] = x[1]+DELTA;  x_p[2] = x[2];
	x_m[0] = x[0];  x_m[1] = x[1]-DELTA;  x_m[2] = x[2];
	dk_dx[1] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) ) * (invD/2.0);
	x_p[0] = x[0];  x_p[1] = x[1];  x_p[2] = x[2]+DELTA;
	x_m[0] = x[0];  x_m[1] = x[1];  x_m[2] = x[2]-DELTA;
	dk_dx[2] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) )* (invD/2.0);

	double dT_dx[3];
	dT_dx[0] =  ( manusol_get_sol(x[0]+DELTA, x[1], x[2], t) -
		manusol_get_sol(x[0]-DELTA, x[1], x[2], t) ) * (invD/2.0);
	dT_dx[1] =  ( manusol_get_sol(x[0], x[1]+DELTA, x[2], t) -
		manusol_get_sol(x[0], x[1]-DELTA, x[2], t) ) * (invD/2.0);
	dT_dx[2] =  ( manusol_get_sol(x[0], x[1], x[2]+DELTA, t) -
		manusol_get_sol(x[0], x[1], x[2]-DELTA, t) ) * (invD/2.0);

	double d2T_dx2[3];
	d2T_dx2[0] = (
			    manusol_get_sol(x[0]+DELTA, x[1], x[2], t) +
			    manusol_get_sol(x[0]-DELTA, x[1], x[2], t) -
			2.0*manusol_get_sol(x[0]      , x[1], x[2], t) ) * invD2;
	d2T_dx2[1] = (
			    manusol_get_sol(x[0], x[1]+DELTA, x[2], t) +
			    manusol_get_sol(x[0], x[1]-DELTA, x[2], t) -
			2.0*manusol_get_sol(x[0], x[1]      , x[2], t) ) * invD2;
	d2T_dx2[2] = (
			    manusol_get_sol(x[0], x[1], x[2]+DELTA, t) +
			    manusol_get_sol(x[0], x[1], x[2]-DELTA, t) -
			2.0*manusol_get_sol(x[0], x[1], x[2]      , t) ) * invD2;

	double val =
		a * (dT_dt + v[0]*dT_dx[0] + v[1]*dT_dx[1] + v[2]*dT_dx[2]) -
		(dk_dx[0]*dT_dx[0] + dk_dx[1]*dT_dx[1] + dk_dx[2]*dT_dx[2]) -
		k * (d2T_dx2[0] + d2T_dx2[1] + d2T_dx2[2]);

	return val;
}
*/


void manusol_set_theo_sol(
		BBFE_DATA* fe,
		double*  theo_sol,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		theo_sol[i] = manusol_get_sol(fe->x[i][0], fe->x[i][1], fe->x[i][2], t);
	}
}


void manusol_set_source(
		BBFE_DATA* fe,
		double*  source,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double a;  double v[3];  double k;
		a = manusol_get_mass_coef(fe->x[i]);
		manusol_get_conv_vel(v, fe->x[i]);
		k = manusol_get_diff_coef(fe->x[i]);
		source[i] = manusol_get_source(fe->x[i], t, a, v, k);
	}
}


void manusol_set_init_value(
		BBFE_DATA* fe,
		double* T)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		T[i] = manusol_get_sol(fe->x[i][0], fe->x[i][1], fe->x[i][2], 0.0);
	}
}


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->T        = BB_std_calloc_1d_double(vals->T,     total_num_nodes);
	vals->error    = BB_std_calloc_1d_double(vals->error, total_num_nodes);
	vals->theo_sol = BB_std_calloc_1d_double(vals->error, total_num_nodes);
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
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "temperature");

	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->theo_sol, vals->T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "abs_error");
	BB_vtk_write_point_vals_scalar(fp, vals->theo_sol, fe->total_num_nodes, "theoretical");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			filename,
			sys->cond.directory,
			t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.T,
			filename,
			sys->cond.directory);

	/**** for manufactured solution ****/
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_sou);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);

	double L2_error = BBFE_elemmat_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis_com),
			t,
			sys->vals.T,
			manusol_get_sol);

	printf("%s L2 error: %e\n", CODENAME, L2_error);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = BBFE_sys_write_add_fopen(fp, "l2_error.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error);
		fclose(fp);
	}

	BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}


void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals)
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

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);

					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
*/
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
*/
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);
/*
					val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis,
						fe->conn[e][i], fe->conn[e][j], 0, 0,
						integ_val);
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


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals,
		double       t)
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
		BBFE_elemmat_set_local_array_scalar(local_T, fe, vals->T, e);

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
				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], a_ip[p]);
/*
				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
*/
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
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


/*for ROM*/
void ROM_read_args(
    int 		argc,
    char* 		argv[],
    ROM_PRM*    rom_prm)
{
	int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_1STDD);
    if(num == -1) {
		printf("\nargs error num_subdomains");
		exit(1);
    }
    else {
        rom_prm->num_subdomains = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_MODES);
    if(num == -1) {
		printf("\nargs error num_modes");
		exit(1);
    }
    else {
		rom_prm->num_modes = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_PADAPTIVE);
    if(num == -1) {
		printf("\nargs error rom_epsilon");
		exit(1);
    }
    else {
        rom_prm->rom_epsilon = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_SOLVER_TYPE);
    if(num == -1) {
		printf("\nargs error solver_type");
		exit(1);
    }
    else {
        rom_prm->solver_type = atof(argv[num+1]);
    }

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);

}

const char* ROM_std_hlpod_get_parted_file_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_2NDD);
        return fname;
    }
    else if(solver_type==4){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}

const char* ROM_std_hlpod_get_metagraph_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_1STDD, INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_METAGRAPH ,INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==4){
        return NULL;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}

void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type)
{
    rom->hlpod_vals.num_1st_subdomains = num_subdomains;
    rom->hlpod_vals.num_2nd_subdomains = num_subdomains;
    rom->hlpod_vals.num_modes_pre = num_modes;
    rom->hlpod_vals.rom_epsilon = 1 - rom_epsilon;

    rom->hlpod_vals.bool_global_mode = false;
    
    if(solver_type == 4){
        rom->hlpod_vals.bool_global_mode = true;
    }

}

void ROM_offline_assign_default_values(
		VALUES*     vals)
{
	vals->snapshot_interval  = POD_DVAL_SNAPSHOT_INTERVAL;
}


void ROM_offline_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, POD_ID_SNAPSHOT_INTERVAL,  vals->snapshot_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_offline_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_offline_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_int_p(
				&(vals->snapshot_interval), filename, POD_ID_SNAPSHOT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_offline_print_all_values(vals);

	printf("\n");
}


void ROM_online_assign_default_values(
		VALUES*     vals)
{
	vals->rom_finish_time  = POD_DVAL_FINISH_TIME;
	vals->rom_output_interval  = POD_DVAL_OUTPUT_INTERVAL;
}


void ROM_online_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %e\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_finish_time);
	printf("%s %s: %d\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_output_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_online_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_online_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_double_p(
				&(vals->rom_finish_time), filename, ROM_ID_FINISH_TIME, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->rom_output_interval), filename, POD_ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_online_print_all_values(vals);

	printf("\n");
}


void ROM_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*   hlpod_vals,
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
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");

	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "abs_error");
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}


void ROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, ROM_OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	ROM_output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			&(sys->rom.hlpod_vals),
			filename,
			sys->cond.directory,
			t);
	
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals_rom.T,
			filename,
			sys->cond.directory);

	double L2_error = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis_com),
			t,
			sys->vals_rom.T,
			sys->vals.T);

	printf("%s L2 error: %e\n", CODENAME, L2_error);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_rom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error);
		fclose(fp);
	}

}


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	manusol_set_init_value(&(sys.fe), sys.vals.T);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));

	monolis_initialize(&(sys.monolis0));

	monolis_com_initialize_by_parted_files(
			&(sys.monolis_com),
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_R(
			&(sys.monolis0),
			sys.fe.total_num_nodes,
			sys.fe.local_num_nodes,
			1,
			sys.fe.total_num_elems,
			sys.fe.conn);

	set_element_mat(
			&(sys.monolis0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));


    /*for ROM *****************************************/
	
    /*for ROM input data*/
	ROM_read_args(argc, argv, &(sys.rom_prm));

    ROM_set_param(
            &(sys.rom),
            sys.rom_prm.num_subdomains,
            sys.rom_prm.num_modes,
            sys.rom_prm.rom_epsilon,
            sys.rom_prm.solver_type);
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            1,
            &(sys.rom.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
    /******************/

    /*only for offline*/
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            1,
			1);
    /******************/

    /**********************************************/

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
		monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
		monolis_clear_mat_value_rhs_R(&(sys.monolis));

		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals),
				t);

		manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
		BBFE_manusol_set_bc_scalar(
				&(sys.fe),
				&(sys.bc),
				sys.vals.theo_sol,
				t);

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				BLOCK_SIZE,
				&(sys.bc),
				sys.monolis.mat.R.B);

		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.monolis_com),
				sys.vals.T,
				MONOLIS_ITER_BICGSTAB,
				MONOLIS_PREC_DIAG,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);


		if(step%sys.vals.snapshot_interval == 0) {
			printf("Set modes T: %d\n", (int)(step/sys.vals.snapshot_interval));

			if(monolis_mpi_get_global_comm_size() == 1){
				ROM_sys_hlpod_fe_set_snap_mat_para(
						sys.vals.T,
						&(sys.rom.hlpod_mat),
						&(sys.bc),
						&(sys.rom.rom_bc),
						sys.fe.total_num_nodes,
						1,
						((int)step/sys.vals.snapshot_interval));
			}
			else{
				ROM_sys_hlpod_fe_set_snap_mat_para(
						sys.vals.T,
						&(sys.rom.hlpod_mat),
						&(sys.bc),
						&(sys.rom.rom_bc),
                        sys.monolis_com.n_internal_vertex,
						1,
						((int)step/sys.vals.snapshot_interval));
			}

		}

		/**********************************************/

		if(step%sys.vals.output_interval == 0) {
			output_files(&sys, file_num, t);
			file_num += 1;
		}

	}

    /*for ROM *****************************************/
    ROM_std_hlpod_set_pod_modes(
		&(sys.rom),
		sys.fe.total_num_nodes,
		sys.monolis_com.n_internal_vertex,
		1,
		"pod_modes",
		sys.cond.directory);
    
	/*only for online*/
    read_calc_conditions(&(sys.vals_rom), sys.cond.directory);                  //set vals

    memory_allocation_nodal_values(&(sys.vals_rom), sys.fe.total_num_nodes);	//set vals

	ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom.hlpod_vals), sys.fe.total_num_nodes, 1);

    ROM_online_read_calc_conditions(&(sys.vals), sys.cond.directory);

	monolis_initialize(&(sys.monolis_rom0));
    monolis_initialize(&(sys.monolis_rom));
	monolis_com_initialize_by_self(&(sys.mono_com_rom));
	monolis_com_initialize_by_self(&(sys.mono_com_rom_solv));

	ROM_std_hlpod_set_monolis_comm(
			&(sys.monolis_com),
			&(sys.mono_com_rom),
			&(sys.mono_com_rom_solv),
			INPUT_DIRECTORYNAME_METAGRAPH,
            INPUT_FILENAME_METAGRAPH,
            sys.rom_prm.solver_type,
			sys.cond.directory);

	if(monolis_mpi_get_global_my_rank() == 0){
		fp = ROM_BB_write_fopen(fp, "l2_error_rom.txt", sys.cond.directory);
		fclose(fp);
	}

	if(monolis_mpi_get_global_my_rank() == 0){
        ROM_std_hlpod_write_solver_prm_fopen("fem_solver_prm", sys.cond.directory);
		ROM_std_hlpod_write_solver_prm_fopen("pod_solver_prm", sys.cond.directory);
	}


    ROM_std_hlpod_read_pod_modes(
            &(sys.rom),
            sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            "pod_modes",
            sys.cond.directory);

	ROM_sys_hlpod_fe_write_pod_modes_vtk(
            &(sys.fe),
            &(sys.rom),
            sys.monolis_com.n_internal_vertex,
            10,		//10本まで可視化
            1,
            "pod_modes_vtk/pod_modes.vtk",
            sys.cond.directory);

	ROM_std_hlpod_online_pre(
            &(sys.monolis_rom0),
            &(sys.monolis_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom),
            sys.fe.total_num_nodes,
            metagraph_name,
            sys.cond.directory);
    
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_rom0), &(sys.monolis_rom));
    monolis_com_initialize_by_self(&(sys.mono_com0));
    /******************/

    /**********************************************/
    ROM_BB_vec_copy(sys.vals.T, sys.vals_rom.T, sys.fe.total_num_nodes);
    
    file_num = 0;
	while (t < sys.vals.rom_finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
        double calctime_fem_t1 = monolis_get_time();
		monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
		monolis_clear_mat_value_rhs_R(&(sys.monolis));

		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals),
				t);

		manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
		BBFE_manusol_set_bc_scalar(
				&(sys.fe),
				&(sys.bc),
				sys.vals.theo_sol,
				t);

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				BLOCK_SIZE,
				&(sys.bc),
				sys.monolis.mat.R.B);

		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.monolis_com),
				sys.vals.T,
				MONOLIS_ITER_CG,
				MONOLIS_PREC_DIAG,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);
        
        double calctime_fem_t2 = monolis_get_time();
		/**********************************************/

        /****************** ROM solver ****************/
        printf("\n%s ----------------- step ROM %d ----------------\n", CODENAME, step);
        double calctime_rom_t1 = monolis_get_time();
		monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
		monolis_clear_mat_value_rhs_R(&(sys.monolis));
    	monolis_clear_mat_value_R(&(sys.monolis_rom));

		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals_rom),
				t);

		manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
		BBFE_manusol_set_bc_scalar(
				&(sys.fe),
				&(sys.bc),
				sys.vals.theo_sol,
				t);

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				BLOCK_SIZE,
				&(sys.bc),
				sys.monolis.mat.R.B);
        
        /*for ROM*/
        ROM_std_hlpod_calc_reduced_mat(
                &(sys.monolis),
                &(sys.monolis_rom),
                &(sys.monolis_com),
                &(sys.mono_com0),
                &(sys.mono_com_rom_solv),
                &(sys.rom),
                sys.fe.total_num_nodes,
                1);
	
        ROM_std_hlpod_solve_ROM(
                &(sys.monolis),
                &(sys.monolis_rom),
                &(sys.mono_com_rom_solv),
                &(sys.rom),
                sys.fe.total_num_nodes,
                1,
                sys.vals.mat_max_iter,
                sys.vals.mat_epsilon,
                MONOLIS_ITER_CG,
                MONOLIS_PREC_DIAG);

        ROM_sys_hlpod_fe_add_Dbc(
                sys.rom.hlpod_vals.sol_vec,
                &(sys.bc),
                sys.fe.total_num_nodes,
                1);

        monolis_mpi_update_R(&(sys.monolis_com), sys.fe.total_num_nodes, 1, sys.rom.hlpod_vals.sol_vec);

        ROM_BB_vec_copy(sys.rom.hlpod_vals.sol_vec, sys.vals_rom.T, sys.fe.total_num_nodes);
        double calctime_rom_t2 = monolis_get_time();
        /************/

		/**********************************************/

		if(step%sys.vals.output_interval == 0) {
			output_files(&sys, file_num, t);
			file_num += 1;
		}

        if(step%sys.vals.output_interval == 0) {
			ROM_output_files(&sys, file_num, t);
                        
            ROM_std_hlpod_write_solver_prm(&(sys.monolis), t, "fem_solver_prm/" , sys.cond.directory);
			ROM_std_hlpod_write_solver_prm(&(sys.monolis_rom), t, "pod_solver_prm/", sys.cond.directory);

            ROM_std_hlpod_output_calc_time(calctime_fem_t1-calctime_fem_t2, t,
					"calctime/time_fem.txt", sys.cond.directory);
            ROM_std_hlpod_output_calc_time(calctime_rom_t1-calctime_rom_t2, t,
					"calctime/time_rom.txt", sys.cond.directory);

			file_num += 1;
		}

	}


	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis0));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
