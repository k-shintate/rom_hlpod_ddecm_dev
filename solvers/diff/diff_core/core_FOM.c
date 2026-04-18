
#include "core_FOM.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
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
const char*         ID_DENSITY  = "#density";
const double      DVAL_DENSITY  = 1000.0;
const char*       ID_VISCOSITY  = "#viscosity";
const double    DVAL_VISCOSITY  = 1.0;

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%06d.vtk";


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z ) * sin( 1.0*t );

	return val;
}

double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z );

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
	//double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));

	double val = (1 + x[0] * x[1]* x[2]);

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
	double val = 0.0;

	return val;
}
*/

double manusol_get_source(
                double x[3],
                double t,
                double a,
                double v[3],
                double k)
{
        double val = 0.0;

        if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 3.70 && x[1] < 3.8 &&x[2] > 3.7 && x[2] < 3.8){
                val = 1000 * (1 + sin(t));
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 2.45 && x[1] < 2.55 &&x[2] > 2.45 && x[2] < 2.55){
                val = 1000 * (1 + sin(t));
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 1.20 && x[1] < 1.3 &&x[2] > 1.2 && x[2] < 1.3){
                val = 1000 * (1 + sin(t));
        }

        return val;
}

double manusol_get_source_without_time(
                double x[3],
                double t,
                double a,
                double v[3],
                double k)
{
        double val = 0.0;

        if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 3.70 && x[1] < 3.8 &&x[2] > 3.7 && x[2] < 3.8){
                val = 1000;
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 2.45 && x[1] < 2.55 &&x[2] > 2.45 && x[2] < 2.55){
                val = 1000;
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 1.20 && x[1] < 1.3 &&x[2] > 1.2 && x[2] < 1.3){
                val = 1000;
        }

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

void manusol_set_theo_sol_without_time(
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
/*
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_sou);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);
*/
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

//	BB_std_free_1d_double(source, sys->fe.total_num_nodes);
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

void set_element_vec_source(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		//VALUES*      vals,
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
		//BBFE_elemmat_set_local_array_scalar(local_T, fe, vals->T, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			//T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source_without_time(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;

				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	//BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	//BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
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

					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);

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
void ROM_set_D_bc_rhs_para(
		//MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*      hlpod_mat,
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

    for(int e = 0; e < fe->total_num_nodes; e++) {
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
                    hlpod_mat->VTf_D_bc[k1] += - val;
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

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
    
    monolis_copy_mat_value_R(&(sys.monolis0), &(sys.monolis));
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
}

//スナップショットの収集
void solver_fom_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

    //monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
    monolis_copy_mat_value_R(&(sys.monolis0), &(sys.monolis));
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

}
