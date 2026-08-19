
#include "core_FOM.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 100000;
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


/*
 * ============================================================
 * Edge-heated transient-diffusion HROM stress benchmark
 * ============================================================
 *
 * Intended component layout (dper=1, --copies 8, --depth-copies 1):
 *
 *   heated face
 *      |
 *      v
 *   x=0                                                x=4
 *    | D0 | D1 | D2 | D3 || D4 | D5 | D6 | D7 |
 *
 * Coefficient distribution:
 *
 *   left  half (D0-D3, x <= 2): q(x) = D_out = O(1)
 *   right half (D4-D7, x >  2): q(x) = D_in
 *
 * with
 *
 *   a(x) = q(x)   (mass coefficient),
 *   k(x) = q(x)   (diffusion coefficient).
 *
 * Governing equation:
 *
 *   a(x) T_t - div(k(x) grad T) = 0.
 *
 * Initial condition:
 *
 *   T(x,0) = T0.
 *
 * Boundary condition REQUIRED for this benchmark:
 *
 *   x = xmin : T = T0 + Q0*(1-exp(-t/tau))
 *   all other external faces : natural zero-flux (adiabatic)
 *
 * Therefore D_bc.dat must contain ONLY the global xmin face.  With the
 * repeat mesh generator use, for example,
 *
 *   --no-auto-bc --bc xmin:1
 *
 * There is no manufactured volumetric source.  Heat enters only from the
 * left end face and diffuses gradually through the component chain.
 *
 * Offline training uses D_out=1 and small D_in, so element weak-form
 * contributions on D4-D7 are small.  The time horizon should be long enough
 * for the FOM temperature front to reach the right half, so POD still sees
 * those states.  Online uses the SAME PDE family with D_in increased toward
 * 1 while the offline POD/HROM data are frozen.
 *
 * manusol_get_sol() is retained only as the time-dependent Dirichlet value
 * provider required by the existing BBFE_manusol_set_bc_scalar() path.  It is
 * NOT an analytical solution in the interior.  FOM output must be used as
 * the reference for ROM/HROM error evaluation.
 * ============================================================
 */

#define HROM_STRESS_SPLIT_X          2.0

void manusol_set_param(const MANUSOL_PARAM* prm)
{
    g_manusol_param = *prm;

    if (!(g_manusol_param.D_out > 0.0) ||
        !(g_manusol_param.D_in  > 0.0) ||
        !(g_manusol_param.tau   > 0.0)) {
        fprintf(stderr,
            "ERROR: edge-heated transient-diffusion stress benchmark requires "
            "D_out>0, D_in>0, tau>0 (D_out=%.6e D_in=%.6e tau=%.6e)\n",
            g_manusol_param.D_out,
            g_manusol_param.D_in,
            g_manusol_param.tau);
        exit(EXIT_FAILURE);
    }

    fprintf(stderr,
        "[TRANSIENT-DIFFUSION-EDGE-STRESS] "
        "left(D0-D3)=%.6e right(D4-D7)=%.6e "
        "T0=%.6e dT=%.6e ramp_tau=%.6e split_x=%.6e\n",
        g_manusol_param.D_out,
        g_manusol_param.D_in,
        g_manusol_param.T0,
        g_manusol_param.Q0,
        g_manusol_param.tau,
        HROM_STRESS_SPLIT_X);
    fprintf(stderr,
        "[TRANSIENT-DIFFUSION-EDGE-STRESS] REQUIRED BC: xmin Dirichlet only; "
        "all other faces natural zero-flux.\n");
}

static double hrom_stress_boundary_temperature(const double t)
{
    return g_manusol_param.T0
         + g_manusol_param.Q0
         * (1.0 - exp(-t/g_manusol_param.tau));
}

static double hrom_stress_q(const double x[3])
{
    if (x[0] <= HROM_STRESS_SPLIT_X) {
        return g_manusol_param.D_out;
    }
    return g_manusol_param.D_in;
}

double manusol_get_sol(
        double x,
        double y,
        double z,
        double t)
{
    /*
     * Compatibility callback for BBFE_manusol_set_bc_scalar().
     * D_bc.dat must contain only xmin, so only boundary nodes consume this
     * value.  This is intentionally NOT the interior analytical solution.
     */
    (void)x;
    (void)y;
    (void)z;
    return hrom_stress_boundary_temperature(t);
}

double manusol_get_diff_coef(double x[3])
{
    return hrom_stress_q(x);
}

double manusol_get_source(
        double x[3],
        double t,
        double a,
        double v[3],
        double k)
{
    (void)x;
    (void)t;
    (void)a;
    (void)v;
    (void)k;
    return 0.0;
}

double manusol_get_sol_without_time(
        double x,
        double y,
        double z,
        double t)
{
    (void)x;
    (void)y;
    (void)z;
    (void)t;
    return g_manusol_param.T0;
}

void manusol_get_conv_vel(
        double v[3],
        double x[3])
{
    (void)x;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}

double manusol_get_mass_coef(double x[3])
{
    return hrom_stress_q(x);
}

/*
double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z ) * sin( 1.0*t );

	return val;
}
*/

/*
double manusol_get_diff_coef(
		double x[3])
{
	double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));

	//double val = (1 + x[0] * x[1]* x[2]);

	return val;
}
*/

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

/*
double manusol_get_source(
                double x[3],
                double t,
                double a,
                double v[3],
                double k)
{
        double val = 0.0;
/*
        if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 3.70 && x[1] < 3.8 &&x[2] > 3.7 && x[2] < 3.8){
                val = 1000 * (1 + sin(t));
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 2.45 && x[1] < 2.55 &&x[2] > 2.45 && x[2] < 2.55){
                val = 1000 * (1 + sin(t));
        }
        else if(x[0] > 2.45 && x[0] < 2.55 && x[1] > 1.20 && x[1] < 1.3 &&x[2] > 1.2 && x[2] < 1.3){
                val = 1000 * (1 + sin(t));
        }
*/
/*
        return val;
}
*/

double manusol_get_source_without_time(
        double x[3],
        double t,
        double a,
        double v[3],
        double k)
{
    return manusol_get_source(x,t,a,v,k);
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
		theo_sol[i] = manusol_get_sol_without_time(fe->x[i][0], fe->x[i][1], fe->x[i][2], t);
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
    for(int i = 0; i < fe->total_num_nodes; ++i) {
        T[i] = g_manusol_param.T0;
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
        VALUES*          vals,
        const char*      filename,
        const char*      directory,
        double           t)
{
    (void)t;

    FILE* fp;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    switch(fe->local_num_nodes) {
        case 4:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
            break;
        case 8:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
            break;
    }

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    BB_vtk_write_point_vals_scalar(
        fp, vals->T, fe->total_num_nodes, "temperature");

    /* Visualize the fixed left/right coefficient distribution as well. */
    double* q = BB_std_calloc_1d_double(q, fe->total_num_nodes);
    for(int i = 0; i < fe->total_num_nodes; ++i) {
        q[i] = manusol_get_diff_coef(fe->x[i]);
    }
    BB_vtk_write_point_vals_scalar(
        fp, q, fe->total_num_nodes, "diffusion_mass_scale");
    BB_std_free_1d_double(q, fe->total_num_nodes);

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

    /*
     * No analytical interior solution is used in the edge-heating benchmark.
     * Compare ROM/HROM against a saved FOM reference in post-processing.
     */
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

void set_element_vec_mass(
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
                

				val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], a_ip[p]);
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

void set_element_mat_mass(
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

    for(int e = 0; e < fe->total_num_elems; e++) {
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


static double dot3(const double* a, const double* b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void set_element_mat_NR(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double*  a_ip;
	double** v_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	for(int e = 0; e < fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p = 0; p < np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int i = 0; i < nl; i++) {
			for(int j = 0; j < nl; j++) {
				for(int p = 0; p < np; p++) {
					val_ip[p] = 0.0;

					/* diffusion */
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i],
							fe->geo[e][p].grad_N[j],
							k_ip[p]);

					/* mass / transient */
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
							basis->N[p][i],
							basis->N[p][j],
							a_ip[p]);

					(void)h_e;
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

	BB_std_free_1d_double(val_ip, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_x, nl, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}


void set_element_vec_NR(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		double      t)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	double*  local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;
	double** v_ip;
	double*  k_ip;
	double*  a_ip;
	double*  T_ip;
	double*  f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int e = 0; e < fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, vals->T, e);

		for(int p = 0; p < np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i = 0; i < nl; i++) {
			for(int p = 0; p < np; p++) {
				val_ip[p] = 0.0;

				double grad_T[3] = {0.0, 0.0, 0.0};
				for(int j = 0; j < nl; j++) {
					grad_T[0] += local_T[j] * fe->geo[e][p].grad_N[j][0];
					grad_T[1] += local_T[j] * fe->geo[e][p].grad_N[j][1];
					grad_T[2] += local_T[j] * fe->geo[e][p].grad_N[j][2];
				}

				/* residual: mass */
				//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
				//		basis->N[p][i], T_ip[p], a_ip[p]);

				/* residual: diffusion */
				val_ip[p] += k_ip[p] * dot3(fe->geo[e][p].grad_N[i], grad_T);

				/* residual: source */
				val_ip[p] -= BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);

				(void)h_e;
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			/* solve: j * delta_u = -r */
			monolis->mat.R.B[fe->conn[e][i]] -= integ_val;
		}
	}

	BB_std_free_1d_double(val_ip, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_1d_double(local_T, nl);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
    
    /* Parameter-dependent a(x),k(x): rebuild the matrix for current D_in/D_out. */
    monolis_clear_mat_value_R(&(sys.monolis));
    monolis_clear_mat_value_rhs_R(&(sys.monolis));
    set_element_mat(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

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


//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom_NR(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
    
    //monolis_copy_mat_value_R(&(sys.monolis0), &(sys.monolis));
    monolis_clear_mat_value_R(&(sys.monolis));
    monolis_clear_mat_value_rhs_R(&(sys.monolis));

    double* T_old = BB_std_calloc_1d_double(T_old, sys.fe.total_num_nodes);

    for(int i = 0; i < sys.fe.total_num_nodes; i++){
        T_old[i] = sys.vals.T[i];
    }

    set_element_mat_NR(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

    set_element_vec_NR(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals),
            t);

    double* theo_sol_old = BB_std_calloc_1d_double(theo_sol_old, sys.fe.total_num_nodes);
    manusol_set_theo_sol(&(sys.fe), theo_sol_old, t - sys.vals.dt);

    manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);

    for(int i = 0; i < sys.fe.total_num_nodes; i++){
        sys.vals.theo_sol[i] -= theo_sol_old[i];
    }

    //manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
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

    for(int i = 0; i < sys.fe.total_num_nodes; i++){
        sys.vals.T[i] += T_old[i];
    }
}


//スナップショットの収集
void solver_fom_collect_snapmat_test(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    double tt = monolis_get_time_global_sync();
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

    //monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
    monolis_copy_mat_value_R(&(sys.monolis0), &(sys.monolis));
    monolis_clear_mat_value_rhs_R(&(sys.monolis));

    double t4 = monolis_get_time_global_sync();
    printf("now4");

    set_element_vec(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals),
            t);
    
    double t3 = monolis_get_time_global_sync();
    printf("now3");

    manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
    BBFE_manusol_set_bc_scalar(
            &(sys.fe),
            &(sys.bc),
            sys.vals.theo_sol,
            t);
    double t2 = monolis_get_time_global_sync();
    printf("now2");
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.monolis),
            sys.fe.total_num_nodes,
            BLOCK_SIZE,
            &(sys.bc),
            sys.monolis.mat.R.B);

    double t1 = monolis_get_time_global_sync();
    printf("now1");

    BBFE_sys_monowrap_solve(
            &(sys.monolis),
            &(sys.monolis_com),
            sys.vals.T,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_DIAG,
            10,
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

void solver_fom_collect_snapmat(
    FE_SYSTEM* sys,
    double t,
    const int step,
    const int iparam)
{
    int rank = monolis_mpi_get_global_my_rank();

#define TRACE(msg) \
    do { \
        fprintf(stderr, \
            "[rank %d param %d step %d] %s\n", \
            rank, iparam, step, msg); \
        fflush(stderr); \
    } while(0)

    /*
     * The caller assembles sys->monolis0 once for the current parameter
     * (D_out,D_in) before entering the time loop.  Coefficients are time
     * independent inside one parameter case, so reuse that matrix at every
     * time step and update only the RHS / Dirichlet contribution.
     */
    monolis_copy_mat_R(
        &(sys->monolis0),
        &(sys->monolis));
    monolis_clear_mat_value_rhs_R(
        &(sys->monolis));

    set_element_vec(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals),
        t);

    manusol_set_theo_sol(
        &(sys->fe),
        sys->vals.theo_sol,
        t);

    BBFE_manusol_set_bc_scalar(
        &(sys->fe),
        &(sys->bc),
        sys->vals.theo_sol,
        t);

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        BLOCK_SIZE,
        &(sys->bc),
        sys->monolis.mat.R.B);

    fprintf(
        stderr,
        "[rank %d param %d step %d] "
        "matrix N=%d NP=%d, internal=%d, total_nodes=%d\n",
        rank,
        iparam,
        step,
        sys->monolis.mat.N,
        sys->monolis.mat.NP,
        sys->monolis_com.n_internal_vertex,
        sys->fe.total_num_nodes);

    fflush(stderr);

    BBFE_sys_monowrap_solve(
        &(sys->monolis),
        &(sys->monolis_com),
        sys->vals.T,
        MONOLIS_ITER_CG,
        MONOLIS_PREC_DIAG,
        sys->vals.mat_max_iter,
        sys->vals.mat_epsilon);

    /*
     * Snapshot
     */
    if(step % sys->vals.snapshot_interval == 0) {

        TRACE("before snapshot");

        /*
         * 1 parameter あたりの time snapshot 数
         */
        int num_time_steps =
            (int)round(
                sys->vals.finish_time /
                sys->vals.dt);

        int num_snap_per_param =
            num_time_steps /
            sys->vals.snapshot_interval;

        /*
         * time方向の snapshot番号
         *
         * 1, 2, 3, ...
         */
        int time_snap_id =
            step /
            sys->vals.snapshot_interval;

        /*
         * parameter方向も含めた global snapshot番号
         *
         * iparam = 0:
         *   1, 2, ..., Ns
         *
         * iparam = 1:
         *   Ns+1, ..., 2Ns
         *
         * ...
         */
        int global_snap_id =
            iparam * num_snap_per_param
            + time_snap_id;

        fprintf(
            stderr,
            "[rank %d] "
            "param=%d step=%d "
            "time_snap=%d global_snap=%d\n",
            rank,
            iparam,
            step,
            time_snap_id,
            global_snap_id);

        fflush(stderr);

        if(monolis_mpi_get_global_comm_size() == 1) {

            ROM_sys_hlpod_fe_set_snap_mat_para(
                sys->vals.T,
                &(sys->rom.hlpod_mat),
                &(sys->bc),
                &(sys->rom.rom_bc),
                sys->fe.total_num_nodes,
                1,
                global_snap_id);
        }
        else {

            ROM_sys_hlpod_fe_set_snap_mat_para(
                sys->vals.T,
                &(sys->rom.hlpod_mat),
                &(sys->bc),
                &(sys->rom.rom_bc),
                sys->monolis_com.n_internal_vertex,
                1,
                global_snap_id);
        }

        TRACE("after snapshot");
    }

    TRACE("LEAVE solver_fom_collect_snapmat");

#undef TRACE
}