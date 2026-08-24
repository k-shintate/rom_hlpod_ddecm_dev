#include "core_FOM.h"
#include "core_NR.h"
#include "fluid_linalg_compat.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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


/* -------------------------------------------------------------------------- */
/* Runtime Nitsche parameters                                                 */
/*                                                                            */
/* Defaults preserve the previous compile-time settings. main() may override  */
/* them once before the time loop. The values are then read-only during solve. */
/* -------------------------------------------------------------------------- */
static double g_t4_nitsche_gamma_n = 100.0;
static double g_t4_nitsche_dt_penalty_coeff = 0.10;

/* -------------------------------------------------------------------------- */
/* TET4/TRI3 Ahmed-body wall-treatment modes                                  */
/*                                                                            */
/* The numeric values are intentionally public-by-convention so main_FOM can  */
/* pass the selected mode without requiring a header change.                  */
/* -------------------------------------------------------------------------- */
#define T4_WALL_STRONG_NOSLIP              0
#define T4_WALL_NITSCHE_NOSLIP             1
#define T4_WALL_NITSCHE_FREESLIP           2
#define T4_WALL_NITSCHE_SPALDING           3
#define T4_WALL_STRONG_NORMAL_FREESLIP     4
#define T4_WALL_STRONG_NORMAL_SPALDING     5

#define T4_WALL_EXCHANGE_WALL_Q            0
#define T4_WALL_EXCHANGE_OPPOSITE_NODE     1

/* Backward-compatible default: the previous Case-C formulation, but main_FOM
 * now selects the mode explicitly and defaults the Spalding exchange point to
 * the owner-TET opposite node. */
static int    g_t4_wall_mode = T4_WALL_NITSCHE_SPALDING;
static int    g_t4_wall_exchange_mode = T4_WALL_EXCHANGE_OPPOSITE_NODE;
static double g_t4_wall_kappa = 0.41;
static double g_t4_wall_B = 5.2;
static double g_t4_wall_ut_eps = 1.0e-12;
static double g_t4_wall_feature_angle_deg = 30.0;

#ifdef __cplusplus
extern "C"
#endif
void BBFE_fluid_set_tet4_wall_runtime_parameters(
    int wall_mode,
    int exchange_mode,
    double gamma_n,
    double dt_penalty_coeff,
    double kappa,
    double B,
    double ut_eps,
    double feature_angle_deg)
{
    if (wall_mode < T4_WALL_STRONG_NOSLIP ||
        wall_mode > T4_WALL_STRONG_NORMAL_SPALDING) {
        fprintf(stderr, "%s ERROR: invalid TET4 wall mode=%d\n", CODENAME, wall_mode);
        exit(EXIT_FAILURE);
    }

    if (exchange_mode != T4_WALL_EXCHANGE_WALL_Q &&
        exchange_mode != T4_WALL_EXCHANGE_OPPOSITE_NODE) {
        fprintf(stderr, "%s ERROR: invalid wall exchange mode=%d\n", CODENAME, exchange_mode);
        exit(EXIT_FAILURE);
    }

    if (!isfinite(gamma_n) || gamma_n <= 0.0 ||
        !isfinite(dt_penalty_coeff) || dt_penalty_coeff < 0.0 ||
        !isfinite(kappa) || kappa <= 0.0 ||
        !isfinite(B) ||
        !isfinite(ut_eps) || ut_eps <= 0.0 ||
        !isfinite(feature_angle_deg) || feature_angle_deg <= 0.0 ||
        feature_angle_deg >= 90.0) {
        fprintf(
            stderr,
            "%s ERROR: invalid wall parameters: mode=%d exch=%d gamma=%.15e "
            "c_dt=%.15e kappa=%.15e B=%.15e ut_eps=%.15e feature=%.15e\n",
            CODENAME,
            wall_mode,
            exchange_mode,
            gamma_n,
            dt_penalty_coeff,
            kappa,
            B,
            ut_eps,
            feature_angle_deg);
        exit(EXIT_FAILURE);
    }

    g_t4_wall_mode = wall_mode;
    g_t4_wall_exchange_mode = exchange_mode;
    g_t4_nitsche_gamma_n = gamma_n;
    g_t4_nitsche_dt_penalty_coeff = dt_penalty_coeff;
    g_t4_wall_kappa = kappa;
    g_t4_wall_B = B;
    g_t4_wall_ut_eps = ut_eps;
    g_t4_wall_feature_angle_deg = feature_angle_deg;
}

static int t4_wall_mode_uses_nitsche_normal(void)
{
    return (
        g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP ||
        g_t4_wall_mode == T4_WALL_NITSCHE_FREESLIP ||
        g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING);
}

static int t4_wall_mode_uses_spalding(void)
{
    return (
        g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING ||
        g_t4_wall_mode == T4_WALL_STRONG_NORMAL_SPALDING);
}

static int t4_wall_mode_uses_projected_strong_normal(void)
{
    return (
        g_t4_wall_mode == T4_WALL_STRONG_NORMAL_FREESLIP ||
        g_t4_wall_mode == T4_WALL_STRONG_NORMAL_SPALDING);
}

#ifdef __cplusplus
extern "C"
#endif
void BBFE_fluid_set_nitsche_runtime_parameters(
    double gamma_n,
    double dt_penalty_coeff)
{
    if (
        !isfinite(gamma_n) || gamma_n <= 0.0 ||
        !isfinite(dt_penalty_coeff) || dt_penalty_coeff < 0.0
    ) {
        fprintf(
            stderr,
            "%s ERROR: invalid Nitsche parameters: gamma_n=%.15e c_dt=%.15e\n",
            CODENAME,
            gamma_n,
            dt_penalty_coeff);
        exit(EXIT_FAILURE);
    }

    g_t4_nitsche_gamma_n = gamma_n;
    g_t4_nitsche_dt_penalty_coeff = dt_penalty_coeff;
}

static void t4n_get_penalty_coefficients(
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double gamma_n,
    double* beta,
    double* beta_mu,
    double* beta_dt)
{
    const double h_eff = fmax(h_n, 1.0e-12);
    const double dt_eff = fmax(dt, 1.0e-30);

    const double b_mu = gamma_n * mu_eff / h_eff;
    const double b_dt =
        gamma_n
        * g_t4_nitsche_dt_penalty_coeff
        * rho
        * h_eff
        / dt_eff;

    if (beta != NULL)    *beta    = b_mu + b_dt;
    if (beta_mu != NULL) *beta_mu = b_mu;
    if (beta_dt != NULL) *beta_dt = b_dt;
}

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

void initialize_velocity_pressure_ahmedbody(
        BBFE_DATA*     fe,
        double** v,
        double* p,
        const int total_num_nodes)
{
    // Initialize velocity array
    for (int i = 0; i < total_num_nodes; i++) {
        v[i][0] = 60.0;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        for(int i = 0; i < fe->local_num_nodes; i++){
            v[fe->conn[e][i]][0] = 0.0;
        }
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
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

void initialize_velocity_pressure_backstep(
    BBFE_DATA*     fe,
        double**    v,
        double*     p,
        const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        if(fe->x[i][0] > 0.5){
            v[i][0] = 1.0;
        }
        else{
            v[i][0] = 0.0;
        }
    }

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

        vals->y_plus = BB_std_calloc_1d_double(vals->y_plus, total_num_nodes);
        vals->y_plus_count = BB_std_calloc_1d_int(vals->y_plus_count, total_num_nodes);
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
        BB_vtk_write_point_vals_scalar(fp, vals->y_plus, fe->total_num_nodes, "yPlus");

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

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

/*
 * z = z0 の節点に、速度3成分のゼロDirichlet条件を追加する。
 *
 * 想定する自由度配置:
 *
 *   block_id = 0 : u_x
 *   block_id = 1 : u_y
 *   block_id = 2 : u_z
 *   block_id = 3 : pressure
 *
 * pressureにはDirichlet条件を追加しない。
 */
int BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
    BBFE_BC* bc,
    const BBFE_DATA* fe,
    const double z0,
    const double tolerance)
{
    if (
        bc == NULL ||
        fe == NULL ||
        bc->D_bc_exists == NULL ||
        bc->imposed_D_val == NULL
    ) {
        fprintf(
            stderr,
            "%s ERROR: invalid argument in "
            "BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane().\n",
            CODENAME
        );

        return -1;
    }

    if (bc->block_size < 3) {
        fprintf(
            stderr,
            "%s ERROR: block_size=%d is smaller than 3.\n",
            CODENAME,
            bc->block_size
        );

        return -1;
    }

    if (
        tolerance <= 0.0 ||
        !isfinite(tolerance)
    ) {
        fprintf(
            stderr,
            "%s ERROR: invalid z-plane tolerance: %.15e\n",
            CODENAME,
            tolerance
        );

        return -1;
    }

    int num_ground_nodes = 0;
    int num_new_constraints = 0;

    for (
        int node_id = 0;
        node_id < bc->total_num_nodes;
        ++node_id
    ) {
        const double z =
            fe->x[node_id][2];

        if (
            !isfinite(z) ||
            fabs(z - z0) > tolerance
        ) {
            continue;
        }

        ++num_ground_nodes;

        /*
         * 速度3成分だけを固定する。
         * pressure成分には触れない。
         */
        for (int block_id = 0; block_id < 3; ++block_id) {
            const int index =
                bc->block_size * node_id
                + block_id;

            if (!bc->D_bc_exists[index]) {
                ++num_new_constraints;
            }

            bc->D_bc_exists[index] = true;
            bc->imposed_D_val[index] = 0.0;
        }
    }

    /*
     * 既存条件との重複を考慮して、
     * D_bc_existsから実際の拘束自由度数を再計算する。
     */
    int num_D_bcs = 0;

    const int num_dofs =
        bc->total_num_nodes
        * bc->block_size;

    for (int i = 0; i < num_dofs; ++i) {
        if (bc->D_bc_exists[i]) {
            ++num_D_bcs;
        }
    }

    bc->num_D_bcs =
        num_D_bcs;

    printf(
        "%s Added zero-velocity Dirichlet condition on z=%.15e:\n"
        "  tolerance             = %.15e\n"
        "  ground nodes          = %d\n"
        "  new velocity DOFs     = %d\n"
        "  total Dirichlet DOFs  = %d\n",
        CODENAME,
        z0,
        tolerance,
        num_ground_nodes,
        num_new_constraints,
        bc->num_D_bcs
    );

    return 0;
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
                                                &(sys.rom_sups.rom_bc), //要変更
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
//                const double tau = tau_e;
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


int read_connectivity_graph_surf(
    BBFE_DATA* fe,
    const char* filename,
    const char* directory,
    int num_integ_points)
{
    FILE* fp = NULL;
    char fname[BUFFER_SIZE];

    const int rank = monolis_mpi_get_global_my_rank();

    snprintf(
        fname,
        BUFFER_SIZE,
        "parted.0/%s.%d",
        filename,
        rank
    );

    /*
     * initialize first
     */
    fe->total_num_elems = 0;
    fe->local_num_nodes = 0;

    fp = BBFE_sys_read_fopen(
        fp,
        fname,
        directory
    );

    /* ========================================================
     * Number of elements
     * ======================================================== */
    if(fscanf(fp, "%d", &(fe->total_num_elems)) != 1) {

        fprintf(
            stderr,
            "%s rank=%d ERROR: "
            "failed to read element count from %s\n",
            CODENAME,
            rank,
            fname
        );

        fclose(fp);
        return -1;
    }

    printf(
        "%s rank=%d: %s Num. elements = %d\n",
        CODENAME,
        rank,
        filename,
        fe->total_num_elems
    );

    if(fe->total_num_elems < 0) {

        fprintf(
            stderr,
            "%s rank=%d ERROR: "
            "negative element count in %s\n",
            CODENAME,
            rank,
            fname
        );

        fclose(fp);
        return -1;
    }


    /* ========================================================
     * Empty element family on this MPI rank
     *
     * file:
     *
     * 0
     *
     * is valid.
     * ======================================================== */
    if(fe->total_num_elems == 0) {

        fe->local_num_nodes = 0;

        fclose(fp);

        printf(
            "%s rank=%d: %s is empty -> skip\n",
            CODENAME,
            rank,
            filename
        );

        return 0;
    }


    /* ========================================================
     * First element
     *
     * format:
     *
     * elem_id NEN node0 node1 ...
     * ======================================================== */
    int elem_id = -1;
    int nen = 0;

    if(fscanf(
        fp,
        "%d %d",
        &elem_id,
        &nen
    ) != 2) {

        fprintf(
            stderr,
            "%s rank=%d ERROR: "
            "failed to read first element header from %s\n",
            CODENAME,
            rank,
            fname
        );

        fclose(fp);
        return -1;
    }

    if(nen <= 0) {

        fprintf(
            stderr,
            "%s rank=%d ERROR: "
            "invalid NEN=%d in %s\n",
            CODENAME,
            rank,
            nen,
            fname
        );

        fclose(fp);
        return -1;
    }

    /*
     * TET -> 4
     * PRI -> 6
     * TRI -> 3
     *
     * automatically detected
     */
    fe->local_num_nodes = nen;


    /* ========================================================
     * Allocation
     *
     * IMPORTANT:
     * local_num_nodes must be set BEFORE this call.
     * ======================================================== */
    BBFE_sys_memory_allocation_elem(
        fe,
        num_integ_points,
        3
    );


    /* ========================================================
     * First element connectivity
     * ======================================================== */
    for(int i = 0; i < fe->local_num_nodes; i++) {

        if(fscanf(
            fp,
            "%d",
            &(fe->conn[0][i])
        ) != 1) {

            fprintf(
                stderr,
                "%s rank=%d ERROR: "
                "failed to read connectivity "
                "e=0 i=%d from %s\n",
                CODENAME,
                rank,
                i,
                fname
            );

            fclose(fp);
            return -1;
        }
    }


    /* ========================================================
     * Remaining elements
     * ======================================================== */
    for(int e = 1; e < fe->total_num_elems; e++) {

        int current_id = -1;
        int current_nen = 0;

        if(fscanf(
            fp,
            "%d %d",
            &current_id,
            &current_nen
        ) != 2) {

            fprintf(
                stderr,
                "%s rank=%d ERROR: "
                "failed element header e=%d in %s\n",
                CODENAME,
                rank,
                e,
                fname
            );

            fclose(fp);
            return -1;
        }

        /*
         * One graph_elem file must contain
         * only one element type.
         */
        if(current_nen != fe->local_num_nodes) {

            fprintf(
                stderr,
                "%s rank=%d ERROR: "
                "different NEN in one file\n"
                "  file      = %s\n"
                "  element   = %d\n"
                "  first NEN = %d\n"
                "  this NEN  = %d\n",
                CODENAME,
                rank,
                fname,
                e,
                fe->local_num_nodes,
                current_nen
            );

            fclose(fp);
            return -1;
        }

        for(int i = 0; i < fe->local_num_nodes; i++) {

            if(fscanf(
                fp,
                "%d",
                &(fe->conn[e][i])
            ) != 1) {

                fprintf(
                    stderr,
                    "%s rank=%d ERROR: "
                    "connectivity read failure "
                    "e=%d i=%d\n",
                    CODENAME,
                    rank,
                    e,
                    i
                );

                fclose(fp);
                return -1;
            }
        }
    }

    fclose(fp);

    printf(
        "%s rank=%d: loaded %s "
        "elements=%d NEN=%d\n",
        CODENAME,
        rank,
        filename,
        fe->total_num_elems,
        fe->local_num_nodes
    );

    return 0;
}

void read_connectivity_graph_surf_internal(
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
    //graph_ndof = fe->total_num_elems;
    fe->total_num_elems = graph_ndof;
    printf("%s Num. elements: %d\n", CODENAME, graph_ndof);
    //fe->local_num_nodes = 4; // 1st-order rectangle
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

void pre_surface_internal(
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

        read_connectivity_graph_surf_internal(
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



void solver_fom_NR2(
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
int max_iter_NR = 20;

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
        //if(i==0){
            //for(int i=0; i<sys.fe.total_num_nodes*4; ++i) rvec[i] = sys.monolis.mat.R.B[i];
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
if (nrm_r / nrm_r_old < 1.0e-6) {
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


}


static int inv3x3(
    const double A[3][3],
    double invA[3][3]);

static double det3x3_local(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int hex8_grad_phys_metric(
    const double xvol[8][3],
    const double dN_dxi[8][3],
    double dN_dx[8][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double J[3][3] = {{0.0}};

    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = det3x3_local(J);
    if (fabs(detJ) < 1.0e-300) return 0;

    if (!inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);

    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void hex8_shape_grad_ref(
    const int SGN[8][3],
    const double xi[3],
    double N[8],
    double dN_dxi[8][3])
{
    for (int a = 0; a < 8; ++a) {
        const double sx = (double)SGN[a][0];
        const double sy = (double)SGN[a][1];
        const double sz = (double)SGN[a][2];

        const double X = 1.0 + sx*xi[0];
        const double Y = 1.0 + sy*xi[1];
        const double Z = 1.0 + sz*xi[2];

        N[a] = 0.125 * X * Y * Z;
        dN_dxi[a][0] = 0.125 * sx * Y * Z;
        dN_dxi[a][1] = 0.125 * X * sy * Z;
        dN_dxi[a][2] = 0.125 * X * Y * sz;
    }
}

static int inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0/det;
    /* invA = adj(A)/det = C^T/det */
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static int hex8_grad_phys(
    const double xvol[8][3],
    const double dN_dxi[8][3],
    double dN_dx[8][3])
{
    double J[3][3] = {{0.0}};

    /* J[i][j] = dx_i / dxi_j */
    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    double invJ[3][3];
    if (!inv3x3(J, invJ)) return 0;

    /* dN/dx_i = sum_j dN/dxi_j * dxi_j/dx_i = dN_dxi[j] * invJ[j][i] */
    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * invJ[j][i];
            }
        }
    }
    return 1;
}

static double eps_nn_from_grad_u(
    const double grad_u[3][3], /* grad_u[a][b] = d u_a / d x_b */
    const double n[3])
{
    double val = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double eps_ab = 0.5*(grad_u[a][b] + grad_u[b][a]);
            val += n[a] * eps_ab * n[b];
        }
    }
    return val;
}


static double spalding_residual_uplus(
    double u_plus,
    double Re_y,
    double kappa,
    double A)
{
    const double x = kappa * u_plus;

    double y_plus;

    /*
      Avoid overflow of exp(x).
      x = 50 is already very large for this purpose.
    */
    double DBL_MAX = 10e13;
    if (x > 50.0) {
        return DBL_MAX * 0.25;
    }

    y_plus =
        u_plus
        + A * (
            exp(x)
            - 1.0
            - x
            - 0.5 * x * x
            - (x * x * x) / 6.0
        );

    return u_plus * y_plus - Re_y;
}


static double compute_utau_all_yplus(
    double ut_norm,
    double y_wall,
    double rho,
    double mu_eff)
{
    const double eps = 1.0e-30;

    if (ut_norm <= eps) return 0.0;
    if (y_wall  <= eps) return 0.0;
    if (rho     <= eps) return 0.0;
    if (mu_eff  <= eps) return 0.0;

    const double nu = mu_eff / rho;
    if (nu <= eps) return 0.0;

    /*
      Re_y = U y / nu
           = u+ y+
    */
    const double Re_y = ut_norm * y_wall / nu;
    if (Re_y <= eps) return 0.0;

    const double kappa = g_t4_wall_kappa;
    const double B     = g_t4_wall_B;
    const double A     = exp(-kappa * B);

    /*
      Solve:
        F(u+) = u+ y+_Spalding(u+) - Re_y = 0
    */
    double lo = 1.0e-12;
    double hi = 1.0;

    while (
        spalding_residual_uplus(hi, Re_y, kappa, A) < 0.0
        && hi < 1.0e6
    ) {
        hi *= 2.0;
    }

    for (int iter = 0; iter < 80; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double fmid =
            spalding_residual_uplus(mid, Re_y, kappa, A);

        if (fmid > 0.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double u_plus = 0.5 * (lo + hi);

    if (u_plus <= eps) return 0.0;

    return ut_norm / u_plus;
}



static void wall_sym_nitsche_local_vecmat(
    double       vec_i[4],
    double       mat_ij[4][4],
    double       Ni,
    double       Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    double       p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double       rho,
    double       mu_eff,       /* VMS 有効粘性: 法線 Nitsche 用 */
    double       mu_wall_law,  /* 分子粘性: all-y+ wall law 用 */
    double       h_n,
    double       y_wall,
    double       dt,
    double       gamma_n)
{
    for (int a = 0; a < 4; ++a) {
        vec_i[a] = 0.0;
        for (int b = 0; b < 4; ++b) {
            mat_ij[a][b] = 0.0;
        }
    }

    const double eps = 1.0e-30;

    const double dt_eff = fmax(dt, eps);
    const double h_eff  = fmax(h_n, 1.0e-12);
    const double y_eff  = fmax(y_wall, 1.0e-12);

    double ur[3];
    for (int a = 0; a < 3; ++a) {
        ur[a] = u_q[a] - Uw[a];
    }

    const double un = dot3(ur, n);

    double ut[3];
    for (int a = 0; a < 3; ++a) {
        ut[a] = ur[a] - un * n[a];
    }

    const double ut_norm = sqrt(dot3(ut, ut));

    /*
      法線方向: no-penetration を symmetric Nitsche で課す。
      接線方向: penalty ではなく all-y+ wall law にする。
    */
    const double beta_base =
        mu_eff / h_eff + rho * h_eff / dt_eff;

    const double beta_n = gamma_n * beta_base;

    /*
      all-y+ wall law:
        tau_w = rho * u_tau^2
        tau_vec = beta_wall * u_t
        beta_wall = rho * u_tau^2 / |u_t|

      |u_t| -> 0 では Spalding の粘性底層極限として
        beta_wall ≈ mu_wall_law / y_wall
      を使う。
    */
    double beta_wall = mu_wall_law / y_eff;

    if (ut_norm > 1.0e-14) {
        const double u_tau =
            compute_utau_all_yplus(
                ut_norm,
                y_eff,
                rho,
                mu_wall_law);

        if (u_tau > 0.0) {
            beta_wall = rho * u_tau * u_tau / ut_norm;
        }
    }

    const double gni = dot3(gradNi, n);
    const double gnj = dot3(gradNj, n);

    const double epsnn_u = eps_nn_from_grad_u(grad_u_q, n);
    const double tnn_u   = -p_q + 2.0 * mu_eff * epsnn_u;

    /*
      Residual:
        法線: symmetric Nitsche
        接線: all-y+ wall law Robin
    */
    for (int a = 0; a < 3; ++a) {
        /* normal consistency */
        vec_i[a] += dt * (-Ni * n[a] * tnn_u);

        /* normal adjoint consistency */
        vec_i[a] += dt * (-(2.0 * mu_eff * n[a] * gni) * un);

        /* normal penalty */
        vec_i[a] += dt * (Ni * beta_n * un * n[a]);

        /* tangential all-y+ wall law */
        vec_i[a] += dt * (Ni * beta_wall * ut[a]);
    }

    /*
      pressure test q = Ni:
      no-penetration Nitsche の pressure consistency
    */
    vec_i[3] += dt * Ni * un;

    /*
      Tangent matrix.
      まずは beta_wall を frozen とする Picard/Newton 混合。
      max_iter_NR = 1 の現状ではこの方が安定。
    */
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double Pt_ab =
                (a == b ? 1.0 : 0.0) - n[a] * n[b];

            /* d/d u_b of normal consistency term */
            mat_ij[a][b] +=
                dt * (-Ni * n[a] *
                      (2.0 * mu_eff * n[b] * gnj));

            /* d/d u_b of normal adjoint-consistency term */
            mat_ij[a][b] +=
                dt * (-(2.0 * mu_eff * n[a] * gni) *
                      Nj * n[b]);

            /* normal penalty */
            mat_ij[a][b] +=
                dt * (Ni * Nj * beta_n * n[a] * n[b]);

            /* tangential all-y+ wall law, frozen beta_wall */
            mat_ij[a][b] +=
                dt * (Ni * Nj * beta_wall * Pt_ab);
        }

        /* d/dp of normal consistency term */
        mat_ij[a][3] +=
            dt * (Ni * Nj * n[a]);
    }

    for (int b = 0; b < 3; ++b) {
        mat_ij[3][b] +=
            dt * (Ni * Nj * n[b]);
    }
}

void set_wall_face_vecmat_symmetric_nitsche(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    /*
     * y+ visualization arrays reset.
     * この関数が呼ばれるたびに、最新の壁面 y+ を作り直す。
     */
    if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            vals->y_plus[i] = 0.0;
            vals->y_plus_count[i] = 0;
        }
    }

    int faces_map[6][4];
    build_hex8_faces_from_sgn(SGN, faces_map);

    int nVolF = 0, nSurfF = 0;
    dl_FaceEntry* volF = dl_build_vol_faces_sorted(fe, faces_map, &nVolF);
    dl_SurfFaceEntry* sF = dl_build_surf_faces_sorted(surf, &nSurfF);

    if (volF == NULL || sF == NULL) {
        fprintf(stderr, "[wall] failed to build face lists\n");
        free(volF);
        free(sF);
        return;
    }

    int* owner = (int*)malloc(sizeof(int) * nSurfF);
    int* lface = (int*)malloc(sizeof(int) * nSurfF);

    if (owner == NULL || lface == NULL) {
        fprintf(stderr, "malloc failed in set_wall_face_vecmat_symmetric_nitsche\n");
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        exit(EXIT_FAILURE);
    }

    for (int e = 0; e < nSurfF; ++e) {
        owner[e] = -1;
        lface[e] = -1;
    }

    int matched = dl_build_owner_map_mergejoin(
        volF,
        nVolF,
        sF,
        nSurfF,
        owner,
        lface);

    if (matched != nSurfF) {
        fprintf(stderr, "[wall] surface-owner mapping failed\n");
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return;
    }

    const int ns = surf->local_num_nodes;      /* normally 4 */
    const int np = basis_surf->num_integ_points;

    if (ns != 4) {
        fprintf(stderr,
            "[wall] set_wall_face_vecmat_symmetric_nitsche assumes QUAD4 surface, ns=%d\n",
            ns);
        free(volF);
        free(sF);
        free(owner);
        free(lface);
        return;
    }

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];

        if (ke < 0) {
            fprintf(stderr,
                "[wall][BUG] owner not found: es=%d owner=%d\n",
                es, ke);
            exit(EXIT_FAILURE);
        }

        {
            int common = 0;

            for (int a = 0; a < ns; ++a) {
                const int sgid = surf->conn[es][a];

                for (int b = 0; b < 8; ++b) {
                    if (fe->conn[ke][b] == sgid) {
                        ++common;
                        break;
                    }
                }
            }

            if (common != ns) {
                fprintf(stderr,
                    "[wall][BUG] surface-owner mismatch: es=%d ke=%d common=%d/%d\n",
                    es, ke, common, ns);

                fprintf(stderr, "  surf nodes:");
                for (int a = 0; a < ns; ++a) {
                    fprintf(stderr, " %d", surf->conn[es][a]);
                }

                fprintf(stderr, "\n  vol nodes:");
                for (int b = 0; b < 8; ++b) {
                    fprintf(stderr, " %d", fe->conn[ke][b]);
                }
                fprintf(stderr, "\n");

                exit(EXIT_FAILURE);
            }
        }

        double xsurf[4][3];
        for (int a = 0; a < ns; ++a) {
            int gid = surf->conn[es][a];
            xsurf[a][0] = fe->x[gid][0];
            xsurf[a][1] = fe->x[gid][1];
            xsurf[a][2] = fe->x[gid][2];
        }

        double xvol[8][3];
        for (int a = 0; a < 8; ++a) {
            int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        double vol_e = 0.0;
        for (int pp = 0; pp < basis_vol->num_integ_points; ++pp) {
            vol_e += basis_vol->integ_weight[pp] * fe->geo[ke][pp].Jacobian;
        }

        double h_e_vms = cbrt(fabs(vol_e));
        if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;

        int s2v[4];
        int ok_s2v = 1;

        for (int a = 0; a < ns; ++a) {
            s2v[a] = -1;

            int sgid = surf->conn[es][a];

            for (int b = 0; b < 8; ++b) {
                if (fe->conn[ke][b] == sgid) {
                    s2v[a] = b;
                    break;
                }
            }

            if (s2v[a] < 0) {
                fprintf(stderr, "[wall] cannot map surface node to volume node\n");
                ok_s2v = 0;
                break;
            }
        }

        if (!ok_s2v) continue;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 8; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 8.0;
        xc[1] /= 8.0;
        xc[2] /= 8.0;

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < ns; ++a) {
            xf[0] += xsurf[a][0];
            xf[1] += xsurf[a][1];
            xf[2] += xsurf[a][2];
        }
        xf[0] /= (double)ns;
        xf[1] /= (double)ns;
        xf[2] /= (double)ns;

        double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < np; ++p) {
            double nrm[3];

            dl_quad4_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nrm);

            const double Jface = sqrt(dot3(nrm, nrm));
            if (Jface <= 1.0e-300) continue;

            double n[3] = {
                nrm[0] / Jface,
                nrm[1] / Jface,
                nrm[2] / Jface
            };

            if (dot3(n, svec) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double xi3[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < ns; ++a) {
                int lv = s2v[a];
                double Na = basis_surf->N[p][a];

                xi3[0] += Na * (double)SGN[lv][0];
                xi3[1] += Na * (double)SGN[lv][1];
                xi3[2] += Na * (double)SGN[lv][2];
            }

            double Nv[8];
            double dN_dxi[8][3];
            double dN_dx[8][3];

            double J_inv_face[3][3];
            double Jacobian_face = 0.0;

            hex8_shape_grad_ref(SGN, xi3, Nv, dN_dxi);

            if (!hex8_grad_phys_metric(
                    xvol,
                    dN_dxi,
                    dN_dx,
                    J_inv_face,
                    &Jacobian_face)) {
                continue;
            }

            double u_q[3]     = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double p_q        = 0.0;

            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3]    = {0.0, 0.0, 0.0};

            for (int a = 0; a < 8; ++a) {
                int gid = fe->conn[ke][a];

                for (int c = 0; c < 3; ++c) {
                    u_q[c]     += Nv[a] * vals->v[gid][c];
                    u_old_q[c] += Nv[a] * vals->v_old[gid][c];

                    for (int d = 0; d < 3; ++d) {
                        grad_u_q[c][d] += dN_dx[a][d] * vals->v[gid][c];
                    }
                }

                p_q += Nv[a] * vals->p[gid];

                for (int d = 0; d < 3; ++d) {
                    grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                }
            }

            double* grad_u_ptr[3] = {
                grad_u_q[0],
                grad_u_q[1],
                grad_u_q[2]
            };

            double mu_eff_face = mu_molecular;
            double tau_face    = 0.0;
            double tau_c_face  = 0.0;

            BBFE_vms_mu_eff_tau(
                &mu_eff_face,
                &tau_face,
                &tau_c_face,
                J_inv_face,
                Jacobian_face,
                h_e_vms,
                u_q,
                u_old_q,
                grad_u_ptr,
                grad_p_q,
                rho,
                mu_molecular,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);

            double h_n =
                fabs((xf[0] - xc[0]) * n[0]
                   + (xf[1] - xc[1]) * n[1]
                   + (xf[2] - xc[2]) * n[2]);

            if (h_n <= 1.0e-12) h_n = 1.0e-12;

            const double y_wall = h_n;
            const double w = basis_surf->integ_weight[p] * Jface;

            /*
             * ------------------------------------------------------------
             * y+ visualization value at this wall quadrature point
             * ------------------------------------------------------------
             */
            if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
                double ur[3];
                for (int a = 0; a < 3; ++a) {
                    ur[a] = u_q[a] - Uw[a];
                }

                const double un = dot3(ur, n);

                double ut[3];
                for (int a = 0; a < 3; ++a) {
                    ut[a] = ur[a] - un * n[a];
                }

                const double ut_norm = sqrt(dot3(ut, ut));

                const double u_tau =
                    compute_utau_all_yplus(
                        ut_norm,
                        y_wall,
                        rho,
                        mu_molecular);

                double yplus_q = 0.0;

                if (u_tau > 0.0 && rho > 0.0 && mu_molecular > 0.0) {
                    yplus_q = rho * u_tau * y_wall / mu_molecular;
                }

                /*
                 * 面上の4節点へ蓄積する。
                 * 求積点重み付き平均にしたい場合は、
                 * y_plus_weight 配列を別途持たせて w で平均する。
                 */
                for (int a = 0; a < ns; ++a) {
                    int gid = surf->conn[es][a];

                    vals->y_plus[gid] += yplus_q;
                    vals->y_plus_count[gid] += 1;
                }
            }

            /*
             * Original Nitsche residual and matrix assembly
             */
            for (int i = 0; i < 8; ++i) {
                int gi = fe->conn[ke][i];

                double vec_i[4];
                double dummy_mat[4][4];
                double zero_grad[3] = {0.0, 0.0, 0.0};

                wall_sym_nitsche_local_vecmat(
                    vec_i,
                    dummy_mat,
                    Nv[i],
                    0.0,
                    dN_dx[i],
                    zero_grad,
                    n,
                    u_q,
                    p_q,
                    grad_u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu_molecular,
                    h_n,
                    y_wall,
                    vals->dt,
                    gamma_n);

                for (int a = 0; a < 4; ++a) {
                    monolis->mat.R.B[4 * gi + a] -= w * vec_i[a];
                }

                for (int j = 0; j < 8; ++j) {
                    int gj = fe->conn[ke][j];

                    double dummy_vec[4];
                    double mat_ij[4][4];

                    wall_sym_nitsche_local_vecmat(
                        dummy_vec,
                        mat_ij,
                        Nv[i],
                        Nv[j],
                        dN_dx[i],
                        dN_dx[j],
                        n,
                        u_q,
                        p_q,
                        grad_u_q,
                        Uw,
                        rho,
                        mu_eff_face,
                        mu_molecular,
                        h_n,
                        y_wall,
                        vals->dt,
                        gamma_n);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                gi,
                                gj,
                                a,
                                b,
                                w * mat_ij[a][b]);
                        }
                    }
                }
            }
        }
    }

    /*
     * Final average at wall nodes.
     * Non-wall nodes are set to 0.0.
     */
    if (vals->y_plus != NULL && vals->y_plus_count != NULL) {
        for (int i = 0; i < fe->total_num_nodes; ++i) {
            if (vals->y_plus_count[i] > 0) {
                vals->y_plus[i] /= (double)vals->y_plus_count[i];
            } else {
                vals->y_plus[i] = 0.0;
            }
        }
    }

    free(volF);
    free(sF);
    free(owner);
    free(lface);
}

/* ---- TET4/TRI3 face-key utilities ---- */
typedef long long t4_node_id_t;

static inline void t4_sort3(t4_node_id_t k[3])
{
    if (k[1] < k[0]) { t4_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
    if (k[2] < k[1]) { t4_node_id_t t = k[1]; k[1] = k[2]; k[2] = t; }
    if (k[1] < k[0]) { t4_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
}


static inline double t4_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double t4_norm3(const double a[3])
{
    return sqrt(t4_dot3(a, a));
}

static inline void t4_normalize3(double a[3])
{
    const double n = t4_norm3(a);
    if (n > 0.0) {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}

static inline void t4_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline int t4_cmp3(const t4_node_id_t a[3], const t4_node_id_t b[3])
{
    if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : +1;
    if (a[1] != b[1]) return (a[1] < b[1]) ? -1 : +1;
    if (a[2] != b[2]) return (a[2] < b[2]) ? -1 : +1;
    return 0;
}

/* Standard outward-oriented faces for a positive-orientation TET4.
 * The orientation is not required for the owner map because keys are sorted;
 * normals are still corrected by centroid direction in the integral loop.
 */
static const int T4_FACE[4][3] = {
    {1, 2, 3},  /* opposite node 0 */
    {0, 3, 2},  /* opposite node 1 */
    {0, 1, 3},  /* opposite node 2 */
    {0, 2, 1}   /* opposite node 3 */
};

typedef struct {
    t4_node_id_t key[3];
    int owner_elem;
    int owner_lface;
} t4_FaceEntry;

typedef struct {
    t4_node_id_t key[3];
    int surf_face_id;
} t4_SurfFaceEntry;

static int t4_cmp_face(const void* A, const void* B)
{
    const t4_FaceEntry* a = (const t4_FaceEntry*)A;
    const t4_FaceEntry* b = (const t4_FaceEntry*)B;
    int c = t4_cmp3(a->key, b->key);
    if (c) return c;
    if (a->owner_elem  != b->owner_elem)  return (a->owner_elem  < b->owner_elem)  ? -1 : +1;
    if (a->owner_lface != b->owner_lface) return (a->owner_lface < b->owner_lface) ? -1 : +1;
    return 0;
}

static int t4_cmp_sface(const void* A, const void* B)
{
    const t4_SurfFaceEntry* a = (const t4_SurfFaceEntry*)A;
    const t4_SurfFaceEntry* b = (const t4_SurfFaceEntry*)B;
    int c = t4_cmp3(a->key, b->key);
    if (c) return c;
    if (a->surf_face_id != b->surf_face_id) return (a->surf_face_id < b->surf_face_id) ? -1 : +1;
    return 0;
}

static t4_FaceEntry* t4_build_vol_faces_sorted(const BBFE_DATA* fe, int* nfaces_out)
{
    const int nfaces = fe->total_num_elems * 4;
    t4_FaceEntry* arr = (t4_FaceEntry*)malloc(sizeof(t4_FaceEntry) * nfaces);
    if (!arr) {
        fprintf(stderr, "[t4][ERR] malloc vol faces failed\n");
        *nfaces_out = 0;
        return NULL;
    }

    int k = 0;
    for (int e = 0; e < fe->total_num_elems; ++e) {
        for (int f = 0; f < 4; ++f) {
            t4_node_id_t key[3] = {
                (t4_node_id_t)fe->conn[e][T4_FACE[f][0]],
                (t4_node_id_t)fe->conn[e][T4_FACE[f][1]],
                (t4_node_id_t)fe->conn[e][T4_FACE[f][2]]
            };
            t4_sort3(key);
            for (int j = 0; j < 3; ++j) arr[k].key[j] = key[j];
            arr[k].owner_elem  = e;
            arr[k].owner_lface = f;
            ++k;
        }
    }

    qsort(arr, nfaces, sizeof(t4_FaceEntry), t4_cmp_face);
    *nfaces_out = nfaces;
    return arr;
}

static t4_SurfFaceEntry* t4_build_surf_faces_sorted(const BBFE_DATA* surf, int* nsurf_out)
{
    const int nsurf = surf->total_num_elems;
    t4_SurfFaceEntry* arr = (t4_SurfFaceEntry*)malloc(sizeof(t4_SurfFaceEntry) * nsurf);
    if (!arr) {
        fprintf(stderr, "[t4][ERR] malloc surf faces failed\n");
        *nsurf_out = 0;
        return NULL;
    }

    for (int e = 0; e < nsurf; ++e) {
        t4_node_id_t key[3] = {
            (t4_node_id_t)surf->conn[e][0],
            (t4_node_id_t)surf->conn[e][1],
            (t4_node_id_t)surf->conn[e][2]
        };
        t4_sort3(key);
        for (int j = 0; j < 3; ++j) arr[e].key[j] = key[j];
        arr[e].surf_face_id = e;
    }

    qsort(arr, nsurf, sizeof(t4_SurfFaceEntry), t4_cmp_sface);
    *nsurf_out = nsurf;
    return arr;
}

static int t4_build_owner_map_mergejoin(
    const t4_FaceEntry* volF,
    int nVolF,
    const t4_SurfFaceEntry* surfF,
    int nSurfF,
    int* owner_elem,
    int* owner_lface)
{
    int i = 0;
    int j = 0;
    int matches = 0;

    while (i < nVolF && j < nSurfF) {
        int c = t4_cmp3(volF[i].key, surfF[j].key);

        if (c == 0) {
            const int es = surfF[j].surf_face_id;
            owner_elem[es]  = volF[i].owner_elem;
            owner_lface[es] = volF[i].owner_lface;
            ++matches;

            t4_node_id_t key0[3] = { volF[i].key[0], volF[i].key[1], volF[i].key[2] };
            do { ++i; } while (i < nVolF && t4_cmp3(key0, volF[i].key) == 0);
            ++j;
        } else if (c < 0) {
            ++i;
        } else {
            ++j;
        }
    }

    return matches;
}

static int t4_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner_out,
    int** lface_out)
{
    if (owner_out == NULL || lface_out == NULL) {
        fprintf(
            stderr,
            "[t4] owner_out or lface_out is NULL\n");

        return -1;
    }

    *owner_out = NULL;
    *lface_out = NULL;

    if (surf == NULL || fe == NULL) {
        fprintf(
            stderr,
            "[t4] surf or fe is NULL\n");

        return -1;
    }

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[t4] invalid number of surface elements: %d\n",
            surf->total_num_elems);

        return -1;
    }

    /*
     * このrankが壁面要素を所有しない場合は正常。
     *
     * owner_out、lface_outはNULLのまま返す。
     * 呼出側の表面要素ループは0回なので問題ない。
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    /*
     * 非空rankに対してのみ要素形式を検査する。
     */
    if (
        fe->local_num_nodes != 4 ||
        surf->local_num_nodes != 3) {

        fprintf(
            stderr,
            "[t4] expected TET4/TRI3 but got "
            "fe->local_num_nodes=%d "
            "surf->local_num_nodes=%d "
            "surf->total_num_elems=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        return -1;
    }

    int nVolF = 0;
    int nSurfF = 0;

    t4_FaceEntry* volF =
        t4_build_vol_faces_sorted(
            fe,
            &nVolF);

    t4_SurfFaceEntry* sF =
        t4_build_surf_faces_sorted(
            surf,
            &nSurfF);

    if (volF == NULL || sF == NULL) {
        free(volF);
        free(sF);

        return -1;
    }

    /*
     * surf->total_num_elems > 0なのでmalloc(0)にはならない。
     */
    int* owner =
        (int*)malloc(
            sizeof(int) * (size_t)nSurfF);

    int* lface =
        (int*)malloc(
            sizeof(int) * (size_t)nSurfF);

    if (owner == NULL || lface == NULL) {
        free(volF);
        free(sF);
        free(owner);
        free(lface);

        return -1;
    }

    for (int e = 0; e < nSurfF; ++e) {
        owner[e] = -1;
        lface[e] = -1;
    }

    const int matched =
        t4_build_owner_map_mergejoin(
            volF,
            nVolF,
            sF,
            nSurfF,
            owner,
            lface);

    if (matched != nSurfF) {
        int cnt = 0;

        for (int e = 0; e < nSurfF; ++e) {
            if (owner[e] < 0) {
                fprintf(
                    stderr,
                    "[t4] surf tri %d has NO owner\n",
                    e);

                ++cnt;
            }
        }

        fprintf(
            stderr,
            "[t4] %d / %d surface triangles unmapped. "
            "Check node IDs / connectivity.\n",
            cnt,
            nSurfF);

        free(volF);
        free(sF);
        free(owner);
        free(lface);

        return -1;
    }

    free(volF);
    free(sF);

    *owner_out = owner;
    *lface_out = lface;

    return 0;
}


/* ========================================================================== */
/* TET4/TRI3 wall mesh and y+ diagnostics                                     */
/* ========================================================================== */

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef T4_WALL_MESH_DIAG_EPS
#define T4_WALL_MESH_DIAG_EPS 1.0e-30
#endif

#define T4_WMD_PI 3.141592653589793238462643383279502884


/* -------------------------------------------------------------------------- */
/* Output-name table                                                          */
/* -------------------------------------------------------------------------- */

static const char*
t4_wmd_stat_name[T4_WALL_MESH_DIAG_NUM_STATS] = {
    "yplus_grad_centroid",
    "yplus_spalding_centroid",
    "yplus_spalding_over_grad",
    "utau",
    "tau_w",
    "Cf",
    "h_n_penalty",
    "owner_tet_face_height",
    "h_t",
    "aspect_ratio_ht_over_hn",
    "nonorthogonality_deg",
    "owner_tet_quality",
    "mu_eff_over_mu",
    "beta_dt_over_beta_mu",
    "penalty_to_shear",
    "owner_CFL",
    "owner_Re_h",
    "owner_Pe_eff"
};


/* -------------------------------------------------------------------------- */
/* Small vector helpers                                                       */
/* -------------------------------------------------------------------------- */

static double t4_wmd_dot3(
    const double a[3],
    const double b[3])
{
    return
        a[0]*b[0]
      + a[1]*b[1]
      + a[2]*b[2];
}


static double t4_wmd_norm3(
    const double a[3])
{
    return sqrt(fmax(t4_wmd_dot3(a, a), 0.0));
}


static void t4_wmd_cross3(
    double c[3],
    const double a[3],
    const double b[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}


static void t4_wmd_normalize3(
    double a[3])
{
    const double norm_a = t4_wmd_norm3(a);

    if (norm_a > T4_WALL_MESH_DIAG_EPS) {
        a[0] /= norm_a;
        a[1] /= norm_a;
        a[2] /= norm_a;
    }
}


/* -------------------------------------------------------------------------- */
/* 3x3 determinant and inverse                                                */
/* -------------------------------------------------------------------------- */

static double t4_wmd_det3(
    const double A[3][3])
{
    return
        A[0][0] * (
            A[1][1]*A[2][2]
          - A[1][2]*A[2][1])
      - A[0][1] * (
            A[1][0]*A[2][2]
          - A[1][2]*A[2][0])
      + A[0][2] * (
            A[1][0]*A[2][1]
          - A[1][1]*A[2][0]);
}


static int t4_wmd_inv3(
    const double A[3][3],
    double invA[3][3],
    double* det_out)
{
    const double detA = t4_wmd_det3(A);

    if (det_out != NULL) {
        *det_out = detA;
    }

    if (fabs(detA) <= 1.0e-300) {
        memset(invA, 0, sizeof(double)*9);
        return 0;
    }

    const double id = 1.0 / detA;

    invA[0][0] =
        (A[1][1]*A[2][2] - A[1][2]*A[2][1]) * id;

    invA[0][1] =
        (A[0][2]*A[2][1] - A[0][1]*A[2][2]) * id;

    invA[0][2] =
        (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * id;

    invA[1][0] =
        (A[1][2]*A[2][0] - A[1][0]*A[2][2]) * id;

    invA[1][1] =
        (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * id;

    invA[1][2] =
        (A[0][2]*A[1][0] - A[0][0]*A[1][2]) * id;

    invA[2][0] =
        (A[1][0]*A[2][1] - A[1][1]*A[2][0]) * id;

    invA[2][1] =
        (A[0][1]*A[2][0] - A[0][0]*A[2][1]) * id;

    invA[2][2] =
        (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * id;

    return 1;
}


/* -------------------------------------------------------------------------- */
/* Exact TET4 geometry                                                        */
/*                                                                            */
/* J[:,0] = x1-x0                                                             */
/* J[:,1] = x2-x0                                                             */
/* J[:,2] = x3-x0                                                             */
/* -------------------------------------------------------------------------- */

static int t4_wmd_tet4_geometry(
    const BBFE_DATA* fe,
    int ke,
    double J[3][3],
    double invJ[3][3],
    double dN_dx[4][3],
    double* detJ_out,
    double* volume_out)
{
    double x[4][3];

    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];

        x[a][0] = fe->x[gid][0];
        x[a][1] = fe->x[gid][1];
        x[a][2] = fe->x[gid][2];
    }

    for (int i = 0; i < 3; ++i) {
        J[i][0] = x[1][i] - x[0][i];
        J[i][1] = x[2][i] - x[0][i];
        J[i][2] = x[3][i] - x[0][i];
    }

    double detJ = 0.0;

    if (!t4_wmd_inv3(J, invJ, &detJ)) {
        if (detJ_out != NULL) {
            *detJ_out = 0.0;
        }

        if (volume_out != NULL) {
            *volume_out = 0.0;
        }

        memset(dN_dx, 0, sizeof(double)*12);
        return 0;
    }

    static const double dN_dxi[4][3] = {
        {-1.0, -1.0, -1.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0}
    };

    /*
     * dN/dx_i = sum_j dN/dxi_j * dxi_j/dx_i
     */
    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;

            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] +=
                    dN_dxi[a][j] * invJ[j][i];
            }
        }
    }

    if (detJ_out != NULL) {
        *detJ_out = detJ;
    }

    if (volume_out != NULL) {
        *volume_out = fabs(detJ) / 6.0;
    }

    return 1;
}


/* -------------------------------------------------------------------------- */
/* Surface nodes to owner TET local nodes                                     */
/* -------------------------------------------------------------------------- */

static int t4_wmd_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3],
    int* opposite_local_node)
{
    int used[4] = {0, 0, 0, 0};

    for (int a = 0; a < 3; ++a) {
        const int surface_gid = surf->conn[es][a];

        s2v[a] = -1;

        for (int b = 0; b < 4; ++b) {
            if (fe->conn[ke][b] == surface_gid) {
                s2v[a] = b;
                used[b] = 1;
                break;
            }
        }

        if (s2v[a] < 0) {
            return 0;
        }
    }

    if (opposite_local_node != NULL) {
        *opposite_local_node = -1;

        for (int b = 0; b < 4; ++b) {
            if (!used[b]) {
                *opposite_local_node = b;
                break;
            }
        }

        if (*opposite_local_node < 0) {
            return 0;
        }
    }

    return 1;
}


/* -------------------------------------------------------------------------- */
/* TET4 velocity and pressure gradients                                       */
/* -------------------------------------------------------------------------- */

static void t4_wmd_tet4_gradients(
    const BBFE_DATA* fe,
    const VALUES* vals,
    int ke,
    const double dN_dx[4][3],
    double grad_u[3][3],
    double grad_p[3])
{
    memset(grad_u, 0, sizeof(double)*9);

    grad_p[0] = 0.0;
    grad_p[1] = 0.0;
    grad_p[2] = 0.0;

    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];

        for (int c = 0; c < 3; ++c) {
            for (int d = 0; d < 3; ++d) {
                grad_u[c][d] +=
                    vals->v[gid][c] * dN_dx[a][d];
            }
        }

        for (int d = 0; d < 3; ++d) {
            grad_p[d] +=
                vals->p[gid] * dN_dx[a][d];
        }
    }
}


/* -------------------------------------------------------------------------- */
/* Mean-ratio TET4 quality                                                    */
/*                                                                            */
/* Regular tetrahedron: quality = 1                                           */
/* Degenerate tetrahedron: quality -> 0                                       */
/* -------------------------------------------------------------------------- */

static double t4_wmd_tet4_quality(
    const BBFE_DATA* fe,
    int ke,
    double volume)
{
    static const int edge[6][2] = {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {1, 3},
        {2, 3}
    };

    double sum_l2 = 0.0;

    for (int e = 0; e < 6; ++e) {
        const int g0 = fe->conn[ke][edge[e][0]];
        const int g1 = fe->conn[ke][edge[e][1]];

        const double dx =
            fe->x[g1][0] - fe->x[g0][0];

        const double dy =
            fe->x[g1][1] - fe->x[g0][1];

        const double dz =
            fe->x[g1][2] - fe->x[g0][2];

        sum_l2 +=
            dx*dx + dy*dy + dz*dz;
    }

    if (
        volume <= T4_WALL_MESH_DIAG_EPS ||
        sum_l2 <= T4_WALL_MESH_DIAG_EPS) {

        return 0.0;
    }

    const double quality =
        12.0 * pow(3.0*volume, 2.0/3.0)
        / sum_l2;

    return fmax(0.0, fmin(1.0, quality));
}


/* -------------------------------------------------------------------------- */
/* Spalding-law utility                                                       */
/* -------------------------------------------------------------------------- */

static double t4_wmd_spalding_yplus(
    double u_plus)
{
    const double kappa = 0.41;
    const double B = 5.2;
    const double A = exp(-kappa * B);

    const double x = kappa * u_plus;

    if (x > 50.0) {
        return DBL_MAX * 0.25;
    }

    return
        u_plus
      + A * (
            exp(x)
          - 1.0
          - x
          - 0.5*x*x
          - x*x*x/6.0);
}


static double t4_wmd_spalding_residual(
    double u_plus,
    double Re_y)
{
    const double y_plus =
        t4_wmd_spalding_yplus(u_plus);

    return u_plus*y_plus - Re_y;
}


static double t4_wmd_compute_utau_spalding(
    double ut_norm,
    double y,
    double rho,
    double mu)
{
    if (
        ut_norm <= T4_WALL_MESH_DIAG_EPS ||
        y <= T4_WALL_MESH_DIAG_EPS ||
        rho <= T4_WALL_MESH_DIAG_EPS ||
        mu <= T4_WALL_MESH_DIAG_EPS) {

        return 0.0;
    }

    const double nu = mu / rho;
    const double Re_y = ut_norm*y/nu;

    if (Re_y <= T4_WALL_MESH_DIAG_EPS) {
        return 0.0;
    }

    /*
     * Solve:
     *
     *   u+ y+(u+) = Re_y
     *
     * by bisection.
     */
    double lo = 1.0e-12;
    double hi = fmax(sqrt(Re_y) + 10.0, 10.0);

    double f_lo =
        t4_wmd_spalding_residual(lo, Re_y);

    double f_hi =
        t4_wmd_spalding_residual(hi, Re_y);

    for (int expand = 0;
         expand < 80 && f_hi < 0.0;
         ++expand) {

        hi *= 2.0;

        f_hi =
            t4_wmd_spalding_residual(hi, Re_y);

        if (hi > 1.0e8) {
            break;
        }
    }

    if (!(f_lo <= 0.0 && f_hi >= 0.0)) {
        /*
         * Linear-law fallback:
         *
         * u+ = y+, Re_y = u+^2
         */
        const double u_plus =
            sqrt(fmax(Re_y, 0.0));

        return ut_norm / fmax(u_plus, 1.0e-30);
    }

    for (int it = 0; it < 100; ++it) {
        const double mid = 0.5*(lo + hi);

        const double f_mid =
            t4_wmd_spalding_residual(mid, Re_y);

        if (f_mid > 0.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const double u_plus =
        0.5*(lo + hi);

    return ut_norm / fmax(u_plus, 1.0e-30);
}


/* -------------------------------------------------------------------------- */
/* Diagnostic structure initialization                                       */
/* -------------------------------------------------------------------------- */

static void t4_wmd_init(
    T4WallMeshDiagnostics* diag)
{
    memset(diag, 0, sizeof(*diag));

    for (int k = 0;
         k < T4_WALL_MESH_DIAG_NUM_STATS;
         ++k) {

        diag->min[k] = 1.0e300;
        diag->max[k] = -1.0e300;
    }
}


/* -------------------------------------------------------------------------- */
/* Add one area-weighted scalar sample                                        */
/* -------------------------------------------------------------------------- */

static void t4_wmd_add_stat(
    T4WallMeshDiagnostics* diag,
    int stat_id,
    double value,
    double area)
{
    if (
        stat_id < 0 ||
        stat_id >= T4_WALL_MESH_DIAG_NUM_STATS ||
        !isfinite(value) ||
        !isfinite(area) ||
        area <= 0.0) {

        return;
    }

    diag->sum[stat_id] +=
        value * area;

    diag->sum2[stat_id] +=
        value * value * area;

    diag->min[stat_id] =
        fmin(diag->min[stat_id], value);

    diag->max[stat_id] =
        fmax(diag->max[stat_id], value);
}


/* -------------------------------------------------------------------------- */
/* Finalize statistics                                                        */
/* -------------------------------------------------------------------------- */

static void t4_wmd_finalize(
    T4WallMeshDiagnostics* diag)
{
    if (diag->area <= 0.0) {
        for (int k = 0;
             k < T4_WALL_MESH_DIAG_NUM_STATS;
             ++k) {

            diag->mean[k] = 0.0;
            diag->rms[k] = 0.0;
            diag->min[k] = 0.0;
            diag->max[k] = 0.0;
        }

        for (int b = 0;
             b < T4_WALL_MESH_DIAG_NUM_YPLUS_BINS;
             ++b) {

            diag->yplus_bin_fraction[b] = 0.0;
        }

        diag->reverse_shear_area_fraction = 0.0;
        diag->nonorth_gt_30_fraction = 0.0;
        diag->nonorth_gt_60_fraction = 0.0;
        diag->quality_lt_0p1_fraction = 0.0;
        diag->quality_lt_0p2_fraction = 0.0;
        diag->mu_ratio_gt_2_fraction = 0.0;
        diag->mu_ratio_gt_10_fraction = 0.0;
        diag->beta_dt_dominant_fraction = 0.0;

        return;
    }

    for (int k = 0;
         k < T4_WALL_MESH_DIAG_NUM_STATS;
         ++k) {

        diag->mean[k] =
            diag->sum[k] / diag->area;

        diag->rms[k] =
            sqrt(fmax(
                diag->sum2[k] / diag->area,
                0.0));

        if (
            diag->min[k] > 0.5e300 ||
            diag->max[k] < -0.5e300) {

            diag->min[k] = 0.0;
            diag->max[k] = 0.0;
        }
    }

    for (int b = 0;
         b < T4_WALL_MESH_DIAG_NUM_YPLUS_BINS;
         ++b) {

        diag->yplus_bin_fraction[b] =
            diag->yplus_bin_area[b] / diag->area;
    }

    if (diag->streamwise_shear_valid_area > 0.0) {
        diag->reverse_shear_area_fraction =
            diag->reverse_shear_area
            / diag->streamwise_shear_valid_area;
    } else {
        diag->reverse_shear_area_fraction = 0.0;
    }

    diag->nonorth_gt_30_fraction =
        diag->nonorth_gt_30_area / diag->area;

    diag->nonorth_gt_60_fraction =
        diag->nonorth_gt_60_area / diag->area;

    diag->quality_lt_0p1_fraction =
        diag->quality_lt_0p1_area / diag->area;

    diag->quality_lt_0p2_fraction =
        diag->quality_lt_0p2_area / diag->area;

    diag->mu_ratio_gt_2_fraction =
        diag->mu_ratio_gt_2_area / diag->area;

    diag->mu_ratio_gt_10_fraction =
        diag->mu_ratio_gt_10_area / diag->area;

    diag->beta_dt_dominant_fraction =
        diag->beta_dt_dominant_area / diag->area;
}


static int t4n_make_owner_map_overlap(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner_out,
    int** lface_out,
    int* num_unmapped_out);



int calc_wall_mesh_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    VALUES* vals,
    double rho,
    double mu_molecular,
    double Uref,
    const double Uw_in[3],
    const double eU_in[3],
    double gamma_n,
    int write_nodal_yplus,
    T4WallMeshDiagnostics* diag)
{
    if (
        surf == NULL ||
        fe == NULL ||
        vals == NULL ||
        diag == NULL) {

        fprintf(
            stderr,
            "[wall_mesh_diag] NULL input\n");

        return -1;
    }

    t4_wmd_init(diag);

    /*
     * Reset the local VTK yPlus field on every rank, including ranks that
     * own no wall faces. This prevents values from a previous time step
     * remaining on an empty-surface rank.
     */
    if (
        write_nodal_yplus &&
        vals->y_plus != NULL &&
        vals->y_plus_count != NULL) {

        for (int i = 0;
             i < fe->total_num_nodes;
             ++i) {

            vals->y_plus[i] = 0.0;
            vals->y_plus_count[i] = 0;
        }
    }

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[wall_mesh_diag] "
            "invalid number of surface elements: %d\n",
            surf->total_num_elems);

        return -1;
    }

    /*
     * A rank that owns no wall faces is valid and contributes zero to all
     * SUM reductions. Keep the neutral MIN/MAX values set by t4_wmd_init().
     *
     * Do not call t4_wmd_finalize() here: converting an empty rank's MIN
     * values to zero before the global reduction would corrupt global minima.
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    /*
     * Element-type checks apply only to ranks with at least one wall face.
     */
    if (
        surf->local_num_nodes != 3 ||
        fe->local_num_nodes != 4) {

        fprintf(
            stderr,
            "[wall_mesh_diag] "
            "TET4/TRI3 expected, got fe=%d surf=%d ne_surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        return -1;
    }

    const double rho_eff =
        fmax(rho, 1.0e-30);

    const double mu_eff_safe =
        fmax(mu_molecular, 1.0e-30);

    const double dt_eff =
        fmax(vals->dt, 1.0e-30);

    const double U_scale =
        fmax(fabs(Uref), 1.0e-30);

    const double Uw[3] = {
        Uw_in ? Uw_in[0] : 0.0,
        Uw_in ? Uw_in[1] : 0.0,
        Uw_in ? Uw_in[2] : 0.0
    };

    double eU[3] = {
        eU_in ? eU_in[0] : 1.0,
        eU_in ? eU_in[1] : 0.0,
        eU_in ? eU_in[2] : 0.0
    };

    t4_wmd_normalize3(eU);

    int* owner = NULL;
    int* lface = NULL;
    int num_unmapped_local = 0;

    if (
        t4n_make_owner_map_overlap(
            surf,
            fe,
            &owner,
            &lface,
            &num_unmapped_local) != 0
    ) {
        fprintf(
            stderr,
            "[wall_mesh_diag] "
            "overlap owner-map construction failed\n");

        return -1;
    }

    for (int es = 0;
         es < surf->total_num_elems;
         ++es) {

        const int ke = owner[es];

        if (ke < 0) {
            /*
            * owner TETを持たない冗長overlap surfaceコピー。
            * 無効形状ではない。
            */
            continue;
        }

        int s2v[3];
        int opposite_local_node = -1;

        if (!t4_wmd_surface_to_volume_nodes(
                surf,
                fe,
                es,
                ke,
                s2v,
                &opposite_local_node)) {

            diag->invalid_face_count += 1.0;
            continue;
        }

        /*
         * Surface coordinates.
         */
        double xtri[3][3];

        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];

            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        const double edge01[3] = {
            xtri[1][0] - xtri[0][0],
            xtri[1][1] - xtri[0][1],
            xtri[1][2] - xtri[0][2]
        };

        const double edge02[3] = {
            xtri[2][0] - xtri[0][0],
            xtri[2][1] - xtri[0][1],
            xtri[2][2] - xtri[0][2]
        };

        double face_cross[3];

        t4_wmd_cross3(
            face_cross,
            edge01,
            edge02);

        const double cross_norm =
            t4_wmd_norm3(face_cross);

        const double face_area =
            0.5 * cross_norm;

        if (
            !isfinite(face_area) ||
            face_area <= 1.0e-300) {

            diag->invalid_face_count += 1.0;
            continue;
        }

        double n[3] = {
            face_cross[0] / cross_norm,
            face_cross[1] / cross_norm,
            face_cross[2] / cross_norm
        };

        /*
         * Face centroid.
         */
        const double xf[3] = {
            (xtri[0][0] + xtri[1][0] + xtri[2][0]) / 3.0,
            (xtri[0][1] + xtri[1][1] + xtri[2][1]) / 3.0,
            (xtri[0][2] + xtri[1][2] + xtri[2][2]) / 3.0
        };

        /*
         * Owner TET centroid.
         */
        double xc[3] = {0.0, 0.0, 0.0};

        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];

            xc[0] += 0.25 * fe->x[gid][0];
            xc[1] += 0.25 * fe->x[gid][1];
            xc[2] += 0.25 * fe->x[gid][2];
        }

        /*
         * Direction from TET centroid toward wall face.
         * This should be outward.
         */
        const double centroid_to_face[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        if (
            t4_wmd_dot3(
                centroid_to_face,
                n) < 0.0) {

            n[0] *= -1.0;
            n[1] *= -1.0;
            n[2] *= -1.0;
        }

        /*
         * Exact TET4 geometry and gradients.
         */
        double J[3][3];
        double invJ[3][3];
        double dN_dx[4][3];
        double detJ = 0.0;
        double volume = 0.0;

        if (!t4_wmd_tet4_geometry(
                fe,
                ke,
                J,
                invJ,
                dN_dx,
                &detJ,
                &volume)) {

            diag->invalid_face_count += 1.0;
            continue;
        }

        if (
            volume <= 1.0e-300 ||
            !isfinite(volume)) {

            diag->invalid_face_count += 1.0;
            continue;
        }

        /*
         * Match the current TET4/TRI3 Nitsche implementation:
         *
         *   h_n = 3V/A = |detJ|/Jface_raw
         *
         * This is the full owner-TET altitude normal to the wall face.
         */
        double h_n =
            3.0 * volume
            / fmax(face_area, 1.0e-300);

        if (
            !isfinite(h_n) ||
            h_n <= 1.0e-12
        ) {
            h_n = 1.0e-12;
        }

        /*
         * Full opposite-node to face-plane height.
         * For a TET centroid, this is nominally 4*h_n.
         */
        const int opposite_gid =
            fe->conn[ke][opposite_local_node];

        const double xopp_to_face[3] = {
            fe->x[opposite_gid][0] - xf[0],
            fe->x[opposite_gid][1] - xf[1],
            fe->x[opposite_gid][2] - xf[2]
        };

        const double face_height =
            fabs(t4_wmd_dot3(
                xopp_to_face,
                n));

        /*
         * Equivalent tangential edge scale for an equal-area equilateral
         * triangle.
         */
        const double h_t =
            sqrt(4.0*face_area/sqrt(3.0));

        const double aspect_ratio =
            h_t / fmax(h_n, 1.0e-30);

        /*
         * Nonorthogonality:
         * angle between centroid-to-face vector and face normal.
         */
        const double d_centroid =
            t4_wmd_norm3(centroid_to_face);

        double nonorthogonality_deg = 0.0;

        if (d_centroid > 1.0e-30) {
            double cos_theta =
                fabs(t4_wmd_dot3(
                    centroid_to_face,
                    n))
                / d_centroid;

            cos_theta =
                fmax(0.0, fmin(1.0, cos_theta));

            nonorthogonality_deg =
                acos(cos_theta)
                * 180.0 / T4_WMD_PI;
        }

        const double tet_quality =
            t4_wmd_tet4_quality(
                fe,
                ke,
                volume);

        /*
         * Exact velocity/pressure gradients in owner TET4.
         */
        double grad_u[3][3];
        double grad_p[3];

        t4_wmd_tet4_gradients(
            fe,
            vals,
            ke,
            dN_dx,
            grad_u,
            grad_p);

        /*
         * Owner centroid velocity.
         */
        double u_c[3] = {0.0, 0.0, 0.0};

        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];

            u_c[0] += 0.25 * vals->v[gid][0];
            u_c[1] += 0.25 * vals->v[gid][1];
            u_c[2] += 0.25 * vals->v[gid][2];
        }

        /*
         * Face-centroid current/old velocity and pressure.
         */
        double u_face[3] = {0.0, 0.0, 0.0};
        double u_old_face[3] = {0.0, 0.0, 0.0};
        double p_face = 0.0;

        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];

            u_face[0] += vals->v[gid][0] / 3.0;
            u_face[1] += vals->v[gid][1] / 3.0;
            u_face[2] += vals->v[gid][2] / 3.0;

            p_face += vals->p[gid] / 3.0;

            if (vals->v_old != NULL) {
                u_old_face[0] += vals->v_old[gid][0] / 3.0;
                u_old_face[1] += vals->v_old[gid][1] / 3.0;
                u_old_face[2] += vals->v_old[gid][2] / 3.0;
            } else {
                u_old_face[0] += vals->v[gid][0] / 3.0;
                u_old_face[1] += vals->v[gid][1] / 3.0;
                u_old_face[2] += vals->v[gid][2] / 3.0;
            }
        }

        /*
         * VMS effective viscosity at face centroid.
         */
        double mu_eff_face = mu_molecular;
        double tau_face = 0.0;
        double tau_c_face = 0.0;

        double* grad_u_ptr[3] = {
            grad_u[0],
            grad_u[1],
            grad_u[2]
        };

        const double h_e_vms =
            cbrt(volume);

        BBFE_vms_mu_eff_tau(
            &mu_eff_face,
            &tau_face,
            &tau_c_face,
            invJ,
            fabs(detJ),
            h_e_vms,
            u_face,
            u_old_face,
            grad_u_ptr,
            grad_p,
            rho,
            mu_molecular,
            vals->dt,
            vals->C_vms,
            vals->vms_cap_coeff);

        if (
            !isfinite(mu_eff_face) ||
            mu_eff_face <= 0.0) {

            mu_eff_face = mu_molecular;
        }

        /*
         * Molecular wall shear:
         *
         * tau_w = P_t mu (grad u + grad u^T) n
         */
        double sym_grad_n[3] = {0.0, 0.0, 0.0};

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                sym_grad_n[a] +=
                    (
                        grad_u[a][b]
                      + grad_u[b][a]
                    ) * n[b];
            }
        }

        const double sym_grad_nn =
            t4_wmd_dot3(
                sym_grad_n,
                n);

        double tau_w_vec[3];

        for (int a = 0; a < 3; ++a) {
            tau_w_vec[a] =
                mu_molecular
                * (
                    sym_grad_n[a]
                  - sym_grad_nn*n[a]
                );
        }

        const double tau_w_mag =
            t4_wmd_norm3(tau_w_vec);

        const double u_tau =
            sqrt(
                tau_w_mag
                / rho_eff);

        /*
         * Gradient-based y+ evaluated at owner TET centroid distance h_n.
         */
        const double yplus_grad =
            rho_eff * u_tau * h_n
            / mu_eff_safe;

        /*
         * Spalding y+ based on owner-centroid tangential velocity.
         */
        const double r_c[3] = {
            u_c[0] - Uw[0],
            u_c[1] - Uw[1],
            u_c[2] - Uw[2]
        };

        const double r_c_n =
            t4_wmd_dot3(
                r_c,
                n);

        const double r_c_t[3] = {
            r_c[0] - r_c_n*n[0],
            r_c[1] - r_c_n*n[1],
            r_c[2] - r_c_n*n[2]
        };

        const double ut_centroid =
            t4_wmd_norm3(r_c_t);

        const double u_tau_spalding =
            t4_wmd_compute_utau_spalding(
                ut_centroid,
                h_n,
                rho_eff,
                mu_eff_safe);

        const double yplus_spalding =
            rho_eff * u_tau_spalding * h_n
            / mu_eff_safe;

        double yplus_ratio = 0.0;

        if (yplus_grad > 1.0e-30) {
            yplus_ratio =
                yplus_spalding / yplus_grad;
        }

        const double Cf =
            2.0 * tau_w_mag
            / fmax(
                rho_eff * U_scale * U_scale,
                1.0e-30);

        /*
         * Nitsche penalty diagnostics.
         */
        double beta_mu = 0.0;
        double beta_dt = 0.0;

        t4n_get_penalty_coefficients(
            rho_eff,
            mu_eff_face,
            h_n,
            dt_eff,
            gamma_n,
            NULL,
            &beta_mu,
            &beta_dt);

        const double beta_dt_over_beta_mu =
            beta_dt
            / fmax(beta_mu, 1.0e-30);

        /*
         * Wall-mode diagnostic: ratio of normal-Nitsche penalty stress scale
         * to the modeled Spalding tangential wall-shear magnitude.
         * This is a stiffness/scale diagnostic only; the two vectors are
         * orthogonal and must not be added as scalar stresses.
         */
        const double beta_n_diag = beta_mu + beta_dt;
        const double normal_penalty_stress =
            beta_n_diag * fabs(r_c_n);
        const double walllaw_shear =
            rho_eff * u_tau_spalding * u_tau_spalding;
        const double penalty_to_shear =
            normal_penalty_stress
            / fmax(walllaw_shear, 1.0e-30);

        const double mu_eff_over_mu =
            mu_eff_face / mu_eff_safe;

        /*
         * Local temporal and convective resolution indicators.
         */
        const double speed_c =
            t4_wmd_norm3(u_c);

        const double CFL =
            speed_c * dt_eff / h_n;

        const double Re_h =
            rho_eff * speed_c * h_n
            / mu_eff_safe;

        const double Pe_eff =
            rho_eff * speed_c * h_n
            / fmax(mu_eff_face, 1.0e-30);

        /*
         * Add area-weighted statistics.
         */
        diag->area += face_area;

        t4_wmd_add_stat(
            diag,
            T4_WMD_YPLUS_GRAD_CENTROID,
            yplus_grad,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_YPLUS_SPALDING_CENTROID,
            yplus_spalding,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_YPLUS_SPALDING_OVER_GRAD,
            yplus_ratio,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_UTAU,
            u_tau,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_TAU_W,
            tau_w_mag,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_CF,
            Cf,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_H_N_PENALTY,
            h_n,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_OWNER_FACE_HEIGHT,
            face_height,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_H_T,
            h_t,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_ASPECT_RATIO,
            aspect_ratio,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_NONORTHOGONALITY_DEG,
            nonorthogonality_deg,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_OWNER_TET_QUALITY,
            tet_quality,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_MU_EFF_OVER_MU,
            mu_eff_over_mu,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_BETA_DT_OVER_BETA_MU,
            beta_dt_over_beta_mu,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_PENALTY_TO_SHEAR,
            penalty_to_shear,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_OWNER_CFL,
            CFL,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_OWNER_RE_H,
            Re_h,
            face_area);

        t4_wmd_add_stat(
            diag,
            T4_WMD_OWNER_PE_EFF,
            Pe_eff,
            face_area);

        /*
         * y+ area bins.
         */
        if (yplus_grad < 1.0) {
            diag->yplus_bin_area[0] += face_area;
        } else if (yplus_grad < 5.0) {
            diag->yplus_bin_area[1] += face_area;
        } else if (yplus_grad < 30.0) {
            diag->yplus_bin_area[2] += face_area;
        } else if (yplus_grad < 100.0) {
            diag->yplus_bin_area[3] += face_area;
        } else {
            diag->yplus_bin_area[4] += face_area;
        }

        /*
         * Reverse streamwise wall shear.
         *
         * n is the outward normal of the fluid domain. Therefore
         *
         *   tau_w_vec = mu P_t (grad u + grad u^T) n
         *
         * is the tangential traction exerted by the wall on the fluid.
         * For attached forward flow this traction points opposite eU_t,
         * so tau_stream < 0 is the normal attached-flow sign. Reverse
         * wall shear is detected by tau_stream > 0.
         */
        const double eU_dot_n =
            t4_wmd_dot3(
                eU,
                n);

        double eU_t[3] = {
            eU[0] - eU_dot_n*n[0],
            eU[1] - eU_dot_n*n[1],
            eU[2] - eU_dot_n*n[2]
        };

        const double eU_t_norm =
            t4_wmd_norm3(eU_t);

        if (eU_t_norm > 1.0e-12) {
            eU_t[0] /= eU_t_norm;
            eU_t[1] /= eU_t_norm;
            eU_t[2] /= eU_t_norm;

            const double tau_stream =
                t4_wmd_dot3(
                    tau_w_vec,
                    eU_t);

            diag->streamwise_shear_valid_area +=
                face_area;

            if (tau_stream > 0.0) {
                diag->reverse_shear_area +=
                    face_area;
            }
        }

        if (nonorthogonality_deg > 30.0) {
            diag->nonorth_gt_30_area += face_area;
        }

        if (nonorthogonality_deg > 60.0) {
            diag->nonorth_gt_60_area += face_area;
        }

        if (tet_quality < 0.1) {
            diag->quality_lt_0p1_area += face_area;
        }

        if (tet_quality < 0.2) {
            diag->quality_lt_0p2_area += face_area;
        }

        if (mu_eff_over_mu > 2.0) {
            diag->mu_ratio_gt_2_area += face_area;
        }

        if (mu_eff_over_mu > 10.0) {
            diag->mu_ratio_gt_10_area += face_area;
        }

        if (beta_dt > beta_mu) {
            diag->beta_dt_dominant_area += face_area;
        }

        /*
         * Store gradient-based y+ for VTK visualization.
         *
         * Maximum over incident wall faces is used.
         */
        if (
            write_nodal_yplus &&
            vals->y_plus != NULL &&
            vals->y_plus_count != NULL) {

            for (int a = 0; a < 3; ++a) {
                const int gid = surf->conn[es][a];

                vals->y_plus[gid] =
                    fmax(
                        vals->y_plus[gid],
                        yplus_grad);

                vals->y_plus_count[gid] = 1;
            }
        }

        (void)p_face;
        (void)tau_face;
        (void)tau_c_face;
    }

    free(owner);
    free(lface);

    t4_wmd_finalize(diag);

    return 0;
}

/* -------------------------------------------------------------------------- */
/* MPI reduction                                                              */
/* -------------------------------------------------------------------------- */

void wall_mesh_diagnostics_allreduce(
    T4WallMeshDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    if (
        diag == NULL ||
        monolis_com == NULL) {

        return;
    }

    enum {
        NSTAT = T4_WALL_MESH_DIAG_NUM_STATS,
        NBIN = T4_WALL_MESH_DIAG_NUM_YPLUS_BINS
    };

    /*
     * SUM pack:
     *
     * area
     * sum[NSTAT]
     * sum2[NSTAT]
     * yplus bins[NBIN]
     * streamwise valid area
     * reverse area
     * nonorth 30/60
     * quality 0.1/0.2
     * mu ratio 2/10
     * beta_dt dominant
     * invalid face count
     */
    double sum_pack[
        1
      + 2*NSTAT
      + NBIN
      + 2
      + 2
      + 2
      + 2
      + 1
      + 1
    ];

    int pos = 0;

    sum_pack[pos++] = diag->area;

    for (int k = 0; k < NSTAT; ++k) {
        sum_pack[pos++] = diag->sum[k];
    }

    for (int k = 0; k < NSTAT; ++k) {
        sum_pack[pos++] = diag->sum2[k];
    }

    for (int b = 0; b < NBIN; ++b) {
        sum_pack[pos++] = diag->yplus_bin_area[b];
    }

    sum_pack[pos++] =
        diag->streamwise_shear_valid_area;

    sum_pack[pos++] =
        diag->reverse_shear_area;

    sum_pack[pos++] =
        diag->nonorth_gt_30_area;

    sum_pack[pos++] =
        diag->nonorth_gt_60_area;

    sum_pack[pos++] =
        diag->quality_lt_0p1_area;

    sum_pack[pos++] =
        diag->quality_lt_0p2_area;

    sum_pack[pos++] =
        diag->mu_ratio_gt_2_area;

    sum_pack[pos++] =
        diag->mu_ratio_gt_10_area;

    sum_pack[pos++] =
        diag->beta_dt_dominant_area;

    sum_pack[pos++] =
        diag->invalid_face_count;

    const int sum_pack_size = pos;

    monolis_allreduce_R(
        sum_pack_size,
        sum_pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    pos = 0;

    diag->area = sum_pack[pos++];

    for (int k = 0; k < NSTAT; ++k) {
        diag->sum[k] = sum_pack[pos++];
    }

    for (int k = 0; k < NSTAT; ++k) {
        diag->sum2[k] = sum_pack[pos++];
    }

    for (int b = 0; b < NBIN; ++b) {
        diag->yplus_bin_area[b] = sum_pack[pos++];
    }

    diag->streamwise_shear_valid_area =
        sum_pack[pos++];

    diag->reverse_shear_area =
        sum_pack[pos++];

    diag->nonorth_gt_30_area =
        sum_pack[pos++];

    diag->nonorth_gt_60_area =
        sum_pack[pos++];

    diag->quality_lt_0p1_area =
        sum_pack[pos++];

    diag->quality_lt_0p2_area =
        sum_pack[pos++];

    diag->mu_ratio_gt_2_area =
        sum_pack[pos++];

    diag->mu_ratio_gt_10_area =
        sum_pack[pos++];

    diag->beta_dt_dominant_area =
        sum_pack[pos++];

    diag->invalid_face_count =
        sum_pack[pos++];

    /*
     * Global maxima.
     */
    double max_pack[NSTAT];

    for (int k = 0; k < NSTAT; ++k) {
        max_pack[k] = diag->max[k];
    }

    monolis_allreduce_R(
        NSTAT,
        max_pack,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    for (int k = 0; k < NSTAT; ++k) {
        diag->max[k] = max_pack[k];
    }

    /*
     * Global minima through max(-min).
     */
    double negative_min_pack[NSTAT];

    for (int k = 0; k < NSTAT; ++k) {
        negative_min_pack[k] = -diag->min[k];
    }

    monolis_allreduce_R(
        NSTAT,
        negative_min_pack,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    for (int k = 0; k < NSTAT; ++k) {
        diag->min[k] = -negative_min_pack[k];
    }

    t4_wmd_finalize(diag);
}


/* -------------------------------------------------------------------------- */
/* Output one statistic                                                       */
/* -------------------------------------------------------------------------- */

static void t4_wmd_output_stat(
    const char* stat_name,
    const T4WallMeshDiagnostics* diag,
    int stat_id,
    double t,
    const char* directory)
{
    char filename[256];

    snprintf(
        filename,
        sizeof(filename),
        "wall_%s_min.txt",
        stat_name);

    ROM_std_hlpod_output_add_calc_time(
        diag->min[stat_id],
        t,
        filename,
        directory);

    snprintf(
        filename,
        sizeof(filename),
        "wall_%s_mean.txt",
        stat_name);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean[stat_id],
        t,
        filename,
        directory);

    snprintf(
        filename,
        sizeof(filename),
        "wall_%s_rms.txt",
        stat_name);

    ROM_std_hlpod_output_add_calc_time(
        diag->rms[stat_id],
        t,
        filename,
        directory);

    snprintf(
        filename,
        sizeof(filename),
        "wall_%s_max.txt",
        stat_name);

    ROM_std_hlpod_output_add_calc_time(
        diag->max[stat_id],
        t,
        filename,
        directory);
}


/* -------------------------------------------------------------------------- */
/* Wall diagnostic output                                                     */
/* -------------------------------------------------------------------------- */

void output_wall_mesh_diagnostics_tet4_tri3(
    const T4WallMeshDiagnostics* diag,
    double t,
    const char* directory)
{
    if (
        diag == NULL ||
        directory == NULL) {

        return;
    }

    ROM_std_hlpod_output_add_calc_time(
        diag->area,
        t,
        "wall_meshdiag_area.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->invalid_face_count,
        t,
        "wall_meshdiag_invalid_face_count.txt",
        directory);

    for (int k = 0;
         k < T4_WALL_MESH_DIAG_NUM_STATS;
         ++k) {

        t4_wmd_output_stat(
            t4_wmd_stat_name[k],
            diag,
            k,
            t,
            directory);
    }

    ROM_std_hlpod_output_add_calc_time(
        diag->yplus_bin_fraction[0],
        t,
        "wall_area_fraction_yplus_lt_1.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->yplus_bin_fraction[1],
        t,
        "wall_area_fraction_yplus_1_5.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->yplus_bin_fraction[2],
        t,
        "wall_area_fraction_yplus_5_30.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->yplus_bin_fraction[3],
        t,
        "wall_area_fraction_yplus_30_100.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->yplus_bin_fraction[4],
        t,
        "wall_area_fraction_yplus_ge_100.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->reverse_shear_area_fraction,
        t,
        "wall_reverse_shear_area_fraction.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->nonorth_gt_30_fraction,
        t,
        "wall_area_fraction_nonorth_gt_30deg.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->nonorth_gt_60_fraction,
        t,
        "wall_area_fraction_nonorth_gt_60deg.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->quality_lt_0p1_fraction,
        t,
        "wall_area_fraction_owner_tet_quality_lt_0p1.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->quality_lt_0p2_fraction,
        t,
        "wall_area_fraction_owner_tet_quality_lt_0p2.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mu_ratio_gt_2_fraction,
        t,
        "wall_area_fraction_mu_eff_over_mu_gt_2.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mu_ratio_gt_10_fraction,
        t,
        "wall_area_fraction_mu_eff_over_mu_gt_10.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->beta_dt_dominant_fraction,
        t,
        "wall_area_fraction_beta_dt_dominant.txt",
        directory);
}


/* ========================================================================== */
/* Domain continuity diagnostics                                              */
/* ========================================================================== */

static void t4_domain_continuity_init(
    T4DomainContinuityDiagnostics* diag)
{
    memset(diag, 0, sizeof(*diag));
}


static void t4_domain_continuity_finalize(
    T4DomainContinuityDiagnostics* diag)
{
    if (diag->volume > 0.0) {
        diag->mean_abs_div_u =
            diag->int_abs_div_u / diag->volume;

        diag->rms_div_u =
            sqrt(fmax(
                diag->int_div_u2 / diag->volume,
                0.0));
    } else {
        diag->mean_abs_div_u = 0.0;
        diag->rms_div_u = 0.0;
        diag->max_abs_div_u = 0.0;
    }
}


int calc_domain_continuity_diagnostics_tet4(
    const BBFE_DATA* fe,
    const VALUES* vals,
    T4DomainContinuityDiagnostics* diag)
{
    if (
        fe == NULL ||
        vals == NULL ||
        diag == NULL) {

        return -1;
    }

    t4_domain_continuity_init(diag);

    if (fe->local_num_nodes != 4) {
        fprintf(
            stderr,
            "[continuity_diag] "
            "TET4 expected, got local_num_nodes=%d\n",
            fe->local_num_nodes);

        return -1;
    }

    for (int ke = 0;
         ke < fe->total_num_elems;
         ++ke) {

        double J[3][3];
        double invJ[3][3];
        double dN_dx[4][3];
        double detJ = 0.0;
        double volume = 0.0;

        if (!t4_wmd_tet4_geometry(
                fe,
                ke,
                J,
                invJ,
                dN_dx,
                &detJ,
                &volume)) {

            diag->invalid_element_count += 1.0;
            continue;
        }

        if (
            volume <= 1.0e-300 ||
            !isfinite(volume)) {

            diag->invalid_element_count += 1.0;
            continue;
        }

        double div_u = 0.0;

        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];

            div_u +=
                vals->v[gid][0] * dN_dx[a][0]
              + vals->v[gid][1] * dN_dx[a][1]
              + vals->v[gid][2] * dN_dx[a][2];
        }

        const double abs_div_u =
            fabs(div_u);

        diag->volume += volume;

        diag->int_abs_div_u +=
            abs_div_u * volume;

        diag->int_div_u2 +=
            div_u * div_u * volume;

        diag->max_abs_div_u =
            fmax(
                diag->max_abs_div_u,
                abs_div_u);
    }

    t4_domain_continuity_finalize(diag);

    return 0;
}


void domain_continuity_diagnostics_allreduce(
    T4DomainContinuityDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    if (
        diag == NULL ||
        monolis_com == NULL) {

        return;
    }

    double sum_pack[4] = {
        diag->volume,
        diag->int_abs_div_u,
        diag->int_div_u2,
        diag->invalid_element_count
    };

    monolis_allreduce_R(
        4,
        sum_pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    diag->volume =
        sum_pack[0];

    diag->int_abs_div_u =
        sum_pack[1];

    diag->int_div_u2 =
        sum_pack[2];

    diag->invalid_element_count =
        sum_pack[3];

    double max_pack[1] = {
        diag->max_abs_div_u
    };

    monolis_allreduce_R(
        1,
        max_pack,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    diag->max_abs_div_u =
        max_pack[0];

    t4_domain_continuity_finalize(diag);
}


void output_domain_continuity_diagnostics_tet4(
    const T4DomainContinuityDiagnostics* diag,
    double t,
    const char* directory)
{
    if (
        diag == NULL ||
        directory == NULL) {

        return;
    }

    ROM_std_hlpod_output_add_calc_time(
        diag->volume,
        t,
        "domain_volume.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_abs_div_u,
        t,
        "domain_div_u_mean_abs.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->rms_div_u,
        t,
        "domain_div_u_rms.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->max_abs_div_u,
        t,
        "domain_div_u_max_abs.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->invalid_element_count,
        t,
        "domain_invalid_tet_count.txt",
        directory);
}

/* normal = dx/dr x dx/ds.  Its norm is 2*area for an affine TRI3. */
static inline void t4_tri3_normal(
    const double xtri[3][3],
    const double* dN_dr,
    const double* dN_ds,
    double normal[3])
{
    double xr[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};

    for (int a = 0; a < 3; ++a) {
        xr[0] += dN_dr[a] * xtri[a][0];
        xr[1] += dN_dr[a] * xtri[a][1];
        xr[2] += dN_dr[a] * xtri[a][2];
        xs[0] += dN_ds[a] * xtri[a][0];
        xs[1] += dN_ds[a] * xtri[a][1];
        xs[2] += dN_ds[a] * xtri[a][2];
    }

    t4_cross3(xr, xs, normal);
}

/* Correct quadrature-weight convention for reference triangle.
 * If sum(weights)=0.5, returns 1.0. If sum(weights)=1.0, returns 0.5.
 */
static double t4_tri_area_factor(const BBFE_BASIS* basis_surf)
{
    double sw = 0.0;
    for (int p = 0; p < basis_surf->num_integ_points; ++p) {
        sw += basis_surf->integ_weight[p];
    }

    if (fabs(sw) <= 1.0e-300) return 1.0;
    return 0.5 / sw;
}

static int t4_inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0 / det;
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static double t4_det3x3(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

/* TET4 reference gradients for N0=1-r-s-t, N1=r, N2=s, N3=t. */
static void t4_shape_grad_ref(double dN_dxi[4][3])
{
    dN_dxi[0][0] = -1.0; dN_dxi[0][1] = -1.0; dN_dxi[0][2] = -1.0;
    dN_dxi[1][0] =  1.0; dN_dxi[1][1] =  0.0; dN_dxi[1][2] =  0.0;
    dN_dxi[2][0] =  0.0; dN_dxi[2][1] =  1.0; dN_dxi[2][2] =  0.0;
    dN_dxi[3][0] =  0.0; dN_dxi[3][1] =  0.0; dN_dxi[3][2] =  1.0;
}

static int t4_grad_phys_metric(
    const double xvol[4][3],
    double dN_dx[4][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double dN_dxi[4][3];
    t4_shape_grad_ref(dN_dxi);

    double J[3][3] = {{0.0}};
    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = t4_det3x3(J);
    if (fabs(detJ) < 1.0e-300) return 0;
    if (!t4_inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);  /* 6 * physical tet volume */

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

static void t4_grad_tet4_exact(int ke, const BBFE_DATA* fe, const VALUES* vals, double G[3][3])
{
    memset(G, 0, sizeof(double) * 9);

    if (fe->local_num_nodes != 4) return;

    double xvol[4][3];
    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];
        xvol[a][0] = fe->x[gid][0];
        xvol[a][1] = fe->x[gid][1];
        xvol[a][2] = fe->x[gid][2];
    }

    double dN_dx[4][3];
    double J_inv[3][3];
    double Jacobian = 0.0;
    if (!t4_grad_phys_metric(xvol, dN_dx, J_inv, &Jacobian)) return;

    for (int a = 0; a < 4; ++a) {
        const int gid = fe->conn[ke][a];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                G[i][j] += vals->v[gid][i] * dN_dx[a][j];
            }
        }
    }
}

static void t4_compute_all_gradients_tet4_exact(
    const BBFE_DATA* fe,
    const VALUES* vals,
    double (*G_all)[3][3])
{
    const int ne = fe->total_num_elems;

    #pragma omp parallel for if(ne > 256) schedule(static) default(none) shared(fe, vals, G_all, ne)
    for (int ke = 0; ke < ne; ++ke) {
        t4_grad_tet4_exact(ke, fe, vals, G_all[ke]);
    }
}

static int t4_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for (int a = 0; a < 3; ++a) {
        s2v[a] = -1;
        const int sgid = surf->conn[es][a];
        for (int b = 0; b < 4; ++b) {
            if (fe->conn[ke][b] == sgid) {
                s2v[a] = b;
                break;
            }
        }
        if (s2v[a] < 0) return 0;
    }
    return 1;
}

/* Pressure contribution on TRI3/TET4 surface: integrate -p n dS. */
void calc_Cd_p_tet4_tri3(
    BBFE_DATA* surf,
    BBFE_DATA* fe,
    BBFE_BASIS* basis_surf,
    VALUES* vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out)
{
    if (D_out) *D_out = 0.0;
    if (L_out) *L_out = 0.0;

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        return;
    }

    double eU[3] = { eU_in[0], eU_in[1], eU_in[2] };
    double eP[3] = { eP_in[0], eP_in[1], eP_in[2] };
    t4_normalize3(eU);

    /* Gram-Schmidt: make eP perpendicular to eU. */
    const double ep_dot_eu = t4_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];
    t4_normalize3(eP);

    const double area_factor = t4_tri_area_factor(basis_surf);
    double Fp[3] = {0.0, 0.0, 0.0};

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0; xf[1] /= 3.0; xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xc[0] += fe->x[gid][0];
            xc[1] += fe->x[gid][1];
            xc[2] += fe->x[gid][2];
        }
        xc[0] /= 4.0; xc[1] /= 4.0; xc[2] /= 4.0;

        const double svec[3] = { xf[0] - xc[0], xf[1] - xc[1], xf[2] - xc[2] };

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(xtri, basis_surf->dN_dxi[p], basis_surf->dN_det[p], nraw);

            const double Jraw = t4_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = { nraw[0]/Jraw, nraw[1]/Jraw, nraw[2]/Jraw };
            if (t4_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double p_q = 0.0;
            for (int a = 0; a < 3; ++a) {
                const int gid = surf->conn[es][a];
                p_q += basis_surf->N[p][a] * vals->p[gid];
            }
//              p_q = 1.0;
            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;
            Fp[0] += (p_q) * n[0] * w;
            Fp[1] += (p_q) * n[1] * w;
            Fp[2] += (p_q) * n[2] * w;
        }
    }

    if (D_out) *D_out = t4_dot3(Fp, eU);
    if (L_out) *L_out = t4_dot3(Fp, eP);

    free(owner);
    free(lface);
}

/* Viscous contribution on TRI3/TET4 surface.
 * Sign convention follows your current calc_Cd_v:
 *   integrate tau*n on the fluid side, then report body force = -integral.
 * The current HEX8 code stores D/L, not D/qA or L/qA, in Cd_out/Cl_out;
 * this function keeps that behavior. Change the final two assignments if you
 * want nondimensional coefficients.
 */


/* -------------------------------------------------------------------------- */
/* Small TET4/TRI3 utilities                                                  */
/* -------------------------------------------------------------------------- */

typedef long long t4n_node_id_t;

static inline double t4n_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double t4n_norm3(const double a[3])
{
    return sqrt(t4n_dot3(a, a));
}

static inline void t4n_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static inline void t4n_sort3(t4n_node_id_t k[3])
{
    if (k[1] < k[0]) { t4n_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
    if (k[2] < k[1]) { t4n_node_id_t t = k[1]; k[1] = k[2]; k[2] = t; }
    if (k[1] < k[0]) { t4n_node_id_t t = k[0]; k[0] = k[1]; k[1] = t; }
}

static inline int t4n_cmp3(const t4n_node_id_t a[3], const t4n_node_id_t b[3])
{
    if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : +1;
    if (a[1] != b[1]) return (a[1] < b[1]) ? -1 : +1;
    if (a[2] != b[2]) return (a[2] < b[2]) ? -1 : +1;
    return 0;
}

/* Faces of TET4. Orientation is not used for matching because the keys are sorted.
 * Outward normal is corrected geometrically by the vector from cell centroid to face centroid.
 */
static const int T4N_FACE[4][3] = {
    {1, 2, 3},
    {0, 3, 2},
    {0, 1, 3},
    {0, 2, 1}
};

typedef struct {
    t4n_node_id_t key[3];
    int owner_elem;
    int owner_lface;
} t4n_FaceEntry;

typedef struct {
    t4n_node_id_t key[3];
    int surf_face_id;
} t4n_SurfFaceEntry;

static int t4n_cmp_face(const void* A, const void* B)
{
    const t4n_FaceEntry* a = (const t4n_FaceEntry*)A;
    const t4n_FaceEntry* b = (const t4n_FaceEntry*)B;
    int c = t4n_cmp3(a->key, b->key);
    if (c) return c;
    if (a->owner_elem != b->owner_elem) return (a->owner_elem < b->owner_elem) ? -1 : +1;
    if (a->owner_lface != b->owner_lface) return (a->owner_lface < b->owner_lface) ? -1 : +1;
    return 0;
}

static int t4n_cmp_sface(const void* A, const void* B)
{
    const t4n_SurfFaceEntry* a = (const t4n_SurfFaceEntry*)A;
    const t4n_SurfFaceEntry* b = (const t4n_SurfFaceEntry*)B;
    int c = t4n_cmp3(a->key, b->key);
    if (c) return c;
    if (a->surf_face_id != b->surf_face_id) return (a->surf_face_id < b->surf_face_id) ? -1 : +1;
    return 0;
}

static t4n_FaceEntry* t4n_build_vol_faces_sorted(const BBFE_DATA* fe, int* nfaces_out)
{
    const int nfaces = fe->total_num_elems * 4;
    t4n_FaceEntry* arr = (t4n_FaceEntry*)malloc(sizeof(t4n_FaceEntry) * nfaces);
    if (arr == NULL) {
        fprintf(stderr, "[t4n][ERR] malloc vol faces failed\n");
        *nfaces_out = 0;
        return NULL;
    }

    int k = 0;
    for (int e = 0; e < fe->total_num_elems; ++e) {
        for (int f = 0; f < 4; ++f) {
            t4n_node_id_t key[3] = {
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][0]],
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][1]],
                (t4n_node_id_t)fe->conn[e][T4N_FACE[f][2]]
            };
            t4n_sort3(key);
            for (int j = 0; j < 3; ++j) arr[k].key[j] = key[j];
            arr[k].owner_elem  = e;
            arr[k].owner_lface = f;
            ++k;
        }
    }

    qsort(arr, nfaces, sizeof(t4n_FaceEntry), t4n_cmp_face);
    *nfaces_out = nfaces;
    return arr;
}

static t4n_SurfFaceEntry* t4n_build_surf_faces_sorted(const BBFE_DATA* surf, int* nsurf_out)
{
    const int nsurf = surf->total_num_elems;
    t4n_SurfFaceEntry* arr = (t4n_SurfFaceEntry*)malloc(sizeof(t4n_SurfFaceEntry) * nsurf);
    if (arr == NULL) {
        fprintf(stderr, "[t4n][ERR] malloc surf faces failed\n");
        *nsurf_out = 0;
        return NULL;
    }

    for (int es = 0; es < nsurf; ++es) {
        t4n_node_id_t key[3] = {
            (t4n_node_id_t)surf->conn[es][0],
            (t4n_node_id_t)surf->conn[es][1],
            (t4n_node_id_t)surf->conn[es][2]
        };
        t4n_sort3(key);
        for (int j = 0; j < 3; ++j) arr[es].key[j] = key[j];
        arr[es].surf_face_id = es;
    }

    qsort(arr, nsurf, sizeof(t4n_SurfFaceEntry), t4n_cmp_sface);
    *nsurf_out = nsurf;
    return arr;
}

static int t4n_build_owner_map_mergejoin(
    const t4n_FaceEntry* volF,
    int nVolF,
    const t4n_SurfFaceEntry* surfF,
    int nSurfF,
    int* owner_elem,
    int* owner_lface)
{
    int i = 0;
    int j = 0;
    int matches = 0;

    if (
        volF == NULL ||
        surfF == NULL ||
        owner_elem == NULL ||
        owner_lface == NULL
    ) {
        return 0;
    }

    while (
        i < nVolF &&
        j < nSurfF
    ) {
        const int c =
            t4n_cmp3(
                volF[i].key,
                surfF[j].key
            );

        if (c < 0) {
            /*
             * volume側のキーが小さい。
             * 同じキーを持つvolume faceをまとめて進める。
             */
            const t4n_node_id_t key0[3] = {
                volF[i].key[0],
                volF[i].key[1],
                volF[i].key[2]
            };

            do {
                ++i;
            } while (
                i < nVolF &&
                t4n_cmp3(
                    key0,
                    volF[i].key
                ) == 0
            );

            continue;
        }

        if (c > 0) {
            /*
             * surface側のキーが小さい。
             *
             * このrankには対応するvolume faceがない。
             * owner[]は初期値-1のままにする。
             */
            const t4n_node_id_t key0[3] = {
                surfF[j].key[0],
                surfF[j].key[1],
                surfF[j].key[2]
            };

            do {
                ++j;
            } while (
                j < nSurfF &&
                t4n_cmp3(
                    key0,
                    surfF[j].key
                ) == 0
            );

            continue;
        }

        /*
         * volumeとsurfaceの面キーが一致。
         */
        const t4n_node_id_t key0[3] = {
            surfF[j].key[0],
            surfF[j].key[1],
            surfF[j].key[2]
        };

        /*
         * 同じキーを持つvolume faceが複数ある場合は、
         * 従来どおり最初のownerを使用する。
         */
        const int owner =
            volF[i].owner_elem;

        const int local_face =
            volF[i].owner_lface;

        /*
         * 重要:
         *
         * 同じキーを持つsurface要素をすべて対応付ける。
         *
         * 従来実装はsurface側を1個しか進めなかったため、
         * 同一キーの2個目以降が未対応になっていた。
         */
        do {
            const int es =
                surfF[j].surf_face_id;

            owner_elem[es] =
                owner;

            owner_lface[es] =
                local_face;

            ++matches;
            ++j;
        } while (
            j < nSurfF &&
            t4n_cmp3(
                key0,
                surfF[j].key
            ) == 0
        );

        /*
         * 同じキーを持つvolume faceをまとめて進める。
         */
        do {
            ++i;
        } while (
            i < nVolF &&
            t4n_cmp3(
                key0,
                volF[i].key
            ) == 0
        );
    }

    return matches;
}

static int t4n_make_owner_map_overlap(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner_out,
    int** lface_out,
    int* num_unmapped_out)
{
    if (
        owner_out == NULL ||
        lface_out == NULL
    ) {
        return -1;
    }

    *owner_out = NULL;
    *lface_out = NULL;

    if (num_unmapped_out != NULL) {
        *num_unmapped_out = 0;
    }

    if (
        surf == NULL ||
        fe == NULL
    ) {
        return -1;
    }

    if (
        fe->local_num_nodes != 4 ||
        surf->local_num_nodes != 3
    ) {
        fprintf(
            stderr,
            "[t4n-overlap] expected TET4/TRI3, "
            "but fe->local_num_nodes=%d "
            "surf->local_num_nodes=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes
        );

        return -1;
    }

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[t4n-overlap] invalid surface element count: %d\n",
            surf->total_num_elems
        );

        return -1;
    }

    /*
     * 表面要素を持たないrank。
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    int nVolF = 0;
    int nSurfF = 0;

    t4n_FaceEntry* volF =
        t4n_build_vol_faces_sorted(
            fe,
            &nVolF
        );

    t4n_SurfFaceEntry* surfF =
        t4n_build_surf_faces_sorted(
            surf,
            &nSurfF
        );

    if (
        volF == NULL ||
        surfF == NULL
    ) {
        free(volF);
        free(surfF);

        return -1;
    }

    if (nSurfF != surf->total_num_elems) {
        fprintf(
            stderr,
            "[t4n-overlap] inconsistent surface count: "
            "built=%d expected=%d\n",
            nSurfF,
            surf->total_num_elems
        );

        free(volF);
        free(surfF);

        return -1;
    }

    int* owner =
        (int*)malloc(
            sizeof(int)
            * (size_t)nSurfF
        );

    int* lface =
        (int*)malloc(
            sizeof(int)
            * (size_t)nSurfF
        );

    if (
        owner == NULL ||
        lface == NULL
    ) {
        free(volF);
        free(surfF);
        free(owner);
        free(lface);

        return -1;
    }

    /*
     * -1は「このrankにowner TETがない」ことを表す。
     */
    for (int es = 0; es < nSurfF; ++es) {
        owner[es] = -1;
        lface[es] = -1;
    }

    const int matched =
        t4n_build_owner_map_mergejoin(
            volF,
            nVolF,
            surfF,
            nSurfF,
            owner,
            lface
        );

    int num_unmapped = 0;
    int num_printed = 0;

    for (int es = 0; es < nSurfF; ++es) {
        if (owner[es] >= 0) {
            continue;
        }

        ++num_unmapped;

        /*
         * ログが過大にならないよう、最初の16面だけ詳細表示。
         */
        if (num_printed < 16) {
            fprintf(
                stderr,
                "[t4n-overlap] surf tri %d has NO local owner: "
                "nodes=(%d,%d,%d)\n",
                es,
                surf->conn[es][0],
                surf->conn[es][1],
                surf->conn[es][2]
            );

            ++num_printed;
        }
    }

    if (num_unmapped > 0) {
        fprintf(
            stderr,
            "[t4n-overlap] %d / %d surface triangles "
            "have no local owner TET; "
            "only these copies will be skipped. "
            "matched=%d\n",
            num_unmapped,
            nSurfF,
            matched
        );
    }

    free(volF);
    free(surfF);

    *owner_out = owner;
    *lface_out = lface;

    if (num_unmapped_out != NULL) {
        *num_unmapped_out =
            num_unmapped;
    }

    /*
     * ownerなし面が存在しても、配列構築自体は成功。
     */
    return 0;
}

static int t4n_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for (int a = 0; a < 3; ++a) {
        s2v[a] = -1;
        const int sgid = surf->conn[es][a];
        for (int b = 0; b < 4; ++b) {
            if (fe->conn[ke][b] == sgid) {
                s2v[a] = b;
                break;
            }
        }
        if (s2v[a] < 0) return 0;
    }
    return 1;
}

/* normal = dx/dr x dx/ds. For affine TRI3, |normal| = 2 * physical area. */
static inline void t4n_tri3_normal(
    const double xtri[3][3],
    const double* dN_dr,
    const double* dN_ds,
    double normal[3])
{
    double xr[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};

    for (int a = 0; a < 3; ++a) {
        xr[0] += dN_dr[a] * xtri[a][0];
        xr[1] += dN_dr[a] * xtri[a][1];
        xr[2] += dN_dr[a] * xtri[a][2];
        xs[0] += dN_ds[a] * xtri[a][0];
        xs[1] += dN_ds[a] * xtri[a][1];
        xs[2] += dN_ds[a] * xtri[a][2];
    }

    t4n_cross3(xr, xs, normal);
}

/* If sum(weights)=0.5, returns 1.0. If sum(weights)=1.0, returns 0.5. */
static double t4n_tri_area_factor(const BBFE_BASIS* basis_surf)
{
    double sw = 0.0;
    for (int p = 0; p < basis_surf->num_integ_points; ++p) {
        sw += basis_surf->integ_weight[p];
    }

    if (fabs(sw) <= 1.0e-300) return 1.0;
    return 0.5 / sw;
}

static double t4n_det3x3(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int t4n_inv3x3(const double A[3][3], double invA[3][3])
{
    const double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    const double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    const double a20=A[2][0], a21=A[2][1], a22=A[2][2];

    const double c00 =  a11*a22 - a12*a21;
    const double c01 =-(a10*a22 - a12*a20);
    const double c02 =  a10*a21 - a11*a20;
    const double c10 =-(a01*a22 - a02*a21);
    const double c11 =  a00*a22 - a02*a20;
    const double c12 =-(a00*a21 - a01*a20);
    const double c20 =  a01*a12 - a02*a11;
    const double c21 =-(a00*a12 - a02*a10);
    const double c22 =  a00*a11 - a01*a10;

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (fabs(det) < 1.0e-300) return 0;

    const double id = 1.0 / det;
    invA[0][0] = c00*id; invA[0][1] = c10*id; invA[0][2] = c20*id;
    invA[1][0] = c01*id; invA[1][1] = c11*id; invA[1][2] = c21*id;
    invA[2][0] = c02*id; invA[2][1] = c12*id; invA[2][2] = c22*id;
    return 1;
}

static void t4n_shape_grad_ref(double dN_dxi[4][3])
{
    dN_dxi[0][0] = -1.0; dN_dxi[0][1] = -1.0; dN_dxi[0][2] = -1.0;
    dN_dxi[1][0] =  1.0; dN_dxi[1][1] =  0.0; dN_dxi[1][2] =  0.0;
    dN_dxi[2][0] =  0.0; dN_dxi[2][1] =  1.0; dN_dxi[2][2] =  0.0;
    dN_dxi[3][0] =  0.0; dN_dxi[3][1] =  0.0; dN_dxi[3][2] =  1.0;
}

static int t4n_grad_phys_metric(
    const double xvol[4][3],
    double dN_dx[4][3],
    double J_inv[3][3],
    double* Jacobian)
{
    double dN_dxi[4][3];
    t4n_shape_grad_ref(dN_dxi);

    double J[3][3] = {{0.0}};
    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                J[i][j] += xvol[a][i] * dN_dxi[a][j];
            }
        }
    }

    const double detJ = t4n_det3x3(J);
    if (fabs(detJ) < 1.0e-300) return 0;
    if (!t4n_inv3x3(J, J_inv)) return 0;

    *Jacobian = fabs(detJ);  /* 6 * physical TET volume */

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            dN_dx[a][i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dN_dx[a][i] += dN_dxi[a][j] * J_inv[j][i];
            }
        }
    }

    return 1;
}

/* -------------------------------------------------------------------------- */
/* Feature-aware nodal projector for the STRONG_NORMAL_* modes                */
/*                                                                            */
/* MONOLIS' existing Dirichlet wrapper is component-wise in global xyz.       */
/* Arbitrary no-penetration constraints u.n=0 therefore cannot be expressed   */
/* directly without changing the global unknown basis / adding MPCs.          */
/*                                                                            */
/* For the TET4-only branch we enforce the admissible nodal velocity subspace  */
/* explicitly: incident face normals are clustered by feature angle, an       */
/* orthonormal constrained-normal basis is built, and both the state and each  */
/* Newton velocity correction are projected onto its orthogonal complement.   */
/*                                                                            */
/* Rank 1 -> smooth wall: two tangential velocity directions remain.           */
/* Rank 2 -> sharp edge : only the edge tangent remains.                       */
/* Rank 3 -> corner     : velocity is fixed to Uw.                             */
/*                                                                            */
/* This is exact as a nodal state/update projection.  It is not an algebraic   */
/* change-of-basis elimination of the sparse matrix, so the code labels this   */
/* implementation "projected strong normal" in run-time messages.             */
/* -------------------------------------------------------------------------- */

typedef struct T4WallNodeFeature {
    int is_wall;
    int ncluster;
    int rank;
    double cluster_sum[3][3];
    double cluster_weight[3];
    double qn[3][3];
    double Pt[3][3];
} T4WallNodeFeature;

static T4WallNodeFeature* g_t4_wall_node_feature = NULL;
static int g_t4_wall_node_feature_nnode = 0;
static const BBFE_DATA* g_t4_wall_projector_surf_ptr = NULL;
static const BBFE_DATA* g_t4_wall_projector_fe_ptr = NULL;
static double g_t4_wall_projector_feature_angle = -1.0;

static void t4_wall_projector_release(void)
{
    free(g_t4_wall_node_feature);
    g_t4_wall_node_feature = NULL;
    g_t4_wall_node_feature_nnode = 0;
    g_t4_wall_projector_surf_ptr = NULL;
    g_t4_wall_projector_fe_ptr = NULL;
    g_t4_wall_projector_feature_angle = -1.0;
}

static double t4_wall_norm3_local(const double a[3])
{
    return sqrt(t4n_dot3(a, a));
}

static void t4_wall_normalize3_local(double a[3])
{
    const double nrm = t4_wall_norm3_local(a);
    if (nrm > 1.0e-30) {
        a[0] /= nrm;
        a[1] /= nrm;
        a[2] /= nrm;
    }
}

static int t4_wall_projector_prepare(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe)
{
    if (surf == NULL || fe == NULL) return -1;
    if (fe->local_num_nodes != 4 || surf->local_num_nodes != 3) return -1;

    if (g_t4_wall_node_feature != NULL &&
        g_t4_wall_node_feature_nnode == fe->total_num_nodes &&
        g_t4_wall_projector_surf_ptr == surf &&
        g_t4_wall_projector_fe_ptr == fe &&
        fabs(g_t4_wall_projector_feature_angle - g_t4_wall_feature_angle_deg) < 1.0e-12) {
        return 0;
    }

    t4_wall_projector_release();

    g_t4_wall_node_feature = (T4WallNodeFeature*)calloc(
        (size_t)fe->total_num_nodes,
        sizeof(T4WallNodeFeature));
    if (g_t4_wall_node_feature == NULL) return -1;

    g_t4_wall_node_feature_nnode = fe->total_num_nodes;
    g_t4_wall_projector_surf_ptr = surf;
    g_t4_wall_projector_fe_ptr = fe;
    g_t4_wall_projector_feature_angle = g_t4_wall_feature_angle_deg;

    const double cos_threshold =
        cos(g_t4_wall_feature_angle_deg * (3.14159265358979323846 / 180.0));

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int g0 = surf->conn[es][0];
        const int g1 = surf->conn[es][1];
        const int g2 = surf->conn[es][2];

        if (g0 < 0 || g0 >= fe->total_num_nodes ||
            g1 < 0 || g1 >= fe->total_num_nodes ||
            g2 < 0 || g2 >= fe->total_num_nodes) {
            t4_wall_projector_release();
            return -1;
        }

        const double e1[3] = {
            fe->x[g1][0] - fe->x[g0][0],
            fe->x[g1][1] - fe->x[g0][1],
            fe->x[g1][2] - fe->x[g0][2]
        };
        const double e2[3] = {
            fe->x[g2][0] - fe->x[g0][0],
            fe->x[g2][1] - fe->x[g0][1],
            fe->x[g2][2] - fe->x[g0][2]
        };

        double n[3];
        t4n_cross3(e1, e2, n);
        const double nraw = t4_wall_norm3_local(n);
        if (!isfinite(nraw) || nraw <= 1.0e-30) continue;
        const double area = 0.5 * nraw;
        n[0] /= nraw;
        n[1] /= nraw;
        n[2] /= nraw;

        const int gids[3] = {g0, g1, g2};
        for (int aa = 0; aa < 3; ++aa) {
            T4WallNodeFeature* f = &g_t4_wall_node_feature[gids[aa]];
            f->is_wall = 1;

            int best = -1;
            double best_abs_dot = -1.0;
            double best_signed_dot = 1.0;

            for (int k = 0; k < f->ncluster; ++k) {
                double mean[3] = {
                    f->cluster_sum[k][0],
                    f->cluster_sum[k][1],
                    f->cluster_sum[k][2]
                };
                t4_wall_normalize3_local(mean);
                const double d = t4n_dot3(mean, n);
                const double ad = fabs(d);
                if (ad > best_abs_dot) {
                    best_abs_dot = ad;
                    best_signed_dot = d;
                    best = k;
                }
            }

            if (best >= 0 && best_abs_dot >= cos_threshold) {
                const double sgn = (best_signed_dot < 0.0 ? -1.0 : 1.0);
                for (int d = 0; d < 3; ++d) {
                    f->cluster_sum[best][d] += area * sgn * n[d];
                }
                f->cluster_weight[best] += area;
            }
            else if (f->ncluster < 3) {
                const int k = f->ncluster++;
                for (int d = 0; d < 3; ++d) {
                    f->cluster_sum[k][d] = area * n[d];
                }
                f->cluster_weight[k] = area;
            }
            else {
                /* More than three geometric clusters cannot add another
                 * independent constraint in R^3; merge into the closest one. */
                if (best < 0) best = 0;
                const double sgn = (best_signed_dot < 0.0 ? -1.0 : 1.0);
                for (int d = 0; d < 3; ++d) {
                    f->cluster_sum[best][d] += area * sgn * n[d];
                }
                f->cluster_weight[best] += area;
            }
        }
    }

    long long rank_count[4] = {0, 0, 0, 0};
    for (int gid = 0; gid < fe->total_num_nodes; ++gid) {
        T4WallNodeFeature* f = &g_t4_wall_node_feature[gid];

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                f->Pt[a][b] = (a == b ? 1.0 : 0.0);
                f->qn[a][b] = 0.0;
            }
        }

        if (!f->is_wall) continue;

        int rank = 0;
        for (int k = 0; k < f->ncluster && rank < 3; ++k) {
            double q[3] = {
                f->cluster_sum[k][0],
                f->cluster_sum[k][1],
                f->cluster_sum[k][2]
            };
            t4_wall_normalize3_local(q);

            for (int j = 0; j < rank; ++j) {
                const double alpha = t4n_dot3(q, f->qn[j]);
                for (int d = 0; d < 3; ++d) q[d] -= alpha * f->qn[j][d];
            }

            const double qnrm = t4_wall_norm3_local(q);
            if (qnrm <= 1.0e-6) continue;
            for (int d = 0; d < 3; ++d) f->qn[rank][d] = q[d] / qnrm;
            ++rank;
        }
        f->rank = rank;
        if (rank >= 0 && rank <= 3) ++rank_count[rank];

        for (int j = 0; j < rank; ++j) {
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    f->Pt[a][b] -= f->qn[j][a] * f->qn[j][b];
                }
            }
        }
    }

    if (monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s projected-strong wall projector built: feature_angle=%.3f deg "
            "rank1=%lld rank2=%lld rank3=%lld\n",
            CODENAME,
            g_t4_wall_feature_angle_deg,
            rank_count[1],
            rank_count[2],
            rank_count[3]);
    }

    return 0;
}

static void t4_wall_project_relative_vector_at_node(int gid, double v[3])
{
    if (g_t4_wall_node_feature == NULL ||
        gid < 0 || gid >= g_t4_wall_node_feature_nnode ||
        !g_t4_wall_node_feature[gid].is_wall) {
        return;
    }

    const T4WallNodeFeature* f = &g_t4_wall_node_feature[gid];
    const double vin[3] = {v[0], v[1], v[2]};
    for (int a = 0; a < 3; ++a) {
        v[a] = 0.0;
        for (int b = 0; b < 3; ++b) v[a] += f->Pt[a][b] * vin[b];
    }
}

static void t4_wall_project_velocity_state(
    double** v,
    const double Uw[3])
{
    if (v == NULL || Uw == NULL || g_t4_wall_node_feature == NULL) return;
    for (int gid = 0; gid < g_t4_wall_node_feature_nnode; ++gid) {
        if (!g_t4_wall_node_feature[gid].is_wall) continue;
        double r[3] = {
            v[gid][0] - Uw[0],
            v[gid][1] - Uw[1],
            v[gid][2] - Uw[2]
        };
        t4_wall_project_relative_vector_at_node(gid, r);
        for (int d = 0; d < 3; ++d) v[gid][d] = Uw[d] + r[d];
    }
}

static void t4_wall_project_velocity_increment(double** dv)
{
    if (dv == NULL || g_t4_wall_node_feature == NULL) return;
    for (int gid = 0; gid < g_t4_wall_node_feature_nnode; ++gid) {
        if (!g_t4_wall_node_feature[gid].is_wall) continue;
        double q[3] = {dv[gid][0], dv[gid][1], dv[gid][2]};
        t4_wall_project_relative_vector_at_node(gid, q);
        dv[gid][0] = q[0];
        dv[gid][1] = q[1];
        dv[gid][2] = q[2];
    }
}

/* Project momentum residual rows onto the admissible test subspace.  This is
 * essential for the projected-strong modes: normal test equations are not part
 * of the constrained variational space and must not enter the Newton residual
 * norm. */
static void t4_wall_project_velocity_rhs(double* rhs)
{
    if (rhs == NULL || g_t4_wall_node_feature == NULL) return;
    for (int gid = 0; gid < g_t4_wall_node_feature_nnode; ++gid) {
        if (!g_t4_wall_node_feature[gid].is_wall) continue;
        double q[3] = {
            rhs[4*gid + 0],
            rhs[4*gid + 1],
            rhs[4*gid + 2]
        };
        t4_wall_project_relative_vector_at_node(gid, q);
        rhs[4*gid + 0] = q[0];
        rhs[4*gid + 1] = q[1];
        rhs[4*gid + 2] = q[2];
    }
}

/* -------------------------------------------------------------------------- */
/* Local wall operators for all selectable wall modes                         */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Case C: normal weak essential BC + tangential all-y+ wall law             */
/*                                                                            */
/* Normal condition:                                                          */
/*                                                                            */
/*   r_n = (u-Uw).n = 0                                                       */
/*                                                                            */
/* is imposed by the same symmetric normal-only Nitsche formulation used in   */
/* Case B.  The normal penalty is therefore still                             */
/*                                                                            */
/*   beta_n = gamma_n * (mu_eff/h + c_dt rho h/dt).                           */
/*                                                                            */
/* Tangentially, no Dirichlet penalty is used.  Instead the wall traction is  */
/* modeled by a Spalding all-y+ Robin law                                     */
/*                                                                            */
/*   t_t,fluid = - beta_w u_t,                                                */
/*   beta_w     = rho u_tau^2 / |u_t|,                                        */
/*                                                                            */
/* with u_tau obtained from the Spalding relation using the molecular         */
/* viscosity.  beta_w is frozen in the tangent matrix (Picard/Newton hybrid). */
/* -------------------------------------------------------------------------- */

static void t4n_spalding_wall_coefficients(
    const double u_exchange[3],
    const double Uw[3],
    const double n[3],
    double rho,
    double mu_wall_law,
    double y_exchange,
    double ut_exchange[3],
    double* ut_norm_out,
    double* u_tau_out,
    double* beta_wall_out)
{
    const double eps = 1.0e-30;
    const double y_eff = fmax(y_exchange, 1.0e-12);
    const double mu_eff = fmax(mu_wall_law, eps);
    const double rho_eff = fmax(rho, eps);

    const double r[3] = {
        u_exchange[0] - Uw[0],
        u_exchange[1] - Uw[1],
        u_exchange[2] - Uw[2]
    };
    const double rn = t4n_dot3(r, n);
    for (int a = 0; a < 3; ++a) {
        ut_exchange[a] = r[a] - rn * n[a];
    }

    const double ut_norm = sqrt(t4n_dot3(ut_exchange, ut_exchange));
    double beta_wall = mu_eff / y_eff;
    double u_tau = 0.0;

    if (ut_norm > g_t4_wall_ut_eps) {
        u_tau = compute_utau_all_yplus(
            ut_norm,
            y_eff,
            rho_eff,
            mu_eff);

        if (isfinite(u_tau) && u_tau > 0.0) {
            beta_wall = rho_eff * u_tau * u_tau / fmax(ut_norm, g_t4_wall_ut_eps);
        }
    }

    if (!isfinite(beta_wall) || beta_wall < 0.0) {
        beta_wall = mu_eff / y_eff;
    }

    if (ut_norm_out != NULL)   *ut_norm_out = ut_norm;
    if (u_tau_out != NULL)     *u_tau_out = u_tau;
    if (beta_wall_out != NULL) *beta_wall_out = beta_wall;
}

static void t4n_wall_fullvec_nitsche_local_vecmat(
    double vec_i[4],
    double mat_ij[4][4],
    double Ni,
    double Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    double p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double gamma_n)
{
    const double eps = 1.0e-30;
    const double dt_eff = fmax(dt, eps);
    const double h_eff = fmax(h_n, 1.0e-12);

    double r[3];
    for (int a = 0; a < 3; ++a) r[a] = u_q[a] - Uw[a];
    const double rn = t4n_dot3(r, n);

    double beta = 0.0;
    t4n_get_penalty_coefficients(
        rho, mu_eff, h_eff, dt_eff, gamma_n, &beta, NULL, NULL);

    const double div_u =
        grad_u_q[0][0] + grad_u_q[1][1] + grad_u_q[2][2];

    double grad_u_n[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) grad_u_n[a] += grad_u_q[a][b] * n[b];
    }

    double traction[3];
    for (int a = 0; a < 3; ++a) {
        traction[a] = -p_q*n[a] + mu_eff*(grad_u_n[a] + n[a]*div_u);
    }

    const double gni = t4n_dot3(gradNi, n);
    const double gnj = t4n_dot3(gradNj, n);

    for (int a = 0; a < 3; ++a) {
        vec_i[a] += dt * (-Ni * traction[a]);
        vec_i[a] += dt * (-mu_eff * (gni*r[a] + gradNi[a]*rn));
        vec_i[a] += dt * ( Ni * beta * r[a]);
    }
    vec_i[3] += dt * Ni * rn;

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double delta_ab = (a == b ? 1.0 : 0.0);
            const double dtraction_du =
                mu_eff * (delta_ab*gnj + n[a]*gradNj[b]);
            mat_ij[a][b] += dt * (-Ni * dtraction_du);
            mat_ij[a][b] += dt * (-mu_eff *
                (gni*Nj*delta_ab + gradNi[a]*Nj*n[b]));
            mat_ij[a][b] += dt * (Ni*Nj*beta*delta_ab);
        }
        mat_ij[a][3] += dt * (Ni*Nj*n[a]);
    }
    for (int b = 0; b < 3; ++b) {
        mat_ij[3][b] += dt * (Ni*Nj*n[b]);
    }
}

static void t4n_wall_normal_nitsche_local_vecmat(
    double vec_i[4],
    double mat_ij[4][4],
    double Ni,
    double Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    double p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double gamma_n)
{
    const double eps = 1.0e-30;
    const double dt_eff = fmax(dt, eps);
    const double h_eff = fmax(h_n, 1.0e-12);

    const double r[3] = {
        u_q[0]-Uw[0], u_q[1]-Uw[1], u_q[2]-Uw[2]
    };
    const double rn = t4n_dot3(r, n);

    double beta = 0.0;
    t4n_get_penalty_coefficients(
        rho, mu_eff, h_eff, dt_eff, gamma_n, &beta, NULL, NULL);

    const double div_u =
        grad_u_q[0][0] + grad_u_q[1][1] + grad_u_q[2][2];
    double grad_u_n[3] = {0.0,0.0,0.0};
    for (int a=0;a<3;++a) for (int b=0;b<3;++b) grad_u_n[a] += grad_u_q[a][b]*n[b];
    const double traction_n =
        -p_q + mu_eff*(t4n_dot3(grad_u_n,n) + div_u);

    const double gni = t4n_dot3(gradNi,n);
    const double gnj = t4n_dot3(gradNj,n);

    for (int a=0;a<3;++a) {
        const double wn = Ni*n[a];
        const double twn = mu_eff*(n[a]*gni + gradNi[a]);
        vec_i[a] += dt*(-wn*traction_n);
        vec_i[a] += dt*(-twn*rn);
        vec_i[a] += dt*( beta*wn*rn);
    }
    vec_i[3] += dt*Ni*rn;

    for (int a=0;a<3;++a) {
        for (int b=0;b<3;++b) {
            const double dtn = mu_eff*(n[b]*gnj + gradNj[b]);
            mat_ij[a][b] += dt*(-Ni*n[a]*dtn);
            mat_ij[a][b] += dt*(-mu_eff*(n[a]*gni + gradNi[a])*Nj*n[b]);
            mat_ij[a][b] += dt*( beta*Ni*n[a]*Nj*n[b]);
        }
        mat_ij[a][3] += dt*(Ni*Nj*n[a]);
    }
    for (int b=0;b<3;++b) mat_ij[3][b] += dt*(Ni*Nj*n[b]);
}

static void t4n_add_spalding_walllaw_local_vecmat(
    double vec_i[4],
    double mat_ij[4][4],
    double Ni,
    double exchange_Nj,
    const double n[3],
    const double u_exchange[3],
    const double Uw[3],
    double rho,
    double mu_wall_law,
    double y_exchange,
    double dt)
{
    double ut[3] = {0.0,0.0,0.0};
    double beta_wall = 0.0;
    t4n_spalding_wall_coefficients(
        u_exchange,
        Uw,
        n,
        rho,
        mu_wall_law,
        y_exchange,
        ut,
        NULL,
        NULL,
        &beta_wall);

    for (int a=0;a<3;++a) {
        vec_i[a] += dt * Ni * beta_wall * ut[a];
    }

    /* beta_wall is frozen.  For opposite-node exchange, exchange_Nj is 1 only
     * for the owner-TET node opposite the wall face; for legacy wall-q exchange
     * it is the usual face-interpolated Nj. */
    for (int a=0;a<3;++a) {
        for (int b=0;b<3;++b) {
            const double Pt_ab = (a==b ? 1.0 : 0.0) - n[a]*n[b];
            mat_ij[a][b] += dt * Ni * beta_wall * Pt_ab * exchange_Nj;
        }
    }
}

static void t4n_wall_local_vecmat(
    double vec_i[4],
    double mat_ij[4][4],
    double Ni,
    double Nj,
    double exchange_Nj,
    const double gradNi[3],
    const double gradNj[3],
    const double n[3],
    const double u_q[3],
    const double u_exchange[3],
    double p_q,
    const double grad_u_q[3][3],
    const double Uw[3],
    double rho,
    double mu_eff,
    double mu_wall_law,
    double h_n,
    double y_exchange,
    double dt,
    double gamma_n)
{
    for (int a=0;a<4;++a) {
        vec_i[a]=0.0;
        for (int b=0;b<4;++b) mat_ij[a][b]=0.0;
    }

    switch (g_t4_wall_mode) {
    case T4_WALL_STRONG_NOSLIP:
        /* Component-wise strong Dirichlet is applied by main_FOM. */
        return;

    case T4_WALL_NITSCHE_NOSLIP:
        t4n_wall_fullvec_nitsche_local_vecmat(
            vec_i, mat_ij, Ni, Nj, gradNi, gradNj, n, u_q, p_q,
            grad_u_q, Uw, rho, mu_eff, h_n, dt, gamma_n);
        return;

    case T4_WALL_NITSCHE_FREESLIP:
        t4n_wall_normal_nitsche_local_vecmat(
            vec_i, mat_ij, Ni, Nj, gradNi, gradNj, n, u_q, p_q,
            grad_u_q, Uw, rho, mu_eff, h_n, dt, gamma_n);
        return;

    case T4_WALL_NITSCHE_SPALDING:
        t4n_wall_normal_nitsche_local_vecmat(
            vec_i, mat_ij, Ni, Nj, gradNi, gradNj, n, u_q, p_q,
            grad_u_q, Uw, rho, mu_eff, h_n, dt, gamma_n);
        t4n_add_spalding_walllaw_local_vecmat(
            vec_i, mat_ij, Ni, exchange_Nj, n, u_exchange, Uw,
            rho, mu_wall_law, y_exchange, dt);
        return;

    case T4_WALL_STRONG_NORMAL_FREESLIP:
        /* Projected strong normal: no additional surface traction. */
        return;

    case T4_WALL_STRONG_NORMAL_SPALDING:
        /* Projected strong normal + modeled tangential traction only. */
        t4n_add_spalding_walllaw_local_vecmat(
            vec_i, mat_ij, Ni, exchange_Nj, n, u_exchange, Uw,
            rho, mu_wall_law, y_exchange, dt);
        return;

    default:
        fprintf(stderr, "%s ERROR: unsupported TET4 wall mode=%d\n", CODENAME, g_t4_wall_mode);
        exit(EXIT_FAILURE);
    }
}

/*
 * TET4 の体積 V と、その境界 TRI3 の面積 A から
 * 法線方向要素長 h_n = 3V/A を求める。
 *
 * TET4: |detJ| = 6V
 * TRI3: Jface_raw = |dx/dr x dx/ds| = 2A
 * よって h_n = |detJ| / Jface_raw = 3V/A
 */
static double t4_tet_face_height_from_jacobians(
    const double tet_detJ,
    const double Jface_raw)
{
    const double eps = 1.0e-30;
    const double detJ_abs = fabs(tet_detJ);
    const double Jface_abs = fabs(Jface_raw);

    if (detJ_abs <= eps || Jface_abs <= eps) {
        return 1.0e-12;
    }

    return fmax(detJ_abs / Jface_abs, 1.0e-12);
}


/*
 * 必要なヘッダ。
 * 既に読み込まれているものは重複して追加しなくてよい。
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define T4N_OWNER_NO_LOCAL_TET (-1)

static void t4n_fatal_exit(
    const char* format,
    ...)
{
    const int rank =
        monolis_mpi_get_global_my_rank();

    fprintf(
        stderr,
        "[rank %d][wall-t4] FATAL: ",
        rank
    );

    va_list args;

    va_start(
        args,
        format
    );

    vfprintf(
        stderr,
        format,
        args
    );

    va_end(args);

    fputc(
        '\n',
        stderr
    );

    fflush(stderr);
    fflush(stdout);

    exit(EXIT_FAILURE);
}



typedef struct {
    long long surf_total;
    long long owner_mapped;
    long long skipped_no_local_owner;
    long long repeated_surface_copies;

    long long face_assembled;
    long long qp_assembled;

    long long rhs_entries;
    long long matrix_entries;

    long long invalid_owner;
    long long invalid_surface_map;
    long long invalid_tet;
    long long invalid_face;
    long long invalid_coefficient;

    double area_sum;

    double h_n_min;
    double h_n_max;

    double Jface_min;
    double Jface_max;

    double mu_eff_min;
    double mu_eff_max;

    double beta_min;
    double beta_max;

    double normal_dot_min;

    double max_wall_residual;

    double rhs_abs_sum;
    double matrix_abs_sum;

    int max_residual_surface_elem;
    int max_residual_owner_tet;
    int max_residual_qp;

    double max_residual_position[3];
    double max_residual_velocity[3];
    double max_residual_normal[3];
} T4NitscheAssemblyDebug;


static void t4n_init_assembly_debug(
    T4NitscheAssemblyDebug* dbg,
    const int surf_total)
{
    if (dbg == NULL) {
        t4n_fatal_exit(
            "NULL debug data in t4n_init_assembly_debug"
        );
    }

    memset(
        dbg,
        0,
        sizeof(*dbg)
    );

    dbg->surf_total =
        (long long)surf_total;

    dbg->h_n_min =
        HUGE_VAL;

    dbg->Jface_min =
        HUGE_VAL;

    dbg->mu_eff_min =
        HUGE_VAL;

    dbg->beta_min =
        HUGE_VAL;

    dbg->normal_dot_min =
        HUGE_VAL;

    dbg->max_residual_surface_elem =
        -1;

    dbg->max_residual_owner_tet =
        -1;

    dbg->max_residual_qp =
        -1;
}


static void t4n_write_assembly_debug(
    const T4NitscheAssemblyDebug* dbg,
    const long long call,
    const int num_integ_points,
    const double gamma_n,
    const double dt)
{
    if (dbg == NULL) {
        t4n_fatal_exit(
            "NULL debug data in t4n_write_assembly_debug"
        );
    }

    const int rank =
        monolis_mpi_get_global_my_rank();

    char filename[256];

    snprintf(
        filename,
        sizeof(filename),
        "nitsche_assembly_rank_%06d.txt",
        rank
    );

    FILE* fp =
        fopen(
            filename,
            call == 0 ? "w" : "a"
        );

    if (fp == NULL) {
        t4n_fatal_exit(
            "cannot open assembly debug file: %s",
            filename
        );
    }

    if (call == 0) {
        fprintf(
            fp,
            "# call rank "
            "surf_total owner_mapped skipped_no_local_owner "
            "repeated_surface_copies "
            "face_assembled qp_assembled qp_expected "
            "rhs_entries matrix_entries "
            "invalid_owner invalid_surface_map invalid_tet "
            "invalid_face invalid_coefficient "
            "area_sum "
            "h_n_min h_n_max "
            "Jface_min Jface_max "
            "mu_eff_min mu_eff_max "
            "beta_min beta_max "
            "normal_dot_min max_abs_normal_residual "
            "rhs_abs_sum matrix_abs_sum "
            "max_es max_owner max_qp "
            "max_x max_y max_z "
            "max_ux max_uy max_uz "
            "max_nx max_ny max_nz "
            "gamma_n dt\n"
        );
    }

    const long long qp_expected =
        dbg->owner_mapped
        * (long long)num_integ_points;

    const double h_n_min =
        isfinite(dbg->h_n_min)
        ? dbg->h_n_min
        : 0.0;

    const double Jface_min =
        isfinite(dbg->Jface_min)
        ? dbg->Jface_min
        : 0.0;

    const double mu_eff_min =
        isfinite(dbg->mu_eff_min)
        ? dbg->mu_eff_min
        : 0.0;

    const double beta_min =
        isfinite(dbg->beta_min)
        ? dbg->beta_min
        : 0.0;

    const double normal_dot_min =
        isfinite(dbg->normal_dot_min)
        ? dbg->normal_dot_min
        : 0.0;

    fprintf(
        fp,
        "%lld %d "
        "%lld %lld %lld %lld "
        "%lld %lld %lld "
        "%lld %lld "
        "%lld %lld %lld %lld %lld "
        "%.15e "
        "%.15e %.15e "
        "%.15e %.15e "
        "%.15e %.15e "
        "%.15e %.15e "
        "%.15e %.15e "
        "%.15e %.15e "
        "%d %d %d "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e\n",

        call,
        rank,

        dbg->surf_total,
        dbg->owner_mapped,
        dbg->skipped_no_local_owner,
        dbg->repeated_surface_copies,

        dbg->face_assembled,
        dbg->qp_assembled,
        qp_expected,

        dbg->rhs_entries,
        dbg->matrix_entries,

        dbg->invalid_owner,
        dbg->invalid_surface_map,
        dbg->invalid_tet,
        dbg->invalid_face,
        dbg->invalid_coefficient,

        dbg->area_sum,

        h_n_min,
        dbg->h_n_max,

        Jface_min,
        dbg->Jface_max,

        mu_eff_min,
        dbg->mu_eff_max,

        beta_min,
        dbg->beta_max,

        normal_dot_min,
        dbg->max_wall_residual,

        dbg->rhs_abs_sum,
        dbg->matrix_abs_sum,

        dbg->max_residual_surface_elem,
        dbg->max_residual_owner_tet,
        dbg->max_residual_qp,

        dbg->max_residual_position[0],
        dbg->max_residual_position[1],
        dbg->max_residual_position[2],

        dbg->max_residual_velocity[0],
        dbg->max_residual_velocity[1],
        dbg->max_residual_velocity[2],

        dbg->max_residual_normal[0],
        dbg->max_residual_normal[1],
        dbg->max_residual_normal[2],

        gamma_n,
        dt
    );

    fclose(fp);
}


void set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    (void)basis_vol;

    static long long call_count = 0;

    const long long call =
        call_count++;

    const int rank =
        monolis_mpi_get_global_my_rank();

    int* owner = NULL;
    int* lface = NULL;

    double* yplus_weight = NULL;

/*
 * This source is compiled as C++, where a forward goto may not skip
 * initialized automatic variables.  Failure is fatal anyway, so release
 * the three function-owned buffers here and terminate directly.
 */
#define T4N_FAIL(fmt, ...) do {                                      \
    fprintf(                                                         \
        stderr,                                                      \
        "[rank %d][wall-t4] " fmt "\n",                              \
        rank,                                                        \
        __VA_ARGS__);                                                \
    fflush(stderr);                                                  \
    free(yplus_weight);                                              \
    free(owner);                                                     \
    free(lface);                                                     \
    exit(EXIT_FAILURE);                                              \
} while (0)

    if (
        monolis == NULL ||
        fe == NULL ||
        surf == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        Uw == NULL
    ) {
        fprintf(
            stderr,
            "[rank %d][wall-t4] NULL input\n",
            rank);

        fflush(stderr);

        exit(EXIT_FAILURE);
    }

    /*
     * 表面要素を持たないrankは正常なno-op。
     */
    if (surf->total_num_elems == 0) {
        return;
    }

    if (
        fe->local_num_nodes != 4 ||
        surf->local_num_nodes != 3
    ) {
        T4N_FAIL(
            "TET4/TRI3 expected: "
            "fe_nnode=%d surf_nnode=%d nsurf=%d",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);
    }

    const int ns = 3;
    const int nv = 4;

    const int np =
        basis_surf->num_integ_points;

    const double area_factor =
        t4n_tri_area_factor(
            basis_surf);

    if (
        np <= 0 ||
        !isfinite(area_factor) ||
        area_factor <= 0.0
    ) {
        T4N_FAIL(
            "invalid TRI3 quadrature: "
            "np=%d area_factor=%.15e",
            np,
            area_factor);
    }

    /*
     * y+節点投影。
     */
    int do_yplus =
        (
            vals->y_plus != NULL &&
            vals->y_plus_count != NULL
        );

    if (do_yplus) {
        for (
            int i = 0;
            i < fe->total_num_nodes;
            ++i
        ) {
            vals->y_plus[i] = 0.0;
            vals->y_plus_count[i] = 0;
        }

        yplus_weight =
            (double*)calloc(
                (size_t)fe->total_num_nodes,
                sizeof(double));

        if (yplus_weight == NULL) {
            fprintf(
                stderr,
                "[rank %d][wall-t4] WARNING: "
                "yPlus output disabled\n",
                rank);

            fflush(stderr);

            do_yplus = 0;
        }
    }

    /*
     * ----------------------------------------------------------------------
     * 修正3:
     *
     * 厳格版t4n_make_owner_map()ではなく、
     * ownerなし表面コピーをowner[es]=-1のまま返す
     * t4n_make_owner_map_overlap()を使用する。
     * ----------------------------------------------------------------------
     */

    int num_unmapped_local = 0;

    if (
        t4n_make_owner_map_overlap(
            surf,
            fe,
            &owner,
            &lface,
            &num_unmapped_local) != 0
    ) {
        T4N_FAIL(
            "overlap owner-map construction failed: "
            "nsurf=%d",
            surf->total_num_elems);
    }

    if (
        owner == NULL ||
        lface == NULL
    ) {
        T4N_FAIL(
            "overlap owner-map returned NULL arrays: "
            "nsurf=%d",
            surf->total_num_elems);
    }

    /*
     * Newton反復ごとに同じ警告を出さないため、
     * 最初の呼び出しだけ表示する。
     */
    if (
        call == 0 &&
        num_unmapped_local > 0
    ) {
        fprintf(
            stderr,
            "[rank %d][wall-t4] overlap owner map: "
            "surf=%d mapped=%d no_local_owner=%d\n",
            rank,
            surf->total_num_elems,
            surf->total_num_elems
                - num_unmapped_local,
            num_unmapped_local);

        fflush(stderr);
    }

    /*
     * owner TET番号を使うため、
     * この関数内ではlfaceは参照しない。
     */
    (void)lface;

    /*
     * ----------------------------------------------------------------------
     * Nitsche組立て診断値
     * ----------------------------------------------------------------------
     */

    long long mapped_faces = 0;
    long long skipped_no_owner = 0;

    long long assembled_faces = 0;
    long long assembled_qp = 0;

    long long rhs_entries = 0;
    long long matrix_entries = 0;

    double area_sum = 0.0;

    double rhs_abs_sum = 0.0;
    double matrix_abs_sum = 0.0;

    double h_n_min = HUGE_VAL;
    double h_n_max = 0.0;

    double beta_min = HUGE_VAL;
    double beta_max = 0.0;

    double max_wall_residual = 0.0;

    int max_es = -1;
    int max_owner = -1;
    int max_qp = -1;

    double max_x[3] = {
        0.0,
        0.0,
        0.0
    };

    double max_u[3] = {
        0.0,
        0.0,
        0.0
    };

    double max_n[3] = {
        0.0,
        0.0,
        0.0
    };

    /*
     * ----------------------------------------------------------------------
     * 修正4:
     *
     * ownerのない表面コピーだけをスキップする。
     * 他の面のNitsche組立ては継続する。
     * ----------------------------------------------------------------------
     */

    for (
        int es = 0;
        es < surf->total_num_elems;
        ++es
    ) {
        const int ke =
            owner[es];

        /*
         * このrankは表面コピーを持っているが、
         * 対応owner TETを持っていない。
         *
         * この面だけをスキップする。
         */
        if (ke < 0) {
            ++skipped_no_owner;
            continue;
        }

        if (ke >= fe->total_num_elems) {
            T4N_FAIL(
                "invalid mapped owner: "
                "es=%d owner=%d ne=%d",
                es,
                ke,
                fe->total_num_elems);
        }

        /*
         * TRI3とowner TET4が3節点を共有することを確認。
         */
        int common = 0;

        for (int a = 0; a < ns; ++a) {
            const int sgid =
                surf->conn[es][a];

            for (int b = 0; b < nv; ++b) {
                if (
                    fe->conn[ke][b]
                    == sgid
                ) {
                    ++common;
                    break;
                }
            }
        }

        if (common != ns) {
            T4N_FAIL(
                "surface-owner mismatch: "
                "es=%d owner=%d common=%d/%d",
                es,
                ke,
                common,
                ns);
        }

        /*
         * TRI3節点座標。
         */
        double xsurf[3][3];

        /*
         * owner TET4節点座標。
         */
        double xvol[4][3];

        for (int a = 0; a < ns; ++a) {
            const int gid =
                surf->conn[es][a];

            if (
                gid < 0 ||
                gid >= fe->total_num_nodes
            ) {
                T4N_FAIL(
                    "invalid surface node: "
                    "es=%d a=%d gid=%d nn=%d",
                    es,
                    a,
                    gid,
                    fe->total_num_nodes);
            }

            for (int c = 0; c < 3; ++c) {
                xsurf[a][c] =
                    fe->x[gid][c];
            }
        }

        for (int a = 0; a < nv; ++a) {
            const int gid =
                fe->conn[ke][a];

            if (
                gid < 0 ||
                gid >= fe->total_num_nodes
            ) {
                T4N_FAIL(
                    "invalid volume node: "
                    "owner=%d a=%d gid=%d nn=%d",
                    ke,
                    a,
                    gid,
                    fe->total_num_nodes);
            }

            for (int c = 0; c < 3; ++c) {
                xvol[a][c] =
                    fe->x[gid][c];
            }
        }

        /*
         * TRI3局所節点からTET4局所節点への対応。
         */
        int s2v[3];

        if (
            !t4n_surface_to_volume_nodes(
                surf,
                fe,
                es,
                ke,
                s2v)
        ) {
            T4N_FAIL(
                "surface-to-volume mapping failed: "
                "es=%d owner=%d",
                es,
                ke);
        }

        /* For a TET4 face the local node opposite the TRI3 is the one not
         * present in s2v[].  0+1+2+3 = 6. */
        const int opposite_local_node =
            6 - s2v[0] - s2v[1] - s2v[2];

        if (opposite_local_node < 0 || opposite_local_node >= nv) {
            T4N_FAIL(
                "invalid opposite local node: es=%d owner=%d opp=%d",
                es, ke, opposite_local_node);
        }

        const int opposite_gid =
            fe->conn[ke][opposite_local_node];

        /*
         * TET4 metricと物理勾配。
         */
        double dN_dx[4][3];
        double J_inv_face[3][3];

        double Jacobian_face = 0.0;

        if (
            !t4n_grad_phys_metric(
                xvol,
                dN_dx,
                J_inv_face,
                &Jacobian_face)
        ) {
            T4N_FAIL(
                "degenerate TET4: "
                "es=%d owner=%d",
                es,
                ke);
        }

        if (
            !isfinite(Jacobian_face) ||
            fabs(Jacobian_face)
                <= 1.0e-300
        ) {
            T4N_FAIL(
                "invalid TET Jacobian: "
                "es=%d owner=%d detJ=%.15e",
                es,
                ke,
                Jacobian_face);
        }

        /*
         * TET4物理体積は|detJ|/6。
         */
        const double h_e_vms =
            cbrt(
                fabs(Jacobian_face)
                / 6.0);

        if (
            !isfinite(h_e_vms) ||
            h_e_vms <= 1.0e-12
        ) {
            T4N_FAIL(
                "invalid VMS length: "
                "es=%d owner=%d h_e=%.15e",
                es,
                ke,
                h_e_vms);
        }

        /*
         * owner TET4重心。
         */
        double xc[3] = {
            0.0,
            0.0,
            0.0
        };

        /*
         * TRI3面重心。
         */
        double xf[3] = {
            0.0,
            0.0,
            0.0
        };

        for (int a = 0; a < nv; ++a) {
            for (int c = 0; c < 3; ++c) {
                xc[c] +=
                    xvol[a][c];
            }
        }

        for (int c = 0; c < 3; ++c) {
            xc[c] /= 4.0;
        }

        for (int a = 0; a < ns; ++a) {
            for (int c = 0; c < 3; ++c) {
                xf[c] +=
                    xsurf[a][c];
            }
        }

        for (int c = 0; c < 3; ++c) {
            xf[c] /= 3.0;
        }

        /*
         * owner重心から壁面重心へ向かうベクトル。
         */
        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        int face_qp_count = 0;

        for (int p = 0; p < np; ++p) {
            double nrm[3];

            t4n_tri3_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nrm);

            const double Jface_raw =
                t4n_norm3(
                    nrm);

            if (
                !isfinite(Jface_raw) ||
                Jface_raw <= 1.0e-300
            ) {
                T4N_FAIL(
                    "degenerate TRI3: "
                    "es=%d owner=%d qp=%d "
                    "Jface=%.15e",
                    es,
                    ke,
                    p,
                    Jface_raw);
            }

            double n[3] = {
                nrm[0] / Jface_raw,
                nrm[1] / Jface_raw,
                nrm[2] / Jface_raw
            };

            /*
             * 法線をowner重心から面へ向かう方向へ統一。
             */
            if (
                t4n_dot3(
                    n,
                    svec) < 0.0
            ) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            const double normal_dot =
                t4n_dot3(
                    n,
                    svec);

            if (
                !isfinite(normal_dot) ||
                normal_dot <= 0.0
            ) {
                T4N_FAIL(
                    "invalid oriented normal: "
                    "es=%d owner=%d qp=%d "
                    "n_dot_s=%.15e",
                    es,
                    ke,
                    p,
                    normal_dot);
            }

            /*
             * TRI3求積点上のTET4形状関数値。
             */
            double Nv[4] = {
                0.0,
                0.0,
                0.0,
                0.0
            };

            for (int a = 0; a < ns; ++a) {
                if (
                    s2v[a] < 0 ||
                    s2v[a] >= nv
                ) {
                    T4N_FAIL(
                        "invalid s2v: "
                        "es=%d owner=%d "
                        "a=%d s2v=%d",
                        es,
                        ke,
                        a,
                        s2v[a]);
                }

                Nv[s2v[a]] =
                    basis_surf->N[p][a];
            }

            double u_q[3] = {
                0.0,
                0.0,
                0.0
            };

            double u_old_q[3] = {
                0.0,
                0.0,
                0.0
            };

            double p_q = 0.0;

            double grad_u_q[3][3] = {
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0}
            };

            double grad_p_q[3] = {
                0.0,
                0.0,
                0.0
            };

            for (int a = 0; a < nv; ++a) {
                const int gid =
                    fe->conn[ke][a];

                for (int c = 0; c < 3; ++c) {
                    u_q[c] +=
                        Nv[a]
                        * vals->v[gid][c];

                    u_old_q[c] +=
                        Nv[a]
                        * vals->v_old[gid][c];

                    for (int d = 0; d < 3; ++d) {
                        grad_u_q[c][d] +=
                            dN_dx[a][d]
                            * vals->v[gid][c];
                    }
                }

                p_q +=
                    Nv[a]
                    * vals->p[gid];

                for (int d = 0; d < 3; ++d) {
                    grad_p_q[d] +=
                        dN_dx[a][d]
                        * vals->p[gid];
                }
            }

            double* grad_u_ptr[3] = {
                grad_u_q[0],
                grad_u_q[1],
                grad_u_q[2]
            };

            double mu_eff_face =
                mu_molecular;

            double tau_face = 0.0;
            double tau_c_face = 0.0;

            BBFE_vms_mu_eff_tau(
                &mu_eff_face,
                &tau_face,
                &tau_c_face,
                J_inv_face,
                Jacobian_face,
                h_e_vms,
                u_q,
                u_old_q,
                grad_u_ptr,
                grad_p_q,
                rho,
                mu_molecular,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);

            (void)tau_face;
            (void)tau_c_face;

            if (
                !isfinite(mu_eff_face) ||
                mu_eff_face <= 0.0
            ) {
                T4N_FAIL(
                    "invalid mu_eff: "
                    "es=%d owner=%d qp=%d "
                    "mu_eff=%.15e",
                    es,
                    ke,
                    p,
                    mu_eff_face);
            }

            /*
             * 現行solverのh_n定義を維持。
             */
            const double h_n =
                t4_tet_face_height_from_jacobians(
                    Jacobian_face,
                    Jface_raw);

            const double y_exchange =
                h_n;

            double u_exchange[3];
            if (g_t4_wall_exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE) {
                u_exchange[0] = vals->v[opposite_gid][0];
                u_exchange[1] = vals->v[opposite_gid][1];
                u_exchange[2] = vals->v[opposite_gid][2];
            }
            else {
                u_exchange[0] = u_q[0];
                u_exchange[1] = u_q[1];
                u_exchange[2] = u_q[2];
            }

            const double w =
                basis_surf->integ_weight[p]
                * area_factor
                * Jface_raw;

            if (
                !isfinite(h_n) ||
                h_n <= 0.0 ||
                !isfinite(w) ||
                w <= 0.0
            ) {
                T4N_FAIL(
                    "invalid wall coefficient: "
                    "es=%d owner=%d qp=%d "
                    "h_n=%.15e w=%.15e",
                    es,
                    ke,
                    p,
                    h_n,
                    w);
            }

            /*
             * Nitsche penaltyを診断用に再計算。
             */
            const double h_eff =
                fmax(
                    h_n,
                    1.0e-12);

            const double dt_eff =
                fmax(
                    vals->dt,
                    1.0e-30);

            double beta = 0.0;

            t4n_get_penalty_coefficients(
                rho,
                mu_eff_face,
                h_eff,
                dt_eff,
                gamma_n,
                &beta,
                NULL,
                NULL);

            if (
                !isfinite(beta) ||
                beta <= 0.0
            ) {
                T4N_FAIL(
                    "invalid beta: "
                    "es=%d owner=%d qp=%d "
                    "beta=%.15e",
                    es,
                    ke,
                    p,
                    beta);
            }

            area_sum +=
                w;

            h_n_min =
                fmin(
                    h_n_min,
                    h_n);

            h_n_max =
                fmax(
                    h_n_max,
                    h_n);

            beta_min =
                fmin(
                    beta_min,
                    beta);

            beta_max =
                fmax(
                    beta_max,
                    beta);

            /*
             * 壁面残差 r=u-Uw。
             */
            const double r[3] = {
                u_q[0] - Uw[0],
                u_q[1] - Uw[1],
                u_q[2] - Uw[2]
            };

            const double r_norm =
                fabs(t4n_dot3(r, n));

            /*
             * 求積点座標。
             */
            double x_q[3] = {
                0.0,
                0.0,
                0.0
            };

            for (int a = 0; a < ns; ++a) {
                for (int c = 0; c < 3; ++c) {
                    x_q[c] +=
                        basis_surf->N[p][a]
                        * xsurf[a][c];
                }
            }

            if (
                r_norm >
                max_wall_residual
            ) {
                max_wall_residual =
                    r_norm;

                max_es =
                    es;

                max_owner =
                    ke;

                max_qp =
                    p;

                for (int c = 0; c < 3; ++c) {
                    max_x[c] =
                        x_q[c];

                    max_u[c] =
                        u_q[c];

                    max_n[c] =
                        n[c];
                }
            }

            /*
             * y+ nodal projection.  Spalding modes use the selected exchange
             * velocity (opposite-node by default), not the wall-face DOF.
             */
            if (do_yplus) {
                double ut_exchange_diag[3] = {0.0, 0.0, 0.0};
                double u_tau = 0.0;
                t4n_spalding_wall_coefficients(
                    u_exchange,
                    Uw,
                    n,
                    rho,
                    mu_molecular,
                    y_exchange,
                    ut_exchange_diag,
                    NULL,
                    &u_tau,
                    NULL);

                double yplus_q = 0.0;
                if (u_tau > 0.0 && rho > 0.0 && mu_molecular > 0.0) {
                    yplus_q = rho * u_tau * y_exchange / mu_molecular;
                }

                for (int a = 0; a < ns; ++a) {
                    const int gid = surf->conn[es][a];
                    const double nodal_weight = w * basis_surf->N[p][a];
                    vals->y_plus[gid] += nodal_weight * yplus_q;
                    yplus_weight[gid] += nodal_weight;
                    vals->y_plus_count[gid] += 1;
                }
            }

            /*
             * Nitsche残差・行列組立て。
             */
            for (int i = 0; i < nv; ++i) {
                const int gi =
                    fe->conn[ke][i];

                double vec_i[4];
                double dummy_mat[4][4];

                const double zero_grad[3] = {
                    0.0,
                    0.0,
                    0.0
                };

                t4n_wall_local_vecmat(
                    vec_i,
                    dummy_mat,
                    Nv[i],
                    0.0,
                    0.0,
                    dN_dx[i],
                    zero_grad,
                    n,
                    u_q,
                    u_exchange,
                    p_q,
                    grad_u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu_molecular,
                    h_n,
                    y_exchange,
                    vals->dt,
                    gamma_n);

                /*
                 * Nitsche残差。
                 */
                for (int a = 0; a < 4; ++a) {
                    const double value =
                        w
                        * vec_i[a];

                    monolis->mat.R.B[
                        4*gi + a
                    ] -= value;

                    rhs_abs_sum +=
                        fabs(value);

                    ++rhs_entries;
                }

                /*
                 * Nitsche行列。
                 */
                for (int j = 0; j < nv; ++j) {
                    const int gj =
                        fe->conn[ke][j];

                    double dummy_vec[4];
                    double mat_ij[4][4];

                    const double exchange_Nj =
                        (g_t4_wall_exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE)
                        ? (j == opposite_local_node ? 1.0 : 0.0)
                        : Nv[j];

                    t4n_wall_local_vecmat(
                        dummy_vec,
                        mat_ij,
                        Nv[i],
                        Nv[j],
                        exchange_Nj,
                        dN_dx[i],
                        dN_dx[j],
                        n,
                        u_q,
                        u_exchange,
                        p_q,
                        grad_u_q,
                        Uw,
                        rho,
                        mu_eff_face,
                        mu_molecular,
                        h_n,
                        y_exchange,
                        vals->dt,
                        gamma_n);

                    for (int a = 0; a < 4; ++a) {
                        for (int b = 0; b < 4; ++b) {
                            const double value =
                                w
                                * mat_ij[a][b];

                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,
                                gi,
                                gj,
                                a,
                                b,
                                value);

                            matrix_abs_sum +=
                                fabs(value);

                            ++matrix_entries;
                        }
                    }
                }
            }

            ++assembled_qp;
            ++face_qp_count;
        }

        if (face_qp_count != np) {
            T4N_FAIL(
                "incomplete face quadrature: "
                "es=%d owner=%d qp=%d/%d",
                es,
                ke,
                face_qp_count,
                np);
        }

        ++assembled_faces;
        ++mapped_faces;
    }

    /*
     * y+節点投影の正規化。
     */
    if (do_yplus) {
        for (
            int i = 0;
            i < fe->total_num_nodes;
            ++i
        ) {
            if (
                yplus_weight[i] >
                1.0e-30
            ) {
                vals->y_plus[i] /=
                    yplus_weight[i];
            } else {
                vals->y_plus[i] = 0.0;
                vals->y_plus_count[i] = 0;
            }
        }
    }

    const long long expected_qp =
        mapped_faces
        * (long long)np;

    /*
     * すべての面が、
     *
     * mappedまたはskipped
     *
     * のどちらかに分類されていることを確認。
     */
    if (
        mapped_faces
        + skipped_no_owner
            != (long long)surf->total_num_elems ||

        skipped_no_owner
            != (long long)num_unmapped_local ||

        assembled_faces
            != mapped_faces ||

        assembled_qp
            != expected_qp ||

        (
            mapped_faces > 0 &&
            (
                rhs_entries <= 0 ||
                matrix_entries <= 0 ||
                matrix_abs_sum <= 0.0
            )
        )
    ) {
        T4N_FAIL(
            "incomplete assembly: "
            "surf=%d mapped=%lld skipped=%lld "
            "face=%lld qp=%lld/%lld "
            "rhs=%lld mat=%lld "
            "mat_abs=%.15e",
            surf->total_num_elems,
            mapped_faces,
            skipped_no_owner,
            assembled_faces,
            assembled_qp,
            expected_qp,
            rhs_entries,
            matrix_entries,
            matrix_abs_sum);
    }

    /*
     * rank別診断ファイル。
     */
    {
        char filename[256];

        snprintf(
            filename,
            sizeof(filename),
            "nitsche_assembly_rank_%06d.txt",
            rank);

        FILE* fp =
            fopen(
                filename,
                call == 0
                    ? "w"
                    : "a");

        if (fp == NULL) {
            T4N_FAIL(
                "cannot open diagnostic file: %s",
                filename);
        }

        if (call == 0) {
            fprintf(
                fp,
                "# call rank "
                "surf mapped skipped "
                "face qp qp_expected "
                "rhs_entries matrix_entries "
                "area_sum "
                "h_min h_max "
                "beta_min beta_max "
                "max_abs_normal_residual "
                "rhs_abs matrix_abs "
                "max_es max_owner max_qp "
                "max_x max_y max_z "
                "max_ux max_uy max_uz "
                "max_nx max_ny max_nz "
                "gamma_n dt\n");
        }

        fprintf(
            fp,
            "%lld %d "
            "%d %lld %lld "
            "%lld %lld %lld "
            "%lld %lld "
            "%.15e "
            "%.15e %.15e "
            "%.15e %.15e "
            "%.15e "
            "%.15e %.15e "
            "%d %d %d "
            "%.15e %.15e %.15e "
            "%.15e %.15e %.15e "
            "%.15e %.15e %.15e "
            "%.15e %.15e\n",

            call,
            rank,

            surf->total_num_elems,
            mapped_faces,
            skipped_no_owner,

            assembled_faces,
            assembled_qp,
            expected_qp,

            rhs_entries,
            matrix_entries,

            area_sum,

            isfinite(h_n_min)
                ? h_n_min
                : 0.0,

            h_n_max,

            isfinite(beta_min)
                ? beta_min
                : 0.0,

            beta_max,

            max_wall_residual,

            rhs_abs_sum,
            matrix_abs_sum,

            max_es,
            max_owner,
            max_qp,

            max_x[0],
            max_x[1],
            max_x[2],

            max_u[0],
            max_u[1],
            max_u[2],

            max_n[0],
            max_n[1],
            max_n[2],

            gamma_n,
            vals->dt);

        fclose(fp);
    }

    free(yplus_weight);
    free(owner);
    free(lface);

#undef T4N_FAIL
}

/* -------------------------------------------------------------------------- */
/* HEX/TET dispatch wrapper. HEX body remains unchanged.                      */
/* -------------------------------------------------------------------------- */
void set_wall_face_vecmat_symmetric_nitsche_hex_or_tet(
    MONOLIS* monolis,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_vol,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    double gamma_n)
{
    if (
        monolis == NULL ||
        fe == NULL ||
        basis_vol == NULL ||
        surf == NULL ||
        basis_surf == NULL ||
        vals == NULL) {

        fprintf(
            stderr,
            "[wall] NULL input in "
            "set_wall_face_vecmat_symmetric_nitsche_hex_or_tet\n");

        return;
    }

    /*
     * 壁面要素を所有しないrank。
     * 局所行列・局所残差への寄与は0なので正常終了する。
     */
    if (surf->total_num_elems == 0) {
        return;
    }

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[wall] invalid number of surface elements: %d\n",
            surf->total_num_elems);

        return;
    }

    if (fe->local_num_nodes == 8) {
        if (surf->local_num_nodes != 4) {
            fprintf(
                stderr,
                "[wall] HEX8 volume requires QUAD4 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return;
        }

        set_wall_face_vecmat_symmetric_nitsche(
            monolis,
            fe,
            basis_vol,
            surf,
            basis_surf,
            vals,
            rho,
            mu_molecular,
            Uw,
            gamma_n);

        return;
    }

    if (fe->local_num_nodes == 4) {
        if (surf->local_num_nodes != 3) {
            fprintf(
                stderr,
                "[wall] TET4 volume requires TRI3 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return;
        }

        /* Strong no-slip is handled by the component-wise Dirichlet wrapper.
         * Projected-strong free-slip has no tangential traction. */
        if (g_t4_wall_mode == T4_WALL_STRONG_NOSLIP ||
            g_t4_wall_mode == T4_WALL_STRONG_NORMAL_FREESLIP) {
            return;
        }

        set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
            monolis,
            fe,
            basis_vol,
            surf,
            basis_surf,
            vals,
            rho,
            mu_molecular,
            Uw,
            gamma_n);

        return;
    }

    fprintf(
        stderr,
        "[wall] unsupported fe->local_num_nodes=%d. "
        "Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);
}

void calc_Cd_p_hex_or_tet(
    BBFE_DATA* surf,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out)
{
    /*
     * 空rankの局所寄与は0。
     */
    if (D_out != NULL) {
        *D_out = 0.0;
    }

    if (L_out != NULL) {
        *L_out = 0.0;
    }

    if (
        surf == NULL ||
        fe == NULL ||
        basis == NULL ||
        vals == NULL ||
        eU_in == NULL ||
        eP_in == NULL) {

        fprintf(
            stderr,
            "[Cd_p] NULL input\n");

        return;
    }

    /*
     * 壁面要素を持たないrankは正常。
     * main側のMPI_SUMには0を渡す。
     */
    if (surf->total_num_elems == 0) {
        return;
    }

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[Cd_p] invalid number of surface elements: %d\n",
            surf->total_num_elems);

        return;
    }

    if (fe->local_num_nodes == 8) {
        if (surf->local_num_nodes != 4) {
            fprintf(
                stderr,
                "[Cd_p] HEX8 volume requires QUAD4 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return;
        }

        calc_Cd_p(
            surf,
            fe,
            basis,
            vals,
            eU_in,
            eP_in,
            D_out,
            L_out);

        return;
    }

    if (fe->local_num_nodes == 4) {
        if (surf->local_num_nodes != 3) {
            fprintf(
                stderr,
                "[Cd_p] TET4 volume requires TRI3 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return;
        }

        calc_Cd_p_tet4_tri3(
            surf,
            fe,
            basis,
            vals,
            eU_in,
            eP_in,
            D_out,
            L_out);

        return;
    }

    fprintf(
        stderr,
        "[Cd_p] unsupported fe->local_num_nodes=%d. "
        "Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);
}

/*
#define T4_WALL_DIAG_NUM_REGIONS 8

enum {
    T4_WALL_DIAG_FRONT_POS_DRAG = 0,
    T4_WALL_DIAG_REAR_NEG_DRAG  = 1,
    T4_WALL_DIAG_TOP_POS_LIFT   = 2,
    T4_WALL_DIAG_BOTTOM_NEG_LIFT= 3,
    T4_WALL_DIAG_SIDE_POS       = 4,
    T4_WALL_DIAG_SIDE_NEG       = 5,
    T4_WALL_DIAG_SLANT_MIXED    = 6,
    T4_WALL_DIAG_OTHER          = 7
};

static const char* t4_wall_diag_region_name[T4_WALL_DIAG_NUM_REGIONS] = {
    "front_pos_drag",
    "rear_neg_drag",
    "top_pos_lift",
    "bottom_neg_lift",
    "side_pos",
    "side_neg",
    "slant_mixed",
    "other"
};

typedef struct {
    double area;
    double int_n[3];

    double max_abs_un_over_Uref;
    double int_abs_un;
    double int_un;
    double int_un2;
    double mean_abs_un_over_Uref;
    double mean_un_over_Uref;
    double rms_un_over_Uref;
    double Uref_scale;

    double p_min;
    double p_max;
    double p_mean;
    double p_rms;
    double p_area_sum;
    double p2_area_sum;

    double cp_min;
    double cp_max;
    double cp_mean;
    double cp_rms;
    double cp_area_sum;
    double cp2_area_sum;

    // Pressure force in the drag direction: int_Gamma p (n.eU) dGamma.
    //   If eU=(1,0,0), this is int_Gamma p n_x dGamma.
    double int_pnU;
    double area_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_pnU_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_abs_un_region[T4_WALL_DIAG_NUM_REGIONS];
} T4WallDiagnostics;
*/
static double t4wd_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double t4wd_norm3(const double a[3])
{
    return sqrt(t4wd_dot3(a, a));
}

static void t4wd_cross3(double c[3], const double a[3], const double b[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static void t4wd_output_gamma_value(
    const char* prefix,
    const char* component,
    int gamma_index,
    double value,
    double t,
    const char* directory)
{
    char fname[256];

    snprintf(
        fname,
        sizeof(fname),
        "wall_nitsche_%s_%s_gamma_%02d.txt",
        prefix,
        component,
        gamma_index);

    ROM_std_hlpod_output_add_calc_time(value, t, fname, directory);
}

/* -------------------------------------------------------------------------- */
/* Normal-only Nitsche penalty force density                                  */
/*                                                                            */
/* For Case C normal Nitsche + tangential Spalding wall law:                          */
/*                                                                            */
/*     f_body_penalty = - beta ((u-Uw).n) n                                 */
/*                                                                            */
/* This diagnostic corresponds only to the penalty part of the discrete        */
/* Nitsche boundary residual.  It does not include the consistency terms       */
/*     -w.t(u,p)  and  -t(w,q).(u-Uw).                                        */
/* -------------------------------------------------------------------------- */


static void t4_cd_normal_nitsche_penalty_force_density(
    const double u_q[3],
    const double Uw_in[3],
    const double n[3],
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double gamma_n,
    double f_body[3],
    double* beta_out,
    double* beta_mu_out,
    double* beta_dt_out)
{
    const double eps = 1.0e-30;
    const double dt_eff = fmax(dt, eps);
    const double h_eff  = fmax(h_n, 1.0e-12);

    const double Uw[3] = {
        Uw_in ? Uw_in[0] : 0.0,
        Uw_in ? Uw_in[1] : 0.0,
        Uw_in ? Uw_in[2] : 0.0
    };

    double r[3];
    for (int a = 0; a < 3; ++a) {
        r[a] = u_q[a] - Uw[a];
    }

    double beta = 0.0;
    double beta_mu = 0.0;
    double beta_dt = 0.0;

    t4n_get_penalty_coefficients(
        rho,
        mu_eff,
        h_eff,
        dt_eff,
        gamma_n,
        &beta,
        &beta_mu,
        &beta_dt);

    const double rn = t4n_dot3(r, n);

    for (int a = 0; a < 3; ++a) {
        if (g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP) {
            f_body[a] = -beta * r[a];
        }
        else if (g_t4_wall_mode == T4_WALL_NITSCHE_FREESLIP ||
                 g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING) {
            f_body[a] = -beta * rn * n[a];
        }
        else {
            f_body[a] = 0.0;
        }
    }

    if (beta_out)    *beta_out    = beta;
    if (beta_mu_out) *beta_mu_out = beta_mu;
    if (beta_dt_out) *beta_dt_out = beta_dt;
}

static int calc_Cd_v_tet4_tri3_core(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    int include_wall_nitsche,
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    (void)monolis_com;
    (void)Uinf;
    (void)Aref;

    if (D_out)           *D_out           = 0.0;
    if (L_out)           *L_out           = 0.0;
    if (Cd_out)          *Cd_out          = 0.0;
    if (Cl_out)          *Cl_out          = 0.0;
    if (D_nitsche_out)   *D_nitsche_out   = 0.0;
    if (L_nitsche_out)   *L_nitsche_out   = 0.0;
    if (Cd_nitsche_out)  *Cd_nitsche_out  = 0.0;
    if (Cl_nitsche_out)  *Cl_nitsche_out  = 0.0;

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        return -1;
    }

    double (*G_all)[3][3] =
        (double(*)[3][3])malloc(sizeof(double) * fe->total_num_elems * 9);
    if (!G_all) {
        free(owner);
        free(lface);
        return -1;
    }
    t4_compute_all_gradients_tet4_exact(fe, vals, G_all);

    const double area_factor = t4_tri_area_factor(basis_surf);

    double F_vis[3] = {0.0, 0.0, 0.0};
    double F_nit[3] = {0.0, 0.0, 0.0};

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xvol[4][3];
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        int s2v[3];
        if (!t4_surface_to_volume_nodes(surf, fe, es, ke, s2v)) {
            fprintf(stderr,
                "[Cd_v_t4] cannot map surface node to volume node: es=%d ke=%d\n",
                es,
                ke);
            continue;
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        const double (*G)[3] = G_all[ke];
        const double divu = G[0][0] + G[1][1] + G[2][2];

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;
        const int have_metric =
            t4_grad_phys_metric(xvol, dN_dx, J_inv_face, &Jacobian_face);

        double h_e_vms = 1.0e-12;
        if (have_metric) {
            h_e_vms = cbrt(fabs(Jacobian_face) / 6.0);
            if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;
        }

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double Nv[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double p_q = 0.0;
            double u_q[3] = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 4; ++a) {
                const int gid = fe->conn[ke][a];

                p_q += Nv[a] * vals->p[gid];

                for (int c = 0; c < 3; ++c) {
                    u_q[c] += Nv[a] * vals->v[gid][c];

                    if (vals->v_old != NULL) {
                        u_old_q[c] += Nv[a] * vals->v_old[gid][c];
                    } else {
                        u_old_q[c] += Nv[a] * vals->v[gid][c];
                    }
                }

                if (have_metric) {
                    for (int c = 0; c < 3; ++c) {
                        for (int d = 0; d < 3; ++d) {
                            grad_u_q[c][d] +=
                                dN_dx[a][d] * vals->v[gid][c];
                        }
                    }

                    for (int d = 0; d < 3; ++d) {
                        grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                    }
                }
            }

            (void)p_q;

            /*
              Viscous traction consistent with the volume bilinear form

                  mu grad(u):grad(w) + mu div(u) div(w).

              The corresponding boundary traction is

                  t_vis = mu grad(u) n + mu n div(u).

              If the physical Cauchy stress is required instead, replace this
              block by mu(grad u + grad u^T)n - 2/3 mu div(u)n.
            */
            double tau_n[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < 3; ++i) {
                double gradun_i = 0.0;
                for (int j = 0; j < 3; ++j) {
                    gradun_i += G[i][j] * n[j];
                }
                tau_n[i] = mu * gradun_i + mu * divu * n[i];
            }

            const double w = basis_surf->integ_weight[p] * area_factor * Jraw;

            F_vis[0] += tau_n[0] * w;
            F_vis[1] += tau_n[1] * w;
            F_vis[2] += tau_n[2] * w;

            if (include_wall_nitsche && have_metric) {
                double* grad_u_ptr[3] = {
                    grad_u_q[0],
                    grad_u_q[1],
                    grad_u_q[2]
                };

                double mu_eff_face = mu;
                double tau_face = 0.0;
                double tau_c_face = 0.0;

                BBFE_vms_mu_eff_tau(
                    &mu_eff_face,
                    &tau_face,
                    &tau_c_face,
                    J_inv_face,
                    Jacobian_face,
                    h_e_vms,
                    u_q,
                    u_old_q,
                    grad_u_ptr,
                    grad_p_q,
                    rho,
                    mu,
                    vals->dt,
                    vals->C_vms,
                    vals->vms_cap_coeff);

                const double h_n =
                    t4_tet_face_height_from_jacobians(
                        Jacobian_face,
                        Jraw);

                double f_nit_q[3];

                t4_cd_normal_nitsche_penalty_force_density(
                    u_q,
                    Uw,
                    n,
                    rho,
                    mu_eff_face,
                    h_n,
                    vals->dt,
                    gamma_n,
                    f_nit_q,
                    NULL,
                    NULL,
                    NULL);

                F_nit[0] += f_nit_q[0] * w;
                F_nit[1] += f_nit_q[1] * w;
                F_nit[2] += f_nit_q[2] * w;
            }
        }
    }

    double eU[3] = { eU_in[0], eU_in[1], eU_in[2] };
    double eP[3] = { eP_in[0], eP_in[1], eP_in[2] };
    t4_normalize3(eU);

    const double ep_dot_eu = t4_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];
    t4_normalize3(eP);

    /*
      Backward-compatible convention:
      D_out/L_out and Cd_out/Cl_out still contain dimensional forces.
      Nondimensionalization should be done once, consistently with pressure
      and Nitsche contributions, at the caller side.
    */
    const double F_body_vis[3] = {
        -F_vis[0],
        -F_vis[1],
        -F_vis[2]
    };

    const double D = t4_dot3(F_body_vis, eU);
    const double L = t4_dot3(F_body_vis, eP);
    const double D_nit = t4_dot3(F_nit, eU);
    const double L_nit = t4_dot3(F_nit, eP);

    if (D_out)          *D_out          = D;
    if (L_out)          *L_out          = L;
    if (Cd_out)         *Cd_out         = D;
    if (Cl_out)         *Cl_out         = L;
    if (D_nitsche_out)  *D_nitsche_out  = D_nit;
    if (L_nitsche_out)  *L_nitsche_out  = L_nit;
    if (Cd_nitsche_out) *Cd_nitsche_out = D_nit;
    if (Cl_nitsche_out) *Cl_nitsche_out = L_nit;

    free(G_all);
    free(owner);
    free(lface);
    return 0;
}


/* ========================================================================== */
/* Ahmed main-body / support-leg force split diagnostics                      */
/* ========================================================================== */

#ifndef T4_AHMED_SPLIT_NUM_REGIONS
#define T4_AHMED_SPLIT_NUM_REGIONS 2
#endif

enum {
    T4_AHMED_MAIN_BODY = 0,
    T4_AHMED_SUPPORTS  = 1
};

typedef struct T4AhmedBodySupportDiagnostics {
    double area[T4_AHMED_SPLIT_NUM_REGIONS];
    double projected_area[T4_AHMED_SPLIT_NUM_REGIONS];
    double int_n[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_pressure[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_molecular_viscous[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_normal_molecular_viscous[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_walllaw[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_normal_penalty_reaction[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double face_count[T4_AHMED_SPLIT_NUM_REGIONS];
    double support_z_cut;
} T4AhmedBodySupportDiagnostics;

static void t4_ahmed_split_zero(T4AhmedBodySupportDiagnostics* diag)
{
    if (diag != NULL) {
        memset(diag, 0, sizeof(*diag));
    }
}

/*
 * Region classification for the supplied geometry:
 *
 *   face-centroid z < support_z_cut  -> support legs
 *   otherwise                         -> main Ahmed body
 *
 * With the current geometry support_z_cut = 0.05 m is the nominal underbody
 * plane / ground-clearance height.  The face-centroid rule intentionally
 * keeps the vertical support-leg side faces in the support set even when their
 * top vertices lie on z = 0.05 m.
 */
static int t4_ahmed_split_region_from_face_centroid(
    const double xf[3],
    double support_z_cut)
{
    const double z_tol = 1.0e-10;
    return (xf[2] < support_z_cut - z_tol)
        ? T4_AHMED_SUPPORTS
        : T4_AHMED_MAIN_BODY;
}

int calc_ahmed_body_support_diagnostics_tet4_tri3_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw_in[3],
    double gamma_n,
    double support_z_cut,
    const double eU_in[3],
    T4AhmedBodySupportDiagnostics* diag)
{
    t4_ahmed_split_zero(diag);

    if (diag != NULL) {
        diag->support_z_cut = support_z_cut;
    }

    if (
        surf == NULL ||
        fe == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        eU_in == NULL ||
        diag == NULL) {

        fprintf(stderr, "[ahmed-split] NULL input\n");
        return -1;
    }

    if (
        fe->local_num_nodes != 4 ||
        (surf->total_num_elems > 0 && surf->local_num_nodes != 3)) {

        fprintf(
            stderr,
            "[ahmed-split] TET4/TRI3 expected: fe=%d surf=%d ne=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);
        return -1;
    }

    if (!isfinite(support_z_cut)) {
        fprintf(stderr, "[ahmed-split] invalid support_z_cut\n");
        return -1;
    }

    if (surf->total_num_elems == 0) {
        return 0;
    }

    const double Uw[3] = {
        Uw_in != NULL ? Uw_in[0] : 0.0,
        Uw_in != NULL ? Uw_in[1] : 0.0,
        Uw_in != NULL ? Uw_in[2] : 0.0
    };

    double eU[3] = {eU_in[0], eU_in[1], eU_in[2]};
    t4_normalize3(eU);

    int* owner = NULL;
    int* lface = NULL;

    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        free(owner);
        free(lface);
        return -1;
    }

    double (*G_all)[3][3] =
        (double(*)[3][3])malloc(sizeof(double) * fe->total_num_elems * 9);

    if (G_all == NULL) {
        free(owner);
        free(lface);
        return -1;
    }

    t4_compute_all_gradients_tet4_exact(fe, vals, G_all);

    const double area_factor = t4_tri_area_factor(basis_surf);

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0 || ke >= fe->total_num_elems) {
            continue;
        }

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xvol[4][3];
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        int s2v[3];
        if (!t4_surface_to_volume_nodes(surf, fe, es, ke, s2v)) {
            continue;
        }

        const int opposite_local_node =
            6 - s2v[0] - s2v[1] - s2v[2];
        if (opposite_local_node < 0 || opposite_local_node >= 4) {
            continue;
        }
        const int opposite_gid = fe->conn[ke][opposite_local_node];

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        const int region =
            t4_ahmed_split_region_from_face_centroid(xf, support_z_cut);

        diag->face_count[region] += 1.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;

        const int have_metric =
            t4_grad_phys_metric(
                xvol,
                dN_dx,
                J_inv_face,
                &Jacobian_face);

        double h_e_vms = 1.0e-12;
        if (have_metric) {
            h_e_vms = cbrt(fabs(Jacobian_face) / 6.0);
            if (h_e_vms <= 1.0e-12) h_e_vms = 1.0e-12;
        }

        const double (*G)[3] = G_all[ke];
        const double divu = G[0][0] + G[1][1] + G[2][2];

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4_norm3(nraw);
            if (Jraw <= 1.0e-300) {
                continue;
            }

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            const double w =
                basis_surf->integ_weight[p] * area_factor * Jraw;

            double Nv[4] = {0.0, 0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double p_q = 0.0;
            double u_q[3] = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 4; ++a) {
                const int gid = fe->conn[ke][a];

                p_q += Nv[a] * vals->p[gid];

                for (int c = 0; c < 3; ++c) {
                    u_q[c] += Nv[a] * vals->v[gid][c];
                    u_old_q[c] += Nv[a]
                        * (vals->v_old != NULL
                            ? vals->v_old[gid][c]
                            : vals->v[gid][c]);
                }

                if (have_metric) {
                    for (int c = 0; c < 3; ++c) {
                        for (int d = 0; d < 3; ++d) {
                            grad_u_q[c][d] +=
                                dN_dx[a][d] * vals->v[gid][c];
                        }
                    }

                    for (int d = 0; d < 3; ++d) {
                        grad_p_q[d] += dN_dx[a][d] * vals->p[gid];
                    }
                }
            }

            diag->area[region] += w;
            diag->projected_area[region] +=
                0.5 * fabs(t4_dot3(n, eU)) * w;

            for (int c = 0; c < 3; ++c) {
                diag->int_n[region][c] += n[c] * w;
                diag->force_pressure[region][c] += p_q * n[c] * w;
            }

            /* Molecular-viscous traction diagnostic, same convention as Cd_v. */
            double tau_n[3] = {0.0, 0.0, 0.0};
            double gradun_dot_n = 0.0;

            for (int i = 0; i < 3; ++i) {
                double gradun_i = 0.0;
                for (int j = 0; j < 3; ++j) {
                    gradun_i += G[i][j] * n[j];
                }

                gradun_dot_n += n[i] * gradun_i;

                tau_n[i] =
                    mu_molecular * gradun_i
                    + mu_molecular * divu * n[i];

                diag->force_molecular_viscous[region][i] +=
                    -tau_n[i] * w;
            }

            /*
             * Physical molecular normal-viscous reaction compatible with the
             * normal Nitsche traction.  The tangential viscous traction is not
             * used in the Case-C physical force because it is replaced by the
             * wall-law traction.
             */
            for (int i = 0; i < 3; ++i) {
                diag->force_normal_molecular_viscous[region][i] +=
                    -mu_molecular
                    * (gradun_dot_n + divu)
                    * n[i]
                    * w;
            }

            if (have_metric) {
                double* grad_u_ptr[3] = {
                    grad_u_q[0],
                    grad_u_q[1],
                    grad_u_q[2]
                };

                double mu_eff_face = mu_molecular;
                double tau_face = 0.0;
                double tau_c_face = 0.0;

                BBFE_vms_mu_eff_tau(
                    &mu_eff_face,
                    &tau_face,
                    &tau_c_face,
                    J_inv_face,
                    Jacobian_face,
                    h_e_vms,
                    u_q,
                    u_old_q,
                    grad_u_ptr,
                    grad_p_q,
                    rho,
                    mu_molecular,
                    vals->dt,
                    vals->C_vms,
                    vals->vms_cap_coeff);

                if (!isfinite(mu_eff_face) || mu_eff_face <= 0.0) {
                    mu_eff_face = mu_molecular;
                }

                const double h_n =
                    t4_tet_face_height_from_jacobians(
                        Jacobian_face,
                        Jraw);

                const double y_exchange = h_n;

                double u_exchange[3];
                if (g_t4_wall_exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE) {
                    u_exchange[0] = vals->v[opposite_gid][0];
                    u_exchange[1] = vals->v[opposite_gid][1];
                    u_exchange[2] = vals->v[opposite_gid][2];
                } else {
                    u_exchange[0] = u_q[0];
                    u_exchange[1] = u_q[1];
                    u_exchange[2] = u_q[2];
                }

                double ut[3] = {0.0, 0.0, 0.0};
                double beta_wall = 0.0;
                if (t4_wall_mode_uses_spalding()) {
                    t4n_spalding_wall_coefficients(
                        u_exchange,
                        Uw,
                        n,
                        rho,
                        mu_molecular,
                        y_exchange,
                        ut,
                        NULL,
                        NULL,
                        &beta_wall);
                }

                const double rwall[3] = {
                    u_q[0]-Uw[0], u_q[1]-Uw[1], u_q[2]-Uw[2]
                };
                const double rn = t4n_dot3(rwall, n);

                double beta_n = 0.0;
                if (t4_wall_mode_uses_nitsche_normal()) {
                    t4n_get_penalty_coefficients(
                        rho,
                        mu_eff_face,
                        fmax(h_n, 1.0e-12),
                        fmax(vals->dt, 1.0e-30),
                        gamma_n,
                        &beta_n,
                        NULL,
                        NULL);
                }

                for (int c = 0; c < 3; ++c) {
                    if (t4_wall_mode_uses_spalding()) {
                        diag->force_walllaw[region][c] +=
                            beta_wall * ut[c] * w;
                    }

                    if (g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP) {
                        diag->force_normal_penalty_reaction[region][c] +=
                            beta_n * rwall[c] * w;
                    }
                    else if (g_t4_wall_mode == T4_WALL_NITSCHE_FREESLIP ||
                             g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING) {
                        diag->force_normal_penalty_reaction[region][c] +=
                            beta_n * rn * n[c] * w;
                    }
                }
            }
        }
    }

    free(G_all);
    free(owner);
    free(lface);

    return 0;
}

void ahmed_body_support_diagnostics_allreduce(
    T4AhmedBodySupportDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    if (diag == NULL || monolis_com == NULL) {
        return;
    }

    monolis_allreduce_R(
        2,
        diag->area,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        2,
        diag->projected_area,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->int_n[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->force_pressure[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->force_molecular_viscous[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->force_normal_molecular_viscous[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->force_walllaw[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        6,
        &diag->force_normal_penalty_reaction[0][0],
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        2,
        diag->face_count,
        MONOLIS_MPI_SUM,
        monolis_com->comm);
}

static void t4_output_split_scalar(
    double value,
    double t,
    const char* filename,
    const char* directory)
{
    ROM_std_hlpod_output_add_calc_time(
        value,
        t,
        filename,
        directory);
}

void output_ahmed_body_support_drag_diagnostics(
    const T4AhmedBodySupportDiagnostics* diag,
    double t,
    double qA,
    double Aref,
    const double eU_in[3],
    const char* directory)
{
    if (
        diag == NULL ||
        eU_in == NULL ||
        directory == NULL ||
        !isfinite(qA) ||
        fabs(qA) <= 1.0e-300) {
        return;
    }

    double eU[3] = {eU_in[0], eU_in[1], eU_in[2]};
    t4_normalize3(eU);

    const char* prefix[T4_AHMED_SPLIT_NUM_REGIONS] = {
        "ahmed_mainbody",
        "ahmed_support"
    };

    double total_D_pressure = 0.0;
    double total_D_walllaw = 0.0;
    double total_D_molecular = 0.0;
    double total_D_normal_molecular = 0.0;
    double total_D_normal_penalty = 0.0;

    for (int r = 0; r < T4_AHMED_SPLIT_NUM_REGIONS; ++r) {
        char fname[256];

        const double D_pressure =
            t4_dot3(diag->force_pressure[r], eU);
        const double D_molecular =
            t4_dot3(diag->force_molecular_viscous[r], eU);
        const double D_normal_molecular =
            t4_dot3(diag->force_normal_molecular_viscous[r], eU);
        const double D_walllaw =
            t4_dot3(diag->force_walllaw[r], eU);
        const double D_normal_penalty =
            t4_dot3(diag->force_normal_penalty_reaction[r], eU);

        total_D_pressure += D_pressure;
        total_D_walllaw += D_walllaw;
        total_D_molecular += D_molecular;
        total_D_normal_molecular += D_normal_molecular;
        total_D_normal_penalty += D_normal_penalty;

        snprintf(fname, sizeof(fname), "%s_area.txt", prefix[r]);
        t4_output_split_scalar(diag->area[r], t, fname, directory);

        snprintf(fname, sizeof(fname), "%s_face_count.txt", prefix[r]);
        t4_output_split_scalar(diag->face_count[r], t, fname, directory);

        snprintf(fname, sizeof(fname), "%s_Aproj.txt", prefix[r]);
        t4_output_split_scalar(diag->projected_area[r], t, fname, directory);

        snprintf(fname, sizeof(fname), "%s_Aproj_over_Aref.txt", prefix[r]);
        t4_output_split_scalar(
            diag->projected_area[r] / fmax(fabs(Aref), 1.0e-300),
            t,
            fname,
            directory);

        const char axis_name[3] = {'x', 'y', 'z'};
        for (int c = 0; c < 3; ++c) {
            snprintf(
                fname,
                sizeof(fname),
                "%s_int_n%c.txt",
                prefix[r],
                axis_name[c]);
            t4_output_split_scalar(
                diag->int_n[r][c],
                t,
                fname,
                directory);
        }

        snprintf(fname, sizeof(fname), "%s_Cd_pressure.txt", prefix[r]);
        t4_output_split_scalar(D_pressure / qA, t, fname, directory);

        snprintf(
            fname,
            sizeof(fname),
            "%s_Cd_molecular_viscous_diagnostic.txt",
            prefix[r]);
        t4_output_split_scalar(D_molecular / qA, t, fname, directory);

        snprintf(
            fname,
            sizeof(fname),
            "%s_Cd_normal_molecular_viscous.txt",
            prefix[r]);
        t4_output_split_scalar(D_normal_molecular / qA, t, fname, directory);

        snprintf(fname, sizeof(fname), "%s_Cd_walllaw.txt", prefix[r]);
        t4_output_split_scalar(D_walllaw / qA, t, fname, directory);

        snprintf(
            fname,
            sizeof(fname),
            "%s_Cd_normal_penalty_reaction.txt",
            prefix[r]);
        t4_output_split_scalar(D_normal_penalty / qA, t, fname, directory);

        snprintf(
            fname,
            sizeof(fname),
            "%s_Cd_pressure_plus_walllaw.txt",
            prefix[r]);
        t4_output_split_scalar(
            (D_pressure + D_walllaw) / qA,
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "%s_Cd_modelled_total.txt",
            prefix[r]);

        double D_modelled = D_pressure;
        if (g_t4_wall_mode == T4_WALL_STRONG_NOSLIP ||
            g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP) {
            D_modelled += D_molecular;
        }
        else if (g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING ||
                 g_t4_wall_mode == T4_WALL_STRONG_NORMAL_SPALDING) {
            D_modelled += D_normal_molecular + D_walllaw;
        }
        else {
            /* free-slip: zero modeled tangential traction */
            D_modelled += D_normal_molecular;
        }

        t4_output_split_scalar(
            D_modelled / qA,
            t,
            fname,
            directory);
    }

    /* Global wall-model quantities over main body + support legs. */
    t4_output_split_scalar(
        total_D_walllaw / qA,
        t,
        "cylinder_drag_coeff_walllaw.txt",
        directory);

    t4_output_split_scalar(
        (total_D_pressure + total_D_walllaw) / qA,
        t,
        "cylinder_drag_coeff_pressure_plus_walllaw.txt",
        directory);

    double total_D_modelled = total_D_pressure;
    if (g_t4_wall_mode == T4_WALL_STRONG_NOSLIP ||
        g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP) {
        total_D_modelled += total_D_molecular;
    }
    else if (g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING ||
             g_t4_wall_mode == T4_WALL_STRONG_NORMAL_SPALDING) {
        total_D_modelled += total_D_normal_molecular + total_D_walllaw;
    }
    else {
        total_D_modelled += total_D_normal_molecular;
    }

    t4_output_split_scalar(
        total_D_modelled / qA,
        t,
        "cylinder_drag_coeff_modelled_total.txt",
        directory);

    t4_output_split_scalar(
        total_D_normal_molecular / qA,
        t,
        "cylinder_drag_coeff_normal_molecular_viscous.txt",
        directory);

    t4_output_split_scalar(
        total_D_molecular / qA,
        t,
        "cylinder_drag_coeff_molecular_viscous_splitcheck.txt",
        directory);

    t4_output_split_scalar(
        total_D_normal_penalty / qA,
        t,
        "cylinder_drag_coeff_normal_penalty_reaction.txt",
        directory);

    t4_output_split_scalar(
        diag->support_z_cut,
        t,
        "ahmed_support_z_cut.txt",
        directory);
}

/*
 * Backward-compatible entry point: original viscous force only.
 * Use calc_Cd_v_tet4_tri3_nitsche() when the wall is imposed by Nitsche
 * in Case C and you want the normal Nitsche penalty reaction returned
 * separately.
 */
int calc_Cd_v_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out)
{
    return calc_Cd_v_tet4_tri3_core(
        surf,
        fe,
        basis_surf,
        monolis_com,
        vals,
        rho,
        mu,
        Uinf,
        Aref,
        eU_in,
        eP_in,
        0,
        NULL,
        0.0,
        D_out,
        L_out,
        Cd_out,
        Cl_out,
        NULL,
        NULL,
        NULL,
        NULL);
}

/*
 * Extended entry point.
 *
 * D_out/L_out and Cd_out/Cl_out are the same viscous-only quantities as the
 * verified calc_Cd_v() path.  D_nitsche_out/L_nitsche_out and
 * Cd_nitsche_out/Cl_nitsche_out are separate diagnostics for the normal
 * normal Nitsche penalty reaction. The physical tangential wall-law reaction
 * is reported separately by the Ahmed body/support split diagnostics.
 */int calc_Cd_v_tet4_tri3_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    /*
     * 空rankの局所寄与はすべて0。
     */
    if (D_out != NULL) {
        *D_out = 0.0;
    }

    if (L_out != NULL) {
        *L_out = 0.0;
    }

    if (Cd_out != NULL) {
        *Cd_out = 0.0;
    }

    if (Cl_out != NULL) {
        *Cl_out = 0.0;
    }

    if (D_nitsche_out != NULL) {
        *D_nitsche_out = 0.0;
    }

    if (L_nitsche_out != NULL) {
        *L_nitsche_out = 0.0;
    }

    if (Cd_nitsche_out != NULL) {
        *Cd_nitsche_out = 0.0;
    }

    if (Cl_nitsche_out != NULL) {
        *Cl_nitsche_out = 0.0;
    }

    if (
        surf == NULL ||
        fe == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        eU_in == NULL ||
        eP_in == NULL ||
        Uw == NULL) {

        fprintf(
            stderr,
            "[Cd_v_t4_nitsche] NULL input\n");

        return -1;
    }

    /*
     * 壁面要素を所有しないrankは正常。
     * main側のMPI_SUMへ0を渡す。
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    if (
        fe->local_num_nodes != 4 ||
        surf->local_num_nodes != 3) {

        fprintf(
            stderr,
            "[Cd_v_t4_nitsche] TET4/TRI3 expected, "
            "got fe=%d surf=%d ne_surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        return -1;
    }

    return calc_Cd_v_tet4_tri3_core(
        surf,
        fe,
        basis_surf,
        monolis_com,
        vals,
        rho,
        mu,
        Uinf,
        Aref,
        eU_in,
        eP_in,
        1,
        Uw,
        gamma_n,
        D_out,
        L_out,
        Cd_out,
        Cl_out,
        D_nitsche_out,
        L_nitsche_out,
        Cd_nitsche_out,
        Cl_nitsche_out);
}


/*
 * Dispatch wrapper including TET4/TRI3 Nitsche wall reaction.
 * For HEX8/QUAD4 this currently falls back to the verified existing calc_Cd_v().
 * For TET4/TRI3 it calls calc_Cd_v_tet4_tri3_nitsche().
 */int calc_Cd_v_hex_or_tet_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out)
{
    if (D_out != NULL) {
        *D_out = 0.0;
    }

    if (L_out != NULL) {
        *L_out = 0.0;
    }

    if (Cd_out != NULL) {
        *Cd_out = 0.0;
    }

    if (Cl_out != NULL) {
        *Cl_out = 0.0;
    }

    if (D_nitsche_out != NULL) {
        *D_nitsche_out = 0.0;
    }

    if (L_nitsche_out != NULL) {
        *L_nitsche_out = 0.0;
    }

    if (Cd_nitsche_out != NULL) {
        *Cd_nitsche_out = 0.0;
    }

    if (Cl_nitsche_out != NULL) {
        *Cl_nitsche_out = 0.0;
    }

    if (
        surf == NULL ||
        fe == NULL ||
        basis == NULL ||
        vals == NULL) {

        fprintf(
            stderr,
            "[Cd_v_nitsche] NULL input\n");

        return -1;
    }

    /*
     * 空rankは正常なゼロ寄与。
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    if (fe->local_num_nodes == 8) {
        if (surf->local_num_nodes != 4) {
            fprintf(
                stderr,
                "[Cd_v_nitsche] HEX8 volume requires QUAD4 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return -1;
        }

        return calc_Cd_v(
            surf,
            fe,
            basis,
            monolis_com,
            vals,
            rho,
            mu,
            Uinf,
            Aref,
            eU_in,
            eP_in,
            D_out,
            L_out,
            Cd_out,
            Cl_out);
    }

    if (fe->local_num_nodes == 4) {
        if (surf->local_num_nodes != 3) {
            fprintf(
                stderr,
                "[Cd_v_nitsche] TET4 volume requires TRI3 surface, "
                "but surf->local_num_nodes=%d "
                "surf->total_num_elems=%d\n",
                surf->local_num_nodes,
                surf->total_num_elems);

            return -1;
        }

        return calc_Cd_v_tet4_tri3_nitsche(
            surf,
            fe,
            basis,
            monolis_com,
            vals,
            rho,
            mu,
            Uinf,
            Aref,
            eU_in,
            eP_in,
            Uw,
            gamma_n,
            D_out,
            L_out,
            Cd_out,
            Cl_out,
            D_nitsche_out,
            L_nitsche_out,
            Cd_nitsche_out,
            Cl_nitsche_out);
    }

    fprintf(
        stderr,
        "[Cd_v_nitsche] unsupported fe->local_num_nodes=%d. "
        "Expected 8 HEX8 or 4 TET4.\n",
        fe->local_num_nodes);

    return -1;
}



static void t4wd_normalize3(double a[3])
{
    const double n = t4wd_norm3(a);
    if (n > 1.0e-300) {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}

static void t4wd_prepare_dirs(
    double eU[3],
    double eP[3],
    double eS[3],
    const double eU_in[3],
    const double eP_in[3])
{
    eU[0] = eU_in[0];
    eU[1] = eU_in[1];
    eU[2] = eU_in[2];
    t4wd_normalize3(eU);

    eP[0] = eP_in[0];
    eP[1] = eP_in[1];
    eP[2] = eP_in[2];

    const double ep_dot_eu = t4wd_dot3(eP, eU);
    eP[0] -= ep_dot_eu * eU[0];
    eP[1] -= ep_dot_eu * eU[1];
    eP[2] -= ep_dot_eu * eU[2];

    if (t4wd_norm3(eP) <= 1.0e-14) {
        eP[0] = 0.0;
        eP[1] = 1.0;
        eP[2] = 0.0;
        const double dot_tmp = t4wd_dot3(eP, eU);
        eP[0] -= dot_tmp * eU[0];
        eP[1] -= dot_tmp * eU[1];
        eP[2] -= dot_tmp * eU[2];
    }
    t4wd_normalize3(eP);

    t4wd_cross3(eS, eU, eP);
    t4wd_normalize3(eS);
}

static void t4_wall_penalty_scale_init(
    T4WallPenaltyScaleDiagnostics* diag)
{
    if (diag == NULL) {
        return;
    }

    diag->area   = 0.0;
    diag->h_sum  = 0.0;
    diag->h_mean = 0.0;
    diag->h_min  = 1.0e300;
    diag->h_max  = -1.0e300;
}

static void t4_wall_penalty_scale_finalize(
    T4WallPenaltyScaleDiagnostics* diag)
{
    if (diag == NULL) {
        return;
    }

    if (diag->area > 0.0) {
        diag->h_mean = diag->h_sum / diag->area;
    } else {
        diag->h_mean = 0.0;
        diag->h_min  = 0.0;
        diag->h_max  = 0.0;
    }
}


int calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    T4WallPenaltyScaleDiagnostics* scale_diag)
{
    for (int g = 0; g < num_gamma; ++g) {
        if (D_penalty_by_gamma != NULL) D_penalty_by_gamma[g] = 0.0;
        if (L_penalty_by_gamma != NULL) L_penalty_by_gamma[g] = 0.0;
        if (D_tangent_by_gamma != NULL) D_tangent_by_gamma[g] = 0.0;
        if (L_tangent_by_gamma != NULL) L_tangent_by_gamma[g] = 0.0;
        if (D_total_by_gamma != NULL) D_total_by_gamma[g] = 0.0;
        if (L_total_by_gamma != NULL) L_total_by_gamma[g] = 0.0;
        if (D_penalty_mu_by_gamma != NULL) D_penalty_mu_by_gamma[g] = 0.0;
        if (L_penalty_mu_by_gamma != NULL) L_penalty_mu_by_gamma[g] = 0.0;
        if (D_penalty_dt_by_gamma != NULL) D_penalty_dt_by_gamma[g] = 0.0;
        if (L_penalty_dt_by_gamma != NULL) L_penalty_dt_by_gamma[g] = 0.0;
    }

    if (scale_diag != NULL) {
        t4_wall_penalty_scale_init(scale_diag);
    }

    if (surf->local_num_nodes != 3 || fe->local_num_nodes != 4) {
        fprintf(stderr,
            "[wall_gamma] TET4/TRI3 expected, got fe=%d surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes);
        return -1;
    }

    int* owner = NULL;
    int* lface = NULL;
    if (t4_make_owner_map(surf, fe, &owner, &lface) != 0) {
        fprintf(stderr, "[wall_gamma] surface-owner mapping failed\n");
        return -1;
    }

    double eU[3], eP[3], eS[3];
    t4wd_prepare_dirs(eU, eP, eS, eU_in, eP_in);
    (void)eS;

    const double Uw_loc[3] = {
        Uw ? Uw[0] : 0.0,
        Uw ? Uw[1] : 0.0,
        Uw ? Uw[2] : 0.0
    };

    const double area_factor = t4_tri_area_factor(basis_surf);

    for (int es = 0; es < surf->total_num_elems; ++es) {
        const int ke = owner[es];
        if (ke < 0) continue;

        double xtri[3][3];
        for (int a = 0; a < 3; ++a) {
            const int gid = surf->conn[es][a];
            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        double xf[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }
        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        double xc[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 4; ++a) {
            const int gid = fe->conn[ke][a];
            xc[0] += fe->x[gid][0];
            xc[1] += fe->x[gid][1];
            xc[2] += fe->x[gid][2];
        }
        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < basis_surf->num_integ_points; ++p) {
            double nraw[3];
            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw = t4wd_norm3(nraw);
            if (Jraw <= 1.0e-300) continue;

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            if (t4wd_dot3(svec, n) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            double u_q[3] = {0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
                const int gid = surf->conn[es][a];
                const double Na = basis_surf->N[p][a];
                u_q[0] += Na * vals->v[gid][0];
                u_q[1] += Na * vals->v[gid][1];
                u_q[2] += Na * vals->v[gid][2];
            }

            double r_q[3] = {
                u_q[0] - Uw_loc[0],
                u_q[1] - Uw_loc[1],
                u_q[2] - Uw_loc[2]
            };

            double h_n = fabs(
                (xf[0] - xc[0]) * n[0]
              + (xf[1] - xc[1]) * n[1]
              + (xf[2] - xc[2]) * n[2]);

            if (h_n <= 1.0e-12) h_n = 1.0e-12;

            const double dt_eff = fmax(vals->dt, 1.0e-30);

            const double beta_mu_base = mu / h_n;
            const double beta_dt_base =
                g_t4_nitsche_dt_penalty_coeff * rho * h_n / dt_eff;

            const double w =
                basis_surf->integ_weight[p] * area_factor * Jraw;

            if (scale_diag != NULL) {
                scale_diag->area += w;
                scale_diag->h_sum += h_n * w;
                scale_diag->h_min = fmin(scale_diag->h_min, h_n);
                scale_diag->h_max = fmax(scale_diag->h_max, h_n);
            }

            /*
             * Case C normal-velocity Nitsche penalty reaction:
             *
             *     f_penalty = - beta ((u-Uw).n) n
             *
             * Here beta = gamma * (mu/h + c_dt rho h/dt).
             */
            const double rn_q = t4wd_dot3(r_q, n);

            for (int g = 0; g < num_gamma; ++g) {
                const double beta_mu = gamma_list[g] * beta_mu_base;
                const double beta_dt = gamma_list[g] * beta_dt_base;
                const double beta    = beta_mu + beta_dt;

                double f_penalty[3];
                double f_penalty_mu[3];
                double f_penalty_dt[3];

                for (int a = 0; a < 3; ++a) {
                    f_penalty_mu[a] = -beta_mu * rn_q * n[a];
                    f_penalty_dt[a] = -beta_dt * rn_q * n[a];
                    f_penalty[a]    = -beta    * rn_q * n[a];
                }

                const double D_penalty =
                    t4wd_dot3(f_penalty, eU) * w;
                const double L_penalty =
                    t4wd_dot3(f_penalty, eP) * w;

                const double D_penalty_mu =
                    t4wd_dot3(f_penalty_mu, eU) * w;
                const double L_penalty_mu =
                    t4wd_dot3(f_penalty_mu, eP) * w;

                const double D_penalty_dt =
                    t4wd_dot3(f_penalty_dt, eU) * w;
                const double L_penalty_dt =
                    t4wd_dot3(f_penalty_dt, eP) * w;

                if (D_penalty_by_gamma != NULL) {
                    D_penalty_by_gamma[g] += D_penalty;
                }
                if (L_penalty_by_gamma != NULL) {
                    L_penalty_by_gamma[g] += L_penalty;
                }

                /*
                 * Case C has no tangential Nitsche penalty; the tangential wall law is a
                 * tangential wall-law reaction. Keep these outputs zero.
                 */
                if (D_tangent_by_gamma != NULL) {
                    D_tangent_by_gamma[g] += 0.0;
                }
                if (L_tangent_by_gamma != NULL) {
                    L_tangent_by_gamma[g] += 0.0;
                }

                if (D_total_by_gamma != NULL) {
                    D_total_by_gamma[g] += D_penalty;
                }
                if (L_total_by_gamma != NULL) {
                    L_total_by_gamma[g] += L_penalty;
                }

                if (D_penalty_mu_by_gamma != NULL) {
                    D_penalty_mu_by_gamma[g] += D_penalty_mu;
                }
                if (L_penalty_mu_by_gamma != NULL) {
                    L_penalty_mu_by_gamma[g] += L_penalty_mu;
                }

                if (D_penalty_dt_by_gamma != NULL) {
                    D_penalty_dt_by_gamma[g] += D_penalty_dt;
                }
                if (L_penalty_dt_by_gamma != NULL) {
                    L_penalty_dt_by_gamma[g] += L_penalty_dt;
                }
            }
        }
    }

    free(owner);
    free(lface);

    if (scale_diag != NULL) {
        t4_wall_penalty_scale_finalize(scale_diag);
    }

    return 0;
}

static int t4wd_classify_region(
    const double n[3],
    const double eU[3],
    const double eP[3],
    const double eS[3])
{
    const double du = t4wd_dot3(n, eU);
    const double dp = t4wd_dot3(n, eP);
    const double ds = t4wd_dot3(n, eS);
    const double au = fabs(du);
    const double ap = fabs(dp);
    const double as = fabs(ds);
    const double axis_threshold = 0.70;

    if (au >= ap && au >= as && au >= axis_threshold) {
        return (du >= 0.0) ? T4_WALL_DIAG_FRONT_POS_DRAG
                           : T4_WALL_DIAG_REAR_NEG_DRAG;
    }
    if (ap >= au && ap >= as && ap >= axis_threshold) {
        return (dp >= 0.0) ? T4_WALL_DIAG_TOP_POS_LIFT
                           : T4_WALL_DIAG_BOTTOM_NEG_LIFT;
    }
    if (as >= au && as >= ap && as >= axis_threshold) {
        return (ds >= 0.0) ? T4_WALL_DIAG_SIDE_POS
                           : T4_WALL_DIAG_SIDE_NEG;
    }

    return T4_WALL_DIAG_SLANT_MIXED;
}

static void t4_wall_diag_init(T4WallDiagnostics* diag)
{
    diag->area = 0.0;
    diag->int_n[0] = 0.0;
    diag->int_n[1] = 0.0;
    diag->int_n[2] = 0.0;

    diag->max_abs_un_over_Uref = 0.0;
    diag->int_abs_un = 0.0;
    diag->int_un = 0.0;
    diag->int_un2 = 0.0;
    diag->mean_abs_un_over_Uref = 0.0;
    diag->mean_un_over_Uref = 0.0;
    diag->rms_un_over_Uref = 0.0;
    diag->Uref_scale = 1.0;

    diag->max_norm_r_over_Uref = 0.0;
    diag->int_norm_r = 0.0;
    diag->mean_norm_r_over_Uref = 0.0;

    diag->p_min = 1.0e300;
    diag->p_max = -1.0e300;
    diag->p_mean = 0.0;
    diag->p_rms = 0.0;
    diag->p_area_sum = 0.0;
    diag->p2_area_sum = 0.0;

    diag->cp_min = 1.0e300;
    diag->cp_max = -1.0e300;
    diag->cp_mean = 0.0;
    diag->cp_rms = 0.0;
    diag->cp_area_sum = 0.0;
    diag->cp2_area_sum = 0.0;

    diag->int_pnU = 0.0;

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->area_region[r] = 0.0;
        diag->int_pnU_region[r] = 0.0;
        diag->int_abs_un_region[r] = 0.0;
    }
}

static void t4_wall_diag_finalize(T4WallDiagnostics* diag)
{
    if (diag->area > 0.0) {
        diag->p_mean = diag->p_area_sum / diag->area;
        diag->p_rms = sqrt(fmax(diag->p2_area_sum / diag->area, 0.0));

        diag->cp_mean = diag->cp_area_sum / diag->area;
        diag->cp_rms = sqrt(fmax(diag->cp2_area_sum / diag->area, 0.0));

        diag->mean_abs_un_over_Uref =
            diag->int_abs_un / (diag->Uref_scale * diag->area);

        diag->mean_un_over_Uref =
            diag->int_un / (diag->Uref_scale * diag->area);

        diag->rms_un_over_Uref =
            sqrt(fmax(diag->int_un2 / diag->area, 0.0))
            / diag->Uref_scale;

        diag->mean_norm_r_over_Uref =
            diag->int_norm_r / (diag->Uref_scale * diag->area);
    } else {
        diag->p_min = 0.0;
        diag->p_max = 0.0;
        diag->cp_min = 0.0;
        diag->cp_max = 0.0;

        diag->mean_abs_un_over_Uref = 0.0;
        diag->mean_un_over_Uref = 0.0;
        diag->rms_un_over_Uref = 0.0;
        diag->mean_norm_r_over_Uref = 0.0;
    }
}
int calc_wall_diagnostics_tet4_tri3_with_Uw(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    T4WallDiagnostics* diag)
{
    if (
        surf == NULL ||
        fe == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        diag == NULL) {

        fprintf(
            stderr,
            "[wall_diag_with_Uw] NULL input\n");

        return -1;
    }

    t4_wall_diag_init(diag);

    /*
     * wall_diagnostics_allreduce()ではUref_scaleを集約しないため、
     * 表面要素を持たないrankを含めて、全rankで同じ値を設定する。
     */
    const double U_scale =
        fmax(
            fabs(Uref),
            1.0e-300);

    diag->Uref_scale = U_scale;

    if (surf->total_num_elems < 0) {
        fprintf(
            stderr,
            "[wall_diag_with_Uw] "
            "invalid number of surface elements: %d\n",
            surf->total_num_elems);

        return -1;
    }

    /*
     * 表面要素を持たないrankは正常なゼロ寄与。
     *
     * p_min、cp_minなどのMIN用中立値を維持するため、
     * この時点ではt4_wall_diag_finalize()を呼ばない。
     */
    if (surf->total_num_elems == 0) {
        return 0;
    }

    if (
        surf->local_num_nodes != 3 ||
        fe->local_num_nodes != 4) {

        fprintf(
            stderr,
            "[wall_diag_with_Uw] "
            "TET4/TRI3 expected, "
            "got fe=%d surf=%d ne_surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        return -1;
    }

    int* owner = NULL;
    int* lface = NULL;
    int num_unmapped_local = 0;

    if (
        t4n_make_owner_map_overlap(
            surf,
            fe,
            &owner,
            &lface,
            &num_unmapped_local) != 0
    ) {
        fprintf(
            stderr,
            "[wall_diag_with_Uw] "
            "overlap owner-map construction failed\n");

        return -1;
    }

    double eU[3];
    double eP[3];
    double eS[3];

    t4wd_prepare_dirs(
        eU,
        eP,
        eS,
        eU_in,
        eP_in);

    /*
     * Uw == NULLの場合は固定壁Uw=0として扱う。
     */
    const double Uw_loc[3] = {
        Uw != NULL ? Uw[0] : 0.0,
        Uw != NULL ? Uw[1] : 0.0,
        Uw != NULL ? Uw[2] : 0.0
    };

    double qref =
        0.5 * rho * Uref * Uref;

    if (fabs(qref) <= 1.0e-300) {
        qref = 1.0;
    }

    const double area_factor =
        t4_tri_area_factor(basis_surf);

    for (int es = 0;
         es < surf->total_num_elems;
         ++es) {

        const int ke = owner[es];

        if (ke < 0) {
            continue;
        }

        /*
         * TRI3の節点座標。
         */
        double xtri[3][3];

        for (int a = 0; a < 3; ++a) {
            const int gid =
                surf->conn[es][a];

            xtri[a][0] = fe->x[gid][0];
            xtri[a][1] = fe->x[gid][1];
            xtri[a][2] = fe->x[gid][2];
        }

        /*
         * 表面三角形の重心。
         */
        double xf[3] = {
            0.0,
            0.0,
            0.0
        };

        for (int a = 0; a < 3; ++a) {
            xf[0] += xtri[a][0];
            xf[1] += xtri[a][1];
            xf[2] += xtri[a][2];
        }

        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        /*
         * owner TET4の重心。
         */
        double xc[3] = {
            0.0,
            0.0,
            0.0
        };

        for (int a = 0; a < 4; ++a) {
            const int gid =
                fe->conn[ke][a];

            xc[0] += fe->x[gid][0];
            xc[1] += fe->x[gid][1];
            xc[2] += fe->x[gid][2];
        }

        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        /*
         * owner TET重心から壁面重心へ向かうベクトル。
         * 表面法線を流体領域外向きへそろえるために使う。
         */
        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0;
             p < basis_surf->num_integ_points;
             ++p) {

            /*
             * TRI3の非正規化法線。
             */
            double nraw[3];

            t4_tri3_normal(
                xtri,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);

            const double Jraw =
                t4wd_norm3(nraw);

            if (Jraw <= 1.0e-300) {
                continue;
            }

            double n[3] = {
                nraw[0] / Jraw,
                nraw[1] / Jraw,
                nraw[2] / Jraw
            };

            /*
             * 法線を流体領域から壁面外側へ向ける。
             */
            if (
                t4wd_dot3(
                    svec,
                    n) < 0.0) {

                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            /*
             * 積分点の圧力と速度。
             */
            double p_q = 0.0;

            double u_q[3] = {
                0.0,
                0.0,
                0.0
            };

            for (int a = 0; a < 3; ++a) {
                const int gid =
                    surf->conn[es][a];

                const double Na =
                    basis_surf->N[p][a];

                p_q +=
                    Na * vals->p[gid];

                u_q[0] +=
                    Na * vals->v[gid][0];

                u_q[1] +=
                    Na * vals->v[gid][1];

                u_q[2] +=
                    Na * vals->v[gid][2];
            }

            /*
             * ============================================================
             * Full-vector weak wall BC error
             *
             *     r_q = u_q - Uw
             * ============================================================
             */
            const double r_q[3] = {
                u_q[0] - Uw_loc[0],
                u_q[1] - Uw_loc[1],
                u_q[2] - Uw_loc[2]
            };

            /*
             * 面積積分重み。
             */
            const double w =
                basis_surf->integ_weight[p]
                * area_factor
                * Jraw;

            /*
             * 法線方向壁面速度誤差：
             *
             *     un = (u-Uw).n
             */
            const double un =
                t4wd_dot3(
                    r_q,
                    n);

            const double abs_un =
                fabs(un);

            /*
             * full-vector壁面速度誤差：
             *
             *     norm_r = ||u-Uw||
             */
            const double norm_r =
                t4wd_norm3(r_q);

            const double cp =
                (p_q - p_ref) / qref;

            const double pnU =
                p_q
                * t4wd_dot3(
                    n,
                    eU);

            const int region =
                t4wd_classify_region(
                    n,
                    eU,
                    eP,
                    eS);

            /*
             * 面積と法線積分。
             */
            diag->area += w;

            diag->int_n[0] +=
                n[0] * w;

            diag->int_n[1] +=
                n[1] * w;

            diag->int_n[2] +=
                n[2] * w;

            /*
             * 法線方向速度誤差。
             */
            diag->int_abs_un +=
                abs_un * w;

            diag->int_un +=
                un * w;

            diag->int_un2 +=
                un * un * w;

            diag->max_abs_un_over_Uref =
                fmax(
                    diag->max_abs_un_over_Uref,
                    abs_un / U_scale);

            /*
             * Full-vector速度誤差。
             */
            diag->int_norm_r +=
                norm_r * w;

            diag->max_norm_r_over_Uref =
                fmax(
                    diag->max_norm_r_over_Uref,
                    norm_r / U_scale);

            /*
             * 圧力統計。
             */
            diag->p_min =
                fmin(
                    diag->p_min,
                    p_q);

            diag->p_max =
                fmax(
                    diag->p_max,
                    p_q);

            diag->p_area_sum +=
                p_q * w;

            diag->p2_area_sum +=
                p_q * p_q * w;

            /*
             * 圧力係数統計。
             */
            diag->cp_min =
                fmin(
                    diag->cp_min,
                    cp);

            diag->cp_max =
                fmax(
                    diag->cp_max,
                    cp);

            diag->cp_area_sum +=
                cp * w;

            diag->cp2_area_sum +=
                cp * cp * w;

            /*
             * 流れ方向圧力力積分。
             */
            diag->int_pnU +=
                pnU * w;

            /*
             * 壁面領域別統計。
             */
            diag->area_region[region] +=
                w;

            diag->int_pnU_region[region] +=
                pnU * w;

            diag->int_abs_un_region[region] +=
                abs_un * w;
        }
    }

    free(owner);
    free(lface);

    /*
     * 非空rankの局所統計を確定する。
     * global allreduce後にも再度finalizeされる。
     */
    t4_wall_diag_finalize(diag);

    return 0;
}

int calc_wall_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    T4WallDiagnostics* diag)
{
    const double Uw_zero[3] = {
        0.0,
        0.0,
        0.0
    };

    return calc_wall_diagnostics_tet4_tri3_with_Uw(
        surf,
        fe,
        basis_surf,
        vals,
        rho,
        Uref,
        p_ref,
        eU_in,
        eP_in,
        Uw_zero,
        diag);
}

void wall_diagnostics_allreduce(
    T4WallDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    double sum_pack[64];
    int k = 0;

    sum_pack[k++] = diag->area;
    sum_pack[k++] = diag->int_abs_un;
    sum_pack[k++] = diag->int_un;
    sum_pack[k++] = diag->int_un2;
    sum_pack[k++] = diag->int_pnU;
    sum_pack[k++] = diag->p_area_sum;
    sum_pack[k++] = diag->p2_area_sum;
    sum_pack[k++] = diag->cp_area_sum;
    sum_pack[k++] = diag->cp2_area_sum;
    sum_pack[k++] = diag->int_n[0];
    sum_pack[k++] = diag->int_n[1];
    sum_pack[k++] = diag->int_n[2];

    /* full-vector weak essential BC diagnostic */
    sum_pack[k++] = diag->int_norm_r;

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->area_region[r];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->int_pnU_region[r];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        sum_pack[k++] = diag->int_abs_un_region[r];
    }

    monolis_allreduce_R(k, sum_pack, MONOLIS_MPI_SUM, monolis_com->comm);

    k = 0;
    diag->area = sum_pack[k++];
    diag->int_abs_un = sum_pack[k++];
    diag->int_un = sum_pack[k++];
    diag->int_un2 = sum_pack[k++];
    diag->int_pnU = sum_pack[k++];
    diag->p_area_sum = sum_pack[k++];
    diag->p2_area_sum = sum_pack[k++];
    diag->cp_area_sum = sum_pack[k++];
    diag->cp2_area_sum = sum_pack[k++];
    diag->int_n[0] = sum_pack[k++];
    diag->int_n[1] = sum_pack[k++];
    diag->int_n[2] = sum_pack[k++];

    diag->int_norm_r = sum_pack[k++];

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->area_region[r] = sum_pack[k++];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->int_pnU_region[r] = sum_pack[k++];
    }
    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        diag->int_abs_un_region[r] = sum_pack[k++];
    }

    {
        double max_pack[4] = {
            diag->max_abs_un_over_Uref,
            diag->max_norm_r_over_Uref,
            diag->p_max,
            diag->cp_max
        };

        monolis_allreduce_R(
            4,
            max_pack,
            MONOLIS_MPI_MAX,
            monolis_com->comm);

        diag->max_abs_un_over_Uref  = max_pack[0];
        diag->max_norm_r_over_Uref = max_pack[1];
        diag->p_max = max_pack[2];
        diag->cp_max = max_pack[3];
    }

    {
        double min_pack_as_negative[2] = {
            -diag->p_min,
            -diag->cp_min
        };

        monolis_allreduce_R(
            2,
            min_pack_as_negative,
            MONOLIS_MPI_MAX,
            monolis_com->comm);

        diag->p_min = -min_pack_as_negative[0];
        diag->cp_min = -min_pack_as_negative[1];
    }

    t4_wall_diag_finalize(diag);
}
void output_wall_diagnostics_tet4_tri3(
    const T4WallDiagnostics* diag,
    double t,
    const char* directory)
{
    ROM_std_hlpod_output_add_calc_time(
        diag->max_abs_un_over_Uref,
        t,
        "wall_max_abs_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_abs_un_over_Uref,
        t,
        "wall_mean_abs_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_un_over_Uref,
        t,
        "wall_mean_un_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->rms_un_over_Uref,
        t,
        "wall_rms_un_over_Uref.txt",
        directory);

    /*
     * Full-vector weak essential BC diagnostics:
     *
     *     r = u - Uw
     */
    ROM_std_hlpod_output_add_calc_time(
        diag->max_norm_r_over_Uref,
        t,
        "wall_max_norm_u_minus_Uw_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->mean_norm_r_over_Uref,
        t,
        "wall_mean_norm_u_minus_Uw_over_Uref.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_norm_r,
        t,
        "wall_int_norm_u_minus_Uw.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->area,
        t,
        "wall_area_total.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_abs_un,
        t,
        "wall_int_abs_un.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_un,
        t,
        "wall_int_un_signed.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_min,
        t,
        "wall_Cp_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_max,
        t,
        "wall_Cp_max.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_rms,
        t,
        "wall_Cp_rms.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->cp_mean,
        t,
        "wall_Cp_mean.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_min,
        t,
        "wall_p_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_max,
        t,
        "wall_p_max.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->p_rms,
        t,
        "wall_p_rms.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_pnU,
        t,
        "wall_pressure_drag_int_pnU.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[0],
        t,
        "wall_int_nx.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[1],
        t,
        "wall_int_ny.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->int_n[2],
        t,
        "wall_int_nz.txt",
        directory);

    for (int r = 0; r < T4_WALL_DIAG_NUM_REGIONS; ++r) {
        char fname[256];

        snprintf(
            fname,
            sizeof(fname),
            "wall_pressure_drag_%s.txt",
            t4_wall_diag_region_name[r]);

        ROM_std_hlpod_output_add_calc_time(
            diag->int_pnU_region[r],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_area_%s.txt",
            t4_wall_diag_region_name[r]);

        ROM_std_hlpod_output_add_calc_time(
            diag->area_region[r],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_int_abs_un_%s.txt",
            t4_wall_diag_region_name[r]);

        ROM_std_hlpod_output_add_calc_time(
            diag->int_abs_un_region[r],
            t,
            fname,
            directory);
    }
}


int calc_wall_nitsche_gamma_sweep_split_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma)
{
    return calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
        surf,
        fe,
        basis_surf,
        vals,
        rho,
        mu,
        Uw,
        eU_in,
        eP_in,
        gamma_list,
        num_gamma,
        D_penalty_by_gamma,
        L_penalty_by_gamma,
        D_tangent_by_gamma,
        L_tangent_by_gamma,
        D_total_by_gamma,
        L_total_by_gamma,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL);
}

int calc_wall_nitsche_gamma_sweep_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_nitsche_by_gamma,
    double* L_nitsche_by_gamma)
{
    return calc_wall_nitsche_gamma_sweep_split_tet4_tri3(
        surf,
        fe,
        basis_surf,
        vals,
        rho,
        mu,
        Uw,
        eU_in,
        eP_in,
        gamma_list,
        num_gamma,
        NULL,
        NULL,
        NULL,
        NULL,
        D_nitsche_by_gamma,
        L_nitsche_by_gamma);
}

void wall_gamma_sweep_allreduce(
    int num_gamma,
    double* D_nitsche_by_gamma,
    double* L_nitsche_by_gamma,
    MONOLIS_COM* monolis_com)
{
    monolis_allreduce_R(
        num_gamma,
        D_nitsche_by_gamma,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        num_gamma,
        L_nitsche_by_gamma,
        MONOLIS_MPI_SUM,
        monolis_com->comm);
}

void wall_gamma_sweep_split_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    MONOLIS_COM* monolis_com)
{
    if (D_penalty_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_tangent_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_tangent_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_tangent_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_tangent_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_total_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_total_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_total_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_total_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
}

void wall_gamma_sweep_decomposed_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    MONOLIS_COM* monolis_com)
{
    wall_gamma_sweep_split_allreduce(
        num_gamma,
        D_penalty_by_gamma,
        L_penalty_by_gamma,
        D_tangent_by_gamma,
        L_tangent_by_gamma,
        D_total_by_gamma,
        L_total_by_gamma,
        monolis_com);

    if (D_penalty_mu_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_mu_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_mu_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_mu_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (D_penalty_dt_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            D_penalty_dt_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
    if (L_penalty_dt_by_gamma != NULL) {
        monolis_allreduce_R(
            num_gamma,
            L_penalty_dt_by_gamma,
            MONOLIS_MPI_SUM,
            monolis_com->comm);
    }
}

void wall_penalty_scale_diagnostics_allreduce(
    T4WallPenaltyScaleDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    double sum_pack[2] = {
        diag->area,
        diag->h_sum
    };

    monolis_allreduce_R(
        2,
        sum_pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    diag->area = sum_pack[0];
    diag->h_sum = sum_pack[1];

    {
        double max_pack[1] = {diag->h_max};
        monolis_allreduce_R(
            1,
            max_pack,
            MONOLIS_MPI_MAX,
            monolis_com->comm);
        diag->h_max = max_pack[0];
    }

    {
        double min_pack_as_negative[1] = {-diag->h_min};
        monolis_allreduce_R(
            1,
            min_pack_as_negative,
            MONOLIS_MPI_MAX,
            monolis_com->comm);
        diag->h_min = -min_pack_as_negative[0];
    }

    t4_wall_penalty_scale_finalize(diag);
}

void output_wall_penalty_scale_diagnostics_tet4_tri3(
    const T4WallPenaltyScaleDiagnostics* diag,
    double t,
    const char* directory)
{
    ROM_std_hlpod_output_add_calc_time(
        g_t4_nitsche_dt_penalty_coeff,
        t,
        "wall_nitsche_dt_penalty_coeff.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_mean,
        t,
        "wall_h_n_mean.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_min,
        t,
        "wall_h_n_min.txt",
        directory);

    ROM_std_hlpod_output_add_calc_time(
        diag->h_max,
        t,
        "wall_h_n_max.txt",
        directory);
}

void output_wall_gamma_sweep_penalty_decomposition_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_mu_by_gamma,
    const double* L_penalty_mu_by_gamma,
    const double* D_penalty_dt_by_gamma,
    const double* L_penalty_dt_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        if (D_penalty_mu_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_mu",
                "drag",
                g,
                D_penalty_mu_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_mu_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_mu",
                "lift",
                g,
                L_penalty_mu_by_gamma[g],
                t,
                directory);
        }
        if (D_penalty_dt_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_dt",
                "drag",
                g,
                D_penalty_dt_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_dt_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty_dt",
                "lift",
                g,
                L_penalty_dt_by_gamma[g],
                t,
                directory);
        }

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}

void output_wall_gamma_sweep_split_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_by_gamma,
    const double* L_penalty_by_gamma,
    const double* D_tangent_by_gamma,
    const double* L_tangent_by_gamma,
    const double* D_total_by_gamma,
    const double* L_total_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        if (D_penalty_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty",
                "drag",
                g,
                D_penalty_by_gamma[g],
                t,
                directory);
        }
        if (L_penalty_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "penalty",
                "lift",
                g,
                L_penalty_by_gamma[g],
                t,
                directory);
        }
        if (D_tangent_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "tangent",
                "drag",
                g,
                D_tangent_by_gamma[g],
                t,
                directory);
        }
        if (L_tangent_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "tangent",
                "lift",
                g,
                L_tangent_by_gamma[g],
                t,
                directory);
        }
        if (D_total_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "total",
                "drag",
                g,
                D_total_by_gamma[g],
                t,
                directory);
        }
        if (L_total_by_gamma != NULL) {
            t4wd_output_gamma_value(
                "total",
                "lift",
                g,
                L_total_by_gamma[g],
                t,
                directory);
        }

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}

void output_wall_gamma_sweep_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_nitsche_by_gamma,
    const double* L_nitsche_by_gamma,
    double t,
    const char* directory)
{
    for (int g = 0; g < num_gamma; ++g) {
        char fname[256];

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_drag_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            D_nitsche_by_gamma[g],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_lift_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            L_nitsche_by_gamma[g],
            t,
            fname,
            directory);

        snprintf(
            fname,
            sizeof(fname),
            "wall_nitsche_gamma_%02d.txt",
            g);
        ROM_std_hlpod_output_add_calc_time(
            gamma_list[g],
            t,
            fname,
            directory);
    }
}




/*
 * Global maximum wall-error and |div u| location diagnostics
 * with local TET4/TRI3 length and quality information.
 *
 * Replace the previously added maximum-location diagnostic block in
 * fluid_sups_core/core_FOM.c with this entire block.
 *
 * Required helpers already present in core_FOM.c:
 *
 *   t4_make_owner_map()
 *   t4_wmd_surface_to_volume_nodes()
 *   t4_tri3_normal()
 *   t4wd_dot3()
 *   t4wd_norm3()
 *   t4_wmd_dot3()
 *   t4_wmd_norm3()
 *   t4_wmd_cross3()
 *   t4_wmd_tet4_geometry()
 *   t4_wmd_tet4_quality()
 *
 * Output element/node IDs are partition-local IDs.  The MPI rank and
 * coordinates are written together with the IDs.
 */

typedef struct {
    double tet_quality;

    double min_edge_length;
    double max_edge_length;
    double edge_length_ratio;

    double min_altitude;
    double max_altitude;
    double altitude_ratio;

    /*
     * Frobenius condition number:
     *
     *     ||J||_F ||J^{-1}||_F
     */
    double jacobian_condition;
} T4TetShapeMetrics;


typedef struct {
    double face_area;
    double tet_volume;
    double detJ;

    /*
     * Nitsche penalty length:
     * distance from owner-TET centroid to wall-face plane.
     */
    double h_n;

    /*
     * Full opposite-node-to-face-plane height.
     * For a linear TET4, nominally face_height = 4*h_n.
     */
    double face_height;

    /*
     * Equal-area equilateral-triangle tangential scale.
     */
    double h_t;

    double aspect_ratio_ht_over_hn;
    double nonorthogonality_deg;

    T4TetShapeMetrics tet_shape;

    double tet_node_x[4][3];
} T4WallFaceMetrics;


typedef struct {
    int valid;
    double value;

    int rank;
    int surf_elem;
    int owner_elem;
    int owner_lface;
    int quadrature_point;

    int surf_nodes[3];
    int tet_nodes[4];

    double x[3];
    double n[3];
    double u[3];
    double r[3];

    T4WallFaceMetrics mesh;
} T4WallErrorLocation;


typedef struct {
    int valid;
    double abs_div_u;
    double div_u;

    int rank;
    int elem;
    int tet_nodes[4];

    double x[3];
    double volume;
    double detJ;

    T4TetShapeMetrics tet_shape;
    double characteristic_length;

    double node_x[4][3];
    double nodal_velocity[4][3];
} T4DivULocation;


static int t4loc_all_ranks_ok(
    int local_status,
    MONOLIS_COM* monolis_com,
    const char* label)
{
    if (monolis_com == NULL) {
        fprintf(
            stderr,
            "[rank %d][%s] monolis_com is NULL\n",
            monolis_mpi_get_global_my_rank(),
            label != NULL ? label : "location diagnostic");

        return 0;
    }

    double failure_count =
        (local_status == 0)
        ? 0.0
        : 1.0;

    monolis_allreduce_R(
        1,
        &failure_count,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    if (failure_count > 0.5) {
        if (
            local_status != 0 &&
            monolis_mpi_get_global_my_rank() >= 0) {

            fprintf(
                stderr,
                "[rank %d][%s] local calculation failed: status=%d\n",
                monolis_mpi_get_global_my_rank(),
                label != NULL ? label : "location diagnostic",
                local_status);
        }

        if (monolis_mpi_get_global_my_rank() == 0) {
            fprintf(
                stderr,
                "[%s] skipped because %.0f rank(s) failed\n",
                label != NULL ? label : "location diagnostic",
                failure_count);
        }

        return 0;
    }

    return 1;
}


static int t4loc_select_winner_rank(
    int local_valid,
    double local_value,
    MONOLIS_COM* monolis_com,
    double* global_value,
    int* winner_rank)
{
    if (
        monolis_com == NULL ||
        global_value == NULL ||
        winner_rank == NULL) {

        return 0;
    }

    double g =
        local_valid
        ? local_value
        : -1.0;

    monolis_allreduce_R(
        1,
        &g,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    if (
        !isfinite(g) ||
        g < 0.0) {

        *global_value = 0.0;
        *winner_rank = -1;
        return 0;
    }

    const double tol =
        1.0e-12
        * fmax(
            1.0,
            fabs(g));

    double rank_candidate = -1.0;

    if (
        local_valid &&
        fabs(local_value - g) <= tol) {

        rank_candidate =
            (double)monolis_mpi_get_global_my_rank();
    }

    /*
     * If equal maxima occur on multiple ranks, select the largest rank.
     * This guarantees that exactly one rank contributes the location pack.
     */
    monolis_allreduce_R(
        1,
        &rank_candidate,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    *global_value = g;
    *winner_rank =
        (int)llround(rank_candidate);

    return
        (*winner_rank >= 0);
}


static FILE* t4loc_open_append_with_header(
    const char* directory,
    const char* filename,
    const char* header)
{
    if (filename == NULL) {
        return NULL;
    }

    char path[4096];

    if (
        directory == NULL ||
        directory[0] == '\0' ||
        strcmp(directory, ".") == 0) {

        snprintf(
            path,
            sizeof(path),
            "%s",
            filename);

    } else {
        const size_t n =
            strlen(directory);

        snprintf(
            path,
            sizeof(path),
            "%s%s%s",
            directory,
            (
                n > 0 &&
                directory[n - 1] == '/'
            )
            ? ""
            : "/",
            filename);
    }

    FILE* fp =
        fopen(
            path,
            "a+");

    if (fp == NULL) {
        fprintf(
            stderr,
            "[location diagnostic] cannot open %s\n",
            path);

        return NULL;
    }

    if (
        fseek(
            fp,
            0,
            SEEK_END) == 0) {

        const long size =
            ftell(fp);

        if (
            size == 0 &&
            header != NULL) {

            fprintf(
                fp,
                "%s\n",
                header);
        }
    }

    return fp;
}



static double t4loc_frobenius_norm3x3(
    const double A[3][3])
{
    double sum = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            sum += A[i][j] * A[i][j];
        }
    }

    return sqrt(sum);
}


static int t4loc_compute_tet_shape_metrics(
    const BBFE_DATA* fe,
    int ke,
    double volume,
    const double J[3][3],
    const double invJ[3][3],
    T4TetShapeMetrics* metrics)
{
    if (
        fe == NULL ||
        metrics == NULL ||
        ke < 0 ||
        ke >= fe->total_num_elems ||
        !isfinite(volume) ||
        volume <= 1.0e-300) {

        return 0;
    }

    memset(
        metrics,
        0,
        sizeof(*metrics));

    metrics->tet_quality =
        t4_wmd_tet4_quality(
            fe,
            ke,
            volume);

    metrics->min_edge_length = 1.0e300;
    metrics->max_edge_length = 0.0;

    for (int a = 0; a < 4; ++a) {
        const int gida =
            fe->conn[ke][a];

        for (int b = a + 1; b < 4; ++b) {
            const int gidb =
                fe->conn[ke][b];

            const double dx =
                fe->x[gidb][0]
              - fe->x[gida][0];

            const double dy =
                fe->x[gidb][1]
              - fe->x[gida][1];

            const double dz =
                fe->x[gidb][2]
              - fe->x[gida][2];

            const double edge =
                sqrt(
                    dx*dx
                  + dy*dy
                  + dz*dz);

            if (
                !isfinite(edge) ||
                edge <= 1.0e-300) {

                return 0;
            }

            metrics->min_edge_length =
                fmin(
                    metrics->min_edge_length,
                    edge);

            metrics->max_edge_length =
                fmax(
                    metrics->max_edge_length,
                    edge);
        }
    }

    metrics->edge_length_ratio =
        metrics->max_edge_length
        / fmax(
            metrics->min_edge_length,
            1.0e-300);

    /*
     * TET4 faces, each row is the face opposite local node i.
     */
    static const int face_nodes[4][3] = {
        {1, 2, 3},
        {0, 3, 2},
        {0, 1, 3},
        {0, 2, 1}
    };

    metrics->min_altitude = 1.0e300;
    metrics->max_altitude = 0.0;

    for (int f = 0; f < 4; ++f) {
        const int g0 =
            fe->conn[ke][face_nodes[f][0]];

        const int g1 =
            fe->conn[ke][face_nodes[f][1]];

        const int g2 =
            fe->conn[ke][face_nodes[f][2]];

        const double e01[3] = {
            fe->x[g1][0] - fe->x[g0][0],
            fe->x[g1][1] - fe->x[g0][1],
            fe->x[g1][2] - fe->x[g0][2]
        };

        const double e02[3] = {
            fe->x[g2][0] - fe->x[g0][0],
            fe->x[g2][1] - fe->x[g0][1],
            fe->x[g2][2] - fe->x[g0][2]
        };

        double cross[3];

        t4_wmd_cross3(
            cross,
            e01,
            e02);

        const double face_area =
            0.5
            * t4_wmd_norm3(
                cross);

        if (
            !isfinite(face_area) ||
            face_area <= 1.0e-300) {

            return 0;
        }

        /*
         * V = A*h/3.
         */
        const double altitude =
            3.0 * volume
            / face_area;

        if (
            !isfinite(altitude) ||
            altitude <= 1.0e-300) {

            return 0;
        }

        metrics->min_altitude =
            fmin(
                metrics->min_altitude,
                altitude);

        metrics->max_altitude =
            fmax(
                metrics->max_altitude,
                altitude);
    }

    metrics->altitude_ratio =
        metrics->max_altitude
        / fmax(
            metrics->min_altitude,
            1.0e-300);

    metrics->jacobian_condition =
        t4loc_frobenius_norm3x3(J)
        * t4loc_frobenius_norm3x3(invJ);

    if (
        !isfinite(metrics->tet_quality) ||
        !isfinite(metrics->edge_length_ratio) ||
        !isfinite(metrics->altitude_ratio) ||
        !isfinite(metrics->jacobian_condition)) {

        return 0;
    }

    return 1;
}


static void t4loc_init_wall_location(
    T4WallErrorLocation* loc)
{
    if (loc == NULL) {
        return;
    }

    memset(
        loc,
        0,
        sizeof(*loc));

    loc->valid = 0;
    loc->value = -1.0;

    loc->rank = -1;
    loc->surf_elem = -1;
    loc->owner_elem = -1;
    loc->owner_lface = -1;
    loc->quadrature_point = -1;

    for (int a = 0; a < 3; ++a) {
        loc->surf_nodes[a] = -1;
    }

    for (int a = 0; a < 4; ++a) {
        loc->tet_nodes[a] = -1;
    }
}


static void t4loc_store_wall_location(
    T4WallErrorLocation* loc,
    double value,
    int es,
    int ke,
    int local_face,
    int p,
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const double x_q[3],
    const double n[3],
    const double u_q[3],
    const double r_q[3],
    const T4WallFaceMetrics* mesh)
{
    if (
        loc == NULL ||
        surf == NULL ||
        fe == NULL ||
        x_q == NULL ||
        n == NULL ||
        u_q == NULL ||
        r_q == NULL ||
        mesh == NULL) {

        return;
    }

    if (
        loc->valid &&
        value <= loc->value) {

        return;
    }

    loc->valid = 1;
    loc->value = value;

    loc->rank =
        monolis_mpi_get_global_my_rank();

    loc->surf_elem = es;
    loc->owner_elem = ke;
    loc->owner_lface = local_face;
    loc->quadrature_point = p;

    for (int a = 0; a < 3; ++a) {
        loc->surf_nodes[a] =
            surf->conn[es][a];

        loc->x[a] = x_q[a];
        loc->n[a] = n[a];
        loc->u[a] = u_q[a];
        loc->r[a] = r_q[a];
    }

    for (int a = 0; a < 4; ++a) {
        loc->tet_nodes[a] =
            fe->conn[ke][a];
    }

    loc->mesh = *mesh;
}


static int t4loc_reduce_wall_location(
    const T4WallErrorLocation* local,
    MONOLIS_COM* monolis_com,
    T4WallErrorLocation* global)
{
    if (
        local == NULL ||
        monolis_com == NULL ||
        global == NULL) {

        return -1;
    }

    t4loc_init_wall_location(global);

    double global_value = 0.0;
    int winner_rank = -1;

    if (!t4loc_select_winner_rank(
            local->valid,
            local->value,
            monolis_com,
            &global_value,
            &winner_rank)) {

        return 0;
    }

    enum {
        T4LOC_WALL_PACK_SIZE = 51
    };

    double pack[T4LOC_WALL_PACK_SIZE];

    for (int i = 0;
         i < T4LOC_WALL_PACK_SIZE;
         ++i) {

        pack[i] = 0.0;
    }

    if (
        local->valid &&
        monolis_mpi_get_global_my_rank() == winner_rank) {

        int k = 0;

        pack[k++] = (double)local->surf_elem;
        pack[k++] = (double)local->owner_elem;
        pack[k++] = (double)local->owner_lface;
        pack[k++] = (double)local->quadrature_point;

        for (int a = 0; a < 3; ++a) {
            pack[k++] =
                (double)local->surf_nodes[a];
        }

        for (int a = 0; a < 4; ++a) {
            pack[k++] =
                (double)local->tet_nodes[a];
        }

        for (int a = 0; a < 3; ++a) {
            pack[k++] = local->x[a];
        }

        for (int a = 0; a < 3; ++a) {
            pack[k++] = local->n[a];
        }

        for (int a = 0; a < 3; ++a) {
            pack[k++] = local->u[a];
        }

        for (int a = 0; a < 3; ++a) {
            pack[k++] = local->r[a];
        }

        pack[k++] = local->mesh.face_area;
        pack[k++] = local->mesh.tet_volume;
        pack[k++] = local->mesh.detJ;
        pack[k++] = local->mesh.h_n;
        pack[k++] = local->mesh.face_height;
        pack[k++] = local->mesh.h_t;
        pack[k++] = local->mesh.aspect_ratio_ht_over_hn;
        pack[k++] = local->mesh.nonorthogonality_deg;

        pack[k++] = local->mesh.tet_shape.tet_quality;
        pack[k++] = local->mesh.tet_shape.min_edge_length;
        pack[k++] = local->mesh.tet_shape.max_edge_length;
        pack[k++] = local->mesh.tet_shape.edge_length_ratio;
        pack[k++] = local->mesh.tet_shape.min_altitude;
        pack[k++] = local->mesh.tet_shape.max_altitude;
        pack[k++] = local->mesh.tet_shape.altitude_ratio;
        pack[k++] = local->mesh.tet_shape.jacobian_condition;

        for (int a = 0; a < 4; ++a) {
            for (int i = 0; i < 3; ++i) {
                pack[k++] =
                    local->mesh.tet_node_x[a][i];
            }
        }

        if (k != T4LOC_WALL_PACK_SIZE) {
            fprintf(
                stderr,
                "[rank %d][wall location] "
                "internal pack-size error: k=%d expected=%d\n",
                monolis_mpi_get_global_my_rank(),
                k,
                T4LOC_WALL_PACK_SIZE);

            return -1;
        }
    }

    monolis_allreduce_R(
        T4LOC_WALL_PACK_SIZE,
        pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    int k = 0;

    global->valid = 1;
    global->value = global_value;
    global->rank = winner_rank;

    global->surf_elem =
        (int)llround(pack[k++]);

    global->owner_elem =
        (int)llround(pack[k++]);

    global->owner_lface =
        (int)llround(pack[k++]);

    global->quadrature_point =
        (int)llround(pack[k++]);

    for (int a = 0; a < 3; ++a) {
        global->surf_nodes[a] =
            (int)llround(pack[k++]);
    }

    for (int a = 0; a < 4; ++a) {
        global->tet_nodes[a] =
            (int)llround(pack[k++]);
    }

    for (int a = 0; a < 3; ++a) {
        global->x[a] = pack[k++];
    }

    for (int a = 0; a < 3; ++a) {
        global->n[a] = pack[k++];
    }

    for (int a = 0; a < 3; ++a) {
        global->u[a] = pack[k++];
    }

    for (int a = 0; a < 3; ++a) {
        global->r[a] = pack[k++];
    }

    global->mesh.face_area = pack[k++];
    global->mesh.tet_volume = pack[k++];
    global->mesh.detJ = pack[k++];
    global->mesh.h_n = pack[k++];
    global->mesh.face_height = pack[k++];
    global->mesh.h_t = pack[k++];
    global->mesh.aspect_ratio_ht_over_hn = pack[k++];
    global->mesh.nonorthogonality_deg = pack[k++];

    global->mesh.tet_shape.tet_quality = pack[k++];
    global->mesh.tet_shape.min_edge_length = pack[k++];
    global->mesh.tet_shape.max_edge_length = pack[k++];
    global->mesh.tet_shape.edge_length_ratio = pack[k++];
    global->mesh.tet_shape.min_altitude = pack[k++];
    global->mesh.tet_shape.max_altitude = pack[k++];
    global->mesh.tet_shape.altitude_ratio = pack[k++];
    global->mesh.tet_shape.jacobian_condition = pack[k++];

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            global->mesh.tet_node_x[a][i] =
                pack[k++];
        }
    }

    if (k != T4LOC_WALL_PACK_SIZE) {
        return -1;
    }

    return 0;
}


static int t4loc_write_wall_location(
    const T4WallErrorLocation* loc,
    double t,
    const char* directory,
    const char* filename)
{
    if (
        loc == NULL ||
        !loc->valid) {

        return 0;
    }

    FILE* fp =
        t4loc_open_append_with_header(
            directory,
            filename,
            "# t value rank surf_elem owner_tet owner_lface qp "
            "x y z nx ny nz ux uy uz rx ry rz "
            "face_area tet_volume detJ h_n face_height h_t "
            "aspect_ht_over_hn nonorth_deg tet_quality "
            "min_edge max_edge edge_ratio "
            "min_altitude max_altitude altitude_ratio jacobian_condition "
            "surf_n0 surf_n1 surf_n2 "
            "tet_n0 tet_n1 tet_n2 tet_n3 "
            "tet0_x tet0_y tet0_z "
            "tet1_x tet1_y tet1_z "
            "tet2_x tet2_y tet2_z "
            "tet3_x tet3_y tet3_z");

    if (fp == NULL) {
        return -1;
    }

    fprintf(
        fp,
        "%.15e %.15e "
        "%d %d %d %d %d "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e %.15e "
        "%d %d %d "
        "%d %d %d %d "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e\n",
        t,
        loc->value,
        loc->rank,
        loc->surf_elem,
        loc->owner_elem,
        loc->owner_lface,
        loc->quadrature_point,
        loc->x[0],
        loc->x[1],
        loc->x[2],
        loc->n[0],
        loc->n[1],
        loc->n[2],
        loc->u[0],
        loc->u[1],
        loc->u[2],
        loc->r[0],
        loc->r[1],
        loc->r[2],
        loc->mesh.face_area,
        loc->mesh.tet_volume,
        loc->mesh.detJ,
        loc->mesh.h_n,
        loc->mesh.face_height,
        loc->mesh.h_t,
        loc->mesh.aspect_ratio_ht_over_hn,
        loc->mesh.nonorthogonality_deg,
        loc->mesh.tet_shape.tet_quality,
        loc->mesh.tet_shape.min_edge_length,
        loc->mesh.tet_shape.max_edge_length,
        loc->mesh.tet_shape.edge_length_ratio,
        loc->mesh.tet_shape.min_altitude,
        loc->mesh.tet_shape.max_altitude,
        loc->mesh.tet_shape.altitude_ratio,
        loc->mesh.tet_shape.jacobian_condition,
        loc->surf_nodes[0],
        loc->surf_nodes[1],
        loc->surf_nodes[2],
        loc->tet_nodes[0],
        loc->tet_nodes[1],
        loc->tet_nodes[2],
        loc->tet_nodes[3],
        loc->mesh.tet_node_x[0][0],
        loc->mesh.tet_node_x[0][1],
        loc->mesh.tet_node_x[0][2],
        loc->mesh.tet_node_x[1][0],
        loc->mesh.tet_node_x[1][1],
        loc->mesh.tet_node_x[1][2],
        loc->mesh.tet_node_x[2][0],
        loc->mesh.tet_node_x[2][1],
        loc->mesh.tet_node_x[2][2],
        loc->mesh.tet_node_x[3][0],
        loc->mesh.tet_node_x[3][1],
        loc->mesh.tet_node_x[3][2]);

    fclose(fp);

    return 0;
}


/*
 * Locate both
 *
 *   max |(u-Uw).n| / Uref
 *   max ||u-Uw|| / Uref
 *
 * and write the owner-TET/face length and quality information.
 */
int output_wall_max_error_locations_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double Uref,
    const double Uw_in[3],
    MONOLIS_COM* monolis_com,
    double t,
    const char* directory)
{
    int local_status = 0;

    if (
        surf == NULL ||
        fe == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        monolis_com == NULL) {

        local_status = -1;
    }

    int* owner = NULL;
    int* lface = NULL;

    if (
        local_status == 0 &&
        surf->total_num_elems < 0) {

        local_status = -1;
    }

    if (
        local_status == 0 &&
        surf->total_num_elems > 0 &&
        (
            surf->local_num_nodes != 3 ||
            fe->local_num_nodes != 4
        )) {

        fprintf(
            stderr,
            "[wall max location] "
            "TET4/TRI3 expected, got fe=%d surf=%d ne_surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        local_status = -1;
    }

    if (
        local_status == 0 &&
        surf->total_num_elems > 0) {

        if (
            t4_make_owner_map(
                surf,
                fe,
                &owner,
                &lface) != 0) {

            local_status = -1;
        }
    }

    if (!t4loc_all_ranks_ok(
            local_status,
            monolis_com,
            "wall max error location setup")) {

        free(owner);
        free(lface);
        return -1;
    }

    T4WallErrorLocation local_abs_un;
    T4WallErrorLocation local_norm_r;

    t4loc_init_wall_location(
        &local_abs_un);

    t4loc_init_wall_location(
        &local_norm_r);

    const double U_scale =
        fmax(
            fabs(Uref),
            1.0e-300);

    const double Uw[3] = {
        Uw_in != NULL ? Uw_in[0] : 0.0,
        Uw_in != NULL ? Uw_in[1] : 0.0,
        Uw_in != NULL ? Uw_in[2] : 0.0
    };

    if (surf->total_num_elems > 0) {
        for (int es = 0;
             es < surf->total_num_elems;
             ++es) {

            const int ke =
                owner[es];

            if (ke < 0) {
                continue;
            }

            int s2v[3];
            int opposite_local_node = -1;

            if (!t4_wmd_surface_to_volume_nodes(
                    surf,
                    fe,
                    es,
                    ke,
                    s2v,
                    &opposite_local_node)) {

                continue;
            }

            double xtri[3][3];

            for (int a = 0; a < 3; ++a) {
                const int gid =
                    surf->conn[es][a];

                xtri[a][0] = fe->x[gid][0];
                xtri[a][1] = fe->x[gid][1];
                xtri[a][2] = fe->x[gid][2];
            }

            const double edge01[3] = {
                xtri[1][0] - xtri[0][0],
                xtri[1][1] - xtri[0][1],
                xtri[1][2] - xtri[0][2]
            };

            const double edge02[3] = {
                xtri[2][0] - xtri[0][0],
                xtri[2][1] - xtri[0][1],
                xtri[2][2] - xtri[0][2]
            };

            double face_cross[3];

            t4_wmd_cross3(
                face_cross,
                edge01,
                edge02);

            const double cross_norm =
                t4_wmd_norm3(
                    face_cross);

            const double face_area =
                0.5 * cross_norm;

            if (
                !isfinite(face_area) ||
                face_area <= 1.0e-300) {

                continue;
            }

            const double xf[3] = {
                (
                    xtri[0][0]
                  + xtri[1][0]
                  + xtri[2][0]
                ) / 3.0,
                (
                    xtri[0][1]
                  + xtri[1][1]
                  + xtri[2][1]
                ) / 3.0,
                (
                    xtri[0][2]
                  + xtri[1][2]
                  + xtri[2][2]
                ) / 3.0
            };

            double xc[3] = {
                0.0,
                0.0,
                0.0
            };

            T4WallFaceMetrics mesh;

            memset(
                &mesh,
                0,
                sizeof(mesh));

            for (int a = 0; a < 4; ++a) {
                const int gid =
                    fe->conn[ke][a];

                xc[0] +=
                    0.25 * fe->x[gid][0];

                xc[1] +=
                    0.25 * fe->x[gid][1];

                xc[2] +=
                    0.25 * fe->x[gid][2];

                mesh.tet_node_x[a][0] =
                    fe->x[gid][0];

                mesh.tet_node_x[a][1] =
                    fe->x[gid][1];

                mesh.tet_node_x[a][2] =
                    fe->x[gid][2];
            }

            const double centroid_to_face[3] = {
                xf[0] - xc[0],
                xf[1] - xc[1],
                xf[2] - xc[2]
            };

            double face_n[3] = {
                face_cross[0] / cross_norm,
                face_cross[1] / cross_norm,
                face_cross[2] / cross_norm
            };

            if (
                t4_wmd_dot3(
                    centroid_to_face,
                    face_n) < 0.0) {

                face_n[0] *= -1.0;
                face_n[1] *= -1.0;
                face_n[2] *= -1.0;
            }

            double J[3][3];
            double invJ[3][3];
            double dN_dx[4][3];
            double detJ = 0.0;
            double volume = 0.0;

            if (!t4_wmd_tet4_geometry(
                    fe,
                    ke,
                    J,
                    invJ,
                    dN_dx,
                    &detJ,
                    &volume)) {

                continue;
            }

            if (
                !isfinite(volume) ||
                volume <= 1.0e-300) {

                continue;
            }

            mesh.face_area = face_area;
            mesh.tet_volume = volume;
            mesh.detJ = detJ;

            mesh.h_n =
                fabs(
                    t4_wmd_dot3(
                        centroid_to_face,
                        face_n));

            if (mesh.h_n <= 1.0e-12) {
                mesh.h_n = 1.0e-12;
            }

            const int opposite_gid =
                fe->conn[ke][opposite_local_node];

            const double xopp_to_face[3] = {
                fe->x[opposite_gid][0] - xf[0],
                fe->x[opposite_gid][1] - xf[1],
                fe->x[opposite_gid][2] - xf[2]
            };

            mesh.face_height =
                fabs(
                    t4_wmd_dot3(
                        xopp_to_face,
                        face_n));

            mesh.h_t =
                sqrt(
                    4.0 * face_area
                    / sqrt(3.0));

            mesh.aspect_ratio_ht_over_hn =
                mesh.h_t
                / fmax(
                    mesh.h_n,
                    1.0e-300);

            const double d_centroid =
                t4_wmd_norm3(
                    centroid_to_face);

            mesh.nonorthogonality_deg = 0.0;

            if (d_centroid > 1.0e-300) {
                double cos_theta =
                    fabs(
                        t4_wmd_dot3(
                            centroid_to_face,
                            face_n))
                    / d_centroid;

                cos_theta =
                    fmax(
                        0.0,
                        fmin(
                            1.0,
                            cos_theta));

                mesh.nonorthogonality_deg =
                    acos(cos_theta)
                    * 180.0
                    / T4_WMD_PI;
            }

            if (!t4loc_compute_tet_shape_metrics(
                    fe,
                    ke,
                    volume,
                    J,
                    invJ,
                    &mesh.tet_shape)) {

                continue;
            }

            for (int p = 0;
                 p < basis_surf->num_integ_points;
                 ++p) {

                double nraw[3];

                t4_tri3_normal(
                    xtri,
                    basis_surf->dN_dxi[p],
                    basis_surf->dN_det[p],
                    nraw);

                const double Jraw =
                    t4wd_norm3(
                        nraw);

                if (
                    !isfinite(Jraw) ||
                    Jraw <= 1.0e-300) {

                    continue;
                }

                double n[3] = {
                    nraw[0] / Jraw,
                    nraw[1] / Jraw,
                    nraw[2] / Jraw
                };

                if (
                    t4wd_dot3(
                        centroid_to_face,
                        n) < 0.0) {

                    n[0] *= -1.0;
                    n[1] *= -1.0;
                    n[2] *= -1.0;
                }

                double x_q[3] = {
                    0.0,
                    0.0,
                    0.0
                };

                double u_q[3] = {
                    0.0,
                    0.0,
                    0.0
                };

                for (int a = 0; a < 3; ++a) {
                    const int gid =
                        surf->conn[es][a];

                    const double Na =
                        basis_surf->N[p][a];

                    x_q[0] +=
                        Na * fe->x[gid][0];

                    x_q[1] +=
                        Na * fe->x[gid][1];

                    x_q[2] +=
                        Na * fe->x[gid][2];

                    u_q[0] +=
                        Na * vals->v[gid][0];

                    u_q[1] +=
                        Na * vals->v[gid][1];

                    u_q[2] +=
                        Na * vals->v[gid][2];
                }

                const double r_q[3] = {
                    u_q[0] - Uw[0],
                    u_q[1] - Uw[1],
                    u_q[2] - Uw[2]
                };

                const double un =
                    t4wd_dot3(
                        r_q,
                        n);

                const double abs_un_over_Uref =
                    fabs(un)
                    / U_scale;

                const double norm_r_over_Uref =
                    t4wd_norm3(
                        r_q)
                    / U_scale;

                t4loc_store_wall_location(
                    &local_abs_un,
                    abs_un_over_Uref,
                    es,
                    ke,
                    lface[es],
                    p,
                    surf,
                    fe,
                    x_q,
                    n,
                    u_q,
                    r_q,
                    &mesh);

                t4loc_store_wall_location(
                    &local_norm_r,
                    norm_r_over_Uref,
                    es,
                    ke,
                    lface[es],
                    p,
                    surf,
                    fe,
                    x_q,
                    n,
                    u_q,
                    r_q,
                    &mesh);
            }
        }
    }

    free(owner);
    free(lface);

    T4WallErrorLocation global_abs_un;
    T4WallErrorLocation global_norm_r;

    if (
        t4loc_reduce_wall_location(
            &local_abs_un,
            monolis_com,
            &global_abs_un) != 0) {

        local_status = -1;
    }

    if (
        t4loc_reduce_wall_location(
            &local_norm_r,
            monolis_com,
            &global_norm_r) != 0) {

        local_status = -1;
    }

    if (!t4loc_all_ranks_ok(
            local_status,
            monolis_com,
            "wall max error location reduction")) {

        return -1;
    }

    int output_status = 0;

    if (
        monolis_mpi_get_global_my_rank() == 0) {

        if (
            t4loc_write_wall_location(
                &global_abs_un,
                t,
                directory,
                "wall_max_abs_un_location.txt") != 0) {

            output_status = -1;
        }

        if (
            t4loc_write_wall_location(
                &global_norm_r,
                t,
                directory,
                "wall_max_norm_u_minus_Uw_location.txt") != 0) {

            output_status = -1;
        }
    }

    if (!t4loc_all_ranks_ok(
            output_status,
            monolis_com,
            "wall max error location output")) {

        return -1;
    }

    return 0;
}


static void t4loc_init_div_location(
    T4DivULocation* loc)
{
    if (loc == NULL) {
        return;
    }

    memset(
        loc,
        0,
        sizeof(*loc));

    loc->valid = 0;
    loc->abs_div_u = -1.0;
    loc->div_u = 0.0;
    loc->rank = -1;
    loc->elem = -1;

    for (int a = 0; a < 4; ++a) {
        loc->tet_nodes[a] = -1;
    }
}


static int t4loc_reduce_div_location(
    const T4DivULocation* local,
    MONOLIS_COM* monolis_com,
    T4DivULocation* global)
{
    if (
        local == NULL ||
        monolis_com == NULL ||
        global == NULL) {

        return -1;
    }

    t4loc_init_div_location(
        global);

    double global_value = 0.0;
    int winner_rank = -1;

    if (!t4loc_select_winner_rank(
            local->valid,
            local->abs_div_u,
            monolis_com,
            &global_value,
            &winner_rank)) {

        return 0;
    }

    enum {
        T4LOC_DIV_PACK_SIZE = 44
    };

    double pack[T4LOC_DIV_PACK_SIZE];

    for (int i = 0;
         i < T4LOC_DIV_PACK_SIZE;
         ++i) {

        pack[i] = 0.0;
    }

    if (
        local->valid &&
        monolis_mpi_get_global_my_rank() == winner_rank) {

        int k = 0;

        pack[k++] =
            local->div_u;

        pack[k++] =
            (double)local->elem;

        for (int a = 0; a < 4; ++a) {
            pack[k++] =
                (double)local->tet_nodes[a];
        }

        for (int a = 0; a < 3; ++a) {
            pack[k++] =
                local->x[a];
        }

        pack[k++] =
            local->volume;

        pack[k++] =
            local->detJ;

        pack[k++] =
            local->tet_shape.tet_quality;

        pack[k++] =
            local->tet_shape.min_edge_length;

        pack[k++] =
            local->tet_shape.max_edge_length;

        pack[k++] =
            local->tet_shape.edge_length_ratio;

        pack[k++] =
            local->tet_shape.min_altitude;

        pack[k++] =
            local->tet_shape.max_altitude;

        pack[k++] =
            local->tet_shape.altitude_ratio;

        pack[k++] =
            local->tet_shape.jacobian_condition;

        pack[k++] =
            local->characteristic_length;

        for (int a = 0; a < 4; ++a) {
            for (int i = 0; i < 3; ++i) {
                pack[k++] =
                    local->node_x[a][i];
            }
        }

        for (int a = 0; a < 4; ++a) {
            for (int i = 0; i < 3; ++i) {
                pack[k++] =
                    local->nodal_velocity[a][i];
            }
        }

        if (k != T4LOC_DIV_PACK_SIZE) {
            fprintf(
                stderr,
                "[rank %d][div location] "
                "internal pack-size error: k=%d expected=%d\n",
                monolis_mpi_get_global_my_rank(),
                k,
                T4LOC_DIV_PACK_SIZE);

            return -1;
        }
    }

    monolis_allreduce_R(
        T4LOC_DIV_PACK_SIZE,
        pack,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    int k = 0;

    global->valid = 1;
    global->abs_div_u = global_value;
    global->div_u = pack[k++];
    global->rank = winner_rank;

    global->elem =
        (int)llround(pack[k++]);

    for (int a = 0; a < 4; ++a) {
        global->tet_nodes[a] =
            (int)llround(pack[k++]);
    }

    for (int a = 0; a < 3; ++a) {
        global->x[a] =
            pack[k++];
    }

    global->volume =
        pack[k++];

    global->detJ =
        pack[k++];

    global->tet_shape.tet_quality =
        pack[k++];

    global->tet_shape.min_edge_length =
        pack[k++];

    global->tet_shape.max_edge_length =
        pack[k++];

    global->tet_shape.edge_length_ratio =
        pack[k++];

    global->tet_shape.min_altitude =
        pack[k++];

    global->tet_shape.max_altitude =
        pack[k++];

    global->tet_shape.altitude_ratio =
        pack[k++];

    global->tet_shape.jacobian_condition =
        pack[k++];

    global->characteristic_length =
        pack[k++];

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            global->node_x[a][i] =
                pack[k++];
        }
    }

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 3; ++i) {
            global->nodal_velocity[a][i] =
                pack[k++];
        }
    }

    if (k != T4LOC_DIV_PACK_SIZE) {
        return -1;
    }

    return 0;
}


static int t4loc_write_div_location(
    const T4DivULocation* loc,
    double t,
    const char* directory)
{
    if (
        loc == NULL ||
        !loc->valid) {

        return 0;
    }

    FILE* fp =
        t4loc_open_append_with_header(
            directory,
            "domain_max_abs_div_u_location.txt",
            "# t abs_div_u div_u rank elem "
            "x y z volume detJ "
            "tet_quality min_edge max_edge edge_ratio "
            "min_altitude max_altitude altitude_ratio "
            "characteristic_length jacobian_condition "
            "n0 n1 n2 n3 "
            "x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3 "
            "u0x u0y u0z u1x u1y u1z "
            "u2x u2y u2z u3x u3y u3z");

    if (fp == NULL) {
        return -1;
    }

    fprintf(
        fp,
        "%.15e %.15e %.15e "
        "%d %d "
        "%.15e %.15e %.15e "
        "%.15e %.15e "
        "%.15e %.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e "
        "%d %d %d %d "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e\n",
        t,
        loc->abs_div_u,
        loc->div_u,
        loc->rank,
        loc->elem,
        loc->x[0],
        loc->x[1],
        loc->x[2],
        loc->volume,
        loc->detJ,
        loc->tet_shape.tet_quality,
        loc->tet_shape.min_edge_length,
        loc->tet_shape.max_edge_length,
        loc->tet_shape.edge_length_ratio,
        loc->tet_shape.min_altitude,
        loc->tet_shape.max_altitude,
        loc->tet_shape.altitude_ratio,
        loc->characteristic_length,
        loc->tet_shape.jacobian_condition,
        loc->tet_nodes[0],
        loc->tet_nodes[1],
        loc->tet_nodes[2],
        loc->tet_nodes[3],
        loc->node_x[0][0],
        loc->node_x[0][1],
        loc->node_x[0][2],
        loc->node_x[1][0],
        loc->node_x[1][1],
        loc->node_x[1][2],
        loc->node_x[2][0],
        loc->node_x[2][1],
        loc->node_x[2][2],
        loc->node_x[3][0],
        loc->node_x[3][1],
        loc->node_x[3][2],
        loc->nodal_velocity[0][0],
        loc->nodal_velocity[0][1],
        loc->nodal_velocity[0][2],
        loc->nodal_velocity[1][0],
        loc->nodal_velocity[1][1],
        loc->nodal_velocity[1][2],
        loc->nodal_velocity[2][0],
        loc->nodal_velocity[2][1],
        loc->nodal_velocity[2][2],
        loc->nodal_velocity[3][0],
        loc->nodal_velocity[3][1],
        loc->nodal_velocity[3][2]);

    fclose(fp);

    return 0;
}


/*
 * Locate the TET4 with the global maximum |div u| and write its
 * length/quality information.
 *
 * For TET4, grad(u) and div(u) are constant inside each element.
 */
int output_domain_max_div_u_location_tet4(
    const BBFE_DATA* fe,
    const VALUES* vals,
    MONOLIS_COM* monolis_com,
    double t,
    const char* directory)
{
    int local_status = 0;

    if (
        fe == NULL ||
        vals == NULL ||
        monolis_com == NULL) {

        local_status = -1;
    }

    if (
        local_status == 0 &&
        fe->local_num_nodes != 4) {

        fprintf(
            stderr,
            "[div-u max location] "
            "TET4 expected, got fe=%d\n",
            fe->local_num_nodes);

        local_status = -1;
    }

    if (!t4loc_all_ranks_ok(
            local_status,
            monolis_com,
            "maximum divergence location setup")) {

        return -1;
    }

    T4DivULocation local;
    T4DivULocation global;

    t4loc_init_div_location(
        &local);

    for (int ke = 0;
         ke < fe->total_num_elems;
         ++ke) {

        double J[3][3];
        double invJ[3][3];
        double dN_dx[4][3];
        double detJ = 0.0;
        double volume = 0.0;

        if (!t4_wmd_tet4_geometry(
                fe,
                ke,
                J,
                invJ,
                dN_dx,
                &detJ,
                &volume)) {

            continue;
        }

        if (
            !isfinite(volume) ||
            volume <= 1.0e-300) {

            continue;
        }

        T4TetShapeMetrics shape;

        if (!t4loc_compute_tet_shape_metrics(
                fe,
                ke,
                volume,
                J,
                invJ,
                &shape)) {

            continue;
        }

        double div_u = 0.0;

        for (int a = 0; a < 4; ++a) {
            const int gid =
                fe->conn[ke][a];

            div_u +=
                vals->v[gid][0]
                * dN_dx[a][0];

            div_u +=
                vals->v[gid][1]
                * dN_dx[a][1];

            div_u +=
                vals->v[gid][2]
                * dN_dx[a][2];
        }

        const double abs_div_u =
            fabs(div_u);

        if (
            local.valid &&
            abs_div_u <= local.abs_div_u) {

            continue;
        }

        local.valid = 1;
        local.abs_div_u = abs_div_u;
        local.div_u = div_u;
        local.rank =
            monolis_mpi_get_global_my_rank();

        local.elem = ke;
        local.volume = volume;
        local.detJ = detJ;
        local.tet_shape = shape;
        local.characteristic_length =
            cbrt(volume);

        local.x[0] = 0.0;
        local.x[1] = 0.0;
        local.x[2] = 0.0;

        for (int a = 0; a < 4; ++a) {
            const int gid =
                fe->conn[ke][a];

            local.tet_nodes[a] =
                gid;

            local.node_x[a][0] =
                fe->x[gid][0];

            local.node_x[a][1] =
                fe->x[gid][1];

            local.node_x[a][2] =
                fe->x[gid][2];

            local.x[0] +=
                0.25 * local.node_x[a][0];

            local.x[1] +=
                0.25 * local.node_x[a][1];

            local.x[2] +=
                0.25 * local.node_x[a][2];

            local.nodal_velocity[a][0] =
                vals->v[gid][0];

            local.nodal_velocity[a][1] =
                vals->v[gid][1];

            local.nodal_velocity[a][2] =
                vals->v[gid][2];
        }
    }

    if (
        t4loc_reduce_div_location(
            &local,
            monolis_com,
            &global) != 0) {

        local_status = -1;
    }

    if (!t4loc_all_ranks_ok(
            local_status,
            monolis_com,
            "maximum divergence location reduction")) {

        return -1;
    }

    int output_status = 0;

    if (
        monolis_mpi_get_global_my_rank() == 0) {

        output_status =
            t4loc_write_div_location(
                &global,
                t,
                directory);
    }

    if (!t4loc_all_ranks_ok(
            output_status,
            monolis_com,
            "maximum divergence location output")) {

        return -1;
    }

    return 0;
}

/*
 * t4n_nitsche_row_sum_diagnostic.c
 *
 * INSERTION BLOCK
 * ----------------
 * Insert this block in the same C translation unit as the existing static
 * TET4/TRI3 Nitsche helpers. Place it after:
 *
 *   t4n_wall_sym_nitsche_local_vecmat()
 *   t4n_tri3_normal()
 *   t4n_surface_to_volume_nodes()
 *   t4n_grad_phys_metric()
 *   t4n_tri_area_factor()
 *   t4_tet_face_height_from_jacobians()
 *
 * It does not modify the MONOLIS matrix or RHS.
 *
 * Exact identity checked at every surface quadrature point:
 *
 *   sum_i vec_i / dt
 *       = p n
 *         - mu_eff n [n.grad(u)n + div(u)]
 *         + beta ((u-Uw).n) n
 *
 * Since the existing Nitsche force file uses
 *
 *   penalty_file = -beta ((u-Uw).n) n,
 *
 * the same identity is
 *
 *   row_sum = pressure_body + weak_viscous_body - penalty_file.
 */

static void t4n_nitsche_row_sum_diag_zero(
    T4NitscheRowSumDiagnostics* diag)
{
    if (diag != NULL) {
        memset(diag, 0, sizeof(*diag));
    }
}


/* Row diagnostics use the same overlap-aware owner map as assembly. */
static int t4n_row_diag_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int** owner,
    int** lface,
    int* num_unmapped)
{
    return t4n_make_owner_map_overlap(
        surf,
        fe,
        owner,
        lface,
        num_unmapped);
}


int calc_nitsche_row_sum_reaction_tet4_tri3_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw_in[3],
    double gamma_n,
    T4NitscheRowSumDiagnostics* diag)
{
    t4n_nitsche_row_sum_diag_zero(diag);

    if (
        surf == NULL ||
        fe == NULL ||
        basis_surf == NULL ||
        vals == NULL ||
        diag == NULL
    ) {
        fprintf(stderr, "[nitsche-row-sum] NULL input\n");
        return -1;
    }

    if (
        fe->local_num_nodes != 4 ||
        (
            surf->total_num_elems > 0 &&
            surf->local_num_nodes != 3
        )
    ) {
        fprintf(
            stderr,
            "[nitsche-row-sum] TET4/TRI3 expected: "
            "fe=%d surf=%d ne_surf=%d\n",
            fe->local_num_nodes,
            surf->local_num_nodes,
            surf->total_num_elems);

        return -1;
    }

    if (
        !isfinite(vals->dt) ||
        vals->dt <= 0.0
    ) {
        fprintf(
            stderr,
            "[nitsche-row-sum] invalid dt=%.17e\n",
            vals->dt);

        return -1;
    }

    if (surf->total_num_elems == 0) {
        return 0;
    }

    /* Row-sum reaction is a Nitsche diagnostic.  Strong/projected-strong
     * modes have no Nitsche boundary row to report. */
    if (!t4_wall_mode_uses_nitsche_normal()) {
        return 0;
    }

    const double Uw[3] = {
        Uw_in != NULL ? Uw_in[0] : 0.0,
        Uw_in != NULL ? Uw_in[1] : 0.0,
        Uw_in != NULL ? Uw_in[2] : 0.0
    };

    int* owner = NULL;
    int* lface = NULL;
    int num_unmapped = 0;

    if (
        t4n_row_diag_make_owner_map(
            surf,
            fe,
            &owner,
            &lface,
            &num_unmapped) != 0
    ) {
        fprintf(
            stderr,
            "[nitsche-row-sum] owner-map construction failed\n");

        free(owner);
        free(lface);
        return -1;
    }

    const int ns = 3;
    const int nv = 4;
    const int np = basis_surf->num_integ_points;

    const double area_factor =
        t4n_tri_area_factor(
            basis_surf);

    for (int es = 0;
         es < surf->total_num_elems;
         ++es) {

        const int ke = owner[es];

        if (ke < 0) {
            diag->skipped_surface_faces += 1.0;
            continue;
        }

        if (ke >= fe->total_num_elems) {
            diag->invalid_faces += 1.0;
            continue;
        }

        diag->mapped_surface_faces += 1.0;

        double xsurf[3][3];
        double xvol[4][3];
        int face_is_valid = 1;

        for (int a = 0; a < ns; ++a) {
            const int gid = surf->conn[es][a];

            if (gid < 0 || gid >= fe->total_num_nodes) {
                diag->invalid_faces += 1.0;
                face_is_valid = 0;
                break;
            }

            xsurf[a][0] = fe->x[gid][0];
            xsurf[a][1] = fe->x[gid][1];
            xsurf[a][2] = fe->x[gid][2];
        }

        if (!face_is_valid) {
            continue;
        }

        for (int a = 0; a < nv; ++a) {
            const int gid = fe->conn[ke][a];

            if (gid < 0 || gid >= fe->total_num_nodes) {
                diag->invalid_faces += 1.0;
                face_is_valid = 0;
                break;
            }

            xvol[a][0] = fe->x[gid][0];
            xvol[a][1] = fe->x[gid][1];
            xvol[a][2] = fe->x[gid][2];
        }

        if (!face_is_valid) {
            continue;
        }

        int s2v[3];

        if (
            !t4n_surface_to_volume_nodes(
                surf,
                fe,
                es,
                ke,
                s2v)
        ) {
            diag->invalid_faces += 1.0;
            continue;
        }

        const int opposite_local_node =
            6 - s2v[0] - s2v[1] - s2v[2];
        if (opposite_local_node < 0 || opposite_local_node >= nv) {
            diag->invalid_faces += 1.0;
            continue;
        }
        const int opposite_gid = fe->conn[ke][opposite_local_node];

        double dN_dx[4][3];
        double J_inv_face[3][3];
        double Jacobian_face = 0.0;

        if (
            !t4n_grad_phys_metric(
                xvol,
                dN_dx,
                J_inv_face,
                &Jacobian_face)
        ) {
            diag->invalid_faces += 1.0;
            continue;
        }

        double h_e_vms =
            cbrt(fabs(Jacobian_face) / 6.0);

        if (h_e_vms <= 1.0e-12) {
            h_e_vms = 1.0e-12;
        }

        double xc[3] = {0.0, 0.0, 0.0};
        double xf[3] = {0.0, 0.0, 0.0};

        for (int a = 0; a < nv; ++a) {
            xc[0] += xvol[a][0];
            xc[1] += xvol[a][1];
            xc[2] += xvol[a][2];
        }

        xc[0] /= 4.0;
        xc[1] /= 4.0;
        xc[2] /= 4.0;

        for (int a = 0; a < ns; ++a) {
            xf[0] += xsurf[a][0];
            xf[1] += xsurf[a][1];
            xf[2] += xsurf[a][2];
        }

        xf[0] /= 3.0;
        xf[1] /= 3.0;
        xf[2] /= 3.0;

        const double svec[3] = {
            xf[0] - xc[0],
            xf[1] - xc[1],
            xf[2] - xc[2]
        };

        for (int p = 0; p < np; ++p) {
            double nrm[3];

            t4n_tri3_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nrm);

            const double Jface_raw =
                t4n_norm3(nrm);

            if (Jface_raw <= 1.0e-300) {
                diag->invalid_faces += 1.0;
                continue;
            }

            double n[3] = {
                nrm[0] / Jface_raw,
                nrm[1] / Jface_raw,
                nrm[2] / Jface_raw
            };

            if (t4n_dot3(n, svec) < 0.0) {
                n[0] *= -1.0;
                n[1] *= -1.0;
                n[2] *= -1.0;
            }

            /*
             * Restriction of the TET4 basis to the TRI3 face.
             */
            double Nv[4] = {0.0, 0.0, 0.0, 0.0};

            for (int a = 0; a < ns; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double sum_N = 0.0;
            double sum_grad_N[3] = {0.0, 0.0, 0.0};

            for (int i = 0; i < nv; ++i) {
                sum_N += Nv[i];

                for (int d = 0; d < 3; ++d) {
                    sum_grad_N[d] += dN_dx[i][d];
                }
            }

            diag->max_partition_unity_abs_error =
                fmax(
                    diag->max_partition_unity_abs_error,
                    fabs(sum_N - 1.0));

            for (int d = 0; d < 3; ++d) {
                diag->max_gradient_sum_abs_error =
                    fmax(
                        diag->max_gradient_sum_abs_error,
                        fabs(sum_grad_N[d]));
            }

            double u_q[3] = {0.0, 0.0, 0.0};
            double u_old_q[3] = {0.0, 0.0, 0.0};
            double p_q = 0.0;
            double grad_u_q[3][3] = {{0.0}};
            double grad_p_q[3] = {0.0, 0.0, 0.0};

            for (int i = 0; i < nv; ++i) {
                const int gid = fe->conn[ke][i];

                p_q += Nv[i] * vals->p[gid];

                for (int a = 0; a < 3; ++a) {
                    u_q[a] += Nv[i] * vals->v[gid][a];

                    if (vals->v_old != NULL) {
                        u_old_q[a] += Nv[i] * vals->v_old[gid][a];
                    } else {
                        u_old_q[a] += Nv[i] * vals->v[gid][a];
                    }

                    for (int d = 0; d < 3; ++d) {
                        grad_u_q[a][d] +=
                            dN_dx[i][d] * vals->v[gid][a];
                    }
                }

                for (int d = 0; d < 3; ++d) {
                    grad_p_q[d] += dN_dx[i][d] * vals->p[gid];
                }
            }

            double* grad_u_ptr[3] = {
                grad_u_q[0],
                grad_u_q[1],
                grad_u_q[2]
            };

            double mu_eff_face = mu_molecular;
            double tau_face = 0.0;
            double tau_c_face = 0.0;

            BBFE_vms_mu_eff_tau(
                &mu_eff_face,
                &tau_face,
                &tau_c_face,
                J_inv_face,
                Jacobian_face,
                h_e_vms,
                u_q,
                u_old_q,
                grad_u_ptr,
                grad_p_q,
                rho,
                mu_molecular,
                vals->dt,
                vals->C_vms,
                vals->vms_cap_coeff);

            if (
                !isfinite(mu_eff_face) ||
                mu_eff_face <= 0.0
            ) {
                mu_eff_face = mu_molecular;
            }

            /*
             * This must be identical to the current Nitsche assembly.
             * Case C uses the same |detJ|/Jface_raw = 3V/A normal scale.
             */
            const double h_n =
                t4_tet_face_height_from_jacobians(
                    Jacobian_face,
                    Jface_raw);

            const double y_exchange = h_n;

            double u_exchange[3];
            if (g_t4_wall_exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE) {
                u_exchange[0] = vals->v[opposite_gid][0];
                u_exchange[1] = vals->v[opposite_gid][1];
                u_exchange[2] = vals->v[opposite_gid][2];
            } else {
                u_exchange[0] = u_q[0];
                u_exchange[1] = u_q[1];
                u_exchange[2] = u_q[2];
            }

            const double w =
                basis_surf->integ_weight[p]
                * area_factor
                * Jface_raw;

            /* Direct sum of the four local momentum residual rows. */
            double row_q[3] = {0.0, 0.0, 0.0};

            for (int i = 0; i < nv; ++i) {
                double vec_i[4];
                double dummy_mat[4][4];
                const double zero_grad[3] = {0.0, 0.0, 0.0};

                t4n_wall_local_vecmat(
                    vec_i,
                    dummy_mat,
                    Nv[i],
                    0.0,
                    0.0,
                    dN_dx[i],
                    zero_grad,
                    n,
                    u_q,
                    u_exchange,
                    p_q,
                    grad_u_q,
                    Uw,
                    rho,
                    mu_eff_face,
                    mu_molecular,
                    h_n,
                    y_exchange,
                    vals->dt,
                    gamma_n);

                for (int a = 0; a < 3; ++a) {
                    row_q[a] += vec_i[a] / vals->dt;
                }
            }

            /*
             * Closed-form decomposition using exactly the same mu_eff,
             * normal, h_n, beta and quadrature state.
             */
            const double div_u =
                grad_u_q[0][0]
                + grad_u_q[1][1]
                + grad_u_q[2][2];

            double grad_u_n[3] = {0.0, 0.0, 0.0};

            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    grad_u_n[a] += grad_u_q[a][b] * n[b];
                }
            }

            double r[3];

            for (int a = 0; a < 3; ++a) {
                r[a] = u_q[a] - Uw[a];
            }

            const double h_eff = fmax(h_n, 1.0e-12);
            const double dt_eff = fmax(vals->dt, 1.0e-30);

            double beta = 0.0;

            t4n_get_penalty_coefficients(
                rho,
                mu_eff_face,
                h_eff,
                dt_eff,
                gamma_n,
                &beta,
                NULL,
                NULL);

            const double rn = t4n_dot3(r, n);
            const double grad_u_n_dot_n = t4n_dot3(grad_u_n, n);

            double ut_row[3] = {0.0, 0.0, 0.0};
            double beta_wall_row = 0.0;
            if (g_t4_wall_mode == T4_WALL_NITSCHE_SPALDING) {
                t4n_spalding_wall_coefficients(
                    u_exchange,
                    Uw,
                    n,
                    rho,
                    mu_molecular,
                    y_exchange,
                    ut_row,
                    NULL,
                    NULL,
                    &beta_wall_row);
            }

            double pressure_body_q[3];
            double weak_viscous_body_q[3];
            double penalty_file_q[3];
            double reference_q[3];

            for (int a = 0; a < 3; ++a) {
                pressure_body_q[a] = p_q * n[a];

                if (g_t4_wall_mode == T4_WALL_NITSCHE_NOSLIP) {
                    /* full-vector Nitsche row sum */
                    weak_viscous_body_q[a] =
                        -mu_eff_face * (grad_u_n[a] + n[a] * div_u);
                    penalty_file_q[a] = -beta * r[a];
                } else {
                    /* normal Nitsche, with optional modeled wall-law shear */
                    weak_viscous_body_q[a] =
                        -mu_eff_face * n[a] * (grad_u_n_dot_n + div_u)
                        + beta_wall_row * ut_row[a];
                    penalty_file_q[a] = -beta * rn * n[a];
                }

                reference_q[a] =
                    pressure_body_q[a]
                    + weak_viscous_body_q[a]
                    - penalty_file_q[a];

                diag->max_qp_row_reference_abs_error =
                    fmax(
                        diag->max_qp_row_reference_abs_error,
                        fabs(row_q[a] - reference_q[a]));
            }

            for (int a = 0; a < 3; ++a) {
                diag->force_row_sum[a] += w * row_q[a];
                diag->force_reference[a] += w * reference_q[a];
                diag->force_pressure[a] += w * pressure_body_q[a];

                diag->force_weak_viscous_eff[a] +=
                    w * weak_viscous_body_q[a];

                diag->force_penalty_file[a] +=
                    w * penalty_file_q[a];
            }

            diag->quadrature_points += 1.0;
        }

    }

    if (
        fabs(
            diag->skipped_surface_faces
            - (double)num_unmapped) > 0.5
    ) {
        fprintf(
            stderr,
            "[nitsche-row-sum] unmapped count mismatch: "
            "map=%d loop=%.0f\n",
            num_unmapped,
            diag->skipped_surface_faces);

        free(owner);
        free(lface);
        return -1;
    }

    free(owner);
    free(lface);

    return diag->invalid_faces > 0.5 ? -1 : 0;
}


/*
 * All ranks must call this in the same order.
 */
void nitsche_row_sum_diagnostics_allreduce(
    T4NitscheRowSumDiagnostics* diag,
    MONOLIS_COM* monolis_com)
{
    if (diag == NULL || monolis_com == NULL) {
        return;
    }

    monolis_allreduce_R(
        3,
        diag->force_row_sum,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        3,
        diag->force_reference,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        3,
        diag->force_pressure,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        3,
        diag->force_weak_viscous_eff,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        3,
        diag->force_penalty_file,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->mapped_surface_faces,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->skipped_surface_faces,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->quadrature_points,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->invalid_faces,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->max_qp_row_reference_abs_error,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->max_partition_unity_abs_error,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    monolis_allreduce_R(
        1,
        &diag->max_gradient_sum_abs_error,
        MONOLIS_MPI_MAX,
        monolis_com->comm);
}



void solver_fom_VMS(
    FE_SYSTEM        sys,
    const BBFE_DATA* surf_wall,
    const BBFE_BASIS* basis_wall,
    double           t,
    const int        step)
{
    const int ndof = sys.fe.total_num_nodes * 4;

    const int myrank    = monolis_mpi_get_global_my_rank();
    const int comm_size = monolis_mpi_get_global_comm_size();


    /*
      main 関数側の収束判定に合わせる。
      serial 実行では mono_com.n_internal_vertex が未設定の場合があるため、
      1 MPI のときは total_num_nodes を内部節点数として扱う。
    */
    int num_internal_vtx;
    if (comm_size == 1) {
        num_internal_vtx = sys.fe.total_num_nodes;
    }
    else {
        num_internal_vtx = sys.mono_com.n_internal_vertex;
    }


    printf("now\n");

    double tt = monolis_get_time_global_sync();

    if (myrank == 0) {
        printf("\n%s ----------------- Time step %d ----------------\n",
               CODENAME,
               step);
    }

    /*
      it = 0 の Newton 残差を保存し、
      更新後残差との比 |r| / |r_old| で収束判定する。
    */
    double* rvec_old = BB_std_calloc_1d_double(rvec_old, ndof);

    const double tiny        = 1.0e-30;
    const double eps_NR      = 1.0e-6;
    const int    max_iter_NR = 20;

    int converged = 0;

    /*
      固定壁なら Uw = 0。
      移動壁ならここを変更する。
    */
    const double Uw[3] = {0.0, 0.0, 0.0};

    /*
      Nitsche 法線ペナルティ係数。
      まずは 10〜100 程度から調整。
    */
    const double gamma_n = g_t4_nitsche_gamma_n;

    if (t4_wall_mode_uses_projected_strong_normal()) {
        if (t4_wall_projector_prepare(surf_wall, &(sys.fe)) != 0) {
            fprintf(stderr, "%s ERROR: failed to build projected-strong wall projector.\n", CODENAME);
            exit(EXIT_FAILURE);
        }

        /* Put hot-start/current and previous-time states on the admissible
         * wall manifold before any volume residual is evaluated. */
        t4_wall_project_velocity_state(sys.vals.v, Uw);
        t4_wall_project_velocity_state(sys.vals.v_old, Uw);
    }

    for (int it = 0; it < max_iter_NR; ++it) {
        if (myrank == 0) {
            printf("\n%s ----------------- Time step %d : NR step %d ----------------\n",
                   CODENAME,
                   step,
                   it);
        }

        /*
          Newton-Raphson 用の接線行列・残差ベクトルを作成。
          B = -F として扱う。
        */
        monolis_clear_mat_value_R(&(sys.monolis));

        for (int i = 0; i < ndof; ++i) {
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        /*
          体積要素：VMS 版 Newton 行列・残差
        */
        set_element_mat_NR_linear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        set_element_vec_NR_linear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        set_element_mat_NR_nonlinear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        set_element_vec_NR_nonlinear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        /*
          Runtime-selectable Ahmed wall terms.
          Depending on g_t4_wall_mode this assembles full-vector Nitsche,
          normal-only Nitsche, Spalding tangential traction, or no surface
          term for strong/projected-strong modes.  Add before Dirichlet BC.
        */

        set_wall_face_vecmat_symmetric_nitsche_hex_or_tet(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            surf_wall,
            basis_wall,
            &(sys.vals),
            sys.vals.density,
            sys.vals.viscosity,
            Uw,
            gamma_n);

        if (t4_wall_mode_uses_projected_strong_normal()) {
            t4_wall_project_velocity_rhs(sys.monolis.mat.R.B);
        }

        /*
          強制 Dirichlet 条件
        */
        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.monolis),
            sys.fe.total_num_nodes,
            4,
            &(sys.bc_NR),
            sys.monolis.mat.R.B);

        /*
          初回 Newton 反復時の残差を保存。
          main 関数側の store_initial_NR_residual と同じ役割。
        */
        if (it == 0) {
            for (int i = 0; i < ndof; ++i) {
                rvec_old[i] = sys.monolis.mat.R.B[i];
            }
        }

        /*
          線形ソルバで Newton 修正量を解く。
        */
        ROM_monowrap_solve(
            &(sys.monolis),
            &(sys.mono_com),
            sys.monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        BBFE_fluid_sups_renew_velocity(
            sys.vals.delta_v,
            sys.monolis.mat.R.X,
            sys.fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys.vals.delta_p,
            sys.monolis.mat.R.X,
            sys.fe.total_num_nodes);

        if (t4_wall_mode_uses_projected_strong_normal()) {
            /* Restrict every Newton correction to the local tangent/edge
             * subspace.  This is the strong-normal enforcement used by the
             * TET4-only projected modes. */
            t4_wall_project_velocity_increment(sys.vals.delta_v);
        }

        update_velocity_pressure_NR(
            sys.vals.v,
            sys.vals.delta_v,
            sys.vals.p,
            sys.vals.delta_p,
            sys.fe.total_num_nodes);

        if (t4_wall_mode_uses_projected_strong_normal()) {
            /* Re-project to remove round-off and keep the iterate exactly on
             * the nodal admissible subspace before residual re-evaluation. */
            t4_wall_project_velocity_state(sys.vals.v, Uw);
        }

        /*
          更新後の残差を再評価する。
          収束判定用なので、接線行列は不要。
          ただし set_wall_face_vecmat_symmetric_nitsche は行列も足す設計なので、
          clear 後に呼び出して、ここでは B だけを評価に使う。
        */
        monolis_clear_mat_value_R(&(sys.monolis));

        for (int i = 0; i < ndof; ++i) {
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        set_element_vec_NR_linear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        set_element_vec_NR_nonlinear_VMS(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));

        set_wall_face_vecmat_symmetric_nitsche_hex_or_tet(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            surf_wall,
            basis_wall,
            &(sys.vals),
            sys.vals.density,
            sys.vals.viscosity,
            Uw,
            gamma_n);

        if (t4_wall_mode_uses_projected_strong_normal()) {
            t4_wall_project_velocity_rhs(sys.monolis.mat.R.B);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.monolis),
            sys.fe.total_num_nodes,
            4,
            &(sys.bc_NR),
            sys.monolis.mat.R.B);

        /*
          main 関数側の evaluate_NR_residual と同じ考え方で、
          内部自由度だけの二乗和を取り、MPI では総和する。
        */
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            num_internal_vtx * 4,
            1);

        double norm_r = calc_internal_norm_1d(
            sys.monolis.mat.R.B,
            num_internal_vtx * 4,
            1);

        if (comm_size != 1) {
            monolis_allreduce_R(
                1,
                &norm_r_old,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);

            monolis_allreduce_R(
                1,
                &norm_r,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);
        }

        const double nrm_r_old = sqrt(norm_r_old);
        const double nrm_r     = sqrt(norm_r);
        const double rel_r     = nrm_r / fmax(nrm_r_old, tiny);

        if (myrank == 0) {
            printf(" residual total    : nrm = %e, old = %e, rel = %e\n",
                   nrm_r,
                   nrm_r_old,
                   rel_r);
        }

        /*
          main 関数側と同じく、初回残差に対する相対残差で収束判定する。
        */
        if (rel_r < eps_NR) {
            double max_du = 0.0;

            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    const double du =
                        fabs(sys.vals.v[ii][d] - sys.vals.v_old[ii][d]);

                    if (du > max_du) {
                        max_du = du;
                    }
                }
            }

            if (myrank == 0) {
                printf("[step %d] NR converged at it = %d: |r|/|r_old| = %.6e\n",
                       step,
                       it,
                       rel_r);

                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n",
                       step,
                       max_du);
            }

            ROM_BB_vec_copy_2d(
                sys.vals.v,
                sys.vals.v_old,
                sys.fe.total_num_nodes,
                3);

            converged = 1;
            break;
        }
    }

    if (!converged) {
        if (myrank == 0) {
            printf("[step %d] WARNING: NR did not converge within %d iterations.\n",
                   step,
                   max_iter_NR);
        }
    }

    BB_std_free_1d_double(rvec_old, ndof);
}


