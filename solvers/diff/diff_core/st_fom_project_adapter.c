#include "st_fom_project_adapter.h"

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef STSB_DEGREE
#define STSB_DEGREE 1
#endif

#ifndef STSB_NUM_WINDOW_SLABS
#define STSB_NUM_WINDOW_SLABS 4
#endif

/*
 * The test_thermal convection-diffusion driver is expected to provide these
 * application helpers, analogous to BBFE_fluid_pre/finalize.  If your local
 * branch uses different names, override them from Makefile_ST:
 *
 *   ST_CONVDIFF_PRE=my_pre ST_CONVDIFF_FINALIZE=my_finalize
 */
#ifndef ST_FOM_CONVDIFF_PRE
#define ST_FOM_CONVDIFF_PRE BBFE_convdiff_pre
#endif

#ifndef ST_FOM_CONVDIFF_FINALIZE
#define ST_FOM_CONVDIFF_FINALIZE BBFE_convdiff_finalize
#endif

/*
 * BBFE_convdiff_pre() final boolean selects the partitioned-mesh path in the
 * existing diffusion driver.  The supplied execution script always creates
 * parted.0 and launches the solver from that case, so keep it enabled by
 * default.  It can be overridden from Makefile_ST if needed.
 */
#ifndef ST_FOM_CONVDIFF_PARTITIONED
#define ST_FOM_CONVDIFF_PARTITIONED 1
#endif

#define ST_FOM_BUFFER_SIZE 10000

static int st_env_int(const char* name, int default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    long value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }

    errno = 0;
    value = strtol(text, &end, 10);
    if(errno != 0 || end == text || *end != '\0' ||
       value < 0 || value > 2147483647L)
    {
        fprintf(stderr, "ERROR: invalid integer %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static const char* st_case_directory(int argc, char* argv[])
{
    if(argc >= 2 && argv[1] != NULL && argv[1][0] != '-'){
        return argv[1];
    }
    return "./";
}

static void st_read_conditions(ST_FOM_PROJECT* project)
{
    FE_SYSTEM* fom = &project->fom;
    char filename[4096];
    project->num_ip_each_axis = 3;

    fom->vals.dt = 0.01;
    fom->vals.mat_epsilon = 1.0e-10;
    fom->vals.mat_max_iter = 10000;
    project->finish_time = 1.0;
    project->output_interval = 1;
    project->time_degree = STSB_DEGREE;
    project->slabs_per_window = STSB_NUM_WINDOW_SLABS;

    snprintf(filename, sizeof(filename), "%s/%s", project->directory, "cond.dat");

    {
        FILE* fp = fopen(filename, "r");
        if(fp != NULL){
            fclose(fp);
            (void)BB_std_read_file_get_val_int_p(
                &project->num_ip_each_axis,
                filename,
                "#num_ip_each_axis",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_double_p(
                &fom->vals.mat_epsilon,
                filename,
                "#mat_epsilon",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_int_p(
                &fom->vals.mat_max_iter,
                filename,
                "#mat_max_iter",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_double_p(
                &fom->vals.dt,
                filename,
                "#time_spacing",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_double_p(
                &project->finish_time,
                filename,
                "#finish_time",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_int_p(
                &project->output_interval,
                filename,
                "#output_interval",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_int_p(
                &project->time_degree,
                filename,
                "#st_time_degree",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
            (void)BB_std_read_file_get_val_int_p(
                &project->slabs_per_window,
                filename,
                "#st_slabs_per_window",
                ST_FOM_BUFFER_SIZE,
                CODENAME);
        }
    }

    project->time_degree = st_env_int(
        "ST_FOM_TIME_DEGREE", project->time_degree);
    project->slabs_per_window = st_env_int(
        "ST_FOM_SLABS_PER_WINDOW", project->slabs_per_window);

    if(project->num_ip_each_axis <= 0 || !isfinite(fom->vals.dt) || fom->vals.dt <= 0.0 ||
       !isfinite(project->finish_time) || project->finish_time <= 0.0 ||
       fom->vals.mat_epsilon <= 0.0 || fom->vals.mat_max_iter <= 0 ||
       project->time_degree < 0 || project->time_degree > 2 ||
       project->slabs_per_window <= 0)
    {
        fprintf(stderr, "ERROR: invalid ST-FOM calculation conditions.\n");
        exit(EXIT_FAILURE);
    }

}

static int st_allocate_values(ST_FOM_PROJECT* project)
{
    FE_SYSTEM* fom = &project->fom;
    const size_t nnode = (size_t)fom->fe.total_num_nodes;

    fom->vals.T = (double*)calloc(nnode, sizeof(double));
    fom->vals.theo_sol = (double*)calloc(nnode, sizeof(double));

    if(fom->vals.T == NULL || fom->vals.theo_sol == NULL){
        return 1;
    }
    project->values_initialized = 1;

    for(int node = 0; node < fom->fe.total_num_nodes; node++){
        const double exact = manusol_get_sol(
            fom->fe.x[node][0],
            fom->fe.x[node][1],
            fom->fe.x[node][2],
            0.0);
        fom->vals.T[node] = exact;
        fom->vals.theo_sol[node] = exact;
    }
    return 0;
}

int ST_fom_project_initialize(
    ST_FOM_PROJECT* project,
    int argc,
    char* argv[])
{
    FE_SYSTEM* fom;
    const char* directory;
    const int rank = monolis_mpi_get_global_my_rank();

    if(project == NULL){
        return 1;
    }
    memset(project, 0, sizeof(*project));

    directory = st_case_directory(argc, argv);
    snprintf(project->directory, sizeof(project->directory), "%s", directory);
    fom = &project->fom;

    st_read_conditions(project);

    /*
     * Use the project's native convection-diffusion preprocessing routine.
     * Its actual signature (from convdiff_core.h) is
     *
     *   BBFE_convdiff_pre(fe, basis, bc, monolis, monolis_com,
     *                     argc, argv, directory, num_ip_each_axis,
     *                     is_partitioned);
     *
     * It owns mesh/basis/BC setup and initializes the ordinary spatial
     * MONOLIS communication object.  Do not allocate/read BCs or initialize
     * the base MONOLIS object a second time here.
     */
    ST_FOM_CONVDIFF_PRE(
        &fom->fe,
        &fom->basis,
        &fom->bc,
        &fom->monolis,
        &fom->monolis_com,
        argc,
        argv,
        project->directory,
        project->num_ip_each_axis,
        (bool)(ST_FOM_CONVDIFF_PARTITIONED != 0));

    project->fe_initialized = 1;
    project->bc_initialized = 1;
    project->base_monolis_initialized = 1;

    /*
     * BBFE_convdiff_pre() reads the mesh/basis/BC and initializes the
     * spatial communication object, but the element geometry used by the
     * FEM integrators must be built explicitly.  The verified fluid ST-FOM
     * lifecycle does these two calls in this order before any element
     * assembly.  Without them fe->geo[e][p].Jacobian/grad_N remain zero,
     * making both the mass and diffusion matrices identically zero.
     */
    BBFE_elemmat_set_Jacobi_mat(&fom->fe, &fom->basis);
    BBFE_elemmat_set_shapefunc_derivative(&fom->fe, &fom->basis);

    /*
     * Fail early if geometry preparation did not populate the Jacobians.
     * Also print representative manufactured coefficients so a bad physics
     * setup can be distinguished from a geometry/setup failure.
     */
    {
        double min_abs_j = DBL_MAX;
        double max_abs_j = 0.0;
        int nonzero_j = 0;

        for(int e = 0; e < fom->fe.total_num_elems; e++){
            for(int ip = 0; ip < fom->basis.num_integ_points; ip++){
                const double aj = fabs(fom->fe.geo[e][ip].Jacobian);
                if(!isfinite(aj)){
                    fprintf(stderr,
                        "ERROR [rank %d]: non-finite element Jacobian at e=%d ip=%d\n",
                        rank, e, ip);
                    ST_fom_project_finalize(project);
                    return 1;
                }
                if(aj > 0.0){
                    nonzero_j++;
                    if(aj < min_abs_j) min_abs_j = aj;
                    if(aj > max_abs_j) max_abs_j = aj;
                }
            }
        }

        if(nonzero_j == 0){
            fprintf(stderr,
                "ERROR [rank %d]: all element Jacobians are zero after "
                "BBFE_elemmat_set_Jacobi_mat()/set_shapefunc_derivative().\n",
                rank);
            ST_fom_project_finalize(project);
            return 1;
        }

        if(fom->fe.total_num_elems > 0 && fom->basis.num_integ_points > 0){
            const int nl = fom->fe.local_num_nodes;
            double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
            double x_ip[3] = {0.0, 0.0, 0.0};
            double mass_coef = 0.0;
            double diff_coef = 0.0;
            int alloc_ok = (local_x != NULL);

            if(alloc_ok){
                for(int i = 0; i < nl; i++){
                    local_x[i] = (double*)calloc(3u, sizeof(double));
                    if(local_x[i] == NULL){
                        alloc_ok = 0;
                        break;
                    }
                }
            }

            if(alloc_ok){
                BBFE_elemmat_set_local_array_vector(
                    local_x, &fom->fe, fom->fe.x, 0, 3);
                BBFE_std_mapping_vector3d(
                    x_ip, nl, local_x, fom->basis.N[0]);
                mass_coef = manusol_get_mass_coef(x_ip);
                diff_coef = manusol_get_diff_coef(x_ip);
            }

            if(local_x != NULL){
                for(int i = 0; i < nl; i++) free(local_x[i]);
                free(local_x);
            }

            printf(
                "ST-FOM rank %d geometry ready: nonzero_J=%d min|J|=%.15e "
                "max|J|=%.15e sample_mass_coef=%.15e sample_diff_coef=%.15e\n",
                rank, nonzero_j, min_abs_j, max_abs_j, mass_coef, diff_coef);
            fflush(stdout);
        }
    }

    if(st_allocate_values(project) != 0){
        ST_fom_project_finalize(project);
        return 1;
    }

    if(rank == 0){
        printf("\n%s ---------- Normal ST-FOM conditions ----------\n", CODENAME);
        printf("%s case directory       : %s\n", CODENAME, project->directory);
        printf("%s integration pts/axis: %d\n", CODENAME, project->num_ip_each_axis);
        printf("%s time spacing         : %.15e\n", CODENAME, fom->vals.dt);
        printf("%s finish time          : %.15e\n", CODENAME, project->finish_time);
        printf("%s matrix epsilon       : %.15e\n", CODENAME, fom->vals.mat_epsilon);
        printf("%s matrix max iter      : %d\n", CODENAME, fom->vals.mat_max_iter);
        printf("%s ST time degree       : %d\n", CODENAME, project->time_degree);
        printf("%s ST slabs/window      : %d\n", CODENAME, project->slabs_per_window);
        printf("%s ROM/POD              : disabled\n", CODENAME);
        printf("%s -----------------------------------------------\n\n", CODENAME);
    }

    return 0;
}

void ST_fom_project_finalize(ST_FOM_PROJECT* project)
{
    if(project == NULL){
        return;
    }

    if(project->values_initialized){
        free(project->fom.vals.T);
        free(project->fom.vals.theo_sol);
    }

    /*
     * BBFE_convdiff_finalize() releases FE, basis and BC storage allocated by
     * BBFE_convdiff_pre().  Do not free the Dirichlet BC separately.
     */
    if(project->fe_initialized){
        ST_FOM_CONVDIFF_FINALIZE(
            &project->fom.fe,
            &project->fom.basis,
            &project->fom.bc);
    }

    if(project->base_monolis_initialized){
        monolis_finalize(&project->fom.monolis);
    }

    memset(project, 0, sizeof(*project));
}

int ST_fom_project_num_windows(const ST_FOM_PROJECT* project)
{
    double steps_real;
    long long steps;
    long long per_window;

    if(project == NULL || project->fom.vals.dt <= 0.0 ||
       project->finish_time <= 0.0 || project->slabs_per_window <= 0)
    {
        return -1;
    }

    steps_real = project->finish_time / project->fom.vals.dt;
    steps = (long long)llround(steps_real);
    per_window = (long long)project->slabs_per_window;

    if(fabs(steps_real - (double)steps) > 1.0e-9 * fmax(1.0, fabs(steps_real))){
        fprintf(
            stderr,
            "ERROR: finish_time/dt must be an integer. got %.17g\n",
            steps_real);
        return -1;
    }
    if(steps <= 0 || (steps % per_window) != 0){
        fprintf(
            stderr,
            "ERROR: total time steps (%lld) must be divisible by "
            "slabs/window (%lld).\n",
            steps,
            per_window);
        return -1;
    }
    if(steps / per_window > 2147483647LL){
        return -1;
    }

    return (int)(steps / per_window);
}
