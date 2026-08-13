#include "fluid_sups_st_model.h"
#include "fluid_sups_nr_api_compat.h"
#include "fluid_sups_jold.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* FLUID_ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
static const char* FLUID_ID_MAT_EPSILON = "#mat_epsilon";
static const char* FLUID_ID_MAT_MAX_ITER = "#mat_max_iter";
static const char* FLUID_ID_DT = "#time_spacing";
static const char* FLUID_ID_FINISH_TIME = "#finish_time";
static const char* FLUID_ID_OUTPUT_INTERVAL = "#output_interval";
static const char* FLUID_ID_DENSITY = "#density";
static const char* FLUID_ID_VISCOSITY = "#viscosity";
static const char* FLUID_ID_NR_EPSILON = "#nr_epsilon";
static const char* FLUID_ID_NR_MAX_ITER = "#nr_max_iter";
static const char* FLUID_ID_ST_TIME_DEGREE = "#st_time_degree";
static const char* FLUID_ID_ST_SLABS = "#st_slabs_per_window";
static const char* FLUID_ID_WRITE_WINDOW = "#st_write_window_binary";
static const char* FLUID_ID_ST_PREVIOUS_JACOBIAN_FD_EPS =
    "#st_previous_jacobian_fd_epsilon";
static const char* FLUID_ID_ST_DG_JACOBIAN_FD_EPS =
    "#st_dg_jacobian_fd_epsilon";

static const int FLUID_DVAL_NUM_IP_EACH_AXIS = 3;
static const double FLUID_DVAL_MAT_EPSILON = 1.0e-8;
static const int FLUID_DVAL_MAT_MAX_ITER = 10000;
static const double FLUID_DVAL_DT = 0.005;
static const double FLUID_DVAL_FINISH_TIME = 50.0;
static const int FLUID_DVAL_OUTPUT_INTERVAL = 100;
static const double FLUID_DVAL_DENSITY = 1000.0;
static const double FLUID_DVAL_VISCOSITY = 1.0;
static const double FLUID_DVAL_NR_EPSILON = 1.0e-6;
static const int FLUID_DVAL_NR_MAX_ITER = 20;
static const int FLUID_DVAL_ST_TIME_DEGREE = 0;
static const int FLUID_DVAL_ST_SLABS = 4;
static const int FLUID_DVAL_WRITE_WINDOW = 1;
static const double FLUID_DVAL_ST_PREVIOUS_JACOBIAN_FD_EPS = 1.0e-7;
static const double FLUID_DVAL_ST_DG_JACOBIAN_FD_EPS = 1.0e-7;

static const char* FLUID_INPUT_FILENAME_COND = "cond.dat";
static const char* FLUID_INPUT_FILENAME_D_BC_V = "D_bc_v.dat";
static const char* FLUID_OUTPUT_FILENAME_VTK = "st_result_%06d.vtk.%d";
static const int FLUID_BUFFER_SIZE = 10000;

static int fluid_env_int(const char* name, int default_value)
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

static double fluid_env_double(const char* name, double default_value)
{
    const char* text = getenv(name);
    char* end = NULL;
    double value;

    if(text == NULL || text[0] == '\0'){
        return default_value;
    }

    errno = 0;
    value = strtod(text, &end);
    if(errno != 0 || end == text || *end != '\0' || !isfinite(value)){
        fprintf(stderr, "ERROR: invalid floating-point %s=%s\n", name, text);
        exit(EXIT_FAILURE);
    }
    return value;
}

static void fluid_assign_default_values(FLUID_SUPS_VALUES* vals)
{
    memset(vals, 0, sizeof(*vals));
    vals->num_ip_each_axis = FLUID_DVAL_NUM_IP_EACH_AXIS;
    vals->mat_epsilon = FLUID_DVAL_MAT_EPSILON;
    vals->mat_max_iter = FLUID_DVAL_MAT_MAX_ITER;
    vals->dt = FLUID_DVAL_DT;
    vals->finish_time = FLUID_DVAL_FINISH_TIME;
    vals->output_interval = FLUID_DVAL_OUTPUT_INTERVAL;
    vals->density = FLUID_DVAL_DENSITY;
    vals->viscosity = FLUID_DVAL_VISCOSITY;
    vals->nr_epsilon = FLUID_DVAL_NR_EPSILON;
    vals->nr_max_iter = FLUID_DVAL_NR_MAX_ITER;
    vals->st_time_degree = FLUID_DVAL_ST_TIME_DEGREE;
    vals->st_slabs_per_window = FLUID_DVAL_ST_SLABS;
    vals->write_window_binary = FLUID_DVAL_WRITE_WINDOW;
    vals->st_previous_jacobian_fd_epsilon =
        FLUID_DVAL_ST_PREVIOUS_JACOBIAN_FD_EPS;
    vals->st_dg_jacobian_fd_epsilon =
        FLUID_DVAL_ST_DG_JACOBIAN_FD_EPS;
}

static void fluid_read_calc_conditions(
    FLUID_SUPS_VALUES* vals,
    const char* directory)
{
    char filename[FLUID_BUFFER_SIZE];
    FILE* fp;

    fluid_assign_default_values(vals);
    snprintf(filename, sizeof(filename), "%s/%s",
        directory, FLUID_INPUT_FILENAME_COND);

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("%s Calc condition file \"%s\" is not found.\n",
            CODENAME, filename);
        printf("%s Default values are used in this calculation.\n", CODENAME);
    }
    else{
        int num;
        printf("%s Reading condition file \"%s\".\n", CODENAME, filename);
        num = BB_std_read_file_get_val_int_p(
            &vals->num_ip_each_axis, filename,
            FLUID_ID_NUM_IP_EACH_AXIS, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->mat_epsilon, filename,
            FLUID_ID_MAT_EPSILON, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->mat_max_iter, filename,
            FLUID_ID_MAT_MAX_ITER, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->dt, filename,
            FLUID_ID_DT, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->finish_time, filename,
            FLUID_ID_FINISH_TIME, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->output_interval, filename,
            FLUID_ID_OUTPUT_INTERVAL, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->density, filename,
            FLUID_ID_DENSITY, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->viscosity, filename,
            FLUID_ID_VISCOSITY, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->nr_epsilon, filename,
            FLUID_ID_NR_EPSILON, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->nr_max_iter, filename,
            FLUID_ID_NR_MAX_ITER, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->st_time_degree, filename,
            FLUID_ID_ST_TIME_DEGREE, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->st_slabs_per_window, filename,
            FLUID_ID_ST_SLABS, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_int_p(
            &vals->write_window_binary, filename,
            FLUID_ID_WRITE_WINDOW, FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->st_previous_jacobian_fd_epsilon, filename,
            FLUID_ID_ST_PREVIOUS_JACOBIAN_FD_EPS,
            FLUID_BUFFER_SIZE, CODENAME);
        num = BB_std_read_file_get_val_double_p(
            &vals->st_dg_jacobian_fd_epsilon, filename,
            FLUID_ID_ST_DG_JACOBIAN_FD_EPS,
            FLUID_BUFFER_SIZE, CODENAME);
        (void)num;
        fclose(fp);
    }

    /* Environment variables are convenient for staged ST validation runs. */
    vals->st_time_degree = fluid_env_int(
        "FLUID_ST_TIME_DEGREE", vals->st_time_degree);
    vals->st_slabs_per_window = fluid_env_int(
        "FLUID_ST_SLABS_PER_WINDOW", vals->st_slabs_per_window);
    vals->write_window_binary = fluid_env_int(
        "FLUID_ST_WRITE_WINDOW_BINARY", vals->write_window_binary);
    vals->st_previous_jacobian_fd_epsilon = fluid_env_double(
        "FLUID_ST_PREVIOUS_JACOBIAN_FD_EPS",
        vals->st_previous_jacobian_fd_epsilon);
    vals->st_dg_jacobian_fd_epsilon = fluid_env_double(
        "FLUID_ST_DG_JACOBIAN_FD_EPS",
        vals->st_dg_jacobian_fd_epsilon);

    if(vals->num_ip_each_axis <= 0 || vals->mat_epsilon <= 0.0 ||
       vals->mat_max_iter <= 0 || vals->dt <= 0.0 ||
       vals->finish_time <= 0.0 || vals->output_interval <= 0 ||
       vals->density <= 0.0 || vals->viscosity <= 0.0 ||
       vals->nr_epsilon <= 0.0 || vals->nr_max_iter <= 0 ||
       vals->st_time_degree < 0 || vals->st_time_degree > 2 ||
       vals->st_slabs_per_window <= 0 ||
       vals->st_previous_jacobian_fd_epsilon <= 0.0 ||
       vals->st_dg_jacobian_fd_epsilon <= 0.0)
    {
        fprintf(stderr, "ERROR: invalid fluid/ST calculation conditions.\n");
        exit(EXIT_FAILURE);
    }
}

void fluid_sups_st_print_conditions(const FLUID_SUPS_SYSTEM* system)
{
    const FLUID_SUPS_VALUES* vals = &system->vals;

    printf("\n%s ---------- Fluid ST-FOM conditions ----------\n", CODENAME);
    printf("%s #num_ip_each_axis: %d\n", CODENAME, vals->num_ip_each_axis);
    printf("%s #mat_epsilon: %.15e\n", CODENAME, vals->mat_epsilon);
    printf("%s #mat_max_iter: %d\n", CODENAME, vals->mat_max_iter);
    printf("%s #time_spacing: %.15e\n", CODENAME, vals->dt);
    printf("%s #finish_time: %.15e\n", CODENAME, vals->finish_time);
    printf("%s #output_interval: %d\n", CODENAME, vals->output_interval);
    printf("%s #density: %.15e\n", CODENAME, vals->density);
    printf("%s #viscosity: %.15e\n", CODENAME, vals->viscosity);
    printf("%s #nr_epsilon: %.15e\n", CODENAME, vals->nr_epsilon);
    printf("%s #nr_max_iter: %d\n", CODENAME, vals->nr_max_iter);
    printf("%s #st_time_degree: %d\n",
        CODENAME, vals->st_time_degree);
    printf("%s #st_slabs_per_window: %d\n",
        CODENAME, vals->st_slabs_per_window);
    printf("%s #st_previous_jacobian_fd_epsilon: %.15e\n",
        CODENAME, vals->st_previous_jacobian_fd_epsilon);
    printf("%s #st_dg_jacobian_fd_epsilon: %.15e\n",
        CODENAME, vals->st_dg_jacobian_fd_epsilon);
    printf("%s J_old backend is selected by FLUID_ST_JOLD_MODE "
        "(verify|analytic|fd; default=verify).\n", CODENAME);
    printf("%s -----------------------------------------------\n\n", CODENAME);
}

static int fluid_allocate_values(
    FLUID_SUPS_VALUES* vals,
    int total_num_nodes)
{
    vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
    vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);
    vals->delta_v = BB_std_calloc_2d_double(
        vals->delta_v, total_num_nodes, 3);
    vals->delta_p = BB_std_calloc_1d_double(
        vals->delta_p, total_num_nodes);
    vals->v_old = BB_std_calloc_2d_double(
        vals->v_old, total_num_nodes, 3);
    vals->p_old = BB_std_calloc_1d_double(
        vals->p_old, total_num_nodes);

    return (vals->v == NULL || vals->p == NULL || vals->delta_v == NULL ||
            vals->delta_p == NULL || vals->v_old == NULL || vals->p_old == NULL)
        ? 1 : 0;
}

static void fluid_free_values(
    FLUID_SUPS_VALUES* vals,
    int total_num_nodes)
{
    if(vals->v != NULL){
        BB_std_free_2d_double(vals->v, total_num_nodes, 3);
    }
    if(vals->p != NULL){
        BB_std_free_1d_double(vals->p, total_num_nodes);
    }
    if(vals->delta_v != NULL){
        BB_std_free_2d_double(vals->delta_v, total_num_nodes, 3);
    }
    if(vals->delta_p != NULL){
        BB_std_free_1d_double(vals->delta_p, total_num_nodes);
    }
    if(vals->v_old != NULL){
        BB_std_free_2d_double(vals->v_old, total_num_nodes, 3);
    }
    if(vals->p_old != NULL){
        BB_std_free_1d_double(vals->p_old, total_num_nodes);
    }
    vals->v = NULL;
    vals->p = NULL;
    vals->delta_v = NULL;
    vals->delta_p = NULL;
    vals->v_old = NULL;
    vals->p_old = NULL;
}

static void fluid_read_dirichlet_bc(
    BBFE_BC* bc,
    const char* filename,
    const char* directory,
    int total_num_nodes,
    int homogeneous_values)
{
    FILE* fp;
    int n;
    int i;
    int tmp;

    bc->total_num_nodes = total_num_nodes;
    bc->block_size = FLUID_SUPS_PHYSICAL_DOF;
    BBFE_sys_memory_allocation_Dirichlet_bc(
        bc, total_num_nodes, bc->block_size);

    n = total_num_nodes * bc->block_size;
    for(i = 0; i < n; i++){
        bc->D_bc_exists[i] = false;
        bc->imposed_D_val[i] = 0.0;
    }

    fp = NULL;
    fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
    if(fp == NULL){
        printf("%s WARNING: Dirichlet B.C. file \"%s\" is not found.\n",
            CODENAME, filename);
        return;
    }

    BB_std_scan_line(&fp, FLUID_BUFFER_SIZE,
        "%d %d", &bc->num_D_bcs, &tmp);
    printf("%s Num. Dirichlet B.C.: %d, block size: %d\n",
        CODENAME, bc->num_D_bcs, tmp);

    for(i = 0; i < bc->num_D_bcs; i++){
        int node_id;
        int block_id;
        double value;
        int index;

        BB_std_scan_line(&fp, FLUID_BUFFER_SIZE,
            "%d %d %lf", &node_id, &block_id, &value);
        index = bc->block_size * node_id + block_id;
        if(index < 0 || index >= n){
            fprintf(stderr,
                "ERROR: invalid fluid Dirichlet entry node=%d block=%d.\n",
                node_id, block_id);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        bc->D_bc_exists[index] = true;
        bc->imposed_D_val[index] = homogeneous_values ? 0.0 : value;
    }

    fclose(fp);
}

static void fluid_pack_velocity_pressure(
    const FLUID_SUPS_VALUES* vals,
    int total_num_nodes,
    double* vector)
{
    int node;
    for(node = 0; node < total_num_nodes; node++){
        vector[4 * node + 0] = vals->v[node][0];
        vector[4 * node + 1] = vals->v[node][1];
        vector[4 * node + 2] = vals->v[node][2];
        vector[4 * node + 3] = vals->p[node];
    }
}

static void fluid_unpack_velocity_pressure(
    FLUID_SUPS_VALUES* vals,
    int total_num_nodes,
    const double* vector)
{
    int node;
    for(node = 0; node < total_num_nodes; node++){
        vals->v[node][0] = vector[4 * node + 0];
        vals->v[node][1] = vector[4 * node + 1];
        vals->v[node][2] = vector[4 * node + 2];
        vals->p[node] = vector[4 * node + 3];
    }
}

static void fluid_copy_current_to_previous(FLUID_SUPS_SYSTEM* system)
{
    BBFE_fluid_vec_copy_2d(
        system->vals.v,
        system->vals.v_old,
        system->fe.total_num_nodes,
        3);
    BBFE_fluid_vec_copy_1d(
        system->vals.p_old,
        system->vals.p,
        system->fe.total_num_nodes);
}

static void fluid_unpack_previous_velocity_pressure(
    FLUID_SUPS_VALUES* vals,
    int total_num_nodes,
    const double* vector)
{
    int node;

    for(node = 0; node < total_num_nodes; node++){
        vals->v_old[node][0] = vector[4 * node + 0];
        vals->v_old[node][1] = vector[4 * node + 1];
        vals->v_old[node][2] = vector[4 * node + 2];
        vals->p_old[node] = vector[4 * node + 3];
    }
}

static void fluid_load_state_pair(
    FLUID_SUPS_SYSTEM* system,
    const double* current_state,
    const double* previous_state)
{
    fluid_unpack_velocity_pressure(
        &system->vals,
        system->fe.total_num_nodes,
        current_state);
    fluid_unpack_previous_velocity_pressure(
        &system->vals,
        system->fe.total_num_nodes,
        previous_state);
}

static int fluid_unpack_window_slab_state(
    const FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    int slab,
    double* state)
{
    return ST_fom_unpack_dg0_slab(
        layout,
        slab,
        window_state,
        state) == ST_FOM_SUCCESS ? 0 : 1;
}

static void fluid_apply_physical_bc_to_window(
    const FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* window_state)
{
    int node;
    int slab;
    int alpha;
    int component;

    for(node = 0; node < system->fe.total_num_nodes; node++){
        for(component = 0;
            component < FLUID_SUPS_PHYSICAL_DOF;
            component++)
        {
            const int physical_index =
                FLUID_SUPS_PHYSICAL_DOF * node + component;

            if(!system->bc.D_bc_exists[physical_index]){
                continue;
            }

            /*
             * The current cavity data are time independent.  Therefore every
             * temporal interpolation coefficient on a Dirichlet node receives
             * the same physical boundary value.  A later time-dependent BC
             * adapter should evaluate the value at dof_tau[alpha].
             */
            for(slab = 0; slab < layout->num_slabs; slab++){
                for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                    const size_t index = ST_fom_vector_index(
                        layout, node, slab, alpha, component);
                    window_state[index] =
                        system->bc.imposed_D_val[physical_index];
                }
            }
        }
    }
}


static int fluid_window_pressure_gauge_enabled(void)
{
    /*
     * A monolithic incompressible window has one spatially constant pressure
     * null mode for each (slab,time-basis) pressure coefficient unless a
     * gauge is selected.  The ordinary one-step solver can often tolerate the
     * single compatible null mode, but the enlarged all-at-once system is much
     * more sensitive near nonlinear convergence.
     */
    return fluid_env_int("FLUID_ST_PRESSURE_GAUGE", 1) ? 1 : 0;
}

static int fluid_is_pressure_gauge_owner(
    const FLUID_SUPS_SYSTEM* system)
{
    if(system == NULL ||
       !fluid_window_pressure_gauge_enabled())
    {
        return 0;
    }

    /*
     * MONOLIS local numbering places owned nodes first.  After ownership has
     * been resolved, local node zero on MPI rank zero is therefore a valid
     * globally unique gauge owner.
     */
    return
        monolis_mpi_get_global_my_rank() == 0 &&
        system->num_owned_space_nodes > 0;
}

static void fluid_apply_pressure_gauge_to_window_state(
    const FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* window_state)
{
    int slab;
    int alpha;

    if(system == NULL || layout == NULL || window_state == NULL ||
       !fluid_is_pressure_gauge_owner(system))
    {
        return;
    }

    for(slab = 0; slab < layout->num_slabs; slab++){
        for(alpha = 0; alpha < layout->time.n_dof; alpha++){
            const size_t index = ST_fom_vector_index(
                layout,
                0,
                slab,
                alpha,
                3);

            window_state[index] = 0.0;
        }
    }
}

static int fluid_initialize_window_bc(FLUID_SUPS_SYSTEM* system)
{
    ST_FOM_TIME_ELEMENT time;
    const int nnode = system->fe.total_num_nodes;
    const int slabs = system->vals.st_slabs_per_window;
    int window_dof;
    int node;
    int slab;
    int alpha;
    int component;
    int count = 0;

    if(ST_fom_time_element_init_dg(
        &time, system->vals.st_time_degree) != ST_FOM_SUCCESS)
    {
        return 1;
    }

    window_dof =
        FLUID_SUPS_PHYSICAL_DOF * slabs * time.n_dof;

    system->bc_NR_window.total_num_nodes = nnode;
    system->bc_NR_window.block_size = window_dof;
    BBFE_sys_memory_allocation_Dirichlet_bc(
        &system->bc_NR_window,
        nnode,
        window_dof);

    for(node = 0; node < nnode; node++){
        for(slab = 0; slab < slabs; slab++){
            for(alpha = 0; alpha < time.n_dof; alpha++){
                for(component = 0;
                    component < FLUID_SUPS_PHYSICAL_DOF;
                    component++)
                {
                    const int src =
                        FLUID_SUPS_PHYSICAL_DOF * node + component;
                    const int dst =
                        window_dof * node
                        + ((slab * time.n_dof + alpha)
                           * FLUID_SUPS_PHYSICAL_DOF)
                        + component;

                    system->bc_NR_window.D_bc_exists[dst] =
                        system->bc_NR.D_bc_exists[src];
                    system->bc_NR_window.imposed_D_val[dst] = 0.0;

                    if(system->bc_NR_window.D_bc_exists[dst]){
                        count++;
                    }
                }
            }
        }
    }

    /*
     * Pressure gauge for the Newton CORRECTION system.
     *
     * The physical Dirichlet file constrains velocity only.  For the
     * all-at-once incompressible system, fix one pressure coefficient at one
     * globally owned spatial node for every (slab,time-basis) pair.
     *
     * Since the imposed correction is zero, this only selects the pressure
     * gauge; it does not change velocity or pressure gradients.
     */
    if(fluid_is_pressure_gauge_owner(system)){
        for(slab = 0; slab < slabs; slab++){
            for(alpha = 0; alpha < time.n_dof; alpha++){
                const int dst =
                    window_dof * 0
                    + ((slab * time.n_dof + alpha)
                       * FLUID_SUPS_PHYSICAL_DOF)
                    + 3;

                if(!system->bc_NR_window.D_bc_exists[dst]){
                    system->bc_NR_window.D_bc_exists[dst] = true;
                    count++;
                }
                system->bc_NR_window.imposed_D_val[dst] = 0.0;
            }
        }

        printf(
            "%s pressure gauge: rank 0 local-owned node 0, "
            "%d pressure correction constraints "
            "(%d slabs x %d time dof).\\n",
            CODENAME,
            slabs * time.n_dof,
            slabs,
            time.n_dof);
    }

    system->bc_NR_window.num_D_bcs = count;
    system->window_dof = window_dof;
    return 0;
}


static int fluid_resolve_owned_space_nodes(
    FLUID_SUPS_SYSTEM* system)
{
    const int local_nodes =
        system != NULL ? system->fe.total_num_nodes : 0;
    const int communicator_size =
        monolis_mpi_get_global_comm_size();
    const int communication_owned =
        system != NULL ? system->mono_com.n_internal_vertex : 0;
    int owned = 0;

    if(system == NULL || local_nodes <= 0){
        return 1;
    }

    /*
     * MONOLIS/BBFE ownership contract used by the validated fluid and
     * convection-diffusion MPI paths:
     *
     *   serial:
     *       owned = all local nodes
     *
     *   MPI:
     *       owned = mono_com.n_internal_vertex
     *       local = fe.total_num_nodes = owned + halo
     *
     * In this BBFE initialization path system->monolis.mat.N may still be the
     * local (owned+halo) node count, so it must NOT be used as the MPI
     * ownership source.
     */
    if(ST_fom_resolve_owned_space_points(
        local_nodes,
        communicator_size,
        communication_owned,
        &owned) != ST_FOM_SUCCESS)
    {
        fprintf(
            stderr,
            "ERROR [rank %d]: failed to resolve spatial ownership: "
            "local=%d comm_size=%d n_internal_vertex=%d\n",
            monolis_mpi_get_global_my_rank(),
            local_nodes,
            communicator_size,
            communication_owned);
        return 1;
    }

    if(owned <= 0 || owned > local_nodes){
        fprintf(
            stderr,
            "ERROR [rank %d]: invalid resolved spatial ownership: "
            "owned=%d local=%d\n",
            monolis_mpi_get_global_my_rank(),
            owned,
            local_nodes);
        return 1;
    }

    system->num_owned_space_nodes = owned;

    printf(
        "%s spatial ownership rank %d/%d: "
        "owned=%d local=%d halo=%d "
        "(mono_com.n_internal_vertex=%d, base monolis.mat.N=%d)\n",
        CODENAME,
        monolis_mpi_get_global_my_rank(),
        communicator_size,
        owned,
        local_nodes,
        local_nodes - owned,
        communication_owned,
        system->monolis.mat.N);
    fflush(stdout);

    return 0;
}


static int fluid_initialize_window_parallel_pattern(
    FLUID_SUPS_SYSTEM* system)
{
    const int rank = monolis_mpi_get_global_my_rank();
    const int size = monolis_mpi_get_global_comm_size();
    const int local_nodes =
        system != NULL ? system->fe.total_num_nodes : 0;
    const int owned =
        system != NULL ? system->num_owned_space_nodes : 0;
    int status;

    if(system == NULL || local_nodes <= 0 ||
       owned <= 0 || owned > local_nodes ||
       system->window_dof <= 0)
    {
        return 1;
    }

    status = ST_fom_spatial_graph_initialize(
        &system->st_spatial_graph);
    if(status != ST_FOM_SUCCESS){
        return 1;
    }

    /*
     * This is deliberately the same construction used by the verified
     * convection-diffusion STSB FOM: reconstruct the local spatial nodal graph
     * from FE connectivity, then lift the graph by the complete window block
     * size.
     */
    status = ST_fom_spatial_graph_build_from_connectivity(
        &system->st_spatial_graph,
        local_nodes,
        system->fe.total_num_elems,
        system->fe.local_num_nodes,
        system->fe.conn);
    if(status != ST_FOM_SUCCESS){
        fprintf(
            stderr,
            "ERROR [rank %d]: failed to build fluid ST spatial graph "
            "(status=%d).\n",
            rank,
            status);
        return 1;
    }

    status = ST_fom_spatial_graph_validate_connectivity(
        &system->st_spatial_graph,
        local_nodes,
        system->fe.total_num_elems,
        system->fe.local_num_nodes,
        system->fe.conn);
    if(status != ST_FOM_SUCCESS){
        fprintf(
            stderr,
            "ERROR [rank %d]: fluid ST spatial graph does not contain "
            "all FE element couplings (status=%d).\n",
            rank,
            status);
        return 1;
    }

    monolis_initialize(&system->monolis_window);

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        &system->monolis_window,
        system->st_spatial_graph.num_nodes,
        system->window_dof,
        system->st_spatial_graph.index,
        system->st_spatial_graph.item);

    /*
     * Convection-diffusion STSB sets N/NP explicitly after lifting the local
     * spatial graph.  Do the same here.  The communicator itself remains
     * system->mono_com, i.e. the validated physical spatial decomposition.
     */
    system->monolis_window.mat.N = owned;
    system->monolis_window.mat.NP = local_nodes;

    if(size > 1 &&
       system->mono_com.n_internal_vertex != owned)
    {
        fprintf(
            stderr,
            "ERROR [rank %d]: window ownership/communicator mismatch: "
            "window N=%d communicator n_internal_vertex=%d.\n",
            rank,
            owned,
            system->mono_com.n_internal_vertex);
        return 1;
    }

    printf(
        "%s ST spatial MPI rank %d/%d: "
        "owned=%d local=%d halo=%d window_dof/node=%d "
        "window(N,NP)=(%d,%d)\n",
        CODENAME,
        rank,
        size,
        owned,
        local_nodes,
        local_nodes - owned,
        system->window_dof,
        system->monolis_window.mat.N,
        system->monolis_window.mat.NP);
    fflush(stdout);

    if(rank == 0){
        printf(
            "%s ST parallel decomposition: spatial MPI ranks=%d, "
            "temporal MPI ranks=1; all temporal DOFs in a window are "
            "stored as one block on each local spatial node.\n",
            CODENAME,
            size);
        printf(
            "%s Window halo exchange/linear solve communicator: "
            "validated base spatial MONOLIS communicator.\n",
            CODENAME);
    }

    return 0;
}


int fluid_sups_st_export_state(
    const FLUID_SUPS_SYSTEM* system,
    double* state)
{
    if(system == NULL || state == NULL){
        return 1;
    }
    fluid_pack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, state);
    return 0;
}

int fluid_sups_st_import_state(
    FLUID_SUPS_SYSTEM* system,
    const double* state,
    int update_previous_state)
{
    if(system == NULL || state == NULL){
        return 1;
    }
    fluid_unpack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, state);
    if(update_previous_state){
        fluid_copy_current_to_previous(system);
    }
    return 0;
}

int fluid_sups_st_synchronize_state(FLUID_SUPS_SYSTEM* system)
{
    double* vector;
    const int count = system != NULL
        ? system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF : 0;

    if(system == NULL || count <= 0){
        return 1;
    }

    vector = (double*)calloc((size_t)count, sizeof(double));
    if(vector == NULL){
        return 1;
    }

    fluid_pack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, vector);
    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        FLUID_SUPS_PHYSICAL_DOF,
        vector);
    fluid_unpack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, vector);
    free(vector);
    return 0;
}

static void fluid_output_result_file_vtk(
    BBFE_DATA* fe,
    FLUID_SUPS_VALUES* vals,
    const char* filename,
    const char* directory)
{
    FILE* fp = NULL;

    fp = BBFE_sys_write_fopen(fp, filename, directory);
    switch(fe->local_num_nodes){
    case 4:
        BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
        break;
    case 8:
        BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
        break;
    default:
        fprintf(stderr, "ERROR: unsupported element with %d nodes.\n",
            fe->local_num_nodes);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    BB_vtk_write_point_vals_vector(
        fp, vals->v, fe->total_num_nodes, "Velocity");
    BB_vtk_write_point_vals_scalar(
        fp, vals->p, fe->total_num_nodes, "Pressure");
    fclose(fp);
}

int fluid_sups_st_write_output(
    FLUID_SUPS_SYSTEM* system,
    int file_number,
    double time)
{
    char filename[FLUID_BUFFER_SIZE];
    const int rank = monolis_mpi_get_global_my_rank();

    (void)time;
    if(system == NULL || file_number < 0){
        return 1;
    }

    snprintf(filename, sizeof(filename),
        FLUID_OUTPUT_FILENAME_VTK, file_number, rank);
    fluid_output_result_file_vtk(
        &system->fe,
        &system->vals,
        filename,
        system->cond.directory);
    return 0;
}

int fluid_sups_st_write_output_named(
    FLUID_SUPS_SYSTEM* system,
    const char* prefix,
    int file_number,
    double time)
{
    char filename[FLUID_BUFFER_SIZE];
    const int rank = monolis_mpi_get_global_my_rank();
    int written;

    (void)time;

    if(system == NULL || prefix == NULL ||
       prefix[0] == '\0' || file_number < 0)
    {
        return 1;
    }

    written = snprintf(
        filename,
        sizeof(filename),
        "%s_%06d.vtk.%d",
        prefix,
        file_number,
        rank);

    if(written < 0 ||
       (size_t)written >= sizeof(filename))
    {
        return 1;
    }

    fluid_output_result_file_vtk(
        &system->fe,
        &system->vals,
        filename,
        system->cond.directory);

    return 0;
}

static void fluid_set_element_jacobian_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    FLUID_SUPS_VALUES*      vals)
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
                (void)h_e;
        
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

static void fluid_set_element_residuals_NR(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    FLUID_SUPS_VALUES*      vals)
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

static void fluid_set_element_residuals_mass_con(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals)
{
    /*
     * Diagnostic decomposition only.
     *
     * The repository contains the validated full NR residual kernel
     * BBFE_elemmat_fluid_sups_vec_NR(), but no definition of the historical
     * helper BBFE_elemmat_fluid_sups_vec_NR_mass_con().
     *
     * Assemble the exact same full residual used by Newton.  The subsequent
     * BBFE_fluid_calc_internal_norm_mass() reads only the continuity/pressure
     * component, so the nonlinear solve algebra is unchanged.
     */
    fluid_set_element_residuals_NR(monolis, fe, basis, vals);
}

static void fluid_set_element_residuals_momentum_eq(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals)
{
    /*
     * Diagnostic decomposition only.
     *
     * Assemble the exact full residual used by Newton.  The subsequent
     * BBFE_fluid_calc_internal_norm_momentum() reads only the three momentum
     * components.  This does not alter the Newton Jacobian, Newton RHS,
     * update, or convergence decision based on the total residual.
     */
    fluid_set_element_residuals_NR(monolis, fe, basis, vals);
}


static void fluid_add_window_current_jacobian_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int offset = FLUID_SUPS_PHYSICAL_DOF * slab;

    double*** val_ip =
        BB_std_calloc_3d_double(NULL, 4, 4, np);
    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double A[4][4];
    double J_inv[3][3];
    int e;
    int i;
    int j;
    int p;
    int a;
    int b;

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(
                nl, local_p, basis->N[p]);
        }

        for(i = 0; i < nl; i++){
            for(j = 0; j < nl; j++){
                for(p = 0; p < np; p++){
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };
                    double tau;
                    double tau_c;

                    BB_calc_mat3d_inverse(
                        fe->geo[e][p].J,
                        fe->geo[e][p].Jacobian,
                        J_inv);

                    tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        vals->density,
                        vals->viscosity,
                        v_ip[p],
                        vals->dt);
                    tau_c =
                        BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                            J_inv,
                            fe->geo[e][p].Jacobian,
                            vals->density,
                            vals->viscosity,
                            v_ip[p],
                            vals->dt);

                    BBFE_elemmat_fluid_sups_mat_NR(
                        A,
                        J_inv,
                        basis->N[p][i],
                        basis->N[p][j],
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        v_ip[p],
                        grad_v_ip[p],
                        grad_p_ip[p],
                        vals->density,
                        vals->viscosity,
                        tau,
                        tau_c,
                        vals->dt,
                        du_time);

                    for(a = 0; a < 4; a++){
                        for(b = 0; b < 4; b++){
                            val_ip[a][b][p] = A[a][b];
                            A[a][b] = 0.0;
                        }
                    }
                }

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 4; b++){
                        const double value =
                            BBFE_std_integ_calc(
                                np,
                                val_ip[a][b],
                                basis->integ_weight,
                                Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            offset + a,
                            offset + b,
                            value);
                    }
                }
            }
        }
    }

    BB_std_free_3d_double(val_ip, 4, 4, np);
    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
}

static void fluid_add_window_residual_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab,
    int window_dof)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int offset = FLUID_SUPS_PHYSICAL_DOF * slab;

    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** val_ip_vec =
        BB_std_calloc_2d_double(NULL, 4, np);
    double* integ_val_vec =
        BB_std_calloc_1d_double(NULL, 4);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double vec[4];
    double J_inv[3][3];
    int e;
    int i;
    int p;
    int d;

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] = BBFE_std_mapping_scalar(
                nl, local_p, basis->N[p]);
        }

        for(i = 0; i < nl; i++){
            for(p = 0; p < np; p++){
                double du_time[3] = {
                    v_ip[p][0] - v_ip_old[p][0],
                    v_ip[p][1] - v_ip_old[p][1],
                    v_ip[p][2] - v_ip_old[p][2]
                };
                double tau;
                double tau_c;

                BB_calc_mat3d_inverse(
                    fe->geo[e][p].J,
                    fe->geo[e][p].Jacobian,
                    J_inv);

                tau = BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);
                tau_c =
                    BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                        J_inv,
                        fe->geo[e][p].Jacobian,
                        vals->density,
                        vals->viscosity,
                        v_ip[p],
                        vals->dt);

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    basis->N[p][i],
                    fe->geo[e][p].grad_N[i],
                    v_ip[p],
                    v_ip_old[p],
                    grad_v_ip[p],
                    p_ip[p],
                    grad_p_ip[p],
                    vals->density,
                    vals->viscosity,
                    tau,
                    tau_c,
                    vals->dt,
                    du_time);

                for(d = 0; d < 4; d++){
                    val_ip_vec[d][p] = vec[d];
                }
            }

            for(d = 0; d < 4; d++){
                integ_val_vec[d] = BBFE_std_integ_calc(
                    np,
                    val_ip_vec[d],
                    basis->integ_weight,
                    Jacobian_ip);

                monolis->mat.R.B[
                    window_dof * fe->conn[e][i]
                    + offset + d] -= integ_val_vec[d];
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(val_ip_vec, 4, np);
    BB_std_free_1d_double(integ_val_vec, 4);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
}

typedef enum {
    FLUID_JOLD_ANALYTIC = 0,
    FLUID_JOLD_VERIFY   = 1,
    FLUID_JOLD_FD       = 2
} FLUID_JOLD_MODE;

typedef struct {
    int initialized;
    FLUID_JOLD_MODE mode;
    int verify_samples;
    int verify_stride;
    long long seen_count;
    int verified_count;
    int pass_announced;
    double verify_tolerance;
} FLUID_JOLD_CONTROL;

static FLUID_JOLD_CONTROL fluid_jold_control = {0};

static FLUID_JOLD_MODE fluid_jold_mode_from_environment(void)
{
    const char* text = getenv("FLUID_ST_JOLD_MODE");

    if(text == NULL || text[0] == '\0' ||
       strcmp(text, "verify") == 0)
    {
        return FLUID_JOLD_VERIFY;
    }

    if(strcmp(text, "analytic") == 0){
        return FLUID_JOLD_ANALYTIC;
    }

    if(strcmp(text, "fd") == 0 ||
       strcmp(text, "finite_difference") == 0)
    {
        return FLUID_JOLD_FD;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_JOLD_MODE must be "
        "verify, analytic, or fd; got '%s'.\n",
        text);
    exit(EXIT_FAILURE);
}

static void fluid_jold_initialize_control(void)
{
    if(fluid_jold_control.initialized){
        return;
    }

    fluid_jold_control.mode =
        fluid_jold_mode_from_environment();

    fluid_jold_control.verify_samples =
        fluid_env_int("FLUID_ST_JOLD_VERIFY_SAMPLES", 64);

    fluid_jold_control.verify_stride =
        fluid_env_int("FLUID_ST_JOLD_VERIFY_STRIDE", 5000);

    fluid_jold_control.verify_tolerance =
        fluid_env_double("FLUID_ST_JOLD_VERIFY_TOL", 5.0e-6);

    if(fluid_jold_control.verify_samples <= 0){
        fluid_jold_control.verify_samples = 1;
    }
    if(fluid_jold_control.verify_stride <= 0){
        fluid_jold_control.verify_stride = 1;
    }
    if(fluid_jold_control.verify_tolerance <= 0.0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_JOLD_VERIFY_TOL must be positive.\n");
        exit(EXIT_FAILURE);
    }

    fluid_jold_control.seen_count = 0;
    fluid_jold_control.verified_count = 0;
    fluid_jold_control.pass_announced = 0;
    fluid_jold_control.initialized = 1;

    if(monolis_mpi_get_global_my_rank() == 0){
        const char* mode_name =
            fluid_jold_control.mode == FLUID_JOLD_ANALYTIC
                ? "analytic"
                : (fluid_jold_control.mode == FLUID_JOLD_FD
                    ? "fd"
                    : "verify");

        printf(
            "%s J_old element kernel mode: %s",
            CODENAME,
            mode_name);

        if(fluid_jold_control.mode == FLUID_JOLD_VERIFY){
            printf(
                " (samples=%d stride=%d tolerance=%.3e; "
                "FD is used for assembly)",
                fluid_jold_control.verify_samples,
                fluid_jold_control.verify_stride,
                fluid_jold_control.verify_tolerance);
        }
        printf("\n");
    }
}

static void fluid_old_jacobian_ip_fd_reference(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    const double* velocity,
    const double* velocity_old,
    double** grad_velocity,
    double pressure,
    const double* grad_pressure,
    double tau,
    double tau_c,
    const double* du_time,
    double fd_epsilon,
    double J_old[4][3])
{
    double base[4];
    int a;
    int b;

    BBFE_elemmat_fluid_sups_vec_NR(
        base,
        basis->N[p][test_node],
        fe->geo[elem][p].grad_N[test_node],
        velocity,
        velocity_old,
        grad_velocity,
        pressure,
        grad_pressure,
        vals->density,
        vals->viscosity,
        tau,
        tau_c,
        vals->dt,
        du_time);

    for(b = 0; b < 3; b++){
        double old_perturbed[3] = {
            velocity_old[0],
            velocity_old[1],
            velocity_old[2]
        };
        double du_perturbed[3] = {
            du_time[0],
            du_time[1],
            du_time[2]
        };
        double perturbed[4];

        const double h =
            fd_epsilon
            * fmax(
                1.0,
                fmax(
                    fabs(velocity_old[b]),
                    fabs(velocity[b])));

        const double nodal_shape =
            basis->N[p][trial_node];

        old_perturbed[b] +=
            h * nodal_shape;
        du_perturbed[b] -=
            h * nodal_shape;

        BBFE_elemmat_fluid_sups_vec_NR(
            perturbed,
            basis->N[p][test_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            old_perturbed,
            grad_velocity,
            pressure,
            grad_pressure,
            vals->density,
            vals->viscosity,
            tau,
            tau_c,
            vals->dt,
            du_perturbed);

        for(a = 0; a < 4; a++){
            J_old[a][b] =
                (perturbed[a] - base[a]) / h;
        }
    }
}

static void fluid_print_jold_matrix(
    const char* label,
    const double J_old[4][3])
{
    int a;

    fprintf(stderr, "%s\n", label);
    for(a = 0; a < 4; a++){
        fprintf(
            stderr,
            "  [% .15e % .15e % .15e]\n",
            J_old[a][0],
            J_old[a][1],
            J_old[a][2]);
    }
}

static void fluid_verify_analytic_jold_sample(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    double tau,
    const double analytic[4][3],
    const double reference[4][3])
{
    double max_error = 0.0;
    double max_abs_difference = 0.0;
    int max_a = -1;
    int max_b = -1;
    int a;
    int b;

    if(fluid_jold_control.mode != FLUID_JOLD_VERIFY ||
       fluid_jold_control.verified_count >=
           fluid_jold_control.verify_samples)
    {
        return;
    }

    /*
     * Spread verification samples across the actual nonlinear run instead of
     * consuming all samples at the first element / first Newton state.
     */
    fluid_jold_control.seen_count++;

    if(((fluid_jold_control.seen_count - 1)
        % fluid_jold_control.verify_stride) != 0)
    {
        return;
    }

    for(a = 0; a < 4; a++){
        for(b = 0; b < 3; b++){
            const double difference =
                analytic[a][b] - reference[a][b];
            const double scale =
                fmax(
                    1.0,
                    fmax(
                        fabs(analytic[a][b]),
                        fabs(reference[a][b])));
            const double error =
                fabs(difference) / scale;

            if(error > max_error){
                max_error = error;
                max_abs_difference = fabs(difference);
                max_a = a;
                max_b = b;
            }
        }
    }

    if(!isfinite(max_error) ||
       max_error > fluid_jold_control.verify_tolerance)
    {
        fprintf(
            stderr,
            "\nERROR: analytic J_old does not match the linked "
            "BBFE_elemmat_fluid_sups_vec_NR provider.\n"
            "  element/test/trial/ip : %d / %d / %d / %d\n"
            "  sampled call          : %lld\n"
            "  max scaled error      : %.15e\n"
            "  max abs difference    : %.15e\n"
            "  offending entry       : (%d,%d)\n"
            "  N_i / N_j             : %.15e / %.15e\n"
            "  tau / dt              : %.15e / %.15e\n"
            "  density                : %.15e\n"
            "The Newton matrix is still using finite-difference J_old in "
            "verify mode, so use FLUID_ST_JOLD_MODE=fd to continue without "
            "analytic checks, or inspect the provider before enabling "
            "analytic mode.\n\n",
            elem,
            test_node,
            trial_node,
            p,
            fluid_jold_control.seen_count,
            max_error,
            max_abs_difference,
            max_a,
            max_b,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            tau,
            vals->dt,
            vals->density);

        fluid_print_jold_matrix(
            "analytic J_old:", analytic);
        fluid_print_jold_matrix(
            "finite-difference J_old:", reference);

        fflush(stderr);
        exit(EXIT_FAILURE);
    }

    fluid_jold_control.verified_count++;

    if(fluid_jold_control.verified_count >=
           fluid_jold_control.verify_samples &&
       !fluid_jold_control.pass_announced)
    {
        fluid_jold_control.pass_announced = 1;

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s analytic J_old verification PASS: "
                "%d sampled element kernels distributed over %lld calls, "
                "tolerance=%.3e. "
                "verify mode continued with FD assembly.\n",
                CODENAME,
                fluid_jold_control.verified_count,
                fluid_jold_control.seen_count,
                fluid_jold_control.verify_tolerance);
        }
    }
}

static void fluid_old_jacobian_ip(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    const double* velocity,
    const double* velocity_old,
    double** grad_velocity,
    double pressure,
    const double* grad_pressure,
    double tau,
    double tau_c,
    const double* du_time,
    double J_old[4][3])
{
    fluid_jold_initialize_control();

    if(fluid_jold_control.mode == FLUID_JOLD_ANALYTIC){
        BBFE_elemmat_fluid_sups_mat_NR_old(
            J_old,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            vals->density,
            tau,
            vals->dt);
        return;
    }

    /*
     * Both `fd` and `verify` use the already validated finite-difference
     * tangent for the ACTUAL Newton matrix.  This guarantees that turning on
     * analytic verification cannot alter convergence.
     */
    fluid_old_jacobian_ip_fd_reference(
        fe,
        basis,
        vals,
        elem,
        test_node,
        trial_node,
        p,
        velocity,
        velocity_old,
        grad_velocity,
        pressure,
        grad_pressure,
        tau,
        tau_c,
        du_time,
        vals->st_dg_jacobian_fd_epsilon,
        J_old);

    if(fluid_jold_control.mode == FLUID_JOLD_VERIFY){
        double analytic[4][3];

        BBFE_elemmat_fluid_sups_mat_NR_old(
            analytic,
            basis->N[p][test_node],
            basis->N[p][trial_node],
            fe->geo[elem][p].grad_N[test_node],
            velocity,
            vals->density,
            tau,
            vals->dt);

        fluid_verify_analytic_jold_sample(
            fe,
            basis,
            vals,
            elem,
            test_node,
            trial_node,
            p,
            tau,
            analytic,
            J_old);
    }
}

static void fluid_add_window_previous_jacobian_block(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int slab)
{
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const int row_offset =
        FLUID_SUPS_PHYSICAL_DOF * slab;
    const int col_offset =
        FLUID_SUPS_PHYSICAL_DOF * (slab - 1);

    double* Jacobian_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double*** grad_v_ip =
        BB_std_calloc_3d_double(NULL, np, 3, 3);
    double** local_v_old =
        BB_std_calloc_2d_double(NULL, nl, 3);
    double** v_ip_old =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* local_p =
        BB_std_calloc_1d_double(NULL, nl);
    double* p_ip =
        BB_std_calloc_1d_double(NULL, np);
    double** grad_p_ip =
        BB_std_calloc_2d_double(NULL, np, 3);
    double* tau =
        BB_std_calloc_1d_double(NULL, np);
    double* tau_c =
        BB_std_calloc_1d_double(NULL, np);
    double*** tangent_ip =
        BB_std_calloc_3d_double(NULL, 4, 3, np);

    int e;
    int i;
    int j;
    int p;
    int a;
    int b;

    if(slab <= 0){
        return;
    }

    for(e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(
            Jacobian_ip, np, e, fe);

        BBFE_elemmat_set_local_array_vector(
            local_v, fe, vals->v, e, 3);
        BBFE_elemmat_set_local_array_vector(
            local_v_old, fe, vals->v_old, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_p, fe, vals->p, e);

        for(p = 0; p < np; p++){
            double J_inv[3][3];

            BBFE_std_mapping_vector3d(
                v_ip[p], nl, local_v, basis->N[p]);
            BBFE_std_mapping_vector3d(
                v_ip_old[p], nl, local_v_old, basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_v_ip[p], nl, local_v,
                fe->geo[e][p].grad_N);
            BBFE_std_mapping_scalar_grad(
                grad_p_ip[p], nl, local_p,
                fe->geo[e][p].grad_N);
            p_ip[p] =
                BBFE_std_mapping_scalar(
                    nl, local_p, basis->N[p]);

            BB_calc_mat3d_inverse(
                fe->geo[e][p].J,
                fe->geo[e][p].Jacobian,
                J_inv);

            tau[p] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);

            tau_c[p] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    J_inv,
                    fe->geo[e][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    v_ip[p],
                    vals->dt);
        }

        for(i = 0; i < nl; i++){
            for(j = 0; j < nl; j++){
                for(p = 0; p < np; p++){
                    double J_old[4][3];
                    double du_time[3] = {
                        v_ip[p][0] - v_ip_old[p][0],
                        v_ip[p][1] - v_ip_old[p][1],
                        v_ip[p][2] - v_ip_old[p][2]
                    };

                    fluid_old_jacobian_ip(
                        fe,
                        basis,
                        vals,
                        e,
                        i,
                        j,
                        p,
                        v_ip[p],
                        v_ip_old[p],
                        grad_v_ip[p],
                        p_ip[p],
                        grad_p_ip[p],
                        tau[p],
                        tau_c[p],
                        du_time,
                        J_old);

                    for(a = 0; a < 4; a++){
                        for(b = 0; b < 3; b++){
                            tangent_ip[a][b][p] =
                                J_old[a][b];
                        }
                    }
                }

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 3; b++){
                        const double value =
                            BBFE_std_integ_calc(
                                np,
                                tangent_ip[a][b],
                                basis->integ_weight,
                                Jacobian_ip);

                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            fe->conn[e][i],
                            fe->conn[e][j],
                            row_offset + a,
                            col_offset + b,
                            value);
                    }
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_2d_double(local_v, nl, 3);
    BB_std_free_2d_double(v_ip, np, 3);
    BB_std_free_3d_double(grad_v_ip, np, 3, 3);
    BB_std_free_2d_double(local_v_old, nl, 3);
    BB_std_free_2d_double(v_ip_old, np, 3);
    BB_std_free_1d_double(local_p, nl);
    BB_std_free_1d_double(p_ip, np);
    BB_std_free_2d_double(grad_p_ip, np, 3);
    BB_std_free_1d_double(tau, np);
    BB_std_free_1d_double(tau_c, np);
    BB_std_free_3d_double(tangent_ip, 4, 3, np);
}

static int fluid_load_window_state_pair(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int slab,
    double* current_state,
    double* previous_state)
{
    if(fluid_unpack_window_slab_state(
        system, layout, window_state, slab, current_state) != 0)
    {
        return 1;
    }

    if(slab == 0){
        memcpy(
            previous_state,
            previous_trace,
            (size_t)system->fe.total_num_nodes
                * FLUID_SUPS_PHYSICAL_DOF
                * sizeof(double));
    }
    else{
        if(fluid_unpack_window_slab_state(
            system,
            layout,
            window_state,
            slab - 1,
            previous_state) != 0)
        {
            return 1;
        }
    }

    fluid_load_state_pair(
        system,
        current_state,
        previous_state);
    return 0;
}

static int fluid_assemble_window_residual_only(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    const int nstate =
        system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF;
    double* current_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    double* previous_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    int slab;

    if(current_state == NULL || previous_state == NULL){
        free(current_state);
        free(previous_state);
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(slab = 0; slab < layout->num_slabs; slab++){
        if(fluid_load_window_state_pair(
            system,
            layout,
            previous_trace,
            window_state,
            slab,
            current_state,
            previous_state) != 0)
        {
            free(current_state);
            free(previous_state);
            return 1;
        }

        fluid_add_window_residual_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab,
            system->window_dof);
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    free(current_state);
    free(previous_state);
    return 0;
}

static int fluid_assemble_window_newton_system(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    const int nstate =
        system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF;
    double* current_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    double* previous_state =
        (double*)calloc((size_t)nstate, sizeof(double));
    int slab;

    if(current_state == NULL || previous_state == NULL){
        free(current_state);
        free(previous_state);
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(slab = 0; slab < layout->num_slabs; slab++){
        if(fluid_load_window_state_pair(
            system,
            layout,
            previous_trace,
            window_state,
            slab,
            current_state,
            previous_state) != 0)
        {
            free(current_state);
            free(previous_state);
            return 1;
        }

        fluid_add_window_current_jacobian_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab);

        fluid_add_window_residual_block(
            &system->monolis_window,
            &system->fe,
            &system->basis,
            &system->vals,
            slab,
            system->window_dof);

        if(slab > 0){
            fluid_add_window_previous_jacobian_block(
                &system->monolis_window,
                &system->fe,
                &system->basis,
                &system->vals,
                slab);
        }
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    free(current_state);
    free(previous_state);
    return 0;
}

static void fluid_window_rhs_norms(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* total_norm,
    double* momentum_norm,
    double* mass_norm,
    double* slab_total_norm,
    double* slab_momentum_norm,
    double* slab_mass_norm)
{
    const int slabs = layout->num_slabs;
    const int owned = system->num_owned_space_nodes;
    double* reduced =
        (double*)calloc((size_t)3 * (size_t)slabs, sizeof(double));
    int node;
    int slab;
    int alpha;
    int component;

    for(node = 0; node < owned; node++){
        for(slab = 0; slab < slabs; slab++){
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                for(component = 0;
                    component < FLUID_SUPS_PHYSICAL_DOF;
                    component++)
                {
                    const int block_dof = ST_fom_block_dof(
                        layout, slab, alpha, component);
                    const int index =
                        system->window_dof * node + block_dof;
                    const double value =
                        system->monolis_window.mat.R.B[index];

                    reduced[3 * slab + 0] += value * value;
                    if(component < 3){
                        reduced[3 * slab + 1] += value * value;
                    }
                    else{
                        reduced[3 * slab + 2] += value * value;
                    }
                }
            }
        }
    }

    monolis_allreduce_R(
        3 * slabs,
        reduced,
        MONOLIS_MPI_SUM,
        system->mono_com.comm);

    *total_norm = 0.0;
    *momentum_norm = 0.0;
    *mass_norm = 0.0;

    for(slab = 0; slab < slabs; slab++){
        const double total2 =
            fmax(reduced[3 * slab + 0], 0.0);
        const double momentum2 =
            fmax(reduced[3 * slab + 1], 0.0);
        const double mass2 =
            fmax(reduced[3 * slab + 2], 0.0);

        *total_norm += total2;
        *momentum_norm += momentum2;
        *mass_norm += mass2;

        if(slab_total_norm != NULL){
            slab_total_norm[slab] = sqrt(total2);
        }
        if(slab_momentum_norm != NULL){
            slab_momentum_norm[slab] = sqrt(momentum2);
        }
        if(slab_mass_norm != NULL){
            slab_mass_norm[slab] = sqrt(mass2);
        }
    }

    *total_norm = sqrt(fmax(*total_norm, 0.0));
    *momentum_norm = sqrt(fmax(*momentum_norm, 0.0));
    *mass_norm = sqrt(fmax(*mass_norm, 0.0));
    free(reduced);
}

static double fluid_window_relative_norm(
    double value,
    double initial)
{
    return initial > 1.0e-300 ? value / initial : value;
}


static int fluid_window_linear_solver_type(void)
{
    const char* value = getenv("FLUID_ST_WINDOW_LINEAR_SOLVER");

    if(value == NULL || value[0] == '\0' ||
       strcmp(value, "bicgsafe") == 0)
    {
        return MONOLIS_ITER_BICGSAFE;
    }

    if(strcmp(value, "bicgstab") == 0){
        return MONOLIS_ITER_BICGSTAB;
    }

    fprintf(
        stderr,
        "ERROR: FLUID_ST_WINDOW_LINEAR_SOLVER must be "
        "bicgsafe or bicgstab; got '%s'.\n",
        value);
    exit(EXIT_FAILURE);
}

static const char* fluid_window_linear_solver_name(int solver)
{
    if(solver == MONOLIS_ITER_BICGSTAB){
        return "BiCGSTAB";
    }
    return "BiCGSAFE";
}

static double fluid_window_linear_epsilon(
    const FLUID_SUPS_SYSTEM* system)
{
    const double value = fluid_env_double(
        "FLUID_ST_WINDOW_LINEAR_EPSILON",
        system->vals.mat_epsilon);

    if(value <= 0.0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_WINDOW_LINEAR_EPSILON must be positive.\n");
        exit(EXIT_FAILURE);
    }

    return value;
}

static int fluid_window_linear_max_iter(
    const FLUID_SUPS_SYSTEM* system)
{
    const int value = fluid_env_int(
        "FLUID_ST_WINDOW_LINEAR_MAX_ITER",
        system->vals.mat_max_iter);

    if(value <= 0){
        fprintf(
            stderr,
            "ERROR: FLUID_ST_WINDOW_LINEAR_MAX_ITER must be positive.\n");
        exit(EXIT_FAILURE);
    }

    return value;
}

static void fluid_zero_window_linear_solution(
    FLUID_SUPS_SYSTEM* system)
{
    const size_t length =
        (size_t)system->fe.total_num_nodes
        * (size_t)system->window_dof;

    if(system != NULL &&
       system->monolis_window.mat.R.X != NULL)
    {
        memset(
            system->monolis_window.mat.R.X,
            0,
            length * sizeof(double));
    }
}

int fluid_sups_st_solve_window_all_at_once_dg0(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report)
{
    const int slabs = layout != NULL ? layout->num_slabs : 0;
    const int nstate = system != NULL
        ? system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF
        : 0;
    const size_t window_length =
        layout != NULL ? ST_fom_window_vector_length(layout) : 0;
    double* initial_slab_total = NULL;
    double* initial_slab_momentum = NULL;
    double* initial_slab_mass = NULL;
    double* final_slab_total = NULL;
    double* final_slab_momentum = NULL;
    double* final_slab_mass = NULL;
    double* final_state = NULL;
    double initial_total = 0.0;
    double initial_momentum = 0.0;
    double initial_mass = 0.0;
    double final_total = INFINITY;
    double final_momentum = INFINITY;
    double final_mass = INFINITY;
    double relative_total = INFINITY;
    double relative_momentum = INFINITY;
    double relative_mass = INFINITY;
    const int linear_solver =
        system != NULL ? fluid_window_linear_solver_type()
                       : MONOLIS_ITER_BICGSAFE;
    const double linear_epsilon =
        system != NULL ? fluid_window_linear_epsilon(system)
                       : 1.0e-8;
    const int linear_max_iter =
        system != NULL ? fluid_window_linear_max_iter(system)
                       : 10000;
    int iteration;
    int slab;
    int converged = 0;
    int status = 1;

    if(system == NULL || layout == NULL ||
       previous_trace == NULL || window_state == NULL ||
       slabs <= 0 || nstate <= 0 ||
       layout->time.degree != 0 ||
       layout->time.n_dof != 1 ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       layout->num_owned_space_points !=
           system->num_owned_space_nodes ||
       system->window_dof !=
           FLUID_SUPS_PHYSICAL_DOF * slabs)
    {
        return 1;
    }

    initial_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_state =
        (double*)calloc((size_t)nstate, sizeof(double));

    if(initial_slab_total == NULL ||
       initial_slab_momentum == NULL ||
       initial_slab_mass == NULL ||
       final_slab_total == NULL ||
       final_slab_momentum == NULL ||
       final_slab_mass == NULL ||
       final_state == NULL)
    {
        goto cleanup;
    }

    if(window_report != NULL){
        memset(window_report, 0, sizeof(*window_report));
        window_report->relative_total_residual = INFINITY;
        window_report->relative_momentum_residual = INFINITY;
        window_report->relative_mass_residual = INFINITY;
    }

    if(slab_reports != NULL){
        memset(
            slab_reports,
            0,
            (size_t)slabs * sizeof(*slab_reports));
    }

    /*
     * Constant predictor across the whole window.  No slab is separately
     * solved before the window Newton iteration starts.
     */
    memset(window_state, 0, window_length * sizeof(double));
    for(slab = 0; slab < slabs; slab++){
        if(ST_fom_pack_dg0_slab(
            layout,
            slab,
            previous_trace,
            window_state) != ST_FOM_SUCCESS)
        {
            goto cleanup;
        }
    }
    fluid_apply_physical_bc_to_window(
        system, layout, window_state);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        window_state);

    for(iteration = 0;
        iteration < system->vals.nr_max_iter;
        iteration++)
    {
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\n%s ================ ALL-AT-ONCE ST window %d "
                "Newton %d ================\n",
                CODENAME,
                window_id,
                iteration + 1);
        }

        /*
         * Assemble the COMPLETE window Jacobian and residual, then call
         * MONOLIS exactly once for this Newton iteration.
         */
        if(fluid_assemble_window_newton_system(
            system,
            layout,
            previous_trace,
            window_state) != 0)
        {
            goto cleanup;
        }

        if(iteration == 0){
            fluid_window_rhs_norms(
                system,
                layout,
                &initial_total,
                &initial_momentum,
                &initial_mass,
                initial_slab_total,
                initial_slab_momentum,
                initial_slab_mass);
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s window Newton linear solve: "
                "%d dof/node = %d physical dof x %d slabs; "
                "ONE MONOLIS solve; solver=%s eps=%.3e max_iter=%d\n",
                CODENAME,
                system->window_dof,
                FLUID_SUPS_PHYSICAL_DOF,
                slabs,
                fluid_window_linear_solver_name(linear_solver),
                linear_epsilon,
                linear_max_iter);
        }

        /*
         * Use a zero correction initial guess every Newton iteration.
         * Reusing the previous Newton correction is especially harmful when
         * the nonlinear residual has already dropped by several orders.
         */
        fluid_zero_window_linear_solution(system);

        BBFE_sys_monowrap_solve(
            &system->monolis_window,
            &system->mono_com,
            system->monolis_window.mat.R.X,
            linear_solver,
            MONOLIS_PREC_DIAG,
            linear_max_iter,
            linear_epsilon);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            system->monolis_window.mat.R.X);

        {
            size_t row;
            for(row = 0; row < window_length; row++){
                window_state[row] +=
                    system->monolis_window.mat.R.X[row];
            }
        }

        fluid_apply_physical_bc_to_window(
            system, layout, window_state);
        fluid_apply_pressure_gauge_to_window_state(
            system, layout, window_state);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            window_state);

        if(fluid_assemble_window_residual_only(
            system,
            layout,
            previous_trace,
            window_state) != 0)
        {
            goto cleanup;
        }

        fluid_window_rhs_norms(
            system,
            layout,
            &final_total,
            &final_momentum,
            &final_mass,
            final_slab_total,
            final_slab_momentum,
            final_slab_mass);

        relative_total =
            fluid_window_relative_norm(
                final_total, initial_total);
        relative_momentum =
            fluid_window_relative_norm(
                final_momentum, initial_momentum);
        relative_mass =
            fluid_window_relative_norm(
                final_mass, initial_mass);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s ALL-AT-ONCE window %d Newton %d residual: "
                "total=%e momentum=%e mass=%e\n",
                CODENAME,
                window_id,
                iteration + 1,
                relative_total,
                relative_momentum,
                relative_mass);
        }

        if(isfinite(relative_total) &&
           relative_total < system->vals.nr_epsilon)
        {
            converged = 1;
            break;
        }
    }

    if(!converged){
        fprintf(
            stderr,
            "ERROR: all-at-once window Newton did not converge: "
            "window=%d rel_total=%.15e\n",
            window_id,
            relative_total);
        goto cleanup;
    }

    if(slab_reports != NULL){
        for(slab = 0; slab < slabs; slab++){
            slab_reports[slab].converged = 1;
            slab_reports[slab].nonlinear_iterations =
                iteration + 1;
            slab_reports[slab].relative_total_residual =
                fluid_window_relative_norm(
                    final_slab_total[slab],
                    initial_slab_total[slab]);
            slab_reports[slab].relative_momentum_residual =
                fluid_window_relative_norm(
                    final_slab_momentum[slab],
                    initial_slab_momentum[slab]);
            slab_reports[slab].relative_mass_residual =
                fluid_window_relative_norm(
                    final_slab_mass[slab],
                    initial_slab_mass[slab]);
        }
    }

    if(window_report != NULL){
        window_report->converged = 1;
        window_report->nonlinear_iterations =
            iteration + 1;
        window_report->relative_total_residual =
            relative_total;
        window_report->relative_momentum_residual =
            relative_momentum;
        window_report->relative_mass_residual =
            relative_mass;
    }

    if(ST_fom_unpack_dg0_slab(
        layout,
        slabs - 1,
        window_state,
        final_state) != ST_FOM_SUCCESS)
    {
        goto cleanup;
    }

    if(fluid_sups_st_import_state(
        system,
        final_state,
        1) != 0 ||
       fluid_sups_st_synchronize_state(system) != 0)
    {
        goto cleanup;
    }

    status = 0;

cleanup:
    free(initial_slab_total);
    free(initial_slab_momentum);
    free(initial_slab_mass);
    free(final_slab_total);
    free(final_slab_momentum);
    free(final_slab_mass);
    free(final_state);
    return status;
}


/* ========================================================================== */
/* Higher-order dG(1)/dG(2) temporal weak form                               */
/* ========================================================================== */

static double fluid_time_basis_derivative(
    const ST_FOM_TIME_ELEMENT* time,
    int beta,
    double tau)
{
    if(time == NULL || beta < 0 || beta >= time->n_dof){
        return 0.0;
    }

    if(time->degree == 0){
        return 0.0;
    }
    if(time->degree == 1){
        return beta == 0 ? -1.0 : 1.0;
    }
    if(time->degree == 2){
        if(beta == 0){
            return 4.0 * tau - 3.0;
        }
        if(beta == 1){
            return -8.0 * tau + 4.0;
        }
        return 4.0 * tau - 1.0;
    }
    return 0.0;
}

static size_t fluid_local_time_index(
    int nt,
    int nl,
    int beta,
    int local_node,
    int component)
{
    (void)nt;
    return (((size_t)beta * (size_t)nl + (size_t)local_node)
        * (size_t)FLUID_SUPS_PHYSICAL_DOF)
        + (size_t)component;
}

static size_t fluid_local_residual_index(
    int nl,
    int alpha,
    int local_node,
    int component)
{
    return (((size_t)alpha * (size_t)nl + (size_t)local_node)
        * (size_t)FLUID_SUPS_PHYSICAL_DOF)
        + (size_t)component;
}

typedef struct {
    int nl;
    int np;
    int nq;
    int nt;

    double* Jacobian_ip;
    double (*J_inv)[3][3];

    /* [time quadrature][space quadrature] */
    double (*v)[3];
    double (*du_time)[3];
    double (*grad_v)[3][3];
    double* pressure;
    double (*grad_pressure)[3];
    double* tau;
    double* tau_c;

    /* left trace and previous right trace, at space quadrature points */
    double (*v_left)[3];
    double (*grad_v_left)[3][3];
    double* pressure_left;
    double (*grad_pressure_left)[3];
    double (*v_previous)[3];
    double* tau_left;
    double* tau_c_left;

    /* reusable spatial nodal work */
    double** local_v;
    double** local_du;
    double* local_p;
    double** local_v_left;
    double* local_p_left;
    double** local_v_previous;
} FLUID_DG_ELEMENT_FIELDS;

static void fluid_dg_element_fields_clear(
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    if(fields != NULL){
        memset(fields, 0, sizeof(*fields));
    }
}

static int fluid_dg_element_fields_allocate(
    FLUID_DG_ELEMENT_FIELDS* fields,
    int nl,
    int np,
    int nq,
    int nt)
{
    const size_t nqp = (size_t)nq * (size_t)np;

    if(fields == NULL || nl <= 0 || np <= 0 ||
       nq <= 0 || nt <= 0)
    {
        return 1;
    }

    fluid_dg_element_fields_clear(fields);
    fields->nl = nl;
    fields->np = np;
    fields->nq = nq;
    fields->nt = nt;

    fields->Jacobian_ip =
        (double*)calloc((size_t)np, sizeof(double));
    fields->J_inv =
        (double (*)[3][3])calloc(
            (size_t)np, sizeof(*fields->J_inv));

    fields->v =
        (double (*)[3])calloc(nqp, sizeof(*fields->v));
    fields->du_time =
        (double (*)[3])calloc(nqp, sizeof(*fields->du_time));
    fields->grad_v =
        (double (*)[3][3])calloc(
            nqp, sizeof(*fields->grad_v));
    fields->pressure =
        (double*)calloc(nqp, sizeof(double));
    fields->grad_pressure =
        (double (*)[3])calloc(
            nqp, sizeof(*fields->grad_pressure));
    fields->tau =
        (double*)calloc(nqp, sizeof(double));
    fields->tau_c =
        (double*)calloc(nqp, sizeof(double));

    fields->v_left =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->v_left));
    fields->grad_v_left =
        (double (*)[3][3])calloc(
            (size_t)np, sizeof(*fields->grad_v_left));
    fields->pressure_left =
        (double*)calloc((size_t)np, sizeof(double));
    fields->grad_pressure_left =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->grad_pressure_left));
    fields->v_previous =
        (double (*)[3])calloc(
            (size_t)np, sizeof(*fields->v_previous));
    fields->tau_left =
        (double*)calloc((size_t)np, sizeof(double));
    fields->tau_c_left =
        (double*)calloc((size_t)np, sizeof(double));

    fields->local_v =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_du =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_p =
        BB_std_calloc_1d_double(NULL, nl);
    fields->local_v_left =
        BB_std_calloc_2d_double(NULL, nl, 3);
    fields->local_p_left =
        BB_std_calloc_1d_double(NULL, nl);
    fields->local_v_previous =
        BB_std_calloc_2d_double(NULL, nl, 3);

    if(fields->Jacobian_ip == NULL ||
       fields->J_inv == NULL ||
       fields->v == NULL ||
       fields->du_time == NULL ||
       fields->grad_v == NULL ||
       fields->pressure == NULL ||
       fields->grad_pressure == NULL ||
       fields->tau == NULL ||
       fields->tau_c == NULL ||
       fields->v_left == NULL ||
       fields->grad_v_left == NULL ||
       fields->pressure_left == NULL ||
       fields->grad_pressure_left == NULL ||
       fields->v_previous == NULL ||
       fields->tau_left == NULL ||
       fields->tau_c_left == NULL ||
       fields->local_v == NULL ||
       fields->local_du == NULL ||
       fields->local_p == NULL ||
       fields->local_v_left == NULL ||
       fields->local_p_left == NULL ||
       fields->local_v_previous == NULL)
    {
        return 1;
    }

    return 0;
}

static void fluid_dg_element_fields_free(
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    if(fields == NULL){
        return;
    }

    free(fields->Jacobian_ip);
    free(fields->J_inv);
    free(fields->v);
    free(fields->du_time);
    free(fields->grad_v);
    free(fields->pressure);
    free(fields->grad_pressure);
    free(fields->tau);
    free(fields->tau_c);
    free(fields->v_left);
    free(fields->grad_v_left);
    free(fields->pressure_left);
    free(fields->grad_pressure_left);
    free(fields->v_previous);
    free(fields->tau_left);
    free(fields->tau_c_left);

    if(fields->local_v != NULL){
        BB_std_free_2d_double(fields->local_v, fields->nl, 3);
    }
    if(fields->local_du != NULL){
        BB_std_free_2d_double(fields->local_du, fields->nl, 3);
    }
    if(fields->local_p != NULL){
        BB_std_free_1d_double(fields->local_p, fields->nl);
    }
    if(fields->local_v_left != NULL){
        BB_std_free_2d_double(
            fields->local_v_left, fields->nl, 3);
    }
    if(fields->local_p_left != NULL){
        BB_std_free_1d_double(
            fields->local_p_left, fields->nl);
    }
    if(fields->local_v_previous != NULL){
        BB_std_free_2d_double(
            fields->local_v_previous, fields->nl, 3);
    }

    fluid_dg_element_fields_clear(fields);
}

static int fluid_dg_build_element_fields(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    const ST_FOM_TIME_ELEMENT* time,
    int elem,
    const double* coefficient,
    const double* previous_trace_local,
    FLUID_DG_ELEMENT_FIELDS* fields)
{
    const int nl = fields->nl;
    const int np = fields->np;
    const int nq = fields->nq;
    const int nt = fields->nt;
    int q;
    int p;
    int j;
    int beta;
    int d;

    if(fe == NULL || basis == NULL || vals == NULL ||
       time == NULL || coefficient == NULL ||
       previous_trace_local == NULL || fields == NULL)
    {
        return 1;
    }

    BBFE_elemmat_set_Jacobian_array(
        fields->Jacobian_ip, np, elem, fe);

    for(p = 0; p < np; p++){
        BB_calc_mat3d_inverse(
            fe->geo[elem][p].J,
            fe->geo[elem][p].Jacobian,
            fields->J_inv[p]);
    }

    /*
     * Temporal quadrature-point state:
     *
     *   v_q       = sum_beta psi_beta v_beta
     *   du_time_q = sum_beta dpsi_beta/dtau v_beta
     *
     * Since physical dv/dt = du_time_q/dt, the verified BE kernel can be
     * evaluated at the polynomial derivative by setting
     *
     *   v_old,eff = v_q - du_time_q.
     */
    for(q = 0; q < nq; q++){
        const double tau_time = time->tau[q];

        for(j = 0; j < nl; j++){
            for(d = 0; d < 3; d++){
                double value = 0.0;
                double derivative = 0.0;

                for(beta = 0; beta < nt; beta++){
                    const double coeff =
                        coefficient[
                            fluid_local_time_index(
                                nt, nl, beta, j, d)];
                    value += time->psi[beta][q] * coeff;
                    derivative +=
                        fluid_time_basis_derivative(
                            time, beta, tau_time) * coeff;
                }

                fields->local_v[j][d] = value;
                fields->local_du[j][d] = derivative;
            }

            fields->local_p[j] = 0.0;
            for(beta = 0; beta < nt; beta++){
                fields->local_p[j] +=
                    time->psi[beta][q]
                    * coefficient[
                        fluid_local_time_index(
                            nt, nl, beta, j, 3)];
            }
        }

        for(p = 0; p < np; p++){
            const size_t qp =
                (size_t)q * (size_t)np + (size_t)p;
            double grad_storage[3][3];
            double* grad_rows[3] = {
                grad_storage[0],
                grad_storage[1],
                grad_storage[2]
            };
            double grad_p[3];

            BBFE_std_mapping_vector3d(
                fields->v[qp],
                nl,
                fields->local_v,
                basis->N[p]);
            BBFE_std_mapping_vector3d(
                fields->du_time[qp],
                nl,
                fields->local_du,
                basis->N[p]);
            BBFE_std_mapping_vector3d_grad(
                grad_rows,
                nl,
                fields->local_v,
                fe->geo[elem][p].grad_N);

            for(d = 0; d < 3; d++){
                int k;
                fields->grad_pressure[qp][d] = 0.0;
                for(k = 0; k < 3; k++){
                    fields->grad_v[qp][d][k] =
                        grad_storage[d][k];
                }
            }

            fields->pressure[qp] =
                BBFE_std_mapping_scalar(
                    nl,
                    fields->local_p,
                    basis->N[p]);
            BBFE_std_mapping_scalar_grad(
                grad_p,
                nl,
                fields->local_p,
                fe->geo[elem][p].grad_N);
            for(d = 0; d < 3; d++){
                fields->grad_pressure[qp][d] = grad_p[d];
            }

            fields->tau[qp] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor(
                    fields->J_inv[p],
                    fe->geo[elem][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    fields->v[qp],
                    vals->dt);
            fields->tau_c[qp] =
                BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                    fields->J_inv[p],
                    fe->geo[elem][p].Jacobian,
                    vals->density,
                    vals->viscosity,
                    fields->v[qp],
                    vals->dt);
        }
    }

    /* Left trace of current slab and previous right trace. */
    for(j = 0; j < nl; j++){
        for(d = 0; d < 3; d++){
            fields->local_v_left[j][d] = 0.0;
            for(beta = 0; beta < nt; beta++){
                fields->local_v_left[j][d] +=
                    time->e_minus[beta]
                    * coefficient[
                        fluid_local_time_index(
                            nt, nl, beta, j, d)];
            }
            fields->local_v_previous[j][d] =
                previous_trace_local[
                    (size_t)j * FLUID_SUPS_PHYSICAL_DOF
                    + (size_t)d];
        }

        fields->local_p_left[j] = 0.0;
        for(beta = 0; beta < nt; beta++){
            fields->local_p_left[j] +=
                time->e_minus[beta]
                * coefficient[
                    fluid_local_time_index(
                        nt, nl, beta, j, 3)];
        }
    }

    for(p = 0; p < np; p++){
        double grad_storage[3][3];
        double* grad_rows[3] = {
            grad_storage[0],
            grad_storage[1],
            grad_storage[2]
        };
        double grad_p[3];

        BBFE_std_mapping_vector3d(
            fields->v_left[p],
            nl,
            fields->local_v_left,
            basis->N[p]);
        BBFE_std_mapping_vector3d_grad(
            grad_rows,
            nl,
            fields->local_v_left,
            fe->geo[elem][p].grad_N);
        for(d = 0; d < 3; d++){
            int k;
            for(k = 0; k < 3; k++){
                fields->grad_v_left[p][d][k] =
                    grad_storage[d][k];
            }
        }

        fields->pressure_left[p] =
            BBFE_std_mapping_scalar(
                nl,
                fields->local_p_left,
                basis->N[p]);
        BBFE_std_mapping_scalar_grad(
            grad_p,
            nl,
            fields->local_p_left,
            fe->geo[elem][p].grad_N);
        for(d = 0; d < 3; d++){
            fields->grad_pressure_left[p][d] = grad_p[d];
        }

        BBFE_std_mapping_vector3d(
            fields->v_previous[p],
            nl,
            fields->local_v_previous,
            basis->N[p]);

        fields->tau_left[p] =
            BBFE_elemmat_fluid_sups_coef_metric_tensor(
                fields->J_inv[p],
                fe->geo[elem][p].Jacobian,
                vals->density,
                vals->viscosity,
                fields->v_left[p],
                vals->dt);
        fields->tau_c_left[p] =
            BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
                fields->J_inv[p],
                fe->geo[elem][p].Jacobian,
                vals->density,
                vals->viscosity,
                fields->v_left[p],
                vals->dt);
    }

    return 0;
}

/*
 * Directional derivative of the verified residual with respect to an old
 * velocity nodal coefficient.  A BE old-state perturbation produces
 *
 *   delta v_old_ip = N_j delta,
 *   delta du_time  = -N_j delta.
 *
 * Hence this derivative includes every transient SUPS/PSPG contribution
 * present in the verified kernel, without replacing it by -M/dt.
 */
static void fluid_dg_old_jacobian_ip(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    FLUID_SUPS_VALUES* vals,
    int elem,
    int test_node,
    int trial_node,
    int p,
    const double* velocity,
    const double* velocity_old,
    double** grad_velocity,
    double pressure,
    const double* grad_pressure,
    double tau,
    double tau_c,
    const double* du_time,
    double fd_epsilon,
    double J_old[4][3])
{
    (void)fd_epsilon;

    fluid_old_jacobian_ip(
        fe,
        basis,
        vals,
        elem,
        test_node,
        trial_node,
        p,
        velocity,
        velocity_old,
        grad_velocity,
        pressure,
        grad_pressure,
        tau,
        tau_c,
        du_time,
        J_old);
}

static int fluid_dg_gather_element_coefficients(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int elem,
    int slab,
    double* coefficient,
    double* previous_trace_local)
{
    const int nl = system->fe.local_num_nodes;
    const int nt = layout->time.n_dof;
    int j;
    int beta;
    int component;

    for(j = 0; j < nl; j++){
        const int node = system->fe.conn[elem][j];

        for(beta = 0; beta < nt; beta++){
            for(component = 0;
                component < FLUID_SUPS_PHYSICAL_DOF;
                component++)
            {
                coefficient[
                    fluid_local_time_index(
                        nt, nl, beta, j, component)] =
                    window_state[
                        ST_fom_vector_index(
                            layout,
                            node,
                            slab,
                            beta,
                            component)];
            }
        }

        for(component = 0;
            component < FLUID_SUPS_PHYSICAL_DOF;
            component++)
        {
            double value = 0.0;

            if(slab == 0){
                value =
                    previous_trace[
                        (size_t)node
                            * FLUID_SUPS_PHYSICAL_DOF
                        + (size_t)component];
            }
            else{
                for(beta = 0; beta < nt; beta++){
                    value +=
                        layout->time.e_plus[beta]
                        * window_state[
                            ST_fom_vector_index(
                                layout,
                                node,
                                slab - 1,
                                beta,
                                component)];
                }
            }

            previous_trace_local[
                (size_t)j * FLUID_SUPS_PHYSICAL_DOF
                + (size_t)component] = value;
        }
    }

    return 0;
}

static int fluid_dg_add_element_slab(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int elem,
    int slab,
    int assemble_jacobian)
{
    const int nl = system->fe.local_num_nodes;
    const int np = system->basis.num_integ_points;
    const int nt = layout->time.n_dof;
    const int nq = layout->time.num_time_ip;
    const int local_dof =
        nt * nl * FLUID_SUPS_PHYSICAL_DOF;
    const size_t residual_ip_size =
        (size_t)nt * FLUID_SUPS_PHYSICAL_DOF * (size_t)np;
    const size_t jac_ip_size =
        (size_t)nt * (size_t)nt
        * FLUID_SUPS_PHYSICAL_DOF
        * FLUID_SUPS_PHYSICAL_DOF
        * (size_t)np;
    const size_t prev_jac_ip_size =
        (size_t)nt * (size_t)nt
        * FLUID_SUPS_PHYSICAL_DOF
        * 3u
        * (size_t)np;

    double* coefficient = NULL;
    double* previous_trace_local = NULL;
    double* residual_ip = NULL;
    double* jump_ip = NULL;
    double* jac_ip = NULL;
    double* prev_jac_ip = NULL;
    FLUID_DG_ELEMENT_FIELDS fields;
    int i;
    int j;
    int p;
    int q;
    int alpha;
    int beta;
    int a;
    int b;
    int status = 1;

    fluid_dg_element_fields_clear(&fields);

    coefficient =
        (double*)calloc((size_t)local_dof, sizeof(double));
    previous_trace_local =
        (double*)calloc(
            (size_t)nl * FLUID_SUPS_PHYSICAL_DOF,
            sizeof(double));
    residual_ip =
        (double*)calloc(residual_ip_size, sizeof(double));
    jump_ip =
        (double*)calloc(
            (size_t)FLUID_SUPS_PHYSICAL_DOF * (size_t)np,
            sizeof(double));

    if(assemble_jacobian){
        jac_ip =
            (double*)calloc(jac_ip_size, sizeof(double));
        if(slab > 0){
            prev_jac_ip =
                (double*)calloc(
                    prev_jac_ip_size, sizeof(double));
        }
    }

    if(coefficient == NULL ||
       previous_trace_local == NULL ||
       residual_ip == NULL ||
       jump_ip == NULL ||
       (assemble_jacobian && jac_ip == NULL) ||
       (assemble_jacobian && slab > 0 && prev_jac_ip == NULL) ||
       fluid_dg_element_fields_allocate(
           &fields, nl, np, nq, nt) != 0)
    {
        goto cleanup;
    }

    fluid_dg_gather_element_coefficients(
        system,
        layout,
        previous_trace,
        window_state,
        elem,
        slab,
        coefficient,
        previous_trace_local);

    if(fluid_dg_build_element_fields(
        &system->fe,
        &system->basis,
        &system->vals,
        &layout->time,
        elem,
        coefficient,
        previous_trace_local,
        &fields) != 0)
    {
        goto cleanup;
    }

#define FLUID_RES_IP(alpha_, a_, p_) \
    residual_ip[(((size_t)(alpha_) * 4u + (size_t)(a_)) \
        * (size_t)np) + (size_t)(p_)]
#define FLUID_JUMP_IP(a_, p_) \
    jump_ip[((size_t)(a_) * (size_t)np) + (size_t)(p_)]
#define FLUID_JAC_IP(alpha_, beta_, a_, b_, p_) \
    jac_ip[((((((size_t)(alpha_) * (size_t)nt \
        + (size_t)(beta_)) * 4u + (size_t)(a_)) * 4u \
        + (size_t)(b_)) * (size_t)np) + (size_t)(p_))]
#define FLUID_PREV_JAC_IP(alpha_, beta_, a_, b_, p_) \
    prev_jac_ip[((((((size_t)(alpha_) * (size_t)nt \
        + (size_t)(beta_)) * 4u + (size_t)(a_)) * 3u \
        + (size_t)(b_)) * (size_t)np) + (size_t)(p_))]

    for(i = 0; i < nl; i++){
        memset(
            residual_ip,
            0,
            residual_ip_size * sizeof(double));
        memset(
            jump_ip,
            0,
            (size_t)FLUID_SUPS_PHYSICAL_DOF
                * (size_t)np
                * sizeof(double));

        /*
         * Interior temporal weak residual:
         *
         *   int_0^1 psi_alpha R_SUPS(U(tau), U_t(tau)) d tau.
         */
        for(q = 0; q < nq; q++){
            for(p = 0; p < np; p++){
                const size_t qp =
                    (size_t)q * (size_t)np + (size_t)p;
                double old_effective[3] = {
                    fields.v[qp][0] - fields.du_time[qp][0],
                    fields.v[qp][1] - fields.du_time[qp][1],
                    fields.v[qp][2] - fields.du_time[qp][2]
                };
                double grad_storage[3][3];
                double* grad_rows[3] = {
                    grad_storage[0],
                    grad_storage[1],
                    grad_storage[2]
                };
                double vec[4];
                int d;

                for(a = 0; a < 3; a++){
                    for(b = 0; b < 3; b++){
                        grad_storage[a][b] =
                            fields.grad_v[qp][a][b];
                    }
                }

                BBFE_elemmat_fluid_sups_vec_NR(
                    vec,
                    system->basis.N[p][i],
                    system->fe.geo[elem][p].grad_N[i],
                    fields.v[qp],
                    old_effective,
                    grad_rows,
                    fields.pressure[qp],
                    fields.grad_pressure[qp],
                    system->vals.density,
                    system->vals.viscosity,
                    fields.tau[qp],
                    fields.tau_c[qp],
                    system->vals.dt,
                    fields.du_time[qp]);

                for(alpha = 0; alpha < nt; alpha++){
                    const double weight =
                        layout->time.weight[q]
                        * layout->time.psi[alpha][q];

                    for(d = 0; d < 4; d++){
                        FLUID_RES_IP(alpha, d, p) +=
                            weight * vec[d];
                    }
                }
            }
        }

        /*
         * Jacobian and dG upwind jump.  The previous-state derivative J_old
         * is used in two ways:
         *
         *  1. temporal chain rule in the slab interior,
         *  2. transient jump operator at tau=0.
         *
         * For a standard mass term J_old=-M/dt, so
         *
         *   -J_old (U^- - U_prev^+) = M/dt (U^- - U_prev^+).
         *
         * Here the same operation also retains SUPS/PSPG transient pieces.
         */
        for(j = 0; j < nl; j++){
            if(assemble_jacobian){
                memset(jac_ip, 0, jac_ip_size * sizeof(double));
                if(slab > 0){
                    memset(
                        prev_jac_ip,
                        0,
                        prev_jac_ip_size * sizeof(double));
                }
            }

            for(p = 0; p < np; p++){
                double grad_left_storage[3][3];
                double* grad_left_rows[3] = {
                    grad_left_storage[0],
                    grad_left_storage[1],
                    grad_left_storage[2]
                };
                double zero_du[3] = {0.0, 0.0, 0.0};
                double J_old_left[4][3];

                for(a = 0; a < 3; a++){
                    for(b = 0; b < 3; b++){
                        grad_left_storage[a][b] =
                            fields.grad_v_left[p][a][b];
                    }
                }

                fluid_dg_old_jacobian_ip(
                    &system->fe,
                    &system->basis,
                    &system->vals,
                    elem,
                    i,
                    j,
                    p,
                    fields.v_left[p],
                    fields.v_left[p],
                    grad_left_rows,
                    fields.pressure_left[p],
                    fields.grad_pressure_left[p],
                    fields.tau_left[p],
                    fields.tau_c_left[p],
                    zero_du,
                    system->vals.st_dg_jacobian_fd_epsilon,
                    J_old_left);

                for(a = 0; a < 4; a++){
                    for(b = 0; b < 3; b++){
                        const double jump_nodal =
                            fields.local_v_left[j][b]
                            - fields.local_v_previous[j][b];
                        FLUID_JUMP_IP(a, p) -=
                            J_old_left[a][b] * jump_nodal;
                    }
                }

                if(assemble_jacobian){
                    /*
                     * Frozen-J_old linearization of the jump operator.  This
                     * is consistent with the modified-Newton treatment already
                     * used by the verified SUPS tangent for stabilization
                     * coefficients.
                     */
                    for(alpha = 0; alpha < nt; alpha++){
                        for(beta = 0; beta < nt; beta++){
                            for(a = 0; a < 4; a++){
                                for(b = 0; b < 3; b++){
                                    FLUID_JAC_IP(
                                        alpha, beta, a, b, p) +=
                                        -layout->time.e_minus[alpha]
                                        * layout->time.e_minus[beta]
                                        * J_old_left[a][b];

                                    if(slab > 0){
                                        FLUID_PREV_JAC_IP(
                                            alpha, beta, a, b, p) +=
                                            layout->time.e_minus[alpha]
                                            * layout->time.e_plus[beta]
                                            * J_old_left[a][b];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(assemble_jacobian){
                /*
                 * Interior temporal Jacobian.
                 *
                 * Existing current-state tangent direction:
                 *   (delta U, delta U_old, delta du) = (1, 0, 1)
                 *
                 * Existing old-state tangent direction:
                 *   (0, 1, -1)
                 *
                 * A temporal coefficient beta generates:
                 *   delta U     = psi_beta
                 *   delta du    = dpsi_beta/dtau
                 *   delta U_old = psi_beta - dpsi_beta/dtau
                 *
                 * so the exact directional combination of the verified
                 * tangents is
                 *
                 *   psi_beta J_current
                 *   + (psi_beta - dpsi_beta) J_old.
                 */
                for(q = 0; q < nq; q++){
                    const double tau_time =
                        layout->time.tau[q];

                    for(p = 0; p < np; p++){
                        const size_t qp =
                            (size_t)q * (size_t)np
                            + (size_t)p;
                        double old_effective[3] = {
                            fields.v[qp][0]
                                - fields.du_time[qp][0],
                            fields.v[qp][1]
                                - fields.du_time[qp][1],
                            fields.v[qp][2]
                                - fields.du_time[qp][2]
                        };
                        double grad_storage[3][3];
                        double* grad_rows[3] = {
                            grad_storage[0],
                            grad_storage[1],
                            grad_storage[2]
                        };
                        double J_current[4][4] = {{0.0}};
                        double J_old[4][3];

                        for(a = 0; a < 3; a++){
                            for(b = 0; b < 3; b++){
                                grad_storage[a][b] =
                                    fields.grad_v[qp][a][b];
                            }
                        }

                        BBFE_elemmat_fluid_sups_mat_NR(
                            J_current,
                            fields.J_inv[p],
                            system->basis.N[p][i],
                            system->basis.N[p][j],
                            system->fe.geo[elem][p].grad_N[i],
                            system->fe.geo[elem][p].grad_N[j],
                            fields.v[qp],
                            grad_rows,
                            fields.grad_pressure[qp],
                            system->vals.density,
                            system->vals.viscosity,
                            fields.tau[qp],
                            fields.tau_c[qp],
                            system->vals.dt,
                            fields.du_time[qp]);

                        fluid_dg_old_jacobian_ip(
                            &system->fe,
                            &system->basis,
                            &system->vals,
                            elem,
                            i,
                            j,
                            p,
                            fields.v[qp],
                            old_effective,
                            grad_rows,
                            fields.pressure[qp],
                            fields.grad_pressure[qp],
                            fields.tau[qp],
                            fields.tau_c[qp],
                            fields.du_time[qp],
                            system->vals.st_dg_jacobian_fd_epsilon,
                            J_old);

                        for(alpha = 0; alpha < nt; alpha++){
                            const double test_weight =
                                layout->time.weight[q]
                                * layout->time.psi[alpha][q];

                            for(beta = 0; beta < nt; beta++){
                                const double psi_beta =
                                    layout->time.psi[beta][q];
                                const double dpsi_beta =
                                    fluid_time_basis_derivative(
                                        &layout->time,
                                        beta,
                                        tau_time);

                                for(a = 0; a < 4; a++){
                                    for(b = 0; b < 4; b++){
                                        double value =
                                            psi_beta
                                            * J_current[a][b];

                                        if(b < 3){
                                            value +=
                                                (psi_beta - dpsi_beta)
                                                * J_old[a][b];
                                        }

                                        FLUID_JAC_IP(
                                            alpha,
                                            beta,
                                            a,
                                            b,
                                            p) +=
                                            test_weight * value;
                                    }
                                }
                            }
                        }
                    }
                }

                for(alpha = 0; alpha < nt; alpha++){
                    for(beta = 0; beta < nt; beta++){
                        for(a = 0; a < 4; a++){
                            for(b = 0; b < 4; b++){
                                const double value =
                                    BBFE_std_integ_calc(
                                        np,
                                        &FLUID_JAC_IP(
                                            alpha,
                                            beta,
                                            a,
                                            b,
                                            0),
                                        system->basis.integ_weight,
                                        fields.Jacobian_ip);
                                const int row_dof =
                                    ST_fom_block_dof(
                                        layout,
                                        slab,
                                        alpha,
                                        a);
                                const int col_dof =
                                    ST_fom_block_dof(
                                        layout,
                                        slab,
                                        beta,
                                        b);

                                monolis_add_scalar_to_sparse_matrix_R(
                                    &system->monolis_window,
                                    system->fe.conn[elem][i],
                                    system->fe.conn[elem][j],
                                    row_dof,
                                    col_dof,
                                    value);
                            }

                            if(slab > 0){
                                for(b = 0; b < 3; b++){
                                    const double value =
                                        BBFE_std_integ_calc(
                                            np,
                                            &FLUID_PREV_JAC_IP(
                                                alpha,
                                                beta,
                                                a,
                                                b,
                                                0),
                                            system->basis.integ_weight,
                                            fields.Jacobian_ip);
                                    const int row_dof =
                                        ST_fom_block_dof(
                                            layout,
                                            slab,
                                            alpha,
                                            a);
                                    const int col_dof =
                                        ST_fom_block_dof(
                                            layout,
                                            slab - 1,
                                            beta,
                                            b);

                                    monolis_add_scalar_to_sparse_matrix_R(
                                        &system->monolis_window,
                                        system->fe.conn[elem][i],
                                        system->fe.conn[elem][j],
                                        row_dof,
                                        col_dof,
                                        value);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(alpha = 0; alpha < nt; alpha++){
            for(a = 0; a < 4; a++){
                double value;

                for(p = 0; p < np; p++){
                    FLUID_RES_IP(alpha, a, p) +=
                        layout->time.e_minus[alpha]
                        * FLUID_JUMP_IP(a, p);
                }

                value = BBFE_std_integ_calc(
                    np,
                    &FLUID_RES_IP(alpha, a, 0),
                    system->basis.integ_weight,
                    fields.Jacobian_ip);

                system->monolis_window.mat.R.B[
                    system->window_dof
                        * system->fe.conn[elem][i]
                    + ST_fom_block_dof(
                        layout,
                        slab,
                        alpha,
                        a)] -= value;
            }
        }
    }

#undef FLUID_RES_IP
#undef FLUID_JUMP_IP
#undef FLUID_JAC_IP
#undef FLUID_PREV_JAC_IP

    status = 0;

cleanup:
    fluid_dg_element_fields_free(&fields);
    free(coefficient);
    free(previous_trace_local);
    free(residual_ip);
    free(jump_ip);
    free(jac_ip);
    free(prev_jac_ip);
    return status;
}

static int fluid_dg_assemble_window(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state,
    int assemble_jacobian)
{
    int elem;
    int slab;

    if(system == NULL || layout == NULL ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    monolis_clear_mat_value_R(&system->monolis_window);

    for(elem = 0; elem < system->fe.total_num_elems; elem++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            if(fluid_dg_add_element_slab(
                system,
                layout,
                previous_trace,
                window_state,
                elem,
                slab,
                assemble_jacobian) != 0)
            {
                return 1;
            }
        }
    }

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &system->monolis_window,
        system->fe.total_num_nodes,
        system->window_dof,
        &system->bc_NR_window,
        system->monolis_window.mat.R.B);

    return 0;
}

static void fluid_dg_constant_predictor(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state)
{
    int node;
    int slab;
    int alpha;
    int component;

    memset(
        window_state,
        0,
        ST_fom_window_vector_length(layout) * sizeof(double));

    for(node = 0; node < layout->num_space_points; node++){
        for(slab = 0; slab < layout->num_slabs; slab++){
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                for(component = 0;
                    component < FLUID_SUPS_PHYSICAL_DOF;
                    component++)
                {
                    window_state[
                        ST_fom_vector_index(
                            layout,
                            node,
                            slab,
                            alpha,
                            component)] =
                        previous_trace[
                            (size_t)node
                                * FLUID_SUPS_PHYSICAL_DOF
                            + (size_t)component];
                }
            }
        }
    }

    fluid_apply_physical_bc_to_window(
        system, layout, window_state);
    fluid_apply_pressure_gauge_to_window_state(
        system, layout, window_state);
}

static int fluid_sups_st_solve_window_all_at_once_high_order(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report)
{
    const int slabs = layout != NULL ? layout->num_slabs : 0;
    const int nstate = system != NULL
        ? system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF
        : 0;
    const size_t window_length =
        layout != NULL ? ST_fom_window_vector_length(layout) : 0;
    double* initial_slab_total = NULL;
    double* initial_slab_momentum = NULL;
    double* initial_slab_mass = NULL;
    double* final_slab_total = NULL;
    double* final_slab_momentum = NULL;
    double* final_slab_mass = NULL;
    double* final_trace = NULL;
    double initial_total = 0.0;
    double initial_momentum = 0.0;
    double initial_mass = 0.0;
    double final_total = INFINITY;
    double final_momentum = INFINITY;
    double final_mass = INFINITY;
    double relative_total = INFINITY;
    double relative_momentum = INFINITY;
    double relative_mass = INFINITY;
    const int linear_solver =
        system != NULL ? fluid_window_linear_solver_type()
                       : MONOLIS_ITER_BICGSAFE;
    const double linear_epsilon =
        system != NULL ? fluid_window_linear_epsilon(system)
                       : 1.0e-8;
    const int linear_max_iter =
        system != NULL ? fluid_window_linear_max_iter(system)
                       : 10000;
    int iteration;
    int slab;
    int converged = 0;
    int status = 1;

    if(system == NULL || layout == NULL ||
       previous_trace == NULL || window_state == NULL ||
       slabs <= 0 || nstate <= 0 ||
       layout->time.degree < 1 ||
       layout->time.degree > 2 ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       layout->num_owned_space_points !=
           system->num_owned_space_nodes ||
       system->window_dof !=
           ST_fom_window_dof_per_space_point(layout))
    {
        return 1;
    }

    initial_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    initial_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_total =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_momentum =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_slab_mass =
        (double*)calloc((size_t)slabs, sizeof(double));
    final_trace =
        (double*)calloc((size_t)nstate, sizeof(double));

    if(initial_slab_total == NULL ||
       initial_slab_momentum == NULL ||
       initial_slab_mass == NULL ||
       final_slab_total == NULL ||
       final_slab_momentum == NULL ||
       final_slab_mass == NULL ||
       final_trace == NULL)
    {
        goto cleanup;
    }

    if(window_report != NULL){
        memset(window_report, 0, sizeof(*window_report));
        window_report->relative_total_residual = INFINITY;
        window_report->relative_momentum_residual = INFINITY;
        window_report->relative_mass_residual = INFINITY;
    }
    if(slab_reports != NULL){
        memset(
            slab_reports,
            0,
            (size_t)slabs * sizeof(*slab_reports));
    }

    fluid_dg_constant_predictor(
        system,
        layout,
        previous_trace,
        window_state);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        window_state);

    for(iteration = 0;
        iteration < system->vals.nr_max_iter;
        iteration++)
    {
        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "\n%s ================ ALL-AT-ONCE dG(%d) window %d "
                "Newton %d ================\n",
                CODENAME,
                layout->time.degree,
                window_id,
                iteration + 1);
        }

        if(fluid_dg_assemble_window(
            system,
            layout,
            previous_trace,
            window_state,
            1) != 0)
        {
            goto cleanup;
        }

        if(iteration == 0){
            fluid_window_rhs_norms(
                system,
                layout,
                &initial_total,
                &initial_momentum,
                &initial_mass,
                initial_slab_total,
                initial_slab_momentum,
                initial_slab_mass);
        }

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s dG(%d) window Newton linear solve: "
                "%d dof/node = %d physical x %d time dof x %d slabs; "
                "ONE MONOLIS solve; solver=%s eps=%.3e max_iter=%d\n",
                CODENAME,
                layout->time.degree,
                system->window_dof,
                FLUID_SUPS_PHYSICAL_DOF,
                layout->time.n_dof,
                slabs,
                fluid_window_linear_solver_name(linear_solver),
                linear_epsilon,
                linear_max_iter);
        }

        /*
         * Use a zero correction initial guess every Newton iteration.
         * Reusing the previous Newton correction is especially harmful when
         * the nonlinear residual has already dropped by several orders.
         */
        fluid_zero_window_linear_solution(system);

        BBFE_sys_monowrap_solve(
            &system->monolis_window,
            &system->mono_com,
            system->monolis_window.mat.R.X,
            linear_solver,
            MONOLIS_PREC_DIAG,
            linear_max_iter,
            linear_epsilon);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            system->monolis_window.mat.R.X);

        {
            size_t row;
            for(row = 0; row < window_length; row++){
                window_state[row] +=
                    system->monolis_window.mat.R.X[row];
            }
        }

        fluid_apply_physical_bc_to_window(
            system, layout, window_state);
        fluid_apply_pressure_gauge_to_window_state(
            system, layout, window_state);

        monolis_mpi_update_R(
            &system->mono_com,
            system->fe.total_num_nodes,
            system->window_dof,
            window_state);

        if(fluid_dg_assemble_window(
            system,
            layout,
            previous_trace,
            window_state,
            0) != 0)
        {
            goto cleanup;
        }

        fluid_window_rhs_norms(
            system,
            layout,
            &final_total,
            &final_momentum,
            &final_mass,
            final_slab_total,
            final_slab_momentum,
            final_slab_mass);

        relative_total =
            fluid_window_relative_norm(
                final_total, initial_total);
        relative_momentum =
            fluid_window_relative_norm(
                final_momentum, initial_momentum);
        relative_mass =
            fluid_window_relative_norm(
                final_mass, initial_mass);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "%s dG(%d) window %d Newton %d residual: "
                "total=%e momentum=%e mass=%e "
                "(abs total=%.15e initial=%.15e)\n",
                CODENAME,
                layout->time.degree,
                window_id,
                iteration + 1,
                relative_total,
                relative_momentum,
                relative_mass,
                final_total,
                initial_total);
        }

        if(isfinite(relative_total) &&
           relative_total < system->vals.nr_epsilon)
        {
            converged = 1;
            break;
        }
    }

    if(!converged){
        fprintf(
            stderr,
            "ERROR: dG(%d) all-at-once window Newton did not converge: "
            "window=%d rel_total=%.15e\n",
            layout->time.degree,
            window_id,
            relative_total);
        goto cleanup;
    }

    if(slab_reports != NULL){
        for(slab = 0; slab < slabs; slab++){
            slab_reports[slab].converged = 1;
            slab_reports[slab].nonlinear_iterations =
                iteration + 1;
            slab_reports[slab].relative_total_residual =
                fluid_window_relative_norm(
                    final_slab_total[slab],
                    initial_slab_total[slab]);
            slab_reports[slab].relative_momentum_residual =
                fluid_window_relative_norm(
                    final_slab_momentum[slab],
                    initial_slab_momentum[slab]);
            slab_reports[slab].relative_mass_residual =
                fluid_window_relative_norm(
                    final_slab_mass[slab],
                    initial_slab_mass[slab]);
        }
    }

    if(window_report != NULL){
        window_report->converged = 1;
        window_report->nonlinear_iterations =
            iteration + 1;
        window_report->relative_total_residual =
            relative_total;
        window_report->relative_momentum_residual =
            relative_momentum;
        window_report->relative_mass_residual =
            relative_mass;
    }

    if(ST_fom_extract_last_right_trace(
        layout,
        window_state,
        final_trace) != ST_FOM_SUCCESS)
    {
        goto cleanup;
    }

    if(fluid_sups_st_import_state(
        system,
        final_trace,
        1) != 0 ||
       fluid_sups_st_synchronize_state(system) != 0)
    {
        goto cleanup;
    }

    status = 0;

cleanup:
    free(initial_slab_total);
    free(initial_slab_momentum);
    free(initial_slab_mass);
    free(final_slab_total);
    free(final_slab_momentum);
    free(final_slab_mass);
    free(final_trace);
    return status;
}

int fluid_sups_st_solve_window_all_at_once_dg(
    FLUID_SUPS_SYSTEM* system,
    int window_id,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* window_state,
    FLUID_SUPS_NR_REPORT* slab_reports,
    FLUID_SUPS_NR_REPORT* window_report)
{
    if(layout == NULL){
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_sups_st_solve_window_all_at_once_dg0(
            system,
            window_id,
            layout,
            previous_trace,
            window_state,
            slab_reports,
            window_report);
    }

    return fluid_sups_st_solve_window_all_at_once_high_order(
        system,
        window_id,
        layout,
        previous_trace,
        window_state,
        slab_reports,
        window_report);
}


/* ========================================================================== */
/* Public nonlinear FOM adapter used by the NS ST-DDROM                       */
/* ========================================================================== */

static int fluid_validate_rom_layout(
    const FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    if(system == NULL || layout == NULL ||
       layout->num_space_points != system->fe.total_num_nodes ||
       layout->num_owned_space_points != system->num_owned_space_nodes ||
       layout->physical_dof != FLUID_SUPS_PHYSICAL_DOF ||
       ST_fom_window_dof_per_space_point(layout) != system->window_dof)
    {
        return 1;
    }

    return 0;
}

int fluid_sups_st_rom_build_window_lift(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* lift)
{
    const size_t length =
        layout != NULL ? ST_fom_window_vector_length(layout) : 0u;

    if(fluid_validate_rom_layout(system, layout) != 0 ||
       lift == NULL || length == 0u)
    {
        return 1;
    }

    memset(lift, 0, length * sizeof(double));

    fluid_apply_physical_bc_to_window(
        system,
        layout,
        lift);
    fluid_apply_pressure_gauge_to_window_state(
        system,
        layout,
        lift);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        lift);

    return 0;
}

int fluid_sups_st_rom_build_window_predictor(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    double* predictor)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || predictor == NULL)
    {
        return 1;
    }

    fluid_dg_constant_predictor(
        system,
        layout,
        previous_trace,
        predictor);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        predictor);

    return 0;
}

int fluid_sups_st_rom_enforce_window_constraints(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       window_state == NULL)
    {
        return 1;
    }

    fluid_apply_physical_bc_to_window(
        system,
        layout,
        window_state);
    fluid_apply_pressure_gauge_to_window_state(
        system,
        layout,
        window_state);

    monolis_mpi_update_R(
        &system->mono_com,
        system->fe.total_num_nodes,
        system->window_dof,
        window_state);

    return 0;
}

int fluid_sups_st_rom_assemble_window_newton(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_assemble_window_newton_system(
            system,
            layout,
            previous_trace,
            window_state);
    }

    return fluid_dg_assemble_window(
        system,
        layout,
        previous_trace,
        window_state,
        1);
}

int fluid_sups_st_rom_assemble_window_residual(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* previous_trace,
    const double* window_state)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       previous_trace == NULL || window_state == NULL)
    {
        return 1;
    }

    if(layout->time.degree == 0){
        return fluid_assemble_window_residual_only(
            system,
            layout,
            previous_trace,
            window_state);
    }

    return fluid_dg_assemble_window(
        system,
        layout,
        previous_trace,
        window_state,
        0);
}

int fluid_sups_st_rom_get_window_rhs_norms(
    FLUID_SUPS_SYSTEM* system,
    const ST_FOM_WINDOW_LAYOUT* layout,
    double* total_norm,
    double* momentum_norm,
    double* mass_norm)
{
    if(fluid_validate_rom_layout(system, layout) != 0 ||
       total_norm == NULL || momentum_norm == NULL || mass_norm == NULL)
    {
        return 1;
    }

    fluid_window_rhs_norms(
        system,
        layout,
        total_norm,
        momentum_norm,
        mass_norm,
        NULL,
        NULL,
        NULL);

    return 0;
}


static void fluid_evaluate_NR_residual(
    FLUID_SUPS_SYSTEM* sys,
    const double* rvec_old,
    const double* rvec_m_old,
    const double* rvec_c_old,
    double* rel_r,
    double* rel_r_m,
    double* rel_r_c)
{

    int num_internal_vtx;

    if(monolis_mpi_get_global_comm_size()==1){
        num_internal_vtx = sys->fe.total_num_nodes;
    }
    else{
        num_internal_vtx = sys->mono_com.n_internal_vertex;
    }

    /** 更新後の全体残差を評価 **/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_NR(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    double norm_r = BBFE_fluid_calc_internal_norm_total(sys->monolis.mat.R.B, num_internal_vtx);

    double norm_r_old = BBFE_fluid_calc_internal_norm_total(rvec_old, num_internal_vtx);

    /** 更新後の運動方程式残差を評価 **/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_momentum_eq(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    double norm_r_m = BBFE_fluid_calc_internal_norm_momentum(sys->monolis.mat.R.B, num_internal_vtx);

    double norm_r_m_old = BBFE_fluid_calc_internal_norm_momentum(rvec_m_old, num_internal_vtx);

    /** 更新後の連続の式残差を評価 **/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_mass_con(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    double norm_r_c = BBFE_fluid_calc_internal_norm_mass(sys->monolis.mat.R.B, num_internal_vtx);

    double norm_r_c_old = BBFE_fluid_calc_internal_norm_mass(rvec_c_old, num_internal_vtx);

    /** MPI 領域間で総和 **/
    if(monolis_mpi_get_global_comm_size()==1){
    }
    else{
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys->mono_com.comm);

        monolis_allreduce_R(1, &norm_r, MONOLIS_MPI_SUM, sys->mono_com.comm);

        monolis_allreduce_R(1, &norm_r_m_old, MONOLIS_MPI_SUM, sys->mono_com.comm);

        monolis_allreduce_R(1, &norm_r_m, MONOLIS_MPI_SUM, sys->mono_com.comm);

        monolis_allreduce_R(1, &norm_r_c_old, MONOLIS_MPI_SUM, sys->mono_com.comm);

        monolis_allreduce_R( 1, &norm_r_c, MONOLIS_MPI_SUM, sys->mono_com.comm);
    }

    /** ノルムと相対残差を計算 **/
    const double nrm_r_old = sqrt(norm_r_old);
    const double nrm_r     = sqrt(norm_r);

    const double nrm_r_m_old = sqrt(norm_r_m_old);
    const double nrm_r_m     = sqrt(norm_r_m);

    const double nrm_r_c_old = sqrt(norm_r_c_old);
    const double nrm_r_c     = sqrt(norm_r_c);

    *rel_r = nrm_r_old > 1.0e-300 ? nrm_r / nrm_r_old : nrm_r;
    *rel_r_m = nrm_r_m_old > 1.0e-300
        ? nrm_r_m / nrm_r_m_old : nrm_r_m;
    *rel_r_c = nrm_r_c_old > 1.0e-300
        ? nrm_r_c / nrm_r_c_old : nrm_r_c;

    if(monolis_mpi_get_global_my_rank()==0){
        printf(" residual total    : nrm = %e, old = %e, rel = %e\n", nrm_r, nrm_r_old, *rel_r);
        printf(" residual momentum : nrm = %e, old = %e, rel = %e\n", nrm_r_m, nrm_r_m_old, *rel_r_m);
        printf(" residual mass     : nrm = %e, old = %e, rel = %e\n", nrm_r_c, nrm_r_c_old, *rel_r_c);
    }
    
}

static void fluid_store_initial_NR_residual(
    FLUID_SUPS_SYSTEM* sys,
    double* rvec_old,
    double* rvec_m_old,
    double* rvec_c_old)
{
    int num_dof;
    
    if(monolis_mpi_get_global_comm_size()==1){
        num_dof = sys->fe.total_num_nodes * 4;
    }
    else{
        num_dof = sys->mono_com.n_internal_vertex * 4;
    }

    /** 初回 Newton 反復時の運動方程式残差を保存 **/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_momentum_eq(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    BBFE_fluid_vec_copy_1d(rvec_m_old, sys->monolis.mat.R.B, num_dof);

    /** 初回 Newton 反復時の連続の式残差を保存**/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_mass_con(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    BBFE_fluid_vec_copy_1d(
        rvec_c_old,
        sys->monolis.mat.R.B,
        num_dof);

    /** 初回 Newton 反復時の全体残差を保存 **/
    monolis_clear_mat_value_R(&(sys->monolis));

    fluid_set_element_residuals_NR(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    BBFE_fluid_vec_copy_1d(rvec_old, sys->monolis.mat.R.B, num_dof);
}

int fluid_sups_st_solve_step(
    FLUID_SUPS_SYSTEM* system,
    int step,
    double time_old,
    double time_new,
    FLUID_SUPS_NR_REPORT* report)
{
    const int num_dof = system != NULL
        ? system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF : 0;
    double* residual_initial = NULL;
    double* momentum_initial = NULL;
    double* mass_initial = NULL;
    int iteration;
    int converged = 0;
    double rel_total = INFINITY;
    double rel_momentum = INFINITY;
    double rel_mass = INFINITY;

    (void)time_old;
    (void)time_new;

    if(system == NULL || num_dof <= 0){
        return 1;
    }

    if(report != NULL){
        memset(report, 0, sizeof(*report));
        report->relative_total_residual = INFINITY;
        report->relative_momentum_residual = INFINITY;
        report->relative_mass_residual = INFINITY;
    }

    residual_initial = (double*)calloc((size_t)num_dof, sizeof(double));
    momentum_initial = (double*)calloc((size_t)num_dof, sizeof(double));
    mass_initial = (double*)calloc((size_t)num_dof, sizeof(double));
    if(residual_initial == NULL || momentum_initial == NULL ||
       mass_initial == NULL)
    {
        free(residual_initial);
        free(momentum_initial);
        free(mass_initial);
        return 1;
    }

    for(iteration = 0; iteration < system->vals.nr_max_iter; iteration++){
        printf(
            "\n%s ----------------- ST dG(0) step %d - NR step %d "
            "----------------\n",
            CODENAME, step, iteration);

        if(iteration == 0){
            fluid_store_initial_NR_residual(
                system,
                residual_initial,
                momentum_initial,
                mass_initial);
        }

        /* Preserve the verified Newton assembly and solve path. */
        monolis_clear_mat_value_R(&system->monolis);
        fluid_set_element_jacobian_NR(
            &system->monolis,
            &system->fe,
            &system->basis,
            &system->vals);
        fluid_set_element_residuals_NR(
            &system->monolis,
            &system->fe,
            &system->basis,
            &system->vals);

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &system->monolis,
            system->fe.total_num_nodes,
            FLUID_SUPS_PHYSICAL_DOF,
            &system->bc_NR,
            system->monolis.mat.R.B);

        BBFE_sys_monowrap_solve(
            &system->monolis,
            &system->mono_com,
            system->monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            system->vals.mat_max_iter,
            system->vals.mat_epsilon);

        BBFE_fluid_sups_renew_velocity(
            system->vals.delta_v,
            system->monolis.mat.R.X,
            system->fe.total_num_nodes);
        BBFE_fluid_sups_renew_pressure(
            system->vals.delta_p,
            system->monolis.mat.R.X,
            system->fe.total_num_nodes);
        BBFE_fluid_update_velocity_pressure_NR(
            system->vals.v,
            system->vals.delta_v,
            system->vals.p,
            system->vals.delta_p,
            system->fe.total_num_nodes);

        fluid_evaluate_NR_residual(
            system,
            residual_initial,
            momentum_initial,
            mass_initial,
            &rel_total,
            &rel_momentum,
            &rel_mass);

        if(isfinite(rel_total) && rel_total < system->vals.nr_epsilon){
            converged = 1;
            fluid_copy_current_to_previous(system);
            break;
        }
    }

    /* Make the converged state a valid input on every local halo. */
    if(converged && fluid_sups_st_synchronize_state(system) != 0){
        converged = 0;
    }
    if(converged){
        fluid_copy_current_to_previous(system);
    }

    if(report != NULL){
        report->converged = converged;
        report->nonlinear_iterations = converged
            ? iteration + 1 : system->vals.nr_max_iter;
        report->relative_total_residual = rel_total;
        report->relative_momentum_residual = rel_momentum;
        report->relative_mass_residual = rel_mass;
    }

    free(residual_initial);
    free(momentum_initial);
    free(mass_initial);

    if(!converged){
        fprintf(stderr,
            "ERROR: fluid Newton solve did not converge at step %d; "
            "last relative residual = %.15e.\n",
            step, rel_total);
        return 1;
    }
    return 0;
}

int fluid_sups_st_initialize(
    FLUID_SUPS_SYSTEM* system,
    int argc,
    char* argv[])
{
    const char* filename;
    double* state = NULL;
    const char* initial_condition;

    if(system == NULL){
        return 1;
    }
    memset(system, 0, sizeof(*system));

    system->cond.directory = BBFE_fluid_get_directory_name(
        argc, argv, CODENAME);
    fluid_read_calc_conditions(&system->vals, system->cond.directory);

    BBFE_fluid_pre(
        &system->fe,
        &system->basis,
        argc,
        argv,
        system->cond.directory,
        system->vals.num_ip_each_axis);

    filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_INPUT_FILENAME_D_BC_V);
    fluid_read_dirichlet_bc(
        &system->bc,
        filename,
        system->cond.directory,
        system->fe.total_num_nodes,
        0);
    fluid_read_dirichlet_bc(
        &system->bc_NR,
        filename,
        system->cond.directory,
        system->fe.total_num_nodes,
        1);

    if(fluid_allocate_values(
        &system->vals, system->fe.total_num_nodes) != 0)
    {
        return 1;
    }

    BBFE_elemmat_set_Jacobi_mat(&system->fe, &system->basis);
    BBFE_elemmat_set_shapefunc_derivative(&system->fe, &system->basis);

    BBFE_sys_monowrap_init_monomat(
        &system->monolis,
        &system->mono_com,
        &system->fe,
        FLUID_SUPS_PHYSICAL_DOF,
        system->cond.directory);

    /*
     * Allocate the true all-at-once window matrix once:
     *
     *   window_dof = 4 * slabs/window * temporal_dofs/slab.
     *
     * With four slabs this gives 16, 32, 48 dof/node for dG(0), dG(1),
     * dG(2), respectively.
     */
    if(fluid_resolve_owned_space_nodes(system) != 0){
        fprintf(stderr, "ERROR: failed to resolve fluid spatial ownership.\n");
        return 1;
    }

    if(fluid_initialize_window_bc(system) != 0){
        fprintf(stderr, "ERROR: failed to initialize fluid ST window BC/gauge.\n");
        return 1;
    }

    if(fluid_initialize_window_parallel_pattern(system) != 0){
        fprintf(stderr,
            "ERROR: failed to initialize distributed fluid ST window pattern.\n");
        return 1;
    }

    initial_condition = getenv("FLUID_ST_INITIAL_CONDITION");
    if(initial_condition != NULL && strcmp(initial_condition, "cavity") == 0){
        BBFE_fluid_initialize_velocity_pressure(
            system->vals.v,
            system->vals.p,
            system->fe.total_num_nodes);
    }
    else{
        BBFE_fluid_initialize_velocity_pressure_karman_vortex(
            system->vals.v,
            system->vals.v_old,
            system->vals.p,
            system->fe.total_num_nodes);
    }

    state = (double*)calloc(
        (size_t)system->fe.total_num_nodes * FLUID_SUPS_PHYSICAL_DOF,
        sizeof(double));
    if(state == NULL){
        return 1;
    }

    fluid_pack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, state);
    BBFE_fluid_add_Dbc(
        state,
        &system->bc,
        system->fe.total_num_nodes,
        FLUID_SUPS_PHYSICAL_DOF);
    fluid_unpack_velocity_pressure(
        &system->vals, system->fe.total_num_nodes, state);
    free(state);

    if(fluid_sups_st_synchronize_state(system) != 0){
        return 1;
    }
    fluid_copy_current_to_previous(system);
    fluid_sups_st_print_conditions(system);
    return 0;
}

void fluid_sups_st_finalize(FLUID_SUPS_SYSTEM* system)
{
    int nnode;

    if(system == NULL){
        return;
    }
    nnode = system->fe.total_num_nodes;

    BBFE_sys_memory_free_Dirichlet_bc(
        &system->bc, nnode, FLUID_SUPS_PHYSICAL_DOF);
    BBFE_sys_memory_free_Dirichlet_bc(
        &system->bc_NR, nnode, FLUID_SUPS_PHYSICAL_DOF);
    BBFE_sys_memory_free_Dirichlet_bc(
        &system->bc_NR_window, nnode, system->window_dof);
    monolis_finalize(&system->monolis);
    monolis_finalize(&system->monolis_window);
    ST_fom_spatial_graph_finalize(&system->st_spatial_graph);
    fluid_free_values(&system->vals, nnode);
    BBFE_fluid_finalize(&system->fe, &system->basis);
    memset(system, 0, sizeof(*system));
}
