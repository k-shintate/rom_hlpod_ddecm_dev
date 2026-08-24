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

