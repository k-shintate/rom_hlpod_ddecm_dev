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
