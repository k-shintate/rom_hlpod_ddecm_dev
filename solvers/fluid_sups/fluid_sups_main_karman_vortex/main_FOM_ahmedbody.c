#include "core_ROM.h"
//#include "core_Nitsche.h"

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

const int BUFFER_SIZE = 1024;

double hot_start_read_initialize_val(
    double*     int_val,
    const char* input_fname,
    const char* directory)
{
    int BUFFER_SIZE = 1024;
    int total_num_nodes;
    int ndof;
    double t = 0.0;
        FILE* fp;
        char fname[BUFFER_SIZE];
        char id[BUFFER_SIZE];

        fp = BBFE_sys_read_fopen(fp, "hot_start/start_time.dat", directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%lf", &(t));
    fclose(fp);

        fp = BBFE_sys_read_fopen(fp, input_fname, directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes), &(ndof));
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fscanf(fp, "%lf", &(int_val[i * ndof + j]));
        }
    }
        fclose(fp);

    return t;
}

void hot_start_write_initialize_val(
    double*         int_val,
    const int       total_num_nodes,
    const int       ndof,
    const double    time,
    const char*     output_fname,
    const char*     directory)
{
        FILE* fp;
        char fname[BUFFER_SIZE];
        char id[BUFFER_SIZE];

        fp = BBFE_sys_write_fopen(fp, output_fname, directory);
        fprintf(fp, "initialization\n");
    fprintf(fp, "%d %d\n", total_num_nodes, ndof);
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fprintf(fp, "%e ", int_val[i * ndof + j]);
        }
        fprintf(fp, "\n");
    }
        fclose(fp);

    if(monolis_mpi_get_global_my_rank()==0){
        fp = BBFE_sys_write_fopen(fp, "hot_start/start_time.dat", directory);
        fprintf(fp, "start_time\n");
        fprintf(fp, "%lf", time);
        fclose(fp);
    }
}

static int t4_diag_trace_is_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;

    if (!initialized) {
        const char* env =
            getenv("T4_DIAG_TRACE");

        enabled =
            (
                env != NULL &&
                env[0] != '\0' &&
                strcmp(env, "0") != 0
            );

        initialized = 1;
    }

    return enabled;
}

static int t4_all_ranks_status_ok(
    int local_status,
    MONOLIS_COM* monolis_com,
    const char* label)
{
    if (monolis_com == NULL) {
        fprintf(
            stderr,
            "[rank %d][%s] monolis_com is NULL\n",
            monolis_mpi_get_global_my_rank(),
            label != NULL ? label : "diagnostic");

        return 0;
    }

    /*
     * local_status == 0 : success
     * local_status != 0 : failure
     */
    double failure_count =
        (local_status == 0)
        ? 0.0
        : 1.0;

    /*
     * 成功・失敗にかかわらず全rankが呼ぶ。
     */
    monolis_allreduce_R(
        1,
        &failure_count,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    /*
     * failure_countは全rankで同じ値になる。
     * したがって、この後の分岐も全rankで一致する。
     */
    if (failure_count > 0.5) {
        if (local_status != 0) {
            fprintf(
                stderr,
                "[rank %d][%s] local calculation failed: status=%d\n",
                monolis_mpi_get_global_my_rank(),
                label != NULL ? label : "diagnostic",
                local_status);
        }

        if (monolis_mpi_get_global_my_rank() == 0) {
            fprintf(
                stderr,
                "[%s] skipped because %.0f rank(s) failed\n",
                label != NULL ? label : "diagnostic",
                failure_count);
        }

        return 0;
    }

    return 1;
}

static void t4_diag_trace(
    const char* label,
    int step)
{
    if (!t4_diag_trace_is_enabled()) {
        return;
    }

    fprintf(
        stderr,
        "[rank %d] step=%d %s\n",
        monolis_mpi_get_global_my_rank(),
        step,
        label != NULL ? label : "(null)");

    fflush(stderr);
}

static int run_wall_and_continuity_diagnostics_parallel_safe(
    FE_SYSTEM* sys,
    int step,
    double t,
    double U_ref,
    const double Uw[3],
    const double eU[3],
    const double eP[3],
    double gamma_n,
    int write_nodal_yplus)
{
    if (sys == NULL) {
        fprintf(
            stderr,
            "[rank %d][diagnostics] sys is NULL\n",
            monolis_mpi_get_global_my_rank());

        return -1;
    }

    const int rank =
        monolis_mpi_get_global_my_rank();

    int wall_mesh_global_ok = 0;
    int wall_velocity_global_ok = 0;
    int continuity_global_ok = 0;

    /*
     * 局所計算が失敗した場合にも、
     * 未初期化の値を使わないようにする。
     */
    T4WallMeshDiagnostics wall_mesh_diag;
    T4WallDiagnostics wall_velocity_diag;
    T4DomainContinuityDiagnostics continuity_diag;

    memset(
        &wall_mesh_diag,
        0,
        sizeof(wall_mesh_diag));

    memset(
        &wall_velocity_diag,
        0,
        sizeof(wall_velocity_diag));

    memset(
        &continuity_diag,
        0,
        sizeof(continuity_diag));

    /* ================================================================ */
    /* 1. Wall mesh / y+ diagnostics                                    */
    /* ================================================================ */

    t4_diag_trace(
        "before wall mesh local calculation",
        step);

    const int wall_mesh_status_local =
        calc_wall_mesh_diagnostics_tet4_tri3(
            &(sys->surf),
            &(sys->fe),
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            U_ref,
            Uw,
            eU,
            gamma_n,
            write_nodal_yplus,
            &wall_mesh_diag);

    t4_diag_trace(
        "after wall mesh local calculation",
        step);

    wall_mesh_global_ok =
        t4_all_ranks_status_ok(
            wall_mesh_status_local,
            &(sys->mono_com),
            "wall mesh diagnostics");

    if (wall_mesh_global_ok) {
        t4_diag_trace(
            "before wall mesh diagnostics allreduce",
            step);

        wall_mesh_diagnostics_allreduce(
            &wall_mesh_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after wall mesh diagnostics allreduce",
            step);
    }

    /* ================================================================ */
    /* 2. Wall velocity / pressure / weak-BC diagnostics                */
    /* ================================================================ */

    const double p_ref_diag = 0.0;

    t4_diag_trace(
        "before wall velocity local calculation",
        step);

    const int wall_velocity_status_local =
        calc_wall_diagnostics_tet4_tri3_with_Uw(
            &(sys->surf),
            &(sys->fe),
            &(sys->basis_surf),
            &(sys->vals),
            sys->vals.density,
            U_ref,
            p_ref_diag,
            eU,
            eP,
            Uw,
            &wall_velocity_diag);

    t4_diag_trace(
        "after wall velocity local calculation",
        step);

    wall_velocity_global_ok =
        t4_all_ranks_status_ok(
            wall_velocity_status_local,
            &(sys->mono_com),
            "wall velocity diagnostics");

    if (wall_velocity_global_ok) {
        t4_diag_trace(
            "before wall velocity diagnostics allreduce",
            step);

        wall_diagnostics_allreduce(
            &wall_velocity_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after wall velocity diagnostics allreduce",
            step);
    }

    /* ================================================================ */
    /* 3. Domain continuity diagnostics                                 */
    /* ================================================================ */

    t4_diag_trace(
        "before continuity local calculation",
        step);

    const int continuity_status_local =
        calc_domain_continuity_diagnostics_tet4(
            &(sys->fe),
            &(sys->vals),
            &continuity_diag);

    t4_diag_trace(
        "after continuity local calculation",
        step);

    continuity_global_ok =
        t4_all_ranks_status_ok(
            continuity_status_local,
            &(sys->mono_com),
            "domain continuity diagnostics");

    if (continuity_global_ok) {
        t4_diag_trace(
            "before continuity diagnostics allreduce",
            step);

        domain_continuity_diagnostics_allreduce(
            &continuity_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after continuity diagnostics allreduce",
            step);
    }

    if (wall_velocity_global_ok) {
        const int wall_location_status =
            output_wall_max_error_locations_tet4_tri3(
                &(sys->surf),
                &(sys->fe),
                &(sys->basis_surf),
                &(sys->vals),
                U_ref,
                Uw,
                &(sys->mono_com),
                t,
                sys->cond.directory);

        if (
            wall_location_status != 0 &&
            rank == 0
        ) {
            /*
             * 最大誤差位置は補助診断であり、壁面統計本体とは独立。
             * 位置診断が失敗しても、すでにallreduce済みの
             * wall_velocity_diagは出力する。
             */
            fprintf(
                stderr,
                "[wall max location] failed; "
                "wall summary diagnostics will still be written.\n");
        }
    }

    if (continuity_global_ok) {
        const int div_location_status =
            output_domain_max_div_u_location_tet4(
                &(sys->fe),
                &(sys->vals),
                &(sys->mono_com),
                t,
                sys->cond.directory);

        if (
            div_location_status != 0 &&
            rank == 0
        ) {
            /*
             * 最大発散位置も補助診断。
             * 位置出力が失敗しても、領域発散統計本体は出力する。
             */
            fprintf(
                stderr,
                "[div-u max location] failed; "
                "continuity summary diagnostics will still be written.\n");
        }
    }

    /* ================================================================ */
    /* 4. Text output: rank 0 only                                      */
    /* ================================================================ */

    if (rank == 0) {
        if (wall_mesh_global_ok) {
            output_wall_mesh_diagnostics_tet4_tri3(
                &wall_mesh_diag,
                t,
                sys->cond.directory);
        }

        if (wall_velocity_global_ok) {
            fprintf(
                stderr,
                "[wall summary] writing at step=%d t=%.15e directory=%s\n",
                step,
                t,
                sys->cond.directory != NULL
                    ? sys->cond.directory
                    : "(null)");

            output_wall_diagnostics_tet4_tri3(
                &wall_velocity_diag,
                t,
                sys->cond.directory);
        } else {
            fprintf(
                stderr,
                "[wall summary] NOT written: "
                "wall_velocity_global_ok=0 at step=%d t=%.15e\n",
                step,
                t);
        }

        if (continuity_global_ok) {
            output_domain_continuity_diagnostics_tet4(
                &continuity_diag,
                t,
                sys->cond.directory);
        }
    }

    /*
     * どれか1つでも失敗した場合は、mainへ失敗を返す。
     * allreduceの順序自体は全rankで一致している。
     */
    if (
        !wall_mesh_global_ok ||
        !wall_velocity_global_ok ||
        !continuity_global_ok) {

        return -1;
    }

    return 0;
}

/*
 * Update the call in main(): add eP between eU and gamma_n.
 *
 *     run_wall_and_continuity_diagnostics_parallel_safe(
 *         &sys,
 *         step,
 *         t,
 *         U_ref,
 *         Uw,
 *         eU,
 *         eP,
 *         gamma_n,
 *         1);
 */



int main (
                int argc,
                char* argv[])
{
        printf("\n");

        FE_SYSTEM sys;

        monolis_global_initialize();

    double t1 = monolis_get_time();
        double FOM_t1 = monolis_get_time();

        sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);

        read_calc_conditions(&(sys.vals), sys.cond.directory);

        BBFE_fluid_pre(
                        &(sys.fe), &(sys.basis),
                        argc, argv, sys.cond.directory,
                        sys.vals.num_ip_each_axis);

        const char* filename;

        memory_allocation_nodal_values(
                        &(sys.vals),
                        sys.fe.total_num_nodes);

        filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
        BBFE_fluid_sups_read_Dirichlet_bc(
                        &(sys.bc),
                        filename,
                        sys.cond.directory,
                        sys.fe.total_num_nodes,
                        4);

const int ground_bc_status =
    BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
        &(sys.bc),
        &(sys.fe),
        0.0,
        1.0e-10
    );

if (ground_bc_status != 0) {
    fprintf(
        stderr,
        "%s ERROR: failed to add ground no-slip condition.\n",
        CODENAME
    );

    exit(EXIT_FAILURE);
}

    BBFE_fluid_sups_read_Dirichlet_bc_NR(
            &(sys.bc_NR),
            filename,
            sys.cond.directory,
            sys.fe.total_num_nodes,
            4);

    const int ground_bc_NR_status =
    BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
        &(sys.bc_NR),
        &(sys.fe),
        0.0,
        1.0e-10
    );

    if (ground_bc_NR_status != 0) {
        fprintf(
            stderr,
            "%s ERROR: failed to add ground no-slip condition "
            "to Newton-increment BC.\n",
            CODENAME
        );

        exit(EXIT_FAILURE);
    }

    pre_surface(
            //&(sys.monolis_com_surf),
            &(sys.surf),
            &(sys.basis_surf),
            sys.cond.directory,
            "surf_graph.dat",
            sys.vals.num_ip_each_axis);
    
    pre_surface_internal(
        //&(sys.monolis_com_surf),
        &(sys.surf_internal),
        &(sys.basis_surf_internal),
        sys.cond.directory,
        "surf_graph.dat",
        sys.vals.num_ip_each_axis);

    FILE* fp;

   if(monolis_mpi_get_global_my_rank() == 0){
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_velocity.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_velocity.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_pressure.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_pressure.txt", sys.cond.directory);
            fclose(fp);

            //fclose(fp);
    }

    BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

        BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

        //intialize for velocity and pressure
        //initialize_velocity_pressure_karman_vortex(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);
    initialize_velocity_pressure_ahmedbody(&(sys.surf), sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);
    double* vec = BB_std_calloc_1d_double(vec, sys.fe.total_num_nodes*4);

    BBFE_fluid_sups_add_velocity_pressure(
            sys.vals.v,
            sys.vals.p,
            vec,
            sys.fe.total_num_nodes);

    ROM_sys_hlpod_fe_add_Dbc(
            vec,
            &(sys.bc),
            sys.fe.total_num_nodes,
            4);

    BBFE_fluid_sups_renew_velocity(
            sys.vals.v,
            vec,
            sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
            sys.vals.p,
            vec,
            sys.fe.total_num_nodes);

    ROM_BB_vec_copy_2d(
            sys.vals.v,
            sys.vals.v_old,
            sys.fe.total_num_nodes,
            3);

    ROM_BB_vec_copy_2d(
            sys.vals.v_old,
            sys.vals.v,
            sys.fe.total_num_nodes,
            3);

output_files(&sys, 0, 0);
        /****************** solver ********************/
        double t = 0.0;
        int step = 0;
        int file_num = 0;
        int count = 0;  //for ROM

        t = 0.0; step = 0; file_num = 0;

/*
        //if(sys.rom_prm_p.hot_start == 1){
            char fname[BUFFER_SIZE];
            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
            double* val = BB_std_calloc_1d_double(val, 4*sys.fe.total_num_nodes);
            double t_hs = hot_start_read_initialize_val(val, fname, sys.cond.directory);
            int step_hs = 0;

            printf("Hot start time: %lf\n", t);
            //printf("Hot start step: %d\n", step_rom);
            printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            //BBFE_fluid_sups_renew_velocity(sys.vals_rom.v, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_pressure(sys.vals_rom.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        //}
*/

const double U_ref = 60.0;
const double A_ref = 0.112032;

const double eU[3] = {
    1.0,
    0.0,
    0.0
};

const double eP[3] = {
    0.0,
    1.0,
    0.0
};

const double Uw[3] = {
    0.0,
    0.0,
    0.0
};

/*
 * solverが実際に使用する値と一致させる。
 */
const double gamma_n = 100.0;

/*
 * デバッグ中は100ステップごと程度にする。
 * 毎ステップowner mapを再構築すると非常に重い。
 */
const int diagnostic_interval = 1;

    while (t < sys.vals.finish_time) {
        t += sys.vals.dt;
        step += 1;

        if (monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "\n%s ----------------- step %d ----------------\n",
                CODENAME,
                step);

            fflush(stdout);
        }

        /* ================================================================ */
        /* Solver                                                           */
        /* ================================================================ */

        t4_diag_trace(
            "before solver_fom_VMS",
            step);

        solver_fom_VMS(
            sys,
            &(sys.surf),
            &(sys.basis_surf),
            t,
            count);

        t4_diag_trace(
            "after solver_fom_VMS",
            step);

        count++;

        /* ================================================================ */
        /* Heavy diagnostics                                                */
        /* ================================================================ */

        const int do_diagnostics =
            (
                step == 1 ||
                step % diagnostic_interval == 0
            );

        if (do_diagnostics) {
            t4_diag_trace(
                "before parallel diagnostics wrapper",
                step);

            const int diagnostic_status =
                run_wall_and_continuity_diagnostics_parallel_safe(
                    &sys,
                    step,
                    t,
                    U_ref,
                    Uw,
                    eU,
                    eP,
                    gamma_n,
                    1);

            t4_diag_trace(
                "after parallel diagnostics wrapper",
                step);

            if (
                diagnostic_status != 0 &&
                monolis_mpi_get_global_my_rank() == 0) {

                fprintf(
                    stderr,
                    "[main] diagnostics failed at step %d\n",
                    step);
            }
        }

        /* ================================================================ */
        /* VTK output                                                       */
        /* ================================================================ */

        /*
        * output_filesはpartitionごとのVTK出力を行うので、
        * 原則として全rankが呼ぶ。
        *
        * yPlusを更新してから呼ぶ。
        */
        if (step % sys.vals.output_interval == 0) {
            t4_diag_trace(
                "before output_files",
                step);

            output_files(
                &sys,
                step,
                t);

            t4_diag_trace(
                "after output_files",
                step);

            file_num++;
        }

        /* ================================================================ */
        /* Existing force / gamma / wall diagnostics                        */
        /* ================================================================ */

        /*
        * 既存の重い診断も毎ステップではなく、
        * diagnostic_intervalごとに変更する。
        */
        if (do_diagnostics) {
            double D_out = 0.0;
            double L_out = 0.0;
            double Cd_out = 0.0;
            double Cl_out = 0.0;

            double D_nitsche_out = 0.0;
            double L_nitsche_out = 0.0;
            double Cd_nitsche_out = 0.0;
            double Cl_nitsche_out = 0.0;

            double Cd_out_p = 0.0;
            double Cl_out_p = 0.0;

            const double qA =
                0.5
                * sys.vals.density
                * U_ref
                * U_ref
                * A_ref;

            /* ------------------------------------------------------------ */
            /* Local force calculation                                      */
            /* ------------------------------------------------------------ */

            const int force_status_local =
                calc_Cd_v_tet4_tri3_nitsche(
                    &(sys.surf_internal),
                    &(sys.fe),
                    &(sys.basis_surf_internal),
                    &(sys.mono_com),
                    &(sys.vals),
                    sys.vals.density,
                    sys.vals.viscosity,
                    U_ref,
                    A_ref,
                    eU,
                    eP,
                    Uw,
                    gamma_n,
                    &D_out,
                    &L_out,
                    &Cd_out,
                    &Cl_out,
                    &D_nitsche_out,
                    &L_nitsche_out,
                    &Cd_nitsche_out,
                    &Cl_nitsche_out);

            /*
            * 全rankでstatusを統一する。
            */
            const int force_global_ok =
                t4_all_ranks_status_ok(
                    force_status_local,
                    &(sys.mono_com),
                    "viscous/Nitsche force diagnostics");

            if (force_global_ok) {
                /*
                * 必ず全rankが同じ順番で呼ぶ。
                */
                monolis_allreduce_R(
                    1,
                    &Cd_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cl_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cd_nitsche_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cl_nitsche_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);
            }

            /* ------------------------------------------------------------ */
            /* Pressure force                                               */
            /* ------------------------------------------------------------ */

            calc_Cd_p_hex_or_tet(
                &(sys.surf_internal),
                &(sys.fe),
                &(sys.basis_surf_internal),
                &(sys.vals),
                eU,
                eP,
                &Cd_out_p,
                &Cl_out_p);

            monolis_allreduce_R(
                1,
                &Cd_out_p,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);

            monolis_allreduce_R(
                1,
                &Cl_out_p,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);
            
            /* ================================================================ */
            /* Nitsche residual-row reaction diagnostic                         */
            /* ここへ t4n_nitsche_row_sum_main_patch.c の内容を追加              */
            /* ================================================================ */

            T4NitscheRowSumDiagnostics row_diag;

            const int row_status_local =
                calc_nitsche_row_sum_reaction_tet4_tri3_local(
                    &(sys.surf_internal),
                    &(sys.fe),
                    &(sys.basis_surf_internal),
                    &(sys.vals),
                    sys.vals.density,
                    sys.vals.viscosity,
                    Uw,
                    gamma_n,
                    &row_diag);

            const int row_global_ok =
                t4_all_ranks_status_ok(
                    row_status_local,
                    &(sys.mono_com),
                    "Nitsche residual-row reaction");

            if (row_global_ok) {
                nitsche_row_sum_diagnostics_allreduce(
                    &row_diag,
                    &(sys.mono_com));
            }

            /* ------------------------------------------------------------ */
            /* Shared text output: rank 0 only                              */
            /* ------------------------------------------------------------ */

            if (monolis_mpi_get_global_my_rank() == 0) {
                if (force_global_ok) {
                    ROM_std_hlpod_output_add_calc_time(
                        Cd_out / qA,
                        t,
                        "cylinder_drag_coeff_velocity.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cl_out / qA,
                        t,
                        "cylinder_lift_coeff_velocity.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cd_nitsche_out / qA,
                        t,
                        "cylinder_drag_coeff_Nitsche.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cl_nitsche_out / qA,
                        t,
                        "cylinder_lift_coeff_Nitsche.txt",
                        sys.cond.directory);
                }

                ROM_std_hlpod_output_add_calc_time(
                    Cd_out_p / qA,
                    t,
                    "cylinder_drag_coeff_pressure.txt",
                    sys.cond.directory);

                ROM_std_hlpod_output_add_calc_time(
                    Cl_out_p / qA,
                    t,
                    "cylinder_lift_coeff_pressure.txt",
                    sys.cond.directory);
            }
        }

        /* ================================================================ */
        /* Hot-start output, when enabled                                   */
        /* ================================================================ */

        if(step%1 == 0){
                char fname[BUFFER_SIZE];
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
                //hot_start_write_initialize_val(sys.monolis.mat.R.X, sys.fe.total_num_nodes, 4, t, fname, sys.cond.directory);
            }
        else{
        }
    }

        BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
        BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), sys.fe.total_num_nodes, 4);
        monolis_finalize(&(sys.monolis));

        double t2 = monolis_get_time();
        int myrank = monolis_mpi_get_global_my_rank();

        if(myrank == 0) {
                printf("** Total time: %f\n", t2 - t1);
        }

        monolis_global_finalize();

        printf("\n");

        return 0;
}