#include "core_FOM_mixed.h"
#include "core_Nitsche_prism.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* INPUT_FILENAME_D_BC_V = "D_bc_v.dat";
static const char* INPUT_FILENAME_SURF_GRAPH = "surf_graph.dat";

/*
 * TET4 + PRI6 volume mesh with TRI3 cylinder wall owned by PRI6.
 *
 * IMPORTANT:
 *   Cylinder_wall velocity DOFs must NOT also appear as strong Dirichlet rows
 *   in D_bc_v.dat when Nitsche is enabled.
 */
int main(int argc, char* argv[])
{
    FE_SYSTEM_FLUID_MIXED sys;
    memset(&sys, 0, sizeof(sys));

    monolis_global_initialize();

    sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);
    read_calc_conditions(&(sys.vals), sys.cond.directory);

    BBFE_fluid_pre_mixed(
        &(sys.fe_tet), &(sys.basis_tet),
        &(sys.fe_pri), &(sys.basis_pri),
        argc, argv,
        sys.cond.directory,
        sys.vals.num_ip_each_axis);

    const int nnode = sys.fe_tet.total_num_nodes;

    memory_allocation_nodal_values_AB2(&(sys.vals), nnode);

    const char* filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        INPUT_FILENAME_D_BC_V);

    BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
        &(sys.bc), &(sys.vals), filename, sys.cond.directory, nnode, 4);

    BBFE_fluid_sups_read_Dirichlet_bc_NR(
        &(sys.bc_NR), filename, sys.cond.directory, nnode, 4);

    /* Verify that PRI6 uy is not accidentally constrained by D_bc_v.dat. */
    BBFE_fluid_mixed_report_prism_bc(&sys);

    /* ------------------------------------------------------------ */
    /* Geometry: each rank can have TET only, PRI only, both, or none. */
    /* ------------------------------------------------------------ */

    if(BBFE_fluid_mixed_has_elements(&(sys.fe_tet))) {
        BBFE_elemmat_set_Jacobi_mat(&(sys.fe_tet), &(sys.basis_tet));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_tet), &(sys.basis_tet));
    }

    if(BBFE_fluid_mixed_has_elements(&(sys.fe_pri))) {
        BBFE_elemmat_set_Jacobi_mat(&(sys.fe_pri), &(sys.basis_pri));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_pri), &(sys.basis_pri));
    }

    /* We already checked this externally; keep the same check in the solver. */
    BBFE_fluid_mixed_check_positive_jacobian(
        &(sys.fe_tet), &(sys.basis_tet), "TET4", 0.0);

    BBFE_fluid_mixed_check_positive_jacobian(
        &(sys.fe_pri), &(sys.basis_pri), "PRI6", 0.0);

    /* ------------------------------------------------------------ */
    /* Unified TET4 + PRI6 Monolis graph.                           */
    /* ------------------------------------------------------------ */

    BBFE_fluid_init_monomat_mixed(
        &(sys.monolis), &(sys.mono_com), &(sys.fe_tet), sys.cond.directory);

    BBFE_fluid_mixed_report_prism_partition(
        &(sys.fe_pri), &(sys.mono_com));

    /* ------------------------------------------------------------ */
    /* TRI3 cylinder wall. Empty surface files on non-PRI ranks are valid. */
    /* ------------------------------------------------------------ */

    BBFE_fluid_pre_surface_mixed(
        &(sys.surf),
        &(sys.basis_surf),
        sys.cond.directory,
        INPUT_FILENAME_SURF_GRAPH,
        sys.vals.num_ip_each_axis);

    int* owner = NULL;
    int* lface = NULL;
    int local_unmapped = 0;

    if(BBFE_pri6_tri3_make_owner_map(
        &(sys.surf),
        &(sys.fe_pri),
        &owner,
        &lface,
        &local_unmapped) != 0) {

        fprintf(stderr,
            "%s ERROR: PRI6/TRI3 owner-map precheck failed on rank %d.\n",
            CODENAME,
            monolis_mpi_get_global_my_rank());

        exit(EXIT_FAILURE);
    }

    printf(
        "%s rank=%d wall owner-map: surf=%d pri=%d mapped=%d unmapped=%d\n",
        CODENAME,
        monolis_mpi_get_global_my_rank(),
        sys.surf.total_num_elems,
        sys.fe_pri.total_num_elems,
        sys.surf.total_num_elems - local_unmapped,
        local_unmapped);

    int global_unmapped = local_unmapped;
    int global_surface = sys.surf.total_num_elems;
    int global_mapped = sys.surf.total_num_elems - local_unmapped;

    monolis_allreduce_I(
        1, &global_unmapped, MONOLIS_MPI_SUM, sys.mono_com.comm);
    monolis_allreduce_I(
        1, &global_surface, MONOLIS_MPI_SUM, sys.mono_com.comm);
    monolis_allreduce_I(
        1, &global_mapped, MONOLIS_MPI_SUM, sys.mono_com.comm);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s global wall owner-map: surf=%d mapped=%d unmapped=%d\n",
            CODENAME,
            global_surface,
            global_mapped,
            global_unmapped);
    }

    free(owner);
    free(lface);

    if(global_unmapped != 0) {
        fprintf(stderr,
            "%s ERROR: at least one TRI3 wall face has no local PRI6 owner.\n",
            CODENAME);
        exit(EXIT_FAILURE);
    }

    /* ------------------------------------------------------------ */
    /* Initial condition + strong BCs other than the cylinder wall. */
    /* ------------------------------------------------------------ */

    initialize_velocity_pressure_karman_vortex(
        sys.vals.v,
        sys.vals.p,
        nnode);

    double* initial_vec = BB_std_calloc_1d_double(NULL, 4*nnode);

    BBFE_fluid_sups_add_velocity_pressure(
        sys.vals.v,
        sys.vals.p,
        initial_vec,
        nnode);

    for(int i = 0; i < 4*nnode; ++i) {
        if(sys.bc.D_bc_exists[i]) {
            initial_vec[i] = sys.bc.imposed_D_val[i];
        }
    }

    /* Synchronize overlap copies before creating v_old. */
    BBFE_fluid_mixed_sync_solution_vector(
        &(sys.mono_com), nnode, initial_vec);

    BBFE_fluid_sups_renew_velocity(
        sys.vals.v, initial_vec, nnode);

    BBFE_fluid_sups_renew_pressure(
        sys.vals.p, initial_vec, nnode);

    BB_std_free_1d_double(initial_vec, 4*nnode);

    BBFE_fluid_mixed_sync_state(
        &(sys.mono_com), &(sys.vals), nnode);

    ROM_BB_vec_copy_2d(
        sys.vals.v,
        sys.vals.v_old,
        nnode,
        3);

    /* ------------------------------------------------------------ */
    /* PRI6 Nitsche wall.                                          */
    /* ------------------------------------------------------------ */

    PRI6WallOptions wall;
    BBFE_pri6_wall_options_default(&wall);

    wall.wall_mode = PRI6_WALL_NITSCHE_NOSLIP;
    wall.exchange_mode = PRI6_WALL_EXCHANGE_OPPOSITE_FACE;
    wall.gamma_n = 10000.0;
    wall.dt_penalty_coeff = 0.10;

    double t = 0.0;
    int step = 0;
    int file_num = 0;

    /* t=0 output: IsPrismNode and ElementFamily make visual debugging easier. */
    output_files_mixed(&sys, file_num, t);
    ++file_num;

    while(t < sys.vals.finish_time) {
        t += sys.vals.dt;
        ++step;

        solver_fom_NR_mixed_prism_nitsche(
            &sys,
            &(sys.surf),
            &(sys.basis_surf),
            &wall,
            t,
            step);

        if(step % sys.vals.output_interval == 0) {
            output_files_mixed(&sys, file_num, t);
            ++file_num;
        }
    }

    /* ------------------------------------------------------------ */
    /* Finalize.                                                    */
    /* ------------------------------------------------------------ */

    BBFE_fluid_finalize_surface_mixed(
        &(sys.surf), &(sys.basis_surf));

    BBFE_fluid_finalize_mixed(
        &(sys.fe_tet), &(sys.basis_tet),
        &(sys.fe_pri), &(sys.basis_pri));

    BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), nnode, 4);
    BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc_NR), nnode, 4);

    monolis_finalize(&(sys.monolis));
    monolis_global_finalize();

    return 0;
}
