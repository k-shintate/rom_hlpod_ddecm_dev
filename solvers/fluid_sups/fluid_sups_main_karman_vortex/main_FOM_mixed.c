#include "core_FOM_mixed.h"
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

static const char* INPUT_FILENAME_D_BC_V = "D_bc_v.dat";

/*
 * Mixed TET4 + PRI6 version of the Karman-vortex FOM driver.
 *
 * Required volume files:
 *   node.dat
 *   elem_tet.dat
 *   elem_prism.dat
 *   graph.dat
 *   D_bc_v.dat
 *
 * The graph must contain the union of tetrahedron and prism nodal adjacency.
 * For MPI runs, the graph-partition preprocessing must create the matching
 * parted.0 files for node/graph/elem_tet/elem_prism/D_bc_v.
 */
int main(int argc, char* argv[])
{
    printf("\n");

    FE_SYSTEM_FLUID_MIXED sys;
    memset(&sys, 0, sizeof(sys));

    monolis_global_initialize();
    const double t1 = monolis_get_time();

    sys.cond.directory =
        BBFE_fluid_get_directory_name(argc, argv, CODENAME);

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

    /* Physical velocity values (used to impose the initial constrained state). */
    BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
        &(sys.bc),
        &(sys.vals),
        filename,
        sys.cond.directory,
        nnode,
        FLUID_MIXED_BLOCK_SIZE);

    /* Newton increment BC: constrained increments are zero. */
    BBFE_fluid_sups_read_Dirichlet_bc_NR(
        &(sys.bc_NR),
        filename,
        sys.cond.directory,
        nnode,
        FLUID_MIXED_BLOCK_SIZE);

    /* Geometry at integration points for both cell types. */
    BBFE_elemmat_set_Jacobi_mat(&(sys.fe_tet), &(sys.basis_tet));
    BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_tet), &(sys.basis_tet));

    BBFE_elemmat_set_Jacobi_mat(&(sys.fe_pri), &(sys.basis_pri));
    BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_pri), &(sys.basis_pri));

    /* One 4-DOF/node sparse matrix from the unified TET4+PRI6 graph. */
    BBFE_fluid_init_monomat_mixed(
        &(sys.monolis),
        &(sys.mono_com),
        &(sys.fe_tet),
        sys.cond.directory);

    /* Karman-vortex initial condition. */
    initialize_velocity_pressure_karman_vortex(
        sys.vals.v,
        sys.vals.p,
        nnode);

    ROM_BB_vec_copy_2d(
        sys.vals.v,
        sys.vals.v_old,
        nnode,
        3);

    /* Put the physical Dirichlet values into the initial field using the same
       path as the existing FOM driver. */
    double* initial_vec = BB_std_calloc_1d_double(NULL, 4*nnode);

    BBFE_fluid_sups_add_velocity_pressure(
        sys.vals.v,
        sys.vals.p,
        initial_vec,
        nnode);

    for(int i = 0; i < nnode * FLUID_MIXED_BLOCK_SIZE; ++i) {
        if(sys.bc.D_bc_exists[i]) {
            initial_vec[i] = sys.bc.imposed_D_val[i];
        }
    }

    BBFE_fluid_sups_renew_velocity(
        sys.vals.v,
        initial_vec,
        nnode);

    BBFE_fluid_sups_renew_pressure(
        sys.vals.p,
        initial_vec,
        nnode);

    ROM_BB_vec_copy_2d(
        sys.vals.v,
        sys.vals.v_old,
        nnode,
        3);

    BB_std_free_1d_double(initial_vec, 4*nnode);

    double t = 0.0;
    int step = 0;
    int file_num = 0;

    while(t < sys.vals.finish_time) {
        t += sys.vals.dt;
        ++step;

        solver_fom_NR_mixed(&sys, t, step);

        if(step % sys.vals.output_interval == 0) {
            output_files_mixed(&sys, file_num, t);
            ++file_num;
        }

        /* Keep monolis.X synchronized with the current physical state for
           hot-start/output routines that use the 4-DOF packed vector. */
        BBFE_fluid_sups_add_velocity_pressure(
            sys.vals.v,
            sys.vals.p,
            sys.monolis.mat.R.X,
            nnode);

        if(step % 1000 == 0) {
            /* This output routine only consumes nodal fields/coordinates, so
               the tetra FE object is a valid representative of shared nodes. */
            output_result_file_karman_vortex_pressure(
                &(sys.fe_tet),
                &(sys.vals),
                t,
                sys.cond.directory);
        }
    }

    BBFE_fluid_finalize_mixed(
        &(sys.fe_tet), &(sys.basis_tet),
        &(sys.fe_pri), &(sys.basis_pri));

    BBFE_sys_memory_free_Dirichlet_bc(
        &(sys.bc), nnode, FLUID_MIXED_BLOCK_SIZE);
    BBFE_sys_memory_free_Dirichlet_bc(
        &(sys.bc_NR), nnode, FLUID_MIXED_BLOCK_SIZE);

    monolis_finalize(&(sys.monolis));

    const double t2 = monolis_get_time();
    if(monolis_mpi_get_global_my_rank() == 0) {
        printf("** Mixed fluid total time: %f\n", t2-t1);
    }

    monolis_global_finalize();

    printf("\n");
    return 0;
}
