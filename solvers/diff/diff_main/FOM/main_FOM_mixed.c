#include "core_FOM_mixed.h"
#include "core_ROM.h"

#include <stdio.h>

int main(
        int   argc,
        char* argv[])
{
    printf("\n");

    FE_SYSTEM_MIXED sys;

    monolis_global_initialize();
    const double t1 = monolis_get_time();

    sys.cond.directory = BBFE_convdiff_get_directory_name(
            argc, argv, CODENAME);

    read_calc_conditions(
            &(sys.vals),
            sys.cond.directory);

    /*
     * Reads partition-aware files:
     *   node.dat
     *   elem_tet.dat
     *   elem_prism.dat
     *   D_bc.dat
     *
     * In a partitioned run monolis_get_global_input_file_name() resolves the
     * rank-specific files under parted.0/.
     */
    BBFE_convdiff_pre_mixed(
            &(sys.fe_tet),
            &(sys.basis_tet),
            &(sys.fe_pri),
            &(sys.basis_pri),
            &(sys.bc),
            &(sys.monolis),
            &(sys.monolis_com),
            sys.cond.directory,
            sys.vals.num_ip_each_axis);

    /* node.dat is common to tetrahedra and prisms */
    memory_allocation_nodal_values(
            &(sys.vals),
            sys.fe_tet.total_num_nodes);

    manusol_set_init_value(
            &(sys.fe_tet),
            sys.vals.T);

    if(monolis_mpi_get_global_my_rank() == 0) {
        FILE* fp = NULL;
        fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "calctime/time_fem.txt", sys.cond.directory);
        fclose(fp);

        ROM_std_hlpod_write_solver_prm_fopen(
                "fem_solver_prm",
                sys.cond.directory);
    }

    /* element geometry for both element families */
    BBFE_elemmat_set_Jacobi_mat(
            &(sys.fe_tet),
            &(sys.basis_tet));

    BBFE_elemmat_set_shapefunc_derivative(
            &(sys.fe_tet),
            &(sys.basis_tet));

    BBFE_elemmat_set_Jacobi_mat(
            &(sys.fe_pri),
            &(sys.basis_pri));

    BBFE_elemmat_set_shapefunc_derivative(
            &(sys.fe_pri),
            &(sys.basis_pri));

    /*
     * A single sparse matrix pattern is built from graph.dat, which already
     * contains the union of tetrahedron and prism nodal adjacency.
     */
    monolis_initialize(&(sys.monolis0));

    BBFE_convdiff_set_nonzero_pattern_mixed(
            &(sys.monolis0),
            &(sys.fe_tet),
            sys.cond.directory);

    /* Add both element families to the same matrix. */
    set_element_mat_mixed(
            &(sys.monolis0),
            &(sys.fe_tet),
            &(sys.basis_tet),
            &(sys.fe_pri),
            &(sys.basis_pri),
            &(sys.vals));

    monolis_copy_mat_R(
            &(sys.monolis0),
            &(sys.monolis));

    /* ---------------- time integration ---------------- */
    double t = 0.0;
    int step = 0;
    int file_num = 0;

    while(t < sys.vals.finish_time) {
        t += sys.vals.dt;
        step += 1;

        const double calctime_fem_t1 = monolis_get_time_global_sync();

        solver_fom_mixed(
                &sys,
                t,
                step);

        const double calctime_fem_t2 = monolis_get_time_global_sync();

        if(step % sys.vals.output_interval == 0) {
            output_files_mixed(
                    &sys,
                    file_num,
                    t);

            if(monolis_mpi_get_global_my_rank() == 0) {
                ROM_std_hlpod_write_solver_prm(
                        &(sys.monolis),
                        t,
                        "fem_solver_prm/",
                        sys.cond.directory);

                ROM_std_hlpod_output_add_calc_time(
                        calctime_fem_t2 - calctime_fem_t1,
                        t,
                        "calctime/time_fem.txt",
                        sys.cond.directory);
            }

            file_num += 1;
        }
    }

    /* ---------------- cleanup ---------------- */
    BBFE_convdiff_finalize_mixed(
            &(sys.fe_tet),
            &(sys.basis_tet),
            &(sys.fe_pri),
            &(sys.basis_pri),
            &(sys.bc));

    monolis_finalize(&(sys.monolis));
    monolis_finalize(&(sys.monolis0));

    monolis_global_finalize();

    const double t2 = monolis_get_time();
    printf("** Total time: %f\n", t2 - t1);
    printf("\n");

    return 0;
}
