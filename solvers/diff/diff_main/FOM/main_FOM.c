
#include "core_ROM.h"

static const char* INPUT_FILENAME_COND    = "cond.dat";

int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	manusol_set_init_value(&(sys.fe), sys.vals.T);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	if(monolis_mpi_get_global_my_rank() == 0){
	    fp = ROM_BB_write_fopen(fp, "calctime/time_fem.txt", sys.cond.directory);
		fclose(fp);
        ROM_std_hlpod_write_solver_prm_fopen("fem_solver_prm", sys.cond.directory);
	}

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));

	monolis_initialize(&(sys.monolis0));

	monolis_com_initialize_by_parted_files(
			&(sys.monolis_com),
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_R(
			&(sys.monolis0),
			sys.fe.total_num_nodes,
			sys.fe.local_num_nodes,
			1,
			sys.fe.total_num_elems,
			sys.fe.conn);

	set_element_mat(
			&(sys.monolis0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));

    monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;
	double calctime_fem_t1 = monolis_get_time_global_sync();
        solver_fom(sys, t, step);
	double calctime_fem_t2 = monolis_get_time_global_sync();

		if(step%sys.vals.output_interval == 0) {
			output_files(&sys, file_num, t);

            if(monolis_mpi_get_global_my_rank()==0){
                ROM_std_hlpod_write_solver_prm(&(sys.monolis), t, "fem_solver_prm/" , sys.cond.directory);
                ROM_std_hlpod_output_add_calc_time(calctime_fem_t2-calctime_fem_t1, t,
                        "calctime/time_fem.txt", sys.cond.directory);
            }

			file_num += 1;
		}
	}
    /**************************************************/

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis0));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
