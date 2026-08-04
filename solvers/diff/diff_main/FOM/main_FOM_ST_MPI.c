#include "core_ROM.h"
#include "core_FOM_ST.h"
#include "core_FOM_ST_MPI.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ST_MPI_DEGREE
#define ST_MPI_DEGREE 1
#endif

#ifndef ST_MPI_WINDOW_SLABS
#define ST_MPI_WINDOW_SLABS 4
#endif


void BBFE_convdiff_pre_global(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution)
{
	BB_calc_void();

	int n_axis = num_integ_points_each_axis;
	const char* filename;

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_NODE);
	BBFE_sys_read_node(
			fe,
			"node.dat",
			directory);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_ELEM);
	BBFE_sys_read_elem(
			fe,
			"elem.dat",
			directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			fe->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
	BBFE_sys_read_Dirichlet_bc(
			bc,
			"D_bc.dat",
			directory,
			fe->total_num_nodes,
			BLOCK_SIZE);

	BBFE_convdiff_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);

	monolis_initialize(monolis);
	
/*
	monolis_com_initialize_by_parted_files(
			monolis_com,
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_R(
			monolis,
			fe->total_num_nodes,
			fe->local_num_nodes,
			1,
			fe->total_num_elems,
			fe->conn);
*/
}

int main(int argc, char* argv[])
{
	printf("\n");

	FE_SYSTEM sys = {0};
	ST_PARTITION st;
	ST_SPACE_MAP space_map;
	ST_partition_initialize(&st);
	ST_space_map_initialize(&space_map);

	monolis_global_initialize();
	const int my_rank = monolis_mpi_get_global_my_rank();

	const double time_start_total = monolis_get_time();

	/* ------------------------------------------------------------------ */
	/* Ordinary spatial finite-element partition used for element          */
	/* integration, boundary flags and output.                              */
	/* ------------------------------------------------------------------ */

	sys.cond.directory = BBFE_convdiff_get_directory_name(
			argc,
			argv,
			CODENAME);

	read_calc_conditions(&sys.vals, sys.cond.directory);

	BBFE_convdiff_pre_global(
			&sys.fe,
			&sys.basis,
			&sys.bc,
			&sys.monolis,
			&sys.monolis_com,
			argc,
			argv,
			sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);
	
	monolis_com_initialize_by_self(
        &sys.monolis_com);

	memory_allocation_nodal_values(
			&sys.vals,
			sys.fe.total_num_nodes);

	ST_space_map_initialize(&space_map);

	ST_space_map_build_identity_from_mesh(
			&space_map,
			&sys.fe);

	if(space_map.num_total_space_nodes != sys.fe.total_num_nodes ||
	   space_map.num_local_space_elements != sys.fe.total_num_elems) {
		fprintf(
				stderr,
				"ERROR: spatial map and BBFE mesh dimensions differ.\n"
				"  map nodes/elements = %d/%d\n"
				"  BBFE nodes/elements = %d/%d\n",
				space_map.num_total_space_nodes,
				space_map.num_local_space_elements,
				sys.fe.total_num_nodes,
				sys.fe.total_num_elems);
		exit(EXIT_FAILURE);
	}

	/* ------------------------------------------------------------------ */
	/* Output initialization.                                               */
	/* ------------------------------------------------------------------ */

	FILE* fp = NULL;
	fp = BBFE_sys_write_fopen(
			fp,
			"l2_error.txt",
			sys.cond.directory);
	fclose(fp);

	if(my_rank == 0) {
		fp = ROM_BB_write_fopen(
				fp,
				"calctime/time_fem.txt",
				sys.cond.directory);
		fclose(fp);
		ROM_std_hlpod_write_solver_prm_fopen(
				"fem_solver_prm",
				sys.cond.directory);
	}

	/* ------------------------------------------------------------------ */
	/* Spatial element geometry.                                            */
	/* ------------------------------------------------------------------ */

	BBFE_elemmat_set_Jacobi_mat(&sys.fe, &sys.basis);
	BBFE_elemmat_set_shapefunc_derivative(&sys.fe, &sys.basis);

	/* ------------------------------------------------------------------ */
	/* Temporal element and partitioned ST data.                            */
	/* ------------------------------------------------------------------ */

	const int degree = ST_MPI_DEGREE;
	const int num_window_slabs = ST_MPI_WINDOW_SLABS;

	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG(&te, degree);

	ST_partition_read(
			&st,
			sys.cond.directory,
			&te,
			num_window_slabs,
			&sys.monolis_com);

	ST_validate_partition_against_space_mesh(
			&st,
			&space_map,
			&sys.fe);

	if(sys.vals.dt <= 0.0) {
		fprintf(stderr, "ERROR: dt must be positive.\n");
		exit(EXIT_FAILURE);
	}
	if(sys.vals.finish_time <= 0.0) {
		fprintf(stderr, "ERROR: finish_time must be positive.\n");
		exit(EXIT_FAILURE);
	}
	if(sys.vals.output_interval <= 0) {
		fprintf(stderr, "ERROR: output_interval must be positive.\n");
		exit(EXIT_FAILURE);
	}

	const int total_num_steps = (int)llround(
			sys.vals.finish_time / sys.vals.dt);
	const double reconstructed_finish_time =
		(double)total_num_steps * sys.vals.dt;
	const double time_tolerance = 1.0e-12 *
		fmax(1.0, fabs(sys.vals.finish_time));

	if(fabs(reconstructed_finish_time - sys.vals.finish_time)
	   > time_tolerance) {
		fprintf(
				stderr,
				"ERROR: finish_time / dt must be an integer.\n"
				"  finish_time = %.16e\n"
				"  dt          = %.16e\n",
				sys.vals.finish_time,
				sys.vals.dt);
		exit(EXIT_FAILURE);
	}

	if(total_num_steps % num_window_slabs != 0) {
		fprintf(
				stderr,
				"ERROR: total_num_steps must be divisible by "
				"num_window_slabs.\n"
				"  total_num_steps = %d\n"
				"  window_slabs    = %d\n",
				total_num_steps,
				num_window_slabs);
		exit(EXIT_FAILURE);
	}

	const int num_windows = total_num_steps / num_window_slabs;

	/* ------------------------------------------------------------------ */
	/* Replace the ordinary spatial MONOLIS matrix with the local           */
	/* partitioned ST matrix.  Keep sys.monolis_com for ordinary spatial    */
	/* output; the ST solver uses st.monolis_com_st.                         */
	/* ------------------------------------------------------------------ */

	monolis_finalize(&sys.monolis);
	monolis_initialize(&sys.monolis);
	monolis_initialize(&sys.monolis0);

	ST_set_nonzero_pattern_from_partition(&sys.monolis0, &st);
	monolis_clear_mat_value_rhs_R(&sys.monolis0);

	set_element_mat_ST_dG_window_partitioned(
			&sys.monolis0,
			&sys.fe,
			&sys.basis,
			&sys.vals,
			&te,
			&st,
			&space_map);

	monolis_copy_mat_R(&sys.monolis0, &sys.monolis);

	/* ------------------------------------------------------------------ */
	/* Solution and global trace vectors.                                   */
	/* ------------------------------------------------------------------ */

	const int num_local_st_dofs = st.num_total_st_nodes;
	const int num_global_space_nodes = st.meta.num_global_space_nodes;

	double* T_window = (double*)calloc(
			(size_t)num_local_st_dofs,
			sizeof(double));
	double* T_prev_trace_global = (double*)calloc(
			(size_t)num_global_space_nodes,
			sizeof(double));
	double* T_trace_global = (double*)calloc(
			(size_t)num_global_space_nodes,
			sizeof(double));

	if(T_window == NULL ||
	   T_prev_trace_global == NULL ||
	   T_trace_global == NULL) {
		fprintf(stderr, "ERROR: solution-vector allocation failed.\n");
		exit(EXIT_FAILURE);
	}

	/* Initial spatial state, gathered once by global spatial node ID. */
	ST_collect_space_nodal_vector_global(
			&space_map,
			&sys,
			sys.vals.T,
			T_prev_trace_global,
			num_global_space_nodes);

	/* ------------------------------------------------------------------ */
	/* Time-window loop.                                                    */
	/* ------------------------------------------------------------------ */

	int file_num = 0;
	for(int window=0; window<num_windows; window++) {
		const double t_window_start =
			(double)(window * num_window_slabs) * sys.vals.dt;

		memset(
				T_window,
				0,
				(size_t)num_local_st_dofs * sizeof(double));

		solver_fom_space_time_partitioned(
				&sys,
				&te,
				&st,
				&space_map,
				t_window_start,
				window + 1,
				T_prev_trace_global,
				T_window);

		for(int local_slab=0;
				local_slab<num_window_slabs;
				local_slab++) {
			const int global_step =
				window * num_window_slabs + local_slab + 1;
			const double t = (double)global_step * sys.vals.dt;

			if(global_step % sys.vals.output_interval == 0) {
				ST_extract_slab_right_trace_global(
						&te,
						&st,
						local_slab,
						T_window,
						T_trace_global);

				ST_copy_global_vector_to_local_space(
						&space_map,
						T_trace_global,
						sys.vals.T);

				if(my_rank == 0) {
					output_files(
							&sys,
							file_num,
							t);
				}

				file_num++;
			}
		}

		/* Right trace of the last slab becomes inflow to the next window. */
		ST_extract_slab_right_trace_global(
				&te,
				&st,
				num_window_slabs - 1,
				T_window,
				T_prev_trace_global);
	}

	ST_copy_global_vector_to_local_space(
			&space_map,
			T_prev_trace_global,
			sys.vals.T);

	/* ------------------------------------------------------------------ */
	/* Finalization.                                                        */
	/* ------------------------------------------------------------------ */

	free(T_window);
	free(T_prev_trace_global);
	free(T_trace_global);

	ST_partition_finalize(&st);
	ST_space_map_finalize(&space_map);

	BBFE_convdiff_finalize(&sys.fe, &sys.basis, &sys.bc);
	monolis_finalize(&sys.monolis);
	monolis_finalize(&sys.monolis0);

	const double time_end_total = monolis_get_time();
	if(my_rank == 0) {
		printf("** Total time: %f\n", time_end_total - time_start_total);
		printf("\n");
	}

	monolis_global_finalize();
	return 0;
}
