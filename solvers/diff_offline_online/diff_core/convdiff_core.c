
#include "convdiff_core.h"


const char* BBFE_convdiff_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename)
{
	const char* dir_name;

	if(argc < 2) { dir_name = "."; }
	else         { dir_name = argv[1]; }

	printf("%s Main directory: %s\n", codename, dir_name);

	return dir_name;
}


void BBFE_convdiff_pre(
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
			filename,
			directory);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_ELEM);
	BBFE_sys_read_elem(
			fe,
			filename,
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
			filename,
			directory,
			fe->total_num_nodes,
			BLOCK_SIZE);

	BBFE_convdiff_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);

	monolis_initialize(monolis);

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
}


void BBFE_convdiff_set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis)
{
	switch( local_num_nodes ) {
		case 4:
			basis->num_integ_points =
				BBFE_std_integ_tet_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_tet1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_tet1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order tetrahedron.\n", CODENAME);
			break;

		case 8:
			basis->num_integ_points =
				BBFE_std_integ_hex_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_hex1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_hex1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order hexahedron.\n", CODENAME);
			break;
	}
	printf("%s The number of integration points: %d\n", CODENAME, basis->num_integ_points);
}


double BBFE_convdiff_equivval_relative_L2_error_scalar(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		MONOLIS_COM*  monolis_com,
		double        t,
		const double* comp_vec, // [total_num_nodes]
		double        (*func)(double, double, double, double)) // scalar function(x, y, z, t)
{
	double L2_abs_error = 0.0;
	double L2_abs_theo  = 0.0;

	double* val_ip_error;
	double* val_ip_theo;
	double* Jacobian_ip;
	bool* is_internal_elem;
	val_ip_error = BB_std_calloc_1d_double(val_ip_error, basis->num_integ_points);
	val_ip_theo  = BB_std_calloc_1d_double(val_ip_theo , basis->num_integ_points);
	Jacobian_ip  = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x  , fe->local_num_nodes, 3);

	double* local_val;
	local_val = BB_std_calloc_1d_double(local_val, fe->local_num_nodes);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_val[i] = comp_vec[ fe->conn[e][i] ];

			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		for(int p=0; p<(basis->num_integ_points); p++) {
			double val_ip;
			val_ip = BBFE_std_mapping_scalar(
					fe->local_num_nodes,
					local_val,
					basis->N[p]);

			double x_ip[3];
			BBFE_std_mapping_vector3d(
					x_ip,
					fe->local_num_nodes,
					local_x,
					basis->N[p]);

			val_ip_error[p] =
				pow( fabs( val_ip - func(x_ip[0], x_ip[1], x_ip[2], t) ), 2 );
			val_ip_theo[p] =
				pow( fabs( func(x_ip[0], x_ip[1], x_ip[2], t) ), 2 );
		}

		double integ_val_error = BBFE_std_integ_calc(
				basis->num_integ_points,
				val_ip_error,
				basis->integ_weight,
				Jacobian_ip);
		double integ_val_theo  = BBFE_std_integ_calc(
				basis->num_integ_points,
				val_ip_theo,
				basis->integ_weight,
				Jacobian_ip);

		L2_abs_error += integ_val_error;
		L2_abs_theo  += integ_val_theo;
	}

	BB_std_free_1d_double(val_ip_error, basis->num_integ_points);
	BB_std_free_1d_double(val_ip_theo,  basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip,  basis->num_integ_points);
	BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_val, fe->local_num_nodes);

	monolis_allreduce_R(
			1,
			&L2_abs_error,
			MONOLIS_MPI_SUM,
			monolis_com->comm);

	monolis_allreduce_R(
			1,
			&L2_abs_theo,
			MONOLIS_MPI_SUM,
			monolis_com->comm);

	return ( sqrt(L2_abs_error)/sqrt(L2_abs_theo) );
}


void BBFE_convdiff_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		BBFE_BC*     bc)
{
	BBFE_sys_memory_free_integ(basis, 3);
	BBFE_sys_memory_free_shapefunc(basis);

	BBFE_sys_memory_free_node(fe, 3);
	BBFE_sys_memory_free_elem(fe, basis->num_integ_points, 3);

	BBFE_sys_memory_free_Dirichlet_bc(bc, fe->total_num_nodes, BLOCK_SIZE);
}
