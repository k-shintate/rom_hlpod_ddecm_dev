#include "fluid_core.h"

const char* CODENAME = "fluid >";

const char* INPUT_FILENAME_NODE = "node.dat";
const char* INPUT_FILENAME_ELEM = "elem.dat";


const char* BBFE_fluid_get_directory_name(
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


void BBFE_fluid_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis)
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

	BBFE_fluid_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);
}


void BBFE_fluid_set_basis(
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


void BBFE_fluid_sups_renew_velocity(
		double**  v,
		double*   ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			v[i][d] = ans_vec[ 4*i + d ];
		}
	}
}


void BBFE_fluid_sups_renew_pressure(
		double*   p,
		double*   ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		p[i] = ans_vec[ 4*i + 3 ];
	}
}


void BBFE_fluid_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis)
{
	BBFE_sys_memory_free_integ(basis, 3);
	BBFE_sys_memory_free_shapefunc(basis);

	BBFE_sys_memory_free_node(fe, 3);
	BBFE_sys_memory_free_elem(fe, basis->num_integ_points, 3);
}


void BBFE_fluid_initialize_velocity_pressure(
	double**    v,
	double*     p,
	const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] = 0.0;
        }
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void BBFE_fluid_initialize_velocity_pressure_karman_vortex(
	double**    v,
    double**    v_old,
	double*     p,
	const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        v[i][0] = 1.0;
        v_old[i][0] = 1.0;
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void BBFE_fluid_sups_add_velocity_pressure(
                double**  v,
                double*   p,
                double*   sol_vec,
                const int total_num_nodes)
{
    for(int i=0; i<total_num_nodes; i++) {
            for(int d=0; d<3; d++) {
                    sol_vec[ 4*i + d ] = v[i][d];
            }
            sol_vec[ 4*i + 3 ] = p[i];
    }
}

void BBFE_fluid_add_Dbc(
    double*         ansvec,
    BBFE_BC*      	bc,
    const int 		total_num_nodes,
    const int       dof)
{
    for(int i = 0; i < total_num_nodes * dof; i++) {
        if( bc->D_bc_exists[i] ) {
            ansvec[i] = bc->imposed_D_val[i];
        }
    }
}

void BBFE_fluid_update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] += delta_v[i][j];
        }
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] += delta_p[i];
    }
}

double BBFE_fluid_calc_internal_norm_2d(
        double**  in,       //input
        const int num1,
        const int num2)
{
    double norm = 0.0;

    for (int i = 0; i < num1; i++) {
        for(int j = 0; j < num2; j++){
            norm += in[i][j] * in[i][j];
        }
    }

    return norm;
}

double BBFE_fluid_calc_internal_norm_1d(
        double*  in,       //input
        const int num1,
        const int num2)
{
    double norm = 0.0;

    for (int i = 0; i < num1; i++) {    
        norm += in[i] * in[i];
    }

    return norm;
}

void BBFE_fluid_vec_copy_2d(
        double**  in,   //input
        double**  out,  //output
        const int num1,
        const int num2)
{
    for (int i = 0; i < num1; i++) {
        for(int j = 0; j < num2; j++){
            out[i][j] = in[i][j];
        }
    }
}

void BBFE_fluid_vec_copy_1d(
    double* dst,
    const double* src,
    int size)
{
    for (int i = 0; i < size; i++) {
        dst[i] = src[i];
    }
}


double BBFE_fluid_calc_internal_norm_total(
    const double* rvec,
    int n_internal_vertex)
{
    double norm = 0.0;

    for (int i = 0; i < n_internal_vertex; i++) {
        for (int d = 0; d < 4; d++) {
            const int index = 4 * i + d;
            norm += rvec[index] * rvec[index];
        }
    }

    return norm;
}


double BBFE_fluid_calc_internal_norm_momentum(
    const double* rvec,
    int n_internal_vertex)
{
    double norm = 0.0;

    for (int i = 0; i < n_internal_vertex; i++) {
        for (int d = 0; d < 3; d++) {
            const int index = 4 * i + d;
            norm += rvec[index] * rvec[index];
        }
    }

    return norm;
}


double BBFE_fluid_calc_internal_norm_mass(
    const double* rvec,
    int n_internal_vertex)
{
    double norm = 0.0;

    for (int i = 0; i < n_internal_vertex; i++) {
        const int index = 4 * i + 3;
        norm += rvec[index] * rvec[index];
    }

    return norm;
}