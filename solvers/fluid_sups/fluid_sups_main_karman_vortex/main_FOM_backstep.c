#include "core_ROM.h"

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

	memory_allocation_nodal_values_AB2(
			&(sys.vals),
			sys.fe.total_num_nodes);
	
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
			&(sys.bc),
			&(sys.vals),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			4);

    BBFE_fluid_sups_read_Dirichlet_bc_NR(
            &(sys.bc_NR),
            filename,
            sys.cond.directory,
            sys.fe.total_num_nodes,
            4);

/*
    pre_surface(
            //&(sys.monolis_com_surf),
            &(sys.surf),
            &(sys.basis_surf),
            sys.cond.directory,
            "surf_graph.dat",
            sys.vals.num_ip_each_axis);
*/
    FILE* fp;
    if(monolis_mpi_get_global_my_rank() == 0){
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_force.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_force.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinderr_lift_coeff.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "out_flow.txt", sys.cond.directory);
            fclose(fp);
    }

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

	//intialize for velocity and pressure
	initialize_velocity_pressure_karman_vortex(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);

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

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	int count = 0;  //for ROM

	t = 0.0; step = 0; file_num = 0;

    while (t < sys.vals.finish_time) {
        t += sys.vals.dt;
        step += 1;

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
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            //BBFE_fluid_sups_renew_velocity(sys.vals_rom.v, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_pressure(sys.vals_rom.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        //}

        printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
        solver_fom_NR(sys, t, count);
        count ++;

        if(step%sys.vals.output_interval == 0) {
            output_files(&sys, step, t);
            file_num += 1;
        }
        
        if(step%10 == 0) {
            double D_out = 0.0; double L_out = 0.0; double Cd_out = 0.0; double Cl_out = 0.0;
            const double eU[3] = {1.0, 0.0, 0.0};
            const double eP[3] = {0.0, 1.0, 0.0};

            calc_Cd_v(&(sys.surf), &(sys.fe), &(sys.basis_surf), &(sys.mono_com), &(sys.vals), sys.vals.density,
                sys.vals.viscosity, 1, 0.08*1, eU, eP, &D_out, &L_out, &Cd_out, &Cl_out);

            Cd_out = 0.0;
            Cl_out = 0.0;

            calc_Cd_p(&(sys.surf), &(sys.fe), &(sys.basis_surf), &(sys.vals), eU, eP, &Cd_out, &Cl_out);

            monolis_allreduce_R(1, &D_out, MONOLIS_MPI_SUM, sys.mono_com.comm);
            monolis_allreduce_R(1, &L_out, MONOLIS_MPI_SUM, sys.mono_com.comm);
            monolis_allreduce_R(1, &Cd_out, MONOLIS_MPI_SUM, sys.mono_com.comm);
            monolis_allreduce_R(1, &Cl_out, MONOLIS_MPI_SUM, sys.mono_com.comm);

            ROM_std_hlpod_output_add_calc_time(D_out, t,
                            "cylinder_drag_force.txt", sys.cond.directory);
            ROM_std_hlpod_output_add_calc_time(L_out, t,
                            "cylinder_lift_force.txt", sys.cond.directory);
            ROM_std_hlpod_output_add_calc_time(Cd_out, t,
                            "cylinder_drag_coeff.txt", sys.cond.directory);
            ROM_std_hlpod_output_add_calc_time(Cl_out, t,
                            "cylinder_lift_coeff.txt", sys.cond.directory);
        }
                            

        BBFE_fluid_sups_add_velocity_pressure(
                        sys.vals.v,
                        sys.vals.p,
                        sys.monolis.mat.R.X,
                        sys.fe.total_num_nodes);
        
        if(step%1000 == 0){
            output_result_file_karman_vortex_pressure(
                &(sys.fe),
                &(sys.vals),
                t,
                sys.cond.directory);
        }

        if(step%1 == 0){
                char fname[BUFFER_SIZE];
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
                hot_start_write_initialize_val(sys.monolis.mat.R.X, sys.fe.total_num_nodes, 4, t, fname, sys.cond.directory);
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
