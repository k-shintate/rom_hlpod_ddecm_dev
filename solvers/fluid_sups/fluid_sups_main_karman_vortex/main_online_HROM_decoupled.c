
#include "core_ROM.h"
#include "core_HROM.h"

#include "core_HROM_decoupled.h"

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";
static const char* OPTION_HOT_START     = "-hs";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

const int BUFFER_SIZE = 1024;

void ROM_read_args(
    int 		argc,
    char* 		argv[],
    ROM_PRM*    rom_prm)
{
	int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_1STDD);
    if(num == -1) {
		printf("\nargs error num_subdomains");
		exit(1);
    }
    else {
        rom_prm->num_subdomains = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_MODES);
    if(num == -1) {
		printf("\nargs error num_modes");
		exit(1);
    }
    else {
		rom_prm->num_modes = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_PADAPTIVE);
    if(num == -1) {
		printf("\nargs error rom_epsilon");
		exit(1);
    }
    else {
        rom_prm->rom_epsilon = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_SOLVER_TYPE);
    if(num == -1) {
		printf("\nargs error solver_type");
		exit(1);
    }
    else {
        rom_prm->solver_type = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_HOT_START);
    if(num == -1) {
        rom_prm->hot_start = 0;
    }
    else {
        rom_prm->hot_start = 1;
        printf("Hot start is enabled.\n");
    }

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);
    printf("hot_start = %d\n", rom_prm->hot_start);

}

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

double hot_start_read_initialize_val2(
    double*     int_val,
    const int   num_2nd_subdomains,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*         hlpod_meta,
    const char* directory)
{
    int BUFFER_SIZE = 1024;
    int total_num_nodes;
    int ndof;
    double t = 0.0;
    int index = 0;
        FILE* fp;
        char fname[BUFFER_SIZE];
        char id[BUFFER_SIZE];

        fp = BBFE_sys_read_fopen(fp, "hot_start/start_time.dat", directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%lf", &(t));
    fclose(fp);


for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        int n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];
        snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat.%d", subdomain_id);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes), &(ndof));
    for(int i = 0; i < n_internal_vertex; i++) {
        for(int j = 0; j < ndof; j++) {
            fscanf(fp, "%lf", &(int_val[i * ndof + j + index]));
        }
    }

    index += n_internal_vertex*4;
        fclose(fp);
}

    return t;
}

int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();

    double t1 = monolis_get_time();

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
    
    //memory_allocation_nodal_values(&(sys.vals_rom), sys.fe.total_num_nodes);        //set vals
    //memory_allocation_nodal_values(&(sys.vals_hrom), sys.fe.total_num_nodes);        //set vals

    memory_allocation_nodal_values(
            &(sys.vals),
            sys.fe.total_num_nodes);

    memory_allocation_nodal_values(
            &(sys.vals_rom),
            sys.fe.total_num_nodes);

    memory_allocation_nodal_values(
            &(sys.vals_hrom),
            sys.fe.total_num_nodes);

	
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);    

        BBFE_fluid_sups_read_Dirichlet_bc(
                        &(sys.bc),
                        //&(sys.vals),
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
        
    BBFE_fluid_sups_renew_velocity(
                            sys.vals_hrom.v,
                            vec,
                            sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
                            sys.vals_hrom.p,
                            vec,
                            sys.fe.total_num_nodes);

    ROM_BB_vec_copy_2d(
                        sys.vals_hrom.v,
                        sys.vals_hrom.v_old,
                        sys.fe.total_num_nodes,
                        3);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

    /*for ROM *****************************************/
	
    /*for ROM input data*/
	ROM_read_args(argc, argv, &(sys.rom_prm_p));
	ROM_read_args(argc, argv, &(sys.rom_prm_v));

    ROM_set_param(
            &(sys.rom_p),
            sys.rom_prm_p.num_subdomains,
            sys.rom_prm_p.num_modes,
            sys.rom_prm_p.rom_epsilon,
            sys.rom_prm_p.solver_type);
    
    ROM_set_param(
            &(sys.rom_v),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);
    
    ROM_set_param(
            &(sys.rom_sups),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_p.num_modes + sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            4,
            &(sys.rom_sups.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm_v.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm_v.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            3,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
	ROM_std_hlpod_pre(
            &(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

	ROM_std_hlpod_pre(
            &(sys.rom_sups),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            4,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

    /******************/

	/*for online*/
    /*
    ROM_std_hlpod_read_pod_modes_diag(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.mono_com.n_internal_vertex,
		3,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);
    */

    ROM_std_hlpod_read_pod_modes_diag_decoupled(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.mono_com.n_internal_vertex,
		3,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);

	FILE* fp;
	if(monolis_mpi_get_global_my_rank() == 0){
		fp = BBFE_sys_write_fopen(fp, "l2_error_pressure.txt", sys.cond.directory);
		fclose(fp);
		fp = BBFE_sys_write_fopen(fp, "l2_error_velosity.txt", sys.cond.directory);
		fclose(fp);
	}

	if(monolis_mpi_get_global_my_rank() == 0){
        ROM_std_hlpod_write_solver_prm_fopen("fem_solver_prm", sys.cond.directory);
		ROM_std_hlpod_write_solver_prm_fopen("pod_solver_prm", sys.cond.directory);
        ROM_std_hlpod_write_solver_prm_fopen("hr_solver_prm", sys.cond.directory);
	}

    ROM_online_read_calc_conditions(&(sys.vals), sys.cond.directory);

	monolis_com_initialize_by_self(&(sys.mono_com_rom));
	monolis_com_initialize_by_self(&(sys.mono_com_rom_solv));

	ROM_std_hlpod_set_monolis_comm(
			&(sys.mono_com),
			&(sys.mono_com_rom),
			&(sys.mono_com_rom_solv),
			INPUT_DIRECTORYNAME_METAGRAPH,
            INPUT_FILENAME_METAGRAPH,
            sys.rom_prm_v.solver_type,
			sys.cond.directory);

    monolis_initialize(&(sys.monolis_rom0));
    monolis_initialize(&(sys.monolis_mass_rom0));

    ROM_std_hlpod_online_pre(
            &(sys.monolis_rom0),
            &(sys.mono_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            metagraph_name,
            sys.cond.directory);
        
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_rom0), &(sys.monolis_rom));
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_rom0), &(sys.monolis_mass_rom0));
    /********************/

    /*for Hyper-reduction*/
/*
    HROM_memory_allocation_online(&sys, &(sys.rom_sups), &(sys.hrom_sups));
    HROM_pre(&sys, &(sys.rom_sups), &(sys.hrom_sups));
    HROM_pre_online(&sys, &(sys.rom_sups), &(sys.hrom_sups));

    HROM_set_bc_id(
        &(sys.fe),
        (&sys.bc),
        &(sys.hrom_sups.hlpod_ddhr),
        4,
        &(sys.rom_sups.hlpod_mat));
*/
        
    HROM_memory_allocation_online(&sys, &(sys.rom_sups), &(sys.hrom_sups));
    HROM_pre_decoupled(&sys, &(sys.rom_sups), &(sys.hrom_sups));
    HROM_pre_online_decoupled2(&sys, &(sys.rom_sups), &(sys.hrom_sups));

    HROM_set_bc_id(
        &(sys.fe),
        (&sys.bc),
        &(sys.hrom_sups.hlpod_ddhr),
        4,
        &(sys.rom_sups.hlpod_mat));
    
    HROM_ddecm_memory_allocation_sol_vec(&(sys.hrom_sups.hr_vals), sys.fe.total_num_nodes, 4);
   	monolis_initialize(&(sys.monolis_hr));
    monolis_copy_mat_R(&(sys.monolis_hr0), &(sys.monolis_hr));
    /************************/


    /**************************************************/

    /*for rom**/
    read_calc_conditions(&(sys.vals_rom), sys.cond.directory);                      //set vals
    ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom_sups.hlpod_vals), sys.fe.total_num_nodes, 4);

	set_target_parameter(&(sys.vals), sys.cond.directory);
	set_target_parameter(&(sys.vals_rom), sys.cond.directory);
    /**********/

    /*for hrom*/
    read_calc_conditions(&(sys.vals_hrom), sys.cond.directory);                      //set vals
	set_target_parameter(&(sys.vals_hrom), sys.cond.directory);
    /******** */

    solver_rom_NR4(&(sys), 0, 0, 0);
    double tt2 = monolis_get_time_global_sync();
	printf("done solver ROM");
    add_reduced_mat_linear(&(sys), 0, 0, 0);
    double tt3 = monolis_get_time_global_sync();
    // printf("done add_mat_linear");

    ROM_BB_vec_copy_2d(
                        sys.vals.v,
                        sys.vals_rom.v,
                        sys.fe.total_num_nodes,
                        3);

    ROM_BB_vec_copy_2d(
                        sys.vals.v,
                        sys.vals_rom.v_old,
                        sys.fe.total_num_nodes,
                        3);

    ROM_BB_vec_copy_2d(
                        sys.vals.v,
                        sys.vals_hrom.v,
                        sys.fe.total_num_nodes,
                        3);

    ROM_BB_vec_copy_2d(
                        sys.vals.v,
                        sys.vals_hrom.v_old,
                        sys.fe.total_num_nodes,
                        3);

	double FOM_t2 = monolis_get_time();
	double ROM_t1 = monolis_get_time();

	double t = 0.0;
	int file_num = 0;
	int step_rom = 0;
    double t_hs = 0.0;
    int step_hs = 0;
/*
	while (t < sys.vals.rom_finish_time) {
		t += sys.vals.dt;
		step_rom += 1;
		*/
/*
    	if(sys.rom_prm_p.hot_start == 1){
            char fname[BUFFER_SIZE];         
            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
            double* val = BB_std_calloc_1d_double(val, 4*sys.fe.total_num_nodes);
            t_hs = hot_start_read_initialize_val(val, fname, sys.cond.directory);
            step_hs = 0;

            printf("Hot start time: %lf\n", t);
            printf("Hot start step: %d\n", step_rom);
            printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        }
*/

        if(sys.rom_prm_p.hot_start == 1) {
            char fname[BUFFER_SIZE];
            snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat.%d", sys.rom_p.hlpod_meta.subdomain_id[0]);

            double* val = BB_std_calloc_1d_double(val, 4 * sys.fe.total_num_nodes);

            t_hs = hot_start_read_initialize_val2(
                    val,
                    sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                    &(sys.rom_sups.hlpod_mat),
                    &(sys.rom_sups.hlpod_meta),
                    sys.cond.directory);

            step_hs = 0;

            monolis_mpi_update_R(&(sys.mono_com), sys.fe.total_num_nodes, 4, val);

            printf("Hot start time: %lf\n", t);
            printf("Hot start step: %d\n", step_rom);
            printf("sys.vals.finish_time - t = %lf\n",
                    ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            BBFE_fluid_sups_renew_velocity(sys.vals_rom.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals_rom.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals_rom.p, val, sys.fe.total_num_nodes);

	        BBFE_fluid_sups_renew_velocity(sys.vals_hrom.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals_hrom.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals_hrom.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4 * sys.fe.total_num_nodes);
        }

    while (t < sys.vals.rom_finish_time) {
                t += sys.vals.dt;
                step_rom += 1;
		printf("\n%s ----------------- step-ROM %d ----------------\n", CODENAME, step_rom);

		/***************FEM***************/
		
		printf("----------------- normal-FEM ----------------\n");
		double calctime_fem_t2 = monolis_get_time_global_sync();
		//solver_fom(sys, t, step_rom);
        solver_fom_NR(sys, t, step_rom);
		double calctime_fem_t1 = monolis_get_time_global_sync();

        /**********************************/
		
		/****************ROM***************/
		printf("----------------- ROM ----------------\n");

		double calctime_rom_t2 = monolis_get_time_global_sync();
		if(sys.rom_sups.hlpod_vals.bool_global_mode==false){
			//solver_rom(&(sys), step_rom, 0, t);
            //solver_rom_NR5(&(sys), t, step_rom, 0);
		}
		else{
			solver_rom_global_para(
						&(sys.monolis),
						&(sys.mono_com),
						&(sys.rom_sups),
						&(sys),
						step_rom,
						0,
						t);
		}
		double calctime_rom_t1 = monolis_get_time_global_sync();

		/**********************************/

		double calctime_hrom_t2 = monolis_get_time_global_sync();

		if(monolis_mpi_get_global_comm_size() == 1){
			//HROM_nonparallel(sys, &(sys.rom_sups), &(sys.hrom_sups), step_rom, 0, t);
        }
        else{
			//HROM_hierarchical_parallel(sys, &(sys.rom_sups), &(sys.hrom_sups), step_rom, 0, t);
            //solver_hrom_NR(&(sys), t, step_rom, 0);
            solver_hrom_NR_decoupled(&(sys), t, step_rom, 0);
        }

		double calctime_hrom_t1 = monolis_get_time_global_sync();

		if(step_rom%sys.vals.output_interval == 0) {
			HROM_output_files(&sys, file_num, t);
            
	    if(monolis_mpi_get_global_my_rank() == 0){
		
            ROM_std_hlpod_write_solver_prm(&(sys.monolis), t, "fem_solver_prm/" , sys.cond.directory);
			ROM_std_hlpod_write_solver_prm(&(sys.monolis_rom), t, "pod_solver_prm/", sys.cond.directory);
            ROM_std_hlpod_write_solver_prm(&(sys.monolis_rom), t, "hr_solver_prm/", sys.cond.directory);

            ROM_std_hlpod_output_calc_time(calctime_fem_t1-calctime_fem_t2, t,
					"calctime/time_fem.txt", sys.cond.directory);
            ROM_std_hlpod_output_calc_time(calctime_rom_t1-calctime_rom_t2, t,
					"calctime/time_rom.txt", sys.cond.directory);
            ROM_std_hlpod_output_calc_time(calctime_hrom_t1-calctime_hrom_t2, t,
					"calctime/time_hrom.txt", sys.cond.directory);
	    }
			file_num += 1;
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


