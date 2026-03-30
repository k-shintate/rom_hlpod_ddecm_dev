
#include "core_ROM.h"

static const char* OPTION_NUM_MODES = "-nm";
static const char* OPTION_NUM_1STDD = "-nd";
static const char* OPTION_PADAPTIVE = "-pa";
static const char* OPTION_SOLVER_TYPE = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";


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

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);

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
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_P);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc_p),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			1);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc_v),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			3);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.mono_pred) , &(sys.mono_com), &(sys.fe), 3, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_ppe)  , &(sys.mono_com), &(sys.fe), 1, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_ppe0) , &(sys.mono_com), &(sys.fe), 1, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_corr) , &(sys.mono_com), &(sys.fe), 3, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_corr0), &(sys.mono_com), &(sys.fe), 3, sys.cond.directory);

	BBFE_elemmat_set_global_mat_Laplacian_const(
			&(sys.mono_ppe0),  &(sys.fe), &(sys.basis), 1.0);
	BBFE_elemmat_set_global_mat_cmass_const(
			&(sys.mono_corr0), &(sys.fe), &(sys.basis), 1.0, 3);


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
        
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc_p),
            sys.fe.total_num_nodes,
            1,
            &(sys.rom_p.rom_bc));
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc_v),
            sys.fe.total_num_nodes,
            3,
            &(sys.rom_v.rom_bc));

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
    /******************/

	/*only for online*/
	ROM_std_hlpod_read_pod_modes(
            &(sys.rom_p),
            sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            1,
            "pod_modes_p",
            sys.cond.directory);

	ROM_std_hlpod_read_pod_modes(
            &(sys.rom_v),
            sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            3,
            "pod_modes_v",
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

	monolis_initialize(&(sys.mono_ppe_rom0));
    monolis_initialize(&(sys.mono_corr_rom0));

	ROM_std_hlpod_online_pre(
            &(sys.mono_ppe_rom0),
            &(sys.mono_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_p),
            sys.fe.total_num_nodes,
            metagraph_name,
            sys.cond.directory);

	ROM_std_hlpod_online_pre(
            &(sys.mono_corr_rom0),
            &(sys.mono_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_v),
            sys.fe.total_num_nodes,
            metagraph_name,
            sys.cond.directory);

    monolis_initialize(&(sys.mono_pred_rom));    
    monolis_copy_mat_nonzero_pattern_R(&(sys.mono_ppe_rom0), &(sys.mono_ppe_rom));
    monolis_copy_mat_nonzero_pattern_R(&(sys.mono_corr_rom0), &(sys.mono_corr_rom));
    monolis_copy_mat_nonzero_pattern_R(&(sys.mono_corr_rom0) , &(sys.mono_pred_rom));
    monolis_com_initialize_by_self(&(sys.mono_com0));
    /********************/

    /**************************************************/

    read_calc_conditions(&(sys.vals_rom), sys.cond.directory);                      //set vals
    memory_allocation_nodal_values(&(sys.vals_rom), sys.fe.total_num_nodes);        //set vals

	initialize_velocity_pressure(sys.vals_rom.v, sys.vals_rom.p, sys.fe.total_num_nodes);
	initialize_velocity_pressure(sys.vals_rom.v, sys.vals_rom.p, sys.fe.total_num_nodes);

	ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom_p.hlpod_vals), sys.fe.total_num_nodes, 1);
    ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom_v.hlpod_vals), sys.fe.total_num_nodes, 3);

	set_target_parameter(&(sys.vals), sys.cond.directory);
	set_target_parameter(&(sys.vals_rom), sys.cond.directory);

	double FOM_t2 = monolis_get_time();
	double ROM_t1 = monolis_get_time();

	double t = 0.0;
	int file_num = 0;
	int step_rom = 0;

	while (t < sys.vals.rom_finish_time) {
		t += sys.vals.dt;
		step_rom += 1;

		printf("\n%s ----------------- step-ROM %d ----------------\n", CODENAME, step_rom);

		/***********normal FEM***********/
		
		printf("----------------- normal-FEM ----------------\n");
		double calctime_fem_t2 = monolis_get_time();
		solver_fom(sys, t, step_rom);	
		double calctime_fem_t1 = monolis_get_time();
		
		/**************for POD*************/
		printf("----------------- ROM ----------------\n");

		double calctime_rom_t2 = monolis_get_time();
		if(sys.rom_v.hlpod_vals.bool_global_mode==false){
			solver_rom(&(sys), step_rom, t);
		}
		else{/*
			solver_rom_global_para(
						&(sys.monolis),
						&(sys.mono_com),
						&(sys.rom_sups),
						&(sys),
						step_rom,
						0,
						t);
                        */
		}
		double calctime_rom_t1 = monolis_get_time();

		/**********************************/
		double calctime_hr_t2 = monolis_get_time();
		double calctime_hr_t1 = monolis_get_time();

		if(step_rom%sys.vals.output_interval == 0) {
			ROM_output_files(&sys, file_num, t);
            
            ROM_std_hlpod_output_calc_time(calctime_fem_t2-calctime_fem_t1, t,
					"calctime/time_fem.txt", sys.cond.directory);
            ROM_std_hlpod_output_calc_time(calctime_rom_t2-calctime_rom_t1, t,
					"calctime/time_rom.txt", sys.cond.directory);

			file_num += 1;
		}
	}

	BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
	BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc_v), sys.fe.total_num_nodes, 3);
	BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc_p), sys.fe.total_num_nodes, 1);
	monolis_finalize(&(sys.mono_pred));
	monolis_finalize(&(sys.mono_ppe));
	monolis_finalize(&(sys.mono_ppe0));
	monolis_finalize(&(sys.mono_corr));
	monolis_finalize(&(sys.mono_corr0));

	double t2 = monolis_get_time();
	int myrank = monolis_mpi_get_global_my_rank();

	if(myrank == 0) {
		printf("** Total time: %f\n", t2 - t1);
	}

	monolis_global_finalize();

	printf("\n");

	return 0;
}

