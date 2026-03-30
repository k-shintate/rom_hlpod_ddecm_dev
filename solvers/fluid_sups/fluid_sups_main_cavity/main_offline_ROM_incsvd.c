
#include "core_ROM.h"
#include "core_HROM.h"

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";

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

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_fluid_sups_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			4);

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

    monolis_initialize(&(sys.monolis_rom0));

    HROM_std_hlpod_online_pre(
            &(sys.monolis_rom0),
            &(sys.mono_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            4,
            metagraph_name,
            sys.cond.directory);
        
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_rom0), &(sys.monolis_rom));
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);
    /********************/


    /*for ROM set vals***/
    read_calc_conditions(&(sys.vals_rom), sys.cond.directory);                      //set vals
    memory_allocation_nodal_values(&(sys.vals_rom), sys.fe.total_num_nodes);        //set vals

	initialize_velocity_pressure(sys.vals_rom.v, sys.vals_rom.p, sys.fe.total_num_nodes);
	initialize_velocity_pressure(sys.vals_rom.v, sys.vals_rom.p, sys.fe.total_num_nodes);

    ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom_sups.hlpod_vals), sys.fe.total_num_nodes, 4);

	set_target_parameter(&(sys.vals), sys.cond.directory);
	set_target_parameter(&(sys.vals_rom), sys.cond.directory);
    ROM_offline_set_reynolds_num_cases(&(sys.vals_rom), sys.cond.directory);
    ROM_offline_set_reynolds_num_cases(&(sys.vals), sys.cond.directory);
    /*********************/

    /*for Hyper-reduction*/
    HROM_offline_read_calc_conditions_inc_svd(&(sys.hrom_sups.hr_vals), sys.cond.directory);

	HROM_std_hlpod_offline_set_num_snapmat_inc_svd(
			&(sys.rom_sups),
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals_rom.num_cases,
            sys.hrom_sups.hr_vals.incsvd_interval);

    HROM_pre(&sys, &(sys.rom_sups), &(sys.hrom_sups));
    HROM_memory_allocation(&sys, &(sys.rom_sups), &(sys.hrom_sups));

    HROM_set_bc_id(
        &(sys.fe),
        (&sys.bc),
        &(sys.hrom_sups.hlpod_ddhr),
        4,
        &(sys.rom_sups.hlpod_mat));
    /************************/

	double FOM_t2 = monolis_get_time();
	double ROM_t1 = monolis_get_time();

	double t = 0.0;
	int file_num = 0;
	int step_rom = 0;
    int step_inc_svd = 0;
    int num_inc_svd = 0;
    int count = 0;

	for(int i = 0; i < sys.vals.num_cases; i++){
		ROM_offline_set_reynolds_number(&(sys.vals), i);
        ROM_offline_set_reynolds_number(&(sys.vals_rom), i);
	
		t = 0.0; step_rom = 0; step_inc_svd = 0; file_num = 0; count = 0;
        
        while (t < sys.vals.rom_finish_time) {
            t += sys.vals.dt;
            step_rom += 1;
            step_inc_svd += 1;

            printf("\n%s ----------------- step-ROM %d ----------------\n", CODENAME, step_rom);

            /***************FEM***************/
            
            printf("----------------- normal-FEM ----------------\n");
            /**********************************/
            
            /****************ROM***************/
            printf("----------------- ROM ----------------\n");

            double calctime_rom_t2 = monolis_get_time();
            if(sys.rom_sups.hlpod_vals.bool_global_mode==false){
                solver_rom(&(sys), step_rom, 0, t);
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
            double calctime_rom_t1 = monolis_get_time();

            /**********************************/

            double calctime_hr_t2 = monolis_get_time();
            HROM_set_matvec(&(sys),&(sys.rom_sups),&(sys.hrom_sups),step_inc_svd,t);
            double calctime_hr_t1 = monolis_get_time();

            if(step_rom%sys.hrom_sups.hr_vals.incsvd_interval == 0 && num_inc_svd==0){
                HROM_pre_offline_inc_svd1(&sys, &(sys.rom_sups), &(sys.hrom_sups), 1,1,1);
                step_inc_svd = 0;
                num_inc_svd ++;
            }
            else if(step_rom%sys.hrom_sups.hr_vals.incsvd_interval == 0 && num_inc_svd != 0){
                HROM_pre_offline_inc_svd2(&sys, &(sys.rom_sups), &(sys.hrom_sups), 1,1,1);
                step_inc_svd = 0;
                num_inc_svd ++;
            }

            double calctime_fem_t2 = monolis_get_time();
            solver_fom(sys, t, step_rom);	
            double calctime_fem_t1 = monolis_get_time();

            if(step_rom%sys.vals.output_interval == 0) {
                ROM_output_files(&sys, file_num, t);
                            
                ROM_std_hlpod_write_solver_prm(&(sys.monolis), t, "fem_solver_prm/" , sys.cond.directory);
                ROM_std_hlpod_write_solver_prm(&(sys.monolis_rom), t, "pod_solver_prm/", sys.cond.directory);

                ROM_std_hlpod_output_calc_time(calctime_fem_t2-calctime_fem_t1, t,
                        "calctime/time_fem.txt", sys.cond.directory);
                ROM_std_hlpod_output_calc_time(calctime_rom_t2-calctime_rom_t1, t,
                        "calctime/time_rom.txt", sys.cond.directory);

                file_num += 1;
            }
        }
    }

    HROM_pre_offline_inc_svd3(&sys, &(sys.rom_sups), &(sys.hrom_sups), 1,1,1);

	BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
	BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), sys.fe.total_num_nodes, 4);
	monolis_finalize(&(sys.monolis));

	double t2 = monolis_get_time();
	int myrank = monolis_mpi_get_global_my_rank();

	if(myrank == 0) {
		printf("** Total time: %f\n", t2 - t1);
        printf("** Number of Incremental SVD : %d\n", num_inc_svd);
	}

	monolis_global_finalize();

	printf("\n");

	return 0;
}

