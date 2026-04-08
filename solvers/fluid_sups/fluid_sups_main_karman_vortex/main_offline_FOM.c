
#include "core_ROM.h"

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
	double FOM_t1 = monolis_get_time();

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
	BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
			&(sys.bc),
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

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

	//intialize for velocity and pressure
	initialize_velocity_pressure_karman_vortex(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);
    initialize_velocity_pressure_cavity(sys.vals.v, sys.vals.v_old, sys.vals.p, sys.fe.total_num_nodes);

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

    /*for offline******/
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);
	ROM_offline_set_reynolds_num_cases(&(sys.vals), sys.cond.directory);

    /*for hot start*/
    ROM_offline_set_reynolds_number(&(sys.vals), 0);
    double t_hotstart = 0.0;
  /*
    //if(sys.rom_prm_p.hot_start == 1){
        char fname[BUFFER_SIZE];         
        //snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
        snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat");
        double* val = BB_std_calloc_1d_double(val, 4*sys.fe.total_num_nodes);
        t_hotstart = hot_start_read_initialize_val(val, fname, sys.cond.directory);
        BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        printf("Hot start time: %lf\n", t_hotstart);
    //}
    //else{
        //t_hotstart = 0.0;
    //}
*/
//	if(sys.rom_prm_p.hot_start == 1) {
            char fname[BUFFER_SIZE];
            //snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat.%d", sys.rom_p.hlpod_meta.subdomain_id[0]);

            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());

            double* val = BB_std_calloc_1d_double(val, 4 * sys.fe.total_num_nodes);

            t_hotstart = hot_start_read_initialize_val2(
                    val,
                    sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                    &(sys.rom_sups.hlpod_mat),
                    &(sys.rom_sups.hlpod_meta),
                    sys.cond.directory);

                BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
                printf("Hot start time: %lf\n", t_hotstart);
/*
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
            */
  //      }

    printf("sys.vals.finish_time - t_hotstart = %lf\n", ((double)sys.vals.finish_time - t_hotstart));
    /***************/

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hotstart),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			3);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hotstart),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			1);

    /******************/

	/**********************************************/

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	int count = 0;  //for ROM
    double t_hs = 0.0; //for hot start
    int step_hs = 0; //for hot start

	for(int i = 0; i < sys.vals.num_cases; i++){
		ROM_offline_set_reynolds_number(&(sys.vals), i);
		
		t = 0.0; step = 0; file_num = 0;
/*
		//if(sys.rom_prm_p.hot_start == 1){
            char fname[BUFFER_SIZE];         
            //snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
            snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat");
            double* val = BB_std_calloc_1d_double(val, 4*sys.fe.total_num_nodes);
            t = hot_start_read_initialize_val(val, fname, sys.cond.directory);
            step_hs = 0;

            printf("Hot start time: %lf\n", t);
            printf("Hot start step: %d\n", step);
            printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        //}

*/

//		if(sys.rom_prm_p.hot_start == 1) {
            char fname[BUFFER_SIZE];
            //snprintf(fname, BUFFER_SIZE, "hot_start/hot_start.dat.%d",
            //        sys.rom_p.hlpod_meta.subdomain_id[0]);
            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());

            double* val = BB_std_calloc_1d_double(val, 4 * sys.fe.total_num_nodes);

            t_hotstart = hot_start_read_initialize_val2(
                    val,
                    sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                    &(sys.rom_sups.hlpod_mat),
                    &(sys.rom_sups.hlpod_meta),
                    sys.cond.directory);

                //BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
                //printf("Hot start time: %lf\n", t_hotstart);
            //step_hs = 0;


            monolis_mpi_update_R(&(sys.mono_com), sys.fe.total_num_nodes, 4, val);

            printf("Hot start time: %lf\n", t);
            //printf("Hot start step: %d\n", step_rom);
            printf("sys.vals.finish_time - t = %lf\n",
                    ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            //BBFE_fluid_sups_renew_velocity(sys.vals_rom.v, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_velocity(sys.vals_rom.v_old, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_pressure(sys.vals_rom.p, val, sys.fe.total_num_nodes);

            //BBFE_fluid_sups_renew_velocity(sys.vals_hrom.v, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_velocity(sys.vals_hrom.v_old, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_pressure(sys.vals_hrom.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4 * sys.fe.total_num_nodes);

  //      }


		while (t < sys.vals.finish_time - t_hotstart) {
			t += sys.vals.dt;
			step += 1;

			printf("\n%s ----------------- step %d ----------------\n", CODENAME, step + step_hs);

			//if(step > (sys.rom_v.hlpod_vals.num_snapshot * sys.vals.snapshot_interval)){
			//	break;
			//}
			//solver_fom_collect_snapmat(sys, t, count);
            //solver_fom_NR(sys, t, count);
            solver_fom_NR_collect_snapmat(sys, t, count);
			count ++;

			//if(step%sys.vals.output_interval == 0) {
            if(step%10 == 0) {
				output_files(&sys, file_num, t);
				file_num += 1;
			}

            /*
            BBFE_fluid_sups_update_vec(
        		sys.vals.v,
                sys.vals.p,
		        sys.monolis.mat.R.X,
		        sys.fe.total_num_nodes);

            char fname[BUFFER_SIZE];         
            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
            hot_start_write_initialize_val(sys.monolis.mat.R.X, sys.fe.total_num_nodes, 4, t, fname, sys.cond.directory);
            */
            
            /*
            if(step == 40000 && sys.rom_prm_p.hot_start == 0){
                char fname[BUFFER_SIZE];         
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
                hot_start_write_initialize_val(sys.monolis.mat.R.X, sys.fe.total_num_nodes, 4, t, fname, sys.cond.directory);
            }

            if(step%10 == 0 && t > 100) {
				output_result_file_karman_vortex(&(sys.fe), &(sys.vals), t, sys.cond.directory);
    			BBFE_fluid_sups_renew_pressure(sys.vals.p, sys.monolis.mat.R.X, sys.fe.total_num_nodes);
				output_result_file_karman_vortex_pressure(&(sys.fe), &(sys.vals), t, sys.cond.directory);
				output_result_file_karman_vortex_pressure_inf(&(sys.fe), &(sys.vals), t, sys.cond.directory);

			}
            */
		}
	}


    /**********************************************/

	ROM_std_hlpod_set_pod_modes_diag(
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

    /*for writing vtk*/
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

	ROM_sys_hlpod_fe_write_pod_modes_vtk_diag(
		&(sys.fe),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		10,
		10,
		3,
		1,
		"pod_modes_vtk/pod_modes_v.vtk",
		"pod_modes_vtk/pod_modes_p.vtk",
		sys.cond.directory);
    /***************/

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

