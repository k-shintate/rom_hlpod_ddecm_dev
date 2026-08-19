
#include "core_ROM.h"

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";


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

static void check_edge_heat_bc(
    const BBFE_DATA* fe,
    const BBFE_BC* bc)
{
    double xmin = fe->x[0][0];
    double xmax = fe->x[0][0];

    for(int i = 1; i < fe->total_num_nodes; i++) {
        if(fe->x[i][0] < xmin) xmin = fe->x[i][0];
        if(fe->x[i][0] > xmax) xmax = fe->x[i][0];
    }

    const double tol =
        1.0e-10 * fmax(1.0, fabs(xmax - xmin));

    int n_bc = 0;
    int n_wrong = 0;

    for(int i = 0; i < fe->total_num_nodes; i++) {

        if(!bc->D_bc_exists[i]) {
            continue;
        }

        n_bc++;

        if(fabs(fe->x[i][0] - xmin) > tol) {
            fprintf(
                stderr,
                "[EDGE-HEAT-BC-ERROR] "
                "node=%d x=(%.16e %.16e %.16e) xmin=%.16e\n",
                i,
                fe->x[i][0],
                fe->x[i][1],
                fe->x[i][2],
                xmin);

            n_wrong++;
        }
    }

    fprintf(
        stderr,
        "[EDGE-HEAT-BC-CHECK] xmin=%.16e "
        "num_D_bc=%d wrong_D_bc=%d\n",
        xmin,
        n_bc,
        n_wrong);

    if(n_wrong > 0) {
        fprintf(
            stderr,
            "ERROR: edge-heating benchmark requires "
            "Dirichlet BC only on x=xmin\n");
        exit(EXIT_FAILURE);
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

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

    check_edge_heat_bc(
        &(sys.fe),
        &(sys.bc));

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	manusol_set_init_value(&(sys.fe), sys.vals.T);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

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


    /*for ROM ***************************************/
	
    /*for ROM input data*/
	ROM_read_args(argc, argv, &(sys.rom_prm));

    ROM_set_param(
            &(sys.rom),
            sys.rom_prm.num_subdomains,
            sys.rom_prm.num_modes,
            sys.rom_prm.rom_epsilon,
            sys.rom_prm.solver_type);
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            1,
            &(sys.rom.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
    /******************/

    /*for offline phase*/
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            5,
			1);
    /******************/

    /**********************************************/

    monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;

    double train_Din[] = {
        1.0e-3,
        1.0e-4,
        1.0e-5,
        1.0e-6,
        1.0e-7
    };

    int num_train_param = sizeof(train_Din) / sizeof(train_Din[0]);

    for(int iparam = 0; iparam < num_train_param; iparam++)
    {
        MANUSOL_PARAM prm = {
            .radius = 0.75,
            .D_out  = 1.0,
            .D_in   = train_Din[iparam],
            .T0     = 0.0,
            .Q0     = 1.0,
            .tau    = 1.0
        };

        manusol_set_param(&prm);

        /*
        * D_IN が変化するので
        * A(mu) を作り直す
        */
        monolis_clear_mat_value_R(
            &(sys.monolis0));

        set_element_mat(
            &(sys.monolis0),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
        
        monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));

        /*
        * 各 parameter case の初期状態
        */
        manusol_set_init_value(
            &(sys.fe),
            sys.vals.T);

        double t = 0.0;
        int step = 0;

        while(t < sys.vals.finish_time) {

            t += sys.vals.dt;
            step++;

            solver_fom_collect_snapmat(
                &sys,
                t,
                step,
                iparam);

            if(step%sys.vals.output_interval == 0) {
                output_files(&sys, file_num, t);
                        
                file_num += 1;
            }
        }

    }
    /**************************************************/

    /*for ROM *****************************************/
    ROM_std_hlpod_set_pod_modes(
		&(sys.rom),
		sys.fe.total_num_nodes,
		sys.monolis_com.n_internal_vertex,
		1,
		"pod_modes",
		sys.cond.directory);
    
    ROM_sys_hlpod_fe_write_pod_modes_vtk(
		&(sys.fe),
		&(sys.rom),
		sys.fe.total_num_nodes,
		10,
		1,
		"pod_modes_vtk/pod_modes.vtk",
		sys.cond.directory);
    /***************************************************/

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis0));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
