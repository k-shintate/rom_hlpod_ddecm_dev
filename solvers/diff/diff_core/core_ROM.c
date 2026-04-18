
#include "core_ROM.h"

/*for ROM*/
const char*     ROM_ID_FINISH_TIME      = "#rom_finish_time";
const double  POD_DVAL_FINISH_TIME      = 1.0;
const char* POD_ID_OUTPUT_INTERVAL      = "#rom_output_interval";
const int POD_DVAL_OUTPUT_INTERVAL      = 1;
const char* POD_ID_SNAPSHOT_INTERVAL    = "#snapshot_interval";
const int POD_DVAL_SNAPSHOT_INTERVAL    = 1;


static const char* POD_INPUT_FILENAME_COND          = "rom_cond.dat";
static const char* ROM_OUTPUT_FILENAME_VTK          = "rom_result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

static const char* INPUT_DIRECTORYNAME_1STDD        = "parted.0/";
static const char* INPUT_DIRECTORYNAME_2NDD         = "parted.1/";
static const char* INPUT_DIRECTORYNAME_METAGRAPH    = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH         = "metagraph.dat";

const int BUFFER_SIZE = 10000;


const char* ROM_std_hlpod_get_parted_file_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_2NDD);
        return fname;
    }
    else if(solver_type==4){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}

const char* ROM_std_hlpod_get_metagraph_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_1STDD, INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_METAGRAPH ,INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==4){
        return NULL;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}

void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type)
{
    rom->hlpod_vals.num_1st_subdomains = num_subdomains;
    rom->hlpod_vals.num_2nd_subdomains = num_subdomains;
    rom->hlpod_vals.num_modes_pre = num_modes;
    rom->hlpod_vals.rom_epsilon = 1 - rom_epsilon;

    rom->hlpod_vals.bool_global_mode = false;
    
    if(solver_type == 4){
        rom->hlpod_vals.bool_global_mode = true;
    }

}

void ROM_offline_assign_default_values(
		VALUES*     vals)
{
	vals->snapshot_interval  = POD_DVAL_SNAPSHOT_INTERVAL;
}


void ROM_offline_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, POD_ID_SNAPSHOT_INTERVAL,  vals->snapshot_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_offline_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_offline_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_int_p(
				&(vals->snapshot_interval), filename, POD_ID_SNAPSHOT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_offline_print_all_values(vals);

	printf("\n");
}


void ROM_online_assign_default_values(
		VALUES*     vals)
{
	vals->rom_finish_time  = POD_DVAL_FINISH_TIME;
	vals->rom_output_interval  = POD_DVAL_OUTPUT_INTERVAL;
}


void ROM_online_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %e\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_finish_time);
	printf("%s %s: %d\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_output_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_online_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_online_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_double_p(
				&(vals->rom_finish_time), filename, ROM_ID_FINISH_TIME, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->rom_output_interval), filename, POD_ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_online_print_all_values(vals);

	printf("\n");
}


void ROM_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*   hlpod_vals,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");

	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "abs_error");
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}


void ROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, ROM_OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
/*
	ROM_output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			&(sys->rom.hlpod_vals),
			filename,
			sys->cond.directory,
			t);
	
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals_rom.T,
			filename,
			sys->cond.directory);
*/
	double L2_error = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis_com),
			t,
			sys->vals_rom.T,
			sys->vals.T);

	printf("%s L2 error: %e\n", CODENAME, L2_error);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_rom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error);
		fclose(fp);
	}

}


void solver_rom(
    FE_SYSTEM* sys,
    const int step,
    const double t)
{
    printf("\n%s ----------------- step ROM %d ----------------\n", CODENAME, step);

    double t1 = monolis_get_time_global_sync();

    monolis_copy_mat_value_R(&(sys->monolis0), &(sys->monolis));
    monolis_clear_mat_value_rhs_R(&(sys->monolis));
    monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
    monolis_clear_mat_value_rhs_R(&(sys->monolis_rom));

   set_element_vec(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom),
            t);

    manusol_set_theo_sol(&(sys->fe), sys->vals.theo_sol, t);
    BBFE_manusol_set_bc_scalar(
            &(sys->fe),
            &(sys->bc),
            sys->vals.theo_sol,
            t);

    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            BLOCK_SIZE,
            &(sys->bc),
            sys->monolis.mat.R.B);

    double t2 = monolis_get_time_global_sync();

    sys->rom.hlpod_vals.time_calc_reduced_matvec = t2 - t1;
    
    /*for ROM*/
    ROM_std_hlpod_solve_ROM(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom),
            sys->fe.total_num_nodes,
            1,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_DIAG);
    
    t1 = monolis_get_time_global_sync();
    ROM_sys_hlpod_fe_add_Dbc(
            sys->rom.hlpod_vals.sol_vec,
            &(sys->bc),
            sys->fe.total_num_nodes,
            1);

    monolis_mpi_update_R(&(sys->monolis_com), sys->fe.total_num_nodes, 1, sys->rom.hlpod_vals.sol_vec);

    t2 = monolis_get_time_global_sync();
    sys->rom.hlpod_vals.time_calc_reduced_matvec += t2 - t1;

    ROM_BB_vec_copy(sys->rom.hlpod_vals.sol_vec, sys->vals_rom.T, sys->fe.total_num_nodes);
}

void solver_rom_NR4(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    monolis_com_initialize_by_self(&(sys->mono_com0));

    ROM_std_hlpod_calc_reduced_mat2(
        &(sys->monolis),
        &(sys->monolis_rom),
        &(sys->mono_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        4);

}