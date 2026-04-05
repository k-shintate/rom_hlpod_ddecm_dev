
#include "core_ROM.h"
#include "core_HROM.h"

/*for ROM*/
const char*     ROM_ID_FINISH_TIME      = "#rom_finish_time";
const double  POD_DVAL_FINISH_TIME      = 1.0;
const char* POD_ID_OUTPUT_INTERVAL      = "#rom_output_interval";
const int POD_DVAL_OUTPUT_INTERVAL      = 1;
const char* POD_ID_SNAPSHOT_INTERVAL    = "#snapshot_interval";
const int POD_DVAL_SNAPSHOT_INTERVAL    = 1;

const char*         ROM_ID_DENSITY  = "#density";
const char*       ROM_ID_VISCOSITY  = "#viscosity";

static const char* POD_INPUT_FILENAME_COND          = "rom_cond.dat";
static const char* ROM_OUTPUT_FILENAME_VTK          = "rom_result_%06d.vtk";

static const char* INPUT_FILENAME_DENSITY           = "density.dat";    //for parametric study
static const char* INPUT_FILENAME_VISCOSITY         = "viscosity.dat";  //for parametric study

static const char* INPUT_DIRECTORYNAME_1STDD        = "parted.0/";
static const char* INPUT_DIRECTORYNAME_2NDD         = "parted.1/";
static const char* INPUT_DIRECTORYNAME_METAGRAPH    = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH         = "metagraph.dat";

const int BUFFER_SIZE = 10000;


void initialize_velocity_pressure(
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

void ROM_offline_set_reynolds_num_cases(
    VALUES*         vals,
    const char*     directory)
{
	FILE* fp;
	char id[BUFFER_SIZE];
	int num_case_density;
	int num_case_viscosity;

	fp = ROM_BB_read_fopen(fp, INPUT_FILENAME_DENSITY, directory);

	fscanf(fp, "%s %d", id, &(vals->num_cases));
	num_case_density = vals->num_cases;
	vals->density_cases = BB_std_calloc_1d_double(vals->density_cases, vals->num_cases);
	for (int i = 0; i < vals->num_cases; i++) {
		fscanf(fp, "%lf", &(vals->density_cases[i]));
	}
	fclose(fp);

	fp = ROM_BB_read_fopen(fp, INPUT_FILENAME_VISCOSITY, directory);

	fscanf(fp, "%s %d", id, &(num_case_viscosity));
	if(num_case_density != num_case_viscosity){
		printf("Error: The number of cases in density.dat and viscosity.dat are different.\n");
		exit(1);
	}
	vals->viscosity_cases = BB_std_calloc_1d_double(vals->viscosity_cases, vals->num_cases);
	for (int i = 0; i < vals->num_cases; i++) {
		fscanf(fp, "%lf", &(vals->viscosity_cases[i]));
	}
	fclose(fp);

}

void ROM_offline_set_reynolds_number( 
        VALUES*     vals,
        const int   case_id)
{
    vals->density   = vals->density_cases[case_id];
    vals->viscosity = vals->viscosity_cases[case_id];

    printf("\n%s ---------- Calculation Conditions: Parametric Study Case %d ----------\n", 
            CODENAME, case_id);

    printf("%s %s: %e\n", CODENAME, ROM_ID_DENSITY,   vals->density);
    printf("%s %s: %e\n", CODENAME, ROM_ID_VISCOSITY, vals->viscosity);
}


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


void set_target_parameter(
		VALUES*         vals,
		const char*     directory)
{
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int tmp;
	double target_vals;

	snprintf(fname, BUFFER_SIZE, "target_density.dat");
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(tmp));
	fscanf(fp, "%lf", &(target_vals));
	vals->density = target_vals;
	fclose(fp);

	snprintf(fname, BUFFER_SIZE, "target_viscosity.dat");
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(tmp));
	fscanf(fp, "%lf", &(target_vals));
	vals->viscosity = target_vals;
	fclose(fp);

	printf("\n%s ---------- Calculation condition - ROM %d  ----------\n", CODENAME);

	printf("%s %s: %e\n", CODENAME, ROM_ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ROM_ID_VISCOSITY,        vals->viscosity);

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

void ROM_output_vtk_shape(
		BBFE_DATA*     fe,
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

	fclose(fp);
}


void ROM_add_output_vtk_pressure(
		BBFE_DATA*     fe,
		double*        fem_pressure,
		double*    	   rom_pressure,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	/*for pressure*/
	BB_vtk_write_point_vals_scalar(fp, rom_pressure, fe->total_num_nodes, "ROM-Pressure");

	double* error_p;
	error_p = BB_std_calloc_1d_double(error_p, fe->total_num_nodes);
	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, rom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_ROM-Pressure");

	BB_std_free_1d_double(error_p, fe->total_num_nodes);

	fclose(fp);
}


void ROM_add_output_vtk_velocity(
		BBFE_DATA*     fe,
		double**       fem_velocity,
		double**	   rom_velocity,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	double** error_v;
	error_v = BB_std_calloc_2d_double(error_v, fe->total_num_nodes, 3);
	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(fem_velocity[i][j] - rom_velocity[i][j]);
		}
	}
	
	/*for velocity*/
	BB_vtk_write_point_vals_vector(fp, rom_velocity, fe->total_num_nodes, "ROM-Velocity");
	BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_ROM-Velosity");

	fclose(fp);

	BB_std_free_2d_double(error_v, fe->total_num_nodes, 3);
}


void ROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, ROM_OUTPUT_FILENAME_VTK, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	ROM_output_vtk_shape(
			&(sys->fe),
			filename,
			sys->cond.directory,
			t);
	ROM_add_output_vtk_pressure(
			&(sys->fe),
			sys->vals.p,
			sys->vals_rom.p,
			filename,
			sys->cond.directory,
			t);
	ROM_add_output_vtk_velocity(
			&(sys->fe),
			sys->vals.v,
			sys->vals_rom.v,
			filename,
			sys->cond.directory,
			t);
	
	double L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals_rom.p,
			(const double*)sys->vals.p);

	printf("%s L2 error: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	double L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals.v,
			(const double**)sys->vals_rom.v);

	printf("%s L2 error: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}
}


void solver_rom(
    FE_SYSTEM* sys,
    const int step_HR,
    const int step_rom,
    const double t)
{
    monolis_clear_mat_value_R(&(sys->monolis));
	monolis_clear_mat_value_R(&(sys->monolis_rom));
	monolis_com_initialize_by_self(&(sys->mono_com0));

	if(monolis_mpi_get_global_comm_size() == 1){
	}
	else{
		monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
	} 

    set_element_mat(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));
        
    set_element_vec(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc),
        sys->monolis.mat.R.B);

	ROM_std_hlpod_calc_reduced_mat(
		&(sys->monolis),
        &(sys->monolis_rom),
        &(sys->mono_com),
		&(sys->mono_com0),
        &(sys->mono_com_rom_solv),
		&(sys->rom_sups),
        sys->fe.total_num_nodes,
		4);

    ROM_std_hlpod_solve_ROM(
        &(sys->monolis),
        &(sys->monolis_rom),
        &(sys->mono_com_rom_solv),
		&(sys->rom_sups),
        sys->fe.total_num_nodes,
		4,
		sys->vals.mat_max_iter,
		sys->vals.mat_epsilon,
		MONOLIS_ITER_BICGSTAB,
		MONOLIS_PREC_DIAG);

	ROM_sys_hlpod_fe_add_Dbc(
        sys->rom_sups.hlpod_vals.sol_vec,
		&(sys->bc),
		sys->fe.total_num_nodes,
		4);
	
	monolis_mpi_update_R(&(sys->mono_com), sys->fe.total_num_nodes, 4, sys->rom_sups.hlpod_vals.sol_vec);

    BBFE_fluid_sups_renew_velocity(
        sys->vals_rom.v, 
        sys->rom_sups.hlpod_vals.sol_vec,
        sys->fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys->vals_rom.p, 
        sys->rom_sups.hlpod_vals.sol_vec,
        sys->fe.total_num_nodes);
        
    output_files(sys, step_rom, t);

}


void solver_rom_global_para(
	MONOLIS*     monolis,
	MONOLIS_COM* monolis_com,
	ROM*		 rom,
    FE_SYSTEM*   sys,
    const int    step_HR,
    const int    step_rom,
    const double t)
{
    int total_num_nodes = sys->fe.total_num_nodes;
    double** mat = BB_std_calloc_2d_double(mat, rom->hlpod_vals.num_modes, rom->hlpod_vals.num_modes);
	double* rhs = BB_std_calloc_1d_double(rhs, rom->hlpod_vals.num_modes);
    double* mode_coef = BB_std_calloc_1d_double(mode_coef, rom->hlpod_vals.num_modes);
    double* ansvec = BB_std_calloc_1d_double(ansvec, total_num_nodes* 4);

	set_D_bc_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
			&(sys->bc),
            rhs,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);
	
	set_reduced_vec_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
            rhs,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);

	set_reduced_mat_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
            mat,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);
	
    allreduce_global_para(
            monolis_com,
            mat,
            rhs,
            rom->hlpod_vals.num_modes);

	ROM_BB_gauss_elimination(
			rom->hlpod_vals.num_modes,
			mat,
			rhs,
			mode_coef);
		
    ROM_std_hlpod_calc_sol_global_para(
            ansvec,
            rom->hlpod_mat.pod_modes,
            mode_coef,
			sys->fe.total_num_nodes,
			rom->hlpod_vals.num_modes,
			4);
	
	ROM_sys_hlpod_fe_add_Dbc(
            ansvec,
			&(sys->bc),
			sys->fe.total_num_nodes,
			4);
				
	/*解ベクトルのupdate*/
	monolis_mpi_update_R(monolis_com, sys->fe.total_num_nodes, 4, ansvec);

	BB_std_free_2d_double(mat, rom->hlpod_vals.num_modes, rom->hlpod_vals.num_modes);
	BB_std_free_1d_double(rhs, rom->hlpod_vals.num_modes);
	BB_std_free_1d_double(mode_coef, rom->hlpod_vals.num_modes);

    BBFE_fluid_sups_renew_velocity(
        sys->vals_rom.v, 
        ansvec,
        sys->fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys->vals_rom.p, 
        ansvec,
        sys->fe.total_num_nodes);
    
    BB_std_free_1d_double(ansvec, total_num_nodes* 4);
}


void solver_rom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step,
    const int step_hrom)
{
	if(monolis_mpi_get_global_my_rank()==0){
    		printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
	}

    double* rvec_m = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));
    double* rvec_c = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys.fe.total_num_nodes*4, sizeof(double));


    monolis_clear_mat_value_R(&(sys.monolis));
        monolis_clear_mat_value_R(&(sys.monolis_rom));
        monolis_com_initialize_by_self(&(sys.mono_com0));

        if(monolis_mpi_get_global_comm_size() == 1){
        }
        else{
                monolis_copy_mat_value_R(&(sys.monolis_rom0), &(sys.monolis_rom));
        }

	const double rel_tol_v = 1.0e-6;
const double abs_tol_v = 1.0e-12;
const double rel_tol_p = 1.0e-6;
const double abs_tol_p = 1.0e-12;
const double tiny      = 1.0e-30;
int max_iter_NR = 5;


monolis_com_initialize_by_self(&(sys.mono_com0));

    for(int it = 0; it < max_iter_NR; it++){
		if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
		}
        
        monolis_clear_mat_value_R(&(sys.monolis));
		monolis_clear_mat_value_R(&(sys.monolis_rom));
        monolis_clear_mat_value_R(&(sys.monolis_rom0));

        monolis_com_initialize_by_self(&(sys.mono_com0));

                set_element_mat_NR_Tezuer(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals_rom));

                BBFE_sys_monowrap_set_Dirichlet_bc(
                                &(sys.monolis),
                                sys.fe.total_num_nodes,
                                4,
                                &(sys.bc_NR),
                                sys.monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
	    for(int i=0; i<sys.fe.total_num_nodes; ++i){
        	for(int j =0; j<3; j++){
                	rvec_m_old[i*4 + j] = sys.monolis.mat.R.B[i*4 + j];
        	}
        		rvec_c_old[i*3] = sys.monolis.mat.R.B[i*3];
		}
            }

                monolis_show_timelog (&(sys.monolis_rom), true);
                monolis_show_iterlog (&(sys.monolis_rom), true);

        ROM_std_hlpod_calc_reduced_mat(
                &(sys.monolis),
        &(sys.monolis_rom),
        &(sys.mono_com),
                &(sys.mono_com0),
        &(sys.mono_com_rom_solv),
                &(sys.rom_sups),
        sys.fe.total_num_nodes,
                4);

    ROM_std_hlpod_solve_ROM(
        &(sys.monolis),
        &(sys.monolis_rom),
        &(sys.mono_com_rom_solv),
                &(sys.rom_sups),
        sys.fe.total_num_nodes,
                4,
                sys.vals.mat_max_iter,
                sys.vals.mat_epsilon,
                MONOLIS_ITER_BICGSAFE,
                MONOLIS_PREC_DIAG);

        monolis_mpi_update_R(
            &(sys.mono_com),
            sys.fe.total_num_nodes,
            4,
            sys.rom_sups.hlpod_vals.sol_vec);


                BBFE_fluid_sups_renew_velocity(
                                sys.vals_rom.delta_v,
                                //sys.monolis.mat.R.X,
				sys.rom_sups.hlpod_vals.sol_vec,
                                sys.fe.total_num_nodes);

                BBFE_fluid_sups_renew_pressure(
                                sys.vals_rom.delta_p,
                                //sys.monolis.mat.R.X,
				sys.rom_sups.hlpod_vals.sol_vec,
                                sys.fe.total_num_nodes);

        update_velocity_pressure_NR(
                sys.vals_rom.v,
                sys.vals_rom.delta_v,
                sys.vals_rom.p,
                sys.vals_rom.delta_p,
                                sys.fe.total_num_nodes);


monolis_clear_mat_value_R(&(sys.monolis));

        for(int i=0; i<sys.fe.total_num_nodes*4; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;

        }


set_element_mat_NR_Tezuer(
                                &(sys.monolis),
                                &(sys.fe),
                                &(sys.basis),
                                &(sys.vals_rom));

                BBFE_sys_monowrap_set_Dirichlet_bc(
                                &(sys.monolis),
                                sys.fe.total_num_nodes,
                                4,
                                &(sys.bc_NR),
                                sys.monolis.mat.R.B);


for(int i=0; i<sys.fe.total_num_nodes; ++i){
	for(int j =0; j<3; j++){ 
		rvec_m[i*4 + j] = sys.monolis.mat.R.B[i*4 + j];
	}
	rvec_c[i*3] = sys.monolis.mat.R.B[i*3];
}
        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys.mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys.mono_com.n_internal_vertex*4,
            1);
	
        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys.mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys.mono_com.n_internal_vertex*4,
            1);


double norm_v = calc_internal_norm_2d(
    sys.vals_rom.v,
    sys.mono_com.n_internal_vertex,
    3);

double norm_delta_v = calc_internal_norm_2d(
    sys.vals_rom.delta_v,
    sys.mono_com.n_internal_vertex,
    3);

/* 圧力の L2 ノルム（内部自由度のみ） */
double norm_p = 0.0, norm_delta_p = 0.0;
for (int ii = 0; ii < sys.mono_com.n_internal_vertex; ++ii) {
    double pv  = sys.vals_rom.p[ii];
    double dpv = sys.vals_rom.delta_p[ii];
    norm_p       += pv  * pv;
    norm_delta_p += dpv * dpv;
}

/* L∞（最大変化量）：速度は3成分、圧力は1成分 */
double linf_delta_v_local = 0.0;
double linf_delta_p_local = 0.0;
for (int i_node = 0; i_node < sys.mono_com.n_internal_vertex; ++i_node) {
    for (int d = 0; d < 3; ++d) {
        double av = fabs(sys.vals_rom.delta_v[i_node][d]);
        if (av > linf_delta_v_local) linf_delta_v_local = av;
    }
    double ap = fabs(sys.vals_rom.delta_p[i_node]);
    if (ap > linf_delta_p_local) linf_delta_p_local = ap;
}

/* MPI で集約（L2 和：SUM、L∞：MAX） */
monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.mono_com.comm);
monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old, MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys.mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old, MONOLIS_MPI_SUM, sys.mono_com.comm);

/* L∞は MAX で集約（Monolis に MAX が無ければ、自前で rank0 に gather→max でもOK） */
double linf_delta_v = linf_delta_v_local;
double linf_delta_p = linf_delta_p_local;
monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.mono_com.comm);
monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.mono_com.comm);

/* ルートを取って実ノルムに */
double nrm_v        = sqrt(norm_v);
double nrm_dv       = sqrt(norm_delta_v);
double nrm_p        = sqrt(norm_p);
double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m        = sqrt(norm_r_m);
        double nrm_r_m_old       = sqrt(norm_r_m_old);
        double nrm_r_c        = sqrt(norm_r_c);
        double nrm_r_c_old       = sqrt(norm_r_c_old);

/* 相対＋絶対の複合判定（ゼロ割回避） */
double denom_v = fmax(nrm_v,  tiny);
double denom_p = fmax(nrm_p,  tiny);

int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

/* ログ出力を見やすく */
if(monolis_mpi_get_global_my_rank()==0){
printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
       it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
       it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old, nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);

}
/* 収束したら後処理 */
//if (conv_v && conv_p) {
if(it == max_iter_NR -1){
    /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
    double max_du = 0.0;
    for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {
        for (int d = 0; d < 3; ++d) {
            double du = fabs(sys.vals_rom.v[ii][d] - sys.vals_rom.v_old[ii][d]);
            if (du > max_du) max_du = du;
        }
    }
    if(monolis_mpi_get_global_my_rank()==0){
    	printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
    }
    ROM_BB_vec_copy_2d(
        sys.vals_rom.v,
        sys.vals_rom.v_old,
        sys.fe.total_num_nodes,
        3);
    break;
}

    }

BB_std_free_1d_double(rvec_m, sys.fe.total_num_nodes*4);
BB_std_free_1d_double(rvec_m_old, sys.fe.total_num_nodes*4);
BB_std_free_1d_double(rvec_c, sys.fe.total_num_nodes*4);
BB_std_free_1d_double(rvec_c_old, sys.fe.total_num_nodes*4);


}




void solver_rom_NR2(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    if(monolis_mpi_get_global_my_rank()==0){
        printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
    }

    double* rvec_m     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));

    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_com_initialize_by_self(&(sys->mono_com0));

    if(monolis_mpi_get_global_comm_size() == 1){
    }
    else{
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
    }

    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 5;

    monolis_com_initialize_by_self(&(sys->mono_com0));

    for(int it = 0; it < max_iter_NR; it++){
        if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
        }

        monolis_clear_mat_value_R(&(sys->monolis));
        monolis_clear_mat_value_R(&(sys->monolis_rom));
        monolis_clear_mat_value_R(&(sys->monolis_rom0));

        monolis_com_initialize_by_self(&(sys->mono_com0));

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys->fe.total_num_nodes; ++i){
                for(int j=0; j<3; j++){
                    rvec_m_old[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
                }
                rvec_c_old[i*3] = sys->monolis.mat.R.B[i*3];
            }
        }

        //monolis_show_timelog(&(sys->monolis_rom), true);
        //monolis_show_iterlog(&(sys->monolis_rom), true);

        ROM_std_hlpod_calc_reduced_mat(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com),
            &(sys->mono_com0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4);

        ROM_std_hlpod_solve_ROM(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->rom_sups.hlpod_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_rom.delta_v,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_rom.delta_p,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_rom.v,
            sys->vals_rom.delta_v,
            sys->vals_rom.p,
            sys->vals_rom.delta_p,
            sys->fe.total_num_nodes);

        monolis_clear_mat_value_R(&(sys->monolis));

        for(int i=0; i<sys->fe.total_num_nodes*4; ++i){
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        for(int i=0; i<sys->fe.total_num_nodes; ++i){
            for(int j=0; j<3; j++){
                rvec_m[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
            }
            rvec_c[i*3] = sys->monolis.mat.R.B[i*3];
        }

        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_v = calc_internal_norm_2d(
            sys->vals_rom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_rom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_rom.p[ii];
            double dpv = sys->vals_rom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_rom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_rom.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        /* MPI で集約（L2 和：SUM、L∞：MAX） */
        monolis_allreduce_R(1, &norm_v,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);

        /* L∞は MAX で集約 */
        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        /* ルートを取って実ノルムに */
        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m      = sqrt(norm_r_m);
        double nrm_r_m_old  = sqrt(norm_r_m_old);
        double nrm_r_c      = sqrt(norm_r_c);
        double nrm_r_c_old  = sqrt(norm_r_c_old);

        /* 相対＋絶対の複合判定（ゼロ割回避） */
        double denom_v = fmax(nrm_v, tiny);
        double denom_p = fmax(nrm_p, tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        /* ログ出力を見やすく */
        if(monolis_mpi_get_global_my_rank()==0){
            printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                   it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);

            printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
                   it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p,
                   nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old,
                   nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);
        }

        /* 収束したら後処理 */
        //if (conv_v && conv_p) {
        if(it == max_iter_NR - 1){
            /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
            double max_du = 0.0;
            for (int ii = 0; ii < sys->fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys->vals_rom.v[ii][d] - sys->vals_rom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_rom.v,
                sys->vals_rom.v_old,
                sys->fe.total_num_nodes,
                3);

            break;
        }
    }

    BB_std_free_1d_double(rvec_m,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_m_old, sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c_old, sys->fe.total_num_nodes*4);
}


void solver_rom_linear(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_clear_mat_value_R(&(sys->monolis_rom0));
    monolis_com_initialize_by_self(&(sys->mono_com0));

    set_element_mat_NR_linear(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    ROM_std_hlpod_calc_reduced_mat(
        &(sys->monolis),
        &(sys->monolis_rom0),
        &(sys->mono_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        4);

}

void add_reduced_mat_linear(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_clear_mat_value_R(&(sys->monolis_rom0));
    monolis_com_initialize_by_self(&(sys->mono_com0));
double t1 = monolis_get_time_global_sync();
printf("test1");
    set_element_mat_NR_linear(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));
double t2 = monolis_get_time_global_sync();
printf("test2");

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

double t3 = monolis_get_time_global_sync();
printf("test3");

    ROM_std_hlpod_calc_reduced_mat(
        &(sys->monolis),
        &(sys->monolis_rom0),
        &(sys->mono_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        4);
/*
    monolis_clear_mat_value_R(&(sys->monolis));
    //monolis_initialize(&(sys->monolis_mass_rom0));
    //monolis_copy_mat_R(&(sys->monolis_rom0), &(sys->monolis_mass_rom0));
    //monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_clear_mat_value_R(&(sys->monolis_mass_rom0));
    //monolis_com_initialize_by_self(&(sys->mono_com0));

    set_element_mat_NR_mass(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    ROM_std_hlpod_calc_reduced_mat(
        &(sys->monolis),
        &(sys->monolis_mass_rom0),
        &(sys->mono_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        4);
*/
}


void ROM_std_hlpod_reduced_rhs_to_monollis_linear(
    MONOLIS*		monolis,
    MONOLIS_COM*    monolis_com,
    HLPOD_MAT*      hlpod_mat,
    double*         mode_coeff,
    const int       num_2nd_subdomains)
{
    int index = 0;
    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            hlpod_mat->VTf[index + i] = 0.0;
        }
        index += hlpod_mat->num_modes_internal[k];
    }

    monolis_matvec_product_R(monolis, monolis_com, mode_coeff, hlpod_mat->VTf);

    for(int k = 0; k < num_2nd_subdomains; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            //printf("%lf ", hlpod_mat->VTf[index + i]);
            //monolis->mat.R.B[index + i] = hlpod_mat->VTf[index + i];
            hlpod_mat->VTf[index + i]  = - hlpod_mat->VTf[index + i];
        }
        index += hlpod_mat->num_modes_internal[k];
    }
    
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

void solver_rom_NR3(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    if(monolis_mpi_get_global_my_rank()==0){
        printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
    }

    double* rvec_m     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));

    /*
    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_com_initialize_by_self(&(sys->mono_com0));

    if(monolis_mpi_get_global_comm_size() == 1){
    }
    else{
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
    }
    */

    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 3;

    monolis_com_initialize_by_self(&(sys->mono_com0));

    for(int it = 0; it < max_iter_NR; it++){
        if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
        }

        monolis_clear_mat_value_R(&(sys->monolis));
        monolis_clear_mat_value_R(&(sys->monolis_rom));
        //monolis_clear_mat_value_R(&(sys->monolis_rom0));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_rom0));
        monolis_com_initialize_by_self(&(sys->mono_com0));
/*
        ROM_std_hlpod_reduced_rhs_to_monollis_linear(
            &(sys->monolis_rom0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_mat.mode_coef,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains);
*/
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));

/*
        set_element_mat_NR_linear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/
        set_element_mat_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        set_element_vec_NR_linear_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

/*
        set_element_vec_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/
        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys->fe.total_num_nodes; ++i){
                for(int j=0; j<3; j++){
                    rvec_m_old[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
                }
                rvec_c_old[i*3] = sys->monolis.mat.R.B[i*3];
            }
        }

        //monolis_show_timelog(&(sys->monolis_rom), true);
        //monolis_show_iterlog(&(sys->monolis_rom), true);

        
        ROM_std_hlpod_calc_reduced_mat(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com),
            &(sys->mono_com0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4);

        ROM_std_hlpod_solve_ROM(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

    int index = 0;
    for(int k = 0; k < sys->rom_sups.hlpod_vals.num_2nd_subdomains; k++){
        for(int i = 0; i < sys->rom_sups.hlpod_mat.num_modes_internal[k]; i++){
           sys->rom_sups.hlpod_mat.mode_coef_old[index + i] = sys->rom_sups.hlpod_mat.mode_coef[index + i];
            //printf("mode_coeff_old = %e", sys->rom_sups.hlpod_mat.mode_coef_old[index + i]);
        }
        index += sys->rom_sups.hlpod_mat.num_modes_internal[k];
    }

    
    //if(step_hrom%2==0 && it == max_iter_NR - 1){
    if(step_hrom%2==0){
    HROM_get_neib_coordinates(
            &(sys->mono_com_rom),
            &(sys->rom_sups.hlpod_vals),
            &(sys->rom_sups.hlpod_mat),
            1 + sys->mono_com_rom_solv.recv_n_neib,
            sys->rom_sups.hlpod_vals.num_modes_max,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->rom_sups.hlpod_vals.num_modes_pre);
/*
    HROM_ddecm_set_residuals_NR_blas2(
    //HROM_ddecm_set_residuals_NR_PSPG(
            &(sys->fe),
            &(sys->basis),
            &(sys->vals),
            &(sys->bc),
            &(sys->rom_sups.hlpod_mat),
            &(sys->rom_sups.hlpod_vals),
            &(sys->hrom_sups.hlpod_ddhr),
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            step_hrom -1 ,   //index 0 start
            sys->rom_sups.hlpod_vals.num_snapshot,
            1 + sys->mono_com.recv_n_neib,
            sys->vals.dt,
            t);
  */  
    //HROM_ddecm_set_residuals_NR_blas(
    HROM_ddecm_set_residuals_NR_vec(
    //HROM_ddecm_set_residuals_NR_PSPG(
            &(sys->fe),
            &(sys->basis),
            &(sys->vals),
            &(sys->bc),
            &(sys->rom_sups.hlpod_mat),
            &(sys->rom_sups.hlpod_vals),
            &(sys->hrom_sups.hlpod_ddhr),
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            step_hrom -1 ,   //index 0 start
            sys->rom_sups.hlpod_vals.num_snapshot,
            1 + sys->mono_com.recv_n_neib,
            sys->vals.dt,
            t);

	    }

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->rom_sups.hlpod_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_rom.delta_v,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_rom.delta_p,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_rom.v,
            sys->vals_rom.delta_v,
            sys->vals_rom.p,
            sys->vals_rom.delta_p,
            sys->fe.total_num_nodes);

        monolis_clear_mat_value_R(&(sys->monolis));

        for(int i=0; i<sys->fe.total_num_nodes*4; ++i){
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        for(int i=0; i<sys->fe.total_num_nodes; ++i){
            for(int j=0; j<3; j++){
                rvec_m[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
            }
            rvec_c[i*3] = sys->monolis.mat.R.B[i*3];
        }

        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_v = calc_internal_norm_2d(
            sys->vals_rom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_rom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_rom.p[ii];
            double dpv = sys->vals_rom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_rom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_rom.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        /* MPI で集約（L2 和：SUM、L∞：MAX） */
        monolis_allreduce_R(1, &norm_v,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);

        /* L∞は MAX で集約 */
        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        /* ルートを取って実ノルムに */
        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m      = sqrt(norm_r_m);
        double nrm_r_m_old  = sqrt(norm_r_m_old);
        double nrm_r_c      = sqrt(norm_r_c);
        double nrm_r_c_old  = sqrt(norm_r_c_old);

        /* 相対＋絶対の複合判定（ゼロ割回避） */
        double denom_v = fmax(nrm_v, tiny);
        double denom_p = fmax(nrm_p, tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        /* ログ出力を見やすく */
        if(monolis_mpi_get_global_my_rank()==0){
            printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                   it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);

            printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
                   it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p,
                   nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old,
                   nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);
        }

        /* 収束したら後処理 */
        //if (conv_v && conv_p) {
        if(it == max_iter_NR - 1){
            /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
            double max_du = 0.0;
            for (int ii = 0; ii < sys->fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys->vals_rom.v[ii][d] - sys->vals_rom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_rom.v,
                sys->vals_rom.v_old,
                sys->fe.total_num_nodes,
                3);

            break;
        }
    }

    BB_std_free_1d_double(rvec_m,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_m_old, sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c_old, sys->fe.total_num_nodes*4);
    
}


void solver_rom_NR5(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    if(monolis_mpi_get_global_my_rank()==0){
        printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
    }

    double* rvec_m     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));


    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 3;

    monolis_com_initialize_by_self(&(sys->mono_com0));

    for(int it = 0; it < max_iter_NR; it++){
        if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
        }

        monolis_clear_mat_value_R(&(sys->monolis));
        monolis_clear_mat_value_R(&(sys->monolis_rom));
        //monolis_clear_mat_value_R(&(sys->monolis_rom0));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_rom0));
        monolis_com_initialize_by_self(&(sys->mono_com0));
/*
        ROM_std_hlpod_reduced_rhs_to_monollis_linear(
            &(sys->monolis_rom0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_mat.mode_coef,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains);
*/
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));

/*
        set_element_mat_NR_linear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/

        set_element_mat_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        set_element_vec_NR_linear_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

/*
        set_element_vec_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/
        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys->fe.total_num_nodes; ++i){
                for(int j=0; j<3; j++){
                    rvec_m_old[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
                }
                rvec_c_old[i*3] = sys->monolis.mat.R.B[i*3];
            }
        }

        //monolis_show_timelog(&(sys->monolis_rom), true);
        //monolis_show_iterlog(&(sys->monolis_rom), true);

  
        ROM_std_hlpod_calc_reduced_mat(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com),
            &(sys->mono_com0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4);

/*
        HROM_set_element_mat_NR(
            &(sys->monolis_rom),
            &(sys->fe),
            &(sys->vals_rom),
            &(sys->basis),
            &(sys->bc_NR),
            &(sys->rom_sups.hlpod_vals),
            &(sys->rom_sups.hlpod_mat),
            &(sys->hrom_sups.hlpod_ddhr),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->vals.dt);
        
        HROM_ddecm_calc_block_mat_bcsr(
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups.hlpod_vals),
            &(sys->rom_sups.hlpod_mat),
            &(sys->hrom_sups.hlpod_ddhr),
            &(sys->rom_sups.hlpod_meta),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->cond.directory);	
*/
        ROM_std_hlpod_solve_ROM(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->rom_sups.hlpod_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_rom.delta_v,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_rom.delta_p,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_rom.v,
            sys->vals_rom.delta_v,
            sys->vals_rom.p,
            sys->vals_rom.delta_p,
            sys->fe.total_num_nodes);

        monolis_clear_mat_value_R(&(sys->monolis));

        for(int i=0; i<sys->fe.total_num_nodes*4; ++i){
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        for(int i=0; i<sys->fe.total_num_nodes; ++i){
            for(int j=0; j<3; j++){
                rvec_m[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
            }
            rvec_c[i*3] = sys->monolis.mat.R.B[i*3];
        }

        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_v = calc_internal_norm_2d(
            sys->vals_rom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_rom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_rom.p[ii];
            double dpv = sys->vals_rom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_rom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_rom.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        /* MPI で集約（L2 和：SUM、L∞：MAX） */
        monolis_allreduce_R(1, &norm_v,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);

        /* L∞は MAX で集約 */
        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        /* ルートを取って実ノルムに */
        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m      = sqrt(norm_r_m);
        double nrm_r_m_old  = sqrt(norm_r_m_old);
        double nrm_r_c      = sqrt(norm_r_c);
        double nrm_r_c_old  = sqrt(norm_r_c_old);

        /* 相対＋絶対の複合判定（ゼロ割回避） */
        double denom_v = fmax(nrm_v, tiny);
        double denom_p = fmax(nrm_p, tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        /* ログ出力を見やすく */
        if(monolis_mpi_get_global_my_rank()==0){
            printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                   it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);

            printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
                   it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p,
                   nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old,
                   nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);
        }

        /* 収束したら後処理 */
        //if (conv_v && conv_p) {
        if(it == max_iter_NR - 1){
            /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
            double max_du = 0.0;
            for (int ii = 0; ii < sys->fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys->vals_rom.v[ii][d] - sys->vals_rom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_rom.v,
                sys->vals_rom.v_old,
                sys->fe.total_num_nodes,
                3);

            break;
        }
    }

    BB_std_free_1d_double(rvec_m,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_m_old, sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c_old, sys->fe.total_num_nodes*4);
}



void solver_rom_NR6(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    if(monolis_mpi_get_global_my_rank()==0){
        printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
    }

    double* rvec_m     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));


    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 1;

    monolis_com_initialize_by_self(&(sys->mono_com0));

    int index = 0;
    for(int k = 0; k < sys->rom_sups.hlpod_vals.num_2nd_subdomains; k++){
        for(int i = 0; i < sys->rom_sups.hlpod_mat.num_modes_internal[k]; i++){
        }
        index += sys->rom_sups.hlpod_mat.num_modes_internal[k];
    }
    //double* mode_coeff_old = (double*)calloc((size_t)index, sizeof(double));

    index = 0;
    for(int k = 0; k < sys->rom_sups.hlpod_vals.num_2nd_subdomains; k++){
        for(int i = 0; i < sys->rom_sups.hlpod_mat.num_modes_internal[k]; i++){
           //sys->rom_sups.hlpod_mat.mode_coef_old[index + i] = sys->rom_sups.hlpod_mat.mode_coef[index + i];
           sys->rom_sups.hlpod_mat.mode_coef[index + i] = sys->rom_sups.hlpod_mat.mode_coef_old[index + i];
            //printf("mode_coeff_old = %e", sys->rom_sups.hlpod_mat.mode_coef_old[index + i]);
        }
        index += sys->rom_sups.hlpod_mat.num_modes_internal[k];
    }

    for(int it = 0; it < max_iter_NR; it++){
        if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
        }

        monolis_clear_mat_value_R(&(sys->monolis));
        monolis_clear_mat_value_R(&(sys->monolis_rom));
        //monolis_clear_mat_value_R(&(sys->monolis_rom0));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_rom0));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_mass_rom0));
        monolis_com_initialize_by_self(&(sys->mono_com0));

        ROM_std_hlpod_reduced_rhs_to_monollis_linear2(
            &(sys->monolis_rom0),
            &(sys->monolis_mass_rom0),
            &(sys->mono_com_rom),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_mat.mode_coef,
            sys->rom_sups.hlpod_mat.mode_coef_old,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains);

        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));

/*
        set_element_mat_NR_linear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/

        set_element_mat_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
/*
        set_element_vec_NR_linear_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/

        set_element_vec_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys->fe.total_num_nodes; ++i){
                for(int j=0; j<3; j++){
                    rvec_m_old[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
                }
                rvec_c_old[i*3] = sys->monolis.mat.R.B[i*3];
            }
        }

        //monolis_show_timelog(&(sys->monolis_rom), true);
        //monolis_show_iterlog(&(sys->monolis_rom), true);

  
        ROM_std_hlpod_calc_reduced_mat(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com),
            &(sys->mono_com0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4);

        ROM_std_hlpod_solve_ROM_NR(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

    index = 0;
    for(int k = 0; k < sys->rom_sups.hlpod_vals.num_2nd_subdomains; k++){
        for(int i = 0; i < sys->rom_sups.hlpod_mat.num_modes_internal[k]; i++){
           sys->rom_sups.hlpod_mat.mode_coef_old[index + i] = sys->rom_sups.hlpod_mat.mode_coef[index + i];
            //printf("mode_coeff_old = %e", sys->rom_sups.hlpod_mat.mode_coef_old[index + i]);
        }
        index += sys->rom_sups.hlpod_mat.num_modes_internal[k];
    }

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->rom_sups.hlpod_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_rom.delta_v,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_rom.delta_p,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_rom.v,
            sys->vals_rom.delta_v,
            sys->vals_rom.p,
            sys->vals_rom.delta_p,
            sys->fe.total_num_nodes);

        monolis_clear_mat_value_R(&(sys->monolis));

        for(int i=0; i<sys->fe.total_num_nodes*4; ++i){
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        for(int i=0; i<sys->fe.total_num_nodes; ++i){
            for(int j=0; j<3; j++){
                rvec_m[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
            }
            rvec_c[i*3] = sys->monolis.mat.R.B[i*3];
        }

        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_v = calc_internal_norm_2d(
            sys->vals_rom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_rom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_rom.p[ii];
            double dpv = sys->vals_rom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_rom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_rom.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        /* MPI で集約（L2 和：SUM、L∞：MAX） */
        monolis_allreduce_R(1, &norm_v,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);

        /* L∞は MAX で集約 */
        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        /* ルートを取って実ノルムに */
        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m      = sqrt(norm_r_m);
        double nrm_r_m_old  = sqrt(norm_r_m_old);
        double nrm_r_c      = sqrt(norm_r_c);
        double nrm_r_c_old  = sqrt(norm_r_c_old);

        /* 相対＋絶対の複合判定（ゼロ割回避） */
        double denom_v = fmax(nrm_v, tiny);
        double denom_p = fmax(nrm_p, tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        /* ログ出力を見やすく */
        if(monolis_mpi_get_global_my_rank()==0){
            printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                   it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);

            printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
                   it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p,
                   nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old,
                   nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);
        }

        /* 収束したら後処理 */
        //if (conv_v && conv_p) {
        if(it == max_iter_NR - 1){
            /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
            double max_du = 0.0;
            for (int ii = 0; ii < sys->fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys->vals_rom.v[ii][d] - sys->vals_rom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_rom.v,
                sys->vals_rom.v_old,
                sys->fe.total_num_nodes,
                3);

            break;
        }
    }

    BB_std_free_1d_double(rvec_m,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_m_old, sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c_old, sys->fe.total_num_nodes*4);
}



void solver_rom_NR3(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    if(monolis_mpi_get_global_my_rank()==0){
        printf("\n%s ----------------- Time step %d ----------------\n", CODENAME, step);
    }

    double* rvec_m     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_m_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c     = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));
    double* rvec_c_old = (double*)calloc((size_t)sys->fe.total_num_nodes*4, sizeof(double));

    /*
    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_com_initialize_by_self(&(sys->mono_com0));

    if(monolis_mpi_get_global_comm_size() == 1){
    }
    else{
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
    }
    */

    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 3;

    monolis_com_initialize_by_self(&(sys->mono_com0));

    for(int it = 0; it < max_iter_NR; it++){
        if(monolis_mpi_get_global_my_rank()==0){
            printf("\n%s ----------------- ROM Time step %d : NR step %d ----------------\n", CODENAME, step, it);
        }

        monolis_clear_mat_value_R(&(sys->monolis));
        monolis_clear_mat_value_R(&(sys->monolis_rom));
        //monolis_clear_mat_value_R(&(sys->monolis_rom0));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_rom0));
        monolis_com_initialize_by_self(&(sys->mono_com0));
/*
        ROM_std_hlpod_reduced_rhs_to_monollis_linear(
            &(sys->monolis_rom0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_mat.mode_coef,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains);
*/
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));

/*
        set_element_mat_NR_linear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/
        set_element_mat_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        set_element_vec_NR_linear_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

/*
        set_element_vec_NR_nonlinear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));
*/
        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<sys->fe.total_num_nodes; ++i){
                for(int j=0; j<3; j++){
                    rvec_m_old[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
                }
                rvec_c_old[i*3] = sys->monolis.mat.R.B[i*3];
            }
        }

        //monolis_show_timelog(&(sys->monolis_rom), true);
        //monolis_show_iterlog(&(sys->monolis_rom), true);

        
        ROM_std_hlpod_calc_reduced_mat(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com),
            &(sys->mono_com0),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4);

        ROM_std_hlpod_solve_ROM(
            &(sys->monolis),
            &(sys->monolis_rom),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups),
            sys->fe.total_num_nodes,
            4,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

        int index = 0;
        for(int k = 0; k < sys->rom_sups.hlpod_vals.num_2nd_subdomains; k++){
            for(int i = 0; i < sys->rom_sups.hlpod_mat.num_modes_internal[k]; i++){
            sys->rom_sups.hlpod_mat.mode_coef_old[index + i] = sys->rom_sups.hlpod_mat.mode_coef[index + i];
                //printf("mode_coeff_old = %e", sys->rom_sups.hlpod_mat.mode_coef_old[index + i]);
            }
            index += sys->rom_sups.hlpod_mat.num_modes_internal[k];
        }

    
        //if(step_hrom%2==0 && it == max_iter_NR - 1){
        if(step_hrom%2==0){
            HROM_get_neib_coordinates(
                &(sys->mono_com_rom),
                &(sys->rom_sups.hlpod_vals),
                &(sys->rom_sups.hlpod_mat),
                1 + sys->mono_com_rom_solv.recv_n_neib,
                sys->rom_sups.hlpod_vals.num_modes_max,
                sys->rom_sups.hlpod_vals.num_2nd_subdomains,
                sys->rom_sups.hlpod_vals.num_modes_pre);

            HROM_ddecm_set_residuals_NR_blas(
                &(sys->fe),
                &(sys->basis),
                &(sys->vals),
                &(sys->bc),
                &(sys->rom_sups.hlpod_mat),
                &(sys->rom_sups.hlpod_vals),
                &(sys->hrom_p.hlpod_ddhr),
                sys->rom_sups.hlpod_vals.num_2nd_subdomains,
                step_hrom -1 ,   //index 0 start
                sys->rom_sups.hlpod_vals.num_snapshot,
                1 + sys->mono_com.recv_n_neib,
                sys->vals.dt,
                t);

            HROM_ddecm_set_residuals_NR_blas(
                &(sys->fe),
                &(sys->basis),
                &(sys->vals),
                &(sys->bc),
                &(sys->rom_sups.hlpod_mat),
                &(sys->rom_sups.hlpod_vals),
                &(sys->hrom_v.hlpod_ddhr),
                sys->rom_sups.hlpod_vals.num_2nd_subdomains,
                step_hrom -1 ,   //index 0 start
                sys->rom_sups.hlpod_vals.num_snapshot,
                1 + sys->mono_com.recv_n_neib,
                sys->vals.dt,
                t);

	    }

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->rom_sups.hlpod_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_rom.delta_v,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_rom.delta_p,
            sys->rom_sups.hlpod_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_rom.v,
            sys->vals_rom.delta_v,
            sys->vals_rom.p,
            sys->vals_rom.delta_p,
            sys->fe.total_num_nodes);

        monolis_clear_mat_value_R(&(sys->monolis));

        for(int i=0; i<sys->fe.total_num_nodes*4; ++i){
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe.total_num_nodes,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        for(int i=0; i<sys->fe.total_num_nodes; ++i){
            for(int j=0; j<3; j++){
                rvec_m[i*4 + j] = sys->monolis.mat.R.B[i*4 + j];
            }
            rvec_c[i*3] = sys->monolis.mat.R.B[i*3];
        }

        double norm_r_m_old = calc_internal_norm_1d(
            rvec_m_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_m = calc_internal_norm_1d(
            rvec_m,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c_old = calc_internal_norm_1d(
            rvec_c_old,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_r_c = calc_internal_norm_1d(
            rvec_c,
            sys->mono_com.n_internal_vertex*4,
            1);

        double norm_v = calc_internal_norm_2d(
            sys->vals_rom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_rom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_rom.p[ii];
            double dpv = sys->vals_rom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_rom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_rom.delta_p[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        /* MPI で集約（L2 和：SUM、L∞：MAX） */
        monolis_allreduce_R(1, &norm_v,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,         MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_m_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_r_c_old,   MONOLIS_MPI_SUM, sys->mono_com.comm);

        /* L∞は MAX で集約 */
        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        /* ルートを取って実ノルムに */
        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r_m      = sqrt(norm_r_m);
        double nrm_r_m_old  = sqrt(norm_r_m_old);
        double nrm_r_c      = sqrt(norm_r_c);
        double nrm_r_c_old  = sqrt(norm_r_c_old);

        /* 相対＋絶対の複合判定（ゼロ割回避） */
        double denom_v = fmax(nrm_v, tiny);
        double denom_p = fmax(nrm_p, tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        /* ログ出力を見やすく */
        if(monolis_mpi_get_global_my_rank()==0){
            printf("ROM : [NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                   it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);

            printf("ROM : [NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e |rm_old| = %.8e |rm| = %.8e  |rm|/|rm_old| = %.8e |rc_old| = %.8e |rc| = %.8e  |rc|/|rc_old| = %.8e\n",
                   it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p,
                   nrm_r_m_old, nrm_r_m, nrm_r_m / nrm_r_m_old,
                   nrm_r_c_old, nrm_r_c, nrm_r_c / nrm_r_c_old);
        }

        /* 収束したら後処理 */
        //if (conv_v && conv_p) {
        if(it == max_iter_NR - 1){
            /* ——— タイムステップ n+1 の NR 収束直後のレポート ——— */
            double max_du = 0.0;
            for (int ii = 0; ii < sys->fe.total_num_nodes; ++ii) {
                for (int d = 0; d < 3; ++d) {
                    double du = fabs(sys->vals_rom.v[ii][d] - sys->vals_rom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_rom.v,
                sys->vals_rom.v_old,
                sys->fe.total_num_nodes,
                3);

            break;
        }
    }

    BB_std_free_1d_double(rvec_m,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_m_old, sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c,     sys->fe.total_num_nodes*4);
    BB_std_free_1d_double(rvec_c_old, sys->fe.total_num_nodes*4);
    
}