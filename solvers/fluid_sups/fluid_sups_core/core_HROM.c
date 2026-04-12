
#include "core_HROM.h"

const int BUFFER_SIZE = 10000;

const char*     HROM_ID_INC_SVD_INTERVAL    = "#inc_svd_interval";
const int       INC_SVD_INTERVAL            = 10;

static const char* HROM_INPUT_FILENAME_COND          = "hrom_cond.dat";

static const char* HROM_OUTPUT_FILENAME_VTK          = "hrom_result_%06d.vtk";


void HROM_offline_assign_default_values_inc_svd(
		HR_VALUES*      hr_vals)
{
	hr_vals->incsvd_interval  = INC_SVD_INTERVAL;
}

void HROM_offline_print_all_values_inc_svd(
		HR_VALUES*     hr_vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, HROM_ID_INC_SVD_INTERVAL,  hr_vals->incsvd_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}

void HROM_offline_read_calc_conditions_inc_svd(
		HR_VALUES*      hr_vals,
		const char* 	directory)
{
	printf("\n");

	HROM_offline_assign_default_values_inc_svd(hr_vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, HROM_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading HROM conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_int_p(
				&(hr_vals->incsvd_interval), filename, HROM_ID_INC_SVD_INTERVAL, BUFFER_SIZE, CODENAME);


		fclose(fp);
	}

	HROM_offline_print_all_values_inc_svd(hr_vals);

	printf("\n");
}

void HROM_std_hlpod_offline_set_num_snapmat_inc_svd(
        ROM*            rom,
        const double    finish_time,
        const double    dt,
        const int       snapshot_interval,
        const int       num_case,
        const int       incsvd_interval)
{
    double quotient = finish_time / dt / snapshot_interval* num_case / incsvd_interval;

    printf("%s: finish_time = %f, dt = %f, snapshot_interval = %d, num_case = %d, incsvd_interval = %d\n",
           CODENAME, finish_time, dt, snapshot_interval, num_case, incsvd_interval);
    
    if (fmod(quotient, 1.0) == 0.0) {
    	rom->hlpod_vals.num_snapshot = incsvd_interval;
    	printf("%s: num_snapshot = %d\n", CODENAME, rom->hlpod_vals.num_snapshot);
    }
    else{
        printf("Error: num_snapshot = %d is not integer\n");
        exit(1);
    }
}

void HROM_output_vtk_shape(
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

void HROM_add_output_vtk_pressure(
		BBFE_DATA*     fe,
		double*        fem_pressure,
		double*    	   rom_pressure,
        double*    	   hrom_pressure,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	/*for pressure*/
    BB_vtk_write_point_vals_scalar(fp, hrom_pressure, fe->total_num_nodes, "HROM-Pressure");

	double* error_p;
	error_p = BB_std_calloc_1d_double(error_p, fe->total_num_nodes);
	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, rom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_ROM-Pressure");

	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, hrom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_HROM-Pressure");

	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, rom_pressure, hrom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-ROM_HROM-Pressure");

	BB_std_free_1d_double(error_p, fe->total_num_nodes);

	fclose(fp);
}


void HROM_add_output_vtk_velocity(
		BBFE_DATA*     fe,
		double**       fem_velocity,
		double**	   rom_velocity,
        double**	   hrom_velocity,
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
    BB_vtk_write_point_vals_vector(fp, hrom_velocity, fe->total_num_nodes, "HROM-Velocity");

	BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_ROM-Velosity");

	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(fem_velocity[i][j] - hrom_velocity[i][j]);
		}
	}

    BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_HROM-Velosity");

	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(rom_velocity[i][j] - hrom_velocity[i][j]);
		}
	}

    BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-ROM_HROM-Velosity");

	fclose(fp);

	BB_std_free_2d_double(error_v, fe->total_num_nodes, 3);
}


void HROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, HROM_OUTPUT_FILENAME_VTK, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	HROM_output_vtk_shape(
			&(sys->fe),
			filename,
			sys->cond.directory,
			t);
	HROM_add_output_vtk_pressure(
			&(sys->fe),
			sys->vals.p,
			sys->vals_rom.p,
            sys->vals_hrom.p,
			filename,
			sys->cond.directory,
			t);
	HROM_add_output_vtk_velocity(
			&(sys->fe),
			sys->vals.v,
			sys->vals_rom.v,
            sys->vals_hrom.v,
			filename,
			sys->cond.directory,
			t);

	double L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals.p,
			(const double*)sys->vals_rom.p);

	printf("%s L2 error pressure FEM-ROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_fem-rom.txt", sys->cond.directory);
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

	printf("%s L2 error velocity FEM-ROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_fem-rom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}

    L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals.p,
			(const double*)sys->vals_hrom.p);

	printf("%s L2 error pressure FEM-HROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_fem-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals.v,
			(const double**)sys->vals_hrom.v);

	printf("%s L2 error velocity FEM-HROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_fem-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}

    L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals_rom.p,
			(const double*)sys->vals_hrom.p);

	printf("%s L2 error pressure ROM-HROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_rom-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals_rom.v,
			(const double**)sys->vals_hrom.v);

	printf("%s L2 error velocity ROM-HROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_rom-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}
}


void HROM_pre_offline(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr0));
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();

	//for arbit dof ddecm
	HROM_ddecm_get_neib_subdomain_id(
		&(sys->mono_com),
		&(rom->hlpod_mat),
		rom->hlpod_vals.num_2nd_subdomains);

    double t1 = monolis_get_time_global_sync();

    HROM_ddecm_get_neib_num_modes(
        &(sys->mono_com_rom),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
        rom->hlpod_vals.num_modes);

    HROM_ddecm_get_neib_coordinates_pre(
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
        rom->hlpod_vals.num_modes_max);

    printf("%d %d\n", rom->hlpod_vals.num_modes_max, rom->hlpod_vals.num_modes_pre);

	//level2領域の最大基底本数の共有
    HROM_ddecm_get_neib_max_num_modes(
		&(sys->mono_com_rom),
        &(rom->hlpod_vals),
		&(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
		rom->hlpod_vals.num_modes_max);

	HROM_ecm_get_meta_neib(
		&(sys->mono_com_rom_solv),
		&(rom->hlpod_meta),
		sys->cond.directory);

	HROM_ddecm_set_neib(
		&(sys->mono_com_rom_solv),
		&(rom->hlpod_mat),
		&(hrom->hlpod_ddhr),
		&(rom->hlpod_meta),
		num_2nd_subdomains,
		rom->hlpod_vals.num_snapshot,
		sys->cond.directory);

}


void HROM_pre_offline2(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
/*
    	ddhr_set_matvec_RH_for_NNLS_para_volume_const(
        &(sys->fe),
        &(sys->vals),
        &(sys->basis),
        &(rom->hlpod_mat),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_2nd_subdomains,
        0 ,   //index 0 start
        rom->hlpod_vals.num_snapshot,
        1 + sys->mono_com.recv_n_neib,
        sys->vals.dt,
        0);
*/
    HROM_ddecm_write_selected_elems_para_arbit_subd(
        &(sys->mono_com_rom_solv),
        &(sys->fe),
        &(sys->bc),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        &(rom->hlpod_meta),
        sys->fe.total_num_elems,
        rom->hlpod_vals.num_snapshot,
        rom->hlpod_vals.num_modes_pre,
        num_2nd_subdomains,
        10000,
        1.0e-8,
        4,
        sys->cond.directory);

    HROM_ddecm_get_selected_elems_int_ovl(
            &(hrom->hlpod_ddhr),
            sys->cond.directory);

    double t_tmp = monolis_get_time_global_sync();
}


void HROM_pre_online(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();

    HROM_ddecm_get_neib_subdomain_id_2nddd(
        &(sys->mono_com),
        &(rom->hlpod_mat),
        rom->hlpod_vals.num_2nd_subdomains);

    HROM_ddecm_read_selected_elems_para(
        num_2nd_subdomains,
        sys->cond.directory);

    HROM_ddecm_get_selected_elema_add(
        &(hrom->hlpod_ddhr),
        monolis_mpi_get_global_comm_size(),
        sys->cond.directory);

    HROM_ddecm_set_podbasis_ovl(
        &(sys->mono_com),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        sys->fe.total_num_nodes,
        4);

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
        &(sys->monolis_hr0),
        rom->hlpod_meta.num_meta_nodes,
        rom->hlpod_meta.n_dof_list,
        rom->hlpod_meta.index,
        rom->hlpod_meta.item);
}

void HROM_hierarchical_parallel(
    FE_SYSTEM   sys,
    ROM*        rom,
    HROM*       hrom,
    const int   step_HR,
    const int   step_POD,
    const double t)
{
    monolis_initialize(&(sys.monolis_hr));
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_hr0), &(sys.monolis_hr));

    HROM_ddecm_set_reduced_mat_para_save_memory(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals_hrom),
        &(sys.basis),
        &(sys.bc),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    HROM_ddecm_calc_block_mat_bcsr(
        &(sys.monolis_hr),
        &(sys.mono_com_rom_solv),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_meta),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.cond.directory);
	
    HROM_ddecm_set_reduced_vec_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals_hrom),
        &(sys.basis),
        &(hrom->hr_vals),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt,
        t);
    
    HROM_ddecm_set_D_bc_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals_hrom),
        &(sys.basis),
        &(sys.bc),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    HROM_ddecm_to_monollis_rhs_para(
        &(sys.monolis_hr),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
		rom->hlpod_vals.num_2nd_subdomains,
        rom->hlpod_vals.num_modes_pre);

    monolis_show_iterlog (&(sys.monolis_hr), false);

    BBFE_sys_monowrap_solve(
        &(sys.monolis_hr),
        &(sys.mono_com_rom_solv),
        rom->hlpod_mat.mode_coef,
        MONOLIS_ITER_BICGSTAB_N128,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

    HROM_ddecm_calc_block_solution(
        &(sys.mono_com),
        &(sys.fe),
        &(hrom->hr_vals),
        &(rom->hlpod_mat),
        rom->hlpod_vals.num_2nd_subdomains,
		4);

	ROM_sys_hlpod_fe_add_Dbc(
        hrom->hr_vals.sol_vec,
		&(sys.bc),
		sys.fe.total_num_nodes,
		4);

	monolis_mpi_update_R(&(sys.mono_com), sys.fe.total_num_nodes, 4, hrom->hr_vals.sol_vec);

    BBFE_fluid_sups_renew_velocity(
        sys.vals_hrom.v, 
        hrom->hr_vals.sol_vec,
        sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys.vals_hrom.p, 
        hrom->hr_vals.sol_vec,
        sys.fe.total_num_nodes);

	output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);

}



void solver_hrom_NR(
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
        monolis_clear_mat_value_R(&(sys->monolis_hr));
        monolis_clear_mat_value_rhs_R(&(sys->monolis_rom0));
        monolis_com_initialize_by_self(&(sys->mono_com0));

        monolis_initialize(&(sys->monolis_hr));
        monolis_copy_mat_nonzero_pattern_R(&(sys->monolis_rom0), &(sys->monolis_hr));
        monolis_clear_mat_value_R(&(sys->monolis_hr));
        monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_hr));

        HROM_set_element_mat_NR(
            &(sys->monolis_hr),
            &(sys->fe),
            &(sys->vals_hrom),
            &(sys->basis),
            &(sys->bc),
            &(sys->rom_sups.hlpod_vals),
            &(sys->rom_sups.hlpod_mat),
            &(sys->hrom_sups.hlpod_ddhr),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->vals.dt);

        HROM_ddecm_calc_block_mat_bcsr(
            &(sys->monolis_hr),
            &(sys->mono_com_rom_solv),
            &(sys->rom_sups.hlpod_vals),
            &(sys->rom_sups.hlpod_mat),
            &(sys->hrom_sups.hlpod_ddhr),
            &(sys->rom_sups.hlpod_meta),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->cond.directory);	

        set_element_vec_NR_linear(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_hrom));

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
        
        ROM_std_hlpod_calc_reduced_rhs(
            &(sys->monolis),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            4);

        ROM_std_hlpod_reduced_rhs_to_monollis(
            &(sys->monolis_hr),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->rom_sups.hlpod_vals.num_modes_pre);
        
        HROM_set_element_vec_NR(
            &(sys->monolis_hr),
            &(sys->fe),
            &(sys->vals_hrom),
            &(sys->basis),
            &(sys->hrom_sups.hr_vals),
            &(sys->rom_sups.hlpod_vals),
            &(sys->hrom_sups.hlpod_ddhr),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_vals.num_modes_pre,
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            sys->vals.dt,
            t);

        ROM_monowrap_solve(
            &(sys->monolis_hr),
            &(sys->mono_com_rom_solv),
            sys->rom_sups.hlpod_mat.mode_coef,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        HROM_ddecm_calc_block_solution(
            &(sys->mono_com),
            &(sys->fe),
            &(sys->hrom_sups.hr_vals),
            &(sys->rom_sups.hlpod_mat),
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            4);

        monolis_mpi_update_R(
            &(sys->mono_com),
            sys->fe.total_num_nodes,
            4,
            sys->hrom_sups.hr_vals.sol_vec);

        BBFE_fluid_sups_renew_velocity(
            sys->vals_hrom.delta_v,
            sys->hrom_sups.hr_vals.sol_vec,
            sys->fe.total_num_nodes);

        BBFE_fluid_sups_renew_pressure(
            sys->vals_hrom.delta_p,
            sys->hrom_sups.hr_vals.sol_vec,
            sys->fe.total_num_nodes);

        update_velocity_pressure_NR(
            sys->vals_hrom.v,
            sys->vals_hrom.delta_v,
            sys->vals_hrom.p,
            sys->vals_hrom.delta_p,
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
            &(sys->vals_hrom));

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
            sys->vals_hrom.v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_delta_v = calc_internal_norm_2d(
            sys->vals_hrom.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys->mono_com.n_internal_vertex; ++ii) {
            double pv  = sys->vals_hrom.p[ii];
            double dpv = sys->vals_hrom.delta_p[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys->mono_com.n_internal_vertex; ++i_node) {
            for (int d = 0; d < 3; ++d) {
                double av = fabs(sys->vals_hrom.delta_v[i_node][d]);
                if (av > linf_delta_v_local) linf_delta_v_local = av;
            }
            double ap = fabs(sys->vals_hrom.delta_p[i_node]);
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
                    double du = fabs(sys->vals_hrom.v[ii][d] - sys->vals_hrom.v_old[ii][d]);
                    if (du > max_du) max_du = du;
                }
            }
            if(monolis_mpi_get_global_my_rank()==0){
                printf("ROM : [step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_2d(
                sys->vals_hrom.v,
                sys->vals_hrom.v_old,
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


void HROM_pre(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){
    }
    else{
        HROM_pre_offline(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
}


void HROM_memory_allocation(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
        HROM_ddecm_set_element_para(
                &(sys->fe),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                sys->cond.directory);

        HROM_ddecm_memory_allocation_para(
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                &(rom->hlpod_mat),
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains);
}

void HROM_memory_allocation_write_data(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
        HROM_ddecm_set_element_para(
                &(sys->fe),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                sys->cond.directory);

        HROM_ddecm_memory_allocation_para_pre(
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                &(rom->hlpod_mat),
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains);
}



void HROM_set_matvec(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom,
        int step,
        double t)
{
    HROM_get_neib_coordinates(
            &(sys->mono_com_rom),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            1 + sys->mono_com_rom_solv.recv_n_neib,
            rom->hlpod_vals.num_modes_max,
            rom->hlpod_vals.num_2nd_subdomains,
            rom->hlpod_vals.num_modes_pre);

    HROM_ddecm_set_RH_for_NNLS_para(
            &(sys->fe),
            &(sys->vals),
            &(sys->basis),
            &(rom->hlpod_mat),
            &(rom->hlpod_vals),
            &(hrom->hlpod_ddhr),
            rom->hlpod_vals.num_2nd_subdomains,
            step -1 ,   //index 0 start
            rom->hlpod_vals.num_snapshot,
            rom->hlpod_vals.num_modes_pre,
            sys->vals.dt,
            t);

    HROM_ddecm_set_residuals_for_NNLS_para(
            &(sys->fe),
            &(sys->vals),
            &(sys->basis),
            &(sys->bc),
            &(rom->hlpod_mat),
            &(rom->hlpod_vals),
            &(hrom->hlpod_ddhr),
            rom->hlpod_vals.num_2nd_subdomains,
            step -1 ,   //index 0 start
            rom->hlpod_vals.num_snapshot,
            1 + sys->mono_com.recv_n_neib,
            sys->vals.dt,
            t);

}


void HROM_pre_offline3(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        HROM_pre_offline2(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
    else{
        HROM_pre_offline2(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
}


void HROM_memory_allocation_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    HROM_ddecm_set_element_para(
        &(sys->fe),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_2nd_subdomains,
        sys->cond.directory);

    HROM_ddecm_memory_allocation_para_online(
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        sys->fe.total_num_nodes);

}


void HROM_pre_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){		
		HROM_pre_online(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
	}
	else{
		HROM_pre_online(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
	}
}

void HROM_std_hlpod_pre_lpod_para(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM*		 rom,
        const int    dof,
        const char*	 directory)
{
    ROM_std_hlpod_get_neib_vec(
            monolis_com,
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            rom->hlpod_vals.num_modes,
            dof);

    ROM_std_hlpod_get_neib_num_modes_para_subd(
            mono_com_rom,
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            1 + monolis_com->recv_n_neib,
            rom->hlpod_vals.num_modes);

    ROM_std_hlpod_get_neib_num_modes_mode_subd(
            mono_com_rom,
            mono_com_rom_solv,
            &(rom->hlpod_mat),
            &(rom->hlpod_meta),
            1 + monolis_com->recv_n_neib,
            directory);

    ROM_std_hlpod_get_n_dof_list(
            mono_com_rom_solv,
            &(rom->hlpod_mat),
            &(rom->hlpod_meta),
            rom->hlpod_vals.num_modes_pre);

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
            monolis_rom0,
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_meta.n_dof_list,
            rom->hlpod_meta.index,
            rom->hlpod_meta.item);
}


void HROM_std_hlpod_online_pre(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* mono_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM* 		 rom_sups,
        const int 	 total_num_nodes,
        const int 	 dof,
        const char*  metagraph,
        const char*  directory)
{
    ROM_std_hlpod_read_metagraph(
			monolis_rom0,
			mono_com_rom_solv,
			rom_sups,
			metagraph,
			directory);

    if(monolis_mpi_get_global_comm_size() == 1){
    
        if(rom_sups->hlpod_vals.num_1st_subdomains==0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else{
            
        }
    }		
    else{
        if(rom_sups->hlpod_vals.bool_global_mode==false){

            HROM_std_hlpod_pre_lpod_para(
                    monolis_rom0,
                    mono_com,
                    mono_com_rom,
                    mono_com_rom_solv,
                    rom_sups,
                    dof,
                    directory);

        }
        else{				
            ROM_std_hlpod_update_global_modes(
                    mono_com,
                    &(rom_sups->hlpod_mat),
                    total_num_nodes,
                    mono_com->n_internal_vertex,
                    rom_sups->hlpod_vals.num_modes_pre,
                    //4);
                    1);
                
        }
    }
}

void HROM_pre_offline_inc_svd1(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
    ddhr_set_matvec_RH_for_NNLS_para_volume_const(
        &(sys->fe),
        &(sys->vals),
        &(sys->basis),
        &(rom->hlpod_mat),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_2nd_subdomains,
        0 ,   //index 0 start
        rom->hlpod_vals.num_snapshot,
        1 + sys->mono_com.recv_n_neib,
        sys->vals.dt,
        0);


    HROM_ddecm_write_selected_elems_inc_svd_init_with_first_block(
        &(sys->mono_com_rom_solv),
        &(sys->fe),
        &(sys->bc),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        &(rom->hlpod_meta),
        sys->fe.total_num_elems,
        rom->hlpod_vals.num_snapshot,
        rom->hlpod_vals.num_modes_pre,
        num_2nd_subdomains,
        4,
        sys->cond.directory);

}

void HROM_pre_offline_inc_svd2(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
    HROM_ddecm_write_selected_elems_inc_svd_update(
        &(sys->mono_com_rom_solv),
        &(sys->fe),
        &(sys->bc),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        &(rom->hlpod_meta),
        sys->fe.total_num_elems,
        rom->hlpod_vals.num_snapshot,
        rom->hlpod_vals.num_modes_pre,
        num_2nd_subdomains,
        4,
        sys->cond.directory);

}

void HROM_pre_offline_inc_svd3(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
    HROM_ddecm_write_selected_elems_inc_svd(
        &(sys->mono_com_rom_solv),
        &(sys->fe),
        &(sys->bc),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        &(rom->hlpod_meta),
        sys->fe.total_num_elems,
        rom->hlpod_vals.num_snapshot,
        rom->hlpod_vals.num_modes_pre,
        num_2nd_subdomains,
        10000,
        1.0e-8,
        4,
        sys->cond.directory);

    HROM_ddecm_get_selected_elems_int_ovl(
            &(hrom->hlpod_ddhr),
            sys->cond.directory);

    double t_tmp = monolis_get_time_global_sync();

}


void read_NNLS_data(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom)
{
    HROM_ddecm_read_NNLS_data(
            &(sys->fe),
            &(sys->basis),
            &(sys->vals_rom),
            &(sys->bc),
            &(sys->rom_sups.hlpod_mat),
            &(sys->rom_sups.hlpod_vals),
            &(sys->hrom_sups.hlpod_ddhr),
            sys->rom_sups.hlpod_vals.num_2nd_subdomains,
            0,
            sys->rom_sups.hlpod_vals.num_snapshot,
            1 + sys->mono_com.recv_n_neib,
            sys->vals.dt,
            t);
}
