
#include "core_HROM.h"

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

/*for Hyper-reduction*/
void HROM_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
        HR_VALUES*      hr_vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_HR*      hlpod_hr,		
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
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hr_vals->sol_vec, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/

/*for Hyper-reduction*/
void HROM_output_result_file_vtk_para(
		BBFE_DATA*     fe,
		VALUES*        vals,
        HR_VALUES*      hr_vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_DDHR*     hlpod_ddhr,		
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
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hr_vals->sol_vec, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/



/*for Hyper-reduction*/
void HROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	if(monolis_mpi_get_global_comm_size() == 1 &&sys->rom.hlpod_vals.num_2nd_subdomains == 1){
        HROM_output_result_file_vtk(
                &(sys->fe),
                &(sys->vals),
                &(sys->hrom.hr_vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_hr),
                filename,
                sys->cond.directory,
                t);
    }
    else{  
        HROM_output_result_file_vtk_para(
                &(sys->fe),
                &(sys->vals),
                &(sys->hrom.hr_vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_ddhr),
                filename,
                sys->cond.directory,
                t);
    }

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->rom.hlpod_vals.sol_vec,
			filename,
			sys->cond.directory);


	/**** for manufactured solution ****/
/*
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_sou);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);
*/
	double L2_error_fem_hr;
	double L2_error_pod_hr;
	if(monolis_mpi_get_global_comm_size() == 1){
		if(sys->rom.hlpod_vals.num_2nd_subdomains==1){
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->hrom.hr_vals.sol_vec,
					sys->vals.T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hr_vals.sol_vec);
		}
		else{
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->vals.T,
					sys->hrom.hr_vals.sol_vec);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hr_vals.sol_vec);
		}
	}
	else{
		L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->vals.T,
				sys->hrom.hr_vals.sol_vec);
				
		L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->rom.hlpod_vals.sol_vec,
				sys->hrom.hr_vals.sol_vec);
	}

	printf("%s L2 error fem-hrom: %e\n", CODENAME, L2_error_fem_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_fem_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_fem_hr);
		fclose(fp);
	}

	printf("%s L2 error rom-hrom: %e\n", CODENAME, L2_error_pod_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_rom_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_pod_hr);
		fclose(fp);
	}

	//BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}
/********/


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
        const int    dof,
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
		&(sys->monolis_com),
		&(rom->hlpod_mat),
		rom->hlpod_vals.num_2nd_subdomains);		//num_2nd_subdomains

    printf("HROM_ddecm_get_neib_subdomain_id done\n");

    double t1 = monolis_get_time_global_sync();

    HROM_ddecm_get_neib_num_modes(
        &(sys->mono_com_rom),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->monolis_com.recv_n_neib,
        rom->hlpod_vals.num_modes);


    HROM_ddecm_get_neib_coordinates_pre(
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->monolis_com.recv_n_neib,
        rom->hlpod_vals.num_modes_max);

    printf("%d %d\n", rom->hlpod_vals.num_modes_max, rom->hlpod_vals.num_modes_pre);

	//level2領域の最大基底本数の共有
    HROM_ddecm_get_neib_max_num_modes(
		&(sys->mono_com_rom),
        &(rom->hlpod_vals),
		&(rom->hlpod_mat),
        1 + sys->monolis_com.recv_n_neib,
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
	if(monolis_mpi_get_global_comm_size() == 1){
		if(num_2nd_subdomains == 1){
			HROM_ecm_get_selected_elems(
				&(sys->fe),
				&(sys->bc),
				sys->fe.total_num_elems,
				rom->hlpod_vals.num_snapshot,
				rom->hlpod_vals.num_modes,
				10000,
				1.0e-6,
				&(hrom->hlpod_hr),
				sys->cond.directory);
            
            HROM_ecm_set_selected_elems(
				&(sys->fe),
				&(hrom->hlpod_hr),
				sys->fe.total_num_nodes,
				sys->cond.directory);
        }
        else{
    		HROM_ddecm_write_selected_elems(
				&(sys->fe),
				&(sys->bc),
				sys->fe.total_num_elems,
				rom->hlpod_vals.num_snapshot,
				rom->hlpod_vals.num_modes_pre,
				num_2nd_subdomains,
				10000,
				1.0e-6,
				&(hrom->hlpod_ddhr),
				sys->cond.directory);

        }
    }
    else{
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
            1,
            sys->cond.directory);

        HROM_ddecm_get_selected_elems_int_ovl(
                &(hrom->hlpod_ddhr),
                sys->cond.directory);

        double t_tmp = monolis_get_time_global_sync();
    }

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

    if(monolis_mpi_get_global_comm_size() == 1){
		if(num_2nd_subdomains == 1){
			monolis_initialize(&(sys->monolis_hr));

			HROM_ecm_read_selected_elems(
				&(hrom->hlpod_hr),
				num_2nd_subdomains,
				sys->cond.directory);

			HROM_ecm_set_reduced_mat(
				&(sys->monolis_hr0),
				&(sys->fe),
				&(sys->basis),
				&(sys->bc),
				&(rom->hlpod_mat),
				&(hrom->hlpod_hr),
				rom->hlpod_vals.num_modes,
				sys->vals.dt);

			HROM_ecm_monolis_set_matrix(
				&(sys->monolis_hr0),
				&(hrom->hlpod_hr),
				rom->hlpod_vals.num_modes);

		}
		else{
			monolis_initialize(&(sys->monolis_hr));        
			HROM_ddecm_get_neib_subdomain_id_nonpara(
				&(rom->hlpod_mat),
				&(hrom->hlpod_ddhr),
				num_2nd_subdomains);
        
			HROM_ddecm_read_selected_elems(
				&(hrom->hlpod_ddhr),
				num_2nd_subdomains,
				sys->cond.directory);

			HROM_ddecm_set_reduced_mat(
				&(sys->monolis_hr0),
				&(sys->fe),
				&(sys->basis),
				&(sys->bc),
				&(rom->hlpod_mat),
				&(hrom->hlpod_ddhr),
				rom->hlpod_vals.num_modes_pre,
				num_2nd_subdomains,
				sys->vals.dt);

			HROM_ddecm_monolis_set_matrix(
				&(sys->monolis_hr0),
                &(rom->hlpod_mat),
				&(hrom->hlpod_ddhr),
                &(rom->hlpod_meta),
				rom->hlpod_vals.num_modes_pre,
				num_2nd_subdomains);
		}
	}
    else{
        HROM_ddecm_get_neib_subdomain_id_2nddd(
            &(sys->monolis_com),
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
            &(sys->monolis_com),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            sys->fe.total_num_nodes,
            1);

        HROM_ddecm_set_reduced_mat_para_save_memory(
            &(sys->monolis_hr0),
            &(sys->fe),
            &(sys->basis),
            &(sys->bc),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            &(hrom->hlpod_ddhr),
            rom->hlpod_vals.num_modes_pre,
            num_2nd_subdomains,
            sys->vals.dt);

        monolis_get_nonzero_pattern_by_nodal_graph_V_R(
            &(sys->monolis_hr0),
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_meta.n_dof_list,
            rom->hlpod_meta.index,
            rom->hlpod_meta.item);

        HROM_ddecm_calc_block_mat_bcsr(
            &(sys->monolis_hr0),
            &(sys->mono_com_rom_solv),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            &(hrom->hlpod_ddhr),
            &(rom->hlpod_meta),
            rom->hlpod_vals.num_modes_pre,
            rom->hlpod_vals.num_2nd_subdomains,
            sys->cond.directory);
    }
	
}



void HROM_nonparallel(
    FE_SYSTEM       sys,
    ROM*            rom,
    HROM*           hrom,
    const int       step_HR,
    const int       step_POD,
    const double    t)
{
    double t1 = monolis_get_time();    

    monolis_copy_mat_value_R(&(sys.monolis_hr0), &(sys.monolis_hr));
	monolis_clear_mat_value_rhs_R(&(sys.monolis_hr));

    double set_bc1 = monolis_get_time();
	HROM_manusol_set_bc(
			&(sys.fe),
			&(sys.bc),
			BLOCK_SIZE,
			t,
			manusol_get_sol,
			&(sys.rom.hlpod_mat));
	double set_bc2 = monolis_get_time();

    if(sys.rom.hlpod_vals.num_2nd_subdomains == 1){

        HROM_ecm_set_reduced_vec(
            &(sys.monolis_hr),
            &(sys.fe),
            &(sys.basis),
            &(sys.hrom.hr_vals),
            &(sys.hrom.hlpod_hr),
            &(sys.rom.hlpod_mat),
            sys.rom.hlpod_vals.num_modes,
            sys.vals.dt,
            t);
        
        HROM_ecm_set_D_bc(
            &(sys.monolis_hr),
            &(sys.fe),
            &(sys.basis),
            &(sys.bc),
            &(sys.rom.hlpod_mat),
            &(sys.hrom.hlpod_hr),
            sys.rom.hlpod_vals.num_modes,
            sys.vals.dt);

        HROM_ecm_to_monollis_rhs(
            &(sys.monolis_hr),
            &(sys.hrom.hlpod_hr),
            sys.rom.hlpod_vals.num_modes);
        
        double t2 = monolis_get_time();

        BBFE_sys_monowrap_solve(
            &(sys.monolis_hr),
            &(sys.mono_com_rom_solv),
            sys.rom.hlpod_mat.mode_coef,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        double t1 = monolis_get_time();

        HROM_ecm_calc_solution(
            &(sys.fe), 
            &(sys.rom.hlpod_mat), 
            sys.hrom.hr_vals.sol_vec, 
            &(sys.bc), 
            sys.rom.hlpod_vals.num_modes);

        ROM_sys_hlpod_fe_add_Dbc(
            sys.hrom.hr_vals.sol_vec,
            &(sys.bc),
            sys.fe.total_num_nodes,
            1);

        t2 = monolis_get_time();

        output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);
    }
    else{
        double t1 = monolis_get_time();

        HROM_ddecm_set_reduced_vec(
            &(sys.monolis_hr),
            &(sys.fe),
            &(sys.basis),
            &(sys.hrom.hr_vals),
            &(sys.hrom.hlpod_ddhr),
            &(sys.rom.hlpod_mat),
            sys.rom.hlpod_vals.num_modes_pre,
            sys.rom.hlpod_vals.num_2nd_subdomains,
            sys.vals.dt,
            t);
        
        HROM_ddecm_set_D_bc(
            &(sys.monolis_hr),
            &(sys.fe),
            &(sys.basis),
            &(sys.bc),
            &(sys.rom.hlpod_mat),
            &(sys.hrom.hlpod_ddhr),
            sys.rom.hlpod_vals.num_modes_pre,
            sys.rom.hlpod_vals.num_2nd_subdomains,
            sys.vals.dt);

        HROM_ddecm_to_monollis_rhs(
            &(sys.monolis_hr),
            &(sys.rom.hlpod_mat),
            &(sys.hrom.hlpod_ddhr),
            sys.rom.hlpod_vals.num_modes_pre,
            sys.rom.hlpod_vals.num_2nd_subdomains);

        double t2 = monolis_get_time();
  
        BBFE_sys_monowrap_solve(
            &(sys.monolis_hr),
            &(sys.mono_com_rom_solv),
            sys.rom.hlpod_mat.mode_coef,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        t1 = monolis_get_time();

        HROM_ddecm_calc_solution(
            &(sys.fe), 
            &(sys.hrom.hr_vals), 
            &(sys.rom.hlpod_mat), 
            &(sys.hrom.hlpod_ddhr), 
            &(sys.bc), 
            sys.rom.hlpod_vals.num_modes_pre,
            sys.rom.hlpod_vals.num_2nd_subdomains,
            1);
        
        ROM_sys_hlpod_fe_add_Dbc(
            sys.hrom.hr_vals.sol_vec,
            &(sys.bc),
            sys.fe.total_num_nodes,
            1);

        t2 = monolis_get_time();

        output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);
    }
}

void HROM_hierarchical_parallel(
    FE_SYSTEM   sys,
    ROM*        rom,
    HROM*       hrom,
    const int   step_HR,
    const int   step_POD,
    const double t)
{
    double t1 = monolis_get_time_global_sync();

	monolis_copy_mat_value_matrix_R(&(sys.monolis_hr0), &(sys.monolis_hr));
	monolis_clear_mat_value_rhs_R(&(sys.monolis_hr));

	//double set_bc1 = monolis_get_time();
	HROM_manusol_set_bc(
        &(sys.fe),
        &(sys.bc),
        BLOCK_SIZE,
        t,
        manusol_get_sol,
        &(sys.rom.hlpod_mat));
	//double set_bc2 = monolis_get_time();
	
    HROM_ddecm_set_reduced_vec_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.basis),
        &(sys.hrom.hr_vals),
        &(sys.rom.hlpod_vals),
        &(sys.hrom.hlpod_ddhr),
        &(sys.rom.hlpod_mat),
        sys.rom.hlpod_vals.num_modes_pre,
        sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.vals.dt,
        t);
    
    HROM_ddecm_set_D_bc_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.basis),
        &(sys.bc),
        &(sys.rom.hlpod_mat),
        &(sys.hrom.hlpod_ddhr),
        sys.rom.hlpod_vals.num_modes_pre,
        sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    HROM_ddecm_to_monollis_rhs_para(
        &(sys.monolis_hr),
        &(sys.hrom.hlpod_ddhr),
        &(sys.rom.hlpod_mat),
		sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.rom.hlpod_vals.num_modes_pre);

    double t2 = monolis_get_time_global_sync();

    sys.hrom.hr_vals.time_calc_reduced_matvec = t2 - t1;

    t1 = monolis_get_time_global_sync();

    BBFE_sys_monowrap_solve(
        &(sys.monolis_rom),
        &(sys.mono_com_rom_solv),
        sys.rom.hlpod_mat.mode_coef,
        MONOLIS_ITER_CG,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

    t2 = monolis_get_time_global_sync();

    sys.hrom.hr_vals.time_linear_solver = t2 - t1;

    HROM_ddecm_calc_block_solution(
        &(sys.monolis_com),
        &(sys.fe),
        &(sys.hrom.hr_vals),
        &(sys.rom.hlpod_mat),
        sys.rom.hlpod_vals.num_2nd_subdomains,
        1);

	ROM_sys_hlpod_fe_add_Dbc(
        sys.hrom.hr_vals.sol_vec,
		&(sys.bc),
		sys.fe.total_num_nodes,
		1);
    
	monolis_mpi_update_R(&(sys.monolis_com), sys.fe.total_num_nodes, 1, sys.hrom.hr_vals.sol_vec);

	output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);

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
    if(monolis_mpi_get_global_comm_size() == 1){	
        if(rom->hlpod_vals.num_2nd_subdomains == 1){
            HROM_ecm_set_element(
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                &(hrom->hlpod_hr));
        }
        else{
            HROM_ddecm_set_element(
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                sys->cond.directory);
                
            HROM_ddecm_memory_allocation(
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes,
                rom->hlpod_vals.num_2nd_subdomains,
                &(hrom->hlpod_ddhr));
        }
    }
    else{
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
}


void HROM_set_matvec(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom,
        int step,
        double t)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_2nd_subdomains == 1){
            HROM_ecm_set_RH_for_NNLS(
                &(sys->fe),
                &(sys->basis),
                &(rom->hlpod_mat),
                &(rom->hlpod_vals), 
                &(hrom->hlpod_hr),
                step -1 ,	//index 0 start
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes,
                sys->vals.dt,
                t);
            
            HROM_ecm_set_residuals_for_NNLS(
                &(sys->fe),
                &(sys->basis),
                &(sys->bc), 
                &(rom->hlpod_mat),
                &(rom->hlpod_vals), 
                &(hrom->hlpod_hr),
                step -1 ,	//index 0 start
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes,
                sys->vals.dt,
                t);
        }
        else{
            HROM_ddecm_set_RH_for_NNLS(
                &(sys->fe),
                &(sys->basis),
                &(sys->bc), 
                &(rom->hlpod_mat),
                &(rom->hlpod_vals), 
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                step -1 ,	//index 0 start
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                sys->vals.dt,
                t);
            
            HROM_ddecm_set_residuals_for_NNLS(			
                &(sys->fe),
                &(sys->basis),
                &(rom->hlpod_mat),
                &(rom->hlpod_vals), 
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                step -1 ,	//index 0 start
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                sys->vals.dt,
                t);

        }
    }
    else{
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
                &(sys->basis),
                &(sys->bc),
                &(rom->hlpod_mat),
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                step -1 ,   //index 0 start
                rom->hlpod_vals.num_snapshot,
                1 + sys->monolis_com.recv_n_neib,
                sys->vals.dt,
                t);
    }
}


void HROM_pre_offline2(
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
	if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_2nd_subdomains == 1){
            HROM_ecm_memory_allocation_online(
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                &(hrom->hlpod_hr));
        }
        else{
            HROM_ddecm_set_element(
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                sys->cond.directory);
                
            HROM_ddecm_memory_allocation2(
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains,
                &(hrom->hlpod_ddhr));
        }
	}
    else{
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
