
#include "core_ROM.h"
#include "3ph_tr_NR.h"
#include "set_matvec.h"
#include "team21c.h"

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
			&(sys->monolis_com),
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
			&(sys->monolis_com),
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

void output_B_node_ROM(
    FE_SYSTEM*    sys,
    const double* Aphi,
    const double* Aphi_rom, const double t,
    const char* filename, const char* directory)
{
    double** B_cell = BB_std_calloc_2d_double(B_cell, sys->fe.total_num_elems, 3);
    double** B_node = BB_std_calloc_2d_double(B_node, sys->fe.total_num_nodes, 3);
    double** B_cell_rom = BB_std_calloc_2d_double(B_cell_rom, sys->fe.total_num_elems, 3);
    double** B_node_rom = BB_std_calloc_2d_double(B_node_rom, sys->fe.total_num_nodes, 3);

    compute_B_cell_average(&(sys->fe), &(sys->basis), &(sys->ned), Aphi, B_cell);
    accumulate_B_cell_to_nodes(&(sys->fe), B_cell, B_node);

    compute_B_cell_average(&(sys->fe), &(sys->basis), &(sys->ned), Aphi_rom, B_cell_rom);
    accumulate_B_cell_to_nodes(&(sys->fe), B_cell_rom, B_node_rom);

    /*
    FILE* fp = BBFE_sys_write_fopen(fp, filename, directory);
    switch (sys->fe.local_num_nodes){
        case 4: BBFE_sys_write_vtk_shape(fp, &(sys->fe), TYPE_VTK_TETRA); break;
        case 8: BBFE_sys_write_vtk_shape(fp, &(sys->fe), TYPE_VTK_HEXAHEDRON); break;
    }
    fprintf(fp, "POINT_DATA %d\n", sys->fe.total_num_nodes);
    BB_vtk_write_point_vals_vector(fp, B_node, sys->fe.total_num_nodes, "B_node");
    BB_vtk_write_point_vals_vector(fp, B_node_rom, sys->fe.total_num_nodes, "B_node_rom");
    fclose(fp);
    */

    double L2_error = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis_com),
			t,
			(const double**)B_node_rom,
			(const double**)B_node);

	printf("\n%s L2 error: %e\n", CODENAME, L2_error);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_B.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error);
		fclose(fp);
	}

    BB_std_free_2d_double(B_cell, sys->fe.total_num_elems, 3);
    BB_std_free_2d_double(B_node, sys->fe.total_num_nodes, 3);
    BB_std_free_2d_double(B_cell_rom, sys->fe.total_num_elems, 3);
    BB_std_free_2d_double(B_node_rom, sys->fe.total_num_nodes, 3);
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
/*
    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
            monolis_rom0,
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_meta.n_dof_list,
            rom->hlpod_meta.index,
            rom->hlpod_meta.item);
*/
        monolis_get_nonzero_pattern_by_nodal_graph_R(
            monolis_rom0,
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_vals.num_modes_pre,
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
                    // /4);
                    1);
                
        }
    }
}



void solver_rom_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    //double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        monolis_clear_mat_value_R(&(sys.monolis_rom));
        monolis_com_initialize_by_self(&(sys.mono_com0));

        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            //dx[i] = 0.0;
        }

        /* 組み立て */
        set_element_mat_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        ROM_std_hlpod_calc_reduced_mat(
            &(sys.monolis),
            &(sys.monolis_rom),
            &(sys.monolis_com),
            &(sys.mono_com0),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            1);

        ROM_std_hlpod_solve_ROM(
            &(sys.monolis),
            &(sys.monolis_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            1,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

        monolis_mpi_update_R(
            &(sys.monolis_com),
            sys.fe.total_num_nodes,
            1,
            sys.rom_sups.hlpod_vals.sol_vec);

        /* Newton update */
        //update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);
        update_Aphi_NR(x_curr, sys.rom_sups.hlpod_vals.sol_vec, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, A_delta, phi_delta, sys.fe.total_num_elems);

        char fnode[BUFFER_SIZE];
        const char* fn1;
        snprintf(fnode, BUFFER_SIZE, "B_node_rom_%06d.vtk", step);
        fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
        output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, fn1, sys.cond.directory);


        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        //int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        //int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        //if (conv_v && conv_p) {
        if(it == max_iter_NR -1 || (conv_v && conv_p)) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }

    //free(dx);
    free(rvec);
    free(rvec_old);


    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
}


void allreduce_global_para(
    MONOLIS_COM*	monolis_com,
    double**        mat,
    double*         rhs,
    const int		num_modes)
{
    double* vec;
	vec = BB_std_calloc_1d_double(vec, num_modes);

	for(int k2 = 0; k2 < num_modes; k2++){
        vec[k2] = rhs[k2];
	}

	monolis_allreduce_R(
		num_modes,
		vec,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k2 = 0; k2 < num_modes; k2++){
        rhs[k2] = vec[k2];
	}


	double* matrix_buffer = BB_std_calloc_1d_double(matrix_buffer, num_modes * num_modes);
	
	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            matrix_buffer[k1 * num_modes + k2] = mat[k1][k2];
		}
	}

	monolis_allreduce_R(
		num_modes * num_modes,
		matrix_buffer,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            mat[k1][k2] = matrix_buffer[k1 * num_modes + k2];
		}
	}
}


void solver_rom_global_para_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    //double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    int total_num_nodes = sys.fe.total_num_nodes;
    double** mat = BB_std_calloc_2d_double(mat, sys.rom_sups.hlpod_vals.num_modes, sys.rom_sups.hlpod_vals.num_modes);
	double* rhs = BB_std_calloc_1d_double(rhs, sys.rom_sups.hlpod_vals.num_modes);
    double* mode_coef = BB_std_calloc_1d_double(mode_coef, sys.rom_sups.hlpod_vals.num_modes);
    double* ansvec = BB_std_calloc_1d_double(ansvec, total_num_nodes);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        monolis_clear_mat_value_R(&(sys.monolis_rom));
        monolis_com_initialize_by_self(&(sys.mono_com0));

        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        set_reduced_vec_global_para(
                &(sys.monolis),
                &(sys.monolis_com),
                &(sys.fe),
                &(sys.basis),
                &(sys.vals_rom),
                &(sys.ned),
                rhs,
                sys.rom_sups.hlpod_mat.pod_modes,
                sys.rom_sups.hlpod_vals.num_modes,
                x_prev,
                x_curr,
                sys.vals.dt,
                t,
                sys.cond.directory);

        set_reduced_mat_global_para(
                &(sys.monolis),
                &(sys.monolis_com),
                &(sys.fe),
                &(sys.basis),
                &(sys.vals_rom),
                &(sys.ned),
                mat,
                sys.rom_sups.hlpod_mat.pod_modes,
                sys.rom_sups.hlpod_vals.num_modes,
                x_curr,
                sys.vals.dt,
                sys.cond.directory);

        allreduce_global_para(
                &(sys.monolis_com),
                mat,
                rhs,
                sys.rom_sups.hlpod_vals.num_modes);

        ROM_BB_gauss_elimination(
                sys.rom_sups.hlpod_vals.num_modes,
                mat,
                rhs,
                mode_coef);
            
        ROM_std_hlpod_calc_sol_global_para(
                ansvec,
                sys.rom_sups.hlpod_mat.pod_modes,
                mode_coef,
                sys.fe.total_num_nodes,
                sys.rom_sups.hlpod_vals.num_modes,
                1);
        
        monolis_mpi_update_R(&(sys.monolis_com), sys.fe.total_num_nodes, 1, ansvec);

        /* Newton update */
        //update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);
        update_Aphi_NR(x_curr, ansvec, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), ansvec, A_delta, phi_delta, sys.fe.total_num_elems);

        char fnode[BUFFER_SIZE];
        const char* fn1;
        snprintf(fnode, BUFFER_SIZE, "B_node_rom_%06d.vtk", step);
        fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
        output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), ansvec, fn1, sys.cond.directory);


        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        //int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        //int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        //if (conv_v && conv_p) {
        if(it == max_iter_NR -1 || (conv_v && conv_p)) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }

    //free(dx);
    free(rvec);
    free(rvec_old);

    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
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
        &(sys->monolis_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        4);

}


void solver_rom_NR_Aphi_team21a2(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    //double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 1;


        //debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        monolis_clear_mat_value_R(&(sys.monolis_rom));
        monolis_com_initialize_by_self(&(sys.mono_com0));

        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            //dx[i] = 0.0;
        }

        /* 組み立て */
        /*
	set_element_mat_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);

        
        set_element_mat_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        */

        monolis_copy_mat_value_R(&(sys.monolis_rom0), &(sys.monolis_rom));
        

	set_element_vec_NR_Aphi_team21c(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi_team21c(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
//        if(it==0){
//            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
 //       }
/*
        ROM_std_hlpod_calc_reduced_mat(
            &(sys.monolis),
            &(sys.monolis_rom),
            &(sys.monolis_com),
            &(sys.mono_com0),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            1);
*/
        ROM_std_hlpod_solve_ROM(
            &(sys.monolis),
            &(sys.monolis_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom_sups),
            sys.fe.total_num_nodes,
            1,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_DIAG);

        monolis_mpi_update_R(
            &(sys.monolis_com),
            sys.fe.total_num_nodes,
            1,
            sys.rom_sups.hlpod_vals.sol_vec);

        /* Newton update */
        //update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);
        update_Aphi_NR(x_curr, sys.rom_sups.hlpod_vals.sol_vec, n_dof_total, relaxation);



        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, A_delta, phi_delta, sys.fe.total_num_elems);

        char fnode[BUFFER_SIZE];
        const char* fn1;
        snprintf(fnode, BUFFER_SIZE, "B_node_rom_%06d.vtk", step);
        fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
        output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, fn1, sys.cond.directory);


            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);


    //free(dx);
    free(rvec);
    free(rvec_old);


    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
}


void calc_reduced_mat_linear_team21a2(
    FE_SYSTEM* sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total)
{
    monolis_clear_mat_value_R(&(sys->monolis));
    monolis_clear_mat_value_R(&(sys->monolis_rom));
    monolis_clear_mat_value_R(&(sys->monolis_rom0));
    monolis_com_initialize_by_self(&(sys->mono_com0));

    double t1 = monolis_get_time_global_sync();
    printf("test1");

    set_element_mat_NR_Aphi_team21c(&(sys->monolis), &(sys->fe), &(sys->basis), &(sys->ned),
                            x_curr, sys->vals.dt);

        double t2 = monolis_get_time_global_sync();
    printf("test2");

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        1,
        &(sys->bc_NR),
        sys->monolis.mat.R.B);

    double t3 = monolis_get_time_global_sync();
    printf("test3");

    ROM_std_hlpod_calc_reduced_mat(
        &(sys->monolis),
        &(sys->monolis_rom0),
        &(sys->monolis_com),
        &(sys->mono_com0),
        &(sys->mono_com_rom_solv),
        &(sys->rom_sups),
        sys->fe.total_num_nodes,
        1);

}

void ROM_std_hlpod_calc_reduced_rhs3(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    HLPOD_MAT* hlpod_mat,
    const int max_num_bases,
    const int num_2nddd,
    const int dof)
{
    int total_modes = 0;
    for(int k = 0; k < num_2nddd; k++){
        total_modes += hlpod_mat->num_modes_internal[k];
    }

    const int ndof_total = fe->total_num_nodes * dof;

    double** mat = BB_std_calloc_2d_double(mat, total_modes, total_modes);
    double* vec  = BB_std_calloc_1d_double(vec, ndof_total);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] = 0.0;
    }

    for(int r = 0; r < ndof_total; r++){
        monolis->mat.R.B[r] = 0.0;
    }

    set_element_vec_NR_Aphi_team21c_source(
        monolis, fe, basis, ned,
        x_prev, x_curr, dt, current_time);

    int index = 0;
    int index_column = 0;
    int sum = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            int a = index + i;

            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
                    int row = hlpod_mat->node_id[j + sum] * dof + l;

                    hlpod_mat->VTf[a] +=
                        hlpod_mat->pod_modes[row][index_column + i]
                      * monolis->mat.R.B[row];
                }
            }
        }

        index_column += hlpod_mat->num_modes_internal[k];
        index        += hlpod_mat->num_modes_internal[k];
        sum          += hlpod_mat->n_internal_vertex_subd[k];
    }

    int col_b = 0;

    for(int kb = 0; kb < num_2nddd; kb++){
        for(int ib = 0; ib < hlpod_mat->num_modes_internal[kb]; ib++){

            int b = col_b + ib;

            for(int r = 0; r < ndof_total; r++){
                vec[r] = 0.0;
                monolis->mat.R.B[r] = 0.0;
            }

            /*
             * 修正点:
             * モード φ_b は所属サブドメイン kb の内部自由度だけに展開する。
             * 全サブドメインへ vec を入れると、局所モードの支持が壊れる。
             */
            int sum_b = 0;
            for(int k = 0; k < kb; k++){
                sum_b += hlpod_mat->n_internal_vertex_subd[k];
            }

            int mode_col_b = col_b + ib;  /* == b */

            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[kb]; j++){
                for(int l = 0; l < dof; l++){
                    int row = hlpod_mat->node_id[j + sum_b] * dof + l;
                    vec[row] = hlpod_mat->pod_modes[row][mode_col_b];
                }
            }

            set_element_vec_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned,
                vec, vec, dt, current_time);

            int col_a = 0;
            int row_offset = 0;

            for(int ka = 0; ka < num_2nddd; ka++){
                for(int ia = 0; ia < hlpod_mat->num_modes_internal[ka]; ia++){

                    int a = col_a + ia;
                    mat[a][b] = 0.0;

                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[ka]; j++){
                        for(int l = 0; l < dof; l++){
                            int row = hlpod_mat->node_id[j + row_offset] * dof + l;

                            mat[a][b] +=
                                hlpod_mat->pod_modes[row][a]
                              * monolis->mat.R.B[row];
                        }
                    }
                }

                col_a      += hlpod_mat->num_modes_internal[ka];
                row_offset += hlpod_mat->n_internal_vertex_subd[ka];
            }
        }

        col_b += hlpod_mat->num_modes_internal[kb];
    }

    for(int a = 0; a < total_modes; a++){
        for(int b = 0; b < total_modes; b++){
            hlpod_mat->VTf[a] += mat[a][b] * hlpod_mat->mode_coef[b];
        }
    }

    BB_std_free_2d_double(mat, total_modes, total_modes);
    BB_std_free_1d_double(vec, ndof_total);
}


void ROM_std_hlpod_calc_reduced_rhs4(
    MONOLIS* monolis,
    MONOLIS* monolis_rom,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* monolis_com_solv,
    MONOLIS_COM* mono_com0,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_META* hlpod_meta,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    HLPOD_MAT* hlpod_mat,
    const int num_2nddd,
    const int       n_neib_vec,
    const int       num_modes,
    const int       max_num_bases,
    const int       total_num_bases,
    const int dof)
{
    int total_modes = 0;
    for(int k = 0; k < num_2nddd; k++){
        total_modes += hlpod_mat->num_modes_internal[k];
    }

    printf("num_modes = %d, total_modes = %d\n", num_modes, total_modes);

    const int ndof_total = fe->total_num_nodes * dof;

    double* vec  = BB_std_calloc_1d_double(vec, ndof_total);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] = 0.0;
    }

    for(int r = 0; r < ndof_total; r++){
        monolis->mat.R.B[r] = 0.0;
    }

    set_element_vec_NR_Aphi_team21c_source(
        monolis, fe, basis, ned,
        x_prev, x_curr, dt, current_time);

    int index = 0;
    int index_column = 0;
    int sum = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            int a = index + i;

            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
                    int row = hlpod_mat->node_id[j + sum] * dof + l;

                    hlpod_mat->VTf[a] +=
                        hlpod_mat->pod_modes[row][index_column + i]
                      * monolis->mat.R.B[row];
                }
            }
        }

        index_column += hlpod_mat->num_modes_internal[k];
        index        += hlpod_mat->num_modes_internal[k];
        sum          += hlpod_mat->n_internal_vertex_subd[k];
    }

    printf("done Vtf source term\n");

    const int NDOF  = fe->total_num_nodes * dof;

    double* monolis_in;  double* monolis_out;  double* monolis_in2;
    monolis_in = BB_std_calloc_1d_double(monolis_in, NDOF);
    monolis_in2 = BB_std_calloc_1d_double(monolis_in2, NDOF);
    monolis_out = BB_std_calloc_1d_double(monolis_out, NDOF);

    double** mat = BB_std_calloc_2d_double(mat, num_modes, n_neib_vec);

	double t = monolis_get_time_global_sync();

    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
        if(l < num_modes){
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
		monolis->mat.R.B[k] = 0.0;
            }

            for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
                monolis_in[j] = hlpod_mat->pod_modes[j][l];
            }

            set_element_vec_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned,
                monolis_in, monolis_in, dt, current_time);

            //monolis_matvec_product_R(monolis, mono_com0, monolis->mat.R.B, monolis_out);

            int index_row = 0;
            int sum = 0;
            int index_column = 0;
            for(int k = 0; k < num_2nddd; k++){
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                        for(int m = 0; m < dof; m++){
                            index_row = hlpod_mat->node_id[j + sum] * dof + m;
                            mat[index_column + i][l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis->mat.R.B[index_row];
                        }
                    }
                }
                index_column += hlpod_mat->num_modes_internal[k];
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
	    printf("done curl curl mat term1\n");
        }
	else{
	}
        if(l < num_modes){
            monolis_mpi_update_R( monolis_com, fe->total_num_nodes, dof, monolis_in);
        }
        else{
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
                monolis->mat.R.B[k] = 0;
            }
            monolis_mpi_update_R( monolis_com, fe->total_num_nodes, dof, monolis_in);
        }
        int index_column2 = num_modes;
        for(int n = 0; n <  monolis_com->recv_n_neib; n++){
            int iS =  monolis_com->recv_index[n];
            int iE =  monolis_com->recv_index[n + 1];

            for(int k = 0; k < NDOF; k++){
                monolis_in2[k] = 0;
		monolis->mat.R.B[k] = 0.0;
            }

            for(int k = iS; k < iE; k++){
                for(int m = 0; m < dof; m++){
                    int index =  monolis_com->recv_item[k] * dof + m;
                    monolis_in2[index] = monolis_in[index];
                }
            }

            set_element_vec_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned,
                monolis_in2, monolis_in2, dt, current_time);

            //monolis_matvec_product_R(monolis,  mono_com0, monolis_in2, monolis_out);

            int sum = 0;
            int index_row = 0;
            int index_column = 0;

            for(int k = 0; k < num_2nddd; k++){
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                        for(int m = 0; m < dof; m++){

                            index_row = hlpod_mat->node_id[j + sum] * dof + m;
                            if(l < hlpod_mat->num_modes_2nddd[n + 1]){
                                mat[index_column + i][index_column2 + l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis->mat.R.B[index_row];
                            }
                        }
                    }
                }
                index_column += hlpod_mat->num_modes_internal[k];
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
            index_column2 += hlpod_mat->num_modes_2nddd[n + 1];
        }
	printf("done curl curl mat term2\n");
    }

    printf("done curl curl mat term\n");

    BB_std_free_1d_double(monolis_in, NDOF);
    BB_std_free_1d_double(monolis_out, NDOF);
    BB_std_free_1d_double(monolis_in2, NDOF);


    const int M = max_num_bases;
    //const int n_neib_vec = hlpod_vals->n_neib_vec;
    int rank = monolis_mpi_get_global_my_rank();

    double** mat2;
    mat2 = BB_std_calloc_2d_double(mat2, M, M);

    int index1 = 0;
    int index2 = 0;
    int num_modes2 = 0;

    /* 対角ブロック */
    for(int k = 0; k < monolis_com_solv->n_internal_vertex; k++){
        int iS = hlpod_mat->num_modes_1stdd[k];
        int iE = hlpod_mat->num_modes_1stdd[k + 1];

        double frob_sq = 0.0;

        index1 = 0;
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                double val = mat[m][n];
                mat2[index1][index2] = val;

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis_rom,
                    k,
                    k,
                    index1,
                    index2,
                    val);

                index2++;
            }
            index1++;
        }
    }

    /* 非対角ブロック */
    for(int k = 0; k < monolis_com_solv->n_internal_vertex; k++){

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com_solv->n_internal_vertex; j++){

                if(hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){

                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k + 1];

                    index1 = 0;
                    for(int m = IS; m < IE; m++){

                        int IIS = hlpod_mat->num_modes_1stdd[j];
                        int IIE = hlpod_mat->num_modes_1stdd[j + 1];
                        num_modes2 += IIE - IIS;

                        index2 = 0;
                        for(int n = IIS; n < IIE; n++){
                            double val = mat[m][n];
                            mat2[index1][index2] = val;

                            //frob_sq += val * val;

                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis_rom,
                                k,
                                hlpod_meta->item[i],
                                index1,
                                index2,
                                val);

                            index2++;
                        }
                        index1++;
                    }
                }
            }
        }
    }

    printf("done setting curl curl mat term\n");

    BB_std_free_2d_double(mat, total_num_bases, n_neib_vec);
    BB_std_free_2d_double(mat2, M, M);

    double* monolis_out2 = BB_std_calloc_1d_double(monolis_out2, n_neib_vec);
    monolis_matvec_product_R(monolis_rom, mono_com0, hlpod_mat->mode_coef_old, monolis_out2);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] += monolis_out2[a];
    }

    //BB_std_free_2d_double(mat, total_modes, total_modes);
    BB_std_free_1d_double(vec, ndof_total);
}


void ROM_std_hlpod_calc_reduced_rhs5(
    MONOLIS* monolis,
    MONOLIS* monolis_rom,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* monolis_com_solv,
    MONOLIS_COM* mono_com0,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_META* hlpod_meta,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    HLPOD_MAT* hlpod_mat,
    const int num_2nddd,
    const int       n_neib_vec,
    const int       num_modes,
    const int       max_num_bases,
    const int       total_num_bases,
    const int dof)
{
    int total_modes = 0;
    for(int k = 0; k < num_2nddd; k++){
        total_modes += hlpod_mat->num_modes_internal[k];
    }

    printf("num_modes = %d, total_modes = %d\n", num_modes, total_modes);

    const int ndof_total = fe->total_num_nodes * dof;

    double* vec  = BB_std_calloc_1d_double(vec, ndof_total);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] = 0.0;
        hlpod_mat->VTf_tmp[a] = 0.0;
    }

    for(int r = 0; r < ndof_total; r++){
        monolis->mat.R.B[r] = 0.0;
    }

    set_element_vec_NR_Aphi_team21c_source(
        monolis, fe, basis, ned,
        x_prev, x_curr, dt, current_time);

    const double FREQ_HZ_team21c = 50.0;
    const double omega = 2.0 * M_PI * FREQ_HZ_team21c;

    int index = 0;
    int index_column = 0;
    int sum = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
            int a = index + i;

            for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                for(int l = 0; l < dof; l++){
                    int row = hlpod_mat->node_id[j + sum] * dof + l;

                    hlpod_mat->VTf[a] +=
                        hlpod_mat->pod_modes[row][index_column + i]
                      * monolis->mat.R.B[row]  * sin(omega * current_time);

                    hlpod_mat->VTf_tmp[a] += hlpod_mat->pod_modes[row][index_column + i]
                      * monolis->mat.R.B[row];
                }
            }
        }

        index_column += hlpod_mat->num_modes_internal[k];
        index        += hlpod_mat->num_modes_internal[k];
        sum          += hlpod_mat->n_internal_vertex_subd[k];
    }

    printf("done Vtf source term\n");

    const int NDOF  = fe->total_num_nodes * dof;

    double* monolis_in;  double* monolis_out;  double* monolis_in2;
    monolis_in = BB_std_calloc_1d_double(monolis_in, NDOF);
    monolis_in2 = BB_std_calloc_1d_double(monolis_in2, NDOF);
    monolis_out = BB_std_calloc_1d_double(monolis_out, NDOF);

    double** mat = BB_std_calloc_2d_double(mat, num_modes, n_neib_vec);
    double t = monolis_get_time_global_sync();


    set_element_mat_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned, current_time);

    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
        if(l < num_modes){
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
                monolis->mat.R.B[k] = 0.0;
            }

            for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
                monolis_in[j] = hlpod_mat->pod_modes[j][l];
            }

            /*
            set_element_vec_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned,
                monolis_in, monolis_in, dt, current_time);
            */

            monolis_matvec_product_R(monolis, mono_com0, monolis_in, monolis_out);

            int index_row = 0;
            int sum = 0;
            int index_column = 0;
            for(int k = 0; k < num_2nddd; k++){
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                        for(int m = 0; m < dof; m++){
                            index_row = hlpod_mat->node_id[j + sum] * dof + m;
                            mat[index_column + i][l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis_out[index_row];
                        }
                    }
                }
                index_column += hlpod_mat->num_modes_internal[k];
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
            printf("done curl curl mat term1\n");
        }
        else{
        }
        if(l < num_modes){
            monolis_mpi_update_R( monolis_com, fe->total_num_nodes, dof, monolis_in);
        }
        else{
            for(int k = 0; k < NDOF; k++){
                monolis_in[k] = 0;
                monolis->mat.R.B[k] = 0;
            }
            monolis_mpi_update_R( monolis_com, fe->total_num_nodes, dof, monolis_in);
        }
        int index_column2 = num_modes;
        for(int n = 0; n <  monolis_com->recv_n_neib; n++){
            int iS =  monolis_com->recv_index[n];
            int iE =  monolis_com->recv_index[n + 1];

            for(int k = 0; k < NDOF; k++){
                monolis_in2[k] = 0;
                monolis->mat.R.B[k] = 0.0;
            }

            for(int k = iS; k < iE; k++){
                for(int m = 0; m < dof; m++){
                    int index =  monolis_com->recv_item[k] * dof + m;
                    monolis_in2[index] = monolis_in[index];
                }
            }
/*
            set_element_vec_NR_Aphi_team21c_curl_curl(
                monolis, fe, basis, ned,
                monolis_in2, monolis_in2, dt, current_time);
*/
            monolis_matvec_product_R(monolis,  mono_com0, monolis_in2, monolis_out);

            int sum = 0;
            int index_row = 0;
            int index_column = 0;

            for(int k = 0; k < num_2nddd; k++){
                for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                    for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
                        for(int m = 0; m < dof; m++){

                            index_row = hlpod_mat->node_id[j + sum] * dof + m;
                            if(l < hlpod_mat->num_modes_2nddd[n + 1]){
                                mat[index_column + i][index_column2 + l] += hlpod_mat->pod_modes[index_row][index_column + i] * monolis_out[index_row];
                            }
                        }
                    }
                }
                index_column += hlpod_mat->num_modes_internal[k];
                sum += hlpod_mat->n_internal_vertex_subd[k];
            }
            index_column2 += hlpod_mat->num_modes_2nddd[n + 1];
        }
        printf("done curl curl mat term2\n");
    }

    printf("done curl curl mat term\n");

    BB_std_free_1d_double(monolis_in, NDOF);
    BB_std_free_1d_double(monolis_out, NDOF);
    BB_std_free_1d_double(monolis_in2, NDOF);


    const int M = max_num_bases;
    //const int n_neib_vec = hlpod_vals->n_neib_vec;
    int rank = monolis_mpi_get_global_my_rank();

    double** mat2;
    mat2 = BB_std_calloc_2d_double(mat2, M, M);

    int index1 = 0;
    int index2 = 0;
    int num_modes2 = 0;

    /* 対角ブロック */
    for(int k = 0; k < monolis_com_solv->n_internal_vertex; k++){
        int iS = hlpod_mat->num_modes_1stdd[k];
        int iE = hlpod_mat->num_modes_1stdd[k + 1];

        double frob_sq = 0.0;

        index1 = 0;
        for(int m = iS; m < iE; m++){
            index2 = 0;
            for(int n = iS; n < iE; n++){
                double val = mat[m][n];
                mat2[index1][index2] = val;

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis_rom,
                    k,
                    k,
                    index1,
                    index2,
                    val);

                index2++;
            }
            index1++;
        }
    }

    /* 非対角ブロック */
    for(int k = 0; k < monolis_com_solv->n_internal_vertex; k++){

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com_solv->n_internal_vertex; j++){

                if(hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){

                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k + 1];

                    index1 = 0;
                    for(int m = IS; m < IE; m++){

                        int IIS = hlpod_mat->num_modes_1stdd[j];
                        int IIE = hlpod_mat->num_modes_1stdd[j + 1];
                        num_modes2 += IIE - IIS;

                        index2 = 0;
                        for(int n = IIS; n < IIE; n++){
                            double val = mat[m][n];
                            mat2[index1][index2] = val;

                            //frob_sq += val * val;

                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis_rom,
                                k,
                                hlpod_meta->item[i],
                                index1,
                                index2,
                                val);

                            index2++;
                        }
                        index1++;
                    }
                }
            }
        }
    }

    printf("done setting curl curl mat term\n");

    BB_std_free_2d_double(mat, total_num_bases, n_neib_vec);
    BB_std_free_2d_double(mat2, M, M);

    double* monolis_out2 = BB_std_calloc_1d_double(monolis_out2, n_neib_vec);
    monolis_matvec_product_R(monolis_rom, mono_com0, hlpod_mat->mode_coef_old, monolis_out2);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] -= monolis_out2[a];
    }

    //BB_std_free_2d_double(mat, total_modes, total_modes);
    BB_std_free_1d_double(vec, ndof_total);
}


void ROM_std_hlpod_calc_reduced_rhs6(
    MONOLIS* monolis,
    MONOLIS* monolis_rom,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* monolis_com_solv,
    MONOLIS_COM* mono_com0,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_META* hlpod_meta,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    HLPOD_MAT* hlpod_mat,
    const int num_2nddd,
    const int       n_neib_vec,
    const int       num_modes,
    const int       max_num_bases,
    const int       total_num_bases,
    const int dof)
{
    int total_modes = 0;
    for(int k = 0; k < num_2nddd; k++){
        total_modes += hlpod_mat->num_modes_internal[k];
    }
    printf("num_modes = %d, total_modes = %d\n", num_modes, total_modes);

    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] = 0.0;
    }

    const double FREQ_HZ_team21c = 50.0;
    const double omega = 2.0 * M_PI * FREQ_HZ_team21c;

    int index = 0;
    int index_column = 0;
    int sum = 0;

    for(int k = 0; k < num_2nddd; k++){
        for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
                int a = index + i;
                hlpod_mat->VTf[a] += hlpod_mat->VTf_tmp[a] * sin(omega * current_time);
            }
        index_column += hlpod_mat->num_modes_internal[k];
        index        += hlpod_mat->num_modes_internal[k];
        sum          += hlpod_mat->n_internal_vertex_subd[k];
    }

    double* monolis_out2 = BB_std_calloc_1d_double(monolis_out2, n_neib_vec);
    monolis_matvec_product_R(monolis_rom, mono_com0, hlpod_mat->mode_coef_old, monolis_out2);
    for(int a = 0; a < total_modes; a++){
        hlpod_mat->VTf[a] -= monolis_out2[a];
    }

    //BB_std_free_2d_double(mat, total_modes, total_modes);
    //BB_std_free_1d_double(vec, ndof_total);
}



void solver_rom_NR_Aphi_team21a2_fast(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    //double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 1;


        //debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);
/*
        monolis_clear_mat_value_R(&(sys.monolis));
        monolis_clear_mat_value_R(&(sys.monolis_rom));
	//monolis_clear_mat_value_R(&(sys.monolis_rom_mass));
        monolis_com_initialize_by_self(&(sys.mono_com0));

        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            //dx[i] = 0.0;
        }
*/
        //monolis_copy_mat_value_R(&(sys.monolis_rom0), &(sys.monolis_rom));
/*
        ROM_std_hlpod_calc_reduced_rhs3(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t, &(sys.rom_sups.hlpod_mat),
                                sys.rom_sups.hlpod_vals.num_modes_pre,
			       sys.rom_sups.hlpod_vals.num_2nd_subdomains,	1);
*/

	/*
        ROM_std_hlpod_calc_reduced_rhs4(&(sys.monolis), &(sys.monolis_rom_mass), &(sys.monolis_comm), &(sys.mono_com_rom_solv),
                                &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t, &(sys.rom_sups.hlpod_mat),
                                sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                                hlpod_vals.n_neib_vec,
                                rom->hlpod_vals.num_modes,
                                rom->hlpod_vals.num_modes_pre,
                                rom->hlpod_vals.num_modes, 1);
*/
	if(step == 0){
		monolis_clear_mat_value_R(&(sys.monolis));
		monolis_clear_mat_value_R(&(sys.monolis_rom_mass));
		monolis_copy_mat_value_R(&(sys.monolis_rom0), &(sys.monolis_rom));
		monolis_com_initialize_by_self(&(sys.mono_com0));

        	ROM_std_hlpod_calc_reduced_rhs5(&(sys.monolis), &(sys.monolis_rom_mass), &(sys.monolis_com), &(sys.mono_com_rom_solv),&(sys.mono_com0),
                                &(sys.fe), &(sys.basis), &(sys.rom_sups.hlpod_vals), &(sys.rom_sups.hlpod_meta), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t, &(sys.rom_sups.hlpod_mat),
                                sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                                sys.rom_sups.hlpod_vals.n_neib_vec,
                                sys.rom_sups.hlpod_vals.num_modes,
                                sys.rom_sups.hlpod_vals.num_modes_pre,
                                sys.rom_sups.hlpod_vals.num_modes, 1);
	}
	else{
		monolis_copy_mat_value_R(&(sys.monolis_rom0), &(sys.monolis_rom));
        	ROM_std_hlpod_calc_reduced_rhs6(&(sys.monolis), &(sys.monolis_rom_mass), &(sys.monolis_com), &(sys.mono_com_rom_solv),&(sys.mono_com0),
                                &(sys.fe), &(sys.basis), &(sys.rom_sups.hlpod_vals), &(sys.rom_sups.hlpod_meta), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t, &(sys.rom_sups.hlpod_mat),
                                sys.rom_sups.hlpod_vals.num_2nd_subdomains,
                                sys.rom_sups.hlpod_vals.n_neib_vec,
                                sys.rom_sups.hlpod_vals.num_modes,
                                sys.rom_sups.hlpod_vals.num_modes_pre,
                                sys.rom_sups.hlpod_vals.num_modes, 1);
	}

        ROM_std_hlpod_reduced_rhs_to_monollis(
            &(sys.monolis_rom),
            &(sys.rom_sups.hlpod_mat),
            sys.rom_sups.hlpod_vals.num_2nd_subdomains,
            sys.rom_sups.hlpod_vals.num_modes_pre);

        ROM_monowrap_solve(
            &(sys.monolis_rom),
            &(sys.mono_com_rom_solv),
            sys.rom_sups.hlpod_mat.mode_coef,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);


	for(int i = 0; i < sys.rom_sups.hlpod_vals.n_neib_vec; i++){
		sys.rom_sups.hlpod_mat.mode_coef_old[i] += sys.rom_sups.hlpod_mat.mode_coef[i];
	}

	/*
        HROM_ddecm_calc_block_solution(
            &(sys.mono_com),
            &(sys.fe),
            &(sys.hrom_sups.hr_vals),
            &(sys.rom_sups.hlpod_mat),
            sys.rom_sups.hlpod_vals.num_2nd_subdomains,
            4);
*/

if(step >= 1){

	for(int i = 0; i < sys.rom_sups.hlpod_vals.n_neib_vec; i++){
                sys.rom_sups.hlpod_mat.mode_coef[i] = sys.rom_sups.hlpod_mat.mode_coef_old[i];
        }

    ROM_std_hlpod_calc_sol(
        &(sys.rom_sups.hlpod_vals),
        &(sys.rom_sups.hlpod_mat),
        n_dof_total,
        sys.rom_sups.hlpod_vals.num_modes_pre,
        sys.rom_sups.hlpod_vals.num_2nd_subdomains,
        1);


        monolis_mpi_update_R(
            &(sys.monolis_com),
            sys.fe.total_num_nodes,
            1,
            sys.rom_sups.hlpod_vals.sol_vec);

	for(int i = 0; i < sys.fe.total_num_nodes; i++){
                x_curr[i] = sys.rom_sups.hlpod_vals.sol_vec[i];
        }
	copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, A, phi, sys.fe.total_num_elems);


        /* Newton update */
        //update_Aphi_NR(x_curr, sys.rom_sups.hlpod_vals.sol_vec, n_dof_total, relaxation);


        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        //copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        //copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, A_delta, phi_delta, sys.fe.total_num_elems);

        char fnode[BUFFER_SIZE];
        const char* fn1;
        snprintf(fnode, BUFFER_SIZE, "B_node_rom_%06d.vtk", step);
        fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
        output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), sys.rom_sups.hlpod_vals.sol_vec, fn1, sys.cond.directory);

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);
}

    //free(dx);
    free(rvec);
    free(rvec_old);


    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
}

