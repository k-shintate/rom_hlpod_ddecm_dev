
#include "core_FOM.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;
const char*         ID_DENSITY  = "#density";
const double      DVAL_DENSITY  = 1000.0;
const char*       ID_VISCOSITY  = "#viscosity";
const double    DVAL_VISCOSITY  = 1.0;

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%d_%06d.vtk";


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
	vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;

	vals->density          = DVAL_DENSITY;
	vals->viscosity        = DVAL_VISCOSITY;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s %s: %e\n", CODENAME, ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ID_VISCOSITY,        vals->viscosity);
	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->density), filename, ID_DENSITY, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->viscosity), filename, ID_VISCOSITY, BUFFER_SIZE, CODENAME);


		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_vector(fp, vals->v, fe->total_num_nodes, "Velocity");
	BB_vtk_write_point_vals_scalar(fp, vals->p, fe->total_num_nodes, "Pressure");

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	int myrank = monolis_mpi_get_global_my_rank();
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, myrank, file_num);

	output_result_file_vtk(
			&(sys->fe), &(sys->vals), fname_vtk, sys->cond.directory, t);
}


void set_element_mat_pred(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					double tau = BBFE_elemmat_fluid_supg_coef(
							vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

					val_ip[p] = BBFE_elemmat_fluid_fs_mat_pred_expl(
							basis->N[p][i], basis->N[p][j], fe->geo[e][p].grad_N[i], v_ip[p], tau);
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis, fe->conn[e][i], fe->conn[e][j], 0, 0, integ_val);
				monolis_add_scalar_to_sparse_matrix_R(
						monolis, fe->conn[e][i], fe->conn[e][j], 1, 1, integ_val);
				monolis_add_scalar_to_sparse_matrix_R(
						monolis, fe->conn[e][i], fe->conn[e][j], 2, 2, integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);

	BB_std_free_2d_double(v_ip, np, 3);
}


void set_element_vec_pred(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 3, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double*** grad_v_ip;
	grad_v_ip = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			double integ_val[3];
	
			for(int p=0; p<np; p++) {
				double tau = BBFE_elemmat_fluid_supg_coef(
						vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

				double vec[3];
				BBFE_elemmat_fluid_fs_vec_pred_expl(
						vec, basis->N[p][i], fe->geo[e][p].grad_N[i],
						v_ip[p], grad_v_ip[p],
						vals->density, vals->viscosity, tau, vals->dt);

				for(int d=0; d<3; d++) {
					val_ip[d][p] = vec[d];
				}
			}

			for(int d=0; d<3; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

				monolis->mat.R.B[ 3*fe->conn[e][i] + d ] += integ_val[d];
			}
		}
	}
	
	BB_std_free_2d_double(val_ip, 3, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_3d_double(grad_v_ip, np, 3, 3);
}


void set_element_vec_corr(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 3, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;  double* local_p;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	local_p = BB_std_calloc_1d_double(local_p, nl);

	double** v_ip;  double** grad_p_ip;
	v_ip      = BB_std_calloc_2d_double(v_ip     , np, 3);
	grad_p_ip = BB_std_calloc_2d_double(grad_p_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				double vec[3];
				
				BBFE_elemmat_fluid_fs_vec_corr(
						vec, basis->N[p][i], grad_p_ip[p],
						v_ip[p], vals->density, vals->dt);

				for(int d=0; d<3; d++) {
					val_ip[d][p] = vec[d];
				}
			}

			double integ_val[3];
			
			for(int d=0; d<3; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

				monolis->mat.R.B[ 3*(fe->conn[e][i]) + d ] += integ_val[d];
			}
		}
	}

	BB_std_free_2d_double(val_ip, 3, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_1d_double(local_p, nl);
	BB_std_free_2d_double(v_ip     , np, 3);
	BB_std_free_2d_double(grad_p_ip, np, 3);
}


void set_element_vec_ppe(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	double* div_v_ip;
	div_v_ip = BB_std_calloc_1d_double(div_v_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		for(int p=0; p<np; p++) {
			div_v_ip[p] = BBFE_std_mapping_vector3d_div(nl, local_v, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = BBFE_elemmat_fluid_fs_vec_ppe(
						basis->N[p][i], div_v_ip[p], vals->density, vals->dt);
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_v, fe->local_num_nodes, 3);
	BB_std_free_1d_double(div_v_ip, np);
}


//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

    BBFE_sys_monowrap_copy_mat(&(sys.mono_ppe0) , &(sys.mono_ppe));
    BBFE_sys_monowrap_copy_mat(&(sys.mono_corr0), &(sys.mono_corr));

    monolis_clear_mat_value_R(&(sys.mono_pred));
    monolis_clear_mat_value_rhs_R(&(sys.mono_ppe));
    monolis_clear_mat_value_rhs_R(&(sys.mono_corr));

    printf("%s --- prediction step ---\n", CODENAME);
    set_element_mat_pred(
            &(sys.mono_pred),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    set_element_vec_pred(
            &(sys.mono_pred),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_pred),
            sys.fe.total_num_nodes,
            3,
            &(sys.bc_v),
            sys.mono_pred.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_pred),
            &(sys.mono_com),
            sys.mono_pred.mat.R.X,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    BBFE_fluid_renew_velocity(
            sys.vals.v, 
            sys.mono_pred.mat.R.X,
            sys.fe.total_num_nodes);

    printf("%s --- pressure Poisson eq. ---\n", CODENAME);
    set_element_vec_ppe(
            &(sys.mono_ppe),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_ppe),
            sys.fe.total_num_nodes,
            1,
            &(sys.bc_p),
            sys.mono_ppe.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_ppe),
            &(sys.mono_com),
            sys.vals.p,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    printf("%s --- Correction step ---\n", CODENAME);
    set_element_vec_corr(
            &(sys.mono_corr),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_corr),
            sys.fe.total_num_nodes,
            3,
            &(sys.bc_v),
            sys.mono_corr.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_corr),
            &(sys.mono_com),
            sys.mono_corr.mat.R.X,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    BBFE_fluid_renew_velocity(
            sys.vals.v, 
            sys.mono_corr.mat.R.X,
            sys.fe.total_num_nodes);

}

//スナップショットの収集
void solver_FOM_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    const int step)
{
    printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

    BBFE_sys_monowrap_copy_mat(&(sys.mono_ppe0) , &(sys.mono_ppe));
    BBFE_sys_monowrap_copy_mat(&(sys.mono_corr0), &(sys.mono_corr));

    monolis_clear_mat_value_R(&(sys.mono_pred));
    monolis_clear_mat_value_rhs_R(&(sys.mono_ppe));
    monolis_clear_mat_value_rhs_R(&(sys.mono_corr));

    printf("%s --- prediction step ---\n", CODENAME);
    set_element_mat_pred(
            &(sys.mono_pred),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    set_element_vec_pred(
            &(sys.mono_pred),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_pred),
            sys.fe.total_num_nodes,
            3,
            &(sys.bc_v),
            sys.mono_pred.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_pred),
            &(sys.mono_com),
            sys.mono_pred.mat.R.X,
            MONOLIS_ITER_BICGSTAB,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    BBFE_fluid_renew_velocity(
            sys.vals.v, 
            sys.mono_pred.mat.R.X,
            sys.fe.total_num_nodes);

    printf("%s --- pressure Poisson eq. ---\n", CODENAME);
    set_element_vec_ppe(
            &(sys.mono_ppe),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_ppe),
            sys.fe.total_num_nodes,
            1,
            &(sys.bc_p),
            sys.mono_ppe.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_ppe),
            &(sys.mono_com),
            sys.vals.p,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);
    
    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_sys_hlpod_fe_set_snap_mat_para(
                    sys.vals.p,
                    &(sys.rom_p.hlpod_mat),
                    &(sys.bc_p),
                    &(sys.rom_p.rom_bc),
                    sys.fe.total_num_nodes,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
        }
        else{
            ROM_sys_hlpod_fe_set_snap_mat_para(
                    sys.vals.p,
                    &(sys.rom_p.hlpod_mat),
                    &(sys.bc_p),
                    &(sys.rom_p.rom_bc),
                    sys.mono_com.n_internal_vertex,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
        }

    }

    printf("%s --- Correction step ---\n", CODENAME);
    set_element_vec_corr(
            &(sys.mono_corr),
            &(sys.fe),
            &(sys.basis),
            &(sys.vals));
    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys.mono_corr),
            sys.fe.total_num_nodes,
            3,
            &(sys.bc_v),
            sys.mono_corr.mat.R.B);
    BBFE_sys_monowrap_solve(
            &(sys.mono_corr),
            &(sys.mono_com),
            sys.mono_corr.mat.R.X,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_SOR,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    BBFE_fluid_renew_velocity(
            sys.vals.v, 
            sys.mono_corr.mat.R.X,
            sys.fe.total_num_nodes);

    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_sys_hlpod_fe_set_snap_mat_para(
                    sys.mono_corr.mat.R.X,
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc_v),
                    &(sys.rom_v.rom_bc),	//要変更
                    sys.fe.total_num_nodes,
                    3,
                    ((int)step/sys.vals.snapshot_interval));
        }
        else{
            ROM_sys_hlpod_fe_set_snap_mat_para(
                    sys.mono_corr.mat.R.X,
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc_v),
                    &(sys.rom_v.rom_bc),
                    sys.mono_com.n_internal_vertex,
                    3,
                    ((int)step/sys.vals.snapshot_interval));
        }

    }

}
