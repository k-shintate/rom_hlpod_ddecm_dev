#pragma once

#include "fluid_core.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

void output_cavity_center_vx(
    VALUES*        vals,
    const char*    method,
    const char*    directory);

void initialize_velocity_pressure_karman_vortex(
	double** v,
	double* p,
	const int total_num_nodes);

void initialize_velocity_pressure_cavity(
	double** v,
    double** v_old,
	double* p,
	const int total_num_nodes);

void BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void output_result_file_karman_vortex(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure_inf(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void memory_allocation_nodal_values(
    VALUES*         vals,
    const int       total_num_nodes);

void assign_default_values(
    VALUES*         vals);

void print_all_values(
    VALUES*         vals);

void read_calc_conditions(
    VALUES*         vals,
    const char*     directory);

void output_result_file_vtk(
    BBFE_DATA*     fe,
    VALUES*        vals,
    const char*    filename,
    const char*    directory,
    double         t);

void output_files(
    FE_SYSTEM*  sys,
    int         file_num,
    double      t);

void set_element_mat_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_corr(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_ppe(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void BBFE_fluid_sups_read_Dirichlet_bc(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
    BBFE_BC*     bc,
    VALUES*      vals,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void set_element_mat(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void solver_fom_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void memory_allocation_nodal_values_AB2(
    VALUES*         vals,
    const int       total_num_nodes);

void ROM_BB_vec_copy_2d(
    double**  in,   //input
    double**  out,  //output
    const int num1,
    const int num2);

void BBFE_fluid_sups_add_velocity_pressure(
    double**  v,
    double*   p,
    double*   sol_vec,
    const int total_num_nodes); 

void solver_fom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void pre_surface(
    BBFE_DATA*    surf,
    BBFE_BASIS*   basis,
    const char*   directory,
    const char*   fname,
    int           num_integ_points_each_axis);

void calc_Cd_p(
    BBFE_DATA*  surf,
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    VALUES*     vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out);

int calc_Cd_v(const BBFE_DATA* surf, const BBFE_DATA* fe, const BBFE_BASIS* basis, MONOLIS_COM* monolis_com,
                       const VALUES* vals,
                       double rho, double mu, double Uinf, double Aref,
                       const double eU_in[3], const double eP_in[3],
                       double* D_out, double* L_out, double* Cd_out, double* Cl_out);

void update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes);

double calc_internal_norm_2d(
    double**  in,       //input
        const int num1,
    const int num2);

    double calc_internal_norm_1d(
        double*  in,       //input
        const int num1,
        const int num2);

void set_element_mat_NR_Tezuer(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

double BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

void solver_fom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
                BBFE_BC*     bc,
                const char*  filename,
                const char*  directory,
                const int    total_num_nodes,
                const int    block_size);

void BBFE_fluid_sups_add_velocity_pressure(
                double**  v,
                double*   p,
                double*   sol_vec,
                const int total_num_nodes);

void solver_fom_NR_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

void ROM_ddecm_set_reduced_mat_para_save_memory(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);


void ROM_ddecm_set_reduced_vec_para(
    MONOLIS*        monolis,
    BBFE_DATA*      fe,
    VALUES*         vals,
    BBFE_BASIS*     basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int       num_modes,
    const int       num_subdomains,
    const double    dt,
    double          t);

void ROM_set_reduced_vec_para_debug(
        MONOLIS*        monolis,
        BBFE_DATA*      fe,
        VALUES*         vals,
        BBFE_BASIS*     basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
        HLPOD_MAT*      hlpod_mat,
        const int       num_modes,
        const int       num_subdomains,
        const double    dt,
        double          t);

void ROM_set_reduced_mat_para_save_memory_debug(
        MONOLIS*        monolis,
        BBFE_DATA*      fe,
        VALUES*         vals,
        BBFE_BASIS*     basis,
        BBFE_BC*        bc,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int       num_modes,
        const int       num_subdomains,
        const double    dt);

void HROM_ddecm_set_residuals_NR_Tezuer(
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals,
    BBFE_BC*        bc,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int               num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_reduced_mat_para_save_memory_NR(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);

void HROM_ddecm_set_reduced_vec_para_NR(
        MONOLIS*         monolis,
        BBFE_DATA*       fe,
        VALUES*          vals,
        BBFE_BASIS*      basis,
        HR_VALUES*       hr_vals,
        HLPOD_VALUES*    hlpod_vals,
        HLPOD_DDHR*      hlpod_ddhr,
        HLPOD_MAT*       hlpod_mat,
        const int        num_modes,
        const int        num_subdomains,
        const double     dt,
        double           t);


void HROM_ddecm_set_reduced_mat_para_save_memory_NR2(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);

void HROM_ddecm_set_reduced_vec_para_NR2(
        MONOLIS*         monolis,
        BBFE_DATA*       fe,
        VALUES*          vals,
        BBFE_BASIS*      basis,
        HR_VALUES*       hr_vals,
        HLPOD_VALUES*    hlpod_vals,
        HLPOD_DDHR*      hlpod_ddhr,
        HLPOD_MAT*       hlpod_mat,
        const int        num_modes,
        const int        num_subdomains,
        const double     dt,
        double           t);

