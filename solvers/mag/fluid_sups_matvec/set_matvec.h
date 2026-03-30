#pragma once

#include "BBFE/sys/FE_dataset.h"

#include "fluid_sups_dataset.h"

#include "set_reduced_matvec.h"


void ddhr_lb_set_reduced_mat_para_save_memory_debug(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
	VALUES*         vals,
    BBFE_BASIS* 	basis,
    BBFE_BC*     	bc,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_modes,
    const int 		num_subdomains,
    const double    dt);

void ddhr_lb_set_D_bc_para_debug(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS* 	basis,
    BBFE_BC*     	bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt);

void ddhr_lb_set_reduced_vec_para_debug(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS*	 	basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt,
    double       	t);

void HROM_ddecm_set_reduced_mat_para_save_memory(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
	VALUES*         vals,
    BBFE_BASIS* 	basis,
    BBFE_BC*     	bc,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_modes,
    const int 		num_subdomains,
    const double    dt);

void HROM_ddecm_set_D_bc_para(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS* 	basis,
    BBFE_BC*     	bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt);

void HROM_ddecm_set_reduced_vec_para(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS*	 	basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt,
    double       	t);

/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
    MONOLIS*        monolis,
    MONOLIS_COM*	monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    double**        mat,
    double**        pod_modes,
    const int 		num_modes);

void set_D_bc_global_para(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*     	bc,
    double*         rhs,
    double**        pod_modes,
    const int 		num_modes);

void set_reduced_vec_global_para(
    MONOLIS*        monolis,        
    MONOLIS_COM*	monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    double*         rhs,
    double**        pod_modes,
    const int		num_modes);

void allreduce_global_para(
    MONOLIS_COM*	monolis_com,
    double**        mat,
    double*         rhs,
    const int		num_modes);
