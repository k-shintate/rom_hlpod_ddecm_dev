
#pragma once

#include "rom_dataset.h"

void ROM_std_hlpod_set_snapmat_nobc(
    double*       	comp_vec,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int		dof,
    const int 		count);

void ROM_std_hlpod_calc_reduced_mat_seq(
    MONOLIS*		monolis,
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*     hlpod_mat,
    const int 		total_num_nodes,
    const int		num_base,
    const int		dof);

void ROM_std_hlpod_calc_reduced_mat_seq_block(
    MONOLIS*		monolis,
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*     hlpod_mat,
    const int 		total_num_nodes,
    const int		num_base,
    const int		num_2nddd,
    const int 		dof);

//省メモリver
void ROM_std_hlpod_calc_reduced_mat_save_memory(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    MONOLIS_COM*    mono_com0,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int       n_neib_vec,
    const int       num_2nddd,
    const int       num_modes,
    const int 		dof);

void ROM_std_hlpod_set_reduced_mat(
    MONOLIS*		monolis,
    HLPOD_MAT*   hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		total_num_nodes,
    const int		num_base,
    const int		num_2nddd,
    const int 		dof);

void ROM_std_hlpod_set_reduced_mat_para(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		num_modes,
    const int		total_num_bases,
    const int		num_2nddd);

void ROM_std_hlpod_reduced_rhs_to_monollis(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int       num_2nd_subdomains,
    const int		num_modes);


void ROM_std_hlpod_calc_reduced_rhs(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		max_num_bases,
    const int		num_2nddd,
    const int 		dof);

void ROM_std_hlpod_calc_sol(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int		num_2nddd,
    const int		max_num_bases,
    const int 		dof);

void ROM_std_hlpod_calc_sol_global_para(
    double*         ansvec,
    double**        pod_modes,
    double*         mode_coef,
    const int 		total_num_nodes,
    const int 		num_base,
    const int		dof);

void ROM_std_hlpod_update_global_modes(
    MONOLIS_COM*	monolis_com,
    HLPOD_MAT*		hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int		ndof);

void ROM_std_hlpod_calloc_mode_coef_rhs(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		max_num_bases,
    const int		total_num_bases,
    const int		num_2nddd);

void ROM_std_hlpod_calc_reduced_rhs_add(
    MONOLIS*		monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		max_num_bases,
    const int		num_2nddd,
    const int 		dof);