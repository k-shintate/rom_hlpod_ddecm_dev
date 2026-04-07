
#pragma once

#include "monolis.h"
#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

#include "rom_dataset.h"
#include "std.h"

#include "write_BB.h"
#include "write_std.h"
#include "read_BB.h"
#include "read_std.h"

void ROM_std_hlpod_set_podmodes(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		num_snapshots,
    const double    rom_epsilon,
    const int 	  	dof,
    const char*  	label,
    const char*  	directory);

void ROM_std_hlpod_free_local_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_2nd_subdomains,
    const int 		num_modes,
    const int 		dof);

void ROM_std_hlpod_free_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		dof);

void ROM_std_hlpod_free_local_podmodes_para(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_2nd_subdomains,
    const int 		num_modes,
    const int 		dof);

void ROM_std_hlpod_free_global_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes,
    const int 		dof);

void ROM_std_hlpod_read_podmodes(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		num_snapshots,
    const double    rom_epsilon,
    const int 	  	dof,
    const char*  	label,
    const char*  	directory);

void ROM_std_hlpod_set_podmodes_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*	    hlpod_mat,
    double** 		v,
    double** 		p,
    const int 		total_num_nodes,
    const int 		num_modes_max_1,
    const int 		num_modes_max_2,
    const int		dof_1,
    const int 		dof_2,
    const char*  	directory);

void ROM_std_hlpod_read_podmodes_local(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int       num_modes,
    const int       num_snapshots,
    const int       num_2nd_subdomains,
    const int 		dof,
    const char*  	label,
    const char*     directory);

void ROM_std_hlpod_read_podmodes_local_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int		num_2nd_subdomains,
    const int 		dof,
    const double    rom_epsilon,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int		num_2nd_subdomains,
    const double    rom_epsilon,
    const int 		dof,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int       num_modes,
    const int       num_snapshots,
    const int       num_2nd_subdomains,
    const double    rom_epsilon,
    const int		dof,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    int* 			n_internal_vertex_subd,
    int* 			node_id_local,
    const int 		num_2nd_subdomains,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int		num_2nd_subdomains,
    const double    rom_epsilon,
    const int 		dof,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para_diag_v(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para_diag_p(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_local_para_Aphi(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_global_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    double**		snapmat,
    const int		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const double    rom_epsilon,
    const int 		dof,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_read_podmodes_global_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    const int       num_modes,
    const int		dof,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_global_para_Aphi(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory);

void ROM_std_hlpod_set_podmodes_global_para_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory);

void HROM_ddecm_set_podbasis_ovl(
    MONOLIS_COM*  	monolis_com,
	HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
	const int		total_num_nodes,
    const int       dof);
