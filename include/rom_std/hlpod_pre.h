
#pragma once

#include "rom_dataset.h"
#include "std.h"
#include "read_BB.h"
#include "read_std.h"

void ROM_std_hlpod_get_meta_neib(
    MONOLIS_COM*  	monolis_com,
    HLPOD_META*		hlpod_meta,
    const char*     metagraph,
    const char*     directory);

void ROM_std_hlpod_get_subdomain_id(
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_META*		hlpod_meta,
    const int       num,
    const char*     metagraph,
    const char*     directory);

void ROM_std_hlpod_read_node_id_pod_subd(
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int		num_2nd_subdomains,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_read_node_id_pod_subd_para(
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int 		N_internal_vertex,
    const int		num_2nd_subdomains,
    const char*     label,
    const char*     directory);

void ROM_std_hlpod_set_nonzero_pattern(
    MONOLIS*     	monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		num_base);

void ROM_std_hlpod_set_nonzero_pattern_bcsr(
    MONOLIS*     	monolis,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		num_modes,
    const char*     label,
    const char*		directory);

void ROM_std_hlpod_set_nonzero_pattern_bcsr_para(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const char*     metagraph,
    const char*		directory);

void ROM_std_hlpod_get_n_dof_list(
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		max_num_bases);
