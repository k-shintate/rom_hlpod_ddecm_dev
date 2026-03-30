#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void HROM_ddecm_memory_allocation_sol_vec(
    HR_VALUES*      hr_vals,
    const int       total_num_nodes,
    const int       dof);

void HROM_ddecm_memory_allocation_para_online(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int       total_num_nodes);

void HROM_ddecm_memory_allocation_para(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*     	hlpod_mat,
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int		num_subdomains);

void HROM_ddecm_memory_free_para(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int		num_subdomains);

/*for visualization*/
void HROM_ddecm_set_selected_elems_para_vis(
    BBFE_DATA*     	fe,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		total_num_nodes,
    const int		num_subdomains,
    const char*     directory);

void HROM_ddecm_to_monollis_rhs_para(
    MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
    HLPOD_MAT*      hlpod_mat,
    const int		num_2nddd,
    const int		max_num_bases);

void HROM_ddecm_calc_block_solution(
	MONOLIS_COM*	monolis_com,
	BBFE_DATA* 		fe,
    HR_VALUES*      hr_vals,
	HLPOD_MAT*      hlpod_mat,
	const int		num_2nddd,
    const int 		dof);

void HROM_ddecm_add_monollis_rhs_para(
	MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
	HLPOD_MAT*    hlpod_mat,
	const int		num_2nddd,
	const int		max_num_bases);