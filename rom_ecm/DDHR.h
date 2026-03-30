#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void HROM_ddecm_memory_allocation(
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int		num_subdomains,
    HLPOD_DDHR*     hlpod_ddhr);

void HROM_ddecm_calc_solution(
	BBFE_DATA* 		fe,
    HR_VALUES*      hr_vals,
	HLPOD_MAT*     hlpod_mat,
	HLPOD_DDHR*     hlpod_ddhr,
	BBFE_BC*     	bc,
    int 			num_base,
	const int		num_subdomains,
	const int		dof);

void HROM_ddecm_monolis_set_matrix(
	MONOLIS*     	monolis,
	HLPOD_MAT*     hlpod_mat,
	HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_META*     hlpod_meta,
    const int 		num_base,
	const int		num_2nddd);

void HROM_ddecm_to_monollis_rhs(
	MONOLIS*		monolis,
	HLPOD_MAT*     hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
	const int 		num_base,
	const int		num_subdomains);

void HROM_ddecm_memory_allocation2(
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int		num_subdomains,
    HLPOD_DDHR*     hlpod_ddhr);

void HROM_ddecm_set_element(
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_subdomains,
    const char*     directory);

void HROM_ddecm_get_neib_subdomain_id_nonpara(
    HLPOD_MAT* 	hlpod_mat,
    HLPOD_DDHR* 	hlpod_ddhr,
    const int       num_subdomains);

void HROM_ddecm_write_selected_elems(
    BBFE_DATA*     	fe,
    BBFE_BC*     	bc,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int 		num_subdomains,
    const int       max_iter, //NNLS
    const double    tol,      //NNLS
    HLPOD_DDHR*     hlpod_ddhr,
    const char*		directory);

void HROM_ddecm_set_selected_elems_vis(
    BBFE_DATA*     	fe,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		total_num_nodes,
    const int		num_subdomains,
    const char*     directory);

void HROM_ddecm_read_selected_elems(
    HLPOD_DDHR*     hlpod_ddhr,
    const int num_subdomains,
    const char* directory);
