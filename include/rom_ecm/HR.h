#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void HROM_ecm_set_element(
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    HLPOD_HR*       hlpod_hr);

void HROM_ecm_memory_allocation_online(
    const int       total_num_nodes,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    HLPOD_HR*       hlpod_hr);

void HROM_ecm_get_selected_elems(
    BBFE_DATA*     	fe,
    BBFE_BC*     	bc,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int       max_iter, //NNLS
    const double    tol,      //NNLS
    HLPOD_HR*       hlpod_hr,
    const char*		directory);

void HROM_ecm_calc_solution(
	BBFE_DATA* 		fe,
	HLPOD_MAT*     hlpod_mat,
	double*       HR_T,
	BBFE_BC*     	bc,
    int 			num_base);

void HROM_ecm_monolis_set_matrix(
    MONOLIS*     	monolis,
    HLPOD_HR*      hlpod_hr,
    const int 		num_basis);

void HROM_ecm_to_monollis_rhs(
    MONOLIS*		monolis,
    HLPOD_HR*       hlpod_rh,
    const int 		k);

void HROM_ecm_set_selected_elems(
    BBFE_DATA*     	fe,
    HLPOD_HR*       hlpod_hr,
    const int		total_num_nodes,
    const char*     directory);

void HROM_ecm_read_selected_elems(
    HLPOD_HR*     hlpod_hr,
    const int num_subdomains,
    const char* directory);