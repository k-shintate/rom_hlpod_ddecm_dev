#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void HROM_ddecm_get_selected_elems_int_ovl_decoupled(
	HLPOD_DDHR*     hlpod_ddhr,
	const char*     directory);

void HROM_ddecm_read_selected_elems_para_decoupled(
	const int num_subdomains,
	const char* directory);

void HROM_ddecm_write_selected_elems_para_arbit_subd_decoupled(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory);

void HROM_ddecm_set_neib_decoupled(
		MONOLIS_COM*  	monolis_com,
		HLPOD_MAT* 	hlpod_mat,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_subdomains,
		const int       num_snapshots,
		const char*     directory);

void HROM_ddecm_get_selected_elema_add_decoupled(
	HLPOD_DDHR*     hlpod_ddhr,
	const int       num_parallel_subdomains,
	const char*     directory);

void ddhr_lb_set_selected_elems_para_decoupled(
		BBFE_DATA*     	fe,
		HLPOD_DDHR*     hlpod_ddhr,
		const int		total_num_nodes,
		const int		num_subdomains,
		const char*     directory);

void HROM_ddecm_write_selected_elems_para_arbit_subd_svd_decoupled(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory);

