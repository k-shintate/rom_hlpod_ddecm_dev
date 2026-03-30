#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void HROM_ddecm_get_selected_elems_int_ovl(
    HLPOD_DDHR*     hlpod_ddhr,
    const char*		directory);

void HROM_ddecm_read_selected_elems_para(
    const int 		num_subdomains,
    const char*		directory);

void HROM_ddecm_write_selected_elems_para(
    MONOLIS_COM*  	monolis_com,
    BBFE_DATA*     	fe,
    BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int 		num_subdomains,
    const int       max_iter, //NNLS
    const double    tol,      //NNLS
    const int       dof,
    const char*		directory);

void HROM_ddecm_write_selected_elems_svd(
	MONOLIS_COM*  	monolis_com,
	BBFE_DATA*     	fe,
	BBFE_BC*     	bc,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*      hlpod_mat,
	HLPOD_META*		hlpod_meta,
	const int       total_num_elem,
	const int       total_num_snapshot,
	const int       total_num_modes,
	const int 		num_subdomains,
	const int       max_iter, //NNLS
	const double    tol,      //NNLS
    const int       dof,
	const char*		directory);

void HROM_ddecm_write_selected_elems_para_arbit_subd(
    MONOLIS_COM*        monolis_com,
    BBFE_DATA*          fe,
    BBFE_BC*            bc,
    HLPOD_VALUES*       hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*         hlpod_meta,
    const int       total_num_elem,
    const int       total_num_snapshot,
    const int       total_num_modes,
    const int           num_subdomains,
    const int       max_iter, //NNLS
    const double    tol,      //NNLS
    const int       dof,
    const char*         directory);

void HROM_ddecm_write_selected_elems_para_arbit_subd_svd(
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


void HROM_ecm_get_meta_neib(
    MONOLIS_COM*  	monolis_com,
    HLPOD_META*		hlpod_meta,
    const char*     directory);

void HROM_ddecm_set_neib(
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_META*		hlpod_meta,
    const int 		num_subdomains,
    const int       num_snapshots,
    const char*     directory);

void HROM_ddecm_set_element_para(
    BBFE_DATA*     	fe,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_subdomains,
    const char*     directory);

void HROM_ddecm_get_selected_elema_add(
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_parallel_subdomains,
    const char*		directory);

void HROM_ddecm_calc_block_mat_bcsr(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_META*		hlpod_meta,
    const int 		num_bases,
    const int		num_2nddd,
    const char*		directory);

void ddhr_lb_set_selected_elems_para(
    BBFE_DATA*     	fe,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		total_num_nodes,
    const int		num_subdomains,
    const char*     directory);

//level1領域の最大基底本数の共有
void HROM_ddecm_get_neib_max_num_modes(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    const int       np,
    const int       num_modes_pre);

//level1領域の選択された基底(p-adaptive)本数の共有
void HROM_ddecm_get_neib_num_modes(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    const int       np,
    const int       num_my_modes);

void HROM_ddecm_get_neib_subdomain_id(
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    const int 		num_modes);

void HROM_ddecm_get_neib_coordinates_pre(
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int       np,				//並列計算領域数
	const int       max_num_basis);

void HROM_get_neib_coordinates(
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int      np,				//並列計算領域数
	const int       max_num_basis,	//level 1の基底本数 (並列計算領域が担当する基底本数の総和)
	const int 		num_subdomains,
	const int		max_num_bases);

void HROM_set_bc_id(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    HLPOD_DDHR* hlpod_ddhr,
    const int   num_dofs_on_node,
    HLPOD_MAT*	hlpod_mat);

void HROM_manusol_set_bc(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    const int   num_dofs_on_node,
    double     	t,
    double      (*func)(double, double, double, double), // scalar function(x, y, z, t)
    HLPOD_MAT*	hlpod_mat);

void HROM_ddecm_get_neib_subdomain_id_2nddd(
	MONOLIS_COM*  	monolis_com,
	HLPOD_MAT* 	    hlpod_mat,
	const int 		num_modes);
