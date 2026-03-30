#pragma once

#include "monolis.h"

#include "BB/std.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


typedef struct
{
	double* sol_vec;
    int    incsvd_interval;

    double time_calc_reduced_matvec;
    double time_linear_solver;

} HR_VALUES;


//Hyper-reduction (ECM) の計算情報
typedef struct
{
	int num_selected_elems;
	int num_selected_elems_D_bc;

	int* id_selected_elems;
	int* id_selected_elems_D_bc;

	double* elem_weight;
	double* elem_weight_D_bc;

	double** matrix;
	double* RH;
	double** reduced_mat;
	double* reduced_RH;

	int num_selected_nodes_D_bc;
	int* D_bc_node_id;
	double* imposed_D_val;
	bool* D_bc_exists;

} HLPOD_HR;


//domain decomposition + Hyper-reduction (ECM) の計算情報
typedef struct
{
	double** reduced_mat;
	double* reduced_RH;

	int* num_selected_elems;
	int* num_selected_elems_D_bc;

	int** id_selected_elems;
	int** id_selected_elems_D_bc;

	double** elem_weight;
	double** elem_weight_D_bc;

	double*** matrix;
	double** RH;

	int* num_selected_nodes_D_bc;
	int** D_bc_node_id;
	double** imposed_D_val;
	bool** D_bc_exists;

//for DDECM
	int** elem_id_local;
	int* num_elems;

//for DDECM_par
	int* total_num_elems;
	int** ovl_elem_global_id;

	int* ovl_id_selected_elems;
	int* ovl_id_selected_elems_D_bc;

	double* ovl_elem_weight;
	double* ovl_elem_weight_D_bc;

	int ovl_num_selected_elems;
	int ovl_num_selected_elems_D_bc;

	int* parallel_elems_id;

	int num_internal_elems;

	int* num_neibs; 
	int* num_neib_modes_1stdd;
	int** neibs_1stdd;
	int* num_modes_1stdd;
	int* num_neib_modes;
	int* num_neib_modes_1stdd_sum;
	int* num_internal_modes_1stdd_sum;

} HLPOD_DDHR;

typedef struct
{
    HR_VALUES       hr_vals;
    HLPOD_HR        hlpod_hr;
    HLPOD_DDHR      hlpod_ddhr;

} HROM;
