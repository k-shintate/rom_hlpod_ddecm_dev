#pragma once

#include "hlpod_dataset.h"
#include "diff_dataset.h"

double manusol_get_sol(
		double x,
		double y,
		double z,
		double t);

void manusol_get_conv_vel(
		double v[3],
		double x[3]);

double manusol_get_mass_coef(
		double x[3]);


double manusol_get_diff_coef(
		double x[3]);

double manusol_get_source(
		double x[3],
		double t,
		double a,
		double v[3],
		double k);

void manusol_set_theo_sol(
		BBFE_DATA* fe,
		double*  theo_sol,
		double   t);

void manusol_set_source(
		BBFE_DATA* fe,
		double*  source,
		double   t);

void manusol_set_init_value(
		BBFE_DATA* fe,
		double* T);

void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes);

void assign_default_values(
		VALUES*     vals);

void print_all_values(
		VALUES*  vals);

void read_calc_conditions(
		VALUES*     vals,
		const char* directory);

void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t);

void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t);

void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals);

void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals,
		double       t);

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
		FE_SYSTEM sys,
		double t,
		const int step);

//スナップショットの収集
void solver_fom_collect_snapmat(
		FE_SYSTEM sys,
		double t,
		const int step);