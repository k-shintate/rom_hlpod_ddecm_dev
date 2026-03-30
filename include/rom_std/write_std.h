#pragma once

#include "BB/vtk.h"
#include "BB/std.h"

#include "monolis.h"

#include "rom_dataset.h"

#include <stdlib.h>

void ROM_std_hlpod_write_num_modes(
    const int       num_modes,
    const int		subdomain_id,
    const char*		label,
    const char*		directory);

int ROM_std_hlpod_read_num_modes(
    const int		subdomain_id,
    const char*		label,
    const char*		directory);

void ROM_std_hlpod_write_pod_modes(
    double**		basis,
    const int 		num_modes,
    const int 		n_internal_vertex,
    const int		subdomain_id,
    const int 		dof,
    const char*		label,
    const char*		directory);

void ROM_std_hlpod_write_singular_values(
    double*			singular_value,
    const int 		num_modes,
    const int 		subdomain_id,
    const char*		label,
    const char*		directory);

void ROM_std_hlpod_write_time_svd(
    double			time_svd,
    const int 		subdomain_id,
    const char*		label,
    const char*		directory);

void hr_write_NNLS_residual(
    const double 	residual,
    const int		myrank,
    const int		num,
    const char*		directory);

void hr_write_NNLS_num_elems(
    const int 		num_elems,
    const int		myrank,
    const int		num,
    const char*		directory);

void ddhr_lb_write_reduced_mat(
    double**		reduced_mat,
    const int		m,
    const int 		n,
    const char*		directory);

void ROM_std_hlpod_output_calc_time(
    double			calctime,
    double 			t,
    const char*  	fname,
    const char*  	directory);

void ROM_std_hlpod_output_add_calc_time(
    double			calctime,
    double 			t,
    const char*  	fname,
    const char*  	directory);

void ROM_std_hlpod_write_solver_prm(
    MONOLIS*  		monolis,
    double 			t,
    const char*     fname,
    const char*  	directory);

void ROM_std_hlpod_write_solver_prm_fopen(
    const char*     output_directory,
    const char*  	directory);

void ROM_std_hlpod_output_fem_solver_prm_fopen(
    const char*  	directory);

void ROM_std_hlpod_fem_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t);

void ROM_std_hlpod_output_rom_solver_prm_fopen(
    const char*  	directory);

void ROM_std_hlpod_rom_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t);

void ROM_std_hlpod_hrom_solver_prm_fopen(
    const char*  	directory);

void ROM_std_hlpod_hrom_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t);
