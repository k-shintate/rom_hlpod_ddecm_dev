#pragma once

#include "core_FOM.h"
#include "hlpod_dataset.h"
#include "fluid_fs_dataset.h"

#include "BBFE/manusol/manusol.h"


void initialize_velocity_pressure(
	double** v,
	double* p,
	const int total_num_nodes);

void initialize_velocity_pressure_karman_vortex(
	double** v,
	double* p,
	const int total_num_nodes);

void ROM_std_hlpod_online_memory_allocation_ansvec(
    HLPOD_VALUES*		hlpod_vals,
    const int		total_num_nodes,
    const int	    dof);

void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type);

void ROM_offline_set_reynolds_num_cases(
    VALUES*     vals,
    const char* directory);

void ROM_offline_set_reynolds_number( 
    VALUES*     vals,
    const int   case_id);

const char* ROM_std_hlpod_get_parted_file_name(
    const int   solver_type);

const char* ROM_std_hlpod_get_metagraph_name(
    const int   solver_type);

void set_target_parameter(
    VALUES*         vals,
    const char*     directory);

void ROM_set_ansvec(
    VALUES*        vals,
    HLPOD_VALUES*     hlpod_vals,
    const int       total_num_nodes);

void ROM_offline_read_calc_conditions(
    VALUES*         vals,
    const char* 	directory);

void ROM_online_read_calc_conditions(
    VALUES*         vals,
    const char* 	directory);

void ROM_output_files(
    FE_SYSTEM*      sys,
    int             file_num,
    double          t);

void solver_rom(
    FE_SYSTEM* sys,
    const int step,
    const double t);

void solver_rom_global_para(
    MONOLIS*     monolis,
    MONOLIS_COM* monolis_com,
    ROM*		 rom,
    FE_SYSTEM* sys,
    const int step_HR,
    const int step_rom,
    const double t);

