#pragma once

#include "core_FOM.h"

//#include "set_matvec.h"
#include "convdiff_core.h"
#include "BBFE/manusol/manusol.h"

#include "nedelec_core.h"
#include "3ph_tr_NR.h"
#include "set_matvec.h"

void initialize_velocity_pressure(
	double**        v,
	double*         p,
	const int       total_num_nodes);

void ROM_std_hlpod_online_memory_allocation_ansvec(
    HLPOD_VALUES*	hlpod_vals,
    const int		total_num_nodes,
    const int	    dof);

void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type);

void ROM_offline_set_reynolds_num_cases(
    VALUES*         vals,
    const char*     directory);

void ROM_offline_set_reynolds_number( 
    VALUES*         vals,
    const int       case_id);

const char* ROM_std_hlpod_get_parted_file_name(
    const int       solver_type);

const char* ROM_std_hlpod_get_metagraph_name(
    const int       solver_type);

void set_target_parameter(
    VALUES*         vals,
    const char*     directory);

void ROM_set_ansvec(
    VALUES*         vals,
    HLPOD_VALUES*   hlpod_vals,
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

void output_B_node_ROM(
    FE_SYSTEM*    sys,
    const double* Aphi,
    const double* Aphi_rom, const double t,
    const char* filename, const char* directory);

void solver_rom(
    FE_SYSTEM*      sys,
    const int       step_HR,
    const int       step_rom,
    const double    t);

void solver_rom_global_para(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    ROM*		    rom,
    FE_SYSTEM*      sys,
    const int       step_HR,
    const int       step_rom,
    const double    t);

void HROM_std_hlpod_online_pre(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* mono_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    ROM* 		 rom_sups,
    const int 	 total_num_nodes,
    const int    dof,
    const char*  metagraph,
    const char*  directory);

void solver_rom_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total);

void solver_rom_global_para_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total);

void solver_rom_NR_Aphi_team21a2(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total);
