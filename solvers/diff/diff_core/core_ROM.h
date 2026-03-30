
#pragma once

#include "core_FOM.h"
#include "hlpod_dataset.h"
#include "diff_dataset.h"

#include "BBFE/manusol/manusol.h"


void ROM_read_args(
    int 		argc,
    char* 		argv[],
    ROM_PRM*    rom_prm);

const char* ROM_std_hlpod_get_parted_file_name(
    const int   solver_type);

const char* ROM_std_hlpod_get_metagraph_name(
    const int   solver_type);

void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type);

void ROM_offline_assign_default_values(
    VALUES*     vals);

void ROM_offline_print_all_values(
    VALUES*     vals);

void ROM_offline_read_calc_conditions(
    VALUES*         vals,
    const char* 	directory);

void ROM_online_assign_default_values(
    VALUES*     vals);

void ROM_online_print_all_values(
    VALUES*     vals);

void ROM_online_read_calc_conditions(
    VALUES*         vals,
    const char* 	directory);

void ROM_output_result_file_vtk(
    BBFE_DATA*     fe,
    VALUES*        vals,
    HLPOD_VALUES*   hlpod_vals,
    const char*    filename,
    const char*    directory,
    double         t);

void ROM_output_files(
    FE_SYSTEM* sys,
    int file_num,
    double t);

void solver_rom(
    FE_SYSTEM* sys,
    const int step,
    const double t);
