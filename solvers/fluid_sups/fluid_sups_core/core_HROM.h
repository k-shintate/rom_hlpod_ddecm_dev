#pragma once

#include "core_FOM.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

#include "BBFE/manusol/manusol.h"

#include "ecm_write.h"
#include "inc_svd.h"

#include "HR.h"
#include "DDHR.h"
#include "DDHR_para.h"
#include "DDHR_para_lb.h"
#include "set_matvec.h"
#include "set_matvec_NNLS.h"
#include "set_modes.h"
#include "monowrap.h"

#include "core_FOM.h"
#include "core_ROM.h"
#include "core_NR.h"


void HROM_offline_read_calc_conditions_inc_svd(
    HR_VALUES*      hr_vals,
    const char* 	directory);

void HROM_std_hlpod_offline_set_num_snapmat_inc_svd(
    ROM*            rom,
    const double    finish_time,
    const double    dt,
    const int       snapshot_interval,
    const int       num_case,
    const int       incsvd_interval);

void HROM_output_files(
    FE_SYSTEM*      sys,
    int             file_num,
    double          t);

void HROM_pre_offline(
    FE_SYSTEM* sys,
    ROM*            rom,
    HROM*           hrom,
    const int num_modes,
    const int num_snapshot,
    const int num_2nd_subdomains);

void HROM_pre_offline2(
    FE_SYSTEM* 		sys,
    ROM*            rom,
    HROM*           hrom,
    const int 		num_modes,
    const int 		num_snapshot,
    const int 		num_2nd_subdomains);

void HROM_pre_online(
    FE_SYSTEM* sys,
    ROM*            rom,
    HROM*           hrom,
    const int num_modes,
    const int num_snapshot,
    const int num_2nd_subdomains);

void HROM_nonparallel(
    FE_SYSTEM       sys,
    ROM*            rom,
    HROM*           hrom,
    const int       step_HR,
    const int       step_POD,
    const double    t);

void HROM_hierarchical_parallel(
    FE_SYSTEM 		sys,
    ROM*            rom,
    HROM*           hrom,
    const int 		step_HR,
    const int 		step_POD,
    const double 	t);

void HROM_hierarchical_parallel_NR(
    FE_SYSTEM   sys,
    ROM*        rom,
    HROM*       hrom,
    const int   step_HR,
    const int   step_POD,
    const double t);

void HROM_pre(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom);

void HROM_memory_allocation(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom);

void HROM_set_matvec(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom,
    int step,
    double t);

void HROM_pre_offline3(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom);

void HROM_memory_allocation_online(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom);

void HROM_pre_online(
    FE_SYSTEM* sys,
    ROM*        rom,
    HROM*       hrom);

void HROM_std_hlpod_pre_lpod_para(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    ROM*		 rom,
    const int    dof,
    const char*	 directory);

void HROM_std_hlpod_online_pre(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* mono_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    ROM* 		 rom_sups,
    const int 	 total_num_nodes,
    const int 	 dof,
    const char*  metagraph,
    const char*  directory);

void HROM_pre_offline_inc_svd1(
    FE_SYSTEM* sys,
    ROM*            rom,
    HROM*           hrom,
    const int num_modes,
    const int num_snapshot,
    const int num_2nd_subdomains);

void HROM_pre_offline_inc_svd2(
    FE_SYSTEM* sys,
    ROM*            rom,
    HROM*           hrom,
    const int num_modes,
    const int num_snapshot,
    const int num_2nd_subdomains);

void HROM_pre_offline_inc_svd3(
    FE_SYSTEM* sys,
    ROM*            rom,
    HROM*           hrom,
    const int num_modes,
    const int num_snapshot,
    const int num_2nd_subdomains);

void solver_hrom_NR(
    FE_SYSTEM *  sys,
    double      t,
    const int   step,
    const int   step_hrom);