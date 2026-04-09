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

#include "DDHR_para_lb_decoupled.h"

void HROM_pre_offline_decoupled(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains);

void HROM_pre_offline2_decoupled(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains);
    
void HROM_pre_online_decoupled(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
        const char* name1,
        const char* name2,
		const int num_2nd_subdomains);

void solver_hrom_NR_decoupled(
        FE_SYSTEM *  sys,
        double      t,
        const int   step,
        const int   step_hrom);

void HROM_pre_decoupled(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_pre_offline3_decoupled(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_pre_online_decoupled2(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_std_hlpod_pre_lpod_para_decoupled(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM*		 rom,
        const int    dof,
        const char*	 directory);

void HROM_std_hlpod_online_pre_decoupled(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* mono_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM* 		 rom_sups,
        const int 	 total_num_nodes,
        const int 	 dof,
        const char*  metagraph,
        const char*  directory);
