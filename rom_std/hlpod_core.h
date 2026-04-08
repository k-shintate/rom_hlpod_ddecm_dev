#pragma once

#include "monolis.h"

#include "BB/std.h"

#include "rom_dataset.h"

#include "hlpod_pre.h"
#include "hlpod_comm.h"
#include "hlpod_matvec.h"

#include "set_modes.h"
#include "monowrap.h"

void ROM_std_hlpod_online_memory_allocation_ansvec(
    HLPOD_VALUES*	hlpod_vals,
    const int		total_num_nodes,
    const int	    dof);

void ROM_std_hlpod_offline_set_num_snapmat(
    ROM*            rom,
    const double    finish_time,
    const double    dt,
    const int       snapshot_interval,
    const int       num_case);

void ROM_std_hlpod_offline_memory_allocation_snapmat(
    ROM*            rom,
    const int		total_num_nodes,
    const int		n_internal_vertex,
    const double    finish_time,
    const double    dt,
    const int       snapshot_interval,
    const int       num_case,
    const int	    dof);

void ROM_std_hlpod_set_pod_modes_diag(
    ROM* 		rom_v,
    ROM* 		rom_p,
    ROM* 		rom_sups,
    const int 	total_num_nodes,
    const int 	n_internal_vertex,
    const int 	ndof1,
    const int 	ndof2,
    const char* label1,
    const char* label2,
    const char* directory);

void ROM_std_hlpod_read_pod_modes_diag(
    ROM* 		rom_v,
    ROM* 		rom_p,
    ROM* 		rom_sups,
    const int 	total_num_nodes,
    const int 	n_internal_vertex,
    const int 	ndof1,
    const int 	ndof2,
    const char* label1,
    const char* label2,
    const char* directory);

void ROM_std_hlpod_read_pod_modes_Aphi(
        ROM* 		rom_v,
        ROM* 		rom_p,
        ROM* 		rom_sups,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	ndof1,
        const int 	ndof2,
        const char* label1,
        const char* label2,
        const char* directory);

void ROM_std_hlpod_set_pod_modes(
    ROM* 			rom,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		ndof,
    const char* 	label,
    const char* 	directory);

void ROM_std_hlpod_read_pod_modes(
    ROM* 		rom,
    const int 	total_num_nodes,
    const int 	n_internal_vertex,
    const int 	ndof,
    const char* label,
    const char* directory);

void ROM_std_hlpod_pre(
    ROM*         rom,
    const int    total_num_nodes,
    const int    n_internal_vertex,
    const int    dof,
    const char*  metagraph,
    const char*  label_pod_subd,
    const char*	 directory);

void ROM_std_hlpod_set_monolis_comm(
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    const char*	 metagraph_parted0,
    const char*  metagraph,
    const int    solver_type,
    const char*	 directory);

void ROM_std_hlpod_read_metagraph(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* mono_com_rom_solv,
    ROM*		 rom,
    const char*  metagraph,
    const char*	 directory);

void ROM_std_hlpod_pre_lpod_para(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    ROM*		 rom,
    const char*	 metagraph_parted0,
    const char*  metagraph,
    const char*	 directory);

void ROM_std_hlpod_online_pre(
    MONOLIS*     monolis_rom0,
    MONOLIS_COM* mono_com,
    MONOLIS_COM* mono_com_rom,
    MONOLIS_COM* mono_com_rom_solv,
    ROM* 		 rom_sups,
    const int 	 total_num_nodes,
    const char*  metagraph,
    const char*  directory);

void ROM_std_hlpod_calc_reduced_mat(
    MONOLIS*     monolis,
    MONOLIS*     monolis_rom,
    MONOLIS_COM* monolis_com,
    MONOLIS_COM* mono_com0,
    MONOLIS_COM* mono_com_rom,
    ROM*		 rom,
    const int	 total_num_nodes,
    const int	 dof);

void ROM_std_hlpod_solve_ROM(
    MONOLIS*     monolis,
    MONOLIS*     monolis_rom,
    MONOLIS_COM* mono_com_rom,
    ROM*		 rom,
    const int	 total_num_nodes,
    const int	 dof,
    const int 	 mat_max_iter,
    const double mat_epsilon,
    const int 	 label_solver,
    const int	 label_prec);

void ROM_std_hlpod_calc_reduced_mat2(
        MONOLIS*     monolis,
        MONOLIS*     monolis_rom,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com0,
        MONOLIS_COM* mono_com_rom,
        ROM*		 rom,
        const int	 total_num_nodes,
        const int	 dof);
    
void ROM_std_hlpod_solve_ROM_NR(
    MONOLIS*     monolis,
    MONOLIS*     monolis_rom,
    MONOLIS_COM* mono_com_rom,
    ROM*		 rom,
    const int	 total_num_nodes,
    const int	 dof,
    const int 	 mat_max_iter,
    const double mat_epsilon,
    const int 	 label_solver,
    const int	 label_prec);


void ROM_std_hlpod_read_pod_modes_diag_decoupled(
        ROM* 		rom_v,
        ROM* 		rom_p,
        ROM* 		rom_sups,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	ndof1,
        const int 	ndof2,
        const char* label1,
        const char* label2,
        const char* directory);

