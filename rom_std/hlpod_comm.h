#pragma once

#include "rom_dataset.h"
#include "std.h"

void ROM_std_hlpod_get_neib_vec(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    const int 		num_modes,
    const int       ndof);

void ROM_std_hlpod_get_neib_vec_save_memory(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*		hlpod_vals,
    const int 		num_modes);

void ROM_std_hlpod_set_comm_table_para_subd(
    MONOLIS_COM*	monolis_com,
    MONOLIS_COM*	mono_com_rom,
    const int 		n_neib);

void ROM_std_hlpod_get_neib_num_modes_para_subd(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    const int       np,
    const int       num_my_modes);

void ROM_std_hlpod_get_neib_num_modes_mode_subd(
    MONOLIS_COM*  	mono_com,
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		np,     //変更
    const char*     directory);
