#pragma once

#include "BBFE/sys/FE_dataset.h"

#include "fluid_fs_dataset.h"
#include "set_reduced_matvec.h"

/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
    MONOLIS*        monolis,
    MONOLIS_COM*	monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    double**        mat,
    double**        pod_modes,
    const int 		num_modes);

void set_D_bc_global_para(
    MONOLIS*        monolis,
    MONOLIS_COM*    monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*     	bc,
    double*         rhs,
    double**        pod_modes,
    const int 		num_modes);

void set_reduced_vec_global_para(
    MONOLIS*        monolis,        
    MONOLIS_COM*	monolis_com,
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    double*         rhs,
    double**        pod_modes,
    const int		num_modes);

void allreduce_global_para(
    MONOLIS_COM*	monolis_com,
    double**        mat,
    double*         rhs,
    const int		num_modes);
