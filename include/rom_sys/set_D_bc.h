#pragma once

#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

#include "BBFE/sys/FE_dataset.h"

#include "hlpod_pre.h"
#include "rom_dataset.h"
#include "write_fe.h"


void ROM_sys_hlpod_fe_set_bc_id(
    BBFE_BC*   	bc,
    const int   total_num_nodes,
    const int   num_dofs_on_node,
    ROM_BC*	    rom_bc);

void ROM_sys_hlpod_fe_manusol_set_bc_id(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    const int   num_dofs_on_node,
    double     	t,
    double      (*func)(double, double, double, double),
    ROM_BC*	    rom_bc);

void ROM_sys_hlpod_fe_monowrap_set_Dirichlet_bc(
    MONOLIS*    monolis,
    const int   num_dofs_on_node,
    BBFE_BC*    bc,
    double*     g_rhs,
    ROM_BC*	    rom_bc);

void HROM_manusol_set_bc(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    const int   num_dofs_on_node,
    double     	t,
    double      (*func)(double, double, double, double), // scalar function(x, y, z, t)
    HLPOD_MAT*	hlpod_mat);
