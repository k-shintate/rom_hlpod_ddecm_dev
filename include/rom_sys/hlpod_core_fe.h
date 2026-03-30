
#pragma once

#include "monolis.h"

#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include "BBFE/std/integ.h"
#include "BBFE/std/shapefunc.h"
#include "BBFE/std/mapping.h"
#include "BBFE/std/surface.h"

#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "BBFE/sys/monowrap.h"

#include "BBFE/elemmat/set.h"
#include "BBFE/elemmat/equivval.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


#include "BBFE/sys/FE_dataset.h"

#include "std.h"
#include "rom_dataset.h"
#include "write_fe.h"
#include "set_D_bc.h"


void ROM_sys_hlpod_fe_set_snap_mat_para(
    double*       	comp_vec,
    HLPOD_MAT*      hlpod_mat,
    BBFE_BC*      	bc,
    ROM_BC*		    rom_bc,
    const int 		total_num_nodes,
    const int		dof,
    const int 		count);

void ROM_sys_hlpod_fe_add_Dbc(
    double*         ansvec,
    BBFE_BC*      	bc,
    const int 		total_num_nodes,
    const int       dof);

void ROM_sys_hlpod_fe_write_pod_modes_vtk(
    BBFE_DATA*		fe,
    ROM*			rom,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int 		ndof,
    const char*     label,
    const char*		directory);

void ROM_sys_hlpod_fe_write_pod_modes_vtk_diag(
    BBFE_DATA*		fe,
    ROM*			rom,
    const int 		n_internal_vertex,
    const int 		num_modes1,
    const int 		num_modes2,
    const int		ndof1,
    const int 		ndof2,
    const char*     label_v,
    const char*     label_p,
    const char*		directory);

double ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
    BBFE_DATA*    	fe,
    BBFE_BASIS*   	basis,
    MONOLIS_COM*  	monolis_com,
    double        	t,
    const double* 	comp_vec,
    const double* 	pod_comp_vec);

double ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    MONOLIS_COM*    monolis_com,
    double          t,
    const double**  fem_vec,
    const double**  rom_vec);
