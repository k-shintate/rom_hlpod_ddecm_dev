#pragma once

#include "BBFE/sys/FE_dataset.h"

#include "mag_dataset.h"

#include "set_reduced_matvec.h"

#include "./../mag_core/convdiff_core.h"
#include "./../mag_core/nedelec_core.h"

#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include <complex.h>

#include "mag_dataset.h"
#include "core_ROM.h"

void set_reduced_mat_global_para(
    MONOLIS* monolis,
    MONOLIS_COM* monolis_com,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES*         vals,
    NEDELEC* ned,
    double**        mat,
    double**        pod_modes,
    const int 		num_modes,
    const double* x_curr,
    double dt,
    const char* directory);

void set_reduced_vec_global_para(
    MONOLIS* monolis,
    MONOLIS_COM* monolis_com,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES*         vals,
    NEDELEC* ned,
    double*         rhs,
    double**        pod_modes,
    const int		num_modes,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time,
    const char* directory);