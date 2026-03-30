#pragma once

#include "convdiff_core.h"
#include "hlpod_dataset.h"
#include "mag_dataset.h"

#include "core_FOM.h"


void solver_fom_NR_Aphi_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total);

void solver_fom_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total);