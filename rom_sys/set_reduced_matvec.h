#pragma once

#include "BBFE/sys/FE_dataset.h"

void set_element_reduced_rhs(
    double*			vec,
    double**		pod_modes,
    const int		index,
    const double	integ_val,
    const int		num_modes);

void set_element_reduced_rhs_Dbc(
    BBFE_BC*     	bc,
    double*			vec,
    double**		pod_modes,
    const int		index_i,
    const int		index_j,
    const double	integ_val,
    const int		num_modes);

void set_element_reduced_mat(
    double**		mat,
    double**		pod_modes,
    const int		index_i,
    const int		index_j,
    const double	integ_val,
    const int		num_modes);
