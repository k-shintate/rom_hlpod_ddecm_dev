#pragma once


#include "3ph_tr_NR.h"

void set_element_mat_nedelec_Aphi_team21a0(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned);

void set_element_vec_nedelec_Aphi_team21a0(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned);

void apply_dirichlet_bc_for_A_and_phi_team21a0(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);

void apply_dirichlet_bc_ned(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);

void set_element_mat_nedelec_Aphi_team21a0_debug(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned);

void set_element_mat_NR_Aphi_team21a02(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt);

void set_element_vec_NR_Aphi_team21a02(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time);

void apply_dirichlet_bc_ned_R(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);

void add_phi_gauge_bc_team21a0(
MONOLIS* monolis,
	    	BBFE_DATA* fe,
    NEDELEC*   ned,
    BBFE_BC*   bc);

void set_element_mat_NR_Aphi_team21a03(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt);

void set_element_vec_NR_Aphi_team21a03(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time);

