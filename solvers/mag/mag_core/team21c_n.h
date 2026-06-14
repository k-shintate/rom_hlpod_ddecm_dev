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