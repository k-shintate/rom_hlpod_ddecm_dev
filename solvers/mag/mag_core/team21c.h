
#include "3ph_tr_NR.h"

void set_element_mat_NR_Aphi_team21c(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt
);

void set_element_vec_NR_Aphi_team21c(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
);

void apply_dirichlet_bc_for_A_and_phi_team21c(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);
