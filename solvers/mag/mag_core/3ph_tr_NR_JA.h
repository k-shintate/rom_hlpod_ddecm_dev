#include "3ph_tr_NR.h"

void ja_init_states_zero(void);

void ja_allocate_states(int n_elem, int n_ip);

void initialize_g_phase_current(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time);

int solver_fom_NR_Aphi_JA(
    FE_SYSTEM* sys,
    double* x_prev,
    double* x_curr,
    int n_dof_total,
    double dt,
    double current_time);

void log_accuracy_metrics_JA(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    int step,
    double t,
    double dt,
    double P_input);