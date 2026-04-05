
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


double calc_copper_shield_loss_EM1(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt);

void log_copper_shield_loss_EM1(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    double shield_loss_inst);


void log_copper_shield_loss_EM1_cycle_average(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    double shield_loss_inst);

typedef struct {
    double loss_total;   /* ∫ sigma |E|^2 dV */
    double loss_A;       /* ∫ sigma |dA/dt|^2 dV */
    double loss_phi;     /* ∫ sigma |grad phi|^2 dV */
    double loss_cross;   /* ∫ 2 sigma (dA/dt)·(grad phi) dV */

    double vol_shield;   /* shield volume */
    double rms_dA_dt;    /* sqrt( (1/V) ∫ |dA/dt|^2 dV ) */
    double rms_grad_phi; /* sqrt( (1/V) ∫ |grad phi|^2 dV ) */
    double rms_E;        /* sqrt( (1/V) ∫ |E|^2 dV ) */

    double max_dA_dt;
    double max_grad_phi;
    double max_E;
} SHIELD_LOSS_DIAG;

SHIELD_LOSS_DIAG calc_copper_shield_loss_EM1_diag(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt
);


void log_copper_shield_loss_EM1_diag(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    const SHIELD_LOSS_DIAG* d);

void log_coil_ampere_turn_diag_team21c(
    FE_SYSTEM* sys,
    int step,
    double current_time);


typedef struct {
    double int_dA2;      /* ∫ |dA/dt|^2 dV */
    double int_gradphi2; /* ∫ |grad phi|^2 dV */
    double int_cross;    /* ∫ (dA/dt · grad phi) dV */
    double vol_shield;   /* ∫ 1 dV */
    double rho;          /* correlation coefficient */
} SHIELD_FIELD_INT_DIAG;

SHIELD_FIELD_INT_DIAG calc_shield_field_integrals_EM1(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    double dt,
    int use_phi_in_shield   /* 1: grad phi を含める, 0: grad phi = 0 */
);

void log_shield_field_integrals_EM1(
    FE_SYSTEM* sys,
    int step,
    double t,
    double dt,
    const SHIELD_FIELD_INT_DIAG* d
);


void apply_dirichlet_bc_for_A_and_phi_team21a2(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned);