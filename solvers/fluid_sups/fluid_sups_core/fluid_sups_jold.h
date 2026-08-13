#pragma once

#include <math.h>

/*
 * Analytic previous-state tangent for the verified SUPS/PSPG transient
 * structure.
 *
 * This kernel is the derivative of the element residual integrand with
 * respect to one previous-time nodal velocity coefficient U_old(j,b).
 *
 * The verified caller constructs
 *
 *   du_time = velocity - velocity_old
 *
 * and tau is evaluated from the CURRENT velocity. Therefore:
 *
 *   d tau / d U_old = 0,
 *   d tau_c / d U_old = 0.
 *
 * For the residual structure used by the current SUPS implementation, the
 * previous velocity enters only through the transient momentum residual.
 *
 * Let
 *
 *   a_i = velocity . grad_N_i.
 *
 * IMPORTANT: the linked verified residual provider receives
 *
 *   du_time = v - v_old
 *
 * as the time-increment quantity.  Its old-state derivative therefore has
 * NO additional 1/dt factor here.  The sampled runtime verification showed
 * exactly the expected factors:
 *
 *   momentum: previous formula / FD = 1/dt
 *   pressure: previous formula / FD = density/dt.
 *
 * Hence the provider-consistent old-state tangent is
 *
 *   momentum rows:
 *     J_old[a][b]
 *       = -density * N_j
 *         * ( N_i + tau*density*a_i ) delta_ab
 *
 *   pressure/continuity row:
 *     J_old[3][b]
 *       = -N_j * tau * grad_N_i[b].
 *
 * The tau_c / LSIC term contains the current divergence only and therefore
 * has no previous-state derivative.
 *
 * IMPORTANT:
 * This analytic expression is accompanied by a runtime finite-difference
 * verification path in fluid_sups_st_model.c.  The first run should use
 *
 *   FLUID_ST_JOLD_MODE=verify
 *
 * so the expression is checked directly against
 * BBFE_elemmat_fluid_sups_vec_NR from the linked verified provider.
 */
static inline void BBFE_elemmat_fluid_sups_mat_NR_old(
    double J_old[4][3],
    double N_i,
    double N_j,
    const double grad_N_i[3],
    const double velocity[3],
    double density,
    double tau,
    double dt)
{
    const double advective_test =
        velocity[0] * grad_N_i[0]
        + velocity[1] * grad_N_i[1]
        + velocity[2] * grad_N_i[2];

    /*
     * Provider-consistent transient scaling.
     *
     * BBFE_elemmat_fluid_sups_vec_NR is called with
     *   du_time = v - v_old
     * and the verified finite-difference derivative demonstrates that the
     * residual kernel does not apply an extra 1/dt to this old-state
     * direction here.
     */
    const double previous_transient_trial_momentum =
        -density * N_j;

    const double previous_transient_trial_pspg =
        -N_j;

    const double momentum_test =
        N_i + tau * density * advective_test;

    (void)dt;

    int a;
    int b;

    for(a = 0; a < 4; a++){
        for(b = 0; b < 3; b++){
            J_old[a][b] = 0.0;
        }
    }

    for(b = 0; b < 3; b++){
        J_old[b][b] =
            previous_transient_trial_momentum
            * momentum_test;

        J_old[3][b] =
            previous_transient_trial_pspg
            * tau
            * grad_N_i[b];
    }
}
