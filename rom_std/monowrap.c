
#include "monowrap.h"

void ROM_monowrap_solve(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon)
{
    monolis_set_method   (monolis, solver_type);
    monolis_set_precond  (monolis, precond_type);
    monolis_set_maxiter  (monolis, num_max_iters);
    monolis_set_tolerance(monolis, epsilon);
    monolis_show_iterlog (monolis, true);

    monolis_solve_R(
            monolis,
            monolis_com,
            monolis->mat.R.B,
            sol_vec);
}

