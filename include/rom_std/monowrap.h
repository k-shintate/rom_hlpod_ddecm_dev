
#pragma once

#include "monolis.h"

void ROM_monowrap_solve(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon);
