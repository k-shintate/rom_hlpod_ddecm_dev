#pragma once

#include "convdiff_core.h"

double BBFE_elemmat_mag_mat_curl(
    const double *curlNi,
    const double *curlNj,
    double k);

double BBFE_elemmat_mag_mat_mass(
    const double *Ni,
    const double *Nj,
    double a);