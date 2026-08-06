
#pragma once

#include "fluid_core.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

void set_wall_face_vecmat_symmetric_nitsche(
    MONOLIS* monolis,
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_eff,
    const double Uw[3],
    double gamma_n);
