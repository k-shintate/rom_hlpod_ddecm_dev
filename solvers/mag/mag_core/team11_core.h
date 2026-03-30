#include "convdiff_core.h"
#include "nedelec_core.h"
#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include <limits.h>


/* ======= 例：t=1,5,10,20 ms を出力 ======= */
void output_B_node_vtk_team11_series(BBFE_DATA* fe, const char* dir);

void apply_nedelec_boundary_conditions_TEAM11_time(
    MONOLIS*    monolis,
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    double      Bz_t,
    double      t_sec);

void apply_nedelec_boundary_conditions_TEAM11_NR(
    MONOLIS*    monolis,
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    double      Bz_t,
    double      t_sec,
    double*      val);