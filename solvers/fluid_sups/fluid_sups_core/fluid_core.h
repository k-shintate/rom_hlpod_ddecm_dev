#pragma once

#include "monolis.h"

#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include "BBFE/std/integ.h"
#include "BBFE/std/shapefunc.h"
#include "BBFE/std/mapping.h"
#include "BBFE/std/surface.h"

#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "BBFE/sys/monowrap.h"

#include "BBFE/elemmat/set.h"
#include "BBFE/elemmat/equivval.h"
#include "BBFE/elemmat/fluid.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


/* ==========================================================================
 * Global code name
 * ========================================================================== */

extern const char* CODENAME;


/* ==========================================================================
 * Basic fluid utilities
 * ========================================================================== */

const char* BBFE_fluid_get_directory_name(
        int         argc,
        char*       argv[],
        const char* codename);


void BBFE_fluid_pre(
        BBFE_DATA*    fe,
        BBFE_BASIS*   basis,
        int           argc,
        char*         argv[],
        const char*   directory,
        int           num_integ_points_each_axis);


void BBFE_fluid_set_basis(
        BBFE_BASIS*   basis,
        int           local_num_nodes,
        int           num_integ_points_each_axis);


void BBFE_fluid_sups_renew_velocity(
        double**  v,
        double*   ans_vec,
        const int total_num_nodes);


void BBFE_fluid_sups_renew_pressure(
        double*   p,
        double*   ans_vec,
        const int total_num_nodes);


void BBFE_fluid_sups_update_vec(
        double**  v,
        double*   p,
        double*   ans_vec,
        const int total_num_nodes);


void BBFE_fluid_finalize(
        BBFE_DATA*   fe,
        BBFE_BASIS*  basis);


/* ==========================================================================
 * Shared FOM / NR / ST-FOM API
 * ========================================================================== */

void BBFE_fluid_initialize_velocity_pressure(
        double**    v,
        double*     p,
        const int   total_num_nodes);


void BBFE_fluid_initialize_velocity_pressure_karman_vortex(
        double**    v,
        double**    v_old,
        double*     p,
        const int   total_num_nodes);


void BBFE_fluid_sups_add_velocity_pressure(
        double**    v,
        double*     p,
        double*     sol_vec,
        const int   total_num_nodes);


void BBFE_fluid_add_Dbc(
        double*         ansvec,
        BBFE_BC*        bc,
        const int       total_num_nodes,
        const int       dof);


void BBFE_fluid_update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes);


double BBFE_fluid_calc_internal_norm_2d(
        double**    in,
        const int   num1,
        const int   num2);


double BBFE_fluid_calc_internal_norm_1d(
        double*     in,
        const int   num1,
        const int   num2);


void BBFE_fluid_vec_copy_2d(
        double**    in,
        double**    out,
        const int   num1,
        const int   num2);


void BBFE_fluid_vec_copy_1d(
        double*       dst,
        const double* src,
        int           size);


double BBFE_fluid_calc_internal_norm_total(
        const double* rvec,
        int           n_internal_vertex);


double BBFE_fluid_calc_internal_norm_momentum(
        const double* rvec,
        int           n_internal_vertex);


double BBFE_fluid_calc_internal_norm_mass(
        const double* rvec,
        int           n_internal_vertex);