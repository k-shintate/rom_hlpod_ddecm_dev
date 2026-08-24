#pragma once

/*
 * Mixed tetrahedron/prism extension for the FOM convection-diffusion solver.
 *
 * This file intentionally does not modify FE_SYSTEM in diff_dataset.h.
 * The existing tetra/hexa FOM remains usable, while FE_SYSTEM_MIXED is used by
 * main_FOM_mixed.c.
 */

#include "core_FOM.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FOM_MIXED_ELEM_TET_FILE   "graph_elem_tet.dat"
#define FOM_MIXED_ELEM_PRI_FILE   "graph_elem_prism.dat"
#define FOM_MIXED_GRAPH_FILE      "graph.dat"

/*
 * node.dat is shared by both element families.  fe_tet and fe_pri each read
 * the file so that both BBFE_DATA objects own valid coordinate arrays and can
 * be finalized independently.
 */
typedef struct
{
    BBFE_BASIS   basis_tet;
    BBFE_DATA    fe_tet;

    BBFE_BASIS   basis_pri;
    BBFE_DATA    fe_pri;

    BBFE_BC      bc;

    MONOLIS      monolis;
    MONOLIS      monolis0;
    MONOLIS_COM  monolis_com;

    CONDITIONS   cond;
    VALUES       vals;

} FE_SYSTEM_MIXED;

/* ---------- prism reference element ---------- */

int FOM_mixed_integ_prism_set_arbitrary_points(
        int      num_points_in_each_axis,
        double** integ_point,
        double*  integ_weight);

void FOM_mixed_shapefunc_prism1st_get_val(
        const double xi[3],
        double*      N);

void FOM_mixed_shapefunc_prism1st_get_derivative(
        const double xi[3],
        double*      dN_dxi,
        double*      dN_deta,
        double*      dN_dzeta);

void BBFE_convdiff_set_basis_prism(
        BBFE_BASIS* basis,
        int         num_integ_points_each_axis);

/* ---------- mixed mesh input / graph ---------- */

void BBFE_convdiff_pre_mixed(
        BBFE_DATA*    fe_tet,
        BBFE_BASIS*   basis_tet,
        BBFE_DATA*    fe_pri,
        BBFE_BASIS*   basis_pri,
        BBFE_BC*      bc,
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        const char*   directory,
        int           num_integ_points_each_axis);

void BBFE_convdiff_set_nonzero_pattern_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_ref,
        const char* directory);

/* ---------- assembly / solve ---------- */

void set_element_mat_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        VALUES*     vals);

void set_element_vec_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        VALUES*     vals,
        double      t);

void solver_fom_mixed(
        FE_SYSTEM_MIXED* sys,
        double           t,
        int              step);

/* ---------- output / error ---------- */

void output_result_file_vtk_mixed(
        BBFE_DATA*  fe_tet,
        BBFE_DATA*  fe_pri,
        VALUES*     vals,
        const char* filename,
        const char* directory,
        double      t);

void output_files_mixed(
        FE_SYSTEM_MIXED* sys,
        int              file_num,
        double           t);

double BBFE_convdiff_equivval_relative_L2_error_scalar_mixed(
        BBFE_DATA*    fe_tet,
        BBFE_BASIS*   basis_tet,
        BBFE_DATA*    fe_pri,
        BBFE_BASIS*   basis_pri,
        MONOLIS_COM*  monolis_com,
        double        t,
        const double* comp_vec,
        double        (*func)(double, double, double, double));

/* ---------- cleanup ---------- */

void BBFE_convdiff_finalize_mixed(
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        BBFE_BC*    bc);

#ifdef __cplusplus
}
#endif
