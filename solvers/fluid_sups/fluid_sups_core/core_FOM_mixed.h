#pragma once

#include "core_FOM.h"

/*
 * Mixed TET4 + PRI6 full-order fluid system.
 *
 * Tetrahedra and prisms share exactly the same node numbering and nodal
 * unknowns (u,v,w,p).  They therefore assemble into one 4-DOF/node Monolis
 * matrix.  The two BBFE_DATA objects own separate element connectivities and
 * geometry arrays but read the same node.dat.
 */
typedef struct
{
    BBFE_BASIS  basis_tet;
    BBFE_DATA   fe_tet;

    BBFE_BASIS  basis_pri;
    BBFE_DATA   fe_pri;

    CONDITIONS  cond;
    VALUES      vals;

    BBFE_BC     bc;
    BBFE_BC     bc_NR;

    MONOLIS     monolis;
    MONOLIS_COM mono_com;

    /* Optional surface data.  It is not initialized by the mixed pre routine;
       keep the existing surface preprocessing path if force diagnostics are
       required. */
    BBFE_BASIS  basis_surf;
    BBFE_DATA   surf;

} FE_SYSTEM_FLUID_MIXED;


/*
 * A partition does not necessarily contain both element types.
 * Treat an element family as active only when it has both a valid local
 * node count and at least one local element.
 */
static inline int BBFE_fluid_mixed_has_elements(const BBFE_DATA* fe)
{
    return (fe != NULL && fe->local_num_nodes > 0 && fe->total_num_elems > 0);
}

#define FLUID_MIXED_BLOCK_SIZE 4
#define FLUID_MIXED_ELEM_TET_FILE "graph_elem_tet.dat"
#define FLUID_MIXED_ELEM_PRI_FILE "graph_elem_prism.dat"
#define FLUID_MIXED_GRAPH_FILE    "graph.dat"

/*
 * Implemented in the existing core_FOM.c, but not declared by the
 * original core_FOM.h.  The mixed surface preprocessor calls it from a
 * separate translation unit, so an explicit prototype is required
 * (especially when the .c files are compiled as C++).
 */
void BBFE_convdiff_set_basis_surface(
    BBFE_BASIS* basis,
    int         local_num_nodes,
    int         num_integ_points_each_axis);


/* MPI overlap synchronization / diagnostics */
void BBFE_fluid_mixed_sync_solution_vector(
    MONOLIS_COM* mono_com,
    int          total_num_nodes,
    double*      vec4);

void BBFE_fluid_mixed_sync_state(
    MONOLIS_COM* mono_com,
    VALUES*      vals,
    int          total_num_nodes);

void BBFE_fluid_mixed_report_prism_bc(
    const FE_SYSTEM_FLUID_MIXED* sys);

void BBFE_fluid_mixed_report_prism_components(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int step,
    int newton_iteration);

void BBFE_fluid_mixed_report_prism_partition(
    const BBFE_DATA*   fe_pri,
    const MONOLIS_COM* mono_com);

void BBFE_fluid_mixed_report_prism_update(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int                          step,
    int                          newton_iteration);

/* Physical-owner Newton increment diagnostics over the full mixed domain. */
void BBFE_fluid_mixed_report_newton_diagnostics(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int                          step,
    int                          newton_iteration);

void BBFE_fluid_mixed_check_positive_jacobian(
    const BBFE_DATA*  fe,
    const BBFE_BASIS* basis,
    const char*       element_name,
    double            tolerance);

/* Empty surf_graph.dat.<rank> files containing only 0 are valid. */
void BBFE_fluid_pre_surface_mixed(
    BBFE_DATA*   surf,
    BBFE_BASIS*  basis,
    const char*  directory,
    const char*  filename,
    int          num_integ_points_each_axis);

/*
 * Read the same surface file as BBFE_fluid_pre_surface_mixed(), but keep only
 * the first N_internal entries described by
 * parted.0/<filename>.n_internal.<rank>.  Use this for physical SUM-type
 * diagnostics/forces where each face must occur once globally.
 */
void BBFE_fluid_pre_surface_internal_mixed(
    BBFE_DATA*   surf,
    BBFE_BASIS*  basis,
    const char*  directory,
    const char*  filename,
    int          num_integ_points_each_axis);

void BBFE_fluid_finalize_surface_mixed(
    BBFE_DATA*  surf,
    BBFE_BASIS* basis);

/* PRI6 reference-element support */
int BBFE_fluid_mixed_integ_prism_set_arbitrary_points(
    int      num_points_in_each_axis,
    double** integ_point,
    double*  integ_weight);

void BBFE_fluid_mixed_shapefunc_prism1st_get_val(
    const double xi[3],
    double*      N);

void BBFE_fluid_mixed_shapefunc_prism1st_get_derivative(
    const double xi[3],
    double*      dN_dxi,
    double*      dN_deta,
    double*      dN_dzeta);

void BBFE_fluid_set_basis_prism(
    BBFE_BASIS* basis,
    int         num_integ_points_each_axis);

/* Mixed mesh input / geometry */
void BBFE_fluid_pre_mixed(
    BBFE_DATA*   fe_tet,
    BBFE_BASIS*  basis_tet,
    BBFE_DATA*   fe_pri,
    BBFE_BASIS*  basis_pri,
    int          argc,
    char*        argv[],
    const char*  directory,
    int          num_integ_points_each_axis);

/* Read the unified tetra+prism graph and create one 4x4-block sparse matrix. */
void BBFE_fluid_init_monomat_mixed(
    MONOLIS*          monolis,
    MONOLIS_COM*      monolis_com,
    const BBFE_DATA*  fe_tet,
    const BBFE_DATA*  fe_pri,
    const char*       directory);

/* Mixed SUPS/NR assembly.  Existing element kernels are called twice and
   accumulate into the same Monolis matrix/RHS. */
void set_element_mat_NR_Tezuer_mixed(
    MONOLIS*    monolis,
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri,
    VALUES*     vals);


/* Mixed VMS/NR assembly. The element kernels themselves are provided by
 * core_NR; these wrappers accumulate TET4 and PRI6 contributions into one
 * unified Monolis matrix/RHS. */
void set_element_mat_NR_linear_VMS_mixed(
    MONOLIS*    monolis,
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri,
    VALUES*     vals);

void set_element_vec_NR_linear_VMS_mixed(
    MONOLIS*    monolis,
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri,
    VALUES*     vals);

void set_element_mat_NR_nonlinear_VMS_mixed(
    MONOLIS*    monolis,
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri,
    VALUES*     vals);

void set_element_vec_NR_nonlinear_VMS_mixed(
    MONOLIS*    monolis,
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri,
    VALUES*     vals);

void solver_fom_NR_mixed(
    FE_SYSTEM_FLUID_MIXED* sys,
    double                 t,
    int                    step);

/* One legacy VTK file containing both VTK_TETRA(10) and VTK_WEDGE(13). */
void output_result_file_vtk_mixed(
    const BBFE_DATA* fe_tet,
    const BBFE_DATA* fe_pri,
    const VALUES*    vals,
    const char*      filename,
    const char*      directory,
    double           t);

void output_files_mixed(
    FE_SYSTEM_FLUID_MIXED* sys,
    int                    file_num,
    double                 t);

void BBFE_fluid_finalize_mixed(
    BBFE_DATA*  fe_tet,
    BBFE_BASIS* basis_tet,
    BBFE_DATA*  fe_pri,
    BBFE_BASIS* basis_pri);

