#ifndef CORE_FOM_ST_SPACEBLOCK_H
#define CORE_FOM_ST_SPACEBLOCK_H

/*
 * Windowed space-time dG implementation with spatial-only MPI partitioning.
 *
 * MONOLIS node:
 *   one ordinary spatial mesh node
 *
 * Degrees of freedom at each MONOLIS node:
 *   q = slab * te->n_dof + alpha
 *
 * Vector layout:
 *   value[space_node * num_dofs_per_node + q]
 *
 * The spatial nodal graph is reconstructed from the partitioned FE mesh.
 */

#include <stddef.h>
#include "core_FOM_ST.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int num_nodes;
    int num_items;
    int* index;
    int* item;
} STSB_SPATIAL_GRAPH;

typedef struct {
    double residual_all2;
    double rhs_all2;
    double residual_free2;
    double rhs_free2;
    double residual_dirichlet2;
    double rhs_dirichlet2;
} STSB_CONSTRAINED_RESIDUAL_SUM;

/*
 * Reference FOM initialization shared by the verified spatial-block FOM and
 * ST-DDROM collect/online paths.  The spatial graph is always reconstructed
 * from the FE connectivity so both executables use exactly the same matrix
 * nonzero pattern and assembly route.
 */
void STSB_initialize_reference_fom(
    FE_SYSTEM* sys,
    ST_TIME_ELEMENT* te,
    STSB_SPATIAL_GRAPH* graph,
    int argc,
    char* argv[],
    int degree,
    int num_window_slabs,
    int* num_windows,
    int* num_internal_space_nodes);

void STSB_finalize_reference_fom(
    FE_SYSTEM* sys,
    STSB_SPATIAL_GRAPH* graph);

void STSB_spatial_graph_build_from_mesh(
    STSB_SPATIAL_GRAPH* graph,
    const BBFE_DATA* fe);

/* Diagnostic-only helpers.  They never change the FOM solve path. */
int STSB_evaluate_constrained_residual(
    FE_SYSTEM* sys,
    int st_dof_per_space_node,
    const double* window_solution,
    STSB_CONSTRAINED_RESIDUAL_SUM* global_sum);

int STSB_compute_window_checksum(
    const FE_SYSTEM* sys,
    int st_dof_per_space_node,
    const double* window_solution,
    const double* right_trace,
    double* solution_l2,
    double* trace_l2);

int STSB_get_num_dofs_per_node(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs);

int STSB_get_block_dof(
    int slab,
    int alpha,
    const ST_TIME_ELEMENT* te);

size_t STSB_get_vector_index(
    int spatial_node,
    int slab,
    int alpha,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs);

void STSB_validate_configuration(
    const BBFE_DATA* fe,
    const VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs);

void STSB_spatial_graph_initialize(
    STSB_SPATIAL_GRAPH* graph);

void STSB_spatial_graph_read(
    STSB_SPATIAL_GRAPH* graph,
    const char* label,
    const char* directory);

void STSB_spatial_graph_finalize(
    STSB_SPATIAL_GRAPH* graph);

void STSB_validate_spatial_graph_against_mesh(
    const STSB_SPATIAL_GRAPH* graph,
    const BBFE_DATA* fe);

void STSB_set_nonzero_pattern_from_spatial_graph(
    MONOLIS* monolis,
    const STSB_SPATIAL_GRAPH* graph,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_internal_space_nodes);

void STSB_compute_element_mass_diffusion(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    int elem,
    double* Me,
    double* Ke);

void STSB_compute_element_source(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    const ST_TIME_ELEMENT* te,
    int elem,
    double slab_start_time,
    double dt,
    double* Fe);

void STSB_assemble_window_matrix(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs);

void STSB_add_previous_trace_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    const double* previous_trace);

void STSB_add_source_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time);

void STSB_assemble_window_rhs(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    const double* previous_trace);

void STSB_apply_dirichlet_bc(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    double* rhs);

void STSB_restore_dirichlet_values(
    BBFE_DATA* fe,
    BBFE_BC* bc,
    VALUES* vals,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    double* window_solution);

void STSB_set_previous_trace_from_nodal_value(
    const BBFE_DATA* fe,
    const double* nodal_value,
    double* previous_trace);

void STSB_extract_slab_right_trace(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int slab,
    int num_space_nodes,
    const double* window_solution,
    double* nodal_trace);

void STSB_extract_last_right_trace(
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    int num_space_nodes,
    const double* window_solution,
    double* previous_trace);

void STSB_solve_window(
    FE_SYSTEM* sys,
    const ST_TIME_ELEMENT* te,
    int num_window_slabs,
    double window_start_time,
    int window_id,
    const double* previous_trace,
    double* window_solution);

int STSB_get_num_windows(
    double finish_time,
    double dt,
    int num_window_slabs);

#ifdef __cplusplus
}
#endif

#endif
