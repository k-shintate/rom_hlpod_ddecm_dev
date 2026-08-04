#ifndef CORE_FOM_ST_MPI_H
#define CORE_FOM_ST_MPI_H

#include "core_FOM_ST.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Partition-file contract used by this implementation.
 *
 * Global metadata:
 *   <case>/<ST_MPI_TOP_DIR>/st_meta.dat
 *
 * Rank-local files resolved with
 * monolis_get_global_input_file_name(ST_MPI_TOP_DIR, ST_MPI_PART_DIR, label):
 *
 * graph.dat
 *   num_total_st_nodes
 *   local_row_id num_adj local_col_id...
 *   ... one row for every local ST node, internal nodes first
 *
 * st_node_map.dat
 *   num_internal_st_nodes num_total_st_nodes
 *   local_st_id global_st_id slab_id alpha global_space_node_id
 *
 * elem.dat
 *   num_local_st_cells num_st_nodes_per_cell
 *   local_st_node_id...
 *   ... same local-cell order as st_cell_map.dat
 *
 * st_cell_map.dat
 *   num_local_st_cells num_space_nodes_per_element
 *   local_cell_id global_st_cell_id slab_id global_space_element_id
 *       global_space_node_id...
 *
 * Standard MONOLIS node/communication files must also exist in the same
 * ST partition directory and be consistent with st_node_map.dat.
 *
 * When BBFE reads a partitioned spatial mesh, the following files are
 * additionally needed:
 *
 * space_node_map.dat
 *   num_internal_space_nodes num_total_space_nodes
 *   local_space_node_id global_space_node_id
 *
 * space_elem_map.dat
 *   num_local_space_elements
 *   local_space_element_id global_space_element_id
 *
 * The current main_FOM_ST_MPI.c instead uses a replicated spatial mesh:
 * every rank reads the complete ordinary node.dat and elem.dat.  In that
 * mode ST_space_map_build_identity_from_mesh() is used and the space-map
 * files are not read.
 */

#ifndef ST_MPI_TOP_DIR
#define ST_MPI_TOP_DIR "."
#endif

#ifndef ST_MPI_PART_DIR
#define ST_MPI_PART_DIR "parted"
#endif

typedef struct {
	int64_t* keys;
	int*     values;
	unsigned char* used;
	size_t capacity;
	size_t size;
} ST_ID_MAP;

typedef struct {
	int degree;
	int num_time_dofs;
	int num_window_slabs;

	int num_global_space_nodes;
	int num_global_space_elements;
	int num_space_nodes_per_element;

	int64_t num_global_st_nodes;
	int64_t num_global_st_cells;
} ST_GLOBAL_META;

typedef struct {
	ST_GLOBAL_META meta;

	int num_internal_st_nodes;
	int num_total_st_nodes;

	int64_t* local_to_global_st;
	int* slab_id;
	int* time_dof;
	int* global_space_node_id;
	ST_ID_MAP global_to_local_st;

	int* graph_index;
	int* graph_item;

	int num_local_st_cells;
	int num_st_nodes_per_cell;
	int* cell_st_conn;

	int64_t* global_st_cell_id;
	int* cell_slab_id;
	int* cell_global_space_elem_id;
	int* cell_global_space_conn;

	MONOLIS_COM monolis_com_st;
} ST_PARTITION;

typedef struct {
	int num_internal_space_nodes;
	int num_total_space_nodes;
	int num_local_space_elements;

	int64_t* local_to_global_space_node;
	int64_t* local_to_global_space_elem;

	ST_ID_MAP global_to_local_space_node;
	ST_ID_MAP global_to_local_space_elem;
} ST_SPACE_MAP;

void ST_ID_MAP_initialize(ST_ID_MAP* map, size_t expected_size);
void ST_ID_MAP_insert(ST_ID_MAP* map, int64_t key, int value);
int  ST_ID_MAP_find(const ST_ID_MAP* map, int64_t key);
void ST_ID_MAP_finalize(ST_ID_MAP* map);

void ST_partition_initialize(ST_PARTITION* st);
void ST_partition_read(
		ST_PARTITION* st,
		const char* directory,
		const ST_TIME_ELEMENT* te,
		int num_window_slabs,
		const MONOLIS_COM* parent_com);
void ST_partition_finalize(ST_PARTITION* st);

void ST_space_map_initialize(ST_SPACE_MAP* map);
void ST_space_map_build_identity_from_mesh(
		ST_SPACE_MAP* map,
		const BBFE_DATA* fe);
void ST_space_map_read(
		ST_SPACE_MAP* map,
		const char* directory,
		const char* top_dir,
		const char* part_dir);
void ST_space_map_finalize(ST_SPACE_MAP* map);

int64_t ST_global_gid_partitioned(
		int slab,
		int alpha,
		int global_space_node,
		int num_time_dofs,
		int num_global_space_nodes);

void ST_set_nonzero_pattern_from_partition(
		MONOLIS* monolis,
		const ST_PARTITION* st);

void set_element_mat_ST_dG_window_partitioned(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map);

void set_element_vec_ST_dG_window_partitioned(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		const double* T_prev_trace_global,
		double t_window_start);

void BBFE_sys_monowrap_set_Dirichlet_bc_ST_partitioned_scalar(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BC* bc,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		double* rhs);

void ST_restore_Dirichlet_values_partitioned_scalar(
		BBFE_DATA* fe,
		BBFE_BC* bc,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		double* T_window);

void solver_fom_space_time_partitioned(
		FE_SYSTEM* sys,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		int window_step,
		const double* T_prev_trace_global,
		double* T_window);

void ST_extract_slab_right_trace_global(
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		int slab,
		const double* T_window,
		double* T_global);

void ST_copy_replicated_space_nodal_vector_global(
		const ST_SPACE_MAP* space_map,
		const double* local_vector,
		double* global_vector,
		int num_global_space_nodes);

void ST_copy_global_vector_to_local_space(
		const ST_SPACE_MAP* space_map,
		const double* global_vector,
		double* local_vector);

void ST_validate_partition_against_space_mesh(
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		const BBFE_DATA* fe);

void ST_collect_space_nodal_vector_global(
		const ST_SPACE_MAP* space_map,
		const FE_SYSTEM* sys,
		const double* local_vector,
		double* global_vector,
		int num_global_space_nodes);

#ifdef __cplusplus
}
#endif

#endif
