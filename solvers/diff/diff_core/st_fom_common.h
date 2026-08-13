#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ST_FOM_MAX_TIME_DOF 3
#define ST_FOM_MAX_TIME_IP  3
#define ST_FOM_WINDOW_MAGIC 0x53544657u /* "STFW" */
#define ST_FOM_WINDOW_VERSION 1u

typedef enum {
    ST_FOM_SUCCESS = 0,
    ST_FOM_ERR_ARGUMENT = 1,
    ST_FOM_ERR_DIMENSION = 2,
    ST_FOM_ERR_ALLOCATION = 3,
    ST_FOM_ERR_IO = 4,
    ST_FOM_ERR_CALLBACK = 5,
    ST_FOM_ERR_UNSUPPORTED = 6
} ST_FOM_STATUS;

typedef enum {
    ST_FOM_TIME_FAMILY_DG = 1
} ST_FOM_TIME_FAMILY;

typedef struct {
    int family;
    int degree;
    int n_dof;
    int is_nodal;
    int num_time_ip;

    double dof_tau[ST_FOM_MAX_TIME_DOF];
    double e_minus[ST_FOM_MAX_TIME_DOF];
    double e_plus[ST_FOM_MAX_TIME_DOF];

    double Mt[ST_FOM_MAX_TIME_DOF][ST_FOM_MAX_TIME_DOF];
    double Dt[ST_FOM_MAX_TIME_DOF][ST_FOM_MAX_TIME_DOF];
    double E_self[ST_FOM_MAX_TIME_DOF][ST_FOM_MAX_TIME_DOF];
    double E_prev[ST_FOM_MAX_TIME_DOF][ST_FOM_MAX_TIME_DOF];

    double tau[ST_FOM_MAX_TIME_IP];
    double weight[ST_FOM_MAX_TIME_IP];
    double psi[ST_FOM_MAX_TIME_DOF][ST_FOM_MAX_TIME_IP];
} ST_FOM_TIME_ELEMENT;

typedef struct {
    int num_space_points;
    int num_owned_space_points;
    int physical_dof;
    int num_slabs;
    ST_FOM_TIME_ELEMENT time;
} ST_FOM_WINDOW_LAYOUT;

/*
 * Node-major ordering used by the distributed ST-FOM/ST-DDROM implementation:
 *
 *   [space point][slab][time basis][physical component].
 */
int ST_fom_block_dof(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    int time_dof,
    int component);

size_t ST_fom_vector_index(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int space_point,
    int slab,
    int time_dof,
    int component);

int ST_fom_window_dof_per_space_point(
    const ST_FOM_WINDOW_LAYOUT* layout);

size_t ST_fom_window_vector_length(
    const ST_FOM_WINDOW_LAYOUT* layout);

int ST_fom_time_element_init_dg(
    ST_FOM_TIME_ELEMENT* time,
    int degree);

int ST_fom_layout_initialize(
    ST_FOM_WINDOW_LAYOUT* layout,
    int num_space_points,
    int num_owned_space_points,
    int physical_dof,
    int num_slabs,
    int degree);

/*
 * Resolve the number of owned spatial points for serial/MPI execution.
 *
 * MONOLIS communication tables may report n_internal_vertex == 0 in a
 * one-process run because no distributed ownership table is required.
 * In that case every local spatial point is owned.  In a distributed run,
 * the communication-table count is authoritative.
 */
int ST_fom_resolve_owned_space_points(
    int num_local_space_points,
    int communicator_size,
    int communication_owned_space_points,
    int* num_owned_space_points);

/* Pack one physical dG(0) slab state into a node-major ST window vector. */
int ST_fom_pack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* state,
    double* window_state);

int ST_fom_unpack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* state);

/* Extract the right trace of an arbitrary nodal dG time element. */
int ST_fom_extract_slab_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* trace);

int ST_fom_extract_last_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    double* trace);

/*
 * Generic exact causal dG(0) window runner.
 *
 * A nonlinear or linear one-step solver is supplied as a callback.  The runner
 * groups consecutive verified one-step solutions into a space-time window.
 * For dG(0), this is algebraically equivalent to forward substitution of the
 * block-lower-triangular monolithic space-time system.
 */
typedef int (*ST_FOM_SOLVE_STEP_CALLBACK)(
    void* user_context,
    int global_step,
    double time_old,
    double time_new);

typedef int (*ST_FOM_EXPORT_STATE_CALLBACK)(
    void* user_context,
    double* state /* [num_space_points][physical_dof] */);

typedef int (*ST_FOM_WRITE_OUTPUT_CALLBACK)(
    void* user_context,
    int file_number,
    int global_step,
    double time);

typedef int (*ST_FOM_WINDOW_COMPLETE_CALLBACK)(
    void* user_context,
    int window_id,
    double time_start,
    double time_end,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    const double* right_trace);

typedef struct {
    void* user_context;
    ST_FOM_SOLVE_STEP_CALLBACK solve_step;
    ST_FOM_EXPORT_STATE_CALLBACK export_state;
    ST_FOM_WRITE_OUTPUT_CALLBACK write_output;
    ST_FOM_WINDOW_COMPLETE_CALLBACK window_complete;
} ST_FOM_CAUSAL_CALLBACKS;

typedef struct {
    ST_FOM_WINDOW_LAYOUT layout;
    int num_windows;
    double dt;
    int output_interval;
    int write_output;
} ST_FOM_CAUSAL_OPTIONS;

int ST_fom_run_causal_dg0_windows(
    const ST_FOM_CAUSAL_OPTIONS* options,
    const ST_FOM_CAUSAL_CALLBACKS* callbacks);

/* Binary format intended for later ST-DDROM snapshot import. */
typedef struct {
    uint32_t magic;
    uint32_t version;
    uint32_t window_id;
    uint32_t degree;
    uint32_t num_slabs;
    uint32_t num_time_dofs;
    uint32_t physical_dof;
    uint32_t num_local_space_points;
    uint32_t num_owned_space_points;
    double time_start;
    double time_end;
} ST_FOM_WINDOW_FILE_HEADER;

int ST_fom_make_directory_recursive(const char* directory);

int ST_fom_write_owned_window_binary(
    const char* filename,
    int window_id,
    double time_start,
    double time_end,
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state);

/* Local sum-of-squares helper; MPI reduction remains the caller's choice. */
int ST_fom_owned_difference_squares(
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* reference,
    const double* approximation,
    double* difference_norm2,
    double* reference_norm2,
    double* approximation_norm2);

#ifdef __cplusplus
}
#endif
