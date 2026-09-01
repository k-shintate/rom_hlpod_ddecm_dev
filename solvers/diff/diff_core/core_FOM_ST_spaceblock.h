#ifndef CORE_FOM_ST_SPACEBLOCK_H
#define CORE_FOM_ST_SPACEBLOCK_H

#include <stdio.h>

#include "core_FOM.h"
#include "st_fom_common.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    ST_FOM_SOLVE_CAUSAL = 0,
    ST_FOM_SOLVE_WINDOW = 1
};

typedef struct {
    FE_SYSTEM* fom;

    ST_FOM_WINDOW_LAYOUT layout;
    ST_FOM_SPATIAL_GRAPH spatial_graph;

    /* Optional monolithic all-at-once window matrix. */
    MONOLIS monolis_window;
    BBFE_BC bc_window;

    /* Default robust path: one dG slab block at a time. */
    MONOLIS monolis_slab;
    BBFE_BC bc_slab;

    int num_owned_space_nodes;
    int window_dof;
    int slab_dof;
    int solve_mode;

    int window_matrix_initialized;
    int window_bc_initialized;
    int slab_matrix_initialized;
    int slab_bc_initialized;
    int graph_initialized;

    double* window_solution;
    double* slab_solution;
    double* previous_trace;
    double* trace_work;

    FILE* accuracy_fp;
    char accuracy_filename[4096];
    char graph_filename[4096];

    /* Optional ParaView/VTK output of each slab right trace. */
    int write_vtk;
    int vtk_interval;
    FILE* vtk_pvd_fp;
    char vtk_pvd_filename[4096];
} ST_DIFFUSION_SPACEBLOCK;

/*
 * Initialize a normal (non-ROM) space-time diffusion solver on an already
 * initialized spatial FE_SYSTEM.  The ordinary spatial MONOLIS communicator
 * in fom->monolis_com is reused; all temporal DOFs remain block DOFs on each
 * spatial node.
 *
 * Runtime selection:
 *   ST_FOM_SOLVE_MODE=causal  (default)
 *       Exact block forward substitution in slab order.  Each linear solve
 *       contains all time-basis coefficients of one dG slab (degree+1 DOFs).
 *       For the lower-triangular linear dG window system this is algebraically
 *       equivalent to solving the complete all-at-once window.
 *
 *   ST_FOM_SOLVE_MODE=window
 *       Solve the complete all-at-once ST window in one MONOLIS call.
 */
int ST_diffusion_spaceblock_initialize(
    ST_DIFFUSION_SPACEBLOCK* solver,
    FE_SYSTEM* fom,
    int time_degree,
    int slabs_per_window,
    const char* accuracy_filename);

void ST_diffusion_spaceblock_finalize(
    ST_DIFFUSION_SPACEBLOCK* solver);

int ST_diffusion_spaceblock_solve_window(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int window_id,
    double time_window_start);

/*
 * Solve num_windows consecutive windows.  Accuracy is written at every slab
 * right endpoint using manusol_get_sol().
 *
 * VTK output is controlled at runtime by:
 *   ST_FOM_WRITE_VTK=1          enable VTK output
 *   ST_FOM_VTK_INTERVAL=N       write every N time steps
 *   ST_FOM_VTK_PVD=name.pvd     collection filename (default st_results.pvd)
 *
 * Each MPI rank writes one legacy VTK piece and rank 0 writes a PVD collection
 * so that ParaView can load all pieces and all output times from one file.
 */
int ST_diffusion_spaceblock_run(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int num_windows);

#ifdef __cplusplus
}
#endif

#endif /* CORE_FOM_ST_SPACEBLOCK_H */
