#ifndef ROM_STD_STROM_MONOLIS_H
#define ROM_STD_STROM_MONOLIS_H

#include "core_FOM_ST.h"
#include "rom_std_strom.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * MONOLIS/BCSR representation of the reduced all-at-once system.
 *
 * One MONOLIS graph node represents one time window.  The block component
 * index represents the POD-mode index inside that window.  Therefore this
 * representation requires a common (uniform) number of POD modes in every
 * window.
 */
typedef struct {
    MONOLIS monolis;
    MONOLIS_COM* monolis_com;

    int num_windows;
    int block_size;
    int total_reduced_dof;

    /* Window-level CSR graph retained for the complete solver lifetime. */
    int* graph_index;
    int* graph_item;
    int num_graph_items;

    double* solution;

    int max_iter;
    double tolerance;
} ROM_STROM_MONOLIS_SOLVER;

/* Check that all windows have the same rank and BCSR-compatible offsets. */
int ROM_std_strom_check_uniform_rank(
    const ROM_STROM_SYSTEM* system,
    int* block_size);

/*
 * Build the lower block-bidiagonal window graph:
 *
 *   row 0 : {0}
 *   row w : {w-1, w}, w > 0.
 *
 * The self entry is explicit because it creates the dense D_w block across
 * all POD components of the same graph node.
 */
int ROM_std_strom_build_bcsr_graph(
    int num_windows,
    int** graph_index,
    int** graph_item,
    int* num_graph_items);

/* Optional human-readable dump of the retained metagraph. */
int ROM_std_strom_write_bcsr_graph(
    const ROM_STROM_MONOLIS_SOLVER* solver,
    const char* filename);

/* Allocate the MONOLIS BCSR pattern and solver work arrays. */
int ROM_std_strom_monolis_initialize(
    ROM_STROM_MONOLIS_SOLVER* solver,
    const ROM_STROM_SYSTEM* system,
    MONOLIS_COM* monolis_com,
    int max_iter,
    double tolerance);

/* Assemble D_w, -L_w and g_w into the retained MONOLIS graph. */
int ROM_std_strom_monolis_assemble(
    ROM_STROM_MONOLIS_SOLVER* solver,
    const ROM_STROM_SYSTEM* system);

/* Solve the nonsymmetric reduced all-at-once system by MONOLIS BiCGSTAB. */
int ROM_std_strom_solve_all_at_once_monolis(
    ROM_STROM_SYSTEM* system,
    ROM_STROM_MONOLIS_SOLVER* solver);

/* Relative algebraic residual of the current reduced coefficients. */
double ROM_std_strom_reduced_residual_norm(
    const ROM_STROM_SYSTEM* system);

void ROM_std_strom_monolis_finalize(
    ROM_STROM_MONOLIS_SOLVER* solver);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STROM_MONOLIS_H */
