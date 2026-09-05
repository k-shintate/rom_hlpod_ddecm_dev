#ifndef HROM_STAGE3_PARALLEL_MODE0_H
#define HROM_STAGE3_PARALLEL_MODE0_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Parallel Stage-3 mode-0 NNLS.
 *
 * Reads DDECM/stage3_rank.<rank>.bin (H3NNLS4 / version 4),
 * constructs the global row space using the stored global row keys,
 * forms the union of Stage-2 selected original-global element IDs,
 * distributes those candidate columns contiguously across `comm`, and
 * solves the sparse NNLS with MONOLIS' column-distributed parallel solver.
 *
 * input_rank_count:
 *   > 0 : read rank files [0, input_rank_count)
 *   <=0 : rank 0 auto-detects a contiguous sequence from rank 0.
 *
 * tol_inner:
 *   > 0 : use as MONOLIS inner/KKT tolerance.
 *   <=0 : use 1e-14 * max(global_rows, global_candidate_columns)
 *          after common Frobenius normalization of A and b.
 *
 * Legacy Stage-3 outputs (authoritative; same responsibility as serial Stage 3):
 *   DDECM/lb_selected_elem.<fine-subdomain>.txt
 *   DDECM/lb_selected_elem_D_bc.<fine-subdomain>.txt
 *
 * DDECM/selected_elem.<rank>.txt and selected_elem_D_bc.<rank>.txt are NOT
 * generated here.  They belong to the existing parallel/online conversion path.
 *
 * output_relative_path:
 *   Optional diagnostic file containing original-global IDs.  Pass NULL or an
 *   empty string to disable it.  Online HROM does NOT depend on this file.
 */
int HROM_stage3_parallel_mode0_run(
    MPI_Comm comm,
    const char* directory,
    int input_rank_count,
    int max_iter,
    double tol_outer,
    double tol_inner,
    double weight_tol,
    const char* output_relative_path);

#ifdef __cplusplus
}
#endif

#endif
