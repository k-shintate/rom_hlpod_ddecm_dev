#pragma once

#include "rom_std_stddrom.h"
#include "set_modes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Import the rank-local, owned-row HLPOD basis into the generic distributed
 * ST-DDROM system.
 *
 * Required row ordering:
 *
 *   [owned spatial node][slab][time basis][ux,uy,uz,p]
 *
 * i.e. exactly ST_fom_vector_index() restricted to the owned-node prefix.
 */
int ROM_std_stdd_import_hlpod_local_basis(
    ROM_STDD_SYSTEM* system,
    const HLPOD_MAT* hlpod_mat,
    int num_modes);

/* Copy one window's local reduced coefficients into a legacy HLPOD array. */
int ROM_std_stdd_export_window_coefficients(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* mode_coef,
    int mode_coef_capacity);

/*
 * Optional compatibility helper for existing HLPOD output/visualization code.
 * hlpod_mat->pod_modes must already be allocated as
 * [internal_st_dof][num_modes_capacity].
 */
int ROM_std_stdd_export_local_basis_to_hlpod(
    const ROM_STDD_SYSTEM* system,
    HLPOD_MAT* hlpod_mat,
    int pod_mode_capacity);

#ifdef __cplusplus
}
#endif
