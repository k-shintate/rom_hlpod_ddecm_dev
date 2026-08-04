#ifndef ROM_STD_STDDROM_HLPOD_BRIDGE_H
#define ROM_STD_STDDROM_HLPOD_BRIDGE_H

#include "rom_std_stddrom.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Import an existing HLPOD local basis into the ST-DDROM rank-local basis.
 * The HLPOD basis must use the same node-major ST block ordering as the
 * spatial-block FOM:
 *
 *   row = local_space_node * st_dof_per_space_node + q.
 *
 * Only owned rows are copied.  This bridge lets the existing
 * ROM_std_hlpod_set/read_podmodes_local_para path be retained while the
 * reduced all-at-once graph and solver are supplied by ST-DDROM.
 */
int ROM_std_stdd_import_hlpod_local_basis(
    ROM_STDD_SYSTEM* system,
    const HLPOD_MAT* hlpod_mat,
    int num_modes);

/* Copy one window's solved coefficients into a legacy HLPOD mode array. */
int ROM_std_stdd_export_window_coefficients(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* mode_coef,
    int mode_coef_capacity);

#ifdef __cplusplus
}
#endif

#endif /* ROM_STD_STDDROM_HLPOD_BRIDGE_H */
