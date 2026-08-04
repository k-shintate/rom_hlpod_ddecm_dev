#include "rom_std_stddrom_hlpod_bridge.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

int ROM_std_stdd_import_hlpod_local_basis(
    ROM_STDD_SYSTEM* system,
    const HLPOD_MAT* hlpod_mat,
    int num_modes)
{
    if(system == NULL || hlpod_mat == NULL || hlpod_mat->pod_modes == NULL ||
       num_modes <= 0 || num_modes > system->num_modes_capacity)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    free(system->basis);
    system->basis = (double*)calloc(
        (size_t)system->internal_st_dof
            * (size_t)system->num_modes_capacity,
        sizeof(double));
    if(system->basis == NULL){
        return ROM_STDD_ERR_ALLOCATION;
    }

    system->num_modes = num_modes;

    for(int row = 0; row < system->internal_st_dof; row++){
        for(int mode = 0; mode < num_modes; mode++){
            system->basis[
                (size_t)row * (size_t)system->num_modes_capacity
                + (size_t)mode] =
                hlpod_mat->pod_modes[row][mode];
        }
    }

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_export_window_coefficients(
    const ROM_STDD_SYSTEM* system,
    int window_id,
    double* mode_coef,
    int mode_coef_capacity)
{
    if(system == NULL || mode_coef == NULL ||
       window_id < 0 || window_id >= system->num_windows ||
       mode_coef_capacity < system->num_modes)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memcpy(
        mode_coef,
        system->reduced_coef
            + (size_t)window_id * (size_t)system->num_modes_capacity,
        (size_t)system->num_modes * sizeof(double));

    return ROM_STDD_SUCCESS;
}
