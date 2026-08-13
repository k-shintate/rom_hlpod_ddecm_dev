#include "rom_std_stddrom_hlpod_bridge.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

static int bridge_validate_common_mode_count(
    int local_modes)
{
    int maximum = local_modes;
    int negative_minimum = -local_modes;

    monolis_allreduce_I(
        1,
        &maximum,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    monolis_allreduce_I(
        1,
        &negative_minimum,
        MONOLIS_MPI_MAX,
        monolis_mpi_get_global_comm());

    return maximum == -negative_minimum
        ? ROM_STDD_SUCCESS
        : ROM_STDD_ERR_DIMENSION;
}

int ROM_std_stdd_import_hlpod_local_basis(
    ROM_STDD_SYSTEM* system,
    const HLPOD_MAT* hlpod_mat,
    int num_modes)
{
    if(system == NULL || hlpod_mat == NULL ||
       hlpod_mat->pod_modes == NULL ||
       system->internal_st_dof <= 0 ||
       system->num_modes_capacity <= 0 ||
       num_modes <= 0 ||
       num_modes > system->num_modes_capacity)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    if(bridge_validate_common_mode_count(num_modes)
       != ROM_STDD_SUCCESS)
    {
        /*
         * The current generic distributed reduced solver uses a fixed block
         * size, so every spatial rank must expose the same number of modes.
         */
        return ROM_STDD_ERR_DIMENSION;
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
                (size_t)row
                    * (size_t)system->num_modes_capacity
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
       system->reduced_coef == NULL ||
       window_id < 0 ||
       window_id >= system->num_windows ||
       system->num_modes <= 0 ||
       mode_coef_capacity < system->num_modes)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    memcpy(
        mode_coef,
        system->reduced_coef
            + (size_t)window_id
                * (size_t)system->num_modes_capacity,
        (size_t)system->num_modes * sizeof(double));

    return ROM_STDD_SUCCESS;
}

int ROM_std_stdd_export_local_basis_to_hlpod(
    const ROM_STDD_SYSTEM* system,
    HLPOD_MAT* hlpod_mat,
    int pod_mode_capacity)
{
    if(system == NULL || hlpod_mat == NULL ||
       system->basis == NULL || hlpod_mat->pod_modes == NULL ||
       system->num_modes <= 0 ||
       pod_mode_capacity < system->num_modes)
    {
        return ROM_STDD_ERR_ARGUMENT;
    }

    for(int row = 0; row < system->internal_st_dof; row++){
        for(int mode = 0; mode < system->num_modes; mode++){
            hlpod_mat->pod_modes[row][mode] =
                system->basis[
                    (size_t)row
                        * (size_t)system->num_modes_capacity
                    + (size_t)mode];
        }
    }

    return ROM_STDD_SUCCESS;
}
