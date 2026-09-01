#ifndef ST_FOM_PROJECT_ADAPTER_H
#define ST_FOM_PROJECT_ADAPTER_H

#include "core_FOM.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Thin adapter between the existing diffusion project and the new normal
 * space-time solver.  It intentionally keeps project-specific mesh/BC setup
 * outside core_FOM_ST_spaceblock.c.
 */
typedef struct {
    FE_SYSTEM fom;
    char directory[4096];

    double finish_time;
    int output_interval;
    int num_ip_each_axis;
    int time_degree;
    int slabs_per_window;

    int base_monolis_initialized;
    int bc_initialized;
    int values_initialized;
    int fe_initialized;
} ST_FOM_PROJECT;

int ST_fom_project_initialize(
    ST_FOM_PROJECT* project,
    int argc,
    char* argv[]);

void ST_fom_project_finalize(ST_FOM_PROJECT* project);

int ST_fom_project_num_windows(const ST_FOM_PROJECT* project);

#ifdef __cplusplus
}
#endif

#endif /* ST_FOM_PROJECT_ADAPTER_H */
