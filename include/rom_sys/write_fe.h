#pragma once

#include "BBFE/sys/write.h"
#include "BB/vtk.h"
#include "BB/std.h"

#include "monolis.h"

#include "rom_dataset.h"

void ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
    const int 		n_internal_vertex,
    BBFE_DATA*		fe,
    double**     	pod_modes,
    const int 		num_modes,
    const char*		output_fname,
    const char*		directory);

void ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
    const int 		n_internal_vertex,
    BBFE_DATA*		fe,
    double**     	pod_modes,
    const int 		num_modes,
    const char*		output_fname,
    const char*		directory);

void ROM_sys_hlpod_fe_write_pod_modes(
    BBFE_DATA*		fe,
    double**		pod_modes,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int 		ndof,
    const char*     label,
    const char*		directory);

void ROM_sys_hlpod_fe_write_local_pod_modes(
    BBFE_DATA*		fe,
    double**		pod_modes,
    const int 		n_internal_vertex,
    const int 		num_modes,
    const int		ndof,
    const char*     label,
    const char*		directory);

void ROM_sys_hlpod_fe_write_pod_modes_diag(
    BBFE_DATA*		fe,
    double**		pod_modes,
    const int 		n_internal_vertex,
    const int 		num_modes_v,
    const int 		num_modes_vi,
    const char*     label_v,
    const char*		label_p,
    const char*		directory);

void ROM_sys_hlpod_fe_write_local_pod_modes_diag(
    BBFE_DATA*		fe,
    double**		pod_modes,
    const int 		n_internal_vertex,
    const int 		num_modes_v,
    const int 		num_modes_vi,
    const char*     label_v,
    const char*     label_p,
    const char*		directory);

void ROM_sys_hlpod_fe_write_local_pod_modes_diag_id(
    BBFE_DATA*		fe,
    double**		pod_modes,
    int*			node_id_local,
    const int 		n_internal_vertex,
    const int 		num_modes_v,
    const int 		num_modes_vi,
    const char*     label_v,
    const char*     label_p,
    const char*		directory);
