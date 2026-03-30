#pragma once

#include "read_BB.h"

int ROM_std_hlpod_read_n_internal(
    FILE*        fp,
    const char*  fname,
    const char*  directory);

void ROM_std_hlpod_read_node_id(
    FILE*        fp,
    int*         node_id,
    const int    n_internal_vertex,
    const char*  fname,
    const char*  directory);

void ROM_std_hlpod_read_singular_values(
    double*                     singular_value,
    const int           num_modes,
    const int           subdomain_id,
    const char*         label,
    const char*         directory);

int ROM_std_hlpod_read_num_modes(
    const int		subdomain_id,
    const char*		label,
    const char*		directory);

void ROM_std_hlpod_read_pod_modes_node(
    FILE*        fp,
    double**     S,
    const int    dof,
    const int    n_internal_vertex_in,
    const int    num_modes,
    const int    subdomain_id,
    const char*  label,
    const char*  directory);
