#pragma once

#include "lb_dataset.h"

int read_total_num_D_bc(
    const int       Ddof,
    int*            mydomain,
    const char*     directory);

int read_total_num_D_bc_p(
    const int       Ddof,
    int*            mydomain,
    const char*     label,
    const char*     directory);

int read_total_num_D_bc_v(
    const int       Ddof,
    int*            mydomain,
    const char*     directory);

