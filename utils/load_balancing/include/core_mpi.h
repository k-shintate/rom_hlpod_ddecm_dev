#pragma once

#include "read_mpi.h"
#include "write_mpi.h"
#include "graph.h"
#include "lb_dataset.h"

#include "std.h"

void lb_set_node_internal(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_neib_list(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

//recv,idの作成
void lb_set_node_overlap(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_node_overlap_without_write(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_node_send(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_D_bc(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_D_bc_p(
    const int       nd1,
    const char*     directory,
    const char*     label,
    LPOD_LB*        lpod_lb);

void lb_set_D_bc_v(
    const int       nd1,
    const char*     directory,
    const int       dof,
    const char*     label,
    LPOD_LB*        lpod_lb);

void lb_set_elem_internal(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);

void lb_set_elem_internal_2(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);