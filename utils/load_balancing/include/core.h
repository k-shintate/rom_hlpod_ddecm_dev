#pragma once

#include "lb_dataset.h"
#include "std.h"

void lb_read_meta_graph(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);

void lb_set_node_global_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);

void lb_write_node_local_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);

void lb_write_sendrecv(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);

void lb_set_elem_global_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);


void lb_write_elem_local_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb);