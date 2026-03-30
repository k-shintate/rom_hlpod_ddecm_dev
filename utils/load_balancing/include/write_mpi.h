#pragma once

#include "lb_dataset.h"

void LB_write_node(
    const int ie,
    LPOD_LB*        lpod_lb,
    int* num_search,
    double** global_x_internal,
    double** global_x_overlap,
    int* overlap_id_for_sort,
    const char* directory);

void LB_write_node_internal(
    const int ndof, 
    const int total_n_internal,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_node_id(
    const int ndof, 
    const int ie,
    int* overlap_id,
    int* num_search,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_node_recv(
    int* recv_index,
    int* num_search,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_node_send(
    int* recv_index,
    int* num_send,
    int** send_id,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_D_bc(
    int* global_D_bc_index,
    const int return_val,
    double* global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_D_bc_p(
    int* global_D_bc_index,
    const int return_val,
    double* global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const char*     label,
    const char* directory);

void LB_write_D_bc_v(
    int*            global_D_bc_index,
    const int       return_val,
    double**         global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const int       dof,
    const char*     directory);

void LB_write_elem(
    const int local_node_dof,
    const int num_internal,
    const int num_total,
    const int total_n_overlap,
    int* elem_id_bef_sort_overlap,
    int* id_overlap,
    int* sum_internal_global_id,
    int** global_conn,
    int* sum_internal_conn_id,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_elem_internal(
    const int ndof, 
    const int num,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_elem_id(
    const int ndof, 
    const int total_num_elem,
    const int total_num_elem_internal,
    const int total_num_elem_overlap,
    int* sum_internal_global_id,
    int* elem_id_bef_sort_overlap,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_elem_id_2(
    const int ndof, 
    const int total_num_elem,
    const int total_num_elem_internal,
    int* sum_internal_global_id,
    LPOD_LB*        lpod_lb,
    const char* directory);

void LB_write_elem_2(
    const int local_node_dof,
    const int total_num_elem,
    const int num_internal,
    int* sum_internal_global_id,
    int** global_conn,
    int* sum_internal_conn_id,
    LPOD_LB*        lpod_lb,
    const char* directory);