#pragma once

#include "hlpod_dataset.h"
#include "fluid_fs_dataset.h"

void memory_allocation_nodal_values(
    VALUES*     vals,
    const int   total_num_nodes);

void assign_default_values(
    VALUES*     vals);

void print_all_values(
    VALUES*     vals);

void read_calc_conditions(
    VALUES*     vals,
    const char* directory);

void output_result_file_vtk(
    BBFE_DATA*     fe,
    VALUES*        vals,
    const char*    filename,
    const char*    directory,
    double         t);

void output_files(
    FE_SYSTEM*  sys,
    int         file_num,
    double      t);

void set_element_mat_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_corr(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_ppe(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void solver_FOM_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);
