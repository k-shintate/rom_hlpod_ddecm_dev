
#pragma once

#include "convdiff_core.h"

#include "monolis_utils.h"
#include "monolis_wrapper_scalapack_c.h"

void ROM_BB_vec_copy_1d(
    double*  in,   //input
    double*  out,  //output
    const int num1,
    const int num2);

void ROM_BB_vec_copy_2d(
    double**  in,   //input
    double**  out,  //output
    const int num1,
    const int num2);

double calc_internal_norm_1d(
    double*  in,       //input
    const int num1,
    const int num2);

double _Complex BBFE_std_integ_calc_C(
    const int num_integ_points,
    const double _Complex *value,
    const double *weight,
    const double *Jacobian);

double calc_nodal_error_vector(
    double** edge_result,
    double** source,
    double** edge_error,
    int N);

void copy_Aphi_to_V_phi(
    double _Complex * Aphi,
    double * V,
    double * phi,
    const int total_num_nodes,
    const int total_num_edges);

void copy_Aphi_to_V_phi_time(
    double * Aphi,
    double * V,
    double * phi,
    const int total_num_nodes,
    const int total_num_edges);

void compute_Js(
    double x,
    double y,
    double z,
    double nu,
    double sigma,
    double omega,
    double _Complex Js[3]);

void compute_Js_time(
    double x,
    double y,
    double z,
    double t,
    double lambda,
    double Nu,
    double Sigma,
    double out[3]);

double**** BB_std_calloc_4d_double(
    double****  array,
    const int  size1,
    const int  size2,
    const int  size3,
    const int  size4);

int cholesky_decomposition(
    double** A,
    double** L,
    int n);

void calc_eigenmodes(
    double** matrix_R_L,
    double** matrix_R_A,
    int size);
