#pragma once

#include "monolis.h"

#include "BB/std.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>

typedef struct
{
    /*User-specified value*/
    int num_modes_pre;		    //Maximum number of basis functions per POD computaion subdomains
    double rom_epsilon;		    //mode acquisition threshold
    int num_1st_subdomains;     //POD computaion subdomains
    /**/

    double* sol_vec;		    //solution vector of ROM
    int num_2nd_subdomains;     //parallel computation subdomains
    int num_modes;              //number of modes
    int num_modes_max;		    //maximum number of modes per parallel computation subdomains
    int num_snapshot;           //number of snapshots
    int n_neib_vec;			    //number of modes including neighboring subdomains in the parallel computation domain

    bool bool_global_mode;

    double time_calc_reduced_matvec;
    double time_linear_solver;
    

} HLPOD_VALUES;


//informations about matrix
typedef struct
{
    double** KV;                    //reduced matrix (KV)
    double** VTKV;                  //reduced matrix (VTKV)
    double* VTf;                    //reduced rhs vector
    double* VTf_tmp;                    //reduced rhs vector
    double* VTf_D_bc;                    //reduced rhs vector
    double* VTf_source;                    //reduced rhs vector

    double* mode_coef;              //coefficients for pod modes
    double* mode_coef_old;              //coefficients for pod modes
    double** pod_modes;             //pod modes
    double** snapmat;               //snapshot matrix
                                            
    double** neib_vec;               //pod_modes updated based on 2nd subdomains comm

    double** pod_modes_decoupled_v;             //pod modes
    double** neib_vec_decoupled_v;               //pod_modes updated based on 2nd subdomains comm
    double** pod_modes_decoupled_p;             //pod modes
    double** neib_vec_decoupled_p;               //pod_modes updated based on 2nd subdomains comm
    double** pod_basis_hr_decoupled_v;
    double** pod_basis_hr_decoupled_p;

    int* n_internal_vertex_subd;     //n_internal_vertex for 1st subdomains
    int* node_id;                    //node id for 1st subdomains
    int* num_modes_internal;         //number of pod modes in internal subdomains

    int* num_modes_2nddd;            //number of pod modes in 2nd subdomains (own)
    int* num_modes_2nddd_max;        //number of pod modes in 2nd subdomains (own, max)

    int* num_modes_1stdd;            //number of pod modes in 1st subdomains (own, index format)
    int* num_modes_1stdd_neib;       //number of pod modes in 1st subdomains (own + neib)
    

    //追加
    int* num_neib_modes_sum;
    double* mode_coef_1stdd;
    int num_metagraph_nodes;
    int* subdomain_id_in_nodes;
    int* subdomain_id_in_nodes_2nddd;
    double** pod_basis_hr;
    int* max_num_neib_modes;
    int** subdomain_id_in_nodes_internal;

    int* hr_D_bc_node_id;
    int num_hr_D_bc_nodes;

} HLPOD_MAT;

//informations about metagraph
typedef struct
{
    int num_meta_nodes;              //number of nodes in neib 1st subdomain

    int* subdomain_id;              //id of 1st subdomains (internal)
    int* subdomain_id_neib;         //id of 1st subdomains (after updated based on 2nd subdomains comm)

    int n_internal_sum;             //sum of n_internal_vertex of 1st subdomains (neib)

    int* item;                      //item used in l-pod calculations in BCSR format
    int* index;                     //index used in l-pod calculations in BCSR format

    int* n_dof_list;                //for arbit_dof_monolis_solver

//追加
    int* n_internal;
    int* global_id;
    int* my_global_id;
    int num_neib;
    int* neib_id;                  //id of 2nd subdomains (neib)

} HLPOD_META;


typedef struct
{
    int* D_bc_node_id;               //node_id for Dirichlet b.c.

} ROM_BC;


typedef struct
{
    HLPOD_VALUES    hlpod_vals;
    HLPOD_MAT       hlpod_mat;
    ROM_BC		    rom_bc;
    HLPOD_META		hlpod_meta;

} ROM;
