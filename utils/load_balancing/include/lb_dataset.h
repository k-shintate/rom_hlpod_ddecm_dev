#pragma once

#include "monolis.h"

#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include "BBFE/std/integ.h"
#include "BBFE/std/shapefunc.h"
#include "BBFE/std/mapping.h"

#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

static const int BUFFER_SIZE = 10000;

static const int BLOCK_SIZE                   = 1;
static const char* CODENAME                   = "convdiff >";

static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";
static const char* INPUT_FILENAME_D_BC        = "D_bc.dat";

static const char* OUTPUT_FILENAME_DOMAIN_GRAPH          = "metagraph.dat";
static const char* OUTPUT_FILENAME_NODAL_WEIGHT          = "node_weight.dat";
static const char* OUTPUT_FILENAME_PARTED_METAGRAPH      = "metagraph_parted.0";


//逐次
typedef struct
{
	int* Ddof;
	int** list_n_neib;
	int** global_id;
	int max_n_neib;

	int total_num_nodes;
	int total_num_elems;

    bool** bool_global_id_internal;
    bool** bool_global_id_overlap;

	int* sum_n_internal;

} POD_LB;

//並列
typedef struct
{
	//一部配列にする必要がない
	/*メタグラフ関連*/
	int myrank;
	int Ddof;
	int max_Ddof;

	int* domain_list;

	int* meta_n_neib;
	int meta_domain_neib;
	bool*** overlap_exists;
	int** meta_list_num_n_neib;
	int** meta_list_n_neib;

	int* n_neib_graph;
	int* meta_list_neib;

	int* meta_domain_list_neib;
	int** domain_list_neib;
	/*******************/

	int* n_internal_vertex;

	int sum_overlap;
	int total_n_internal;

	int* node_id_internal;
	int* node_id_internal_for_sort;

	int* node_id_local_overlap;
	int** recv_id;

	int* total_num_nodes_local;

	/*globalとlocalの対応関係*/
	int* local_id;
	int* local_id_for_sort;
	/************************/
	int** meta_list_num_n_neib_send;

	int total_num_nodes_lb;
	int total_num_D_bc;
	
} LPOD_LB;