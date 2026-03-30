#pragma once

#include "monolis.h"

#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include "BBFE/std/integ.h"
#include "BBFE/std/shapefunc.h"
#include "BBFE/std/mapping.h"
#include "BBFE/std/surface.h"

#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "BBFE/sys/monowrap.h"

#include "BBFE/elemmat/set.h"
#include "BBFE/elemmat/equivval.h"
#include "BBFE/elemmat/convdiff.h"

#include "BBFE/manusol/manusol.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

static const int BLOCK_SIZE                   = 1;
static const char* CODENAME                   = "convdiff >";

static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";
static const char* INPUT_FILENAME_D_BC        = "D_bc.dat";


const char* BBFE_convdiff_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename);

void BBFE_convdiff_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution);

void BBFE_convdiff_set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis);

double BBFE_convdiff_equivval_relative_L2_error_scalar(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		MONOLIS_COM*  monolis_com,
		double        t,
		const double* comp_vec, // [total_num_nodes]
		double        (*func)(double, double, double, double)); // scalar function(x, y, z, t)

void BBFE_convdiff_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		BBFE_BC*     bc);
