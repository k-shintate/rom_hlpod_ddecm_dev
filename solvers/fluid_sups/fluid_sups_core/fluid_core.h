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
#include "BBFE/elemmat/fluid.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

static const char* CODENAME                   = "fluid >";

static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";


const char* BBFE_fluid_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename);

void BBFE_fluid_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis);

void BBFE_fluid_set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis);

void BBFE_fluid_sups_renew_velocity(
		double**  v,
		double*   ans_vec,
		const int total_num_nodes);

void BBFE_fluid_sups_renew_pressure(
		double*   p,
		double*   ans_vec,
		const int total_num_nodes);

void BBFE_fluid_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis);
