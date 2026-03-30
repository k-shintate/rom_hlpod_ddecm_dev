
#pragma once

#include "BBFE/sys/write.h"
#include "BB/vtk.h"
#include "BB/std.h"

#include "monolis.h"

#include "write_BB.h"
#include "rom_dataset.h"

void output_hr_monolis_solver_prm(
    MONOLIS*  		monolis,
	const char*  	directory,
	double 			t);
