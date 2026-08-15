#include "fluid_sups_st_model.h"
#include "fluid_sups_nr_api_compat.h"
#include "fluid_sups_jold.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fluid_sups_st_model/system.inc"
#include "fluid_sups_st_model/element_nr.inc"
#include "fluid_sups_st_model/window_dg0.inc"
#include "fluid_sups_st_model/window_dg_high_order.inc"
#include "fluid_sups_st_model/rom_adapter.inc"
#include "fluid_sups_st_model/step_lifecycle.inc"
