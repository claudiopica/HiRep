/***************************************************************************\
* Copyright (c) 2022, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "Utils/background_field.h"
#include "Utils/boundary_conditions.h"
#include "Utils/clover_exp.h"
#include "Utils/data_storage.h"
#include "Utils/eva_deflation.h"
#include "Utils/gauge_anisotropy.h"
#include "Utils/gaugefix.h"
#include "Utils/HYP_smearing.h"
#include "Utils/mat_utils.h"
#include "Utils/max_eig.h"
#include "Utils/print_compile_options.h"
#include "Utils/shift_fields.h"
#include "Utils/single_double_utils.h"
#include "Utils/spatial_transformations.h"
#include "Utils/suN_mat_utils.h"
#include "Utils/test_utils.h"
#include "Utils/timing.h"
#include "Utils/wilsonflow.h"
#include "Utils/work_space.h"


// void cross_prod(suNg_vector *v1, suNg_vector *v2, suNg_vector *v3); //TOOD: this is not defined in libhr
// void cross_prod_flt(suNg_vector_flt *v1, suNg_vector_flt *v2, suNg_vector_flt *v3); //TOOD: this is not defined in libhr


#endif
