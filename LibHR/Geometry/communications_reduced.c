/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "io.h"
#include "random.h"
#include "memory.h"

#define _GEOM_TYPE spinor
#define _IS_SPINOR_LIKE 1

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications_reduced.c.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/communications_reduced.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE
