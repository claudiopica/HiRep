/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "io.h"
#include "random.h"

// TODO: put gpu as last suffix
// TODO: fill buffers needs gpu suffix
// TODO: MPI error management
    
#define random_double ranlxd
#define random_float ranlxs

static void *gf_sendrecv_guard=NULL;

static inline void zeroes_double(double* dbl, int n) {
    for (int i = 0; i < n; ++i) { dbl[i] = 0.0; }
}

static inline void zeroes_float(float* flt, int n) {
    for (int i = 0; i < n; ++i) { flt[i] = 0.0f; }
}

/* Spinor-like fields */
#define _GEOM_TYPE spinor
#define _IS_SPINOR_LIKE 1

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE

/* Gauge fields */
#define _GEOM_TYPE gauge
#define _IS_SPINOR_LIKE 0

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/communications.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE