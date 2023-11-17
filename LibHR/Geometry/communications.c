/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "io.h"
#include "random.h"

#define random_double ranlxd
#define random_float ranlxs

void probe_mpi(void) {
#ifdef WITH_MPI
    int flag;
    MPI_Status status[1];
    // MPI_Testall(nreq, f->comm_req, &flag, status);
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, cart_comm, &flag, status);
#endif
}

/* Spinor-like fields */
#define _GEOM_TYPE spinor
#define _IS_SPINOR_LIKE 1

#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE

/* Gauge fields */
#define _GEOM_TYPE gauge
#define _IS_SPINOR_LIKE 0

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _MPI_REAL MPI_FLOAT
#define _REAL float
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _MPI_REAL MPI_DOUBLE
#define _REAL double
#include "TMPL/lock_communications.c.tmpl"
#include "TMPL/communications.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE
