/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef SPINOR_FIELD_H
#define SPINOR_FIELD_H

#include "geometry.h"
#include "suN_types.h"
#include "error.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* MPI data */
#define _MPI_FIELD_DATA
#ifdef WITH_MPI
#undef _MPI_FIELD_DATA
#define _MPI_FIELD_DATA MPI_Request *comm_req;
#endif //WITH_MPI

/* GPU data */
#define _GPU_FIELD_DATA(_type)
#ifdef WITH_GPU
#undef _GPU_FIELD_DATA
#define _GPU_FIELD_DATA(_type) _type *gpu_ptr;
#endif //WITH_MPI

typedef struct {
	complex up[NF*(2*NF+1)];
	complex dn[NF*(2*NF+1)];
} ldl_t;

#define _DECLARE_FIELD_STRUCT(_name,_type) \
typedef struct _##_name { \
_type *ptr; \
geometry_descriptor *type;\
_MPI_FIELD_DATA \
_GPU_FIELD_DATA(_type) \
} _name


_DECLARE_FIELD_STRUCT(suNg_field, suNg);
_DECLARE_FIELD_STRUCT(suNg_scalar_field, suNg_vector);
_DECLARE_FIELD_STRUCT(suNg_field_flt, suNg_flt);
_DECLARE_FIELD_STRUCT(suNf_field, suNf);
_DECLARE_FIELD_STRUCT(suNf_field_flt, suNf_flt);
_DECLARE_FIELD_STRUCT(spinor_field, suNf_spinor);
_DECLARE_FIELD_STRUCT(spinor_field_flt, suNf_spinor_flt);
_DECLARE_FIELD_STRUCT(suNg_av_field, suNg_algebra_vector);
_DECLARE_FIELD_STRUCT(scalar_field, double);
_DECLARE_FIELD_STRUCT(ldl_field, ldl_t);
_DECLARE_FIELD_STRUCT(suNfc_field, suNfc);


/* LOOPING MACRO */

#ifdef CHECK_SPINOR_MATCHING

#define _TWO_SPINORS_MATCHING(s1,s2) \
	error((s1)->type!=(s2)->type,1,__FILE__ ": ", "Spinors don't match!");

#define _ARRAY_SPINOR_MATCHING(s,n) \
	for(int i=0; i<n; i++) \
		error((s)->type!=((s)+i)->type,1,__FILE__ ": ", "Spinors don't match!");

#else /* CHECK_SPINOR_MATCHING */

#define _TWO_SPINORS_MATCHING(s1,s2)

#define _ARRAY_SPINOR_MATCHING(s,n)

#endif /* CHECK_SPINOR_MATCHING */


#define _ONE_SPINOR_FOR_RED(s,redop1,redop2) _MASTER_FOR_RED((s)->type,_spinor_for_is,redop1,redop2)
#define _ONE_SPINOR_FOR(s) _ONE_SPINOR_FOR_RED(s,,)
#define _ONE_SPINOR_FOR_SUM(s,...) _ONE_SPINOR_FOR_RED(s,_omp_sum(__VA_ARGS__),)

#define _TWO_SPINORS_FOR_RED(s1,s2,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)
#define _TWO_SPINORS_FOR(s1,s2) _TWO_SPINORS_FOR_RED(s1,s2,,)
#define _TWO_SPINORS_FOR_SUM(s1,s2,...) _TWO_SPINORS_FOR_RED(s1,s2,_omp_sum(__VA_ARGS__),)

#define _THREE_SPINORS_FOR_RED(s1,s2,s3,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _TWO_SPINORS_MATCHING(s1,s3); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)
#define _THREE_SPINORS_FOR(s1,s2,s3) _THREE_SPINORS_FOR_RED(s1,s2,s3,,)

#include "field_ordering.h"

#define _FIELD_AT(s,i) (((s)->ptr)+i-(s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr)+coord_to_index(i-(s)->type->master_shift,mu))
#define _6FIELD_AT(s,i,mu) (((s)->ptr)+((i-(s)->type->master_shift)*6+mu))

#define _SPINOR_PTR(s) _FIELD_AT(s,_spinor_for_is)


#endif
