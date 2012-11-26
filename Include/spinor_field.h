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

#define _DECLARE_FIELD_STRUCT(_name,_type) \
typedef struct _##_name { \
_type *ptr; \
geometry_descriptor *type;\
_MPI_FIELD_DATA \
_GPU_FIELD_DATA(_type) \
} _name


_DECLARE_FIELD_STRUCT(suNg_field, suNg);
_DECLARE_FIELD_STRUCT(suNg_field_flt, suNg_flt);
_DECLARE_FIELD_STRUCT(suNf_field, suNf);
_DECLARE_FIELD_STRUCT(suNf_field_flt, suNf_flt);
_DECLARE_FIELD_STRUCT(spinor_field, suNf_spinor);
_DECLARE_FIELD_STRUCT(spinor_field_flt, suNf_spinor_flt);
_DECLARE_FIELD_STRUCT(suNg_av_field, suNg_algebra_vector);
_DECLARE_FIELD_STRUCT(scalar_field, double);



/* LOOPING MACRO */

#ifndef _PIECE_INDEX
#define _PIECE_INDEX(i) i##_pindex
#endif
#define _SPINOR_PTR(s) s##_ptr

#ifndef _DECLARE_INT_ITERATOR
#define _DECLARE_INT_ITERATOR(i) int i, _PIECE_INDEX(i)
#endif
#define _DECLARE_SPINOR_ITERATOR(s) suNf_spinor* _SPINOR_PTR(s)
#define _DECLARE_SPINOR_FLT_ITERATOR(s) suNf_spinor_flt* _SPINOR_PTR(s)


#ifdef CHECK_SPINOR_MATCHING

#define _TWO_SPINORS_MATCHING(s1,s2) \
	error((s1)->type!=(s2)->type,1,__FILE__ ": ", "Spinors don't match!");

#define _ARRAY_SPINOR_MATCHING(s,i,n) \
	for(i=0; i<n; i++) \
		error((s)->type!=((s)+i)->type,1,__FILE__ ": ", "Spinors don't match!");

#else /* CHECK_SPINOR_MATCHING */

#define _TWO_SPINORS_MATCHING(s1,s2)

#define _ARRAY_SPINOR_MATCHING(s,i,n)

#endif /* CHECK_SPINOR_MATCHING */



#define _ONE_SPINOR_FOR(s,i) \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(s)->type->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(s)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s)=(s)->ptr+(s)->type->master_start[_PIECE_INDEX(i)]-(s)->type->master_shift; \
	    i<=(s)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s)++ \
	   )

#define _TWO_SPINORS_FOR(s1,s2,i) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(s1)->type->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(s1)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s1)=(s1)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]-(s1)->type->master_shift, \
	              _SPINOR_PTR(s2)=(s2)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]-(s1)->type->master_shift; \
	    i<=(s1)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++ \
	   )

#define _THREE_SPINORS_FOR(s1,s2,s3,i) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	_TWO_SPINORS_MATCHING(s1,s3); \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(s1)->type->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(s1)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s1)=(s1)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]-(s1)->type->master_shift, \
	              _SPINOR_PTR(s2)=(s2)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]-(s1)->type->master_shift, \
	              _SPINOR_PTR(s3)=(s3)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]-(s1)->type->master_shift; \
	    i<=(s1)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++, _SPINOR_PTR(s3)++ \
	   )

#include "field_ordering.h"

#define _FIELD_AT(s,i) (((s)->ptr)+i-(s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr)+coord_to_index(i-(s)->type->master_shift,mu))


#endif
