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

typedef struct _suNg_field {
  suNg *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} suNg_field;

typedef struct _suNg_field_flt {
  suNg_flt *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} suNg_field_flt;


typedef struct _suNf_field {
  suNf *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} suNf_field;

typedef struct _suNf_field_flt {
  suNf_flt *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} suNf_field_flt;

typedef struct _spinor_field {
  suNf_spinor* ptr;
  geometry_descriptor* type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} spinor_field;

typedef struct _spinor_field_flt {
  suNf_spinor_flt *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} spinor_field_flt;

typedef struct _suNg_av_field {
  suNg_algebra_vector *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} suNg_av_field;

typedef struct _scalar_field {
  double *ptr;
  geometry_descriptor *type;
#ifdef WITH_MPI
  MPI_Request *comm_req; 
#endif
} scalar_field;


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
	for(i=(s)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s)=(s)->ptr+(s)->type->master_start[_PIECE_INDEX(i)]; \
	    i<=(s)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s)++ \
	   )

#define _TWO_SPINORS_FOR(s1,s2,i) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(s1)->type->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(s1)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s1)=(s1)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)], \
	              _SPINOR_PTR(s2)=(s2)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]; \
	    i<=(s1)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++ \
	   )

#define _THREE_SPINORS_FOR(s1,s2,s3,i) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	_TWO_SPINORS_MATCHING(s1,s3); \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(s1)->type->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(s1)->type->master_start[_PIECE_INDEX(i)], _SPINOR_PTR(s1)=(s1)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)], \
	              _SPINOR_PTR(s2)=(s2)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)], \
	              _SPINOR_PTR(s3)=(s3)->ptr+(s1)->type->master_start[_PIECE_INDEX(i)]; \
	    i<=(s1)->type->master_end[_PIECE_INDEX(i)]; \
	    i++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++, _SPINOR_PTR(s3)++ \
	   )

#include "field_ordering.h"

#define _FIELD_AT(s,i) (((s)->ptr)+i)
#define _4FIELD_AT(s,i,mu) (((s)->ptr)+coord_to_index(i,mu))
/*
#define _SPINOR_ADDR(s) ((s)->ptr)
*/

#endif
