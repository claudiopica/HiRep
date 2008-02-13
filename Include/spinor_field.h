#ifndef SPINOR_FIELD_H
#define SPINOR_FIELD_H

#include "geometry.h"
#include "error.h"

typedef struct {
	spinor_descriptor* type;
	suNf_spinor* ptr;
} spinor_field;

typedef struct {
	spinor_descriptor* type;
	suNf_spinor_flt* ptr;
} spinor_field_flt;


#define _LAT_INDEX(i) i##_latindex
#define _SPINOR_PTR(s) s##_ptr

#define _DECLARE_INT_ITERATOR(i) int i, _LAT_INDEX(i)
#define _DECLARE_SPINOR_ITERATOR(s) suNf_spinor _SPINOR_PTR(s)

#define _TWO_SPINORS_MATCHING(s1,s2) \
	error((s1)->type!=(s2)->type,1,"not available", "Spinors don't match!");

#define _ONE_SPINOR_FOR(s,i,type) \
	for(i=0,type=(s)->type;i<type->local_pieces;i++) \
	for(_LAT_INDEX(i)=type->start[i], _SPINOR_PTR(s)=(s)->ptr+type->start[i]; \
	    _LAT_INDEX(i)<=type->end[i]; \
	    _LAT_INDEX(i)++, _SPINOR_PTR(s)++ \
	   )

#define _TWO_SPINORS_FOR(s1,s2,i,type) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	for(i=0,type=(s1)->type;i<type->local_pieces;i++) \
	for(_LAT_INDEX(i)=type->start[i], _SPINOR_PTR(s1)=(s1)->ptr+type->start[i], \
	              _SPINOR_PTR(s2)=(s2)->ptr+type->start[i]; \
	    _LAT_INDEX(i)<=type->end[i]; \
	    _LAT_INDEX(i)++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++ \
	   )

#define _THREE_SPINORS_FOR(s1,s2,s3,i,type) \
	_TWO_SPINORS_MATCHING(s1,s2); \
	_TWO_SPINORS_MATCHING(s1,s3); \
	for(i=0,type=(s1)->type;i<type->local_pieces;i++) \
	for(_LAT_INDEX(i)=type->start[i], _SPINOR_PTR(s1)=(s1)->ptr+type->start[i], \
	              _SPINOR_PTR(s2)=(s2)->ptr+type->start[i], \
	              _SPINOR_PTR(s3)=(s3)->ptr+type->start[i]; \
	    _LAT_INDEX(i)<=type->end[i]; \
	    _LAT_INDEX(i)++, _SPINOR_PTR(s1)++, _SPINOR_PTR(s2)++, _SPINOR_PTR(s3)++ \
	   )

#define _SPINOR_AT(s,i) (((s)->ptr)+i)
#define _SPINOR_ADDR(s) ((s)->ptr)


#endif
