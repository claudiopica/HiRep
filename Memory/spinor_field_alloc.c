#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "geometry.h"
#include "linear_algebra.h"

spinor_field* alloc_spinor_field_f(unsigned int n, geometry_descriptor* type)
{
	suNf_spinor *p;
	spinor_field *s;
	unsigned int i;


	s=amalloc(n*sizeof(spinor_field),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*type->gsize*sizeof(suNf_spinor),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");
	
	for(i=0; i<n; i++) {
		s[i].ptr=p+i*type->gsize;
		s[i].type=type;
	}

	return s;
}

spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n, geometry_descriptor *type)
{
	suNf_spinor_flt *p;
	spinor_field_flt *s;
	unsigned int i;

	s=amalloc(n*sizeof(spinor_field_flt),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*type->gsize*sizeof(suNf_spinor_flt),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");
	
	for(i=0; i<n; i++) {
		s[i].ptr=p+i*type->gsize;
		s[i].type=type;
	}
	
	return s;
}

void free_spinor_field(spinor_field *s)
{
	afree(s->ptr);
	afree(s);
}

void free_spinor_field_flt(spinor_field_flt *s)
{
	afree(s->ptr);
	afree(s);
}


