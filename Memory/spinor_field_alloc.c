#include <stdlib.h>
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "memory.h"
#include "global.h"


/*****  spinor_field = suNf_spinor*  *****/
spinor_field* alloc_spinor_field_f(unsigned int n)
{
	suNf_spinor *p;
	unsigned int len;
	spinor_field *s;
	unsigned int i;

	get_spinor_len(&len);

	s=amalloc(n*sizeof(spinor_field),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*len*sizeof(suNf_spinor),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");
	
	for(i=0; i<n; i++)
		s[i]=p+i*len;

	return s;
}

spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n)
{
	suNf_spinor_flt *p;
	unsigned int len;
	spinor_field_flt *s;
	unsigned int i;

	get_spinor_len(&len);

	s=amalloc(n*sizeof(spinor_field_flt),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*len*sizeof(suNf_spinor_flt),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");
	
	for(i=0; i<n; i++)
		s[i]=p+i*len;

	return s;
}

void free_spinor_field(spinor_field *s)
{
	afree(*s);
	afree(s);
}

void free_spinor_field_flt(spinor_field_flt *s)
{
	afree(*s);
	afree(s);
}


/*****  spinor_field = suNf_spinor
spinor_field* alloc_spinor_field_f(unsigned int n)
{
	suNf_spinor *p;
	unsigned int len;

	get_spinor_len(&len);
	p=amalloc(n*len*sizeof(suNf_spinor),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field");

	return p;
}

spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n)
{
	suNf_spinor_flt *p;
	unsigned int len;

	get_spinor_len(&len);
	p=amalloc(n*len*sizeof(suNf_spinor_flt),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		   "Could not allocate memory space for the spinor field");

	return p;
}
*/


