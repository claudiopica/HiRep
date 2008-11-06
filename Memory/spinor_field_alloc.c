/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "geometry.h"
#include "linear_algebra.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif


spinor_field* alloc_spinor_field_f(unsigned int n, geometry_descriptor* type)
{
	suNf_spinor *p;
	spinor_field *s;
	unsigned int i;
#ifdef WITH_MPI /* MPI variables */
	MPI_Request *r;
#endif 
	s=amalloc(n*sizeof(spinor_field),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*type->gsize*sizeof(suNf_spinor),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");

#ifdef WITH_MPI
	if ((type->nbuffers)>0) {
	  r=amalloc(n*2*type->nbuffers*sizeof(MPI_Request),ALIGN);
	  error(r==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
		"Could not allocate memory space for the spinor field (3)");
	  for (i=0; i<n*2*type->nbuffers; ++i)
	    r[i]=MPI_REQUEST_NULL;
	} else {
	  r=NULL;
	}	
#endif

	for(i=0; i<n; ++i) {
	  s[i].ptr=p+i*type->gsize;
	  s[i].type=type;
#ifdef WITH_MPI
	  if (r==NULL) {
	    s[i].comm_req=NULL;
	  } else {
	    s[i].comm_req=r+i*2*type->nbuffers;
	  }
#endif
	}

	return s;
}

spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n, geometry_descriptor *type)
{
	suNf_spinor_flt *p;
	spinor_field_flt *s;
	unsigned int i;
#ifdef WITH_MPI /* MPI variables */
	MPI_Request *r;
#endif 

	s=amalloc(n*sizeof(spinor_field_flt),ALIGN);
	error(s==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (1)");
	p=amalloc(n*type->gsize*sizeof(suNf_spinor_flt),ALIGN);
	error(p==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		  "Could not allocate memory space for the spinor field (2)");
	
#ifdef WITH_MPI
	if (type->nbuffers>0) {
	  r=amalloc(n*2*type->nbuffers*sizeof(MPI_Request),ALIGN);
	  error(r==NULL,1,"alloc_spinor_field_f_flt [field_alloc.c]",
		"Could not allocate memory space for the spinor field (3)");
	  for (i=0; i<n*2*type->nbuffers; ++i)
	    r[i]=MPI_REQUEST_NULL;
	} else {
	  r=NULL;
	}
#endif

	for(i=0; i<n; ++i) {
	  s[i].ptr=p+i*type->gsize;
	  s[i].type=type;
#ifdef WITH_MPI
	  if (r==NULL) {
	    s[i].comm_req=NULL;
	  } else {
	    s[i].comm_req=r+i*2*type->nbuffers;
	  }	
#endif
	}
	
	return s;
}

void free_spinor_field(spinor_field *s)
{
	afree(s->ptr);
#ifdef WITH_MPI
	afree(s->comm_req);
#endif
	afree(s);
}

void free_spinor_field_flt(spinor_field_flt *s)
{
	afree(s->ptr);
#ifdef WITH_MPI
	afree(s->comm_req);
#endif
	afree(s);
}


