/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random_fields.c
*
* Pseudorandom generation of fields
*
*******************************************************************************/

#include "random.h"
#include "libhr_core.h"
#include "memory.h"

void random_u(suNg_field *gf)
{
  error(gf==NULL,1,"random_u [random_fields.c]",
	"Attempt to access unallocated memory space");   
  
  _MASTER_FOR(gf->type,ix) {
    /* unroll 4 directions */
    suNg *ptr=(gf->ptr)+coord_to_index(ix,0);
    random_suNg(ptr++);
    random_suNg(ptr++);
    random_suNg(ptr++);
    random_suNg(ptr);
  }

  #ifdef WITH_GPU
    copy_to_gpu_gfield(gf);
  #endif

  start_sendrecv_gfield(gf);
}

void random_u_f(suNf_field *gf) 
{
  error(gf==NULL, 1, "random_u_f [random_fields.c]", 
  "Attempt to access unallocated memory space");

  _MASTER_FOR(gf->type,ix) 
  {
    suNf *ptr = (gf->ptr)+coord_to_index(ix, 0);
    random_suNf(ptr++);
    random_suNf(ptr++);
    random_suNf(ptr++);
    random_suNf(ptr);
  }

  #ifdef WITH_GPU
    copy_to_gpu_gfield_f(gf);
  #endif

  //start_sendrecv_gfield(gf);
}

void unit_u(suNg_field *gf)
{
  suNg unity;

  error(gf==NULL,1,"unit_u [random_fields.c]",
	"Attempt to access unallocated memory space");   
  
   _suNg_unit(unity);
  _MASTER_FOR(gf->type,ix) {
    /* unroll 4 directions */
    suNg *ptr=(gf->ptr)+coord_to_index(ix,0);
    *(ptr++)=unity;
    *(ptr++)=unity;
    *(ptr++)=unity;
    *(ptr)=unity;
  }

  #ifdef WITH_GPU
    copy_to_gpu_gfield(gf);
  #endif

  start_sendrecv_gfield(gf);

}

void random_s(suNg_scalar_field *sf)
{
	error(sf==NULL,1,"random_s [random_fields.c]","Attempt to access unallocated memory space");

	_MASTER_FOR(sf->type,ix)
	{
		suNg_vector *ptr = _FIELD_AT(sf,ix);
		gaussian_suNg_vector(ptr);
	}

  #ifdef WITH_GPU
    copy_to_gpu_suNg_scalar_field(sf);
  #endif

	start_sendrecv_suNg_scalar_field(sf);
	complete_sendrecv_suNg_scalar_field(sf);
}

void zero_s(suNg_scalar_field *sf)
{
	error(sf==NULL,1,"zero_s [random_fields.c]","Attempt to access unallocated memory space");

	suNg_vector zero_vector;
	_vector_zero_g(zero_vector);

	_MASTER_FOR(sf->type,ix)
	{
		suNg_vector *ptr = _FIELD_AT(sf,ix);
		*ptr = zero_vector;
	}

  #ifdef WITH_GPU
    copy_to_gpu_suNg_scalar_field(sf);
  #endif

	start_sendrecv_suNg_scalar_field(sf);
	complete_sendrecv_suNg_scalar_field(sf);
}
