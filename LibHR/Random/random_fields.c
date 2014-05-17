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

#include <stdlib.h>

#include "global.h"
#include "random.h"
#include "error.h"
#include "field_ordering.h"
#include "spinor_field.h"
#include "communications.h"
#include <math.h>

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

  start_gf_sendrecv(gf);

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

  start_gf_sendrecv(gf);

}

