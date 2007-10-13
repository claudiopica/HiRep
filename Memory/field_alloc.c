/*******************************************************************************
*
* File field_alloc.c
*
* Functions for fields allocation
*
*******************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "memory.h"
#include "global.h"


void free_field(void *u)
{
   afree(u);
}

suNg* alloc_gfield()
{
   int ix;
   suNg unity,*p;

   p=amalloc(4*VOLUME*sizeof(suNg),ALIGN);
   error(p==NULL,1,"alloc_gfield [field_alloc.c]",
         "Could not allocate memory space for the gauge field");

   _suNg_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}

suNf* alloc_gfield_f()
{
   int ix;
   suNf unity,*p;

   p=amalloc(4*VOLUME*sizeof(suNf),ALIGN);
   error(p==NULL,1,"alloc_gfield_f [field_alloc.c]",
         "Could not allocate memory space for the gauge field");

   _suNf_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}


suNg_flt* alloc_gfield_flt()
{
   int ix;
   suNg_flt unity,*p;

   p=amalloc(4*VOLUME*sizeof(suNg_flt),ALIGN);
   error(p==NULL,1,"alloc_gfield_flt [field_alloc.c]",
         "Could not allocate memory space for the gauge field");

   _suNg_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}

suNf_flt* alloc_gfield_f_flt()
{
   int ix;
   suNf_flt unity,*p;
 
   p=amalloc(4*VOLUME*sizeof(suNf_flt),ALIGN);
   error(p==NULL,1,"alloc_gfield_f_flt [field_alloc.c]",
         "Could not allocate memory space for the gauge field");

   _suNf_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}

suNg_algebra_vector* alloc_momenta()
{
   suNg_algebra_vector *p;

   p=amalloc(4*VOLUME*sizeof(suNg_algebra_vector),ALIGN);
   error(p==NULL,1,"alloc_momenta [field_alloc.c]",
         "Could not allocate memory space for the momenta field");
   
   return p;
}
