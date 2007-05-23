/*******************************************************************************
*
* File field_alloc.c
*
* Functions for fields allocation
*
*******************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "global.h"

#ifdef P4
#  define ALIGN 7
#else
#  define ALIGN 5 
#endif  

void free_field(void *u)
{
   afree(u);
}

suNg* alloc_gfield()
{
   int ix;
   suNg unity,*p;

   p=amalloc(4*VOLUME*sizeof(suNg),ALIGN);
   error(p==NULL,1,"alloc_u [start.c]",
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
   error(p==NULL,1,"alloc_u [start.c]",
         "Could not allocate memory space for the gauge field");

   _suNf_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}


suNg_dble* alloc_gfield_dble()
{
   int ix;
   suNg_dble unity,*p;

   p=amalloc(4*VOLUME*sizeof(suNg_dble),ALIGN);
   error(p==NULL,1,"alloc_ud [start.c]",
         "Could not allocate memory space for the gauge field");

   _suNg_unit(unity);

   for (ix=0;ix<4*VOLUME;ix+=2){
       *(p+ix)=unity;
       *(p+ix+1)=unity;
   }
   
   return p;
}

suNf_spinor* alloc_spinor_field_f()
{
   suNf_spinor *p;

   p=amalloc(VOLUME*sizeof(suNf_spinor),ALIGN);
   error(p==NULL,1,"alloc_u [start.c]",
         "Could not allocate memory space for the spinor field");
   
   return p;
}

suNg_algebra_vector* alloc_momenta()
{
   suNg_algebra_vector *p;

   p=amalloc(4*VOLUME*sizeof(suNg_algebra_vector),ALIGN);
   error(p==NULL,1,"alloc_u [start.c]",
         "Could not allocate memory space for the momenta field");
   
   return p;
}
