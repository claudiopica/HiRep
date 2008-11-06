/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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
#include "spinor_field.h"
#include "geometry.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

void free_gfield(suNg_field *u)
{
   afree(u->ptr);
#ifdef WITH_MPI
   afree(u->comm_req);   
#endif
   afree(u);
}

suNg_field *alloc_gfield(geometry_descriptor* type)
{
   int ix;
   suNg unity;
   suNg_field *gf;

   gf=amalloc(sizeof(suNg_field),ALIGN);
   error(gf==NULL,1,"alloc_gfield [field_alloc.c]",
         "Could not allocate memory space for the gauge field (1)");
   gf->ptr=amalloc(4*type->gsize*sizeof(suNg),ALIGN);
   error((gf->ptr)==NULL,1,"alloc_gfield [field_alloc.c]",
         "Could not allocate memory space for the gauge field (2)");

   gf->type=type;

#ifdef WITH_MPI
   if (type->nbuffers>0) {
     gf->comm_req=amalloc(2*type->nbuffers*sizeof(MPI_Request),ALIGN);
     error((gf->comm_req)==NULL,1,"alloc_gfield [field_alloc.c]",
	   "Could not allocate memory space for the gauge field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       gf->comm_req[ix]=MPI_REQUEST_NULL;
   } else {
     gf->comm_req=NULL;
   }
#endif

   /* set gauge field to unity */
   _suNg_unit(unity);
   for (ix=0;ix<4*type->gsize;++ix)
     *((gf->ptr)+ix)=unity;

   return gf;
}

void free_gfield_f(suNf_field *u)
{
  afree(u->ptr);
#ifdef WITH_MPI
  afree(u->comm_req);   
#endif
  afree(u);
}

suNf_field* alloc_gfield_f(geometry_descriptor* type)
{
  int ix;
  suNf unity;
  suNf_field *gf;

  gf=amalloc(sizeof(suNf_field),ALIGN);
  error(gf==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (1)");
  gf->ptr=amalloc(4*type->gsize*sizeof(suNf),ALIGN);
  error((gf->ptr)==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (2)");

  gf->type=type;

#ifdef WITH_MPI
  if (type->nbuffers>0) {
    gf->comm_req=amalloc(2*type->nbuffers*sizeof(MPI_Request),ALIGN);
    error((gf->comm_req)==NULL,1,"alloc_gfield [field_alloc.c]",
	  "Could not allocate memory space for the gauge field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       gf->comm_req[ix]=MPI_REQUEST_NULL;
  } else {
    gf->comm_req=NULL;
  }
#endif

  /* set gauge field to unity */
  _suNf_unit(unity);
  for (ix=0;ix<4*type->gsize;++ix)
    *((gf->ptr)+ix)=unity;

  return gf;
}

void free_gfield_flt(suNg_field_flt *u)
{
  afree(u->ptr);
#ifdef WITH_MPI
  afree(u->comm_req);   
#endif
  afree(u);
}

suNg_field_flt* alloc_gfield_flt(geometry_descriptor* type)
{
  int ix;
  suNg_flt unity;
  suNg_field_flt *gf;

  gf=amalloc(sizeof(suNg_field_flt),ALIGN);
  error(gf==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (1)");
  gf->ptr=amalloc(4*type->gsize*sizeof(suNg_flt),ALIGN);
  error((gf->ptr)==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (2)");

  gf->type=type;

#ifdef WITH_MPI
  if (type->nbuffers>0) {
    gf->comm_req=amalloc(2*type->nbuffers*sizeof(MPI_Request),ALIGN);
    error((gf->comm_req)==NULL,1,"alloc_gfield [field_alloc.c]",
	  "Could not allocate memory space for the gauge field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       gf->comm_req[ix]=MPI_REQUEST_NULL;
  } else {
    gf->comm_req=NULL;
  }
#endif


  /* set gauge field to unity */
  _suNg_unit(unity);
  for (ix=0;ix<4*type->gsize;++ix)
    *((gf->ptr)+ix)=unity;

  return gf;
}

void free_gfield_f_flt(suNf_field_flt *u)
{
  afree(u->ptr);
#ifdef WITH_MPI
  afree(u->comm_req);   
#endif
  afree(u);
}

suNf_field_flt* alloc_gfield_f_flt(geometry_descriptor* type)
{
  int ix;
  suNf_flt unity;
  suNf_field_flt *gf;

  gf=amalloc(sizeof(suNf_field_flt),ALIGN);
  error(gf==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (1)");
  gf->ptr=amalloc(4*type->gsize*sizeof(suNf_flt),ALIGN);
  error((gf->ptr)==NULL,1,"alloc_gfield [field_alloc.c]",
	"Could not allocate memory space for the gauge field (2)");

  gf->type=type;

#ifdef WITH_MPI
  if (type->nbuffers>0) {
    gf->comm_req=amalloc(2*type->nbuffers*sizeof(MPI_Request),ALIGN);
    error((gf->comm_req)==NULL,1,"alloc_gfield [field_alloc.c]",
	  "Could not allocate memory space for the gauge field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       gf->comm_req[ix]=MPI_REQUEST_NULL;
  } else {
    gf->comm_req=NULL;
  }
#endif

  /* set gauge field to unity */
  _suNf_unit(unity);
  for (ix=0;ix<4*type->gsize;++ix)
    *((gf->ptr)+ix)=unity;

  return gf;
}

void free_avfield(suNg_av_field *u)
{
   afree(u->ptr);
#ifdef WITH_MPI
   afree(u->comm_req);   
#endif
   afree(u);
}

suNg_av_field *alloc_avfield(geometry_descriptor* type)
{
   int ix;
   suNg_av_field *af;

   af=amalloc(sizeof(*af),ALIGN);
   error(af==NULL,1,"alloc_avfield [field_alloc.c]",
         "Could not allocate memory space for the av field (1)");
   af->ptr=amalloc(4*type->gsize*sizeof(*(af->ptr)),ALIGN);
   error((af->ptr)==NULL,1,"alloc_avfield [field_alloc.c]",
         "Could not allocate memory space for the av field (2)");

   af->type=type;

#ifdef WITH_MPI
   if (type->nbuffers>0) {
     af->comm_req=amalloc(2*type->nbuffers*sizeof(*(af->comm_req)),ALIGN);
     error((af->comm_req)==NULL,1,"alloc_avfield [field_alloc.c]",
	   "Could not allocate memory space for the av field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       af->comm_req[ix]=MPI_REQUEST_NULL;
   } else {
     af->comm_req=NULL;
   }
#endif

   return af;
}

void free_sfield(scalar_field *u)
{
   afree(u->ptr);
#ifdef WITH_MPI
   afree(u->comm_req);   
#endif
   afree(u);
}

scalar_field *alloc_sfield(geometry_descriptor* type)
{
   int ix;
   scalar_field *af;

   af=amalloc(sizeof(*af),ALIGN);
   error(af==NULL,1,"alloc_avfield [field_alloc.c]",
         "Could not allocate memory space for the av field (1)");
   af->ptr=amalloc(type->gsize*sizeof(*(af->ptr)),ALIGN);
   error((af->ptr)==NULL,1,"alloc_avfield [field_alloc.c]",
         "Could not allocate memory space for the av field (2)");

   af->type=type;

#ifdef WITH_MPI
   if (type->nbuffers>0) {
     af->comm_req=amalloc(2*type->nbuffers*sizeof(*(af->comm_req)),ALIGN);
     error((af->comm_req)==NULL,1,"alloc_avfield [field_alloc.c]",
	   "Could not allocate memory space for the av field (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       af->comm_req[ix]=MPI_REQUEST_NULL;
   } else {
     af->comm_req=NULL;
   }
#endif

   return af;
}


suNg_field *alloc_gtransf(geometry_descriptor* type)
{
   int ix;
   suNg unity;
   suNg_field *gf;

   gf=amalloc(sizeof(suNg_field),ALIGN);
   error(gf==NULL,1,"alloc_gtransf [field_alloc.c]",
         "Could not allocate memory space for the gauge transformation (1)");
   gf->ptr=amalloc(type->gsize*sizeof(suNg),ALIGN);
   error((gf->ptr)==NULL,1,"alloc_gtransf [field_alloc.c]",
         "Could not allocate memory space for the gauge transformation (2)");

   gf->type=type;

#ifdef WITH_MPI
   if (type->nbuffers>0) {
     gf->comm_req=amalloc(2*type->nbuffers*sizeof(MPI_Request),ALIGN);
     error((gf->comm_req)==NULL,1,"alloc_gtransf [field_alloc.c]",
	   "Could not allocate memory space for the gauge transformation (3)");
     for (ix=0; ix<2*type->nbuffers; ++ix)
       gf->comm_req[ix]=MPI_REQUEST_NULL;
   } else {
     gf->comm_req=NULL;
   }
#endif

   /* set gauge field to unity */
   _suNg_unit(unity);
   for (ix=0;ix<type->gsize;++ix)
     *((gf->ptr)+ix)=unity;

   return gf;
}

void free_gtransf(suNg_field *u)
{
  afree(u->ptr);
#ifdef WITH_MPI
  afree(u->comm_req);   
#endif
  afree(u);
}

