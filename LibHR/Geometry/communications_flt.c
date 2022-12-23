/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "communications.h"
#include "logger.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include "suN_types.h"
#include "global.h"
#include "utils.h"
#include <string.h>
#ifdef WITH_MPI
#include <hr_mpi.h>
#endif
#include "logger.h"

#ifdef MPI_TIMING
struct timeval gfstart, gfend, gfetime,sfstart, sfend, sfetime;
int gf_control=0,sf_control=0;
#endif


/* functions for filling sending buffer */
#ifdef WITH_MPI
#ifdef WITH_NEW_GEOMETRY
static void sync_gauge_field_flt(suNg_field_flt *gf) {
  sync_field(gf->type, 4*sizeof(*gf->ptr), 0, gf->ptr, gf->sendbuf_ptr);
}
#else
static void sync_gauge_field_flt(suNg_field_flt *gf) {
  geometry_descriptor *gd=gf->type;
  for(int i=0; i<gd->ncopies_gauge; ++i) {
    /* this assumes that the 4 directions are contiguous in memory !!! */
    memcpy(((gf->ptr)+4*gd->copy_to[i]),((gf->ptr)+4*gd->copy_from[i]),4*(gd->copy_len[i])*sizeof(*(gf->ptr)));
  }
}
#endif

#ifdef WITH_NEW_GEOMETRY
static void sync_spinor_field_flt(spinor_field_flt *p) {
  sync_field(p->type, sizeof(*p->ptr), 1, p->ptr, p->sendbuf_ptr);
}
#else
static void sync_spinor_field_flt(spinor_field_flt *p) {
  geometry_descriptor *gd = p->type;
  for(int i=0; i<gd->ncopies_spinor; ++i) {
    memcpy((p->ptr+gd->copy_to[i]-gd->master_shift),(p->ptr+gd->copy_from[i]-gd->master_shift),(gd->copy_len[i])*sizeof(*(p->ptr)));
  }
}
#endif

#endif /* WITH_MPI */

#if defined(WITH_NEW_GEOMETRY) && defined(WITH_MPI)
static void *gf_sendrecv_guard=NULL;
#endif

void complete_gf_sendrecv_flt(suNg_field_flt *gf) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int nreq=2*gf->type->nbuffers_gauge;

#ifdef WITH_NEW_GEOMETRY
  if (gf_sendrecv_guard!=NULL && gf_sendrecv_guard != gf->comm_req)  
    error(1, 1, "complete_gf_sendrecv " __FILE__, "More simultaneous communication attempted. Existing...\n");
  gf_sendrecv_guard = NULL;
#endif

  if(nreq>0) {
    MPI_Status status[nreq];

    mpiret=MPI_Waitall(nreq, gf->comm_req, status);

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<nreq; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k, 
              status[k].MPI_SOURCE, 
              status[k].MPI_TAG, 
              mesg);
        }
      }
      error(1,1,"complete_gf_sendrecv_flt" __FILE__,"Cannot complete communications");
    }
#endif
  }

#ifdef MPI_TIMING
    if(gf_control>0)
      {
	gettimeofday(&gfend,0);
	timeval_subtract(&gfetime,&gfend,&gfstart);
	lprintf("MPI TIMING",0,"complete_gf_sendrecv" __FILE__ " %ld sec %ld usec\n",gfetime.tv_sec,gfetime.tv_usec);
	gf_control=0;
      }
#endif

#endif /* WITH_MPI */
}

void start_gf_sendrecv_flt(suNg_field_flt *gf) {
#ifdef WITH_MPI
  int i, mpiret; (void)mpiret; // Remove warning of variable set but not used
  geometry_descriptor *gd=gf->type;

  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_gf_sendrecv_flt(gf);

#ifdef WITH_NEW_GEOMETRY
  gf_sendrecv_guard=(void*)gf->comm_req;
#endif

  /* fill send buffers */
  sync_gauge_field_flt(gf);

#ifdef MPI_TIMING
  error(gf_control>0,1,"start_gf_sendrecv_flt" __FILE__,"Multiple send without receive");
  gettimeofday(&gfstart,0);  
  gf_control=1;
#endif

  for (i=0; i<(gd->nbuffers_gauge); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((float*)((gf->sendbuf_ptr)+4*gd->sbuf_start[i]), /* buffer */
        (gd->sbuf_len[i])*sizeof(suNg_flt)/sizeof(float)*4,   /* lenght in units of flaots */
        MPI_FLOAT,                                            /* basic datatype */
        gd->sbuf_to_proc[i],                                  /* cid of destination */
        i,                                                    /* tag of communication */
        cart_comm,                                            /* use the cartesian communicator */
        &(gf->comm_req[2*i])                                  /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv_flt" __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((float*)((gf->ptr)+4*gd->rbuf_start[i]), /* buffer */
        (gd->rbuf_len[i])*sizeof(suNg_flt)/sizeof(float)*4,   /* lenght in units of floats */
        MPI_FLOAT,                                            /* basic datatype */
        gd->rbuf_from_proc[i],                                /* cid of origin */
        i,                                                    /* tag of communication */
        cart_comm,                                            /* use the cartesian communicator */
        &(gf->comm_req[2*i+1])                                /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv_flt" __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}

#if defined(WITH_NEW_GEOMETRY) && defined(WITH_MPI)
static void *sf_sendrecv_guard=NULL;
#endif

void complete_sf_sendrecv_flt(spinor_field_flt *sf) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int nreq=2*sf->type->nbuffers_spinor;

#ifdef WITH_NEW_GEOMETRY
  if (sf_sendrecv_guard!=NULL && sf_sendrecv_guard != sf->comm_req) {
    print_trace();
    error(1, 1, "complete_sf_sendrecv " __FILE__, "More simultaneous communication attempted. Existing...\n");
  }
  sf_sendrecv_guard = NULL;
#endif


if(nreq>0) {
    MPI_Status status[nreq];

    mpiret=MPI_Waitall(nreq, sf->comm_req, status);

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<nreq; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k, 
              status[k].MPI_SOURCE, 
              status[k].MPI_TAG, 
              mesg);
        }
      }
      error(1,1,"complete_sf_sendrecv_flt" __FILE__,"Cannot complete communications");
    }
#endif
  }

#ifdef MPI_TIMING
    if(sf_control>0)
      {
	gettimeofday(&sfend,0);
	timeval_subtract(&sfetime,&sfend,&sfstart);
	lprintf("MPI TIMING",0,"complete_sf_sendrecv_flt" __FILE__ " %ld sec %ld usec\n",sfetime.tv_sec,sfetime.tv_usec);
	sf_control=0;
      }
#endif

#endif /* WITH_MPI */
}

void start_sf_sendrecv_flt(spinor_field_flt *sf) {
#ifdef WITH_MPI
  int i, mpiret; (void)mpiret; // Remove warning of variable set but not used
  geometry_descriptor *gd=sf->type;


  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_sf_sendrecv_flt(sf);

#ifdef WITH_NEW_GEOMETRY
  sf_sendrecv_guard=(void*)sf->comm_req;
#endif

  /* fill send buffers */
  sync_spinor_field_flt(sf);
#ifdef MPI_TIMING
  error(sf_control>0,1,"start_sf_sendrecv_flt" __FILE__,"Multiple send without receive");
  gettimeofday(&sfstart,0);  
  sf_control=1;
#endif

  for (i=0; i<(gd->nbuffers_spinor); ++i) {
    mpiret=MPI_Irecv((float*)((sf->ptr)+(gd->rbuf_start[i])-(gd->master_shift)),  /* buffer */
        (gd->rbuf_len[i])*(sizeof(suNf_spinor_flt)/sizeof(float)),                /* lenght in units of float */
        MPI_FLOAT,                                                                /* basic datatype */
        gd->rbuf_from_proc[i],                                                    /* cid of origin */
        i,                                                                        /* tag of communication */
        cart_comm,                                                                /* use the cartesian communicator */
        &(sf->comm_req[2*i+1])                                                     /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv_flt" __FILE__,"Cannot start receive buffer");
    }
#endif
  }

  int shift = gd->master_shift;
#ifdef WITH_NEW_GEOMETRY
  shift = 0;
#endif

  for (i=0; i<(gd->nbuffers_spinor); ++i) {
    mpiret=MPI_Isend((float*)((sf->sendbuf_ptr)+(gd->sbuf_start[i])-shift),  /* buffer */
        (gd->sbuf_len[i])*(sizeof(suNf_spinor_flt)/sizeof(float)),                /* lenght in units of floats */
        MPI_FLOAT,                                                                /* basic datatype */
        gd->sbuf_to_proc[i],                                                      /* cid of destination */
        i,                                                                        /* tag of communication */
        cart_comm,                                                                /* use the cartesian communicator */
        &(sf->comm_req[2*i])                                                     /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv_flt" __FILE__,"Cannot start send buffer");
    }
#endif
  }

#endif /* WITH_MPI */
}

