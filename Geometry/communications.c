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
#include <string.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include "logger.h"

void global_sum(double *d, int n) {
#ifdef WITH_MPI
  int mpiret;
  double pres[n];

  mpiret=MPI_Allreduce(d,pres,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"global_sum " __FILE__,"Cannot perform global_sum");
  }
#endif
  while(n>0) { 
    --n;
    d[n]=pres[n];	
  }
#else
  /* for non mpi do nothing */
  return;
#endif
}

void bcast(double *d, int n) {
#ifdef WITH_MPI
  int mpiret;

  mpiret=MPI_Bcast(d, n, MPI_DOUBLE, 0,MPI_COMM_WORLD);
#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"bcast " __FILE__,"Cannot perform global_sum");
  }
#endif

#else
  /* for non mpi do nothing */
  return;
#endif
}

/* functions for filling sending buffer */
static void sync_gauge_field(suNg_field *gf) {
  int i;
  geometry_descriptor *gd=gf->type;
  /* int j, mu, x, y; */

  for(i=0; i<gd->ncopies; ++i) {
    /* this assumes that the 4 directions are contiguous in memory !!! */
    memcpy(((gf->ptr)+4*gd->copy_to[i]),((gf->ptr)+4*gd->copy_from[i]),4*(gd->copy_len[i])*sizeof(*(gf->ptr)));
    /*   
         for(j=0; j<gd->copy_len[i]; j++) {
         x=gd->copy_from[i]+j;
         y=gd->copy_to[i]+j;
         for(mu=0; mu<4; mu++)
     *pu_gauge(y,mu) = *pu_gauge(x,mu);
     }
     */
  }
}

static void sync_spinor_field(spinor_field *p) {
  int i;
  /* int j, x, y; */
  geometry_descriptor *gd = p->type;

  for(i=0; i<gd->ncopies; ++i) {
    memcpy((p->ptr+gd->copy_to[i]),(p->ptr+gd->copy_from[i]),(gd->copy_len[i])*sizeof(*(p->ptr)));
    /*
       for(j=0; j<gd->copy_len[i]; j++) {
       x=gd->copy_from[i]+j;
       y=gd->copy_to[i]+j;
       p->ptr[y] = p->ptr[x];
       }
       */
  }
}

static void sync_gauge_transf(suNg_field *gf) {
  int i;
  geometry_descriptor *gd=gf->type;

  for(i=0; i<gd->ncopies; ++i) {
    memcpy(((gf->ptr)+gd->copy_to[i]),((gf->ptr)+gd->copy_from[i]),(gd->copy_len[i])*sizeof(*(gf->ptr)));
  }
}

static void sync_spinor_field_flt(spinor_field_flt *p) {
  int i;
  /* int j, x, y; */
  geometry_descriptor *gd = p->type;

  for(i=0; i<gd->ncopies; ++i) {
    memcpy((p->ptr+gd->copy_to[i]),(p->ptr+gd->copy_from[i]),(gd->copy_len[i])*sizeof(*(p->ptr)));
    /*
       for(j=0; j<gd->copy_len[i]; j++) {
       x=gd->copy_from[i]+j;
       y=gd->copy_to[i]+j;
       p->ptr[y] = p->ptr[x];
       }
       */
  }
}

/* This variable contains the information of the current status of communications
 * Values:
 * 0 => No communications pending
 *
 */ 
/* static unsigned int comm_status=0; */

void complete_gf_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int mpiret;
  int nreq=2*gf->type->nbuffers;

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
      error(1,1,"complete_gf_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#endif /* WITH_MPI */
}

void start_gf_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int i, mpiret;
  geometry_descriptor *gd=gf->type;

  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_gf_sendrecv(gf);

  /* fill send buffers */
  sync_gauge_field(gf);

  for (i=0; i<(gd->nbuffers); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((gf->ptr)+4*gd->sbuf_start[i], /* buffer */
        (gd->sbuf_len[i])*sizeof(suNg)/sizeof(double)*4, /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((gf->ptr)+4*gd->rbuf_start[i], /* buffer */
        (gd->rbuf_len[i])*sizeof(suNg)/sizeof(double)*4, /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}

void complete_sf_sendrecv(spinor_field *sf) {
#ifdef WITH_MPI
  int mpiret;
  int nreq=2*sf->type->nbuffers;

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
      error(1,1,"complete_gf_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#endif /* WITH_MPI */
}

void start_sf_sendrecv(spinor_field *sf) {
#ifdef WITH_MPI
  int i, mpiret;
  geometry_descriptor *gd=sf->type;


  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_sf_sendrecv(sf);

  /* fill send buffers */
  sync_spinor_field(sf);

  for (i=0; i<(gd->nbuffers); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((sf->ptr)+(gd->sbuf_start[i]), /* buffer */
        (gd->sbuf_len[i])*(sizeof(suNf_spinor)/sizeof(double)), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(sf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((sf->ptr)+(gd->rbuf_start[i]), /* buffer */
        (gd->rbuf_len[i])*(sizeof(suNf_spinor)/sizeof(double)), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(sf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}


void complete_gt_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int mpiret;
  int nreq=2*gf->type->nbuffers;

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
      error(1,1,"complete_gt_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#endif /* WITH_MPI */
}

void start_gt_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int i, mpiret;
  geometry_descriptor *gd=gf->type;

  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_gt_sendrecv(gf);

  /* fill send buffers */
  sync_gauge_transf(gf);

  for (i=0; i<(gd->nbuffers); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((gf->ptr)+gd->sbuf_start[i], /* buffer */
        (gd->sbuf_len[i])*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gt_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((gf->ptr)+gd->rbuf_start[i], /* buffer */
        (gd->rbuf_len[i])*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gt_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}

