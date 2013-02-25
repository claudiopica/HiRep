/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"
#include <assert.h>


#define _PROJ_BIT (1<<4) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_gfield(c) if((c)&_PROJ_BIT) \
                          {(c)=0;project_gauge_field();} \
                        else \
                          {++(c); start_gf_sendrecv(u_gauge);} \
                        complete_gf_sendrecv(u_gauge);

const double lambda = 0.1931833275037836; /* Omelyan et al */ 
                                         /* lamdba = 1./6. Sexton Weingarten */ 

static unsigned int count=0; /*used to project gfield */



void gauge_integrator(suNg_av_field *momenta, double tlen, integrator_par *int_par){
  #ifdef TIMING
  struct timeval start, end;
  struct timeval start1, end1;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  _DECLARE_INT_ITERATOR(i);

  int n;
  double dt=tlen/((double)int_par->nsteps);

  if(int_par->nsteps == 0)  return;

  lprintf("MD_INT",10+int_par->level*10,"Starting new MD trajectory with gauge_integrator, level %d.\n",int_par->level);
  lprintf("MD_INT",10+int_par->level*10,"MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",
    int_par->level,tlen,int_par->nsteps,dt);
  for(n=1;n<=int_par->nsteps;++n) {

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&start1,0);
    #endif

    _MASTER_FOR(&glattice,i) {
      ExpX(dt,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
      ExpX(dt,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
      ExpX(dt,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
      ExpX(dt,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
    }

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&end1,0);
    timeval_subtract(&etime,&end1,&start1);
    lprintf("TIMING",0,"exp in gauge_integrator %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
    #endif

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&start1,0);
    #endif

    _proj_gfield(count);

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&end1,0);
    timeval_subtract(&etime,&end1,&start1);
    lprintf("TIMING",0,"proj_gfield in gauge_integrator %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
    #endif

    lprintf("MD_INT",20+int_par->level*10,"MD step (level %d): %d/%d\n",int_par->level,n,int_par->nsteps);
   }
  lprintf("MD_INT",10+int_par->level*10,"MD trajectory completed, level %d.\n",int_par->level);

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start1,0);
  #endif

  represent_gauge_field();

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end1,0);
  timeval_subtract(&etime,&end1,&start1);
  lprintf("TIMING",0,"represent_gauge_field in gauge_integrator %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"gauge_integrator %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}




void O2MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par){
  #ifdef TIMING
  struct timeval start, end;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
   _TWO_SPINORS_MATCHING(u_gauge,momenta);
#endif

  if(int_par->nsteps == 0)  return;
  
  int n;
  double dt=tlen/((double)int_par->nsteps);
    
  lprintf("MD_INT",10+int_par->level*10,"Starting new MD trajectory with O2MN_multistep, level %d.\n",int_par->level);
  lprintf("MD_INT",10+int_par->level*10,"MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",
    int_par->level,tlen,int_par->nsteps,dt);

  (*int_par->force)(lambda*dt,momenta,int_par->force_par);
  
  for(n=1;n<=int_par->nsteps;++n) {
    
    (*int_par->next->integrator)(momenta, dt/2., int_par->next);
    
    (*int_par->force)((1.-2.*lambda)*dt,momenta,int_par->force_par);
    
    (*int_par->next->integrator)(momenta, dt/2., int_par->next);
    
    /* Update of momenta */
    if(n<int_par->nsteps) {
      (*int_par->force)(2.*lambda*dt,momenta,int_par->force_par);
      lprintf("MD_INT",20+int_par->level*10,"MD step (level %d): %d/%d\n",int_par->level,n,int_par->nsteps);
    } else {
      (*int_par->force)(lambda*dt,momenta,int_par->force_par);
      lprintf("MD_INT",10+int_par->level*10,"MD trajectory completed, level %d.\n",int_par->level);
    }
  }

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"O2MN_multistep[%d] %.6f s\n",int_par->level,1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}




#undef _PROJ_BIT
#undef _proj_gfield
