/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica, Ari Hietanen          *   
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

static unsigned int count = 0; /*used to project gfield */
#define _PROJ_BIT (1<<4) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_gfield(c) if((c)&_PROJ_BIT) \
                          {(c)=0;project_gauge_field();} \
                        else \
                          {++(c); start_gf_sendrecv(u_gauge);} \
                        complete_gf_sendrecv(u_gauge);

void update_gauge_field(suNg_av_field *momenta, double dt) {
#ifdef TIMING
	struct timeval start, end;
	struct timeval start1, end1;
	struct timeval etime;

#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
#endif
	gettimeofday(&start,0);
#endif

	_MASTER_FOR(&glattice,i)
	{
		ExpX(dt,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
		ExpX(dt,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
		ExpX(dt,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
		ExpX(dt,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
	}

	_proj_gfield(count);
	represent_gauge_field();

#ifdef TIMING
#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
#endif
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"gauge_integrator %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
#endif
}

void gforce(suNg_av_field *force, double dt, int nmon, const monomial **mon_list)
{
  for(int n = 0; n < nmon; n++)
  {
    const monomial *m = mon_list[n];
    switch(m->type)
    {
      case PureGauge:
        force0(dt, force, &((mon_pg_par*)(m->par))->beta);
        break;
      case HMC:
        force_hmc(dt, force, &((mon_hmc_par*)m->par)->fpar);
        break;
      case RHMC:
        force_rhmc(dt, force, &((mon_rhmc_par*)m->par)->fpar);
        break;
      case TM:
        force_hmc_tm(dt,force, &((mon_tm_par*)m->par)->fpar);
        break;
      case TM_alt:
        force_hmc(dt,force, &((mon_tm_par*)m->par)->fpar);
        break;
      case Hasenbusch:
        force_hmc(dt, force, &((mon_hasenbusch_par*)m->par)->fpar);
        break;
      case Hasenbusch_tm:
        force_hmc_tm(dt,force, &((mon_hasenbusch_tm_par*)m->par)->fpar);
        break;
      case Hasenbusch_tm_alt:
        force_hmc(dt, force, &((mon_hasenbusch_tm_par*)m->par)->fpar);
      default:
        lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
        break;
    }
  }
}

void leapfrog_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par)
{
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
	
	error(int_par->nsteps==0,1,"leapfrog","Error nsteps 0\n");
	double dt = tlen / int_par->nsteps;

	lprintf("MD_INT",10+int_par->level*10,"Starting new MD trajectory using leapfrog.\n");
	lprintf("MD_INT",10+int_par->level*10,"MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",
	  int_par->level,tlen,int_par->nsteps,dt);
	
	integrator_par* ipn=int_par->next;
	if(ipn)
	{
		ipn->integrator(momenta, dt/2, ipn);
	}
	else
	{
		update_gauge_field(momenta, dt/2);
	}

	for(int n = 1; n < int_par->nsteps; n++)
	{
		gforce(momenta, dt, int_par->nmon, int_par->mon_list);
		if(ipn)
		{
			ipn->integrator(momenta, dt, ipn);
		}
		else
		{
			update_gauge_field(momenta, dt);
		}
	}

	gforce(momenta, dt, int_par->nmon, int_par->mon_list);
	if(ipn)
	{
		ipn->integrator(momenta, dt/2, ipn);
	}
	else
	{
		update_gauge_field(momenta, dt/2);
	}
}


void O2MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par)
{
	const double lambda = 0.1931833275037836;

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

	error(int_par->nsteps==0,1,"O2MN_multistep","Error nsteps 0\n");
	double dt = tlen / int_par->nsteps;

	lprintf("MD_INT",10+int_par->level*10,"Starting new MD trajectory with O2MN_multistep, level %d.\n",int_par->level);
	lprintf("MD_INT",10+int_par->level*10,"MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",
	  int_par->level,tlen,int_par->nsteps,dt);
 
	integrator_par* ipn=int_par->next;
	gforce(momenta, lambda*dt, int_par->nmon, int_par->mon_list);

	for(int n = 1; n <= int_par->nsteps; n++)
	{
		if(ipn)
		{
			ipn->integrator(momenta, dt/2., ipn);
			gforce(momenta, (1.-2.*lambda)*dt, int_par->nmon, int_par->mon_list);
			ipn->integrator(momenta, dt/2., ipn);
		}
		else
		{
			update_gauge_field(momenta, dt/2);
			gforce(momenta, (1.-2.*lambda)*dt, int_par->nmon, int_par->mon_list);
			update_gauge_field(momenta, dt/2);
		}

		/* Update of momenta */
		if(n < int_par->nsteps)
		{
			gforce(momenta, 2.*lambda*dt, int_par->nmon, int_par->mon_list);
			lprintf("MD_INT",20+int_par->level*10,"MD step (level %d): %d/%d\n",int_par->level,n,int_par->nsteps);
		}
		else
		{
			gforce(momenta, lambda*dt, int_par->nmon, int_par->mon_list);
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


/* 4th order  I.P. Omelyan, I.M. Mryglod, R. Folk, computer Physics Communications 151 (2003) 272-314 */
void O4MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par)
{
	const double r1 = 0.08398315262876693;
	const double r2 = 0.2539785108410595;
	const double r3 = 0.6822365335719091;
	const double r4 = -0.03230286765269967;

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
	
	error(int_par->nsteps==0,1,"O4MN_multistep","Error nsteps 0\n");
	double dt = tlen / int_par->nsteps;

	lprintf("MD_INT",10+int_par->level*10,"Starting new MD trajectory with O4MN_multistep, level %d.\n",int_par->level);
	lprintf("MD_INT",10+int_par->level*10,"MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",
			int_par->level,tlen,int_par->nsteps,dt);

	integrator_par* ipn=int_par->next;
	gforce(momenta, r1*dt, int_par->nmon, int_par->mon_list);

	for(int n = 1; n <= int_par->nsteps; n++)
	{
		
		if(ipn)
		{
			ipn->integrator(momenta, r2*dt, ipn);
			gforce(momenta, r3*dt, int_par->nmon, int_par->mon_list);
			ipn->integrator(momenta, r4*dt, ipn);
			gforce(momenta, (0.5-r1-r3)*dt, int_par->nmon, int_par->mon_list);
			ipn->integrator(momenta, (1.-2.*(r2+r4))*dt, ipn);
			gforce(momenta, (0.5-r1-r3)*dt, int_par->nmon, int_par->mon_list);
			ipn->integrator(momenta, r4*dt, ipn);
			gforce(momenta, r3*dt, int_par->nmon, int_par->mon_list);
			ipn->integrator(momenta, r2*dt, ipn);
		}
		else
		{
			update_gauge_field(momenta, r2*dt);
			gforce(momenta, r3*dt, int_par->nmon, int_par->mon_list);
			update_gauge_field(momenta, r4*dt);
			gforce(momenta, (0.5-r1-r3)*dt, int_par->nmon, int_par->mon_list);
			update_gauge_field(momenta, (1.-2.*(r2+r4))*dt);
			gforce(momenta, (0.5-r1-r3)*dt, int_par->nmon, int_par->mon_list);
			update_gauge_field(momenta, r4*dt);
			gforce(momenta, r3*dt, int_par->nmon, int_par->mon_list);
			update_gauge_field(momenta, r2*dt);
		}
		
		/* Update of momenta */
		if(n < int_par->nsteps)
		{
			gforce(momenta, 2*r1*dt, int_par->nmon, int_par->mon_list);
			lprintf("MD_INT",20+int_par->level*10,"MD step (level %d): %d/%d\n",int_par->level,n,int_par->nsteps);
		}
		else
		{
			gforce(momenta, r1*dt, int_par->nmon, int_par->mon_list);
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
