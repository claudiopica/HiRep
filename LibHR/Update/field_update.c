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

// Project gauge field every 2^_PROJ_BIT changes
static unsigned int count = 0;
#define _PROJ_BIT (1<<4)

void update_gauge_field(double dt, void *vpar)
{
#ifdef TIMING
	struct timeval start, end;
	struct timeval start1, end1;
	struct timeval etime;

#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
#endif
	gettimeofday(&start,0);
#endif

	field_gauge_par *par = (field_gauge_par*)vpar;
	suNg_field *gfield = *par->field;
	suNg_av_field *force = *par->momenta;

	_MASTER_FOR(&glattice,ix)
	{
		ExpX(dt, _4FIELD_AT(force,ix,0), _4FIELD_AT(gfield,ix,0));
		ExpX(dt, _4FIELD_AT(force,ix,1), _4FIELD_AT(gfield,ix,1));
		ExpX(dt, _4FIELD_AT(force,ix,2), _4FIELD_AT(gfield,ix,2));
		ExpX(dt, _4FIELD_AT(force,ix,3), _4FIELD_AT(gfield,ix,3));
	}

	if(count & _PROJ_BIT)
	{
		count = 0;
		project_gauge_field();
		complete_gf_sendrecv(u_gauge);
	}
	else
	{
		count++;
		start_gf_sendrecv(u_gauge);
		complete_gf_sendrecv(u_gauge);
	}

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
