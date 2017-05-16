/***************************************************************************\
 * Copyright (c) 2017                                                        *
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
	lprintf("TIMING",0,"update_gauge_field %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
#endif
}

void update_scalar_field(double dt, void *vpar)
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

	field_scalar_par *par = (field_scalar_par*)vpar;
	suNg_scalar_field *s_field = *par->field;
	suNg_scalar_field *force = *par->momenta;
	suNg_vector mom_star;

	_MASTER_FOR(&glattice,ix)
	{
		vector_star(&mom_star, _FIELD_AT(force,ix));
		_vector_mul_add_assign_g(*_FIELD_AT(s_field,ix), dt, mom_star);
	}

	start_sc_sendrecv(s_field);
	complete_sc_sendrecv(s_field);

#ifdef TIMING
#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
#endif
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"update_scalar_field %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
#endif
}

