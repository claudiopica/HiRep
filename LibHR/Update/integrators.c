/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica, Ari Hietanen          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"

void monomial_force(double dt, integrator_par *par)
{
	for(int n = 0; n < par->nmon; n++)
	{
		const monomial *m = par->mon_list[n];
		m->update_force(dt, m->force_par);
	}
}

void monomial_field(double dt, integrator_par *par)
{
	for(int n = 0; n < par->nmon; n++)
	{
		const monomial *m = par->mon_list[n];
		if(m->update_field)
		{
			m->update_field(dt, m->field_par);
		}
	}
	if(par->next)
	{
		par->next->integrator(dt, par->next);
	}
}

void leapfrog_multistep(double tlen, integrator_par *par)
{
	double dt = tlen / par->nsteps;
	int level = 10+par->level*10;

	if(par->nsteps == 0)
	{
		return;
	}

	lprintf("MD_INT", level, "Starting new MD trajectory with Leapfrog\n");
	lprintf("MD_INT", level, "MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n", par->level, tlen, par->nsteps, dt);
	
	for(int n = 0; n < par->nsteps; n++)
	{
		if(n == 0)
		{
			monomial_force(dt/2, par);
		}
		else
		{
			monomial_force(dt, par);
		}
		monomial_field(dt, par);
	}
	monomial_force(dt/2, par);
}


void O2MN_multistep(double tlen, integrator_par *par)
{
	const double lambda = 0.1931833275037836;
	double dt = tlen / par->nsteps;
	int level = 10+par->level*10;

	if(par->nsteps == 0)
	{
		return;
	}

	lprintf("MD_INT", level, "Starting new MD trajectory with O2MN_multistep\n");
	lprintf("MD_INT", level, "MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n", par->level, tlen, par->nsteps, dt);
 
	for(int n = 0; n < par->nsteps; n++)
	{
		if(n == 0)
		{
			monomial_force(lambda*dt, par);
		}
		else
		{
			monomial_force(2*lambda*dt, par);
		}
		monomial_field(dt/2, par);
		monomial_force((1-2*lambda)*dt, par);
		monomial_field(dt/2, par);
	}
	monomial_force(lambda*dt, par);
}

/* 4th order  I.P. Omelyan, I.M. Mryglod, R. Folk, computer Physics Communications 151 (2003) 272-314 */
void O4MN_multistep(double tlen, integrator_par *par)
{
	const double r1 = 0.08398315262876693;
	const double r2 = 0.2539785108410595;
	const double r3 = 0.6822365335719091;
	const double r4 = -0.03230286765269967;
	double dt = tlen / par->nsteps;
	int level = 10+par->level*10;

	if(par->nsteps == 0)
	{
		return;
	}

	lprintf("MD_INT", level, "Starting new MD trajectory with O4MN_multistep\n");
	lprintf("MD_INT", level, "MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n", par->level, tlen, par->nsteps, dt);

	for(int n = 0; n < par->nsteps; n++)
	{
		if(n == 0)
		{
			monomial_force(r1*dt, par);
		}
		else
		{
			monomial_force(2*r1*dt, par);
		}
		monomial_field(r2*dt, par);
		monomial_force(r3*dt, par);
		monomial_field(r4*dt, par);
		monomial_force((0.5-r1-r3)*dt, par);
		monomial_field((1-2*(r2+r4))*dt, par);
		monomial_force((0.5-r1-r3)*dt, par);
		monomial_field(r4*dt, par);
		monomial_force(r3*dt, par);
		monomial_field(r2*dt, par);
	}
	monomial_force(r1*dt, par);
}
