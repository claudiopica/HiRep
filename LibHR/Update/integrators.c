/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica, Ari Hietanen          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include <stdio.h>

static double total_action(scalar_field* la){
	double tot_action=0.0;
	_MASTER_FOR(&glattice,i){
		tot_action += *_FIELD_AT(la,i);
	}

	return tot_action;	
}
void monomial_force(double dt, int nmon, const monomial **mon_list)
{
	scalar_field *la = alloc_sfield(1,&glattice); 
	local_hmc_action(NEW,la,suN_momenta,scalar_momenta);
//	printf("Force: mon type: action = %f \n",total_action(la));
	for(int n = 0; n < nmon; n++)
	{
		const monomial *m = mon_list[n];
//		local_hmc_action(NEW,la,suN_momenta,scalar_momenta);
		m->update_force(dt, m->force_par);
	}

}

void monomial_field(double dt, int nmon, const monomial **mon_list, integrator_par *ipn)
{
	scalar_field *la = alloc_sfield(1,&glattice); 
	for(int n = 0; n < nmon; n++)
	{
		const monomial *m = mon_list[n];
		if(m->update_field)
		{
//			local_hmc_action(NEW,la,suN_momenta,scalar_momenta);
			m->update_field(dt, m->field_par);
		}
	}
	local_hmc_action(NEW,la,suN_momenta,scalar_momenta);
	printf("Field: mon type: action = %f \n",total_action(la));
	if(ipn)
	{
		ipn->integrator(dt, ipn);
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
	lprintf("MD_INT", level, "MD parameters: level=%d tlen=%1.6f nsteps=%d => dt=%1.6f\n",par->level, tlen, par->nsteps, dt);
	
	for(int n = 0; n < par->nsteps; n++)
	{
		if(n == 0)
		{
			monomial_force(dt/2, par->nmon, par->mon_list);
		}
		else
		{
			monomial_force(dt, par->nmon, par->mon_list);
		}
		monomial_field(dt, par->nmon, par->mon_list, par->next);
	}
	monomial_force(dt/2, par->nmon, par->mon_list);
}


void O2MN_multistep(double tlen, integrator_par *par)
{
	const double lambda = 0.1931833275037836;
	const double dt = tlen / par->nsteps;
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
			monomial_force(lambda*dt, par->nmon, par->mon_list);
		}
		else
		{
			monomial_force(2*lambda*dt, par->nmon, par->mon_list);
		}
		monomial_field(dt/2, par->nmon, par->mon_list, par->next);
		monomial_force((1-2*lambda)*dt, par->nmon, par->mon_list);
		monomial_field(dt/2, par->nmon, par->mon_list, par->next);
		
	}
	monomial_force(lambda*dt, par->nmon, par->mon_list);
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
			monomial_force(r1*dt, par->nmon, par->mon_list);
		}
		else
		{
			monomial_force(2*r1*dt, par->nmon, par->mon_list);
		}
		monomial_field(r2*dt, par->nmon, par->mon_list, par->next);
		monomial_force(r3*dt, par->nmon, par->mon_list);
		monomial_field(r4*dt, par->nmon, par->mon_list, par->next);
		monomial_force((0.5-r1-r3)*dt, par->nmon, par->mon_list);
		monomial_field((1-2*(r2+r4))*dt, par->nmon, par->mon_list, par->next);
		monomial_force((0.5-r1-r3)*dt, par->nmon, par->mon_list);
		monomial_field(r4*dt, par->nmon, par->mon_list, par->next);
		monomial_force(r3*dt, par->nmon, par->mon_list);
		monomial_field(r2*dt, par->nmon, par->mon_list, par->next);
	}
	monomial_force(r1*dt, par->nmon, par->mon_list);
}

