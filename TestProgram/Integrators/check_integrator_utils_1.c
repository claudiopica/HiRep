/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "logger.h"
#include "check_integrator_1.h"
#include "random.h"
#include "io.h"
#include "representation.h"
#include "update.h"
#include "memory.h"
#include "utils.h"
#include "communications.h"
#include "clover_tools.h"
#include "linear_algebra.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

input_hmc hmc_var = init_input_hmc(hmc_var);

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif

/* add dirname to filename and return it */
static char *add_dirname(char *dirname, char *filename)
{
	static char buf[256];
	strcpy(buf, dirname);
	return strcat(buf, filename);
}

/* convert string to lowercase */
static void slower(char *str)
{
	while (*str)
	{
		*str = (char)(tolower(*str));
		++str;
	}
}

/* read g_start string and decide what the initial config should be.
 * It also fill the content of the variable start in hmc_flow with the
 * id of the first config to be used.
 * returns:
 * 0 => g_start is a valid file name which should be used as starting config
 * 1 => use unit gauge config
 * 2 => use a random config
 */
static int parse_gstart(hmc_flow *rf)
{
	int ret = 0;
	int len = 0;
	char buf[256];
	char *ptr;

	ptr = strrchr(rf->g_start, 'n');
	ret = (int)(ptr - rf->g_start);
	len = strlen(rf->run_name);
	if (ptr && (ret == len || strncmp(rf->g_start, rf->run_name, len) == 0))
	{
		rf->start = atoi(ptr + 1) + 1;
		return 0;
	}

	rf->start = 1; /* reset rf->start */

	/* try other matches */
	strcpy(buf, rf->g_start);
	slower(buf);

	ret = strcmp(buf, "unit");
	if (ret == 0)
	{
		lprintf("FLOW", 0, "Starting a new run from a unit conf!\n");
		return 1;
	}

	ret = strcmp(buf, "random");
	if (ret == 0)
	{
		lprintf("FLOW", 0, "Starting a new run from a random conf!\n");
		return 2;
	}

	lprintf("ERROR", 0, "Invalid starting gauge conf specified [%s]\n", rf->g_start);
	error(1, 1, "parse_gstart " __FILE__, "invalid config name");

	return -1;
}

/* Initialize the Monte Carlo.
 * This performs the following operations:
 * 1) read from the specified input file the flow variables 
 *    and the hmc parameters;
 * 2) set the starting gauge field
 * 3) init the hmc update
 */
int init_mc_ghmc(hmc_flow *rf, char *ifile)
{

	int start_t;

	/* flow defaults */
	strcpy(rf->g_start, "invalid");
	strcpy(rf->run_name, "run_name");
	strcpy(rf->last_conf, "invalid");
	strcpy(rf->conf_dir, "./");
	rf->save_freq = 0;
	rf->meas_freq = 0;
	rf->hmc_v = &hmc_var;

	read_input(hmc_var.read, ifile);
	read_input(rf->read, ifile);

	/* Read the action and initialize fields */
	read_action(ifile, &hmc_var.hmc_p.integrator);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
	set_csw(&hmc_var.hmc_p.csw);
#endif

	/* initialize boundary conditions */
	BCs_pars_t BCs_pars = {
		.fermion_twisting_theta = {0., 0., 0., 0.},
		.gauge_boundary_improvement_cs = 1.,
		.gauge_boundary_improvement_ct = 1.,
		.chiSF_boundary_improvement_ds = 1.,
		.SF_BCs = 0};
#ifdef FERMION_THETA
	BCs_pars.fermion_twisting_theta[0] = hmc_var.hmc_p.theta[0];
	BCs_pars.fermion_twisting_theta[1] = hmc_var.hmc_p.theta[1];
	BCs_pars.fermion_twisting_theta[2] = hmc_var.hmc_p.theta[2];
	BCs_pars.fermion_twisting_theta[3] = hmc_var.hmc_p.theta[3];
#endif
#ifdef ROTATED_SF
	BCs_pars.chiSF_boundary_improvement_ds = hmc_var.hmc_p.SF_ds;
#endif
#if defined(BASIC_SF) || defined(ROTATED_SF)
	BCs_pars.gauge_boundary_improvement_ct = hmc_var.hmc_p.SF_ct;
	error(hmc_var.hmc_p.SF_background != 0 && hmc_var.hmc_p.SF_background != 1, 0, "init_mc_ghmc" __FILE__, "Wrong value of SF_background\n");
	BCs_pars.SF_BCs = hmc_var.hmc_p.SF_background;
#endif

	init_BCs(&BCs_pars);

	/* fix conf_dir name: put a / at the end of string */
	start_t = strlen(rf->conf_dir);
	if (rf->conf_dir[start_t - 1] != '/')
	{
		rf->conf_dir[start_t] = '/';
		rf->conf_dir[start_t + 1] = '\0';
	}

	/* set initial configuration */
	start_t = parse_gstart(rf);

	/* init gauge field */
	switch (start_t)
	{
	case 0:
		read_gauge_field(add_dirname(rf->conf_dir, rf->g_start));
		if (u_scalar)
		{
			char configname[256] = "scalar_";
			strcat(configname, rf->g_start);
			read_scalar_field(add_dirname(rf->conf_dir, configname));
		}
		break;
	case 1:
		unit_u(u_gauge);
		if (u_scalar)
		{
			zero_s(u_scalar);
		}
		break;
	case 2:
		random_u(u_gauge);
		if (u_scalar)
		{
			random_s(u_scalar);
		}

		break;
	}

	start_gf_sendrecv(u_gauge);

	apply_BCs_on_fundamental_gauge_field();

	represent_gauge_field();

	/* init HMC */
	init_ghmc(&hmc_var.hmc_p);

	return 0;
}

/* clean up memory */
int end_mc()
{
	free_ghmc();
	free_BCs();

	/* free memory */
	free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
	free_gfield_f(u_gauge_f);
#endif
	if (u_scalar != NULL)
	{
		free_scalar_field(u_scalar);
	}

	return 0;
}
static void suNg_av_field_copy(suNg_av_field *g1, suNg_av_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
	_TWO_SPINORS_MATCHING(g1, g2);
#endif
	memcpy(g1->ptr, g2->ptr, 4 * g1->type->gsize_gauge * sizeof(*(g1->ptr)));
}

double integrate_ghmc(int regenerate, ghmc_par *update_par)
{
	static suNg_av_field *suN_momenta_copy;
	static double deltaH;
	static suNg_field *u_gauge_copy = NULL;
	static scalar_field *la = NULL; /* local action field for Metropolis test */
	static suNg_scalar_field *u_scalar_copy = NULL;
	static spinor_field **pf_copy = NULL;

	if (suN_momenta_copy == NULL)
	{
		suN_momenta_copy = alloc_avfield(&glattice);

		u_gauge_copy = alloc_gfield(&glattice);
		suNg_field_copy(u_gauge_copy, u_gauge);

		if (u_scalar != NULL)
		{
			u_scalar_copy = alloc_scalar_field(&glattice);
			suNg_scalar_field_copy(u_scalar_copy, u_scalar);
		}
		pf_copy = (spinor_field **)malloc(num_mon() * sizeof(spinor_field *));
		for (int i = 0; i < num_mon(); i++)
			pf_copy[i] = NULL;
	}

	if (la == NULL)
		la = alloc_sfield(1, &glattice);

	if (regenerate == 0)
	{
		/* generate new momenta */
		lprintf("HMC", 30, "Generating gaussian momenta and pseudofermions...\n");

		gaussian_momenta(suN_momenta);
		suNg_av_field_copy(suN_momenta_copy, suN_momenta);

		if (u_scalar != NULL)
		{
			gaussian_scalar_momenta(scalar_momenta);
			//suNg_av_field_copy(scalar_momenta_copy, scalar_momenta);
		}
	}

	/* generate new pseudofermions */
	spinor_field *msf;
	if (regenerate == 0)
	{
		for (int i = 0; i < num_mon(); ++i)
		{
			const monomial *m = mon_n(i);
			msf = (spinor_field *)m->pseudofermion(m);

			m->gaussian_pf(m);

			if (msf != NULL)
			{
				if (pf_copy[i] == NULL)
					pf_copy[i] = alloc_spinor_field_f(1, &glattice);

				spinor_field_copy_f(pf_copy[i], msf);
			}
		}
	}

	/* compute starting action */
	lprintf("HMC", 30, "Computing action density...\n");
	local_hmc_action(NEW, la, suN_momenta, scalar_momenta);

	/* correct pseudofermion distribution */
	for (int i = 0; i < num_mon(); ++i)
	{
		const monomial *m = mon_n(i);
		m->correct_pf(m);
	}

	/* integrate molecular dynamics */
	lprintf("HMC", 30, "MD integration...\n");
	update_par->integrator->integrator(update_par->tlen, update_par->integrator);

	/* project and represent gauge field */
	project_gauge_field();
	represent_gauge_field();

	/* compute new action */
	lprintf("HMC", 30, "Computing new action density...\n");
	for (int i = 0; i < num_mon(); ++i)
	{
		const monomial *m = mon_n(i);
		m->correct_la_pf(m);
	}
	local_hmc_action(DELTA, la, suN_momenta, scalar_momenta);

	/* Metropolis test */
	_OMP_PRAGMA(single)
	{
		deltaH = 0.0;
	}
	_MASTER_FOR_SUM(la->type, i, deltaH)
	{
		deltaH += *_FIELD_AT(la, i);
	}

	global_sum(&deltaH, 1);

	//Restore the copies

	suNg_av_field_copy(suN_momenta, suN_momenta_copy);
	//we need to wreite the function for the scalar momenta Antonio
	suNg_field_copy(u_gauge, u_gauge_copy);
	if (u_scalar != NULL)
	{
		suNg_scalar_field_copy(u_scalar, u_scalar_copy);
	}

	for (int i = 0; i < num_mon(); ++i)
	{
		const monomial *m = mon_n(i);
		msf = (spinor_field *)m->pseudofermion(m);

		m->gaussian_pf(m);
		if (msf != NULL)
		{
			error(pf_copy[i] == NULL, 0, "integrate_ghmc", "The pseudo fermion must be allocated to be copied, wrong sequence in regenarate");

			spinor_field_copy_f(msf, pf_copy[i]);
		}
	}

	return deltaH;
}

void set_integrator_nsteps(ghmc_par *gpar, int nsteps)
{
	lprintf("set_integrator_nsteps", 0, "Setting all the integrator nsteps to %d\n", nsteps);
	integrator_par *pint = gpar->integrator;
	do
	{
		pint->nsteps = nsteps;
		pint = pint->next;
	} while (pint != NULL);
}

void set_integrator_type(ghmc_par *gpar, int type)
{
	if (type == LEAPFROG)
		lprintf("set_integrator_type", 0, "Setting all the integrator type to leapfrog");
	else if (type == O2MN)
		lprintf("set_integrator_type", 0, "Setting all the integrator type to O2MN");
	else if (type == O4MN)
		lprintf("set_integrator_type", 0, "Setting all the integrator type to O4MN");
	else
		error(0 == 0, 0, "set_integrator_type", "Wrong integrator identifier");

	integrator_par *pint = gpar->integrator;
	do
	{
		if (type == LEAPFROG)
			pint->integrator = &leapfrog_multistep;
		else if (type == O2MN)
			pint->integrator = &O2MN_multistep;
		else if (type == O4MN)
			pint->integrator = &O4MN_multistep;

		pint = pint->next;
	} while (pint != NULL);
}