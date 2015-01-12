/***************************************************************************
 * Copyright (c) 2014, Martin Hansen                                       *
 * All rights reserved.                                                    *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "update.h"
#include "geometry.h"
#include "logger.h"
#include "memory.h"
#include "global.h"
#include "linear_algebra.h"

#define MAX_SECTIONS 8
#define MAX_VALUES   16
#define MAX_LENGTH   32

typedef struct {
	int num;
	char key[MAX_VALUES][MAX_LENGTH];
	char value[MAX_VALUES][MAX_LENGTH];
} section;

static integrator_par *ip = 0;
static section sec_int[MAX_SECTIONS];
static section sec_mon[MAX_SECTIONS];
static int n_int = 0;
static int n_mon = 0;
static int last_error = 0;

static void check(int test, const char* format, ...)
{
	char buf[128];
	if(test)
	{
		va_list args;
		va_start(args, format);
		vsprintf(buf, format, args);
		va_end(args);
		lprintf("ACTION", 0, buf);

		finalize_process();
		exit(1);
	}
}

static char* find_string(section *sec, char *key)
{
	for(int i = 0; i < sec->num; i++)
	{
		if(strcmp(key, sec->key[i]) == 0)
		{
			last_error = 0;
			return sec->value[i];
		}
	}

	last_error = 1;
	return 0;
}

static double find_double(section *sec, char *key)
{
	for(int i = 0; i < sec->num; i++)
	{
		if(strcmp(key, sec->key[i]) == 0)
		{
			last_error = 0;
			return atof(sec->value[i]);
		}
	}

	last_error = 1;
	return 0;
}

static void add_monomial_to_integrator(const monomial *m, int level)
{
	integrator_par *iter = ip;
	while(iter->level != level) iter = iter->next;
	iter->mon_list[iter->nmon] = m;
	iter->nmon++;
}

static void setup_monomials()
{
	spinor_field *stmp;
	monomial mon;
	const monomial *mret;
	section *cur;
	char *type;
	int level;

	#ifdef UPDATE_EO
	stmp = alloc_spinor_field_f(1, &glat_even);
	#else
	stmp = alloc_spinor_field_f(1, &glattice);
	#endif

	for(int i = 0; i < n_mon; i++)
	{
		cur = &sec_mon[i];
		memset(&mon, 0, sizeof(mon));
		mon.id = i;

		type = find_string(cur, "type");
		check(last_error, "Monomial type not specified\n");

		level = find_double(cur, "level");
		check(last_error, "Unable to find 'level' in monomial");
		check(level<0||level>=n_int, "Invalid integrator level %d in monomial\n", level);

		if(strcmp(type, "gauge") == 0)
		{
			mon_pg_par par;
			mon.par = &par;
			mon.type = PureGauge;

			// Find parameters
			par.beta = find_double(cur, "beta");
			check(last_error, "Unable to find 'beta' in monomial of type 'gauge'\n");

			// Add monomial
			mret = add_mon(&mon);

			// Monomial information
			lprintf("ACTION", 10, "Monomial %d: level = %d, type = gauge, beta = %1.6f\n", i, level, par.beta);
		}
		else if(strcmp(type, "hmc") == 0)
		{
			mon_hmc_par par;
			mon.par = &par;
			mon.type = HMC;
			par.pf = stmp;

			// Find parameters
			par.mass = find_double(cur, "mass");
			check(last_error, "Unable to find 'mass' in monomial of type 'hmc'\n");

			mon.MT_prec = find_double(cur, "mt_prec");
			check(last_error, "Unable to find 'mt_prec' in monomial of type 'hmc'\n");

			mon.force_prec = find_double(cur, "force_prec");
			check(last_error, "Unable to find 'force_prec' in monomial of type 'hmc'\n");

			par.mre_past = find_double(cur, "mre_past");
			check(last_error, "Unable to find 'mre_past' in monomial of type 'hmc'\n");

			// Setup force parameters
			par.fpar.n_pf = 1;
			par.fpar.pf = par.pf;
			par.fpar.inv_err2 = mon.force_prec;
			par.fpar.inv_err2_flt = 1e-6;
			par.fpar.mass = par.mass;
			par.fpar.b = 0;
			par.fpar.hasenbusch = 0;

			// Setup chronological inverter
			mre_init(&par.fpar.mpar, par.mre_past, mon.force_prec);

			// Add monomial
			mret = add_mon(&mon);
			
			// Monomial information
			lprintf("ACTION", 10, "Monomial %d: level = %d, type = hmc, mass = %1.6f, force_prec = %1.2e, mt_prec = %1.2e\n",
					i, level, par.mass, mon.force_prec, mon.MT_prec);
		}
		else if(strcmp(type, "hasenbusch") == 0)
		{
			mon_hasenbusch_par par;
			mon.par = &par;
			mon.type = Hasenbusch;
			par.pf = stmp;

			// Find parameters
			par.mass = find_double(cur, "mass");
			check(last_error, "Unable to find 'mass' in monomial of type 'hasenbusch'\n");

			par.dm = find_double(cur, "dm");
			check(last_error, "Unable to find 'dm' in monomial of type 'hasenbusch'\n");

			mon.MT_prec = find_double(cur, "mt_prec");
			check(last_error, "Unable to find 'mt_prec' in monomial of type 'hasenbusch'\n");

			mon.force_prec = find_double(cur, "force_prec");
			check(last_error, "Unable to find 'force_prec' in monomial of type 'hasenbusch'\n");

			par.mre_past = find_double(cur, "mre_past");
			check(last_error, "Unable to find 'mre_past' in monomial of type 'hasenbusch'\n");

			// Setup force parameters
			par.fpar.n_pf = 1;
			par.fpar.pf = par.pf;
			par.fpar.inv_err2 = mon.force_prec;
			par.fpar.inv_err2_flt = 1e-6;
			par.fpar.mass = par.mass;
			#ifdef UPDATE_EO
			par.fpar.b = (4.+par.mass+par.dm)*(4.+par.mass+par.dm)-(4.+par.mass)*(4.+par.mass);
			#else
			par.fpar.b = par.dm;
			#endif
			par.fpar.hasenbusch = 2;
			
			// Setup chronological inverter
			mre_init(&par.fpar.mpar, par.mre_past, mon.force_prec);

			// Add monomial
			mret = add_mon(&mon);

			// Monomial information
			lprintf("ACTION", 10, "Monomial %d: level = %d, type = hasenbusch, mass = %1.6f, dm = %1.6f, force_prec = %1.2e, mt_prec = %1.2e\n",
					i, level, par.mass, par.dm, mon.force_prec, mon.MT_prec);
		}
		else if(strcmp(type, "rhmc") == 0)
		{
			mon_rhmc_par par;
			mon.par = &par;
			mon.type = RHMC;
			par.pf = stmp;
			par.ratio.order = 16;

			// Find parameters
			par.mass = find_double(cur, "mass");
			check(last_error, "Unable to find 'mass' in monomial of type 'rhmc'\n");

			par.ratio.n = -1*find_double(cur, "n");
			check(last_error, "Unable to find 'n' in monomial of type 'rhmc'\n");

			par.ratio.d = find_double(cur, "d");
			check(last_error, "Unable to find 'd' in monomial of type 'rhmc'\n");

			mon.MT_prec = find_double(cur, "mt_prec");
			check(last_error, "Unable to find 'mt_prec' in monomial of type 'rhmc'\n");

			mon.MD_prec = find_double(cur, "md_prec");
			check(last_error, "Unable to find 'md_prec' in monomial of type 'rhmc'\n");

			mon.force_prec = find_double(cur, "force_prec");
			check(last_error, "Unable to find 'force_prec' in monomial of type 'rhmc'\n");

			// Setup force parameters
			par.fpar.n_pf = 1;
			par.fpar.pf = par.pf;
			par.fpar.inv_err2 = mon.force_prec;
			par.fpar.mass = par.mass;
			par.ratio.rel_error = mon.MD_prec;
			par.fpar.ratio = &par.ratio;

			// Add monomial
			mret = add_mon(&mon);

			// Monomial information
			lprintf("ACTION", 10, "Monomial %d: level = %d, type = rhmc, mass = %1.6f, force_prec = %1.2e, mt_prec = %1.2e, md_prec = %1.2e\n",
					i, level, par.mass, mon.force_prec, mon.MT_prec, mon.MD_prec);
		}
		else
		{
			check(1, "Unknown monomial type\n");
		}

		add_monomial_to_integrator(mret, level);
	}

	free_spinor_field_f(stmp);
}

static void setup_integrator()
{
	integrator_par tmp;
	integrator_par *iter = 0;
	section *cur;
	char *type;
	int map[MAX_SECTIONS];
	int id;

	// Determine the ordering of the integrator levels and check for inconsistencies
	for(int i = 0; i < n_int; i++)
	{
		map[i] = -1;
	}

	for(int i = 0; i < n_int; i++)
	{
		id = find_double(&sec_int[i], "level");

		check(last_error, "Integrator level not specified\n");
		check(id<0||id>=n_int, "Invalid integrator level %d\n", id);
		check(map[id]>=0, "Integrator level %d already exists\n", id);

		map[id] = i;
	}

	// Construct integrator structure
	for(int i = 0; i < n_int; i++)
	{
		cur = &sec_int[map[i]];

		type = find_string(cur, "type");
		check(last_error, "Unable to find 'type' in integrator\n");

		tmp.nsteps = find_double(cur, "steps");
		check(last_error, "Unable to find 'steps' in integrator\n");

		if(strcmp(type, "o2mn") == 0)
		{
			tmp.integrator = &O2MN_multistep;
			lprintf("INTEGRATOR", 10, "Level %d: type = o2mn, steps = %d\n", i, tmp.nsteps);
		}
		else if(strcmp(type, "o4mn") == 0)
		{
			tmp.integrator = &O4MN_multistep;
			lprintf("INTEGRATOR", 10, "Level %d: type = o4mn, steps = %d\n", i, tmp.nsteps);
		}
		else if(strcmp(type, "lf") == 0)
		{
			tmp.integrator = &leapfrog_multistep;
			lprintf("INTEGRATOR", 10, "Level %d: type = lf, steps = %d\n", i, tmp.nsteps);
		}
		else
		{
			check(1, "Unknown integrator type\n");
		}

		if(iter == 0)
		{
			iter = malloc(sizeof(integrator_par));
			ip = iter;
		}
		else
		{
			iter->next = malloc(sizeof(integrator_par));
			iter = iter->next;
		}

		// Fill structure
		iter->nsteps = tmp.nsteps;
		iter->integrator = tmp.integrator;
		iter->mon_list = calloc(n_mon, sizeof(iter->mon_list));
		iter->nmon = 0;
		iter->level = i;
		iter->next = 0;
	}
}

static void parse_lines(section *sec, char *lines)
{
	char key[MAX_LENGTH];
	char val[MAX_LENGTH];
	int count;
	int id;

	while(sscanf(lines, "%s = %s%n", key, val, &count) == 2)
	{
		id = sec->num++;
		strcpy(sec->key[id], key);
		strcpy(sec->value[id], val);
		lines += count;
	}
}

void read_action(char *filename, integrator_par **ip_ptr)
{
	FILE *fp;
	char *content;
	char *lines;
	char type[16];
	int fsz = 0;
	int count = 0;
	int pos = 0;

	// Open file
	fp = fopen(filename, "r");

	if(fp == NULL)
	{
		check(1, "Unable to open file\n");
	}

	// Clear memory
	memset(sec_int, 0, sizeof(sec_int));
	memset(sec_mon, 0, sizeof(sec_mon));

	// Determine file size
	fseek(fp, 0, SEEK_END);
	fsz = ftell(fp);
	rewind(fp);

	// Read file content
	content = malloc(fsz);
	lines = malloc(fsz);
	fread(content, fsz, 1, fp);
	fclose(fp);

	// Find all sections in the file
	while(pos < fsz)
	{
		if(sscanf(content + pos, "%s { %[^}]%n", type, lines, &count) == 2)
		{
			if(strcmp(type, "monomial") == 0)
			{
				parse_lines(&sec_mon[n_mon], lines);
				n_mon++;
			}

			if(strcmp(type, "integrator") == 0)
			{
				parse_lines(&sec_int[n_int], lines);
				n_int++;
			}

			pos += count + 1;
		}
		else
		{
			pos++;
		}
	}

	// Free memory
	free(content);
	free(lines);

	// Print log info
	if(n_mon < n_int)
	{
		check(1, "There should be at least as many monomials as integrator levels\n");
	}
	else
	{
		lprintf("ACTION", 10, "Found %d monomials and %d integrator levels\n", n_mon, n_int);
	}

	// Setup structures
	setup_integrator();
	setup_monomials();

	// Check that all integrator levels have monomials
	integrator_par *iter = ip;
	while(iter)
	{
		check(iter->nmon == 0, "No monomials on integrator level %d\n", iter->level);
		iter = iter->next;
	}

	// Return integrator structure
	*ip_ptr = ip;
}