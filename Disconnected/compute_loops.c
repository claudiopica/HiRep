/*******************************************************************************
 *
 * Compute some disconnected loops
 * Copyright (c) 2014, R. Arthur, V. Drach, A. Hietanen 
 * All rights reserved.
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "disconnected.h"
#include "clover_tools.h"
#include "setup.h"
#include "data_storage.h"

#define PI 3.141592653589793238462643383279502884197

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif
/* Disonnected parameters */
typedef struct _input_loops
{
	char mstring[256];
	double precision;
	int nhits;
	int source_type;
	int n_mom;
	double csw;
	char configlist[256];
	/* for the reading function */
	input_record_t read[8];

} input_loops;

#define init_input_loops(varname)                                                             \
	{                                                                                         \
		.read = {                                                                             \
			{"Fermion mass", "disc:mass = %s", STRING_T, (varname).mstring},                  \
			{"inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision},   \
			{"number of inversions per cnfg", "disc:nhits = %d", INT_T, &(varname).nhits},    \
			{"Source type ", "disc:source_type = %d", INT_T, &(varname).source_type},         \
			{"maximum component of momentum", "disc:n_mom = %d", INT_T, &(varname).n_mom},    \
			{"Configuration list:", "disc:configlist = %s", STRING_T, &(varname).configlist}, \
			{"csw", "disc:csw = %f", DOUBLE_T, &(varname).csw},                               \
			{NULL, NULL, INT_T, NULL}                                                         \
		}                                                                                     \
	}

char input_filename[256] = "input_file";
int Nsource;
double M;

enum
{
	UNKNOWN_CNFG,
	DYNAMICAL_CNFG,
	QUENCHED_CNFG
};

input_loops disc_var = init_input_loops(disc_var);

typedef struct
{
	char string[256];
	int t, x, y, z;
	int nc, nf;
	double b, m;
	int n;
	int type;
} filename_t;

int main(int argc, char *argv[])
{
	int i;
	FILE *list;

	double m[256];
	char list_filename[256] = "";
	char cnfg_filename[256] = "";

	struct timeval start, end, etime;
	gettimeofday(&start, 0);
	/* setup process id and communications */

	setup_process(&argc, &argv);
	setup_gauge_fields();

	read_input(glb_var.read, get_input_filename());
	read_input(disc_var.read, get_input_filename());
	read_input(rlx_var.read, get_input_filename());
	strcpy(list_filename, disc_var.configlist);
	lprintf("MAIN", 0, "list_filename = %s \n", list_filename, disc_var.configlist);
	if (strcmp(list_filename, "") != 0)
	{
		error((list = fopen(list_filename, "r")) == NULL, 1, "main [measure_spectrum.c]",
			  "Failed to open list file\n");
	}
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
	set_csw(&disc_var.csw);
#endif

	init_BCs(NULL);

	m[0] = atof(disc_var.mstring); //

	lprintf("MAIN", 0, "Inverter precision = %e\n", disc_var.precision);
	lprintf("MAIN", 0, "Mass[%d] = %f\n", 0, m[0]);
	lprintf("MAIN", 0, "nhits = %d\n", 0, disc_var.nhits);
	i = 0;
	while (++i)
	{

		if (list != NULL)
			if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
				break;

		lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
		/* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
		read_gauge_field(cnfg_filename);
		represent_gauge_field();

		lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
		full_plaquette();

		lprintf("CORR", 0, "Number of noise vector : nhits = %i \n", disc_var.nhits);

		measure_loops(m, disc_var.nhits, i, disc_var.precision, disc_var.source_type, disc_var.n_mom, DONTSTORE, NULL);

		if (list == NULL)
			break;
	}

	if (list != NULL)
		fclose(list);

	finalize_process();

	free_BCs();

	free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
	free_gfield_f(u_gauge_f);
#endif

	/* close communications */
	finalize_process();

	gettimeofday(&end, 0);
	timeval_subtract(&etime, &end, &start);
	lprintf("TIMING", 0, "Inversions and contractions for configuration  [%s] done [%ld sec %ld usec]\n", cnfg_filename, etime.tv_sec, etime.tv_usec);

	return 0;
}
