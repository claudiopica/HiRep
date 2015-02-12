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
#include "update.h"
#include "global.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "logger.h"
#include "reweight.h"

// Variables
static char  input_filename[64] = "input_file";
static char output_filename[64] = "out_reweight";
static char config_filename[64] = "";

static double mass = 0;
static int cnfg_number = 0;

// Mass parameters
static double mass_old = 0;
static double mass_new = 0;

// Theta angle
static double theta_old[4] = {0.0, 0.0, 0.0, 0.0};
static double theta_new[4] = {0.0, 0.0, 0.0, 0.0};

// Wrappers
void Dirac_operator(spinor_field *out, spinor_field *in)
{
	Dphi_eopre(mass, out, in);
}

void load_cnfg(char *filename)
{
	int len;
	char *basename;
	int x, y, z, t, nc;
	char repr[4];

	len = strlen(filename);
	basename  = strrchr(filename, '/');
	filename[len - 1] = '\0';

	if(basename == NULL)
	{
		basename = filename;
	}
	else
	{
		basename++;
	}

	sscanf(basename, "%*[^_]_%dx%dx%dx%dnc%dr%3snf%*db%*fm%lfn%d", &t, &x, &y, &z, &nc, repr, &mass_old, &cnfg_number);
	mass_old = -mass_old;
	mass = mass_old;

	if(t != GLB_T || x != GLB_X || y != GLB_Y || z != GLB_Z)
	{
		lprintf("WARNING", 0, "Size read from config name (%d,%d,%d,%d) is different from the lattice size!\n", t, x, y, z);
	}

	if(nc != NG)
	{
		lprintf("WARNING", 0, "Gauge group read from config name (NG=%d) is not the one used in this code!\n", nc);
	}

	if(strcmp(repr, repr_name) != 0)
	{
		lprintf("WARNING", 0, "Representation (repr=%s) is not the one used in this code!\n", repr);
	}

	read_gauge_field(filename);
	represent_gauge_field();
	apply_BCs_on_fundamental_gauge_field();
}

void set_theta(double *a)
{
	#ifdef FERMION_THETA
	for(int i = 0; i < 4; i++)
	{
		eitheta[i].re = cos(a[i]);
		eitheta[i].im = sin(a[i]);
	}
	#endif
}

double apply_operator(double *old, double *new, spinor_field *in)
{
	double shift, val;
	spinor_field *tmp, *res;

	// Inverter
	mshift_par mpar;
	mpar.err2 = rw_var.precision;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;

	// Allocate spinor fields
	tmp = alloc_spinor_field_f(1, in->type);
	res = alloc_spinor_field_f(1, in->type);
	spinor_field_zero_f(tmp);

	// Apply dirac operators
	#ifdef REWEIGHT_THETA
	set_theta(old);
	spinor_field_g5_assign_f(in);
	g5QMR_mshift(&mpar, &Dirac_operator, in, tmp);
	spinor_field_g5_assign_f(in);

	set_theta(new);
	Dirac_operator(res, tmp);
	spinor_field_g5_assign_f(res);
	val = spinor_field_prod_re_f(res, res);
	#endif

	#ifdef REWEIGHT_MASS
	mass = old[0];
	g5QMR_mshift(&mpar, &Dirac_operator, in, tmp);

	mass = new[0];
	Dirac_operator(res, tmp);
	val = spinor_field_prod_re_f(res, res);
	#endif

	// Free spinor fields
	free_spinor_field_f(tmp);
	free_spinor_field_f(res);

	// Return
	return val;
}

void reweight(double steps, double *result)
{
	double denominator;
	double numerator;
	double old[4];
	double new[4];
	spinor_field *eta;

	// Allocate spinor
	eta = alloc_spinor_field_f(1, &glat_even);

	// Perform calculation
	for(int n = 0; n < steps; n++)
	{
		#ifdef REWEIGHT_THETA
		for(int i = 0; i < 4; i++)
		{
			old[i] = theta_old[i] + (theta_new[i] - theta_old[i]) * (n / steps);
			new[i] = theta_old[i] + (theta_new[i] - theta_old[i]) * ((n + 1) / steps);
		}
		#endif

		#ifdef REWEIGHT_MASS
		old[0] = mass_old + (mass_new - mass_old) * (n / steps);
		new[0] = mass_old + (mass_new - mass_old) * ((n + 1) / steps);
		#endif

		gaussian_spinor_field(eta);
		denominator = spinor_field_prod_re_f(eta, eta);
		numerator = apply_operator(old, new, eta);
		result[n] = (denominator - numerator);
	}

	// Free spinor
	free_spinor_field_f(eta);
}

void read_cmdline(int argc, char *argv[])
{
	int op = 0;
	int ip = 0;
	int cp = 0;

	for(int i = 1; i < argc; i++)
	{
		if(strcmp("-o", argv[i]) == 0)
		{
			op = ++i;
			continue;
		}
		if(strcmp("-i", argv[i]) == 0)
		{
			ip = ++i;
			continue;
		}
		if(strcmp("-l", argv[i]) == 0)
		{
			cp = ++i;
			continue;
		}
	}

	if(op > 0 && op < argc)
	{
		strcpy(output_filename, argv[op]);
	}

	if(ip > 0 && ip < argc)
	{
		strcpy(input_filename, argv[ip]);
	}

	if(cp > 0 && cp < argc)
	{
		strcpy(config_filename, argv[cp]);
	}
}

int main(int argc, char *argv[])
{
	char sbuf[512];
	char rbuf[512];

	// Read commandline
	read_cmdline(argc, argv);

	// Setup process
	setup_process(&argc, &argv);
	logger_setlevel(0, 10);

	if(PID != 0)
	{
		logger_disable();
	}
	else
	{
		sprintf(sbuf, ">>%s", output_filename);
		logger_stdout(sbuf);
		sprintf(sbuf, "err_%d", PID);
		freopen(sbuf, "w", stderr);
	}

	// Read settings
	read_input(glb_var.read, input_filename);
	read_input(rlx_var.read, input_filename);
	read_input(rw_var.read,  input_filename);

	// Initialize stuff
	rlxd_init(rlx_var.rlxd_level, rlx_var.rlxd_seed+PID);

	if(geometry_init() == 1)
	{
		finalize_process();
		return 0;
	}

	geometry_mpi_eo();

	// Allocate gauge field
	u_gauge = alloc_gfield(&glattice);
	u_gauge_f = alloc_gfield_f(&glattice);
	unit_u(u_gauge);
	represent_gauge_field();

	///////////////////////////////////////////////////////////////////////////////////////////
	// START OF ACTUAL COMPUTATION                                                           //
	///////////////////////////////////////////////////////////////////////////////////////////

	// Handle input parameters
	theta_old[0] = rw_var.old_theta_t / GLB_T;
	theta_old[1] = rw_var.old_theta_x / GLB_X;
	theta_old[2] = rw_var.old_theta_y / GLB_Y;
	theta_old[3] = rw_var.old_theta_z / GLB_Z;

	theta_new[0] = rw_var.new_theta_t / GLB_T;
	theta_new[1] = rw_var.new_theta_x / GLB_X;
	theta_new[2] = rw_var.new_theta_y / GLB_Y;
	theta_new[3] = rw_var.new_theta_z / GLB_Z;

	mass_old = rw_var.old_mass;
	mass_new = rw_var.new_mass;

	// Print log info
	lprintf("REWEIGHT", 10, "Reweightning parameters\n");
	lprintf("REWEIGHT", 10, "hits = %d\n", rw_var.hits);
	lprintf("REWEIGHT", 10, "substeps = %d\n", rw_var.steps);
	lprintf("REWEIGHT", 10, "inverter precision = %1.4e\n", rw_var.precision);

	#ifdef REWEIGHT_MASS
	lprintf("REWEIGHT", 10, "old mass = %1.6f\n", rw_var.old_mass);
	lprintf("REWEIGHT", 10, "new mass = %1.6f\n", rw_var.new_mass);
	#endif

	#ifdef REWEIGHT_THETA
	lprintf("REWEIGHT", 10, "old twisting angles = %1.4f %1.4f %1.4f %1.4f\n", rw_var.old_theta_t, rw_var.old_theta_x, rw_var.old_theta_y, rw_var.old_theta_z);
	lprintf("REWEIGHT", 10, "new twisting angles = %1.4f %1.4f %1.4f %1.4f\n", rw_var.new_theta_t, rw_var.new_theta_x, rw_var.new_theta_y, rw_var.new_theta_z);
	#endif

	// Variables
	FILE *fp;
	char filename[512];
	double result[rw_var.steps];
	double summed[rw_var.steps];
	double weight;

	// Open config file
	fp = fopen(config_filename, "r");

	if(fp == NULL)
	{
		lprintf("WARNING", 0, "Unable to open config file\n");
		exit(0);
	}

	while(fgets(filename, 512, fp) != NULL)
	{
		load_cnfg(filename);
		memset(summed, 0, sizeof(summed));

		for(int i = 0; i < rw_var.hits; i++)
		{
			reweight(rw_var.steps, result);
			memset(rbuf, 0, sizeof(rbuf));
			weight = 1;

			for(int k = 0; k < rw_var.steps; k++)
			{
				summed[k] += exp(result[k]);
				sprintf(sbuf, "%1.6e ", result[k]);
				strcat(rbuf, sbuf);
				weight *= summed[k] / (i+1);
			}

			lprintf("REWEIGHT", 10, "cnfg: %d, hit: %d, weight: %1.6e, exponents: %s\n", cnfg_number, i+1, sqrt(1.0/weight), rbuf);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	// END OF ACTUAL COMPUTATION                                                             //
	///////////////////////////////////////////////////////////////////////////////////////////

	free_gfield(u_gauge);
	free_gfield_f(u_gauge_f);
	finalize_process();
}
