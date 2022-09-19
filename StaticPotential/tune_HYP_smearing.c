/*******************************************************************************
 *
 * Search the best parameters for the HYP smearing
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

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional
#endif

#ifdef BC_XYZ_TWISTED
#error This code does not work with the twisted BCs
#endif

char cnfg_filename[256] = "";
char list_filename[256] = "";
char input_filename[256] = "input_file";
char output_filename[256] = "HYPsmearing.out";
enum
{
  UNKNOWN_CNFG,
  DYNAMICAL_CNFG,
  QUENCHED_CNFG
};

typedef struct
{
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;

int parse_cnfg_filename(char *filename, filename_t *fn)
{
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL)
  {
    basename = tmp + 1;
  }

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif
  hm = sscanf(basename, "%*[^_]_%dx%dx%dx%d%*[Nn]c%dr" repr_name "%*[Nn]f%db%lfm%lfn%d",
              &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->nf), &(fn->b), &(fn->m), &(fn->n));
  if (hm == 9)
  {
    fn->m = -fn->m; /* invert sign of mass */
    fn->type = DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }
#undef repr_name

  double kappa;
  hm = sscanf(basename, "%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
              &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->nf), &(fn->b), &kappa, &(fn->n));
  if (hm == 9)
  {
    fn->m = .5 / kappa - 4.;
    fn->type = DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }

  hm = sscanf(basename, "%dx%dx%dx%d%*[Nn]c%db%lfn%d",
              &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->b), &(fn->n));
  if (hm == 7)
  {
    fn->type = QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  hm = sscanf(basename, "%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
              &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->b), &(fn->n));
  if (hm == 7)
  {
    fn->type = QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  fn->type = UNKNOWN_CNFG;
  return UNKNOWN_CNFG;
}

void read_cmdline(int argc, char *argv[])
{
  int i, ai = 0, ao = 0, ac = 0, al = 0, am = 0;
  FILE *list = NULL;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-i") == 0)
      ai = i + 1;
    else if (strcmp(argv[i], "-o") == 0)
      ao = i + 1;
    else if (strcmp(argv[i], "-c") == 0)
      ac = i + 1;
    else if (strcmp(argv[i], "-l") == 0)
      al = i + 1;
    else if (strcmp(argv[i], "-m") == 0)
      am = i;
  }

  if (am != 0)
  {
    print_compiling_info();
    exit(0);
  }

  if (ao != 0)
    strcpy(output_filename, argv[ao]);
  if (ai != 0)
    strcpy(input_filename, argv[ai]);

  error((ac == 0 && al == 0) || (ac != 0 && al != 0), 1, "parse_cmdline [tune_HYP_smearing.c]",
        "Syntax: tune_HYP_smearing { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

  if (ac != 0)
  {
    strcpy(cnfg_filename, argv[ac]);
    strcpy(list_filename, "");
  }
  else if (al != 0)
  {
    strcpy(list_filename, argv[al]);
    error((list = fopen(list_filename, "r")) == NULL, 1, "parse_cmdline [tune_HYP_smearing.c]",
          "Failed to open list file\n");
    error(fscanf(list, "%s", cnfg_filename) == 0, 1, "parse_cmdline [tune_HYP_smearing.c]",
          "Empty list file\n");
    fclose(list);
  }
}

int main(int argc, char *argv[])
{
  int i;
  char tmp[256];
  FILE *list;
  filename_t fpars;
  double mtp[6859], mtp0, mtp1;
  double weight[3];
  int ibest;

  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc, &argv);

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0, 10);
  if (PID != 0)
  {
    logger_disable();
  }
  if (PID == 0)
  {
    sprintf(tmp, ">%s", output_filename);
    logger_stdout(tmp);
    sprintf(tmp, "err_%d", PID);
    freopen(tmp, "w", stderr);
  }

  lprintf("MAIN", 0, "Compiled with macros: %s\n", MACROS);
  lprintf("MAIN", 0, "PId =  %d [world_size: %d]\n\n", PID, WORLD_SIZE);
  lprintf("MAIN", 0, "input file [%s]\n", input_filename);
  lprintf("MAIN", 0, "output file [%s]\n", output_filename);
  if (list_filename != NULL)
    lprintf("MAIN", 0, "list file [%s]\n", list_filename);
  else
    lprintf("MAIN", 0, "cnfg file [%s]\n", cnfg_filename);

  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename, &fpars);

  read_input(glb_var.read, input_filename);
  GLB_T = fpars.t;
  GLB_X = fpars.x;
  GLB_Y = fpars.y;
  GLB_Z = fpars.z;
  error(fpars.type == UNKNOWN_CNFG, 1, "tune_HYP_smearing.c", "Bad name for a configuration file");
  error(fpars.nc != NG, 1, "tune_HYP_smearing.c", "Bad NG");

  lprintf("MAIN", 0, "Gauge group: SU(%d)\n", NG);
  lprintf("MAIN", 0, "Fermion representation: " REPR_NAME " [dim=%d]\n", NF);

  /* setup communication geometry */
  if (geometry_init() == 1)
  {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* setup random numbers */
  read_input(rlx_var.read, input_filename);
  lprintf("MAIN", 0, "RLXD [%d,%d]\n", rlx_var.rlxd_level, rlx_var.rlxd_seed + MPI_PID);
  rlxd_init(rlx_var.rlxd_level, rlx_var.rlxd_seed);

  init_BCs(NULL);

  /* alloc global gauge fields */
  u_gauge = alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f = alloc_gfield_f(&glattice);
#endif

  list = NULL;
  if (strcmp(list_filename, "") != 0)
  {
    error((list = fopen(list_filename, "r")) == NULL, 1, "main [mk_mesons.c]",
          "Failed to open list file\n");
  }

  for (i = 0; i < 6859; i++)
    mtp[i] = 0.;
  mtp0 = 0.;

  i = 0;
  while (1)
  {

    if (list != NULL)
      if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
        break;

    i++;

    lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());

    HYP_span_parameters(mtp);
    mtp0 += min_tplaq(u_gauge);

    ibest = HYP_best_parameters(mtp, weight);
    lprintf("MAIN", 0, "Best weights for HYP smearing: %.2f %.2f %.2f\n", weight[0], weight[1], weight[2]);
    lprintf("MAIN", 0, "The smallest plaquette has been increased from %e to %e\n", mtp0 / i, mtp[ibest] / i);

    if (list == NULL)
      break;
  }

  if (list != NULL)
    fclose(list);

  ibest = HYP_best_parameters(mtp, weight);
  mtp0 /= i;
  mtp1 = mtp[ibest] / i;

  lprintf("MAIN", 0, "Best weights for HYP smearing: %.2f %.2f %.2f\n", weight[0], weight[1], weight[2]);
  lprintf("MAIN", 0, "The smallest plaquette has been increased from %e to %e\n", mtp0, mtp1);

  free_BCs();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}
