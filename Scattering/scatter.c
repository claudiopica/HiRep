/******************************************************************************
 *
 *
 * File scatter.c
 * This code  contains all the contractions necessary for rho to pi pi scattering
 * Checks of the rho to pi pi calculations (free case) with point source?
 *
 * Author: Tadeusz Janowski
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
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
#include "clover_tools.h"
#include "setup.h"

#include "cinfo.c"
//#include "IOroutines.c"
//#include "scatter_functions.h"

#define PI M_PI

/**
 * @brief Structure containing data from the input file relevant to scattering.
 */
typedef struct _input_scatt
{
  char mstring[256];
  double csw;
  double precision;
  int nhits;
  int tsrc;
  char outdir[256], bc[16], p[256], configlist[256];

  /* for the reading function */
  input_record_t read[11];

} input_scatt;

#define init_input_scatt(varname)                                                      \
  {                                                                                    \
    .read = {                                                                          \
      {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},       \
      {"csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw},                              \
      {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},   \
      {"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},    \
      {"Source time:", "mes:tsrc = %d", INT_T, &(varname).tsrc},                       \
      {"Output directory:", "mes:outdir = %s", STRING_T, &(varname).outdir},           \
      {"Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist}, \
      {"Boundary conditions:", "mes:bc = %s", STRING_T, &(varname).bc},                \
      {"Momenta:", "mes:p = %s", STRING_T, &(varname).p},                              \
      {NULL, NULL, INT_T, NULL}                                                        \
    }                                                                                  \
  }

char cnfg_filename[256] = "";
char list_filename[256] = "";
char prop_filename[256] = "";
char source_filename[256] = "";
char input_filename[256] = "input_file";
char output_filename[256] = "meson_scattering.out";
int Nsource;
double M;

enum
{
  UNKNOWN_CNFG,
  DYNAMICAL_CNFG,
  QUENCHED_CNFG
};

input_scatt mes_var = init_input_scatt(mes_var);

typedef struct
{
  char string[256];
  char configlist[256];
  char outdir[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;

int main(int argc, char *argv[])
{
  FILE *list = NULL;
  int tau = 0;
  double m[256];

  setup_process(&argc, &argv);

  setup_gauge_fields();

  read_input(glb_var.read, get_input_filename());
  read_input(mes_var.read, get_input_filename());
  read_input(rlx_var.read, get_input_filename());

#ifdef WITH_CLOVER
  clover_init(mes_var.csw);
#endif

  m[0] = atof(mes_var.mstring);
  init_propagator_eo(1, m, mes_var.precision);
  strcpy(list_filename, mes_var.configlist);
  lprintf("MAIN", 0, "outdir %s\n", mes_var.outdir);
  lprintf("MAIN", 0, "%s %s\n", list_filename, mes_var.configlist);

  if (strcmp(list_filename, "") != 0)
  {
    error((list = fopen(list_filename, "r")) == NULL, 1, "main [mk_mesons.c]",
          "Failed to open list file\n");
  }

  int numsources = mes_var.nhits;
  char path[256];
  strcpy(path, mes_var.outdir);
  int Nmom;
  int **p = getmomlist(mes_var.p, &Nmom);

  lprintf("MAIN", 0, "Boundary conditions: %s\n", mes_var.bc);
  lprintf("MAIN", 0, "The momenta are: %s\n", mes_var.p);
  lprintf("MAIN", 0, "mass is : %s\n", mes_var.mstring);
  lprintf("MAIN", 0, "Number of momenta: %d\n", Nmom);
  lprintf("MAIN", 0, "The momenta are:\n");

  for (int i = 0; i < Nmom; i++)
  {
    lprintf("MAIN", 0, "p%d = (%d, %d, %d)\n", i + 1, p[i][0], p[i][1], p[i][2]);
  }

  while (1)
  {
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
    if (list != NULL)
      if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
        break;

    lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    struct mo_0 *mo_p0[numsources];
    struct mo_p *mo_p[Nmom][numsources];
    for (int i = 0; i < numsources; i++)
    {
      mo_p0[i] = (struct mo_0 *)malloc(sizeof(struct mo_0));
      for (int j = 0; j < Nmom; j++)
      {
        mo_p[j][i] = (struct mo_p *)malloc(sizeof(struct mo_p));
      }
      lprintf("MAIN", 0, "Initiating mo, source = %d\n", i);
      init_mo_0(mo_p0[i]);
      for (int j = 0; j < Nmom; j++)
        init_mo_p(mo_p[j][i], p[j][0], p[j][1], p[j][2]);
    }

    for (int src = 0; src < numsources; ++src)
    {
      struct src_common src0;
      struct src_p *src_pn = (struct src_p *)malloc(Nmom * sizeof(struct src_p));
      struct prop_common prop0;
      struct prop_p *p_p = (struct prop_p *)malloc(Nmom * sizeof(struct prop_p));

      init_src_common(&src0, tau);
      make_prop_common(&prop0, &src0, 4, tau, mes_var.bc);
      gen_mo_0(mo_p0[src], &prop0, &src0, tau);

      for (int i = 0; i < Nmom; i++)
      {
        init_src_p(src_pn + i, &src0, p[i][0], p[i][1], p[i][2]);
        make_prop_p(p_p + i, src_pn + i, &src0, 4, tau, mes_var.bc);
        gen_mo_p(mo_p[i][src], &prop0, p_p + i, &src0, tau);
      }

      free_src_common(&src0);
      free_prop_common(&prop0);
      for (int i = 0; i < Nmom; i++)
      {
        free_src_p(src_pn + i);
        free_prop_p(p_p + i);
      }
    }
    lprintf("MAIN", 0, "num sources: %d, path: %s\n", numsources, path);
    IOold_0(mo_p0, numsources, path, cnfg_filename);
    //IO_json_0(mo_p0, numsources, path,cnfg_filename);
    for (int i = 0; i < Nmom; i++)
    {
      IOold_p(mo_p[i], numsources, path, cnfg_filename);
      //IO_json_p(mo_p[i], numsources, path,cnfg_filename);
    }

    for (int src = 0; src < numsources; src++)
    {
      free_mo_0(mo_p0[src]);
      for (int i = 0; i < Nmom; i++)
      {
        free_mo_p(mo_p[i][src]);
      }
    }

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Configuration : analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

    if (list == NULL)
      break;
  }
  lprintf("DEBUG", 0, "ALL done, deallocating\n");

  if (list != NULL)
    fclose(list);
  finalize_process();
  free_BCs();
  free_gfield(u_gauge);
  free_propagator_eo();

  return 0;
}
