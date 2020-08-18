/*******************************************************************************
*  
* Converter from and to openQCD format.
* Code modified by Fernando Romero-Lopez
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
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
#include "moreio.h"
#include "setup.h"

#ifdef WITH_MPI
#error Please compile without MPI!
#endif
#if !(defined(BC_X_PERIODIC)) || !(defined(BC_Y_PERIODIC)) || !(defined(BC_Z_PERIODIC)) || !(defined(BC_T_PERIODIC))
#error Please compile with Periodic Boundary Conditions
#endif

static void slower(char *str)
{
  while (*str)
  {
    *str = (char)(tolower(*str));
    ++str;
  }
}
typedef struct _format_type
{
  char name[256];
  void (*read)(char *);
  void (*write)(char *);
} format_type;

#define nformats 4

format_type format[nformats] = {
    {.name = "OBC openQCD to HiRep", .read = read_gauge_field_openQCD, .write = write_gauge_field_hirep_pbc_to_obc},
    {.name = "SF  openQCD to HiRep", .read = read_gauge_field_openQCD_SF, .write = write_gauge_field_hirep_pbc_to_sf},
    {.name = "PBC openQCD to HiRep", .read = read_gauge_field_openQCD, .write = write_gauge_field},
    {.name = "HiRep PBC to openQCD", .read = read_gauge_field, .write = write_gauge_field_openQCD}};

format_type *conf_format;

typedef struct _conf_details
{
  char cnfgin[256];
  char cnfgout[256];
  char type[128];
  input_record_t read[3];
} conf_details;

#define init_conf_details(varname)                                      \
  {                                                                     \
    .read = {                                                           \
      {"conf name", "conf name= %s", STRING_T, &((varname).cnfgin[0])}, \
      {"conf type", "conf type= %s", STRING_T, &((varname).type[0])},   \
      {NULL, NULL, 0, NULL}                                             \
    }                                                                   \
  }

int main(int argc, char *argv[])
{

  //Read general command line arguments and input file
  setup_process(&argc, &argv);

  setup_gauge_fields();

  conf_details conf_info = init_conf_details(conf_info);

  read_input(conf_info.read, get_input_filename());

  slower(conf_info.type);

  if (strcmp(conf_info.type, "openqcd_obc") == 0)
  {
    conf_format = format;
    sprintf(conf_info.cnfgout, "HiRep_obc_%s", conf_info.cnfgin);
  }
  else if (strcmp(conf_info.type, "openqcd_sf") == 0)
  {
    conf_format = format + 1;
    sprintf(conf_info.cnfgout, "HiRep_sf_%s", conf_info.cnfgin);
  }
  else if (strcmp(conf_info.type, "openqcd_pbc") == 0)
  {
    conf_format = format + 2;
    sprintf(conf_info.cnfgout, "HiRep_pbc_%s", conf_info.cnfgin);
  }
  else if (strcmp(conf_info.type, "hirep") == 0)
  {
    conf_format = format + 3;
    sprintf(conf_info.cnfgout, "openQCD_%s", conf_info.cnfgin);
  }
  else
  {
    error(0 == 0, 0, "MAIN", "Unknonw configuration type (openqcd_pbc openqcd_obc openqcd_sf hirep) are the only allowed.");
  }

  lprintf("MAIN", 0, "Converting the configuration from %s\n", conf_format->name);

  //Use reader from HiRep or OpenQCD and it writes in the other way around
  conf_format->read(conf_info.cnfgin);

  conf_format->write(conf_info.cnfgout);

  finalize_process();

  return 0;
}
