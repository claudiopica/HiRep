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


typedef struct _format_type
{
  char name[256];
  void (*read)(char *);
  void (*write)(char *);
} format_type;

#define nformats 2

format_type format[nformats] = {
    {.name = "openQCD to HiRep", .read = read_gauge_field_openQCD, .write = write_gauge_field},
    {.name = "HiRep to openQCD", .read = read_gauge_field, .write = write_gauge_field_openQCD}};

char input_filename[1024];
char output_filename[1024];

format_type *conf_format;

void print_cmdline_info()
{
  error(1, 1, "parse_cmdline [converter_openQCD.c]", "\nSyntax (1): converter_openQCD -m\n* Show compilation information.\n\nSyntax (2): converter_openQCD -i <input file> [-o <log_file>] -H <HiRep_config>\n* Convert HiRep_config to openQCD format \n\nSyntax (3): converter_openQCD -i <input file> [-o <log_file>] -O <openQCD_config>\n* Convert openQCD_config to HiRep format \n\n ");
}

static void converter_read_cmdline(int argc, char *argv[])
{
  int i;
  int aO = 0, aH = 0;
  int aopenqcd = 0, ahirep = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-O") == 0)
    {
      aopenqcd = 1;
      aO = i;
    }
    else if (strcmp(argv[i], "-H") == 0)
    {
      ahirep = 1;
      aH = i;
    }
  }

  if (aopenqcd + ahirep != 1)
  {
    lprintf("ERROR", 0, "Incompatible options -O and -H.\n");
    print_cmdline_info();
  }

  if (aO != 0)
  {
    strcpy(input_filename, argv[aO + 1]);
    sprintf(output_filename, "HiRep_%s", input_filename);
  }
  if (aH != 0)
  {
    strcpy(input_filename, argv[aH + 1]);
    sprintf(output_filename, "openQCD_%s", input_filename);
  }

  if (aopenqcd == 1)
  {
    conf_format = format;
  }
  else
  {
    conf_format = format + 1;
  }
}

int main(int argc, char *argv[])
{

  //Read general command line arguments and input file
  setup_process(&argc, &argv);

  //Read specific command line arguments
  converter_read_cmdline(argc, argv);

  setup_gauge_fields();

  //Use reader from HiRep or OpenQCD and it writes in the other way around
  conf_format->read(input_filename);

  conf_format->write(output_filename);

  finalize_process();

  return 0;
}
