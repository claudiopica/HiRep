/*******************************************************************************
*
* Testing geometry
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "geometry.h"
#include "logger.h"
#include "random.h"
#include "setup.h"


int main(int argc,char *argv[])
{
   /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  test_geometry_mpi_eo();
  
  finalize_process();

  return 0;
}
