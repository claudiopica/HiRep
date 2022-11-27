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

  // sync_field(u_gauge->type, 4*sizeof(*u_gauge->ptr), 0, u_gauge->ptr);
#ifdef WITH_NEW_GEOMETRY
  test_define_geometry();  
#else
  test_geometry_mpi_eo();
#endif

  finalize_process();

  return 0;
}
