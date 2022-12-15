/*******************************************************************************
*
* NOCOMPILE = !WITH_GPU
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
#include "new_geom_gpu.h"
#include "representation.h"
#include "memory.h"
#include "communications.h"


int main(int argc,char *argv[])
{
   /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  lprintf("test",1,"GLATTICE GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glattice.gsize_gauge,glattice.gsize_spinor, glattice.master_shift);
  lprintf("test",1,"GLATEVEN GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glat_even.gsize_gauge,glat_even.gsize_spinor, glat_even.master_shift);
  lprintf("test",1,"GLAT_ODD GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glat_odd.gsize_gauge,glat_odd.gsize_spinor, glat_odd.master_shift);

  setup_gauge_fields();
  random_u(u_gauge);
  represent_gauge_field();
  copy_to_gpu_gfield_f(u_gauge_f);

  sync_gpu_gfield_f(u_gauge_f);
  #ifdef WITH_NEW_GEOMETRY
    test_define_geometry();  
  #else
    test_geometry_mpi_eo();
  #endif

  sendbuf_report();

  finalize_process();

  return 0;
}
