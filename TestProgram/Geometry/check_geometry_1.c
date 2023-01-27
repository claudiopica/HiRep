/*******************************************************************************
*
* Testing geometry
*
*******************************************************************************/

#include "libhr.h"

int main(int argc,char *argv[])
{
   /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  lprintf("test",1,"GLATTICE GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glattice.gsize_gauge,glattice.gsize_spinor, glattice.master_shift);
  lprintf("test",1,"GLATEVEN GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glat_even.gsize_gauge,glat_even.gsize_spinor, glat_even.master_shift);
  lprintf("test",1,"GLAT_ODD GSIZE=%d SPINORSIZE=%d SHIFT=%d\n",glat_odd.gsize_gauge,glat_odd.gsize_spinor, glat_odd.master_shift);

  setup_gauge_fields();

  #ifdef WITH_NEW_GEOMETRY
    int errors = test_define_geometry();  
  #else
    int errors = 0; //the test will exit with error if failure occurs
    test_geometry_mpi_eo();
  #endif

  sendbuf_report();

  finalize_process();

  return errors;
}
