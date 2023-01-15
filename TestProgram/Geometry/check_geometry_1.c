/*******************************************************************************
*
* NOCOMPILE= WITH_GPU
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

  //sync_field(u_gauge->type, 4*sizeof(*u_gauge->ptr), 0, u_gauge->ptr);
  #ifdef WITH_NEW_GEOMETRY
    test_define_geometry();  
  #else
    test_geometry_mpi_eo();
  #endif

  sendbuf_report();

  finalize_process();

  return 0;
}
