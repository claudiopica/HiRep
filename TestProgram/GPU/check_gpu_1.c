/*******************************************************************************
*
*  NOCOMPILE= WITH_MPI
*
* Test of modules
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
#include "global.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "setup.h"
#include "logger.h"

static double EPSILON=1.e-12;

int main(int argc, char *argv[])
{

      spinor_field *test,*original;
      int return_value = 0;
      /* setup process id and communications */
      logger_map("DEBUG", "debug");

      setup_process(&argc, &argv);

      test=alloc_spinor_field_f(2, &glattice);
      original=alloc_spinor_field_f(2,&glattice);
      gaussian_spinor_field(&original[0]);
      spinor_field_copy_f(&test[0],&original[0]);

      /*field copy check*/
      //spinor_field_sub_assign_f(&test[0],&original[0]);
      //double d0=spinor_field_sqnorm_f(&test[0])/spinor_field_sqnorm_f(&original[0]);
      //lprintf("GPU",0,"copy difference: %.2e\n",d0);
      /*end field copy check*/

      spinor_field_togpuformat(&test[1], &test[0]);
      spinor_field_tocpuformat(&test[0], &test[1]);

      spinor_field_sub_assign_f(&test[0],&original[0]);
      double d=spinor_field_sqnorm_f(&test[0])/spinor_field_sqnorm_f(&original[0]);
      lprintf("GPU TEST",0,"difference: %.2e\n",d);
      if (d > EPSILON){
            lprintf("GPU",0,"Test failed ?\n");
            return_value +=1;
      }
      
      //setup_gauge_fields();

      free_spinor_field_f(test);
      free_spinor_field_f(original);

      finalize_process();
      return return_value;
}
