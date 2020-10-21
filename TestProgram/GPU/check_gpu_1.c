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
void random_spinor_field_cpu(spinor_field *s)
{
	int i;
	geometry_descriptor *type = s->type;
	for (i = 0; i < type->local_master_pieces; i++)
		gauss((double *)(s->ptr + (type->master_start[i] - type->master_shift)), (type->master_end[i] - type->master_start[i] + 1) * sizeof(suNf_spinor) / sizeof(double));
}
void unit_array(double *a, int len)
{
      for(int i=0;i<len;i++)
      {
            a[i]=1.;
      }
}
void unit_spinor_field_cpu(spinor_field *s)
{
	int i;
	geometry_descriptor *type = s->type;
	for (i = 0; i < type->local_master_pieces; i++)
		unit_array((double *)(s->ptr + (type->master_start[i] - type->master_shift)), (type->master_end[i] - type->master_start[i] + 1) * sizeof(suNf_spinor) / sizeof(double));
}

int main(int argc, char *argv[])
{

      spinor_field *test,*original;
      int return_value = 0;
      /* setup process id and communications */
      logger_map("DEBUG", "debug");

      setup_process(&argc, &argv);

      test=alloc_spinor_field_f(2, &glattice);
      original=alloc_spinor_field_f(2,&glattice);
      random_spinor_field_cpu(&original[0]);
      //gaussian_spinor_field(&original[0]);
      spinor_field_copy_f_cpu(&test[0],&original[0]);

      /*field copy check*/
      //spinor_field_sub_assign_f(&test[0],&original[0]);
      //double d0=spinor_field_sqnorm_f(&test[0])/spinor_field_sqnorm_f(&original[0]);
      //lprintf("GPU",0,"copy difference: %.2e\n",d0);
      /*end field copy check*/

      spinor_field_togpuformat(&test[1], &test[0]);
      spinor_field_tocpuformat(&test[0], &test[1]);

      spinor_field_sub_assign_f_cpu(&test[0],&original[0]);
      double d=spinor_field_sqnorm_f_cpu(&test[0])/spinor_field_sqnorm_f_cpu(&original[0]);
      lprintf("GPU TEST",0,"norm: %.2e\n",spinor_field_sqnorm_f_cpu(&original[0]));
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
