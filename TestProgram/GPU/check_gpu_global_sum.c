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

      double *vector,sum_cpu,sum_gpu,difference;
      int n=1024;
      int return_value = 0;
      /* setup process id and communications */
      logger_map("DEBUG", "debug");

      setup_process(&argc, &argv);

      cudaMallocManaged((void **) &vector,n*sizeof(double),cudaMemAttachGlobal);
      unit_array(vector,n);
      sum_cpu=0.;
      sum_gpu=0.;
      difference=0.;
      /******SUM ON CPU******/
      for(int i=0;i<n;i++)
      {
            //printf("%d  -  %.2e\n",i,vector[i]);
            sum_cpu+=vector[i];
      }
      lprintf("GPU SUM TEST",0,"global sum cpu: %.3e\n",sum_cpu);
      /******SUM ON GPU******/
      sum_gpu=global_sum_gpu(vector,n);
      lprintf("GPU SUM TEST",0,"global sum gpu: %.3e\n",sum_gpu);
      /******DIFFERENCE TEST******/
      difference=sum_cpu-sum_gpu;
      if(difference>EPSILON)
      {
            lprintf("GPU SUM TEST",0,"difference: %.3e - test failed\n",difference);
            return_value +=1;
      }
      else
      {
            lprintf("GPU SUM TEST",0,"difference: %.3e - test passed\n",difference);
      }
      
      cudaFree(vector);
      finalize_process();
      return return_value;
}
