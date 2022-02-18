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

        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){
             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             gauss((double*)(_FIELD_AT(s, start)), N * sizeof(suNf_spinor) / sizeof(double));
        }
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
        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){
             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             unit_array((double*)(_FIELD_AT(s, start)), N * sizeof(suNf_spinor) / sizeof(double));
        }
}

int main(int argc, char *argv[])
{

      spinor_field *sf1, *sf2, *sf3;
      int sfsize = 1;
      double norm_cpu, norm_cpu2, norm_cpu3;

      /* setup process id and communications */
      //logger_setlevel(0,10000);
      logger_map("DEBUG", "debug");
      
      setup_process(&argc, &argv);

      // Allocate memory for CPU and GPU spinor fields
      sf1=alloc_spinor_field_f(sfsize, &glattice);
      sf2=alloc_spinor_field_f(sfsize, &glattice);
      sf3=alloc_spinor_field_f(sfsize, &glattice);
    
      printf("gsize_spinor = %d\n", sf3->type->gsize_gauge );     
      printf("sf1->ptr = %p", sf1->ptr );     
      printf("  sf1->gpu_ptr = %p\n", sf1->gpu_ptr );     
      printf("sf2->ptr = %p", sf2->ptr );     
      printf("  sf2->gpu_ptr = %p\n", sf2->gpu_ptr );     
      printf("sf3->ptr = %p", sf3->ptr );     
      printf("  sf3->gpu_ptr = %p\n", sf3->gpu_ptr );     
      
      // Assign random field to CPU
      for (int i=0;i<sfsize;i++){
        random_spinor_field_cpu(&sf1[i]);
        unit_spinor_field_cpu(&sf2[i]);
        //unit_spinor_field_cpu(&sf3[i]);
        //spinor_field_zero_f_cpu(&sf2[i]);         
        spinor_field_zero_f_cpu(&sf3[i]);         
      }

      norm_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
      printf("norm 1 = %f\n", norm_cpu);      


      for (int i=0;i<sfsize;i++){
        spinor_field_copy_f_cpu(&sf2[i],&sf1[i]);
        spinor_field_copy_to_gpu_f(&sf1[i]); //
        spinor_field_zero_f_cpu(&sf1[i]);    //
        spinor_field_copy_to_gpu_f(&sf2[i]);
        spinor_field_zero_f_cpu(&sf2[i]);
        spinor_field_copy_to_gpu_f(&sf3[i]);
        spinor_field_zero_f_cpu(&sf3[i]);
        /* s1+=r*s2 */
        double r = 2.0;
        //spinor_field_mul_f(&sf3[i],r,&sf2[i]);
        spinor_field_copy_from_gpu_f(&sf1[i]);
        spinor_field_copy_from_gpu_f(&sf2[i]);
        spinor_field_copy_from_gpu_f(&sf3[i]);
      }

      printf("sf1->ptr = %p", sf1->ptr );     
      printf("  sf1->gpu_ptr = %p\n", sf1->gpu_ptr );     
      printf("sf2->ptr = %p", sf2->ptr );     
      printf("  sf2->gpu_ptr = %p\n", sf2->gpu_ptr );     
      printf("sf3->ptr = %p", sf3->ptr );     
      printf("  sf3->gpu_ptr = %p\n", sf3->gpu_ptr );     
      
      // Calculate norm on CPU
      norm_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
      norm_cpu2 = spinor_field_sqnorm_f_cpu(&sf2[0]); 
      norm_cpu3 = spinor_field_sqnorm_f_cpu(&sf3[0]); 
      
      lprintf("GPU TEST",2,"Norm CPU 1: %.6e\n",norm_cpu); 
      lprintf("GPU TEST",2,"Norm CPU 2: %.6e\n",norm_cpu2); 
      lprintf("GPU TEST",2,"Norm CPU 3: %.6e\n",norm_cpu3); 

}
