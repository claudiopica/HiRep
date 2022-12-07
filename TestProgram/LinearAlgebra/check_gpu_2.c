/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
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

#define _TEST_LIN_ALG(sfsize,s1,s2,s3,s4,i,opt_gpu,opt_cpu,check)        \
  do {                                                                   \
      check = 0.0;                                                       \
      for (int i=0;i<(sfsize);i++){                                      \
        gaussian_spinor_field(&s1[i]);                                   \
        copy_to_gpu_spinor_field_f(&s1[i]);                              \
        spinor_field_zero_f_cpu(&s1[i]);                                 \
        gaussian_spinor_field(&s2[i]);                                   \
        copy_to_gpu_spinor_field_f(&s2[i]);                              \
        spinor_field_zero_f_cpu(&s2[i]);                                 \
        spinor_field_zero_f_cpu(&s3[i]);                                 \
        copy_to_gpu_spinor_field_f(&s3[i]);                              \
        opt_gpu;                                                         \
        copy_from_gpu_spinor_field_f(&s1[i]);                            \
        copy_from_gpu_spinor_field_f(&s2[i]);                            \
        copy_from_gpu_spinor_field_f(&s3[i]);                            \
        spinor_field_zero_f_cpu(&s4[i]);                                 \
        opt_cpu;                                                         \
        spinor_field_sub_assign_f_cpu(&s4[i],&s3[i]);                    \
        check += spinor_field_sqnorm_f_cpu(&s4[i]);                      \
      }                                                                  \
  } while(0)

#define _TEST_RED_SUM(sfsize,s1,s2,i,opt_gpu,opt_cpu,check)              \
  do {                                                                   \
      check = 0.0;                                                       \
      for (int i=0;i<(sfsize);i++){                                      \
        gaussian_spinor_field(&s1[i]);                                   \
        copy_to_gpu_spinor_field_f(&s1[i]);                              \
        spinor_field_zero_f_cpu(&s1[i]);                                 \
        gaussian_spinor_field(&s2[i]);                                   \
        copy_to_gpu_spinor_field_f(&s2[i]);                              \
        spinor_field_zero_f_cpu(&s2[i]);                                 \
        check += opt_gpu;                                                \
        copy_from_gpu_spinor_field_f(&s1[i]);                            \
        copy_from_gpu_spinor_field_f(&s2[i]);                            \
        check -= opt_cpu;                                                \
      }                                                                  \
  } while(0)

#define _TEST_RED_SUM2(sfsize,s1,s2,i,opt_gpu,opt_cpu,check)             \
  do {                                                                   \
      check = 0.0;                                                       \
      for (int i=0;i<(sfsize);i++){                                      \
        gaussian_spinor_field(&s1[i]);                                   \
        copy_to_gpu_spinor_field_f(&s1[i]);                              \
        spinor_field_zero_f_cpu(&s1[i]);                                 \
        gaussian_spinor_field(&s2[i]);                                   \
        copy_to_gpu_spinor_field_f(&s2[i]);                              \
        spinor_field_zero_f_cpu(&s2[i]);                                 \
        check += opt_gpu;                                                \
        copy_from_gpu_spinor_field_f(&s1[i]);                            \
        copy_from_gpu_spinor_field_f(&s2[i]);                            \
        check -= opt_cpu;                                                \
      }                                                                  \
  } while(0)

static double EPSILON=1.e-12;

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

void print_spinor_field_cpu(spinor_field *s)
{
        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){

             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             hr_complex *r = (hr_complex*)(_FIELD_AT(s, start));

             for(int i=0;i<N*sizeof(suNf_spinor)/sizeof(hr_complex);i++) {

                 printf("spinor[%d,%d] = %f, %f\n", ixp, i, creal(r[i]), cimag(r[i]));

             }
        }

        printf("\n\nDone.\n\n");
}


int main(int argc, char *argv[])
{
      spinor_field *sf1, *sf2, *sf3, *sf4;
      int sfsize = 2;
      double r, s, check, sc_cpu, sc_gpu;
      hr_complex c, d;

      /* setup process id and communications */
      //logger_setlevel(0,10000);
      logger_map("DEBUG", "debug");
      setup_process(&argc, &argv);

      // Allocate memory for CPU and GPU spinor fields
      sf1=alloc_spinor_field_f(sfsize, &glattice);
      sf2=alloc_spinor_field_f(sfsize, &glattice);
      sf3=alloc_spinor_field_f(sfsize, &glattice);
      sf4=alloc_spinor_field_f(sfsize, &glattice);

      // sf3=r*sf1
      r = 10.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_mul_f(&sf3[i],r,&sf1[i]),
             spinor_field_mul_f_cpu(&sf4[i],r,&sf1[i]),
                    check
      );
      //printf("sqnorm sf3 = %f\n", spinor_field_sqnorm_f_cpu(&sf3[0]));
      //printf("sqnorm sf4 = %f\n", spinor_field_sqnorm_f_cpu(&sf4[0]));
      lprintf("GPU TEST",2,"Kernel Check (s2=r*s1):\t\t %.10e\n", check);

      // sf3=c*sf1
      c = 2.0+1.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_mulc_f(&sf3[i],c,&sf1[i]),
             spinor_field_mulc_f_cpu(&sf4[i],c,&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2=c*s1):\t\t %.10e\n", check);

      // sf3+=r*sf1
      r = 2.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_mul_add_assign_f(&sf3[i],r,&sf1[i]),
             spinor_field_mul_add_assign_f_cpu(&sf4[i],r,&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2+=r*s1):\t\t %.10e\n", check);

      // sf3+=c*sf1
      c = 2.0+1.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_mulc_add_assign_f(&sf3[i],c,&sf1[i]),
             spinor_field_mulc_add_assign_f_cpu(&sf4[i],c,&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2+=c*s1):\t\t %.10e\n", check);

      // sf3=sf1+sf2
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_add_f(&sf3[i],&sf2[i],&sf1[i]),
             spinor_field_add_f_cpu(&sf4[i],&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3=s1+s2):\t\t %.10e\n", check);

      // sf3=sf1-sf2
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_sub_f(&sf3[i],&sf2[i],&sf1[i]),
             spinor_field_sub_f_cpu(&sf4[i],&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3=s1-s2):\t\t %.10e\n", check);

      // sf3+=sf1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_add_assign_f(&sf3[i],&sf1[i]),
             spinor_field_add_assign_f_cpu(&sf4[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2+=s1):\t\t %.10e\n", check);

      // sf3-=sf1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_sub_assign_f(&sf3[i],&sf1[i]),
             spinor_field_sub_assign_f_cpu(&sf4[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2-=s1):\t\t %.10e\n", check);

      // sf3=0
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_zero_f(&sf3[i]),
             spinor_field_zero_f_cpu(&sf4[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2=0):\t\t %.10e\n", check);

      // sf3=-sf1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_minus_f(&sf3[i],&sf1[i]),
             spinor_field_minus_f_cpu(&sf4[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2=-s1):\t\t %.10e\n", check);

      // sf3=r*sf1+s*sf2
      r = 2.0; s = 3.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_lc_f(&sf3[i],r,&sf1[i],s,&sf2[i]),
             spinor_field_lc_f_cpu(&sf4[i],r,&sf1[i],s,&sf2[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3=r*s1+s*s2):\t %.10e\n", check);

      // sf3+=r*sf1+s*sf2
      r = 2.0; s = 3.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_lc_add_assign_f(&sf3[i],r,&sf1[i],s,&sf2[i]),
             spinor_field_lc_add_assign_f_cpu(&sf4[i],r,&sf1[i],s,&sf2[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3+=r*s1+s*s2):\t %.10e\n", check);

      // sf3=c*sf1+d*sf2
      c = 2.0+3.0*I; s = 3.0+4.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_clc_f(&sf3[i],c,&sf1[i],d,&sf2[i]),
             spinor_field_clc_f_cpu(&sf4[i],c,&sf1[i],d,&sf2[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3=c*s1+d*s2):\t %.10e\n", check);

      // sf3+=c*sf1+d*sf2
      c = 2.0+3.0*I; s = 3.0+4.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_clc_add_assign_f(&sf3[i],c,&sf1[i],d,&sf2[i]),
             spinor_field_clc_add_assign_f_cpu(&sf4[i],c,&sf1[i],d,&sf2[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s3+=c*s1+d*s2):\t %.10e\n", check);

      // sf3=g5*sf1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_g5_f(&sf3[i],&sf1[i]),
             spinor_field_g5_f_cpu(&sf4[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2=g5*s1):\t\t %.10e\n", check);

      // sf3+=c*g5*sf1
      c = 2.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,sf4,i,
             spinor_field_g5_mulc_add_assign_f(&sf3[i],c,&sf1[i]),
             spinor_field_g5_mulc_add_assign_f_cpu(&sf4[i],c,&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (s2+=c*g5*s1):\t %.10e\n", check);

      // <sf1,sf1>
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_sqnorm_f(&sf1[i]),
                    spinor_field_sqnorm_f_cpu(&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (<s1,s1>):\t\t %.10e\n", check);

      // Re<sf2,sf1>
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_prod_re_f(&sf2[i],&sf1[i]),
                    spinor_field_prod_re_f_cpu(&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check2 (Re<s1,s2>):\t\t %.10e\n", check);
      _TEST_RED_SUM2(sfsize,sf1,sf2,i,
                    spinor_field_sqnorm_f(&sf1[i]),
                    spinor_field_sqnorm_f_cpu(&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (<s1,s1>):\t\t %.10e\n", check);

      // Im<sf2,sf1>
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_prod_im_f(&sf2[i],&sf1[i]),
                    spinor_field_prod_im_f_cpu(&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (Im<s1,s2>):\t\t %.10e\n", check);

      // <sf2,sf1>
      hr_complex res;
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_prod_f(&sf2[i],&sf1[i]),
                    spinor_field_prod_f_cpu(&sf2[i],&sf1[i]),
                    res
      );
      lprintf("GPU TEST",2,"Kernel Check (<s1,s2>):\t\t %.2e, %.2e\n", creal(res), cimag(res));

      // Re<g5*sf2,sf1>
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_g5_prod_re_f(&sf2[i],&sf1[i]),
                    spinor_field_g5_prod_re_f_cpu(&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (Re<g5*s1,s2>):\t %.10e\n", check);

      // Im<g5*sf2,sf1>
      _TEST_RED_SUM(sfsize,sf1,sf2,i,
                    spinor_field_g5_prod_im_f(&sf2[i],&sf1[i]),
                    spinor_field_g5_prod_im_f_cpu(&sf2[i],&sf1[i]),
                    check
      );
      lprintf("GPU TEST",2,"Kernel Check (Im<g5*s1,s2>):\t %.10e\n", check);

      free_spinor_field_f(sf1);
      free_spinor_field_f(sf2);
      free_spinor_field_f(sf3);
      free_spinor_field_f(sf4);

       finalize_process();
}
