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

#define _TEST_LIN_ALG(sfsize,s1,s2,s3,i,kernel,sqn_s1,sqn_s2,sqn_s3)		\
  do {										\
      sqn_s1 = sqn_s2 = sqn_s3 = 0.0;	 					\
      for (int i=0;i<(sfsize);i++){						\
        random_spinor_field_cpu(&s1[i]); 					\
        spinor_field_copy_f_cpu(&s2[i],&s1[i]);					\
        spinor_field_copy_f_cpu(&s3[i],&s1[i]);					\
        spinor_field_copy_to_gpu_f(&s1[i]); 					\
        spinor_field_zero_f_cpu(&s1[i]);    					\
        spinor_field_copy_to_gpu_f(&s2[i]);					\
        spinor_field_zero_f_cpu(&s2[i]);					\
        spinor_field_copy_to_gpu_f(&s3[i]);					\
        spinor_field_zero_f_cpu(&s3[i]);					\
	kernel;        								\
        spinor_field_copy_from_gpu_f(&s1[i]);					\
        spinor_field_copy_from_gpu_f(&s2[i]);					\
        spinor_field_copy_from_gpu_f(&s3[i]);					\
        sqn_s1 += spinor_field_sqnorm_f_cpu(&s1[i]);  				\
        sqn_s2 += spinor_field_sqnorm_f_cpu(&s2[i]);  				\
        sqn_s3 += spinor_field_sqnorm_f_cpu(&s3[i]);  				\
      }      									\
  } while(0)

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

void print_spinor_field_cpu(spinor_field *s)
{
        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){

             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             double *r = (double*)(_FIELD_AT(s, start));

             for(int i=0;i<N*sizeof(suNf_spinor)/sizeof(double);i++) {

                 printf("spinor[%d,%d] = %f\n", ixp, i, r[i]);

             }
        }

        printf("\n\nDone.\n\n");

}


int main(int argc, char *argv[])
{

      spinor_field *sf1, *sf2, *sf3;
      double sqnorm_sf1, sqnorm_sf2, sqnorm_sf3, check;
      int sfsize = 1;
      
      double r, s;
      hr_complex c, d;
      
      /* setup process id and communications */
      //logger_setlevel(0,10000);
      logger_map("DEBUG", "debug");
      
      setup_process(&argc, &argv);

      // Allocate memory for CPU and GPU spinor fields
      sf1=alloc_spinor_field_f(sfsize, &glattice);
      sf2=alloc_spinor_field_f(sfsize, &glattice);
      sf3=alloc_spinor_field_f(sfsize, &glattice);
       
      //unit_spinor_field_cpu(&sf1[0]);       
      //printf("sqnorm sf1 = %e\n", spinor_field_sqnorm_f_cpu(&sf1[0]));
      //print_spinor_field_cpu(&sf1[0]);

        random_spinor_field_cpu(&sf1[0]); 					
        spinor_field_copy_f_cpu(&sf2[0],&sf1[0]);					
        spinor_field_copy_to_gpu_f(&sf1[0]); 					
        spinor_field_zero_f_cpu(&sf1[0]);    					
        spinor_field_copy_to_gpu_f(&sf2[0]);					
        spinor_field_zero_f_cpu(&sf2[0]);					
	spinor_field_mul_f(&sf2[0],10.0,&sf1[0]);
        spinor_field_copy_from_gpu_f(&sf1[0]);
        spinor_field_copy_from_gpu_f(&sf2[0]);					
        printf("sqnorm sf1 = %e\n", spinor_field_sqnorm_f_cpu(&sf1[0]));
        printf("sqnorm sf2 = %e\n", spinor_field_sqnorm_f_cpu(&sf2[0]));

 
      // s2=r*s1 
      r = 10.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_mul_f(&sf2[i],r,&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      printf("sqnorm sf2 = %f\n", sqnorm_sf2);
      printf("sqnorm sf1 = %f\n", sqnorm_sf1*pow(r,2));
      check = fabs(sqnorm_sf1*pow(r,2)-sqnorm_sf2);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2=r*s1):\t\t %e\n", check);	
      
      // s2=c*s1 
      c = 2.0+1.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_mulc_f(&sf2[i],c,&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(cabs(c),2)-sqnorm_sf2);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2=c*s1):\t\t %e\n", check);	

      // s2+=r*s1 
      r = 2.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_mul_add_assign_f(&sf2[i],r,&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(1+r,2)-sqnorm_sf2);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2+=r*s1):\t\t %e\n", check);	

      // s2+=c*s1 
      c = 2.0+1.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_mulc_add_assign_f(&sf2[i],c,&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(cabs(1+c),2)-sqnorm_sf2);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2+=c*s1):\t\t %e\n", check);	

      // s3=s1+s2 
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_add_f(&sf3[i],&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(4*sqnorm_sf1-sqnorm_sf3);	 		
      lprintf("GPU TEST",2,"Kernel Check (s3=s1+s2):\t\t %e\n", check);	
      
      // s3=s1-s2 
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_sub_f(&sf3[i],&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = sqnorm_sf3;	 		
      lprintf("GPU TEST",2,"Kernel Check (s3=s1-s2):\t\t %e\n", check);	

      // s2+=s1 
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_add_assign_f(&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(4*sqnorm_sf1-sqnorm_sf2);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2+=s1):\t\t %e\n", check);	

      // s2-=s1 
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_sub_assign_f(&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = sqnorm_sf2;	 		
      lprintf("GPU TEST",2,"Kernel Check (s2-=s1):\t\t %e\n", check);	

      // s2=0
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_zero_f(&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = sqnorm_sf2;	 		
      lprintf("GPU TEST",2,"Kernel Check (s2=0):\t\t %e\n", check);	

      // s2=-s1 
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_minus_f(&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf2-sqnorm_sf1);	 		
      lprintf("GPU TEST",2,"Kernel Check (s2=-s1):\t\t %e\n", check);	

      // s3=r*s1+s*s2 
      r = 2.0; s = 3.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_lc_f(&sf3[i],r,&sf1[i],s,&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(r+s,2)-sqnorm_sf3);	 		
      lprintf("GPU TEST",2,"Kernel Check (s3=r*s1+s*s2):\t %e\n", check);	

      // s3+=r*s1+s*s2 
      r = 2.0; s = 3.0;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_lc_add_assign_f(&sf3[i],r,&sf1[i],s,&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(1+r+s,2)-sqnorm_sf3);	 		
      lprintf("GPU TEST",2,"Kernel Check (s3+=r*s1+s*s2):\t %e\n", check);	

      // s3=c*s1+d*s2 
      c = 2.0+3.0*I; s = 3.0+4.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_clc_f(&sf3[i],c,&sf1[i],d,&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(cabs(c+d),2)-sqnorm_sf3);	 		
      lprintf("GPU TEST",2,"Kernel Check (s3=c*s1+d*s2):\t %e\n", check);	

      // s3+=c*s1+d*s2 
      c = 2.0+3.0*I; s = 3.0+4.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_clc_add_assign_f(&sf3[i],c,&sf1[i],d,&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = fabs(sqnorm_sf1*pow(cabs(1+c+d),2)-sqnorm_sf3);	 		
      lprintf("GPU TEST",2,"Kernel Check (s3+=c*s1+d*s2):\t %e\n", check);	

      // s2=g5*s1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_g5_f(&sf2[i],&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = 0.0;
      for (int i=0;i<(sfsize);i++){
          spinor_field_g5_f_cpu(&sf3[i],&sf1[i]);
          spinor_field_sub_f_cpu(&sf1[i],&sf3[i],&sf2[i]);
      	  check += spinor_field_sqnorm_f_cpu(&sf1[i]);      
      }
      lprintf("GPU TEST",2,"Kernel Check (s2=g5*s1):\t\t %e\n", check);	

      // s1=g5*s1
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_g5_assign_f(&sf2[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = 0.0;
      for (int i=0;i<(sfsize);i++){
          spinor_field_g5_assign_f_cpu(&sf1[i]);
          spinor_field_sub_f_cpu(&sf3[i],&sf1[i],&sf2[i]);
      	  check += spinor_field_sqnorm_f_cpu(&sf3[i]);      
      }
      lprintf("GPU TEST",2,"Kernel Check (s1=g5*s1):\t\t %e\n", check);

      // s2+=c*g5*s1
      c = 2.0+3.0*I;
      _TEST_LIN_ALG(sfsize,sf1,sf2,sf3,i,
 		    spinor_field_g5_mulc_add_assign_f(&sf2[i],c,&sf1[i]),
                    sqnorm_sf1,sqnorm_sf2,sqnorm_sf3
      );
      check = 0.0;
      for (int i=0;i<(sfsize);i++){
          spinor_field_g5_mulc_add_assign_f_cpu(&sf3[i],c,&sf1[i]);
          spinor_field_sub_f_cpu(&sf1[i],&sf3[i],&sf2[i]);
      	  check += spinor_field_sqnorm_f_cpu(&sf1[i]);      
      }
      lprintf("GPU TEST",2,"Kernel Check (s2+=c*g5*s1):\t %e\n", check);

      free_spinor_field_f(sf1);
      free_spinor_field_f(sf2);
      free_spinor_field_f(sf3);

}
