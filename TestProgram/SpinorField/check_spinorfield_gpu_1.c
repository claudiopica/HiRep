
/******************************************************************************
 *
 * File check11.c
 *
 * Consistency checks on the programs in the module linalg
 *
 * Author: luigi del debbio <luigi.del.debbio@ed.ac.uk>
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "linear_algebra.h"
#include "representation.h"
#include "global.h"
#include "logger.h"
#include "utils.h"
#include "io.h"

#include "communications.h"

#define MAX_ROTATE 50

static complex v[25];
static double EPSILON=1.e-12;
static spinor_field *ppk[5];

double sfdiff (spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, &glattice);
  alloc_spinor_field_f_gpu(1, tmp);
  spinor_field_copy_f_cpu(tmp,sf);
  spinor_field_copy_to_gpu_f(tmp);
  spinor_field_sub_f(tmp,tmp,sf);
  res= spinor_field_sqnorm_f(tmp);
  free_spinor_field_gpu(tmp);
  free_spinor_field(tmp);
  return res;
}

int main(int argc,char *argv[])
{
  
  spinor_field *sf1,*sf2;
  int sfsize = 10;
  double norm_cpu;
  double norm_gpu;
  double res_cpu,res_gpu;
  complex c1,c2;
  int i;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");

#ifdef WITH_MPI
  sprintf(pame,">out_%d",PID); logger_stdout(pame);
  sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
#endif

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input(glb_var.read,"test_input");
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed);


  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */


  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("CPTEST",0,"gsize=%d\n",glattice.gsize);
  lprintf("CPTEST",0,"nbuffers=%d\n",glattice.nbuffers);
  lprintf("CPTEST",0,"lmp=%d\n",glattice.local_master_pieces);
  lprintf("CPTEST",0,"ncopies=%d\n",glattice.ncopies);
  
//	Allocates memory for cpu spinor field. GPU pointer = NULL
  sf1=alloc_spinor_field_f(sfsize, &glattice);
  sf2=alloc_spinor_field_f(sfsize, &glattice);
  alloc_spinor_field_f_gpu(sfsize, sf1);
  alloc_spinor_field_f_gpu(sfsize,sf2);
	
// CPU part set to gaussian
  for (i=0;i<sfsize;i++){
	  
    gaussian_spinor_field(&sf1[i]);
  }
	

// Copy content of CPU field to GPU field
  for (i=0;i<sfsize;i++){
    spinor_field_copy_to_gpu_f(&sf1[i]);
  }

// Copying to the other spinor field... Now all fields are equal,    Copy from 2nd arg to 1st
  for (i=0;i<sfsize;i++){
    spinor_field_copy_f(&sf2[i],&sf1[i]);
    spinor_field_copy_f_cpu(&sf2[i],&sf1[i]);
  }  

  //Check spinor_field_prod_im
  res_gpu = spinor_field_prod_im_f(&sf1[0],&sf1[1]);
  res_cpu = spinor_field_prod_im_f_cpu(&sf1[0],&sf1[1]);
  lprintf("LA TEST",0,"Check spinor_field_prod_im\n gpu=%1.10g, cpu=%1.10g, \n gpu-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu-res_cpu);

  //Check spinor_field_g5_prod_re
  res_gpu = spinor_field_g5_prod_re_f(&sf1[0],&sf1[1]);
  res_cpu = spinor_field_g5_prod_re_f_cpu(&sf1[0],&sf1[1]);
  lprintf("LA TEST",0,"Check spinor_field_g5_prod_re\n gpu=%1.10g, cpu=%1.10g, \n gpu-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu-res_cpu);

  //Check spinor_field_sqnorm
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  lprintf("LA TEST",0,"Check spinor_field_sqnorm\n gpu=%1.10g, cpu=%1.10g, \n gpu-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu-res_cpu);

  //Check spinor_field_mul_add_assign
  spinor_field_mul_add_assign_f(&sf1[0],2.5,&sf1[1]);
  spinor_field_mul_add_assign_f_cpu(&sf1[0],2.5,&sf1[1]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_mul_add_assign\n sqnorm(gpu)=%1.10g, sqnorm(cpu)=%1.10g, \n sqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));

  //Check spinor_field_mulc_add_assign
  c1.re = 2.34;
  c1.im = -1.11;
  
  spinor_field_mulc_add_assign_f(&sf1[0],c1,&sf1[1]);
  spinor_field_mulc_add_assign_f_cpu(&sf1[0],c1,&sf1[1]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_mulc_add_assign\n sqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \n sqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));


  //Check spinor_field_mulc
  c1.re = -2.3334;
  c1.im = 1.241;
  
  spinor_field_mulc_f(&sf1[0],c1,&sf1[1]);
  spinor_field_mulc_f_cpu(&sf1[0],c1,&sf1[1]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_mulc\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));

  //Check spinor_field_sub
  spinor_field_sub_f(&sf1[0],&sf1[1],&sf1[2]);
  spinor_field_sub_f_cpu(&sf1[0],&sf1[1],&sf1[2]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_sub\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));

  //Check spinor_field_minus
  spinor_field_minus_f(&sf1[0],&sf1[1]);
  spinor_field_minus_f_cpu(&sf1[0],&sf1[1]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_minus\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));

  //Check spinor_field_lc_add_assing_f
  spinor_field_lc_add_assign_f(&sf1[0],1.1,&sf1[1],-1.5,&sf2[2]);
  spinor_field_lc_add_assign_f_cpu(&sf1[0],1.1,&sf1[1],-1.5,&sf2[2]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_lc_add_assing_f\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));


  //Check spinor_field_clc_add_assing_f
  c1.re = 5.34;c1.im = -1.11;
  c2.re = 0.34;c1.im = -0.11;
  
  spinor_field_clc_add_assign_f(&sf1[0],c1,&sf1[1],c2,&sf2[2]);
  spinor_field_clc_add_assign_f_cpu(&sf1[0],c1,&sf1[1],c2,&sf2[2]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_clc_add_assing_f\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));


  //Check spinor_field_g5_assign
  spinor_field_g5_assign_f(&sf1[0]);
  spinor_field_g5_assign_f_cpu(&sf1[0]);
  
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_g5_assign\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));

  //Check spinor_field_lc2

  
  spinor_field_lc2_f(1.1,-5.5,&sf1[0],&sf1[1]);
  spinor_field_lc2_f_cpu(1.1,-5.5, &sf1[0],&sf1[1]);

  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  
  lprintf("LA TEST",0,"Check spinor_field_clc_add_assing_f\nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g, \nsqnorm(gpu-cpu)= %1.10g\n\n",res_gpu,res_cpu,sfdiff(&sf1[0]));
  
  
  lprintf("LA TEST",0,"DONE!\n");

  free_spinor_field_gpu(sf1);
  free_spinor_field_gpu(sf2);
  free_spinor_field(sf1);
  free_spinor_field(sf2);
	
  finalize_process();

}
