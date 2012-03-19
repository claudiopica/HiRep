
/******************************************************************************
 *
 * File global_sum_gpu.c
 *
 * Speed and validity of improved global sum
 *
 * Author: Ulrik Ishøj Søndergaard
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

/*extern "C" {
double spinor_field_prod_re_old_f(spinor_field* s1, spinor_field* s2);
double spinor_field_prod_im_old_f(spinor_field* s1, spinor_field* s2);
double spinor_field_g5_prod_im_old_f(spinor_field* s1, spinor_field* s2);
double spinor_field_g5_prod_re_old_f(spinor_field* s1, spinor_field* s2);
complex spinor_field_prod_old_f(spinor_field* s1, spinor_field* s2);
double spinor_field_sqnorm_old_f(spinor_field* s1);
}*/

double sfdiff_gpu (spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, sf->type);
  spinor_field_copy_f_cpu(tmp,sf);
  spinor_field_copy_to_gpu_f(tmp);
  spinor_field_sub_f(tmp,tmp,sf);
  res= spinor_field_sqnorm_f(tmp);
  free_spinor_field_f(tmp);
  return res;
}

double sfdiff (spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, sf->type);

  spinor_field_copy_f(tmp,sf);
  spinor_field_copy_from_gpu_f(tmp);
  spinor_field_sub_f_cpu(tmp,tmp,sf);

  res=spinor_field_sqnorm_f_cpu(tmp);

  free_spinor_field_f(tmp);
  return res;
}

int main(int argc,char *argv[])
{
  
  spinor_field *sf1,*sf2;
	int sfsize = 2;
//  double norm_cpu;
	// double norm_gpu;
	double res_cpu,res_gpu, res_gpu_opt;
	complex c_res_cpu,c_res_gpu, c_res_gpu_opt;
  complex c1,c2;
  int i;
	
	gpu_timer t1;
	float elapsed;
	

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
  
 t1 = gpuTimerStart();
  //	Allocates memory for cpu & gpu spinor field. 
 sf1=alloc_spinor_field_f(sfsize, &glattice);
//  sf2=alloc_spinor_field_f(sfsize, &glattice);
 elapsed = gpuTimerStop(t1);
 lprintf("TEST",0,"Time for allocation: %1.10gms\n",elapsed);
 
 t1 = gpuTimerStart();
// CPU part set to gaussian
  for (i=0;i<sfsize;i++){
    gaussian_spinor_field(&sf1[i]);
  }


// Copy content of CPU field to GPU field
  for (i=0;i<sfsize;i++){
    spinor_field_copy_to_gpu_f(&sf1[i]);
  }

// Copying to the other spinor field... Now all fields are equal,    Copy from 2nd arg to 1st
  //for (i=0;i<sfsize;i++){
  //  spinor_field_copy_f(&sf2[i],&sf1[i]);
//    spinor_field_copy_f_cpu(&sf2[i],&sf1[i]);
  //}  
  elapsed = gpuTimerStop(t1);
  lprintf("TEST",0,"Time for Gaussian+cpy: %1.10gms\n",elapsed);
	
	
  
  //Check Gaussian Spinor Field
  /*gaussian_spinor_field_cpu(&sf1[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
  res_gpu = spinor_field_sqnorm_f(&sf1[0]);
  lprintf("TEST",0,"Check gaussian_spinor_field\n gpu=%1.10g, cpu=%1.10g, \n gpu-cpu= %1.10g\n\n", res_gpu, res_cpu, res_gpu-res_cpu);*/
	
  //Check spinor_field_prod_re
  //	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
  //  res_gpu = spinor_field_prod_re_old_f(&sf1[0],&sf1[1]);
  res_cpu = spinor_field_prod_re_f_cpu(&sf1[0],&sf1[1]);
  t1 = gpuTimerStart();
  res_gpu_opt = spinor_field_prod_re_f(&sf1[0],&sf1[1]);
  elapsed = gpuTimerStop(t1);
  lprintf("TEST",0,"Check spinor_field_prod_re\n old gpu=%1.10g, cpu=%1.10g, opt gpu=%1.10g, \n opt-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu_opt,res_gpu_opt-res_cpu);
  lprintf("TEST",0,"Time: =%1.10g \n\n",elapsed);
  	
  //Check spinor_field_prod_im
  //	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
  //  res_gpu = spinor_field_prod_im_old_f(&sf1[0],&sf1[1]);
  res_cpu = spinor_field_prod_im_f_cpu(&sf1[0],&sf1[1]);
  t1 = gpuTimerStart();
  res_gpu_opt = spinor_field_prod_im_f(&sf1[0],&sf1[1]);
  elapsed = gpuTimerStop(t1);
  lprintf("TEST",0,"Check spinor_field_prod_im\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu_opt,res_gpu_opt-res_cpu);
  lprintf("TEST",0,"Time: =%1.10g \n\n",elapsed);
  
	//Check spinor_field_prod
//	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
  //    c_res_gpu = spinor_field_prod_old_f(&sf1[0],&sf1[1]);
	c_res_cpu = spinor_field_prod_f_cpu(&sf1[0],&sf1[1]);
	c_res_gpu_opt = spinor_field_prod_f(&sf1[0],&sf1[1]);
	
	lprintf("TEST",0,"Im: Check spinor_field_prod\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n",c_res_gpu.im, c_res_cpu.im,c_res_gpu_opt.im,c_res_gpu_opt.im-c_res_cpu.im);
	lprintf("TEST",0,"Re: Check spinor_field_prod\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n\n",c_res_gpu.re,c_res_cpu.re,c_res_gpu_opt.re,c_res_gpu_opt.re-c_res_cpu.re);
	
	//Check spinor_field_g5_prod_re
//	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
//	res_gpu = spinor_field_g5_prod_re_old_f(&sf1[0],&sf1[1]);
	res_cpu = spinor_field_g5_prod_re_f_cpu(&sf1[0],&sf1[1]);
	res_gpu_opt = spinor_field_g5_prod_re_f(&sf1[0],&sf1[1]);
	lprintf("TEST",0,"Check spinor_field_g5_prod_re\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu_opt,res_gpu_opt-res_cpu);
	
	//Check spinor_g5_field_prod_im
//	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
//	res_gpu = spinor_field_g5_prod_im_old_f(&sf1[0],&sf1[1]);
	res_cpu = spinor_field_g5_prod_im_f_cpu(&sf1[0],&sf1[1]);
	res_gpu_opt = spinor_field_g5_prod_im_f(&sf1[0],&sf1[1]);
	lprintf("TEST",0,"Check spinor_field_g5_prod_im\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu_opt,res_gpu_opt-res_cpu);

	//Check spinor_field_sqnorm
//	for (i=0;i<sfsize;i++){ spinor_field_copy_to_gpu_f(&sf1[i]); }
//	res_gpu = spinor_field_sqnorm_old_f(&sf1[0]);
	res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
	res_gpu_opt = spinor_field_sqnorm_f(&sf1[0]);
	lprintf("TEST",0,"Check spinor_field_sqnorm\n gpu=%1.10g, cpu=%1.10g, opt=%1.10g, \n opt-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu_opt,res_gpu_opt-res_cpu);
	
	int NumberOfRuns = 500;
	

	lprintf("TEST",0,"::::::::::::spinor_field_prod_re:::::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu = spinor_field_prod_re_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu_opt = spinor_field_prod_re_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n",elapsed);
	
	lprintf("TEST",0,"::::::::::::spinor_field_prod_im:::::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu = spinor_field_prod_im_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu_opt = spinor_field_prod_im_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n\n",elapsed);
	
	
	lprintf("TEST",0,"::::::::::::spinor_field_prod::::::::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		c_res_gpu = spinor_field_prod_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		c_res_gpu_opt = spinor_field_prod_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n\n",elapsed);
	
	
	lprintf("TEST",0,"::::::::::::spinor_field_g5_prod_re::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu = spinor_field_g5_prod_re_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu_opt = spinor_field_g5_prod_re_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n\n",elapsed);
	
	
	lprintf("TEST",0,"::::::::::::spinor_field_g5_prod_im::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu = spinor_field_g5_prod_im_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu_opt = spinor_field_g5_prod_im_f(&sf1[0],&sf1[1]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n\n",elapsed);
	
	
	lprintf("TEST",0,"::::::::::::spinor_field_qsnorm::::::::::::::::\n\n");
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu = spinor_field_sqnorm_f(&sf1[0]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for original: %1.10gms\n\n",elapsed);
	t1 = gpuTimerStart();
	for (int i=0; i<NumberOfRuns; i++) {	
		res_gpu_opt = spinor_field_sqnorm_f(&sf1[0]);
	}
	elapsed = gpuTimerStop(t1);
	lprintf("TEST",0,"Time for optimized: %1.10gms\n\n",elapsed);
	
	
	
	
	lprintf("TEST",0,"DONE!\n");


  free_spinor_field_f(sf1);
  //free_spinor_field_f(sf2);
	
  finalize_process();

}
