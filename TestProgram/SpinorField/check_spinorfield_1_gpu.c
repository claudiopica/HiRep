
/******************************************************************************
 *
 *
 *  NOCOMPILE= !WITH_GPU
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
#include "setup.h"
#include "communications.h"

#define MAX_ROTATE 50

static complex v[25];
static double EPSILON=1.e-12;
static spinor_field *ppk[5];

double sfdiff_gpu(spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, sf->type);
  spinor_field_copy_f_cpu(tmp,sf);
  copy_to_gpu_spinor_field_f(tmp);
  spinor_field_sub_f(tmp,tmp,sf);
  res=spinor_field_sqnorm_f(tmp);
  free_spinor_field_f(tmp);
  return res;
}

double sfdiff(spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, sf->type);

  spinor_field_copy_f(tmp,sf);
  copy_from_gpu_spinor_field_f(tmp);
  spinor_field_sub_f_cpu(tmp,tmp,sf);

  res=spinor_field_sqnorm_f_cpu(tmp);

  free_spinor_field_f(tmp);
  return res;
}

int main(int argc,char *argv[])
{

  spinor_field *sf1,*sf2, *sf3;
  int sfsize = 3;
  double norm_cpu;
  double norm_gpu;
  double res_cpu,res_cpu2,res_gpu;
  complex c1,c2;
  int i;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");

#ifdef WITH_MPI
  // Throws errors during compilation. 
  // For logging, lprintf should be used. (SAM)
  //sprintf(pame,">out_%d",PID); logger_stdout(pame);
  //sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
#endif

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE);

  /* read input file */
  read_input(glb_var.read,"test_input");
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed);


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

  lprintf("CPTEST",0,"gsize=%d\n",glattice.gsize_spinor);
  lprintf("CPTEST",0,"nbuffers=%d\n",glattice.nbuffers_spinor);
  lprintf("CPTEST",0,"lmp=%d\n",glattice.local_master_pieces);
  lprintf("CPTEST",0,"ncopies=%d\n",glattice.ncopies_spinor);


  input_gpu gpu_var = init_input_gpu(gpu_var);
  read_input(gpu_var.read,"test_input");
  init_gpu(gpu_var);

  // Allocates memory for cpu & gpu spinor field.
  sf1=alloc_spinor_field_f(sfsize, &glattice);
  sf2=alloc_spinor_field_f(sfsize, &glattice);
  sf3=alloc_spinor_field_f(sfsize, &glattice);

  lprintf("LA TEST",0,"Maximum grid size: %d\n", grid_size_max_gpu);

// CPU part set to gaussian
  for (i=0;i<sfsize;i++){
    gaussian_spinor_field(&sf1[i]);
  }

// Copy content of CPU field to GPU field
  for (i=0;i<sfsize;i++){
    copy_to_gpu_spinor_field_f(&sf1[i]);
  }

// Copying to the other spinor field... Now all fields are equal, Copy from 2nd arg to 1st
  for (i=0;i<sfsize;i++){
    spinor_field_copy_f(&sf2[i],&sf1[i]);
    spinor_field_copy_f_cpu(&sf2[i],&sf1[i]);
    spinor_field_copy_f_cpu(&sf3[i],&sf1[i]);
  }


  // Check copy to and from GPU, copy sf3 to GPU and back to CPU.
  // Check that it is still equal to sf1 on CPU.
  copy_to_gpu_spinor_field_f(&sf3[0]);
  copy_from_gpu_spinor_field_f(&sf3[0]);
  res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  res_cpu2 = spinor_field_sqnorm_f_cpu(&sf3[0]);

  lprintf("LA TEST",0,"Check copy_tofrom_spinor_field\n copied=%1.10g, ref=%1.10g, \n copied-ref= %1.10g\n\n",res_cpu2,res_cpu,res_cpu2-res_cpu);

  // gaussian_spinor_field_cpu(&sf1[0]);
  // res_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
  // spinor_field_copy_to_gpu_f(&sf1[0]);
  // res_gpu = spinor_field_sqnorm_f(&sf1[0]);

  // lprintf("LA TEST",0,"Check gaussian_spinor_field\n gpu=%1.10g, cpu=%1.10g, \n gpu-cpu= %1.10g\n\n",res_gpu,res_cpu,res_gpu-res_cpu);


  lprintf("LA TEST",0,"DONE!\n");


  free_spinor_field_f(sf1);
  free_spinor_field_f(sf2);

  finalize_process();

}
