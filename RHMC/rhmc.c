/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main RHMC program
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
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"

#include "communications.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

spinor_field *stmp=NULL;

void D(spinor_field *out, spinor_field *in) {
  const double m=-1.2;
  Dphi(m,out,in);
}

void H(spinor_field *out, spinor_field *in) {
  const double m=-1.2;
  g5Dphi(m,out,in);
}

void M(spinor_field *out, spinor_field *in) {
  const double m=-1.2;
  
  g5Dphi(m,stmp,in);
  g5Dphi(m,out,stmp);

}

int main(int argc,char *argv[])
{
  int i, acc;
  char tmp[256];

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,40);
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input("input_file");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",input_p.rlxd_level,input_p.rlxd_seed);
  rlxd_init(input_p.rlxd_level,input_p.rlxd_seed);

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
  test_geometry_mpi_eo();

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("CPTEST",0,"PSIGN=%d\n",PSIGN);

  lprintf("CPTEST",0,"Global Lattice\n");
  lprintf("CPTEST",0,"gsize=%d ",glattice.gsize);
  lprintf("CPTEST",0,"nbuffers=%d ",glattice.nbuffers);
  lprintf("CPTEST",0,"lmp=%d ",glattice.local_master_pieces);
  lprintf("CPTEST",0,"ncopies=%d\n",glattice.ncopies);

  lprintf("CPTEST",0,"Even Lattice\n");
  lprintf("CPTEST",0,"gsize=%d ",glat_even.gsize);
  lprintf("CPTEST",0,"nbuffers=%d ",glat_even.nbuffers);
  lprintf("CPTEST",0,"lmp=%d ",glat_even.local_master_pieces);
  lprintf("CPTEST",0,"ncopies=%d\n",glat_even.ncopies);

  lprintf("CPTEST",0,"Odd Lattice\n");
  lprintf("CPTEST",0,"gsize=%d ",glat_odd.gsize);
  lprintf("CPTEST",0,"nbuffers=%d ",glat_odd.nbuffers);
  lprintf("CPTEST",0,"lmp=%d ",glat_odd.local_master_pieces);
  lprintf("CPTEST",0,"ncopies=%d\n",glat_odd.ncopies);

  /* disable logger for MPI processes != 0 */
  if (PID!=0) { logger_disable(); }

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  random_u(u_gauge);
  represent_gauge_field();


  /* init RHMC */
  input_p.rhmc_p.integrator=&O2MN_multistep;
  input_p.rhmc_p.mshift_solver=&cg_mshift; /* this is not used in the code now */
  input_p.rhmc_p.MD_par=&input_p.int_p;
  init_rhmc(&input_p.rhmc_p);
  
  lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());
  
  acc=0;
  for(i=1;i<(10+1);++i) {
    int rr;
    double perc;
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    rr=update_rhmc();
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else {
      acc+=rr;
    }
    perc=(acc==0)?0.:(float)(100*acc)/(float)(i);
    
    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM = %ld\n",i,acc,i,perc,getMVM());
    lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
  }
  
  
   free_rhmc();


  /* free memory */
  free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
  free_gfield_f(u_gauge_f);
#endif
  /* close communications */
  finalize_process();

  return 0;

}
