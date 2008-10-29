/******************************************************************************
*
* Test of hermiticity
*
******************************************************************************/

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

static double hmass=0.1;


void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void H(spinor_field *out, spinor_field *in){
   g5Dphi(-hmass,out,in);
}

void M(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
   g5Dphi_eopre_sq(-hmass, out, in);
#else
   g5Dphi_sq(-hmass, out, in);
#endif
}

void test_herm(spinor_operator S, char *name){
   spinor_field *s1, *s2, *s3, *s4;
   double tau;

#ifdef UPDATE_EO
   s1=alloc_spinor_field_f(4,&glat_even);
#else
   s1=alloc_spinor_field_f(4,&glattice);
#endif
   s2=s1+1;
   s3=s2+1;
   s4=s3+1;
   
   lprintf("RESULT",0,"Test if %s is hermitean: ",name);

   gaussian_spinor_field(s1);
   gaussian_spinor_field(s2);
   S(s3,s1);
   S(s4,s2);

   tau=spinor_field_prod_re_f(s2,s3);
   tau-=spinor_field_prod_re_f(s4,s1);
   tau+=spinor_field_prod_im_f(s2,s3);
   tau-=spinor_field_prod_im_f(s4,s1);
   tau/=sqrt(spinor_field_sqnorm_f(s1));
   tau/=sqrt(spinor_field_sqnorm_f(s2));
   if (fabs(tau)>1.e-7) 
     lprintf("RESULT",0,"FAILED ");
   else 
     lprintf("RESULT",0,"OK ");
   lprintf("RESULT",0,"[norm = %e]\n",tau);
   
   free_spinor_field(s1);

}


int main(int argc,char *argv[])
{
  char tmp[256];

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");
#ifdef WITH_MPI
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
#endif

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input("test_input");
  rlxd_init(1,12345);

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

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  lprintf("MAIN",0,"Generating a random gauge field... ");
  random_u(u_gauge);
  /* questo va rimosso quando la geometria viene sistemata! */
  lprintf("MAIN",0,"done.\n");
  represent_gauge_field();
   
	 /*
   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
   
   tau = 1./sqrt(spinor_field_sqnorm_f(s1));
   spinor_field_mul_f(s1,tau,s1);
   tau = 1./sqrt(spinor_field_sqnorm_f(s2));
   spinor_field_mul_f(s2,tau,s2);
*/

  
   test_herm(&M,"M");
   /* test_herm(&H,"H"); */

  finalize_process();
}
