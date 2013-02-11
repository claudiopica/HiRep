/******************************************************************************
*
* Test of modules
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

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;


void D_dbl(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void H_dbl(spinor_field *out, spinor_field *in){
   g5Dphi(hmass,out,in);
}

static spinor_field *tmp;

void M_dbl(spinor_field *out, spinor_field *in){
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}

spinor_operator D={&D_dbl,NULL}; 
spinor_operator H={&H_dbl,NULL}; 
spinor_operator M={&M_dbl,NULL}; 


int main(int argc,char *argv[])
{
   int i;
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;


   MINRES_par par;

   int cgiters;

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
   
   u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
#endif
   
   lprintf("MAIN",0,"Generating a random gauge field... ");
   fflush(stdout);
   random_u(u_gauge);
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);
   lprintf("MAIN",0,"done.\n");
   represent_gauge_field();
   
   par.err2=1.e-28;
   par.max_iter=0;
   res=alloc_spinor_field_f(4,&glattice);
   s1=res+1;
   s2=s1+1;
   tmp=s2+1;
   
   gaussian_spinor_field(s1);
   
   
   /* TEST CG_M */

   lprintf("CG TEST",0,"Testing  multishift\n");
   lprintf("CG TEST",100,"---------------------\n");

   cgiters = GMRES(&par, M, s1, res,NULL);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   M.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("CG TEST",0,"test  = %e (req. %e)\n",tau,par.err2);

   free_spinor_field_f(res);

   finalize_process();

   exit(0);
}
