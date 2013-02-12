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
static spinor_field_flt *tmp_flt;

void M_dbl(spinor_field *out, spinor_field *in){
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}

void M_flt(spinor_field_flt *out, spinor_field_flt *in){
   g5Dphi_flt(-hmass,tmp_flt,in); 
   g5Dphi_flt(-hmass,out,tmp_flt);
}

spinor_operator D={&D_dbl,NULL}; 
spinor_operator H={&H_dbl,NULL}; 
spinor_operator M={&M_dbl,&M_flt}; 


int main(int argc,char *argv[])
{
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;
   spinor_field_flt *s1_flt, *s2_flt;
   spinor_field_flt *res_flt;


   inverter_par par;

   int cgiters;

   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */
   //   logger_setlevel(0,10000); /* log all */
   logger_setlevel(0,100); /* log all */
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
   u_gauge_flt=alloc_gfield_flt(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
   u_gauge_f_flt=alloc_gfield_f_flt(&glattice);
#endif
   
   lprintf("MAIN",0,"Generating a random gauge field... ");
   fflush(stdout);
   random_u(u_gauge);
   assign_ud2u();
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);
   lprintf("MAIN",0,"done.\n");
   represent_gauge_field();
   represent_gauge_field_flt();
   //   assign_ud2u_f();

   par.err2=1.e-25;
   par.err2_flt=1.e-10;
   par.max_iter=0;
   par.max_iter_flt =0;
   par.kry_dim = 40;

   res=alloc_spinor_field_f(4,&glattice);
   s1=res+1;
   s2=s1+1;
   tmp=s2+1;

   res_flt=alloc_spinor_field_f_flt(4,&glattice);
   s1_flt=res_flt+1;
   s2_flt=s1_flt+1;
   tmp_flt = s2_flt+1;
   
   gaussian_spinor_field(s1);
   assign_sd2s(s1_flt,s1);   
   
   /* TEST GMRES_M */

   lprintf("GMRES TEST",0,"Testing  GMRES\n");
   lprintf("GMRES TEST",100,"---------------------\n");

   cgiters = GMRES(&par, M, s1, res,NULL);
   lprintf("GMRES TEST",0,"Converged in %d iterations\n",cgiters);
   M.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("GMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2);

   free_spinor_field_f(res);



   lprintf("GMRES TEST",0,"Testing  GMRES flt\n");
   lprintf("GMRES TEST",100,"---------------------\n");

   cgiters = GMRES_flt(&par, M, s1_flt, res_flt,NULL);
   lprintf("GMRES TEST",0,"Converged in %d iterations\n",cgiters);
   M.flt(s2_flt,res_flt);
   spinor_field_sub_assign_f_flt(s2_flt,s1_flt);
   tau=spinor_field_sqnorm_f_flt(s2_flt)/spinor_field_sqnorm_f_flt(s1_flt);
   lprintf("GMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2_flt);

   free_spinor_field_f(res);




   finalize_process();

   exit(0);
}
