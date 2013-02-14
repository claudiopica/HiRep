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

static double hmass=-0.949000;
//static double hmass=0.11;


void D_dbl(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void D_flt(spinor_field_flt *out, spinor_field_flt *in){
   Dphi_flt(hmass,out,in);
}

void H_dbl(spinor_field *out, spinor_field *in){
   g5Dphi(hmass,out,in);
}

void H_flt(spinor_field_flt *out, spinor_field_flt *in){
   g5Dphi_flt(hmass,out,in);
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

void Mg5_dbl(spinor_field *out, spinor_field *in){
   Dphi(-hmass,tmp,in); 
   Dphi(-hmass,out,tmp);
}

void Mg5_flt(spinor_field_flt *out, spinor_field_flt *in){
   Dphi_flt(-hmass,tmp_flt,in); 
   Dphi_flt(-hmass,out,tmp_flt);
}



spinor_operator D={&D_dbl,&D_flt}; 
spinor_operator H={&H_dbl,&H_flt}; 
spinor_operator M={&M_dbl,&M_flt}; 
spinor_operator Mg5={&M_dbl,&M_flt}; 


inverter_par par_precon;
g5QMR_fltacc_par parg5;
spinor_field_flt *s1_flt, *s2_flt;

void precon_dbl(spinor_field *out, spinor_field *in){
  //  int cgiters;
  assign_sd2s(s1_flt,in);   
  (void) GMRES_flt(&par_precon, H, s1_flt, s2_flt, NULL);
  assign_s2sd(out,s2_flt);
}

void preconQMR_dbl(spinor_field *out, spinor_field *in){
  int cgiters;
  assign_sd2s(s1_flt,in);
  spinor_field_g5_f_flt(s1_flt,s1_flt);
  (void) g5QMR_flt(&par_precon, D, s1_flt, s2_flt);
  spinor_field_g5_f_flt(s2_flt,s2_flt);
  assign_s2sd(out,s2_flt);
}

void precontrivial_dbl(spinor_field *out, spinor_field *in){
  spinor_field_copy_f(out,in);
}

void preconQMRflacc_dbl(spinor_field *out, spinor_field *in){
  int cgiters;
  //  assign_sd2s(s1_flt,in);
  spinor_field_g5_f(tmp,in);
  (void) g5QMR_fltacc(&parg5, D, tmp, out);
  spinor_field_g5_f(out,out);
  //  assign_s2sd(out,s2_flt);
}


spinor_operator precon_GMRES={&precon_dbl,NULL}; 
spinor_operator precon_QMR={&preconQMR_dbl,NULL}; 
spinor_operator precon_QMR_flacc={&preconQMRflacc_dbl,NULL}; 
spinor_operator precon_trivial={&precontrivial_dbl,NULL}; 

int main(int argc,char *argv[])
{
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;
   spinor_field_flt *res_flt;
   gpu_timer t1;
   float elapsed;

   inverter_par par;

   int cgiters;

   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */
   //   logger_setlevel(0,10000); /* log all */
   logger_setlevel(0,100); 
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

   input_gpu gpu_var = init_input_gpu(gpu_var);
   read_input(gpu_var.read,"test_input");
   init_gpu(gpu_var);
   
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

   //   random_u(u_gauge);
   read_gauge_field("conf");

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
   par.kry_dim = 12;

   par_precon.err2=1.e-25;
   par_precon.err2_flt=1.e-6;
   par_precon.max_iter_flt =0;
   par_precon.kry_dim = 10;

   parg5.err2=1.e-25;
   parg5.err2_flt=1.e-4;
   parg5.max_iter=0;
   parg5.max_iter_flt =0;


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
   
   /* TEST FGMRES */

   lprintf("FGMRES TEST",0,"Testing  FGMRES\n");
   lprintf("FGMRES TEST",100,"---------------------\n");

   /*   cgiters=0;
   lprintf("FGMRES TEST",0,"Using flt QMR as preconditioner\n");
   t1 = gpuTimerStart();   
   cgiters = FGMRES(&par, H, s1, res,NULL,precon_QMR_flacc);
   elapsed = gpuTimerStop(t1);
   lprintf("FGMRES TEST",0,"Converged in %d iterations\n",cgiters);
   lprintf("GMRES TEST",0,"Time: %1.10g s\n",elapsed/1000);
   H.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("GMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2);
   */
   /*   cgiters=0;
   lprintf("GMRES TEST",0,"Testing  GMRES\n");
   lprintf("GMRES TEST",100,"---------------------\n");
   t1 = gpuTimerStart();   
   cgiters = GMRES(&par, H, s1, res,NULL);
   elapsed = gpuTimerStop(t1);
   lprintf("GMRES TEST",0,"Converged in %d iterations\n",cgiters);
   lprintf("GMRES TEST",0,"Time: %1.10g s\n",elapsed/1000);
   H.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("GMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2);*/

   lprintf("GMRES TEST",0,"Testing  g5QMR flt\n");
   lprintf("GMRES TEST",100,"---------------------\n");
   t1 = gpuTimerStart();   
   cgiters = g5QMR_fltacc(&parg5, D, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("GMRES TEST",0,"g5QMR converged in %d iterations\n",cgiters);
   lprintf("GMRES TEST",0,"Time: %1.10g s\n",elapsed/1000);
   D.flt(s2_flt,res_flt);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("GMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2);

   cgiters=0;
   lprintf("FGMRES TEST",0,"Using flt GMRES as preconditioner\n");
   t1 = gpuTimerStart();   
   cgiters = FGMRES(&par, H, s1, res,NULL,precon_GMRES);
   elapsed = gpuTimerStop(t1);
   lprintf("FGMRES TEST",0,"Converged in %d iterations\n",cgiters);
   lprintf("FGMRES TEST",0,"Time: %1.10g s\n",elapsed/1000);
   H.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("FGMRES TEST",0,"test  = %e (req. %e)\n",tau,par.err2);



   free_spinor_field_f(res);

   free_spinor_field_f_flt(res_flt);

   finalize_process();

   exit(0);
}
