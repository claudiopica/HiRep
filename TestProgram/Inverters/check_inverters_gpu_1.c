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
#include "gpu.h"

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

static spinor_field_flt *tmp_flt;

void M_flt(spinor_field_flt *out, spinor_field_flt *in){
   g5Dphi_flt(-hmass,tmp_flt,in); 
   g5Dphi_flt(-hmass,out,tmp_flt);
}


spinor_operator D={&D_dbl,NULL}; 
spinor_operator H={&H_dbl,NULL}; 
spinor_operator M={&M_dbl,&M_flt}; 


int main(int argc,char *argv[])
{
   char pame[256];
   int i;
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;

   mshift_par par;
   int cgiters;
   float elapsed;
   gpu_timer t1;
   
   
   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */
   logger_setlevel(0,10000); /* log all */
   //   logger_setlevel(0,10); /* log all */
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
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
   u_gauge_f=alloc_gfield_f(&glattice);
#endif

   u_gauge_flt=alloc_gfield_flt(&glattice);
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
   u_gauge_f_flt=alloc_gfield_f_flt(&glattice);
#endif

   
   lprintf("MAIN",0,"Generating a random gauge field... ");
   fflush(stdout);
   random_u(u_gauge);
   gfield_copy_to_gpu(u_gauge);
   assign_ud2u_cpu();
   gfield_copy_to_gpu_flt(u_gauge_flt);
   lprintf("MAIN",0,"done.\n");
   
   lprintf("MAIN",0,"Representing gauge field... ");
   represent_gauge_field();
   //gfield_copy_to_gpu_f(u_gauge_f);
   represent_gauge_field_flt();
   //gfield_copy_to_gpu_f_flt(u_gauge_f_flt);
   lprintf("MAIN",0,"done.\n");

   par.n = 6;
   par.shift=(double*)malloc(sizeof(double)*(par.n));
   par.err2=1.e-28;
   par.max_iter=0;
   res=alloc_spinor_field_f(par.n+3,&glattice);
   s1=res+par.n;
   s2=s1+1;
   tmp=s2+1;
   
   tmp_flt = alloc_spinor_field_f_flt(1,&glattice);
   
   par.shift[0]=+0.1;
   par.shift[1]=-0.21;
   par.shift[2]=+0.05;
   par.shift[3]=-0.01;
   par.shift[4]=-0.15;
   par.shift[5]=-0.05;
   
   par.n = 6;

   gaussian_spinor_field(s1);
   spinor_field_copy_to_gpu_f(s1);
   
   /* TEST CG_M */

   lprintf("CG TEST",0,"Testing CG multishift\n");
   lprintf("CG TEST",0,"---------------------\n");

   t1 = gpuTimerStart();   
   cgiters = cg_mshift(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);

   lprintf("CG TEST",0,"\n\nTesting CG multishift with single precision acceleration\n");
   lprintf("CG TEST",0,"------------------------------------------------------------\n");
   
   t1 = gpuTimerStart();
   cgiters = cg_mshift_flt(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
     spinor_field_zero_f(&res[i]);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);


   lprintf("CG TEST",0,"\n\nTesting CG multishift version 2 with single precision acceleration\n");
   lprintf("CG TEST",0,"------------------------------------------------------------\n");

   t1 = gpuTimerStart();
   cgiters = cg_mshift_flt2(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
     spinor_field_zero_f(&res[i]);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);

   lprintf("CG TEST",0,"DONE!\n");

   free_spinor_field_f(res);

   free_spinor_field_f_flt(tmp_flt);


   free(par.shift);

   finalize_process();

   exit(0);
}
