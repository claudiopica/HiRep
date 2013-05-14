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
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "update.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=-0.949000;
//static double hmass=-1.0;

// ------------------ Building operators ----------
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
   tmp->type=&glattice;
   empty_buffers(tmp);
   tmp->type=in->type;
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}

void M_flt(spinor_field_flt *out, spinor_field_flt *in){
   tmp->type=&glattice;
   empty_buffers(tmp);
   tmp->type=in->type;
   g5Dphi_flt(-hmass,tmp_flt,in); 
   g5Dphi_flt(-hmass,out,tmp_flt);
}

void Mg5_dbl(spinor_field *out, spinor_field *in){
   tmp->type=&glattice;
   empty_buffers(tmp);
   tmp->type=in->type;
   Dphi(-hmass,tmp,in); 
   Dphi(-hmass,out,tmp);
}

void Mg5_flt(spinor_field_flt *out, spinor_field_flt *in){
   tmp->type=&glattice;
   empty_buffers(tmp);
   tmp->type=in->type;
   Dphi_flt(-hmass,tmp_flt,in); 
   Dphi_flt(-hmass,out,tmp_flt);
}

spinor_operator D={&D_dbl,&D_flt}; 
spinor_operator Hop={&H_dbl,&H_flt}; 
spinor_operator M={&M_dbl,&M_flt}; 
spinor_operator Mg5={&M_dbl,&M_flt}; 


// ------------------ Defining preconditioners ----------

mshift_par par_precon;
mshift_par par_cg;

void precon_dbl(spinor_field *out, spinor_field *in){
  DDalphaAMG(&par_precon, Hop, in, out);
}
void precon_SAP_dbl(spinor_field *out, spinor_field *in){
  SAP_prec(3,&cg_mshift,&par_precon, Hop, in, out);
}
void precontrivial_dbl(spinor_field *out, spinor_field *in){
  spinor_field_copy_f(out,in);
}

spinor_operator precon_DDalphaAMG={&precon_dbl,NULL}; 
spinor_operator precon_SAP={&precon_SAP_dbl,NULL}; 
spinor_operator precon_trivial={&precontrivial_dbl,NULL}; 



// ------------------ Main function --------------
int main(int argc,char *argv[])
{
// ------------------ Standard stuff ----------
   char pame[256];
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;
   struct timeval start, end, etime;

   inverter_par par;

   int cgiters;

   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */
  logger_setlevel(0,5); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
//    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
  }
   
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
   lprintf("MAIN",100,"Logger test\n");
   
   /* setup lattice geometry */
   geometry_mpi_eo();
   /* test_geometry_mpi_eo(); */
   
   
// ------------------ Initializing new geometries ----------
   
    init_geometry_SAP();
    init_geometry_nocomm();
  
  
  
   lprintf("MAIN",100,"Logger test 0.1\n");
   
   u_gauge=alloc_gfield(&glattice);
   u_gauge_flt=alloc_gfield_flt(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
   u_gauge_f_flt=alloc_gfield_f_flt(&glattice);
#endif

   lprintf("MAIN",0,"Generating a random gauge field... ");
   fflush(stdout);

     random_u(u_gauge);
//   read_gauge_field("run1_12x12x12x12nc3rSYMnf2b6.000000m1.150000n523");

   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);
   lprintf("MAIN",0,"done.\n");
   represent_gauge_field();
   
   
   lprintf("MAIN",100,"Logger test 1\n");
   
  // ------------------ Inverter parameters (For the FGMRES) ---------- 
   par.err2=1.e-28;
   par.err2_flt=1.e-10;
   par.max_iter=0;
   par.max_iter_flt =0;
   par.kry_dim = 30;

  // ------------------ Inverter parameters (For the cg inside) ---------- 
   par_precon.n=1;
   par_precon.shift=(double*)malloc(sizeof(double)*(par_precon.n));
   par_precon.err2=1.e-2;
   par_precon.max_iter=0;
   par_precon.shift[0]=0.0;
   
  // ------------------ Inverter parameters (For CG comparison) ---------- 
   par_cg.n=1;
   par_cg.shift=(double*)malloc(sizeof(double)*(par_precon.n));
   par_cg.err2=1.e-28;
   par_cg.max_iter=0;
   par_cg.shift[0]=0.0;
   
   res=alloc_spinor_field_f(4,&glattice);
   s1=res+1;
   s2=s1+1;
   tmp=s2+1;


      
   gaussian_spinor_field(s1);
   
   
   
   
// ------------------ FGMRES with SAP test ---------- 

   spinor_field_zero_f(res);
   
   lprintf("FGMRES with SAP",0,"\n Testing  FGMRES with SAP preconditioner\n");
   lprintf("FGMRES with SAP",0,"---------------------\n");

   cgiters=0;
   gettimeofday(&start,0);
   cgiters = FGMRES(&par, Hop, s1, res,NULL,precon_SAP);
   gettimeofday(&end,0);
   timeval_subtract(&etime,&end,&start);
   
   lprintf("FGMRES with SAP",0," Converged in %d iterations\n",cgiters);
   Hop.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("FGMRES with SAP",0," Normalized residual  = %e (req. %e)\n",tau,par.err2);
   lprintf("FGMRES with SAP",0," Inversion took: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
   lprintf("FGMRES with SAP",0,"---------------------\n");   
   

   
// ------------------ FGMRES with DDalphaAMG test ---------- 

   spinor_field_zero_f(res);
   
   DDalphaAMG_setup(&par_precon, Hop, 5, 3, 6);
   // Arguments are (cg_parameters, operator, Number of vectors to span near kernel, (not used), n_inv)
   
   lprintf("FGMRES with DDalphaAMG",0,"\n Testing  FGMRES with DDalphaAMG preconditioner\n");
   lprintf("FGMRES with DDalphaAMG",0,"---------------------\n");

   cgiters=0;
   gettimeofday(&start,0);   
   cgiters = FGMRES(&par, Hop, s1, res,NULL,precon_DDalphaAMG);
   gettimeofday(&end,0);
   timeval_subtract(&etime,&end,&start);
   
   lprintf("FGMRES with DDalphaAMG",0," Converged in %d iterations\n",cgiters);
   Hop.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("FGMRES with DDalphaAMG",0," Normalized residual  = %e (req. %e)\n",tau,par.err2);
   lprintf("FGMRES with DDalphaAMG",0," Inversion took: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
   lprintf("FGMRES with DDalphaAMG",0,"---------------------\n");   

   

// ------------------ FGMRES with trivial test ---------- 

   spinor_field_zero_f(res);
   
   lprintf("FGMRES with trivial",0,"\n Testing  FGMRES with trivial preconditioner\n");
   lprintf("FGMRES with trivial",0,"---------------------\n");

   cgiters=0;
   gettimeofday(&start,0);   
   cgiters = FGMRES(&par, Hop, s1, res,NULL,precon_trivial);
   gettimeofday(&end,0);
   timeval_subtract(&etime,&end,&start);
   
 // cgiters = FGMRES(&par, Hop, s1, res,NULL,precon_trivial);
   lprintf("FGMRES with trivial",0," Converged in %d iterations\n",cgiters);
   Hop.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("FGMRES with trivial",0," Normalized residual  = %e (req. %e)\n",tau,par.err2);
   lprintf("FGMRES with trivial",0," Inversion took: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
   lprintf("FGMRES with trivial",0,"---------------------\n");
   
      

// ------------------ CG inverter ---------- 

   spinor_field_zero_f(res);
   
   lprintf("CG inverter",0,"\n Testing  CG inverter\n");
   lprintf("CG inverter",0,"---------------------\n");

   cgiters=0;
   gettimeofday(&start,0);   
   cgiters = cg_mshift(&par_cg, Hop, s1, res);
   gettimeofday(&end,0);
   timeval_subtract(&etime,&end,&start);
   
 // cgiters = FGMRES(&par, Hop, s1, res,NULL,precon_trivial);
   lprintf("CG inverter",0," Converged in %d iterations\n",cgiters);
   Hop.dbl(s2,res);
   spinor_field_sub_assign_f(s2,s1);
   tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
   lprintf("CG inverter",0," Normalized residual  = %e (req. %e)\n",tau,par.err2);
   lprintf("CG inverter",0," Inversion took: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
   lprintf("CG inverter",0,"---------------------\n");






   free_spinor_field_f(res);


   finalize_process();

   exit(0);
}
