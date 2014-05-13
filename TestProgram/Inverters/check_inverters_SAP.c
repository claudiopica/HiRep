/******************************************************************************
*
* Test of Schwarz Alternating Procedure
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


void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void H(spinor_field *out, spinor_field *in){
   g5Dphi(hmass,out,in);
}

static spinor_field *tmp;

void M(spinor_field *out, spinor_field *in){
   empty_buffers(tmp);
   tmp->type=in->type;
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}


int main(int argc,char *argv[])
{
   char pame[256];
   int i;
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;

   mshift_par par;

   int cgiters;

   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */

  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
  }
      
   lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
   
   /* read input file */
   read_input(glb_var.read,"test_input");
   
   
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
   
   init_geometry_SAP();
   /* test_geometry_mpi_eo(); */
    
    /* setup random numbers */
    read_input(rlx_var.read,"test_input");
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

   
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

   lprintf("REV TEST",0,"Initial plaquette: %1.8e\n",avr_plaquette());

   par.n = 1;
   par.shift=(double*)malloc(sizeof(double)*(par.n));
   par.err2=1.e-10;
   par.max_iter=0;
   res=alloc_spinor_field_f(par.n+3,&glattice);
   s1=res+par.n;
   s2=s1+1;
   tmp=s2+1;
   
   par.shift[0]=0.0;
   /*
   lprintf("res->type",0,"\ninner_master_pieces: %d\n",res->type->inner_master_pieces);
   lprintf("res->type",0,"local_master_pieces: %d\n",res->type->local_master_pieces);
   lprintf("res->type",0,"total_spinor_master_pieces: %d\n",res->type->total_spinor_master_pieces);
   lprintf("res->type",0,"total_gauge_master_pieces: %d\n",res->type->total_gauge_master_pieces);
   lprintf("res->type",0,"master_start[0]: %d\n",res->type->master_start[0]);
   lprintf("res->type",0,"master_start[1]: %d\n\n",res->type->master_start[1]);
   */
   
   
   gaussian_spinor_field(s1);
   	
 	g5Dphi(0.1,s2,s1);
   
   /* TEST CG_M */

	
	s1->type=&glat_red;
	s2->type=&glat_red;
	
	
	spinor_field_zero_f(s1);
	
       lprintf("CGTEST",0,"spinor_field_sqnorm_f(s1)=%e\n",spinor_field_sqnorm_f(s1));	// ULRIK
	
 	g5Dphi(0.1,s2,s1);
 	
       lprintf("CGTEST",0,"spinor_field_sqnorm_f(s1)=%e\n",spinor_field_sqnorm_f(s1));	// ULRIK
       lprintf("CGTEST",0,"spinor_field_sqnorm_f(s2)=%e\n",spinor_field_sqnorm_f(s2));	// ULRIK


	s1->type=&glattice;
	s2->type=&glattice;


   lprintf("SAP TEST",0,"Testing SAP\n");
   lprintf("SAP TEST",0,"---------------------\n");

   //cgiters = cg_mshift(&par, &M, s1, res);
   cgiters=0;
   spinor_field_zero_f(res);
   SAP_prec(5,&cg_mshift,&par, &M, s1, res);
   
 //  lprintf("SAP TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("SAP TEST",0,"test SAP = %e (req. %e)\n",tau,par.err2);
   }

   

   free_spinor_field_f(res);
   free(par.shift);

   finalize_process();

   return 0;
}
