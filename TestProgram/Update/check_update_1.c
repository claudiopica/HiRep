/******************************************************************************* 
* Check that the molecular dynamics evolution is reversible
* 
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

extern suNg_av_field *momenta;

static void flip_mom()
{

  suNg_algebra_vector  *dptr;
 
  _DECLARE_INT_ITERATOR(ix);

  geometry_descriptor *gd=momenta->type;
  
  int dx;

   _MASTER_FOR(gd,ix) {
     for(dx=0;dx <4 ; dx++)
       {
	 dptr=(suNg_algebra_vector*)(momenta->ptr+4*ix+dx);

	_algebra_vector_mul_g(*dptr,-1.0,*dptr);
	
       }
   }
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
   rlxd_init(1,12345+10*PID);
   
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
   
   
   int rr;
   rhmc_par rpar;
   int_par t_par;
   
   u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
#endif


   rpar.beta = 4.9;
   rpar.mass = -1.3045822102425876010; /* k=0.1855 */
   rpar.mass = -0.8; /* k=0.1855 */
   rpar.nf = 2;
   rpar.MT_prec = 1.e-10;
   rpar.MD_prec = 1.e-10;
   rpar.HB_prec = 1.e-10;
   rpar.force_prec = 1.e-20;
   rpar.n_pf = 2;
   rpar.integrator=&O2MN_multistep;
   rpar.MD_par=&t_par;
   rpar.mshift_solver=&cg_mshift; /* this is not used in the code now */
   
   t_par.tlen = 1.;
   t_par.nsteps = 2;
   t_par.gsteps = 5;
   
   lprintf("MAIN",0,"Generating a random gauge field... ");fflush(stdout);

   random_u(u_gauge);

   lprintf("MAIN",0,"done.\n");

  
   
   project_gauge_field();
   
   represent_gauge_field();
   
   
   init_rhmc(&rpar);

   lprintf("REV TEST",0,"MVM during RHMC initialzation: %ld\n",getMVM());
   lprintf("REV TEST",0,"Initial plaquette: %1.8e\n",avr_plaquette());
   
   rr=update_rhmc_o();
   if(rr<0) {
     lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
     return 1;
   }
   lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());
   
   flip_mom();

   rr=update_rhmc_o();
   if(rr<0) {
     lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
     return 1;
   }
   lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());
   
   free_rhmc();
   
   free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_gfield_f(u_gauge_f);
#endif
   
   finalize_process();
   exit(0);
}
