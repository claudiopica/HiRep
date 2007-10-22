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

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static suNg oldgauge[4*VOLUME];

int main(int argc,char *argv[])
{
   int rr;
   rhmc_par rpar;
   int_par t_par;

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   
   level=1;
   seed=666;
   rlxd_init(level,seed);

   printf("ranlux: level = %d, seed = %d\n",level,seed); 

   geometry_eo_lexi();
   test_geometry();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif

   rpar.beta = 4.9;
   rpar.mass = -1.3045822102425876010; /* k=0.1855 */
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
   
   printf("Generating a random gauge field... ");fflush(stdout);
   random_u();
   printf("done.\n");
   project_gauge_field();
   represent_gauge_field();

   suNg_field_copy(oldgauge,u_gauge);
   
   init_rhmc(&rpar);

   printf("MVM during RHMC initialzation: %ld\n",getMVM());
   printf("Initial plaquette: %1.8e\n",avr_plaquette());

	 rr=update_rhmc_o();
	 if(rr<0) {
		 printf("Error in updating the gauge field!!\n");
		 return 1;
	 }
   printf("Plaquette: %1.8e\n",avr_plaquette());

   flip_mom();
   
   rr=update_rhmc_o();
   if(rr<0) {
      printf("Error in updating the gauge field!!\n");
      return 1;
   }
   printf("Plaquette: %1.8e\n",avr_plaquette());

   free_rhmc();
   
   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif
   
   exit(0);
}
