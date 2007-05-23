#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "rational_functions.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* State quantities for RHMC */
suNg *u_gauge_old;
suNg_algebra_vector *momenta;
suNf_spinor *pf;
rhmc_par _update_par;
rational_app r_S;  /* used for computing the action S in the metropolis test */
rational_app r_MD; /* used in the action MD evolution */
rational_app r_HB;  /* used in pseudofermions heatbath */
double minev, maxev; /* min and max eigenvalue of H^2 */
/* END of State */

static float *la;

/* this is the basic operator used in the update */
void H2(suNf_spinor *out, suNf_spinor *in){
  suNf_spinor *tmp=(suNf_spinor*)malloc(sizeof(suNf_spinor)*VOLUME);
  g5Dphi(_update_par.mass, tmp, in);
  g5Dphi(_update_par.mass, out, tmp);
  free(tmp);
}

/* this is the basic operator used in the update */
/*
void H(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(_update_par.mass, out, in);
}
*/

void init_rhmc(rhmc_par *par){
   /* fare test su input par e copiare in _update_par */
   _update_par=*par;

   /* allocate space for the backup copy of gfield */
   u_gauge_old=alloc_gfield();
   suNg_field_copy(u_gauge_old,u_gauge);

   /* allocate momenta */
   momenta = alloc_momenta();

   /* if (par->nf!=0) { pure gauge */ /* not the pure gauge */
     /* allocate pseudofermions */
     pf=alloc_spinor_field_f();

     /* allocate memory for the local action */
     la=(float*)malloc(sizeof(float)*VOLUME);
     
     /* represent gauge field and find min and max eigenvalue of H^2 */
     represent_gauge_field();
     find_spec_H2(&maxev,&minev, par->mass); /* find spectral interval of H^2 */
     
     /* set up rational approx needed for RHMC */
     r_S.order=15;
     r_S.a=(double*)malloc(sizeof(double)*(r_S.order+1));
     r_S.b=(double*)malloc(sizeof(double)*r_S.order);
     inv_fourrt_coef(&r_S, minev, maxev); /* for nf==1 */
     
     r_MD.order=7;
     r_MD.a=(double*)malloc(sizeof(double)*(r_MD.order+1));
     r_MD.b=(double*)malloc(sizeof(double)*r_MD.order);
     inv_sqrt_coef(&r_MD, minev, maxev); /* for nf==1 */
     
     r_HB.order=15;
     r_HB.a=(double*)malloc(sizeof(double)*(r_HB.order+1));
     r_HB.b=(double*)malloc(sizeof(double)*r_HB.order);
     fourrt_coef(&r_HB, minev, maxev); /* for nf==1 */
     /*} pure gauge */
}

void free_rhmc(){
   /* free momenta */
  free_field(u_gauge_old);
  free_field(momenta);
  /*  if (_update_par.nf!=0) { pure gauge */
  free_field(pf);

  free(la);
  
  free(r_S.a); free(r_S.b);
  free(r_MD.a); free(r_MD.b);
  free(r_HB.a); free(r_HB.b);
  /*} pure gauge*/
}

int update_rhmc(){

   double deltaH;
   double oldmax,oldmin;
   int i;
   
   /* generate new momenta and pseudofermions */
   gaussian_momenta(momenta);
   gaussian_spinor_field(pf);

   /* compute starting action */
   local_hmc_action(NEW, la, momenta, pf, pf);
   
   /* compute H2^{a/2}*pf */
   rational_func(&r_HB, &H2, pf, pf);

   /* integrate molecular dynamics */
   /* leapfrog(momenta, _update_par.tlen, _update_par.nsteps); */
   O2MN_multistep(momenta, _update_par.tlen, _update_par.nsteps, 3); 

	 /* project gauge field */
	 project_gauge_field();
	 represent_gauge_field();

   /* test min and max eigenvalue of H2 and update approx if necessary */
   /* now it just tests the approx !!! */
   oldmax = maxev; /* save old max */
   oldmin = minev; /* save old min */
   find_spec_H2(&maxev,&minev, _update_par.mass); /* find spectral interval of H^2 */
   inv_fourrt_coef(&r_S, minev, maxev); /* for nf==1 */
   inv_sqrt_coef(&r_MD, minev, maxev); /* for nf==1 */
   fourrt_coef(&r_HB, minev, maxev); /* for nf==1 */

   /* compute H2^{-a/2}*pf oppure H2^{-a}*pf */
   /* here we choose the first strategy which is more symmetric */
   rational_func(&r_S, &H2, pf, pf);

   /* compute new action */
   local_hmc_action(DELTA, la, momenta, pf, pf);

   /* Metropolis test */
   deltaH=0.;
   for(i=0; i<VOLUME; ++i) {
      deltaH+=la[i];
   }
   printf("[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));

   if(deltaH<0.) {
      suNg_field_copy(u_gauge_old,u_gauge);
   } else {
      float r;
      ranlxs(&r,1);
      if(r<exp(-deltaH)) {
				suNg_field_copy(u_gauge_old,u_gauge);
      } else {
				printf("NON Accettato\n");
				suNg_field_copy(u_gauge,u_gauge_old);
				represent_gauge_field();

				/* revert the approx to the old one */
				maxev=oldmax;
				minev=oldmin;
				inv_fourrt_coef(&r_S, minev, maxev); /* for nf==1 */
				inv_sqrt_coef(&r_MD, minev, maxev); /* for nf==1 */
				fourrt_coef(&r_HB, minev, maxev); /* for nf==1 */

				return 0;
      }
   }

	 return 1;
}


