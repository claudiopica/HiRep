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
static suNg *u_gauge_old;
suNg_algebra_vector *momenta;
suNf_spinor **pf;
rhmc_par _update_par;
rational_app r_S;  /* used for computing the action S in the metropolis test */
rational_app r_MD; /* used in the action MD evolution */
rational_app r_HB;  /* used in pseudofermions heatbath */
double minev, maxev; /* min and max eigenvalue of H^2 */
/* END of State */


/* local action array for metropolis test */
static double *la=0;

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

static int gcd(int a, int b) {
	while (b!=0){
		int t=b;
		b=a%t;
		a=t;
	}
	return a;
}

static void reduce_fraction(int *a, int *b){
	int f;
	f=gcd(abs(*a),abs(*b));
	if (*b!=0 && f!=1){
		*a/=f;
		*b/=f;
	}
}

void init_rhmc(rhmc_par *par){
	int i;

	/* fare test su input par e copiare in _update_par */
	_update_par=*par;

	/* allocate space for the backup copy of gfield */
	u_gauge_old=alloc_gfield();
	suNg_field_copy(u_gauge_old,u_gauge);

	/* allocate momenta */
	momenta = alloc_momenta();

	/* allocate pseudofermions */
	pf=malloc(_update_par.n_pf*sizeof(*pf));
	for(i=0;i<_update_par.n_pf;++i) {
		pf[i]=alloc_spinor_field_f();
	}

	/* allocate memory for the local action */
	if(la==0)
		la=malloc(sizeof(*la)*VOLUME);

	/* represent gauge field and find min and max eigenvalue of H^2 */
	represent_gauge_field();
	find_spec_H2(&maxev,&minev, par->mass); /* find spectral interval of H^2 */

	/* set up rational approx needed for RHMC */
	/* r_S = x^{-Nf/(4*NPf)} is used in the metropolis test */
	r_S.order=1;
	r_S.n=-_update_par.nf;
	r_S.d=4*_update_par.n_pf; 
	reduce_fraction(&r_S.n,&r_S.d);
	r_S.rel_error=_update_par.MT_prec;
	r_app_alloc(&r_S);
	r_app_set(&r_S,minev,maxev);
	/* r_D = x^{-Nf/(2*NPf)} is used in the molecula dynamics */
	r_MD.order=1;
	r_MD.n=-_update_par.nf;
	r_MD.d=2*_update_par.n_pf; 
	reduce_fraction(&r_MD.n,&r_MD.d);
	r_MD.rel_error=_update_par.MD_prec;
	r_app_alloc(&r_MD);
	r_app_set(&r_MD,minev,maxev);
	/* r_D = x^{+Nf/(4*NPf)} is used in the heat bath for pseudofermions */
	r_HB.order=1;
	r_HB.n=_update_par.nf;
	r_HB.d=4*_update_par.n_pf; 
	reduce_fraction(&r_HB.n,&r_HB.d);
	r_HB.rel_error=_update_par.HB_prec;
	r_app_alloc(&r_HB);
	r_app_set(&r_HB,minev,maxev);

}

void free_rhmc(){
	int i;
   /* free momenta */
  free_field(u_gauge_old);
  free_field(momenta);
	for (i=0;i<_update_par.n_pf;++i)
		free_field(pf[i]);
	free(pf);

  if(la!=0) free(la);
  
  r_app_free(&r_S);
  r_app_free(&r_MD);
  r_app_free(&r_HB);
}

int update_rhmc(){

   double deltaH;
   double oldmax,oldmin;
   int i;
   
   /* generate new momenta and pseudofermions */
   gaussian_momenta(momenta);
	 for (i=0;i<_update_par.n_pf;++i)
		 gaussian_spinor_field(pf[i]);

   /* compute starting action */
   local_hmc_action(NEW, la, momenta, pf, pf);
   
   /* compute H2^{a/2}*pf */
	 for (i=0;i<_update_par.n_pf;++i)
		 rational_func(&r_HB, &H2, pf[i], pf[i]);

   /* integrate molecular dynamics */
	 _update_par.integrator(momenta,_update_par.MD_par);
	 /*leapfrog(momenta, _update_par.tlen, _update_par.nsteps);
	 O2MN_multistep(momenta, _update_par.tlen, _update_par.nsteps, 3);*/

	 /* project gauge field */
	 project_gauge_field();
	 represent_gauge_field();

   /* test min and max eigenvalue of H2 and update approx if necessary */
   /* now it just tests the approx !!! */
   oldmax = maxev; /* save old max */
   oldmin = minev; /* save old min */
   find_spec_H2(&maxev,&minev, _update_par.mass); /* find spectral interval of H^2 */
	 r_app_set(&r_S,minev,maxev);
	 r_app_set(&r_MD,minev,maxev);
	 r_app_set(&r_HB,minev,maxev);

   /* compute H2^{-a/2}*pf or H2^{-a}*pf */
   /* here we choose the first strategy which is more symmetric */
	 for (i=0;i<_update_par.n_pf;++i)
		 rational_func(&r_S, &H2, pf[i], pf[i]);

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
				r_app_set(&r_S,minev,maxev);
				r_app_set(&r_MD,minev,maxev);
				r_app_set(&r_HB,minev,maxev);

				return 0;
      }
   }

	 return 1;
}


