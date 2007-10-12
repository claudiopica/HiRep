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
#include "logger.h"

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

static short int init=0;

/* local action array for metropolis test */
static double *la=0;

/* this is the basic operator used in the update */
static suNf_spinor *h2tmp;
void H2(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(_update_par.mass, h2tmp, in);
  g5Dphi(_update_par.mass, out, h2tmp);
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
	unsigned int len;

	lprintf("RHMC",0,"Initializing...\n");

	/* fare test su input par e copiare in _update_par */
	_update_par=*par;

	lprintf("RHMC",10,
			"Number of Flavors = %d\n"
			"beta = %.8f\n"
			"Mass = %.8f\n"
			"Metropolis test precision = %.8e\n"
			"MD precision = %.8e\n"
			"PF heat-bath precision = %.8e\n"
			"RHMC force precision = %.8e\n"
			"Number of pseudofermions = %d\n"
			"MD trajectory length = %.8f\n"
			"MD steps = %d\n"
			"MD gauge substeps = %d\n"
			,_update_par.nf
			,_update_par.beta
			,_update_par.mass
			,_update_par.MT_prec
			,_update_par.MD_prec
			,_update_par.HB_prec
			,_update_par.force_prec
			,_update_par.n_pf
			,_update_par.MD_par->tlen
			,_update_par.MD_par->nsteps
			,_update_par.MD_par->gsteps
			);

	/* allocate space for the backup copy of gfield */
	u_gauge_old=alloc_gfield();
	suNg_field_copy(u_gauge_old,u_gauge);

	/* allocate h2tmp for H2 */
  h2tmp=alloc_spinor_field_f(1);

	/* allocate momenta */
	momenta = alloc_momenta();

	/* allocate pseudofermions */
	get_spinor_len(&len);
	pf=malloc(_update_par.n_pf*sizeof(*pf));
	pf[0]=alloc_spinor_field_f(_update_par.n_pf);
	for(i=1;i<_update_par.n_pf;++i) {
		pf[i]=pf[i-1]+len;
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

	init = 1;

	lprintf("RHMC",0,"Initialization done.\n");

}

void free_rhmc(){
	/* free momenta */
  free_field(u_gauge_old);
  free_field(momenta);
	free_field(h2tmp);
	free_field(pf[0]);
	free(pf);

  if(la!=0) free(la);
  
  r_app_free(&r_S);
  r_app_free(&r_MD);
  r_app_free(&r_HB);

	init = 0;

	lprintf("RHMC",0,"Memory deallocated.\n");

}

int update_rhmc(){

   double deltaH;
   double oldmax,oldmin;
   int i;
  
	 if(!init)
		 return -1;

   /* generate new momenta and pseudofermions */
	 lprintf("RHMC",30,"Generating gaussian momenta and pseudofermions...\n");
   gaussian_momenta(momenta);
	 for (i=0;i<_update_par.n_pf;++i)
		 gaussian_spinor_field(pf[i]);

   /* compute starting action */
	 lprintf("RHMC",30,"Computing action density...\n");
   local_hmc_action(NEW, la, momenta, pf, pf);
   
   /* compute H2^{a/2}*pf */
	 lprintf("RHMC",30,"Correcting pseudofermions distribution...\n");
	 for (i=0;i<_update_par.n_pf;++i)
		 rational_func(&r_HB, &H2, pf[i], pf[i]);

   /* integrate molecular dynamics */
	 lprintf("RHMC",30,"MD integration...\n");
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

	 lprintf("RHMC",30,"Computing new action density...\n");
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
   lprintf("RHMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));

   if(deltaH<0.) {
      suNg_field_copy(u_gauge_old,u_gauge);
   } else {
      double r;
      ranlxd(&r,1);
      if(r<exp(-deltaH)) {
				suNg_field_copy(u_gauge_old,u_gauge);
      } else {
				lprintf("RHMC",10,"Configuration rejected.\n");
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

	 lprintf("RHMC",10,"Configuration accepted.\n");

	 return 1;
}


