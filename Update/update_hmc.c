#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "error.h"
#include "representation.h"
#include "rational_functions.h"
#include "linear_algebra.h"
#include "inverters.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"

/* State quantities for RHMC declared in fermions_update_common.c */
extern suNg_algebra_vector *momenta;
extern suNf_spinor **pf;
extern rhmc_par _update_par;
extern rational_app r_S;  /* used for computing the action S in the metropolis test */
extern rational_app r_MD; /* used in the action MD evolution */
extern rational_app r_HB;  /* used in pseudofermions heatbath */
extern double minev, maxev; /* min and max eigenvalue of H^2 */
/* END of State */

static short int init=0;

/* local action array for metropolis test */
static double *la=NULL;
/* gauge field copy */
static suNg *u_gauge_old=NULL;
static MINRES_par pfa;

void init_hmc(rhmc_par *par){
	int i;
	unsigned int len;

	if (init) /* already initialized */
		return;

	init_fermions_common();

	lprintf("HMC",0,"Initializing...\n");

#ifdef UPDATE_EO
	set_spinor_len(VOLUME/2);
#else
	set_spinor_len(VOLUME);
#endif

	/* fare test su input par e copiare in _update_par */
	_update_par=*par;
	if (_update_par.nf!=2*_update_par.n_pf) {
		lprintf("HMC",0,"The number of fermions is not twice the number of pseudofermions.\nTry with the RHMC algorithm\n");
		error(1,1,"init_hmc",
				"The HMC algorithm is not suitable for the parameters specified as input\n");
	}

	lprintf("HMC",10,
			"Number of Flavors = %d\n"
			"beta = %.8f\n"
			"Mass = %.8f\n"
			"Metropolis test precision = %.8e\n"
			"RHMC force precision = %.8e\n"
			"Number of pseudofermions = %d\n"
			"MD trajectory length = %.8f\n"
			"MD steps = %d\n"
			"MD gauge substeps = %d\n"
			,_update_par.nf
			,_update_par.beta
			,_update_par.mass
			,_update_par.MT_prec
			,_update_par.force_prec
			,_update_par.n_pf
			,_update_par.MD_par->tlen
			,_update_par.MD_par->nsteps
			,_update_par.MD_par->gsteps
			);

	/* allocate space for the backup copy of gfield */
	u_gauge_old=alloc_gfield();
	suNg_field_copy(u_gauge_old,u_gauge);

	/* allocate momenta */
	momenta = alloc_momenta();

	/* allocate pseudofermions */
	/* we allocate one more pseudofermion for the computation 
	 * of the final action density 
	 */
	get_spinor_len(&len);
	pf=malloc((_update_par.n_pf+1)*sizeof(*pf));
	pf[0]=alloc_spinor_field_f(_update_par.n_pf+1);
	for(i=1;i<_update_par.n_pf+1;++i) {
		pf[i]=pf[i-1]+len;
	}

	/* allocate memory for the local action */
	if(la==NULL)
		la=malloc(sizeof(*la)*VOLUME);

	/* represent gauge field and find min and max eigenvalue of H^2 */
	represent_gauge_field();

	/* set up rational approx needed for HMC */
	/* r_S = x^{-Nf/(4*NPf)} = x^-1/2 is used in the metropolis test */
	/* since H2^-1/2 = H^-1 we use the MINRES inverter */
	pfa.err2=_update_par.MT_prec;
	pfa.err2*=pfa.err2;
	pfa.max_iter=0;
	/* r_MD = x^{-Nf/(2*NPf)} = x^-1 is used in the molecular dynamics in the force*/
	r_MD.order=1;
	r_MD.n=-1;
	r_MD.d=1; 
	r_MD.rel_error=0.; /* not used */
	r_app_alloc(&r_MD);
	r_MD.a[0]=r_MD.b[0]=0.;
	r_MD.a[1]=1.;
	/* r_HB = x^{+Nf/(4*NPf)} = x^1/2 is used in the heat bath for pseudofermions */
	/* we don't need rational approximation for this since H2^1/2 = H */

	init = 1;

	lprintf("RHMC",0,"Initialization done.\n");

}

void free_hmc(){

	if (!init) /* not initialized */
		return;

	/* free momenta */
  free_field(u_gauge_old); u_gauge_old=NULL;
  free_field(momenta); momenta=NULL;
	free_field(pf[0]);
	free(pf); pf=NULL;

  if(la!=NULL) free(la); la=NULL;
  
  r_app_free(&r_S);
  r_app_free(&r_MD);
  r_app_free(&r_HB);

	init = 0;

	lprintf("RHMC",0,"Memory deallocated.\n");

}

int update_hmc(){

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
   
   /* compute H2^{1/2}*pf = H*pf */
	 lprintf("RHMC",30,"Correcting pseudofermions distribution...\n");
	 for (i=0;i<_update_par.n_pf;++i)
		 H(pf[i], pf[i]);

   /* integrate molecular dynamics */
	 lprintf("RHMC",30,"MD integration...\n");
	 _update_par.integrator(momenta,_update_par.MD_par);

	 /* project gauge field */
	 project_gauge_field();
	 represent_gauge_field();

	 lprintf("RHMC",30,"Computing new action density...\n");
   /* compute H2^{-1/2}*pf or H2^{-1}*pf */
   /* here we choose the first strategy which is more symmetric */
	 /* for the HMC H2^-1/2 = H^-1 and we use MINRES for this inversion */
	 for (i=0;i<_update_par.n_pf;++i) {
		 spinor_field_copy_f(pf[_update_par.n_pf],pf[i]);
		 MINRES(&pfa,&H,pf[_update_par.n_pf],pf[i],0);
	 }

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

				return 0;
      }
   }

	 lprintf("RHMC",10,"Configuration accepted.\n");

	 return 1;
}


