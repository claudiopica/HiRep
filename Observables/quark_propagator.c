/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"


spinor_operator loc_H;
static spinor_field *tmpspinor;

static double hmass;
static double hmass_eo;
static void H(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}

static void D(spinor_field *out, spinor_field *in){
  Dphi(hmass,out,in);
}


/*prende spinori lunghi VOLUME/2 !!! */
static void D_pre(spinor_field *out, spinor_field *in){
  Dphi_(OE,tmpspinor,in);
  Dphi_(EO,out,tmpspinor);
  spinor_field_mul_add_assign_f(out,-hmass_eo,in);
}

/*
 * Computes the matrix elements (H^-1)_{x,0}
 * using the QMR inverter with even/odd preconditioning
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void quark_propagator_QMR_eo(FILE *propfile, unsigned int ssite, int nm, double *mass, double acc) {
  mshift_par QMR_par;
  int i;
	int source;
  double *shift;
  spinor_field *in=0;
  spinor_field *resdn=0;
  spinor_field *resd=0;
	spinor_field *res=0;
	int cgiter=0;
	double norm;
  spinor_field *test=0; /*doppia o singola? */

  /* allocate input spinor field */
	set_spinor_len(VOLUME);
  res = alloc_spinor_field_f(1);
  test = alloc_spinor_field_f(1);

  set_spinor_len(VOLUME/2);
	tmpspinor=alloc_spinor_field_f(1);
	resdn=alloc_spinor_field_f(nm);
	resd=alloc_spinor_field_f(nm);
	in = alloc_spinor_field_f(1);

	/* set up inverters parameters */
  shift=(double*)malloc(sizeof(double)*(nm));
  hmass_eo=(4.+mass[0])*(4.+mass[0]); /* we can put any number here!!! */
  for(i=0;i<nm;++i){
    shift[i]=(4.+mass[i])*(4.+mass[i])-hmass_eo;
  }
  QMR_par.n = nm;
  QMR_par.shift = shift;
  QMR_par.err2 = .5*acc;
  QMR_par.max_iter = 0;

  /* noisy background */
	assert(ssite<VOLUME/2);
	gaussian_spinor_field(in);
	for (source=0;source<NF*4*2;++source) {
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);
	
	/* invert noise */
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, in, resdn);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source on an EVEN site */
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 1.;
		cgiter+=g5QMR_mshift(&QMR_par, &D_pre, in, resd);

		for(i=0;i<QMR_par.n;++i){
			spinor_field_sub_f(&resd[i],&resd[i],&resdn[i]); /* compute difference */

			/* compute solution */
			spinor_field_mul_f(res,-(4.+mass[i]),&resd[i]);
			Dphi_(OE,res,&resd[i]);
			if(source&1) ++cgiter; /* count only half of calls. works because the number of sources is even */

			/* this is a test of the solution */
			set_spinor_len(VOLUME);
			hmass=mass[i];
			D(test,res);
			++cgiter;
			( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] -= 1.;
			norm=spinor_field_sqnorm_f(test);
			if(norm>acc)
				lprintf("PROPAGATOR",0,"g5QMR_oe residuum of source [%d] = %e\n",i,norm);

			/* write propagator on file */
			/* multiply by g_5 to match the MINRES version */
			spinor_field_g5_f(res,res);
			/* convert res to single precision */
			assign_sd2s(VOLUME,(suNf_spinor_flt*)_SPINOR_ADDR(res),_SPINOR_ADDR(res));
			error(fwrite(_SPINOR_ADDR(res),(size_t) sizeof(suNf_spinor_flt),(size_t)(VOLUME),propfile)!=(VOLUME),1,"quark_propagator_QMR_eo",
					"Failed to write quark propagator to file");
			set_spinor_len(VOLUME/2);
		}

		/* remove source */
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 0.;
	}

  set_spinor_len(VOLUME);

	lprintf("PROPAGATOR",10,"QMR_eo MVM = %d\n",cgiter);
  
  /* free memory */
  free_spinor_field(in);
	free_spinor_field(tmpspinor);
	free_spinor_field(resdn);
	free_spinor_field(resd);
  free(shift);
	free_spinor_field(test);
	free_spinor_field(res);

}

/*
 * Computes the matrix elements (H^-1)_{x,0}
 * using the QMR inverter
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void quark_propagator_QMR(FILE *propfile, unsigned int ssite, int nm, double *mass, double acc) {
  mshift_par QMR_par;
  int i;
	int source;
  double *shift;
  spinor_field *in=0;
  spinor_field *resdn=0;
  spinor_field *resd=0;
	int cgiter=0;
	double norm;
  spinor_field *test=0;

  /* allocate input spinor field */
	set_spinor_len(VOLUME);
  test = alloc_spinor_field_f(1);
	resdn=alloc_spinor_field_f(nm);
	resd=alloc_spinor_field_f(nm);
	in = alloc_spinor_field_f(1);


	/* set up inverters parameters */
  shift=(double*)malloc(sizeof(double)*(nm));
  hmass=mass[0]; /* we can put any number here!!! */
  for(i=0;i<nm;++i){
    shift[i]=hmass-mass[i];
  }
  QMR_par.n = nm;
  QMR_par.shift = shift;
  QMR_par.err2 = .5*acc;
  QMR_par.max_iter = 0;

  /* noisy background */
	gaussian_spinor_field(in);
	for (source=0;source<NF*4*2;++source) {
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);
	
	/* invert noise */
  cgiter+=g5QMR_mshift(&QMR_par, &D, in, resdn);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source */
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 1.;
		cgiter+=g5QMR_mshift(&QMR_par, &D, in, resd);

		for(i=0;i<QMR_par.n;++i){
			spinor_field_sub_f(&resd[i],&resd[i],&resdn[i]); /* compute difference */

			/* this is a test of the solution */
			D(test,&resd[i]);
			++cgiter;
			spinor_field_mul_add_assign_f(test,-QMR_par.shift[i],&resd[i]);
			( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] -= 1.;
			norm=spinor_field_sqnorm_f(test);
			if(norm>acc)
				lprintf("PROPAGATOR",0,"g5QMR residuum of source [%d] = %e\n",i,norm);

			/* write propagator on file */
			/* multiply by g_5 to match the MINRES version */
			spinor_field_g5_f(&resd[i],&resd[i]);
			/* convert to single precision */
			assign_sd2s(VOLUME,(suNf_spinor_flt*)_SPINOR_ADDR(&resd[i]),_SPINOR_ADDR(&resd[i]));
			error(fwrite(_SPINOR_ADDR(&resd[i]),(size_t) sizeof(suNf_spinor_flt),(size_t)(VOLUME),propfile)!=(VOLUME),1,"quark_propagator_QMR",
					"Failed to write quark propagator to file");
		}

		/* remove source */
		( (double *)_SPINOR_AT_SITE(in,ssite) )[2*source] = 0.;
	}

	lprintf("PROPAGATOR",10,"QMR MVM = %d\n",cgiter);
  
  /* free memory */
	free_spinor_field(in);
	free_spinor_field(resdn);
	free_spinor_field(resd);
  free(shift);
	free_spinor_field(test);

}


/*
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of nm spinor fields
 */
void quark_propagator(unsigned int source, int nm, double *mass, spinor_field *out, double acc) {
  static MINRES_par MINRESpar;
  int i,cgiter;
  spinor_field *in;

  /* allocate input spinor field */
  set_spinor_len(VOLUME);
  in = alloc_spinor_field_f(1);


  /* the source is on the first even site */
  spinor_field_zero_f(in);
  ( (double *)_SPINOR_AT_SITE(in,0) )[2*source] = 1.; /* put in source */

  hmass=mass[0];

  MINRESpar.err2 = acc;
  MINRESpar.max_iter = 0;

	cgiter=0;
  cgiter+=MINRES(&MINRESpar, &H, in, &out[0],0);
  for(i=1;i<nm;++i){
    hmass=mass[i];
    cgiter+=MINRES(&MINRESpar, &H, in, &out[i],&out[i-1]);
  }
	lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);

  /* free input spinor field */
  free_spinor_field(in);

}
